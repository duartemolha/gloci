//! Nearest-neighbour candidates (Task 6): for each A and B file, the `k`
//! nearest non-overlapping B on each side, plus ties at the k-th distance,
//! within `max_distance`, under the strand rule (D55). Checked against a
//! brute-force reference over BED inputs and over BED/VCF inputs whose VCF
//! spans start after their key.

use gloci_core::{ChromTable, Interval, Strand, StrandRule, VcfSpan, Window, distance};
use gloci_engine::{
    ARec, Hit, Inputs, NearestConfig, Neighbours, SweepConfig, SweepSummary, Visitor, sweep_sorted,
};
use gloci_io::{InvalidPolicy, IoError, ReaderOptions};
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, Sides, naive_hits, naive_neighbours};
use gloci_testkit::strategies::{
    CHROMS, Feature, HUGE, MixedInputs, RelRec, VCF_HEADER, VcfShape, relation_inputs, vcf_rec,
    vcf_relation_inputs,
};
use proptest::prelude::*;
use proptest::test_runner::TestCaseError;
use std::num::NonZeroUsize;

fn opts() -> ReaderOptions {
    ReaderOptions {
        chrom_map: None,
        threads: NonZeroUsize::new(1).unwrap(),
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn text(b: &[u8]) -> String {
    String::from_utf8(b.to_vec()).unwrap()
}

/// A neighbour as the visitor saw it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Nb {
    raw: String,
    span: Interval,
    strand: Strand,
    loc: u64,
}

/// Each A (raw line and span), its hits as `(file, B line)`, and for each
/// of the `files` B inputs its `(left, right)` neighbours.
#[derive(Debug, Default)]
struct Near {
    files: u16,
    a: Vec<String>,
    a_span: Vec<Interval>,
    hits: Vec<Vec<(u16, String)>>,
    near: Vec<Vec<(Vec<Nb>, Vec<Nb>)>>,
}

impl Visitor for Near {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        let side = |f: u16, hs: &[Hit<'_>]| -> Vec<Nb> {
            hs.iter()
                .map(|h| {
                    assert_eq!(h.file, f, "a neighbour of file {f} tagged {}", h.file);
                    Nb {
                        raw: text(h.raw),
                        span: h.span,
                        strand: h.strand,
                        loc: h.loc,
                    }
                })
                .collect()
        };
        self.a.push(text(a.raw));
        self.a_span.push(a.meta.span);
        self.hits
            .push(hits.iter().map(|h| (h.file, text(h.raw))).collect());
        self.near.push(
            (0..self.files)
                .map(|f| (side(f, near.left(f)), side(f, near.right(f))))
                .collect(),
        );
        Ok(())
    }
}

/// Sweep in-memory A text against B texts (named `a`, `b0`, `b1`, …) with
/// a karyotype table.
fn run(a: &str, b: &[&str], cfg: &SweepConfig) -> Result<(Near, SweepSummary), IoError> {
    let o = opts();
    let inputs = Inputs {
        a: mem_reader("a", a, &o),
        b: b.iter()
            .enumerate()
            .map(|(i, t)| mem_reader(&format!("b{i}"), t, &o))
            .collect(),
    };
    let mut chroms = ChromTable::karyotype();
    let mut near = Near {
        files: b.len() as u16,
        ..Near::default()
    };
    let summary = sweep_sorted(inputs, &mut chroms, cfg, &mut near)?;
    Ok((near, summary))
}

fn nearest(k: usize, max_distance: Option<u64>, strand: StrandRule) -> SweepConfig {
    SweepConfig {
        window: None,
        inclusive_points: false,
        nearest: Some(NearestConfig {
            k,
            max_distance,
            strand,
            split_by_b_strand: false,
        }),
    }
}

/// `nearest` with `split_by_b_strand`.
fn split(k: usize, strand: StrandRule) -> SweepConfig {
    let mut cfg = nearest(k, None, strand);
    if let Some(n) = &mut cfg.nearest {
        n.split_by_b_strand = true;
    }
    cfg
}

fn raws(v: &[Nb]) -> Vec<&str> {
    v.iter().map(|n| n.raw.as_str()).collect()
}

fn dists(a: Interval, v: &[Nb]) -> Vec<u64> {
    v.iter().map(|n| distance(a, n.span, false)).collect()
}

fn pad() -> impl Strategy<Value = u64> {
    prop_oneof![
        4 => 0u64..4,
        2 => 4u64..100,
        1 => 100u64..100_000,
        1 => Just(HUGE),
        1 => Just(u64::MAX),
    ]
}

fn window() -> impl Strategy<Value = Option<Window>> {
    prop_oneof![
        2 => Just(None),
        1 => (pad(), pad(), any::<bool>()).prop_map(|(left, right, strand_aware)| Some(Window {
            left,
            right,
            strand_aware,
        })),
    ]
}

fn max_distance() -> impl Strategy<Value = u64> {
    prop_oneof![
        2 => 0u64..4,
        2 => 4u64..100,
        1 => 100u64..100_000,
        1 => Just(HUGE),
    ]
}

/// Sort `input`, then sweep it with `k`, every strand rule with
/// `split_by_b_strand` off and on, `--inclusive-points` off and on, and
/// `max_distance` off and at `max_d`.
/// Each A's neighbours must be exactly the naive ones (the same records,
/// in the same order, bytes intact), and its hits unchanged.
fn check_neighbours(
    mut input: MixedInputs,
    k: usize,
    max_d: u64,
    window: Option<Window>,
) -> Result<(), TestCaseError> {
    let rank = naive::karyotype_rank(CHROMS);
    input.a.sort(&rank);
    for f in &mut input.b {
        f.sort(&rank);
    }
    let a = &input.a.recs;
    let b: Vec<Vec<RelRec>> = input.b.iter().map(|f| f.recs.clone()).collect();
    let a_text = input.a.text();
    let a_lines: Vec<String> = a.iter().map(Feature::line).collect();
    let b_texts: Vec<String> = input.b.iter().map(|f| f.text()).collect();
    let b_refs: Vec<&str> = b_texts.iter().map(String::as_str).collect();
    let rules = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite];
    for (strand, split_by_b_strand) in rules.into_iter().flat_map(|r| [(r, false), (r, true)]) {
        for inclusive_points in [false, true] {
            for max_distance in [None, Some(max_d)] {
                let cfg = SweepConfig {
                    window,
                    inclusive_points,
                    nearest: Some(NearestConfig {
                        k,
                        max_distance,
                        strand,
                        split_by_b_strand,
                    }),
                };
                let (seen, _) = run(&a_text, &b_refs, &cfg).unwrap();
                prop_assert_eq!(&seen.a, &a_lines);
                prop_assert_eq!(&seen.hits, &naive_hits(a, &b, window, inclusive_points));
                let mut got: Vec<Vec<Sides>> = Vec::new();
                for files in &seen.near {
                    let mut per_file = Vec::new();
                    for (f, (left, right)) in files.iter().enumerate() {
                        let index = |v: &[Nb]| -> Result<Vec<usize>, TestCaseError> {
                            let mut out = Vec::new();
                            for nb in v {
                                let i = nb.loc - 1 - input.b[f].header_lines();
                                let rec = &b[f][i as usize];
                                prop_assert_eq!(&nb.raw, &rec.line());
                                prop_assert_eq!(nb.span, rec.span());
                                prop_assert_eq!(nb.strand, rec.strand());
                                out.push(i as usize);
                            }
                            Ok(out)
                        };
                        per_file.push((index(left)?, index(right)?));
                    }
                    got.push(per_file);
                }
                prop_assert_eq!(
                    &got,
                    &naive_neighbours(
                        a,
                        &b,
                        k,
                        max_distance,
                        strand,
                        split_by_b_strand,
                        inclusive_points
                    ),
                    "{:?}",
                    cfg
                );
            }
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// BED inputs, where every key is the span start.
    #[test]
    fn neighbours_match_naive(
        input in relation_inputs(),
        k in 1usize..=3,
        max_d in max_distance(),
        window in window(),
    ) {
        check_neighbours(MixedInputs::from(input), k, max_d, window)?;
    }

    /// BED and VCF inputs, where VCF spans (of A and of B) can start after
    /// the key: the left pool must go by A's key, and reading ahead must
    /// not stop at a B whose start is beyond the k-th while a later key
    /// can still start before it.
    #[test]
    fn neighbours_match_naive_vcf_style(
        input in vcf_relation_inputs(),
        k in 1usize..=3,
        max_d in max_distance(),
        window in window(),
    ) {
        check_neighbours(input, k, max_d, window)?;
    }
}

#[test]
fn edge_cases() {
    let a_iv = Interval {
        start: 100,
        end: 110,
    };
    let a = "chr1\t100\t110\ta\t0\t+\n";
    let k1 = nearest(1, None, StrandRule::Any);

    // No B left of A: the right side still has the nearest.
    let b = "chr1\t120\t130\nchr1\t150\t160\n";
    let (seen, _) = run(a, &[b], &k1).unwrap();
    let (left, right) = &seen.near[0][0];
    assert!(left.is_empty());
    assert_eq!(raws(right), ["chr1\t120\t130"]);
    assert_eq!(dists(a_iv, right), [11]);

    // Only B on another chromosome: both sides are empty.
    let (seen, _) = run(a, &["chr2\t10\t20\nchr2\t200\t210\n"], &k1).unwrap();
    assert_eq!(seen.near[0][0], (vec![], vec![]));

    // Five B at the same distance on each side with k = 1: all five, in
    // input order; the farther ones are not there.
    let b = "chr1\t40\t89\tfar\n\
             chr1\t50\t90\tl1\nchr1\t60\t90\tl2\nchr1\t70\t90\tl3\nchr1\t80\t90\tl4\nchr1\t85\t90\tl5\n\
             chr1\t120\t125\tr1\nchr1\t120\t126\tr2\nchr1\t120\t127\tr3\nchr1\t120\t128\tr4\nchr1\t120\t129\tr5\n\
             chr1\t121\t122\tfar\n";
    let (seen, _) = run(a, &[b], &k1).unwrap();
    let (left, right) = &seen.near[0][0];
    let names = |v: &[Nb]| -> Vec<String> {
        v.iter()
            .map(|n| n.raw.rsplit('\t').next().unwrap().to_string())
            .collect()
    };
    assert_eq!(names(left), ["l1", "l2", "l3", "l4", "l5"]);
    assert_eq!(dists(a_iv, left), [11; 5]);
    assert_eq!(names(right), ["r1", "r2", "r3", "r4", "r5"]);
    assert_eq!(dists(a_iv, right), [11; 5]);

    // k larger than the number of B: every B, nearest first.
    let b = "chr1\t10\t20\nchr1\t30\t99\nchr1\t50\t60\nchr1\t111\t112\nchr1\t300\t301\n";
    let (seen, _) = run(a, &[b], &nearest(10, None, StrandRule::Any)).unwrap();
    let (left, right) = &seen.near[0][0];
    assert_eq!(raws(left), ["chr1\t30\t99", "chr1\t50\t60", "chr1\t10\t20"]);
    assert_eq!(dists(a_iv, left), [2, 41, 81]);
    assert_eq!(raws(right), ["chr1\t111\t112", "chr1\t300\t301"]);
    assert_eq!(dists(a_iv, right), [2, 191]);

    // A point on either edge of A is at distance 1, unless
    // --inclusive-points makes it an overlap (a hit, not a neighbour).
    let b = "chr1\t100\t100\tp_start\nchr1\t110\t110\tp_end\n";
    let (seen, _) = run(a, &[b], &k1).unwrap();
    let (left, right) = &seen.near[0][0];
    assert_eq!(raws(left), ["chr1\t100\t100\tp_start"]);
    assert_eq!(raws(right), ["chr1\t110\t110\tp_end"]);
    assert_eq!(dists(a_iv, left), [1]);
    assert_eq!(dists(a_iv, right), [1]);
    assert!(seen.hits[0].is_empty());
    let inclusive = SweepConfig {
        inclusive_points: true,
        ..k1
    };
    let (seen, _) = run(a, &[b], &inclusive).unwrap();
    assert_eq!(seen.near[0][0], (vec![], vec![]));
    assert_eq!(seen.hits[0].len(), 2);

    // A point A (a VCF insertion: a point after POS 101) and B bookended
    // on either side of it are at distance 1.
    let a_ins = format!("{VCF_HEADER}chr1\t101\t.\tA\tAT\t.\t.\t.\n");
    let b = "chr1\t95\t101\nchr1\t101\t105\n";
    let (seen, _) = run(&a_ins, &[b], &k1).unwrap();
    assert_eq!(seen.a_span[0], Interval::point(101));
    let (left, right) = &seen.near[0][0];
    assert_eq!(raws(left), ["chr1\t95\t101"]);
    assert_eq!(raws(right), ["chr1\t101\t105"]);
    assert_eq!(dists(Interval::point(101), left), [1]);
    assert_eq!(dists(Interval::point(101), right), [1]);

    // Without a nearest request (or with k = 0) there are no neighbours.
    for cfg in [SweepConfig::default(), nearest(0, None, StrandRule::Any)] {
        let (seen, _) = run(a, &["chr1\t10\t20\nchr1\t200\t210\n"], &cfg).unwrap();
        assert_eq!(seen.near[0][0], (vec![], vec![]));
    }
}

#[test]
fn max_distance_bounds_both_sides() {
    // B at distance 11 on each side.
    let a = "chr1\t100\t110\n";
    let b = "chr1\t80\t90\nchr1\t120\t130\n";
    for (max, found) in [(0, false), (10, false), (11, true), (HUGE, true)] {
        let (seen, _) = run(a, &[b], &nearest(1, Some(max), StrandRule::Any)).unwrap();
        let (left, right) = &seen.near[0][0];
        assert_eq!(
            (left.len(), right.len()),
            if found { (1, 1) } else { (0, 0) },
            "max distance {max}"
        );
    }
}

#[test]
fn strand_rules_keep_a_pool_per_strand() {
    // Left of A (+): `-` at distance 41, `+` at 61, `-` at 81. Right:
    // unknown at 1, `-` at 2, `+` at 11.
    let a = "chr1\t100\t110\ta\t0\t+\n";
    let b = "chr1\t10\t20\tm1\t0\t-\nchr1\t30\t40\tp1\t0\t+\nchr1\t50\t60\tm2\t0\t-\n\
             chr1\t110\t111\tu\t0\t.\nchr1\t111\t112\tm3\t0\t-\nchr1\t120\t121\tp2\t0\t+\n";
    let name = |v: &[Nb]| -> Vec<String> {
        v.iter()
            .map(|n| n.raw.split('\t').nth(3).unwrap().to_string())
            .collect()
    };
    for (rule, left, right) in [
        (StrandRule::Any, ["m2"], ["u"]),
        (StrandRule::Same, ["p1"], ["p2"]),
        (StrandRule::Opposite, ["m2"], ["m3"]),
    ] {
        let (seen, _) = run(a, &[b], &nearest(1, None, rule)).unwrap();
        let (l, r) = &seen.near[0][0];
        assert_eq!(
            (name(l), name(r)),
            (
                left.map(String::from).to_vec(),
                right.map(String::from).to_vec()
            ),
            "{rule:?}"
        );
    }
    // A with unknown strand (BED without column 6, or any VCF record) has
    // no neighbours under a strand rule.
    for a in [
        "chr1\t100\t110\n".to_string(),
        format!("{VCF_HEADER}chr1\t101\t.\tA\tG\t.\t.\t.\n"),
    ] {
        for rule in [StrandRule::Same, StrandRule::Opposite] {
            let (seen, _) = run(&a, &[b], &nearest(1, None, rule)).unwrap();
            assert_eq!(seen.near[0][0], (vec![], vec![]), "{rule:?}");
        }
    }
}

#[test]
fn split_by_b_strand_keeps_k_per_group() {
    // A (+). Left: `-` m1 at 81, `+` p1 at 41, `.` u1 at 31. Right: `+`
    // p2, p3, p4 at 2, 6 and 9, then `-` m2 at 21 past them. With the
    // split, `-` B are a group of their own and `.` B join `+` (D6): each
    // side holds the k nearest of each group, by distance. One pool would
    // drop m1 behind p1 and u1; one look-ahead count would stop before m2.
    let a = "chr1\t100\t110\ta\t0\t+\n";
    let b = "chr1\t10\t20\tm1\t0\t-\nchr1\t50\t60\tp1\t0\t+\nchr1\t60\t70\tu1\t0\t.\n\
             chr1\t111\t112\tp2\t0\t+\nchr1\t115\t116\tp3\t0\t+\nchr1\t118\t119\tp4\t0\t+\n\
             chr1\t130\t131\tm2\t0\t-\n";
    let name = |v: &[Nb]| -> Vec<String> {
        v.iter()
            .map(|n| n.raw.split('\t').nth(3).unwrap().to_string())
            .collect()
    };
    let names = |v: &[&str]| -> Vec<String> { v.iter().map(|s| s.to_string()).collect() };
    let cases: [(SweepConfig, &[&str], &[&str]); 5] = [
        (nearest(1, None, StrandRule::Any), &["u1"], &["p2"]),
        (split(1, StrandRule::Any), &["u1", "m1"], &["p2", "m2"]),
        (
            split(2, StrandRule::Any),
            &["u1", "p1", "m1"],
            &["p2", "p3", "m2"],
        ),
        // Under a strand rule the flag changes nothing.
        (split(1, StrandRule::Same), &["p1"], &["p2"]),
        (split(1, StrandRule::Opposite), &["m1"], &["m2"]),
    ];
    for (cfg, left, right) in cases {
        let (seen, _) = run(a, &[b], &cfg).unwrap();
        let (l, r) = &seen.near[0][0];
        assert_eq!((name(l), name(r)), (names(left), names(right)), "{cfg:?}");
    }
}

#[test]
fn left_pool_goes_by_a_key() {
    // A1 `chr1 101 . AAATTT AAA` has key 100 and span [103,106); A2
    // `chr1 102 . A G` has key 101 and span [101,102), before A1's span.
    // For A1, B2 [99,102) is the nearest left B (distance 2) and B1 is
    // farther (9). For A2, B2 overlaps, so B1 is the nearest left B. B2
    // ends after A1's key: treating it as gone for good at A1 (it ends
    // before A1's span) would drop B1 behind it. With a maximum distance
    // of 7, B1 is just in reach of A2 (7) though out of reach of A1's span
    // (9): dropping it by A1's span start instead of its key loses it.
    let a =
        format!("{VCF_HEADER}chr1\t101\t.\tAAATTT\tAAA\t.\t.\t.\nchr1\t102\t.\tA\tG\t.\t.\t.\n");
    let b = "chr1\t90\t95\tB1\nchr1\t99\t102\tB2\n";
    for max in [None, Some(7)] {
        let (seen, _) = run(&a, &[b], &nearest(1, max, StrandRule::Any)).unwrap();
        assert_eq!(raws(&seen.near[0][0].0), ["chr1\t99\t102\tB2"], "{max:?}");
        assert_eq!(seen.hits[1], [(0, "chr1\t99\t102\tB2".to_string())]);
        assert_eq!(raws(&seen.near[1][0].0), ["chr1\t90\t95\tB1"], "{max:?}");
    }
    let (seen, _) = run(&a, &[b], &nearest(1, None, StrandRule::Any)).unwrap();
    assert_eq!(
        seen.a_span,
        [
            Interval {
                start: 103,
                end: 106
            },
            Interval {
                start: 101,
                end: 102
            }
        ]
    );
}

#[test]
fn look_ahead_goes_by_b_key() {
    // B is VCF, sorted by key, but its spans start after the key: B1 (key
    // 110) starts at 130 and B2 (key 111) at 141, yet B3 (key 112) starts
    // at 112 and is the nearest right of A [100,105). All three lie past
    // the load limit (105). Stopping once a start passes the k-th start
    // (141 > 130) would miss B3; only a key beyond the k-th start can end
    // the look-ahead, since every later B starts at or after its key.
    let line = |key, shape, id: &str| vcf_rec("chr1", key, shape, id.to_string()).line();
    let b = format!(
        "{VCF_HEADER}{}\n{}\n{}\n{}\n",
        line(110, VcfShape::Del { prefix: 20, len: 1 }, "B1"),
        line(111, VcfShape::Del { prefix: 30, len: 1 }, "B2"),
        line(112, VcfShape::Snv, "B3"),
        line(200, VcfShape::Snv, "B4"),
    );
    let (seen, _) = run(
        "chr1\t100\t105\n",
        &[&b],
        &nearest(1, None, StrandRule::Any),
    )
    .unwrap();
    let right = &seen.near[0][0].1;
    assert_eq!(right.len(), 1);
    assert!(right[0].raw.contains("\tB3\t"), "{right:?}");
    assert_eq!(
        right[0].span,
        Interval {
            start: 112,
            end: 113
        }
    );
}

#[test]
fn look_ahead_stops_at_the_chromosome_end_and_keeps_order_checks() {
    // Reading ahead for the chr1 A meets B's chr2 record: chr1's right
    // side is empty, and the chr2 A still finds that record.
    let a = "chr1\t100\t110\nchr2\t5\t6\n";
    let b = "chr1\t10\t20\nchr2\t50\t60\n";
    let (seen, _) = run(a, &[b], &nearest(1, None, StrandRule::Any)).unwrap();
    assert_eq!(seen.near[0][0].0.len(), 1);
    assert!(seen.near[0][0].1.is_empty());
    assert!(seen.near[1][0].0.is_empty());
    assert_eq!(raws(&seen.near[1][0].1), ["chr2\t50\t60"]);

    // A B record read only by looking ahead is order-checked like any
    // other: the unsorted B is an error naming it.
    let b = "chr1\t120\t130\nchr1\t115\t116\n";
    let cfg = nearest(1, None, StrandRule::Any);
    let err = run("chr1\t100\t110\n", &[b], &cfg).unwrap_err();
    let IoError::Order(e) = err else {
        panic!("expected an order error, got {err}");
    };
    assert_eq!(e.path, "b0");
}

#[test]
fn retention_stays_near_k() {
    // 100,000 B and 10 sparse A points, k = 2: the left pool keeps the two
    // latest ends, and reading ahead stops after a few B.
    let b: String = (0..100_000u64)
        .map(|i| format!("chr1\t{}\t{}\n", i * 10, i * 10 + 5))
        .collect();
    let a: String = (1..=10u64)
        .map(|j| format!("chr1\t{p}\t{p}\n", p = j * 99_997))
        .collect();
    let cfg = nearest(2, None, StrandRule::Any);
    let (seen, summary) = run(&a, &[&b], &cfg).unwrap();
    assert!(
        summary.stats.peak_pool <= 2,
        "peak pool {}",
        summary.stats.peak_pool
    );
    assert!(
        summary.stats.peak_active <= 4,
        "peak active {}",
        summary.stats.peak_active
    );
    for (i, files) in seen.near.iter().enumerate() {
        let (left, right) = &files[0];
        assert_eq!((left.len(), right.len()), (2, 2), "A {i}");
    }

    // 1,000 B ending at the same base are all ties: every one is kept.
    let b: String = (0..1000u64).map(|i| format!("chr1\t{i}\t1000\n")).collect();
    let (seen, summary) = run(
        "chr1\t5000\t5001\n",
        &[&b],
        &nearest(1, None, StrandRule::Any),
    )
    .unwrap();
    assert_eq!(seen.near[0][0].0.len(), 1000);
    assert_eq!(summary.stats.peak_pool, 1000);

    // With k = 1,000 and a maximum distance of 50, the pool keeps only
    // the B in reach: those ending within 50 of the current A's key.
    let b: String = (0..100_000u64)
        .map(|i| format!("chr1\t{}\t{}\n", i * 10, i * 10 + 5))
        .collect();
    // Nine A points well inside B's range, each with 4 to 6 B in reach on
    // either side.
    let a: String = (1..=9u64)
        .map(|j| format!("chr1\t{p}\t{p}\n", p = j * 99_997))
        .collect();
    let (seen, summary) = run(&a, &[&b], &nearest(1000, Some(50), StrandRule::Any)).unwrap();
    assert!(
        summary.stats.peak_pool <= 6,
        "peak pool {}",
        summary.stats.peak_pool
    );
    for (i, files) in seen.near.iter().enumerate() {
        let (left, right) = &files[0];
        assert!(
            (4..=6).contains(&left.len()) && (4..=6).contains(&right.len()),
            "A {i}"
        );
    }

    // With `split_by_b_strand`, reading ahead waits for k B of each group:
    // B alternate `+`/`-` here, so it still stops after a few. Under a
    // strand rule the split changes nothing and there is no second group
    // to wait for: B all `+` stop the look-ahead as soon as without it.
    let alternating: String = (0..100_000u64)
        .map(|i| {
            let strand = if i % 2 == 0 { '+' } else { '-' };
            format!("chr1\t{}\t{}\tb\t0\t{strand}\n", i * 10, i * 10 + 5)
        })
        .collect();
    let plus: String = (0..100_000u64)
        .map(|i| format!("chr1\t{}\t{}\tb\t0\t+\n", i * 10, i * 10 + 5))
        .collect();
    let a: String = (1..=9u64)
        .map(|j| format!("chr1\t{p}\t{p}\ta\t0\t+\n", p = j * 99_997))
        .collect();
    for (b, cfg, sides) in [
        (&alternating, split(1, StrandRule::Any), 2),
        (&plus, split(1, StrandRule::Same), 1),
    ] {
        let (seen, summary) = run(&a, &[b], &cfg).unwrap();
        assert!(
            summary.stats.peak_active <= 6,
            "{cfg:?}: peak active {}",
            summary.stats.peak_active
        );
        assert!(
            summary.stats.peak_pool <= 2,
            "{cfg:?}: peak pool {}",
            summary.stats.peak_pool
        );
        for (i, files) in seen.near.iter().enumerate() {
            let (left, right) = &files[0];
            assert_eq!((left.len(), right.len()), (sides, sides), "{cfg:?}, A {i}");
        }
    }
}

#[test]
fn pool_bytes_stay_bounded_across_a_gap() {
    // Two A 10 Mbp apart with 100,000 B between them: at the second A the
    // whole gap is skipped into the left pool in one step. The pool must
    // reclaim dropped members' bytes as it goes rather than buffer the
    // gap (spec §8.2: memory bounded by depth and k, not file size).
    let b: String = (0..100_000u64)
        .map(|i| {
            let s = 1_000 + i * 99;
            format!("chr1\t{s}\t{}\tgap_feature_{i:06}\t0\t+\n", s + 50)
        })
        .collect();
    let longest = b.lines().map(str::len).max().unwrap();
    let last = b.lines().last().unwrap();
    let a = "chr1\t0\t10\ta1\t0\t+\nchr1\t10000000\t10000010\ta2\t0\t+\n";
    for rule in [StrandRule::Any, StrandRule::Same] {
        let (seen, summary) = run(a, &[&b], &nearest(1, None, rule)).unwrap();
        assert_eq!(raws(&seen.near[1][0].0), [last], "{rule:?}");
        // One member (k = 1, no ties): the arena holds it, the bytes
        // dropped since the last compaction (at most as many) and the
        // record being inserted.
        assert!(
            summary.stats.peak_pool_bytes <= 3 * longest,
            "{rule:?}: peak pool bytes {} for {longest}-byte lines",
            summary.stats.peak_pool_bytes
        );
    }
}
