//! Unsorted mode (D39, D40): every B is loaded into per-chromosome indexes
//! and A is streamed in any order. The visitor must see, for each A, what
//! the sorted sweep shows it: the same hits in the same order (D38) and the
//! same neighbours. Checked against `sweep_sorted` over shuffled A, and
//! against the brute-force references over shuffled A and B, for BED inputs
//! and for BED/VCF inputs whose VCF spans start after their key.

use gloci_core::{ChromTable, Interval, Strand, StrandRule, VcfSpan, Window};
use gloci_engine::{
    ARec, Hit, ITree64, Inputs, NearestConfig, Neighbours, SweepConfig, SweepSummary, Visitor,
    sweep_sorted, sweep_unsorted,
};
use gloci_io::{InvalidPolicy, IoError, ReaderOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, Sides, naive_hits, naive_neighbours};
use gloci_testkit::strategies::{
    CHROMS, Feature, HUGE, MixedInputs, RelRec, TestFile, relation_inputs, vcf_relation_inputs,
};
use proptest::prelude::*;
use proptest::test_runner::TestCaseError;
use std::num::NonZeroUsize;

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn opts() -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn text(b: &[u8]) -> String {
    String::from_utf8(b.to_vec()).unwrap()
}

/// Every field of a hit or neighbour: file, span, strand, raw line, line
/// number.
type FullHit = (u16, Interval, Strand, String, u64);

fn full(h: &Hit<'_>) -> FullHit {
    (h.file, h.span, h.strand, text(h.raw), h.loc)
}

/// One A as the visitor saw it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SeenA {
    line: String,
    loc: u64,
    hits: Vec<FullHit>,
    /// Per B file, its `(left, right)` neighbours.
    near: Vec<(Vec<FullHit>, Vec<FullHit>)>,
}

struct Record {
    files: u16,
    seen: Vec<SeenA>,
}

impl Visitor for Record {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        let side = |hs: &[Hit<'_>]| hs.iter().map(full).collect::<Vec<_>>();
        self.seen.push(SeenA {
            line: text(a.raw),
            loc: a.meta.loc,
            hits: hits.iter().map(full).collect(),
            near: (0..self.files)
                .map(|f| (side(near.left(f)), side(near.right(f))))
                .collect(),
        });
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Sorted,
    Unsorted,
}

/// Sweep in-memory A text against B texts (named `a`, `b0`, `b1`, …) with
/// a karyotype table.
fn run(
    mode: Mode,
    a: &str,
    b: &[&str],
    cfg: &SweepConfig,
) -> Result<(Vec<SeenA>, SweepSummary), IoError> {
    let o = opts();
    let inputs = Inputs {
        a: mem_reader("a", a, &o),
        b: b.iter()
            .enumerate()
            .map(|(i, t)| mem_reader(&format!("b{i}"), t, &o))
            .collect(),
    };
    let mut chroms = ChromTable::karyotype();
    let mut rec = Record {
        files: b.len() as u16,
        seen: Vec::new(),
    };
    let summary = match mode {
        Mode::Sorted => sweep_sorted(inputs, &mut chroms, cfg, &mut rec)?,
        Mode::Unsorted => sweep_unsorted(inputs, &mut chroms, cfg, &mut rec)?,
    };
    if cfg.nearest.is_none() {
        let none = |s: &SeenA| s.near.iter().all(|(l, r)| l.is_empty() && r.is_empty());
        assert!(
            rec.seen.iter().all(none),
            "neighbours without a nearest request"
        );
    }
    Ok((rec.seen, summary))
}

/// Each A's hits as `(file, B line)`, as `naive_hits` gives them.
fn hit_lines(seen: &[SeenA]) -> Vec<Vec<(u16, String)>> {
    seen.iter()
        .map(|s| s.hits.iter().map(|h| (h.0, h.3.clone())).collect())
        .collect()
}

/// Each A's neighbours as indices into each B file's records, as
/// `naive_neighbours` gives them. Record `i` of file `f` is on line
/// `heads[f] + i + 1`, and its line must be the neighbour's raw bytes.
fn neighbour_indices(seen: &[SeenA], heads: &[u64], b: &[Vec<impl Feature>]) -> Vec<Vec<Sides>> {
    seen.iter()
        .map(|s| {
            s.near
                .iter()
                .enumerate()
                .map(|(f, (left, right))| {
                    let index = |v: &[FullHit]| -> Vec<usize> {
                        v.iter()
                            .map(|h| {
                                assert_eq!(h.0, f as u16, "a neighbour of file {f} tagged {}", h.0);
                                let i = (h.4 - 1 - heads[f]) as usize;
                                assert_eq!(h.3, b[f][i].line());
                                assert_eq!(h.1, b[f][i].span());
                                assert_eq!(h.2, b[f][i].strand());
                                i
                            })
                            .collect()
                    };
                    (index(left), index(right))
                })
                .collect()
        })
        .collect()
}

/// Check one unsorted run over `a` and `b` against the brute-force
/// references, for `cfg`.
fn check_naive<R: Feature>(
    seen: &[SeenA],
    a: &[R],
    b: &[Vec<R>],
    heads: &[u64],
    cfg: &SweepConfig,
) {
    assert_eq!(
        hit_lines(seen),
        naive_hits(a, b, cfg.window, cfg.inclusive_points),
        "hits, {cfg:?}"
    );
    if let Some(n) = cfg.nearest {
        assert_eq!(
            neighbour_indices(seen, heads, b),
            naive_neighbours(
                a,
                b,
                n.k,
                n.max_distance,
                n.strand,
                n.split_by_b_strand,
                cfg.inclusive_points
            ),
            "neighbours, {cfg:?}"
        );
    }
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

fn max_distance() -> impl Strategy<Value = Option<u64>> {
    prop::option::of(prop_oneof![
        2 => 0u64..4,
        2 => 4u64..100,
        1 => 100u64..100_000,
        1 => Just(HUGE),
    ])
}

/// `input` with a random order for A's records and for each B's (as
/// positions in the sorted files).
fn shuffled(
    input: MixedInputs,
) -> impl Strategy<Value = (MixedInputs, Vec<usize>, Vec<Vec<usize>>)> {
    let order = |n: usize| Just((0..n).collect::<Vec<usize>>()).prop_shuffle();
    let b: Vec<_> = input.b.iter().map(|f| order(f.recs.len())).collect();
    (order(input.a.recs.len()), b).prop_map(move |(a, b)| (input.clone(), a, b))
}

/// `f`'s records in the order `order` (positions in `f`).
fn permuted(f: &TestFile, order: &[usize]) -> TestFile {
    TestFile {
        format: f.format,
        recs: order.iter().map(|&i| f.recs[i].clone()).collect(),
    }
}

/// No nearest request, then every strand rule with `split_by_b_strand`
/// off and on, each with `--inclusive-points` off and on.
fn configs(window: Option<Window>, k: usize, max_distance: Option<u64>) -> Vec<SweepConfig> {
    let mut out = Vec::new();
    for inclusive_points in [false, true] {
        out.push(SweepConfig {
            window,
            inclusive_points,
            nearest: None,
        });
        for strand in [StrandRule::Any, StrandRule::Same, StrandRule::Opposite] {
            for split_by_b_strand in [false, true] {
                out.push(SweepConfig {
                    window,
                    inclusive_points,
                    nearest: Some(NearestConfig {
                        k,
                        max_distance,
                        strand,
                        split_by_b_strand,
                    }),
                });
            }
        }
    }
    out
}

/// Sort `input`, then for every config of [`configs`]:
/// 1. A in sorted order through the unsorted path (sorted querents
///    throughout) gives exactly what `sweep_sorted` gives.
/// 2. A in the order `a_order` is visited in that order, and each A gets
///    `sweep_sorted`'s hits and neighbours for that A (same line).
/// 3. A and every B shuffled give the brute-force answer over B in its
///    shuffled order: B needs no order.
fn check_unsorted(
    mut input: MixedInputs,
    a_order: Vec<usize>,
    b_orders: Vec<Vec<usize>>,
    k: usize,
    max_distance: Option<u64>,
    window: Option<Window>,
) -> Result<(), TestCaseError> {
    let rank = naive::karyotype_rank(CHROMS);
    input.a.sort(&rank);
    for f in &mut input.b {
        f.sort(&rank);
    }
    let shuffled_a = permuted(&input.a, &a_order);
    let shuffled_b: Vec<TestFile> = input
        .b
        .iter()
        .zip(&b_orders)
        .map(|(f, o)| permuted(f, o))
        .collect();
    let a_text = input.a.text();
    let shuffled_a_text = shuffled_a.text();
    let b_texts: Vec<String> = input.b.iter().map(TestFile::text).collect();
    let b_refs: Vec<&str> = b_texts.iter().map(String::as_str).collect();
    let shuffled_b_texts: Vec<String> = shuffled_b.iter().map(TestFile::text).collect();
    let shuffled_b_refs: Vec<&str> = shuffled_b_texts.iter().map(String::as_str).collect();
    let shuffled_b_recs: Vec<Vec<RelRec>> = shuffled_b.iter().map(|f| f.recs.clone()).collect();
    let heads: Vec<u64> = input.b.iter().map(TestFile::header_lines).collect();
    let b_lens: Vec<u64> = input.b.iter().map(|f| f.recs.len() as u64).collect();
    let a_head = input.a.header_lines();
    for cfg in configs(window, k, max_distance) {
        let (sorted, s_sum) = run(Mode::Sorted, &a_text, &b_refs, &cfg).unwrap();

        let (same, u_sum) = run(Mode::Unsorted, &a_text, &b_refs, &cfg).unwrap();
        prop_assert_eq!(&same, &sorted, "sorted A, {:?}", cfg);
        prop_assert_eq!(u_sum.a_records, s_sum.a_records);
        prop_assert_eq!(u_sum.shared_chroms, s_sum.shared_chroms);
        prop_assert_eq!(&u_sum.b_records, &b_lens, "every B record is loaded");

        let (shuf, _) = run(Mode::Unsorted, &shuffled_a_text, &b_refs, &cfg).unwrap();
        prop_assert_eq!(shuf.len(), sorted.len());
        for (j, (got, &i)) in shuf.iter().zip(&a_order).enumerate() {
            let want = &sorted[i];
            prop_assert_eq!(&got.line, &want.line);
            prop_assert_eq!(got.loc, a_head + j as u64 + 1);
            prop_assert_eq!(&got.hits, &want.hits, "hits of {}, {:?}", got.line, cfg);
            prop_assert_eq!(
                &got.near,
                &want.near,
                "neighbours of {}, {:?}",
                got.line,
                cfg
            );
        }

        let (both, _) = run(Mode::Unsorted, &shuffled_a_text, &shuffled_b_refs, &cfg).unwrap();
        check_naive(&both, &shuffled_a.recs, &shuffled_b_recs, &heads, &cfg);
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// BED inputs, where every key is the span start.
    #[test]
    fn unsorted_equals_sorted(
        (input, a_order, b_orders) in relation_inputs().prop_map(MixedInputs::from).prop_flat_map(shuffled),
        k in 1usize..=3,
        max_d in max_distance(),
        window in window(),
    ) {
        check_unsorted(input, a_order, b_orders, k, max_d, window)?;
    }

    /// BED and VCF inputs, where VCF spans (of A and of B) can start after
    /// the key: sorted A's queries are not in span order.
    #[test]
    fn unsorted_equals_sorted_vcf_style(
        (input, a_order, b_orders) in vcf_relation_inputs().prop_flat_map(shuffled),
        k in 1usize..=3,
        max_d in max_distance(),
        window in window(),
    ) {
        check_unsorted(input, a_order, b_orders, k, max_d, window)?;
    }
}

fn coord() -> impl Strategy<Value = u64> {
    prop_oneof![
        4 => 0u64..64,
        2 => 0u64..100_000,
        1 => (1u64 << 31) - 8..(1u64 << 31) + 8,
        1 => u64::MAX - 64..=u64::MAX,
        1 => any::<u64>(),
    ]
}

fn item() -> impl Strategy<Value = Interval> {
    let len = prop_oneof![
        2 => Just(0u64),
        4 => 0u64..16,
        2 => 0u64..100_000,
        1 => any::<u64>(),
    ];
    (coord(), len).prop_map(|(start, l)| Interval {
        start,
        end: start.saturating_add(l),
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// `ITree64::query(lo, hi)` returns exactly the items whose closed
    /// range `[start, end]` meets `[lo, hi]`, each once.
    #[test]
    fn itree64_matches_brute_force(
        items in prop::collection::vec(item(), 0..300),
        queries in prop::collection::vec((coord(), 0u64..1000), 1..20),
    ) {
        let tree = ITree64::new(items.iter().enumerate().map(|(i, &iv)| (iv, i as u32)).collect());
        for (lo, len) in queries {
            let hi = lo.saturating_add(len);
            let mut got = vec![u32::MAX];
            tree.query(lo, hi, &mut got);
            prop_assert_eq!(got.remove(0), u32::MAX, "query appends to out");
            got.sort_unstable();
            let want: Vec<u32> = items
                .iter()
                .enumerate()
                .filter(|(_, iv)| iv.start <= hi && iv.end >= lo)
                .map(|(i, _)| i as u32)
                .collect();
            prop_assert_eq!(got, want, "[{}, {}]", lo, hi);
        }
    }
}

fn rec(chrom: &str, start: u64, end: u64, strand: char, name: &str) -> TestRec {
    TestRec {
        chrom: chrom.to_string(),
        start,
        end,
        strand,
        name: name.to_string(),
    }
}

fn bed(recs: &[TestRec]) -> String {
    recs.iter().map(|r| r.bed_line() + "\n").collect()
}

fn nearest(k: usize, max_distance: Option<u64>, strand: StrandRule, split: bool) -> SweepConfig {
    SweepConfig {
        window: None,
        inclusive_points: false,
        nearest: Some(NearestConfig {
            k,
            max_distance,
            strand,
            split_by_b_strand: split,
        }),
    }
}

/// Plain hits, `--inclusive-points`, a strand-aware window, and nearest
/// requests (several `k`, a strand rule, the split, a `max_distance`).
fn fixed_configs() -> Vec<SweepConfig> {
    let plain = SweepConfig::default();
    vec![
        plain,
        SweepConfig {
            inclusive_points: true,
            ..plain
        },
        SweepConfig {
            window: Some(Window {
                left: 5,
                right: 20,
                strand_aware: true,
            }),
            ..plain
        },
        nearest(1, None, StrandRule::Any, false),
        nearest(2, None, StrandRule::Any, true),
        nearest(1, None, StrandRule::Same, false),
        nearest(3, Some(60), StrandRule::Opposite, false),
    ]
}

/// The unsorted sweep of `a` against the single B file `b`, for every
/// fixed config, equals the brute-force answer; returns the summaries.
fn check_fixed(a: &[TestRec], b: &[TestRec]) -> Vec<SweepSummary> {
    let (a_text, b_text) = (bed(a), bed(b));
    fixed_configs()
        .iter()
        .map(|cfg| {
            let (seen, summary) = run(Mode::Unsorted, &a_text, &[&b_text], cfg).unwrap();
            check_naive(&seen, a, &[b.to_vec()], &[0], cfg);
            summary
        })
        .collect()
}

const BIG: u64 = 3_000_000_000;

#[test]
fn itree64_fallback() {
    // chr1 is long: some of its B lie past 2^31, so its index is an
    // ITree64. chr2's B fit coitrees' i32. Both files are unsorted.
    let b = [
        rec("chr1", BIG, BIG + 100, '-', "far"),
        rec("chr2", 100, 200, '+', "x"),
        rec("chr1", 10, 20, '+', "small"),
        rec("chr1", BIG - 10, BIG, '+', "book"),
        rec("chr1", BIG + 50, BIG + 50, '.', "pt"),
        rec("chr2", 300, 300, '-', "p2"),
        rec("chr1", BIG + 100, BIG + 300, '+', "after"),
        rec("chr1", BIG + 400, BIG + 410, '-', "later"),
    ];
    let a = [
        // Bookended with far, overlaps after.
        rec("chr1", BIG + 100, BIG + 110, '+', "a1"),
        rec("chr2", 150, 160, '-', "a2"),
        // A point on far's start and book's end.
        rec("chr1", BIG, BIG, '+', "a3"),
        // Equal to pt, inside far.
        rec("chr1", BIG + 50, BIG + 50, '-', "a4"),
        // Past i32 on chr2, whose index is a coitree: no hits, and every
        // chr2 B is a left neighbour.
        rec("chr2", BIG, BIG + 1, '+', "a5"),
        rec("chr1", 0, 15, '+', "a6"),
        // Bookended on book's left.
        rec("chr1", BIG - 20, BIG - 10, '-', "a7"),
        rec("chr1", BIG + 320, BIG + 330, '.', "a8"),
    ];
    for summary in check_fixed(&a, &b) {
        assert_eq!(summary.stats.itree64_indexes, 1, "chr1 is an ITree64");
        assert_eq!(summary.stats.coitrees_indexes, 1, "chr2 is a coitree");
        assert_eq!(summary.b_records, [8]);
        assert!(summary.shared_chroms);
    }
    // The fallback finds real hits: a1 overlaps after, a4 is inside far and
    // equal to pt.
    let (seen, _) = run(
        Mode::Unsorted,
        &bed(&a),
        &[&bed(&b)],
        &SweepConfig::default(),
    )
    .unwrap();
    let names = |s: &SeenA| -> Vec<String> {
        s.hits
            .iter()
            .map(|h| h.3.split('\t').nth(3).unwrap().to_string())
            .collect()
    };
    assert_eq!(names(&seen[0]), ["after"]);
    assert_eq!(names(&seen[3]), ["far", "pt"]);
}

#[test]
fn index_choice_at_the_i32_limit() {
    // coitrees is used while every padded B end is below i32::MAX - 2.
    let lim = (i32::MAX - 2) as u64;
    let b = [
        rec("chr1", lim - 10, lim - 1, '+', "c1"),
        rec("chr1", lim - 1, lim - 1, '-', "c2"),
        rec("chr1", 5, 9, '+', "c3"),
        rec("chr2", lim - 10, lim, '+', "w1"),
        rec("chr2", lim, lim, '-', "w2"),
        rec("chr2", 5, 9, '+', "w3"),
    ];
    let a: Vec<TestRec> = ["chr1", "chr2"]
        .iter()
        .flat_map(|c| {
            [
                rec(c, lim - 1, lim - 1, '+', "p1"),
                rec(c, lim, lim, '-', "p2"),
                rec(c, lim - 5, lim + 5, '+', "span"),
                rec(c, lim + 1, lim + 2, '-', "past"),
                rec(c, u64::MAX - 1, u64::MAX, '+', "end"),
                rec(c, 0, 7, '+', "low"),
            ]
        })
        .collect();
    for (cfg, summary) in fixed_configs().iter().zip(check_fixed(&a, &b)) {
        // chr1 is a coitree and chr2 an ITree64, unless the window (20
        // bases) pads chr1's B past the limit too.
        let want = if cfg.window.is_some() { (0, 2) } else { (1, 1) };
        let got = (
            summary.stats.coitrees_indexes,
            summary.stats.itree64_indexes,
        );
        assert_eq!(got, want, "{cfg:?}");
    }
    let cfg = SweepConfig {
        window: Some(Window {
            left: 0,
            right: 1,
            strand_aware: false,
        }),
        ..SweepConfig::default()
    };
    let (seen, summary) = run(Mode::Unsorted, &bed(&a), &[&bed(&b)], &cfg).unwrap();
    assert_eq!(
        summary.stats.itree64_indexes, 2,
        "a 1-base pad reaches the limit"
    );
    check_naive(&seen, &a, &[b.to_vec()], &[0], &cfg);
}

#[test]
fn sorted_querent_switch() {
    // Overlapping, nested and bookended B on chr1 and chr2, unsorted.
    let b: Vec<TestRec> = (0..60u64)
        .rev()
        .map(|i| {
            let chrom = if i % 3 == 0 { "chr2" } else { "chr1" };
            let start = i * 7 % 200;
            let strand = ['+', '-', '.'][(i % 3) as usize];
            rec(chrom, start, start + 5 + i % 40, strand, &format!("b{i}"))
        })
        .collect();
    let at = |chrom: &str, start: u64, name: &str| rec(chrom, start, start + 6, '+', name);
    // Sorted on chr1, then on chr2; then back on chr1 (out of order), then
    // chr2 again in key order: plain queries from the switch on.
    let mut a: Vec<TestRec> = (0..10u64).map(|i| at("chr1", i * 20, "s1")).collect();
    a.extend((0..6u64).map(|i| at("chr2", i * 30, "s2")));
    a.push(at("chr1", 50, "back"));
    a.extend((0..4u64).map(|i| at("chr2", i * 40, "late")));
    for summary in check_fixed(&a, &b) {
        // One B file, a coitree on each chromosome: one query per A.
        assert_eq!(summary.stats.querent_queries, 16);
        assert_eq!(summary.stats.coitrees_indexes, 2);
    }

    // A smaller key on the same chromosome switches at that record.
    let a = [
        at("chr1", 0, "x"),
        at("chr1", 10, "x"),
        at("chr1", 20, "x"),
        at("chr1", 15, "x"),
        at("chr1", 30, "x"),
    ];
    for summary in check_fixed(&a, &b) {
        assert_eq!(summary.stats.querent_queries, 3);
    }

    // The same A sorted keeps the querent to the end, with the same
    // results for every A (keyed by line) as the sorted sweep.
    let mut sorted = a.to_vec();
    sorted.sort_by_key(|r| r.start);
    let (a_text, b_text) = (
        bed(&sorted),
        bed(&naive::sort(b, naive::karyotype_rank(CHROMS))),
    );
    for cfg in fixed_configs() {
        let (want, _) = run(Mode::Sorted, &a_text, &[&b_text], &cfg).unwrap();
        let (got, summary) = run(Mode::Unsorted, &a_text, &[&b_text], &cfg).unwrap();
        assert_eq!(got, want, "{cfg:?}");
        assert_eq!(summary.stats.querent_queries, 5);
    }
}

#[test]
fn edges_and_empties() {
    let cfg = SweepConfig::default();
    let iv = |start, end| Interval { start, end };

    // B unsorted, its chromosomes interleaved: hits come in B's input
    // order (D38), not by position.
    let b =
        "chr1\t50\t60\tb1\nchr2\t0\t10\tc1\nchr1\t10\t100\tb2\nchr10\t0\t5\td1\nchr1\t55\t56\tb3\n";
    let (seen, summary) = run(Mode::Unsorted, "chr1\t52\t58\n", &[b], &cfg).unwrap();
    let lines: Vec<&str> = seen[0].hits.iter().map(|h| h.3.as_str()).collect();
    assert_eq!(
        lines,
        ["chr1\t50\t60\tb1", "chr1\t10\t100\tb2", "chr1\t55\t56\tb3"]
    );
    assert_eq!(
        seen[0].hits.iter().map(|h| (h.1, h.4)).collect::<Vec<_>>(),
        [(iv(50, 60), 1), (iv(10, 100), 3), (iv(55, 56), 5)]
    );
    assert_eq!(summary.b_records, [5], "every B record is read");
    assert!(summary.shared_chroms);

    // Several B files: hits by file, then input order.
    let (seen, _) = run(
        Mode::Unsorted,
        "chr1\t52\t58\n",
        &["chr1\t57\t70\tx\nchr1\t0\t53\ty\n", b],
        &cfg,
    )
    .unwrap();
    let tags: Vec<(u16, u64)> = seen[0].hits.iter().map(|h| (h.0, h.4)).collect();
    assert_eq!(tags, [(0, 1), (0, 2), (1, 1), (1, 3), (1, 5)]);

    // A on a chromosome B lacks, then on one it has: A keeps its order.
    let a = "chrX\t0\t100\nchr2\t5\t6\nchr1\t0\t1\n";
    let (seen, summary) = run(Mode::Unsorted, a, &[b], &cfg).unwrap();
    assert_eq!(
        seen.iter().map(|s| s.line.as_str()).collect::<Vec<_>>(),
        a.lines().collect::<Vec<_>>()
    );
    assert_eq!(
        seen.iter().map(|s| s.hits.len()).collect::<Vec<_>>(),
        [0, 1, 0]
    );
    assert_eq!(summary.a_records, 3);
    assert!(summary.shared_chroms);

    // No chromosome in common.
    let (seen, summary) = run(Mode::Unsorted, "chrX\t0\t10\n", &[b], &cfg).unwrap();
    assert!(seen[0].hits.is_empty());
    assert!(!summary.shared_chroms);

    // Empty A; empty B; header-only VCF A.
    let (seen, summary) = run(Mode::Unsorted, "", &[b], &cfg).unwrap();
    assert!(seen.is_empty());
    assert_eq!(
        (summary.a_records, summary.b_records.as_slice()),
        (0, &[5][..])
    );
    let (seen, summary) = run(
        Mode::Unsorted,
        a,
        &[""],
        &nearest(1, None, StrandRule::Any, false),
    )
    .unwrap();
    assert!(
        seen.iter()
            .all(|s| s.hits.is_empty() && s.near == [(vec![], vec![])])
    );
    assert_eq!(summary.b_records, [0]);
    assert!(!summary.shared_chroms);
    let (seen, summary) = run(Mode::Unsorted, VCF_HEAD, &[b], &cfg).unwrap();
    assert!(seen.is_empty());
    assert_eq!(summary.a_records, 0);

    // Points on B's edges hit only with --inclusive-points; inside always.
    let a = "chr1\t20\t20\nchr1\t10\t10\nchr1\t15\t15\n";
    let b = "chr1\t10\t20\n";
    let (seen, _) = run(Mode::Unsorted, a, &[b], &cfg).unwrap();
    assert_eq!(
        seen.iter().map(|s| s.hits.len()).collect::<Vec<_>>(),
        [0, 0, 1]
    );
    let inclusive = SweepConfig {
        inclusive_points: true,
        ..cfg
    };
    let (seen, _) = run(Mode::Unsorted, a, &[b], &inclusive).unwrap();
    assert_eq!(
        seen.iter().map(|s| s.hits.len()).collect::<Vec<_>>(),
        [1, 1, 1]
    );
}
