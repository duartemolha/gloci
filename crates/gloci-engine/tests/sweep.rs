//! The sorted sweep (spec §8.2, D36, D38): for each A, in A's input order,
//! every B whose window-padded span overlaps A under §4.3, by file and then
//! B input order. Checked against a brute-force reference, over BED inputs
//! and over BED/VCF inputs whose VCF spans start after their key.

use gloci_core::{ChromTable, Interval, Strand, StrandRule, VcfSpan, Window};
use gloci_engine::{
    ARec, Hit, Inputs, NearestConfig, Neighbours, SweepConfig, SweepSummary, Visitor, sweep_sorted,
};
use gloci_io::{InvalidPolicy, IoError, LEX_SORTED_NOTE, ReaderOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, naive_hits};
use gloci_testkit::strategies::{
    CHROMS, Feature, HUGE, MixedInputs, RelRec, relation_inputs, vcf_relation_inputs,
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

/// What the visitor saw: each A (raw line, key, span) and its hits as
/// `(file, B line)`, and whether any A had neighbours in any of the
/// `files` B inputs.
#[derive(Debug, Default)]
struct Seen {
    files: u16,
    a: Vec<String>,
    a_pos: Vec<(u64, Interval)>,
    hits: Vec<Vec<(u16, String)>>,
    near_seen: bool,
}

fn text(b: &[u8]) -> String {
    String::from_utf8(b.to_vec()).unwrap()
}

impl Visitor for Seen {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        for f in 0..self.files {
            self.near_seen |= !near.left(f).is_empty() || !near.right(f).is_empty();
        }
        self.a.push(text(a.raw));
        self.a_pos.push((a.meta.key, a.meta.span));
        self.hits
            .push(hits.iter().map(|h| (h.file, text(h.raw))).collect());
        Ok(())
    }
}

/// Run the sweep over in-memory inputs `(name, text)` with a karyotype table.
fn sweep(
    a: (&str, &str),
    b: &[(&str, &str)],
    cfg: &SweepConfig,
) -> Result<(Seen, SweepSummary), IoError> {
    let o = opts();
    let inputs = Inputs {
        a: mem_reader(a.0, a.1, &o),
        b: b.iter().map(|(n, t)| mem_reader(n, t, &o)).collect(),
    };
    let mut chroms = ChromTable::karyotype();
    let mut seen = Seen {
        files: b.len() as u16,
        ..Seen::default()
    };
    let summary = sweep_sorted(inputs, &mut chroms, cfg, &mut seen)?;
    if cfg.nearest.is_none() {
        assert!(!seen.near_seen, "neighbours without a nearest request");
    }
    Ok((seen, summary))
}

fn bed(recs: &[TestRec]) -> String {
    recs.iter().map(|r| r.bed_line() + "\n").collect()
}

fn rec(chrom: &str, start: u64, end: u64, name: String) -> TestRec {
    TestRec {
        chrom: chrom.to_string(),
        start,
        end,
        strand: '+',
        name,
    }
}

fn hit(file: u16, line: &str) -> (u16, String) {
    (file, line.to_string())
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

/// Sort `input`, then sweep it without a window, with a window and with a
/// strand-aware window, each with `--inclusive-points` off and on, and
/// compare every run with the naive answer.
fn check_hits(mut input: MixedInputs, left: u64, right: u64) -> Result<(), TestCaseError> {
    let rank = naive::karyotype_rank(CHROMS);
    input.a.sort(&rank);
    for f in &mut input.b {
        f.sort(&rank);
    }
    let a = &input.a.recs;
    let b: Vec<Vec<RelRec>> = input.b.iter().map(|f| f.recs.clone()).collect();
    let a_text = input.a.text();
    let a_lines: Vec<String> = a.iter().map(Feature::line).collect();
    let b_names: Vec<String> = (0..b.len()).map(|i| format!("b{i}")).collect();
    let b_texts: Vec<String> = input.b.iter().map(|f| f.text()).collect();
    let b_inputs: Vec<(&str, &str)> = b_names
        .iter()
        .zip(&b_texts)
        .map(|(n, t)| (n.as_str(), t.as_str()))
        .collect();
    let shared = a
        .iter()
        .any(|ra| b.iter().flatten().any(|rb| rb.chrom() == ra.chrom()));
    let windows = [
        None,
        Some(Window {
            left,
            right,
            strand_aware: false,
        }),
        Some(Window {
            left,
            right,
            strand_aware: true,
        }),
    ];
    for window in windows {
        for inclusive_points in [false, true] {
            let cfg = SweepConfig {
                window,
                inclusive_points,
                nearest: None,
            };
            let (seen, summary) = sweep(("a", &a_text), &b_inputs, &cfg).unwrap();
            prop_assert_eq!(&seen.a, &a_lines);
            prop_assert_eq!(
                &seen.hits,
                &naive_hits(a, &b, window, inclusive_points),
                "window {:?}, inclusive points {}",
                window,
                inclusive_points
            );
            prop_assert_eq!(summary.a_records, a.len() as u64);
            prop_assert_eq!(summary.b_records.len(), b.len());
            for (n, f) in summary.b_records.iter().zip(&b) {
                prop_assert!(*n <= f.len() as u64);
            }
            prop_assert_eq!(summary.shared_chroms, shared);
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// BED inputs, where every key is the span start.
    #[test]
    fn sweep_matches_naive(input in relation_inputs(), left in pad(), right in pad()) {
        check_hits(MixedInputs::from(input), left, right)?;
    }

    /// BED and VCF inputs, where VCF spans (of A and of B) can start after
    /// the key: loading and eviction must go by A's key, not its span.
    #[test]
    fn sweep_matches_naive_vcf_style(input in vcf_relation_inputs(), left in pad(), right in pad()) {
        check_hits(input, left, right)?;
    }
}

/// Sweep `a` against one B and compare with the naive answer; the peak
/// active set must be the overlap depth (plus one), not the file size.
fn check_deep(a: &[TestRec], b: &[TestRec], depth: usize) {
    let (seen, summary) = sweep(
        ("a.bed", &bed(a)),
        &[("b.bed", &bed(b))],
        &SweepConfig::default(),
    )
    .unwrap();
    assert_eq!(seen.a.len(), a.len());
    assert!(
        seen.hits == naive_hits(a, &[b.to_vec()], None, false),
        "hits differ from the naive answer"
    );
    assert!(
        summary.peak_active <= depth + 1,
        "peak active set {} for depth {depth} ({} B records)",
        summary.peak_active,
        b.len()
    );
}

#[test]
fn deep_overlap_whole_chromosome() {
    // One whole-chromosome B against 100,000 A (some of them points, the
    // first on B's left edge).
    let whole = [rec("chr1", 0, 1_000_000_000, "whole".into())];
    let a: Vec<TestRec> = (0..100_000u64)
        .map(|i| rec("chr1", i * 10_000, i * 10_000 + i % 7, format!("a{i}")))
        .collect();
    check_deep(&a, &whole, 1);

    // 1,000 B nested 1,000 deep against 1,000 A.
    let nest = |off: u64, c: u64| -> Vec<TestRec> {
        (0..1000u64)
            .map(|i| {
                rec(
                    "chr1",
                    off + i * 10,
                    off + 100_000 - i * 10,
                    format!("b{c}_{i}"),
                )
            })
            .collect()
    };
    let a: Vec<TestRec> = (0..1000u64)
        .map(|j| rec("chr1", j * 100 + 5, j * 100 + 6, format!("a{j}")))
        .collect();
    check_deep(&a, &nest(0, 0), 1000);

    // That nest ten times along the chromosome: 10,000 B, still 1,000 deep.
    let mut b = Vec::new();
    let mut a = Vec::new();
    for c in 0..10u64 {
        let off = c * 200_000;
        b.extend(nest(off, c));
        a.extend((0..50u64).map(|j| {
            rec(
                "chr1",
                off + j * 2000 + 5,
                off + j * 2000 + 6,
                format!("a{c}_{j}"),
            )
        }));
    }
    check_deep(&a, &b, 1000);
}

#[test]
fn eviction_uses_a_key_not_span() {
    // A1 is a trimmed deletion (key 100, span [101,105)); A2 is an SNV
    // (key 101, span [101,102)). Both overlap B [100,103).
    let a = format!("{VCF_HEAD}chr1\t101\t.\tATTTT\tA\t.\t.\t.\nchr1\t102\t.\tA\tG\t.\t.\t.\n");
    let (seen, _) = sweep(
        ("a.vcf", &a),
        &[("b.bed", "chr1\t100\t103\n")],
        &SweepConfig::default(),
    )
    .unwrap();
    assert_eq!(
        seen.a_pos,
        [
            (
                100,
                Interval {
                    start: 101,
                    end: 105
                }
            ),
            (
                101,
                Interval {
                    start: 101,
                    end: 102
                }
            )
        ]
    );
    assert_eq!(seen.hits, vec![vec![hit(0, "chr1\t100\t103")]; 2]);

    // A later A's span can start before the current A's span: A1 (key 100,
    // span [103,106)) misses B [100,102), A2 (key 101, span [101,102))
    // overlaps it. Evicting by A1's span start would lose B.
    let a = format!("{VCF_HEAD}chr1\t101\t.\tAAATTT\tAAA\t.\t.\t.\nchr1\t102\t.\tA\tG\t.\t.\t.\n");
    let (seen, _) = sweep(
        ("a.vcf", &a),
        &[("b.bed", "chr1\t100\t102\n")],
        &SweepConfig::default(),
    )
    .unwrap();
    assert_eq!(
        seen.a_pos,
        [
            (
                100,
                Interval {
                    start: 103,
                    end: 106
                }
            ),
            (
                101,
                Interval {
                    start: 101,
                    end: 102
                }
            )
        ]
    );
    assert_eq!(seen.hits, vec![vec![], vec![hit(0, "chr1\t100\t102")]]);
}

#[test]
fn b_lex_sorted_is_an_error() {
    // A in karyotype order; the second B in lexicographic order (chr10
    // before chr2). The error names that B, with the lex note.
    let a = "chr1\t10\t20\nchr2\t10\t20\nchr10\t10\t20\n";
    let b1 = "chr1\t15\t25\nchr10\t0\t5\n";
    let b2 = "chr1\t15\t25\nchr10\t15\t25\nchr2\t15\t25\n";
    let err = sweep(
        ("a.bed", a),
        &[("b1.bed", b1), ("b2.bed", b2)],
        &SweepConfig::default(),
    )
    .unwrap_err();
    let IoError::Order(e) = err else {
        panic!("expected an order error, got {err}");
    };
    assert_eq!(e.path, "b2.bed");
    assert!(e.hint.starts_with("run: gloci sort b2.bed"), "{}", e.hint);
    assert!(e.hint.contains(LEX_SORTED_NOTE), "{}", e.hint);

    // A is checked the same way: an unsorted A names A.
    let err = sweep(("a.bed", b2), &[("b1.bed", b1)], &SweepConfig::default()).unwrap_err();
    let IoError::Order(e) = err else {
        panic!("expected an order error, got {err}");
    };
    assert_eq!(e.path, "a.bed");
}

#[test]
fn edges_and_empties() {
    let cfg = SweepConfig::default();

    // A on a chromosome B lacks (chr2); B on chromosomes A lacks (chr10
    // between A's chromosomes, chrM after them).
    let a = "chr1\t10\t20\nchr2\t10\t20\nchrX\t5\t6\n";
    let b = "chr1\t15\t25\nchr10\t0\t100\nchrX\t0\t10\nchrM\t0\t10\n";
    let (seen, summary) = sweep(("a.bed", a), &[("b.bed", b)], &cfg).unwrap();
    assert_eq!(
        seen.hits,
        vec![
            vec![hit(0, "chr1\t15\t25")],
            vec![],
            vec![hit(0, "chrX\t0\t10")]
        ]
    );
    assert!(summary.shared_chroms);

    // A nearest request leaves the hits unchanged; the neighbours
    // themselves are tested in tests/nearest.rs.
    let nearest = SweepConfig {
        nearest: Some(NearestConfig {
            k: 1,
            max_distance: None,
            strand: StrandRule::Any,
        }),
        ..cfg
    };
    let (with_nearest, _) = sweep(("a.bed", a), &[("b.bed", b)], &nearest).unwrap();
    assert_eq!(with_nearest.hits, seen.hits);

    // No chromosome in common.
    let (seen, summary) = sweep(
        ("a.bed", "chr2\t0\t10\n"),
        &[("b.bed", "chr1\t0\t10\nchr10\t0\t10\n")],
        &cfg,
    )
    .unwrap();
    assert_eq!(seen.hits, vec![Vec::<(u16, String)>::new()]);
    assert!(!summary.shared_chroms);

    // Empty A: nothing to visit.
    let (seen, summary) = sweep(("a.bed", ""), &[("b.bed", b)], &cfg).unwrap();
    assert!(seen.a.is_empty());
    assert_eq!(summary.a_records, 0);
    assert!(!summary.shared_chroms);

    // Empty B: every A, without hits.
    let (seen, summary) = sweep(("a.bed", a), &[("b.bed", "")], &cfg).unwrap();
    assert_eq!(seen.hits, vec![Vec::<(u16, String)>::new(); 3]);
    assert_eq!(summary.a_records, 3);
    assert_eq!(summary.b_records, [0]);
    assert!(!summary.shared_chroms);

    // Header-only VCF A.
    let (seen, summary) = sweep(("a.vcf", VCF_HEAD), &[("b.bed", b)], &cfg).unwrap();
    assert!(seen.a.is_empty());
    assert_eq!(summary.a_records, 0);

    // A points on B's left and right edges hit only with
    // --inclusive-points; the point inside always does.
    let a = "chr1\t10\t10\nchr1\t15\t15\nchr1\t20\t20\n";
    let b = "chr1\t10\t20\n";
    let (seen, _) = sweep(("a.bed", a), &[("b.bed", b)], &cfg).unwrap();
    assert_eq!(seen.hits, vec![vec![], vec![hit(0, b.trim_end())], vec![]]);
    let inclusive = SweepConfig {
        inclusive_points: true,
        ..cfg
    };
    let (seen, _) = sweep(("a.bed", a), &[("b.bed", b)], &inclusive).unwrap();
    assert_eq!(seen.hits, vec![vec![hit(0, b.trim_end())]; 3]);
}

/// Every field of a hit: file, span, strand, raw line, line number.
type FullHit = (u16, Interval, Strand, String, u64);

/// Each A's raw line with its full hits, for `in_memory_streams`.
#[derive(Default)]
struct Full(Vec<(String, Vec<FullHit>)>);

impl Visitor for Full {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        _near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        let hits = hits
            .iter()
            .map(|h| (h.file, h.span, h.strand, text(h.raw), h.loc))
            .collect();
        self.0.push((text(a.raw), hits));
        Ok(())
    }
}

#[test]
fn in_memory_streams() {
    // Spec §8.1: the sweep takes record streams, here built from text in
    // memory; hits carry their file, span, strand, raw line and line number.
    let o = opts();
    let inputs = Inputs {
        a: mem_reader(
            "a.bed",
            "chr1\t100\t200\ta1\t0\t+\nchr1\t150\t160\ta2\t0\t-\n",
            &o,
        ),
        b: vec![
            mem_reader(
                "x.bed",
                "chr1\t90\t110\tx1\t0\t-\nchr1\t120\t130\tx2\t0\t+\n",
                &o,
            ),
            mem_reader(
                "y.bed",
                "#comment\nchr1\t100\t100\ty1\nchr1\t155\t300\ty2\t0\t.\n",
                &o,
            ),
        ],
    };
    let mut chroms = ChromTable::karyotype();
    let mut full = Full::default();
    let summary = sweep_sorted(inputs, &mut chroms, &SweepConfig::default(), &mut full).unwrap();
    let iv = |start, end| Interval { start, end };
    let y2 = (
        1,
        iv(155, 300),
        Strand::Unknown,
        "chr1\t155\t300\ty2\t0\t.".to_string(),
        3,
    );
    assert_eq!(
        full.0,
        vec![
            (
                "chr1\t100\t200\ta1\t0\t+".to_string(),
                vec![
                    (
                        0,
                        iv(90, 110),
                        Strand::Minus,
                        "chr1\t90\t110\tx1\t0\t-".to_string(),
                        1
                    ),
                    (
                        0,
                        iv(120, 130),
                        Strand::Plus,
                        "chr1\t120\t130\tx2\t0\t+".to_string(),
                        2
                    ),
                    y2.clone(),
                ]
            ),
            ("chr1\t150\t160\ta2\t0\t-".to_string(), vec![y2]),
        ]
    );
    assert_eq!(summary.a_records, 2);
    assert_eq!(summary.b_records, [2, 2]);
    assert!(summary.shared_chroms);
}
