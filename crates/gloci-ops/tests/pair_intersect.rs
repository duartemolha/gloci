//! `pair-intersect` (D154): the nine types, absent ends, strands, span edges,
//! padding, VCF footprints, unsorted B, thread-count independence and
//! agreement with a brute-force reference.

use gloci_core::{ChromTable, Fraction, Interval, OverlapCriteria, Strand, StrandRule, VcfSpan};
use gloci_io::{BedpeReader, InvalidPolicy, ReaderOptions};
use gloci_ops::pairs::{PairIntersectOptions, PairSummary, PairType, pair_intersect};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{NaivePair, naive_pair_intersect};
use proptest::prelude::*;
use std::io::Cursor;
use std::num::NonZeroUsize;

fn opts() -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::new(1).unwrap(),
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

const TYPES: [(&str, PairType); 9] = [
    ("either", PairType::Either),
    ("neither", PairType::Neither),
    ("both", PairType::Both),
    ("xor", PairType::Xor),
    ("notboth", PairType::NotBoth),
    ("ispan", PairType::InnerSpan),
    ("ospan", PairType::OuterSpan),
    ("notispan", PairType::NotInnerSpan),
    ("notospan", PairType::NotOuterSpan),
];

fn kind(name: &str) -> PairType {
    TYPES.iter().find(|(n, _)| *n == name).unwrap().1
}

fn run_with(
    a: &str,
    b: &str,
    kind: PairType,
    criteria: OverlapCriteria,
    threads: usize,
) -> Result<(String, PairSummary), gloci_io::IoError> {
    let a = BedpeReader::open_from(
        "a.bedpe",
        Box::new(Cursor::new(a.as_bytes().to_vec())),
        &opts(),
    )?;
    let b = mem_reader("b", b, &opts());
    let mut chroms = ChromTable::first_appearance();
    let mut out = Vec::new();
    let o = PairIntersectOptions { kind, criteria };
    let s = pair_intersect(
        a,
        b,
        &mut chroms,
        &o,
        NonZeroUsize::new(threads).unwrap(),
        &mut out,
    )?;
    Ok((String::from_utf8(out).unwrap(), s))
}

fn run(a: &str, b: &str, kind: &str) -> (String, PairSummary) {
    run_with(a, b, self::kind(kind), OverlapCriteria::default(), 1).unwrap()
}

const B: &str = "chr1\t10\t20\tb1\nchr2\t5\t25\tb2\n";
// p1 hits b1 with end 1 only, p2 hits b2 with end 2 only, p3 hits both, p4 neither,
// p5 is on one chromosome and touches nothing.
const A: &str = "\
chr1\t12\t18\tchr1\t1000\t1010\tp1
chr1\t500\t510\tchr2\t10\t20\tp2
chr1\t11\t13\tchr2\t8\t9\tp3
chr1\t500\t510\tchr2\t500\t510\tp4
chr1\t100\t110\tchr1\t200\t210\tp5
";

fn rows(lines: &[&str]) -> String {
    lines.iter().map(|l| format!("{l}\n")).collect()
}

#[test]
fn pair_type_table() {
    let p = |n: &str| match n {
        "p1" => "chr1\t12\t18\tchr1\t1000\t1010\tp1",
        "p2" => "chr1\t500\t510\tchr2\t10\t20\tp2",
        "p3" => "chr1\t11\t13\tchr2\t8\t9\tp3",
        "p4" => "chr1\t500\t510\tchr2\t500\t510\tp4",
        _ => "chr1\t100\t110\tchr1\t200\t210\tp5",
    };
    let (b1, b2) = ("chr1\t10\t20\tb1", "chr2\t5\t25\tb2");
    let ab = |n: &str, b: &str| format!("{}\t{b}", p(n));
    let want: [(&str, Vec<String>); 9] = [
        (
            "either",
            vec![ab("p1", b1), ab("p2", b2), ab("p3", b1), ab("p3", b2)],
        ),
        ("neither", vec![p("p4").into(), p("p5").into()]),
        ("both", vec![ab("p3", b1), ab("p3", b2)]),
        ("xor", vec![ab("p1", b1), ab("p2", b2)]),
        (
            "notboth",
            vec![
                ab("p1", b1),
                ab("p2", b2),
                ab("p4", ".\t-1\t-1\t."),
                ab("p5", ".\t-1\t-1\t."),
            ],
        ),
        // p1: gap [18, 1000) and span [12, 1010) both meet b1. The pairs on two
        // chromosomes (p2, p3, p4) are left out.
        ("ispan", vec![ab("p1", b1)]),
        ("ospan", vec![ab("p1", b1)]),
        ("notispan", vec![p("p5").into()]),
        ("notospan", vec![p("p5").into()]),
    ];
    for (name, lines) in want {
        let (got, s) = run(A, B, name);
        let want = rows(&lines.iter().map(String::as_str).collect::<Vec<_>>());
        assert_eq!(got, want, "--type {name}");
        assert_eq!(s.a_pairs, 5, "{name}");
        let span = name.contains("span");
        assert_eq!(s.span_skipped, if span { 3 } else { 0 }, "{name}");
    }
}

#[test]
fn bedpe_absent_ends() {
    // Pair a: end 1 absent, end 2 hits b1. Pair b: end 1 absent, end 2 misses.
    let a = ".\t-1\t-1\tchr1\t12\t18\ta\n.\t-1\t-1\tchr1\t500\t510\tb\n";
    let b = "chr1\t10\t20\tb1\n";
    let (got, s) = run(a, b, "either");
    assert_eq!(got, ".\t-1\t-1\tchr1\t12\t18\ta\tchr1\t10\t20\tb1\n");
    assert_eq!(s.a_pairs, 2);
    assert_eq!(run(a, b, "both").0, "");
    assert_eq!(
        run(a, b, "xor").0,
        ".\t-1\t-1\tchr1\t12\t18\ta\tchr1\t10\t20\tb1\n"
    );
    assert_eq!(run(a, b, "neither").0, ".\t-1\t-1\tchr1\t500\t510\tb\n");
    assert_eq!(
        run(a, b, "notboth").0,
        ".\t-1\t-1\tchr1\t12\t18\ta\tchr1\t10\t20\tb1\n\
         .\t-1\t-1\tchr1\t500\t510\tb\t.\t-1\t-1\t.\n"
    );
    // Span types skip both, and count them.
    for t in ["ispan", "ospan", "notispan", "notospan"] {
        let (got, s) = run(a, b, t);
        assert_eq!(got, "", "{t}");
        assert_eq!(s.span_skipped, 2, "{t}");
    }
    // Both ends absent: neither reports it, notboth pads it.
    let none = ".\t-1\t-1\t.\t-1\t-1\n";
    assert_eq!(run(none, b, "neither").0, none);
    assert_eq!(run(none, b, "either").0, "");
}

fn strands(rule: StrandRule) -> OverlapCriteria {
    OverlapCriteria {
        strand: rule,
        ..OverlapCriteria::default()
    }
}

#[test]
fn pair_strands_use_each_ends_own() {
    // Strand 1 is +, strand 2 is -. B: b1 + near end 1, b2 + near end 2.
    let a = "chr1\t10\t20\tchr1\t100\t110\tp\t0\t+\t-\n";
    let b = "chr1\t12\t14\tb1\t0\t+\nchr1\t102\t104\tb2\t0\t+\n";
    let line = "chr1\t10\t20\tchr1\t100\t110\tp\t0\t+\t-";
    let go = |rule| {
        run_with(a, b, PairType::Either, strands(rule), 1)
            .unwrap()
            .0
    };
    assert_eq!(
        go(StrandRule::Same),
        format!("{line}\tchr1\t12\t14\tb1\t0\t+\n")
    );
    assert_eq!(
        go(StrandRule::Opposite),
        format!("{line}\tchr1\t102\t104\tb2\t0\t+\n")
    );
}

#[test]
fn pair_span_types_edges() {
    let b = "chr1\t1\t2\tk1\nchr1\t5\t6\tk5\nchr1\t7\t8\tk7\nchr1\t9\t10\tk9\nchr1\t12\t13\tk12\n";
    let go = |a: &str, t: &str| run(a, b, t).0;
    // End 1 left of end 2: the gap is [5, 10).
    let left = "chr1\t2\t5\tchr1\t10\t14\tl\n";
    // End 1 right of end 2: the same gap.
    let right = "chr1\t10\t14\tchr1\t2\t5\tr\n";
    for (a, name) in [(left, "l"), (right, "r")] {
        let got = go(a, "ispan");
        assert_eq!(got.lines().count(), 3, "{name}: {got}");
        assert!(
            got.contains("k5") && got.contains("k7") && got.contains("k9"),
            "{name}"
        );
        assert_eq!(go(a, "notispan"), "");
        // The outer span is [2, 14): all but k1.
        let outer = go(a, "ospan");
        assert_eq!(outer.lines().count(), 4, "{name}");
        assert!(!outer.contains("k1\n"), "{name}");
    }
    // Touching ends: no gap, nothing hits it, and notispan reports the pair.
    let touch = "chr1\t2\t5\tchr1\t5\t8\tt\n";
    assert_eq!(go(touch, "ispan"), "");
    assert_eq!(go(touch, "notispan"), touch);
    // Overlapping ends have no gap either.
    let over = "chr1\t2\t8\tchr1\t5\t12\to\n";
    assert_eq!(go(over, "ispan"), "");
    assert_eq!(go(over, "notispan"), over);
    // Nested ends: the outer span is the container, [2, 12).
    let nested = "chr1\t2\t12\tchr1\t5\t8\tn\n";
    assert_eq!(go(nested, "ospan").lines().count(), 3);
    // Two chromosomes: skipped and counted, in either direction.
    let two = "chr1\t2\t5\tchr2\t10\t14\tx\n";
    for t in ["ispan", "ospan", "notispan", "notospan"] {
        let (got, s) = run(two, b, t);
        assert_eq!((got.as_str(), s.span_skipped), ("", 1), "{t}");
    }
    // -s with a span type is a usage error.
    let e = run_with(left, b, PairType::InnerSpan, strands(StrandRule::Same), 1).unwrap_err();
    assert!(matches!(e, gloci_io::IoError::Usage(_)), "{e}");
    assert!(e.to_string().contains("ispan"), "{e}");
}

#[test]
fn notboth_padding_follows_b() {
    let a = "chr1\t500\t510\tchr1\t600\t610\tp\n";
    // A BED6 B pads with a score of -1; a GFF B with five fixed columns.
    let bed6 = "chr1\t1\t2\tn\t0\t+\n";
    assert_eq!(
        run(a, bed6, "notboth").0,
        "chr1\t500\t510\tchr1\t600\t610\tp\t.\t-1\t-1\t.\t-1\t.\n"
    );
    let gff = "chr1\tsrc\tgene\t1\t2\t.\t+\t.\tID=g\n";
    assert_eq!(
        run(a, gff, "notboth").0,
        "chr1\t500\t510\tchr1\t600\t610\tp\t.\t.\t.\t-1\t-1\t.\t.\t.\t.\n"
    );
    // An empty B: three columns.
    assert_eq!(
        run(a, "", "notboth").0,
        "chr1\t500\t510\tchr1\t600\t610\tp\t.\t-1\t-1\n"
    );
}

#[test]
fn vcf_b_hits_its_footprint() {
    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t101\tv1\tAT\tA\t.\t.\t.\n\
chr1\t300\tv2\t<DEL>\t<DEL>\t.\t.\tEND=320\n";
    // v1 deletes base 102 (0-based 101..102); v2 spans 299..320.
    let a = "chr1\t100\t105\tchr1\t310\t312\tp\nchr1\t0\t50\tchr1\t60\t70\tq\n";
    let (got, _) = run_with(a, vcf, PairType::Either, OverlapCriteria::default(), 1).unwrap();
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines.len(), 2, "{got}");
    assert!(lines[0].starts_with("chr1\t100\t105\tchr1\t310\t312\tp\tchr1\t101\tv1\tAT\tA"));
    assert!(lines[1].contains("\tv2\t<DEL>"));
}

#[test]
fn unsorted_b_is_accepted_and_hits_keep_b_order() {
    let b = "chr2\t5\t25\tz\nchr1\t100\t200\ty\nchr1\t10\t20\tx\nchr1\t12\t30\tw\n";
    let a = "chr1\t15\t16\tchr1\t150\t151\tp\n";
    let (got, _) = run(a, b, "either");
    let names: Vec<&str> = got
        .lines()
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    // End 1 hits x then w (B order), end 2 hits y.
    assert_eq!(names, ["x", "w", "y"]);
}

#[test]
fn fraction_applies_to_each_end() {
    let a = "chr1\t10\t20\tchr1\t100\t110\tp\n";
    let b = "chr1\t10\t13\tb1\nchr1\t100\t106\tb2\n";
    let c = OverlapCriteria {
        min_frac_a: Some(Fraction::parse("0.5").unwrap()),
        ..OverlapCriteria::default()
    };
    // b1 covers 30% of end 1, b2 covers 60% of end 2.
    let (got, _) = run_with(a, b, PairType::Either, c, 1).unwrap();
    assert_eq!(got, "chr1\t10\t20\tchr1\t100\t110\tp\tchr1\t100\t106\tb2\n");
}

fn big_inputs(n: usize) -> (String, String) {
    let mut a = String::new();
    for i in 0..n {
        let s = (i * 37) % 9000;
        a.push_str(&format!(
            "chr{}\t{s}\t{}\tchr{}\t{}\t{}\tp{i}\n",
            i % 3 + 1,
            s + 50,
            (i * 7) % 3 + 1,
            (i * 53) % 9000,
            (i * 53) % 9000 + 40
        ));
    }
    let mut b = String::new();
    for i in 0..300 {
        let s = (i * 31) % 9000;
        b.push_str(&format!("chr{}\t{s}\t{}\tb{i}\n", i % 3 + 1, s + 100));
    }
    (a, b)
}

#[test]
fn threads_do_not_change_the_output() {
    // Enough pairs for several batches.
    let (a, b) = big_inputs(60_000);
    for t in ["either", "notboth", "ospan", "neither"] {
        let one = run_with(&a, &b, kind(t), OverlapCriteria::default(), 1).unwrap();
        for threads in [2, 3, 8] {
            let many = run_with(&a, &b, kind(t), OverlapCriteria::default(), threads).unwrap();
            assert!(one.0 == many.0, "--type {t} at {threads} threads");
            assert_eq!(one.1, many.1);
        }
        assert!(t == "neither" || !one.0.is_empty(), "{t}");
    }
}

#[test]
fn pair_huge_contig_count() {
    // 200,000 contigs, one pair and one B record on each.
    let n = 200_000;
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..n {
        a.push_str(&format!("c{i}\t10\t20\tc{i}\t100\t110\tp{i}\n"));
        b.push_str(&format!("c{i}\t15\t16\tb{i}\n"));
    }
    let (got, s) = run_with(&a, &b, PairType::Both, OverlapCriteria::default(), 1).unwrap();
    // End 1 hits b; end 2 does not.
    assert_eq!(got, "");
    assert_eq!(s.a_pairs, n);
    let (got, _) = run_with(&a, &b, PairType::Xor, OverlapCriteria::default(), 2).unwrap();
    assert_eq!(got.lines().count(), n as usize);
    assert!(got.starts_with("c0\t10\t20\tc0\t100\t110\tp0\tc0\t15\t16\tb0\n"));
}

// --- property test against the brute-force reference -------------------

fn chrom() -> impl Strategy<Value = &'static str> {
    prop::sample::select(vec!["chr1", "chr2", "chr3"])
}

fn strand() -> impl Strategy<Value = char> {
    prop::sample::select(vec!['+', '-', '.'])
}

fn end() -> impl Strategy<Value = Option<(&'static str, u64, u64, char)>> {
    prop_oneof![
        1 => Just(None),
        9 => (chrom(), 0u64..40, 0u64..10, strand()).prop_map(|(c, s, l, st)| Some((c, s, s + l, st))),
    ]
}

fn sgn(c: char) -> Strand {
    match c {
        '+' => Strand::Plus,
        '-' => Strand::Minus,
        _ => Strand::Unknown,
    }
}

proptest! {
    #[test]
    fn pair_intersect_matches_naive(
        pairs in proptest::collection::vec((end(), end()), 0..25),
        recs in proptest::collection::vec((chrom(), 0u64..40, 0u64..10, strand()), 0..25),
        t in 0usize..9,
        rule in 0usize..3,
        inclusive in any::<bool>(),
        frac in proptest::option::of(prop::sample::select(vec!["0.5", "0.2", "1"])),
        threads in 1usize..4,
    ) {
        let (name, ptype) = TYPES[t];
        let rule = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite][rule];
        // The span types refuse a strand rule.
        let rule = if ptype.is_span() { StrandRule::Any } else { rule };
        let criteria = OverlapCriteria {
            strand: rule,
            inclusive_points: inclusive,
            min_frac_a: frac.map(|f| Fraction::parse(f).unwrap()),
            ..OverlapCriteria::default()
        };
        let mut a_text = String::new();
        let mut naive_a = Vec::new();
        for (i, (e1, e2)) in pairs.iter().enumerate() {
            let col = |e: &Option<(&str, u64, u64, char)>| match e {
                Some((c, s, en, _)) => format!("{c}\t{s}\t{en}"),
                None => ".\t-1\t-1".to_string(),
            };
            let st = |e: &Option<(&str, u64, u64, char)>| e.map_or('.', |x| x.3);
            let line = format!("{}\t{}\tp{i}\t0\t{}\t{}", col(e1), col(e2), st(e1), st(e2));
            a_text.push_str(&line);
            a_text.push('\n');
            let conv = |e: &Option<(&str, u64, u64, char)>| {
                e.map(|(c, s, en, st)| (c.to_string(), Interval { start: s, end: en }, sgn(st)))
            };
            naive_a.push(NaivePair { ends: [conv(e1), conv(e2)], line });
        }
        let b_recs: Vec<TestRec> = recs
            .iter()
            .enumerate()
            .map(|(i, (c, s, l, st))| TestRec {
                chrom: c.to_string(),
                start: *s,
                end: s + l,
                strand: *st,
                name: format!("b{i}"),
            })
            .collect();
        let b_text: String = b_recs.iter().map(|r| format!("{}\n", r.bed_line())).collect();
        let (want, skipped) = naive_pair_intersect(&naive_a, &b_recs, name, &criteria);
        let (got, s) = run_with(&a_text, &b_text, ptype, criteria, threads).unwrap();
        let want_text: String = want.iter().map(|l| format!("{l}\n")).collect();
        prop_assert_eq!(got, want_text);
        prop_assert_eq!(s.span_skipped, skipped);
        prop_assert_eq!(s.a_pairs, pairs.len() as u64);
    }
}
