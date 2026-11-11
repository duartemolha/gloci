//! `pair-pair` (D155): orientations, the four types, strands, padding,
//! names, absent ends, thread-count independence and agreement with a
//! brute-force reference.

use gloci_core::{
    ChromTable, Fraction, Interval, OverlapCriteria, Strand, StrandRule, VcfSpan, Window,
};
use gloci_io::{BedpeReader, InvalidPolicy, IoError, ReaderOptions};
use gloci_ops::pairs::{PairPairOptions, PairPairType, PairSummary, pair_pair};
use gloci_testkit::naive::{NaivePair, naive_pair_pair};
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

fn kind(name: &str) -> PairPairType {
    match name {
        "neither" => PairPairType::Neither,
        "either" => PairPairType::Either,
        "both" => PairPairType::Both,
        _ => PairPairType::NotBoth,
    }
}

fn reader(name: &str, text: &str) -> Result<BedpeReader, IoError> {
    BedpeReader::open_from(
        name,
        Box::new(Cursor::new(text.as_bytes().to_vec())),
        &opts(),
    )
}

fn run_with(
    a: &str,
    b: &str,
    o: &PairPairOptions,
    threads: usize,
) -> Result<(String, PairSummary), IoError> {
    let mut chroms = ChromTable::first_appearance();
    let mut out = Vec::new();
    let s = pair_pair(
        reader("a.bedpe", a)?,
        reader("b.bedpe", b)?,
        &mut chroms,
        o,
        NonZeroUsize::new(threads).unwrap(),
        &mut out,
    )?;
    Ok((String::from_utf8(out).unwrap(), s))
}

fn options(k: &str) -> PairPairOptions {
    PairPairOptions {
        kind: kind(k),
        criteria: OverlapCriteria::default(),
        window: None,
        different_names: false,
    }
}

fn run(a: &str, b: &str, k: &str) -> String {
    run_with(a, b, &options(k), 1).unwrap().0
}

const A: &str = "chr1\t100\t200\tchr1\t1000\t1100\ta1\t5\t+\t-\n";

#[test]
fn pair_pair_orientations() {
    // b1 matches in orientation 1, b2 only in orientation 2, b3 on end 1 only.
    let b = "\
chr1\t150\t160\tchr1\t1050\t1060\tb1\t1\t+\t-
chr1\t1050\t1060\tchr1\t150\t160\tb2\t1\t-\t+
chr1\t150\t160\tchr1\t5000\t5010\tb3\t1\t+\t-
";
    let got = run(A, b, "both");
    let want = format!(
        "{a}\t{b1}\n{a}\t{b2}\n",
        a = A.trim_end(),
        b1 = b.lines().next().unwrap(),
        b2 = b.lines().nth(1).unwrap()
    );
    assert_eq!(got, want);
    // A B that matches in both orientations is written once, as it reads.
    let twice = "chr1\t150\t160\tchr1\t150\t160\tu\t1\t+\t+\n";
    let a = "chr1\t100\t200\tchr1\t150\t250\tt\t5\t+\t+\n";
    for k in ["both", "either"] {
        assert_eq!(run(a, twice, k), format!("{}\t{}", a.trim_end(), twice));
    }
}

#[test]
fn pair_pair_types() {
    let a = "\
chr1\t100\t200\tchr1\t1000\t1100\tfull
chr1\t100\t200\tchr2\t500\t600\tone
chr3\t1\t2\tchr3\t5\t6\tnone
";
    let b = "chr1\t150\t160\tchr1\t1050\t1060\tx\nchr1\t180\t190\tchr9\t1\t2\ty\n";
    let l = |s: &str| s.to_string();
    let (full, one, none) = (
        l("chr1\t100\t200\tchr1\t1000\t1100\tfull"),
        l("chr1\t100\t200\tchr2\t500\t600\tone"),
        l("chr3\t1\t2\tchr3\t5\t6\tnone"),
    );
    let (x, y) = (
        l("chr1\t150\t160\tchr1\t1050\t1060\tx"),
        l("chr1\t180\t190\tchr9\t1\t2\ty"),
    );
    let rows = |v: &[String]| v.iter().map(|s| format!("{s}\n")).collect::<String>();
    assert_eq!(run(a, b, "both"), rows(&[format!("{full}\t{x}")]));
    assert_eq!(
        run(a, b, "either"),
        rows(&[
            format!("{full}\t{x}"),
            format!("{full}\t{y}"),
            format!("{one}\t{x}"),
            format!("{one}\t{y}")
        ])
    );
    assert_eq!(run(a, b, "neither"), rows(std::slice::from_ref(&none)));
    assert_eq!(run(a, b, "notboth"), rows(&[one, none]));
}

#[test]
fn pair_pair_window() {
    // x is 10 bases from A's end 1 and 10 from end 2: bookended once padded by 10, so a hit
    // needs 11.
    let b = "chr1\t210\t220\tchr1\t1110\t1120\tx\n";
    let a = "chr1\t100\t200\tchr1\t1000\t1100\tp\n";
    let win = |w: u64, left: u64, right: u64, strand_aware| PairPairOptions {
        window: Some(Window {
            left: left.max(w),
            right: right.max(w),
            strand_aware,
        }),
        ..options("both")
    };
    assert_eq!(run_with(a, b, &win(10, 0, 0, false), 1).unwrap().0, "");
    let got = run_with(a, b, &win(11, 0, 0, false), 1).unwrap().0;
    assert_eq!(got, format!("{}\t{}", a.trim_end(), b));
    // --left 10 --strand-window: a `+` end extends B to the left, a `-` end to the right.
    let sw = PairPairOptions {
        window: Some(Window {
            left: 11,
            right: 0,
            strand_aware: true,
        }),
        ..options("both")
    };
    let plus = "chr1\t100\t200\tchr1\t1000\t1100\tp\t0\t+\t+\n";
    let minus = "chr1\t100\t200\tchr1\t1000\t1100\tp\t0\t-\t-\n";
    assert_eq!(run_with(plus, b, &sw, 1).unwrap().0.lines().count(), 1);
    assert_eq!(run_with(minus, b, &sw, 1).unwrap().0, "");
}

#[test]
fn pair_pair_strands() {
    let a = "chr1\t100\t200\tchr1\t1000\t1100\tp\t0\t+\t-\n";
    let b = "\
chr1\t150\t160\tchr1\t1050\t1060\tsame\t0\t+\t-
chr1\t150\t160\tchr1\t1050\t1060\tflip\t0\t-\t+
chr1\t150\t160\tchr1\t1050\t1060\tunk\t0\t.\t.
";
    let with = |s| PairPairOptions {
        criteria: OverlapCriteria {
            strand: s,
            ..OverlapCriteria::default()
        },
        ..options("both")
    };
    let names = |o: &PairPairOptions| {
        run_with(a, b, o, 1)
            .unwrap()
            .0
            .lines()
            .map(|l| l.split('\t').nth(16).unwrap().to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(names(&options("both")), ["same", "flip", "unk"]);
    assert_eq!(names(&with(StrandRule::Same)), ["same"]);
    assert_eq!(names(&with(StrandRule::Opposite)), ["flip"]);
    // `.` strands never match -s, even `.` with `.`.
    let dots = "chr1\t100\t200\tchr1\t1000\t1100\tp\n";
    let bd = "chr1\t150\t160\tchr1\t1050\t1060\tq\n";
    assert_eq!(
        run_with(dots, bd, &with(StrandRule::Same), 1).unwrap().0,
        ""
    );
    assert_eq!(run(dots, bd, "both").lines().count(), 1);
}

#[test]
fn pair_pair_names() {
    let a = "chr1\t100\t200\tchr1\t1000\t1100\tp\t0\t+\t-\n";
    let b = "\
chr1\t150\t160\tchr1\t1050\t1060\tp\t0\t+\t-
chr1\t150\t160\tchr1\t1050\t1060\tq\t0\t+\t-
";
    let o = PairPairOptions {
        different_names: true,
        ..options("either")
    };
    let got = run_with(a, b, &o, 1).unwrap().0;
    assert_eq!(got.lines().count(), 1);
    assert!(got.contains("\tq\t"));
    // A hit only a same-named B: `neither` writes A.
    let o = PairPairOptions {
        different_names: true,
        ..options("neither")
    };
    let only = "chr1\t150\t160\tchr1\t1050\t1060\tp\t0\t+\t-\n";
    assert_eq!(run_with(a, only, &o, 1).unwrap().0, a);
    // A pair without a name is a data error naming its line.
    let o = PairPairOptions {
        different_names: true,
        ..options("both")
    };
    let e = run_with("chr1\t1\t2\tchr1\t3\t4\n", b, &o, 1).unwrap_err();
    assert!(e.to_string().contains("a.bedpe line 1"), "{e}");
    let e = run_with(a, "chr1\t1\t2\tchr1\t3\t4\n", &o, 1).unwrap_err();
    assert!(e.to_string().contains("b.bedpe line 1"), "{e}");
}

#[test]
fn pair_pair_absent_ends() {
    let a = "chr1\t100\t200\t.\t-1\t-1\ta1\n.\t-1\t-1\t.\t-1\t-1\ta2\n";
    let b = "chr1\t150\t160\t.\t-1\t-1\tb1\n.\t-1\t-1\t.\t-1\t-1\tb2\n";
    // An absent end never hits, not even another absent end.
    assert_eq!(run(a, b, "both"), "");
    assert_eq!(
        run(a, b, "either"),
        "chr1\t100\t200\t.\t-1\t-1\ta1\tchr1\t150\t160\t.\t-1\t-1\tb1\n"
    );
    assert_eq!(run(a, b, "neither"), ".\t-1\t-1\t.\t-1\t-1\ta2\n");
    assert_eq!(run(a, b, "notboth").lines().count(), 2);
}

#[test]
fn pair_pair_threads_match() {
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..3000u64 {
        let s = (i * 7) % 900;
        a.push_str(&format!(
            "chr1\t{s}\t{}\tchr2\t{s}\t{}\tp{i}\n",
            s + 20,
            s + 20
        ));
        if i % 3 == 0 {
            b.push_str(&format!(
                "chr1\t{s}\t{}\tchr2\t{s}\t{}\tq{i}\n",
                s + 5,
                s + 5
            ));
        }
    }
    for k in ["both", "either", "neither", "notboth"] {
        let one = run_with(&a, &b, &options(k), 1).unwrap();
        for t in [2, 3, 4] {
            assert_eq!(run_with(&a, &b, &options(k), t).unwrap(), one, "{k} -t {t}");
        }
    }
}

// --- property test against the brute-force reference -------------------

fn chrom() -> impl Strategy<Value = &'static str> {
    prop::sample::select(vec!["chr1", "chr2"])
}

fn strand() -> impl Strategy<Value = char> {
    prop::sample::select(vec!['+', '-', '.'])
}

type RawEnd = Option<(&'static str, u64, u64, char)>;

fn end() -> impl Strategy<Value = RawEnd> {
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

fn build(prefix: &str, pairs: &[(RawEnd, RawEnd)], names: usize) -> (String, Vec<NaivePair>) {
    let (mut text, mut naive) = (String::new(), Vec::new());
    for (i, (e1, e2)) in pairs.iter().enumerate() {
        let col = |e: &RawEnd| match e {
            Some((c, s, en, _)) => format!("{c}\t{s}\t{en}"),
            None => ".\t-1\t-1".to_string(),
        };
        let st = |e: &RawEnd| e.map_or('.', |x| x.3);
        let line = format!(
            "{}\t{}\t{prefix}{}\t0\t{}\t{}",
            col(e1),
            col(e2),
            i % names,
            st(e1),
            st(e2)
        );
        text.push_str(&line);
        text.push('\n');
        let conv = |e: &RawEnd| {
            e.map(|(c, s, en, st)| (c.to_string(), Interval { start: s, end: en }, sgn(st)))
        };
        naive.push(NaivePair {
            ends: [conv(e1), conv(e2)],
            line,
        });
    }
    (text, naive)
}

proptest! {
    #[test]
    fn pair_pair_matches_naive(
        pairs in proptest::collection::vec((end(), end()), 0..20),
        recs in proptest::collection::vec((end(), end()), 0..20),
        t in 0usize..4,
        rule in 0usize..3,
        inclusive in any::<bool>(),
        pad in proptest::option::of((0u64..6, 0u64..6, any::<bool>())),
        frac in proptest::option::of(prop::sample::select(vec!["0.5", "1"])),
        rdn in any::<bool>(),
        threads in 1usize..4,
    ) {
        let name = ["both", "either", "neither", "notboth"][t];
        let rule = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite][rule];
        // Padding cannot be combined with thresholds (D37).
        let frac = if pad.is_some() { None } else { frac };
        let criteria = OverlapCriteria {
            strand: rule,
            inclusive_points: inclusive,
            min_frac_a: frac.map(|f| Fraction::parse(f).unwrap()),
            ..OverlapCriteria::default()
        };
        let window = pad.map(|(left, right, strand_aware)| Window { left, right, strand_aware });
        let (a_text, a) = build("n", &pairs, 3);
        let (b_text, b) = build("n", &recs, 3);
        let want = naive_pair_pair(&a, &b, name, &criteria, window, rdn);
        let o = PairPairOptions {
            kind: kind(name),
            criteria,
            window,
            different_names: rdn,
        };
        let (got, s) = run_with(&a_text, &b_text, &o, threads).unwrap();
        let want_text: String = want.iter().map(|l| format!("{l}\n")).collect();
        prop_assert_eq!(got, want_text);
        prop_assert_eq!(s.a_pairs, pairs.len() as u64);
    }
}
