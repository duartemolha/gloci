//! `coverage` against the brute-force `naive_depth`: every report, over
//! sorted inputs and shuffled ones with `--unsorted`, with thresholds,
//! strand rules, points, and several `-b` files.

use gloci_core::{ChromOrder, Fraction, Interval, OverlapCriteria, StrandRule, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::coverage::{CoverageOptions, CoverageReport, coverage, write_hist_totals};
use gloci_ops::{RelateInputs, RelateOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, NaiveBlocks, NaiveDepth, naive_depth};
use gloci_testkit::strategies::{CHROMS, relation_inputs};
use proptest::prelude::*;
use std::collections::BTreeMap;
use std::io;
use std::num::NonZeroUsize;

#[derive(Default)]
struct Sink(Vec<String>);

impl RecordWriter for Sink {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        self.0.push(String::from_utf8(raw.to_vec()).unwrap());
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), gloci_io::IoError> {
        Ok(())
    }
}

fn reader_opts() -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

/// The blocks of a record, derived from its coordinates: two blocks with a
/// gap for some records of three or more bases, else the span as one block.
fn blocks_of(r: &TestRec) -> Vec<Interval> {
    let len = r.end - r.start;
    if len >= 3 && !(r.start + r.end).is_multiple_of(3) {
        let l1 = 1 + r.start % (len - 2);
        let l2 = 1 + r.end % (len - l1 - 1);
        vec![
            Interval {
                start: r.start,
                end: r.start + l1,
            },
            Interval {
                start: r.end - l2,
                end: r.end,
            },
        ]
    } else {
        vec![Interval {
            start: r.start,
            end: r.end,
        }]
    }
}

/// The records as BED12 lines (a point as BED6, which has no blocks).
fn text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| format!("{}\n", line(r))).collect()
}

fn line(r: &TestRec) -> String {
    if r.start == r.end {
        return r.bed_line();
    }
    let b = blocks_of(r);
    let sizes: Vec<String> = b.iter().map(|x| (x.end - x.start).to_string()).collect();
    let starts: Vec<String> = b.iter().map(|x| (x.start - r.start).to_string()).collect();
    format!(
        "{}\t{}\t{}\t0\t{}\t{}\t{}",
        r.bed_line(),
        r.start,
        r.end,
        b.len(),
        sizes.join(","),
        starts.join(",")
    )
}

fn shuffle(v: Vec<TestRec>, seed: u64) -> Vec<TestRec> {
    let key = |i: usize| {
        let mut z = seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z ^ (z >> 27)
    };
    let mut keyed: Vec<(u64, TestRec)> = v
        .into_iter()
        .enumerate()
        .map(|(i, r)| (key(i), r))
        .collect();
    keyed.sort_by_key(|k| k.0);
    keyed.into_iter().map(|k| k.1).collect()
}

fn ratio(n: u64, d: u64) -> String {
    (n as f64 / d as f64).to_string()
}

/// The rows gloci should write for `report`, from the naive depth profiles
/// (single `-b` file), and the `--hist` totals.
fn expected(a: &[TestRec], depths: &[NaiveDepth], report: CoverageReport) -> Vec<String> {
    let mut rows = Vec::new();
    for (ra, nd) in a.iter().zip(depths) {
        let line = line(ra);
        let len = nd.len();
        let covered = nd.covered();
        match report {
            CoverageReport::Summary => {
                let frac = if len == 0 {
                    ".".into()
                } else {
                    ratio(covered, len)
                };
                rows.push(format!("{line}\t{}\t{covered}\t{len}\t{frac}", nd.count));
            }
            CoverageReport::Mean => {
                let v = if len == 0 {
                    ".".into()
                } else {
                    ratio(nd.total_depth(), len)
                };
                rows.push(format!("{line}\t{v}"));
            }
            CoverageReport::PerBase => {
                for (i, d) in nd.per_base().iter().enumerate() {
                    rows.push(format!("{line}\t{}\t{d}", i + 1));
                }
            }
            CoverageReport::Hist => {
                if len == 0 {
                    continue;
                }
                let mut here: BTreeMap<u64, u64> = BTreeMap::new();
                for (s, e, d) in &nd.runs {
                    *here.entry(*d).or_default() += e - s;
                }
                for (d, b) in here {
                    rows.push(format!("{line}\t{d}\t{b}\t{len}\t{}", ratio(b, len)));
                }
            }
        }
    }
    rows
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn coverage_matches_naive(
        input in relation_inputs(),
        strand in 0u8..3,
        frac in prop::option::of(1u32..=100),
        reciprocal in any::<bool>(),
        either in any::<bool>(),
        min_bases in prop::option::of(0u64..40),
        split in any::<bool>(),
        inclusive_points in any::<bool>(),
        seed in any::<u64>(),
    ) {
        let strand = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite][usize::from(strand)];
        let criteria = OverlapCriteria {
            min_frac_a: frac.map(|p| Fraction::parse(&format!("{p}%")).unwrap()),
            min_frac_b: None,
            reciprocal,
            either,
            min_bases,
            strand,
            inclusive_points,
        };
        let rank = naive::karyotype_rank(CHROMS);
        let a_sorted = naive::sort(input.a.clone(), &rank);
        let b_sorted: Vec<Vec<TestRec>> =
            input.b.iter().map(|f| naive::sort(f.clone(), &rank)).collect();
        for unsorted in [false, true] {
            let (a, b) = if unsorted {
                (
                    shuffle(a_sorted.clone(), seed),
                    b_sorted.iter().map(|f| shuffle(f.clone(), seed ^ 7)).collect::<Vec<_>>(),
                )
            } else {
                (a_sorted.clone(), b_sorted.clone())
            };
            let run = |report: CoverageReport, files: &[Vec<TestRec>]| {
                let o = CoverageOptions {
                    report,
                    relate: RelateOptions {
                        criteria,
                        window: None,
                        split,
                        unsorted,
                        labels: Vec::new(),
                        genome: None,
                        chrom_order: ChromOrder::Karyotype,
                        header: false,
                    },
                };
                let ro = reader_opts();
                let inputs = RelateInputs {
                    a: mem_reader("a.bed", &text(&a), &ro),
                    b: files.iter().enumerate()
                        .map(|(i, f)| mem_reader(&format!("b{i}.bed"), &text(f), &ro))
                        .collect(),
                };
                let mut sink = Sink::default();
                let (_, totals) = coverage(inputs, &o, &mut sink).unwrap();
                (sink.0, totals)
            };
            let nb: Option<NaiveBlocks> = if split { Some(&blocks_of) } else { None };
            // every report over the first file alone
            let one = &b[..1];
            let depths = naive_depth(&a, one, &criteria, nb);
            for report in [
                CoverageReport::Summary,
                CoverageReport::Hist,
                CoverageReport::PerBase,
                CoverageReport::Mean,
            ] {
                // per-base writes a row per base: only for short A
                if report == CoverageReport::PerBase && depths.iter().any(|d| d.len() > 300) {
                    continue;
                }
                let (mut got, totals) = run(report, one);
                if report == CoverageReport::Hist {
                    let mut w = Sink::default();
                    write_hist_totals(&totals, &mut w).unwrap();
                    got.extend(w.0);
                }
                let mut want = expected(&a, &depths, report);
                if report == CoverageReport::Hist {
                    let mut by_depth: BTreeMap<u64, u64> = BTreeMap::new();
                    let mut len = 0u64;
                    for nd in depths.iter().filter(|d| !d.is_empty()) {
                        len += nd.len();
                        for (s, e, d) in &nd.runs {
                            *by_depth.entry(*d).or_default() += e - s;
                        }
                    }
                    if len > 0 {
                        for (d, n) in by_depth {
                            want.push(format!("all\t{d}\t{n}\t{len}\t{}", ratio(n, len)));
                        }
                    }
                }
                prop_assert_eq!(got, want, "{:?} unsorted={}", report, unsorted);
            }
            // several files: one covered fraction per file
            let (got, _) = run(CoverageReport::Summary, &b);
            let per_file: Vec<Vec<NaiveDepth>> =
                b.iter().map(|f| naive_depth(&a, std::slice::from_ref(f), &criteria, nb)).collect();
            if b.len() > 1 {
                let want: Vec<String> = a.iter().enumerate().map(|(i, ra)| {
                    let mut row = line(ra);
                    for pf in &per_file {
                        let d = &pf[i];
                        row.push('\t');
                        if d.is_empty() {
                            row.push('.');
                        } else {
                            row.push_str(&ratio(d.covered(), d.len()));
                        }
                    }
                    row
                }).collect();
                prop_assert_eq!(got, want, "multi unsorted={}", unsorted);
            }
        }
    }
}

#[test]
fn totals_merge_sums_chunks() {
    let one = |a: &str, b: &str| {
        let ro = reader_opts();
        let o = CoverageOptions {
            report: CoverageReport::Hist,
            relate: RelateOptions {
                criteria: OverlapCriteria::default(),
                window: None,
                split: false,
                unsorted: false,
                labels: Vec::new(),
                genome: None,
                chrom_order: ChromOrder::Karyotype,
                header: false,
            },
        };
        let inputs = RelateInputs {
            a: mem_reader("a.bed", a, &ro),
            b: vec![mem_reader("b.bed", b, &ro)],
        };
        coverage(inputs, &o, &mut Sink::default()).unwrap().1
    };
    let t1 = one("chr1\t0\t10\n", "chr1\t0\t4\n");
    let t2 = one("chr2\t0\t10\n", "chr2\t2\t4\nchr2\t3\t8\n");
    let whole = one(
        "chr1\t0\t10\nchr2\t0\t10\n",
        "chr1\t0\t4\nchr2\t2\t4\nchr2\t3\t8\n",
    );
    let mut m = t1.clone();
    m.merge(&t2);
    assert_eq!(m, whole);
    assert_eq!(
        (m.bases_at(0), m.bases_at(1), m.bases_at(2), m.length()),
        (10, 9, 1, 20)
    );
}
