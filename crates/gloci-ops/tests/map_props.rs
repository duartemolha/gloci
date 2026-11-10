//! `map` against the brute-force `naive_map`: sum, count, min, max, bases
//! and indicator, over sorted inputs and shuffled ones with `--unsorted`,
//! with thresholds, strand rules and windows.

use gloci_core::{ChromOrder, Fraction, OverlapCriteria, StrandRule, VcfSpan, Window};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::agg::parse_aggs;
use gloci_ops::map::{MapOptions, map};
use gloci_ops::{RelateInputs, RelateOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, NaiveAgg, naive_map, scored_line};
use gloci_testkit::strategies::{CHROMS, relation_inputs};
use proptest::prelude::*;
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
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

/// B's score column: missing for a quarter of the records, and repeated
/// values, so `distinct` has something to drop.
fn score(r: &TestRec) -> Option<u64> {
    match r.end % 4 {
        0 => None,
        k => Some(k * 10 + r.start % 3),
    }
}

fn text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| format!("{}\n", r.bed_line())).collect()
}

fn b_text(recs: &[TestRec]) -> String {
    recs.iter()
        .map(|r| format!("{}\n", scored_line(r, score(r))))
        .collect()
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

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn map_matches_naive(
        input in relation_inputs(),
        strand in 0u8..3,
        frac in prop::option::of(1u32..=100),
        reciprocal in any::<bool>(),
        either in any::<bool>(),
        min_bases in prop::option::of(0u64..40),
        split in any::<bool>(),
        window in prop::option::of((0u64..30, 0u64..30, any::<bool>())),
        inclusive_points in any::<bool>(),
        skip_unmatched in any::<bool>(),
        seed in any::<u64>(),
    ) {
        let strand = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite][usize::from(strand)];
        // D37: a window and thresholds do not mix
        let (frac, reciprocal, either, min_bases) = if window.is_some() {
            (None, false, false, None)
        } else {
            (frac, reciprocal, either, min_bases)
        };
        let criteria = OverlapCriteria {
            min_frac_a: frac.map(|p| Fraction::parse(&format!("{p}%")).unwrap()),
            min_frac_b: None,
            reciprocal,
            either,
            min_bases,
            strand,
            inclusive_points,
        };
        let window = window.map(|(left, right, strand_aware)| Window { left, right, strand_aware });
        let cols = [
            "0:count", "3:sum", "2:min", "3:max", "0:bases", "0:indicator", "2:sum",
            "5:sum", "5:min", "5:max", "5:collapse", "5:distinct", "5:first", "5:last",
            "4:collapse", "4:last",
        ];
        let naive_aggs = [
            NaiveAgg::Count, NaiveAgg::Sum(3), NaiveAgg::Min(2), NaiveAgg::Max(3),
            NaiveAgg::Bases, NaiveAgg::Indicator, NaiveAgg::Sum(2),
            NaiveAgg::Sum(5), NaiveAgg::Min(5), NaiveAgg::Max(5), NaiveAgg::Collapse(5),
            NaiveAgg::Distinct(5), NaiveAgg::First(5), NaiveAgg::Last(5),
            NaiveAgg::Collapse(4), NaiveAgg::Last(4),
        ];
        let aggs = parse_aggs(&cols.map(String::from)).unwrap();
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
            let o = MapOptions {
                aggs: aggs.clone(),
                skip_unmatched,
                empty: ".".into(),
                delim: ",".into(),
                relate: RelateOptions {
                    criteria,
                    window,
                    split,
                    unsorted,
                    sort: None,
                    labels: Vec::new(),
                    genome: None,
                    chrom_order: ChromOrder::Karyotype,
                    header: false,
                },
            };
            let ro = reader_opts();
            let inputs = RelateInputs {
                a: mem_reader("a.bed", &text(&a), &ro),
                b: b.iter().enumerate()
                    .map(|(i, f)| mem_reader(&format!("b{i}.bed"), &b_text(f), &ro))
                    .collect(),
            };
            let mut sink = Sink::default();
            map(inputs, &o, &mut sink).unwrap();
            let want = naive_map(&a, &b, &criteria, window, &naive_aggs, skip_unmatched, &score);
            prop_assert_eq!(sink.0, want, "unsorted={}", unsorted);
        }
    }
}
