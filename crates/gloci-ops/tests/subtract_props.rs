//! `subtract` against the brute-force `naive_subtract`: pieces, `--whole`
//! and `--sum-overlaps`, over sorted inputs and shuffled ones with
//! `--unsorted`.

use gloci_core::{ChromOrder, Fraction, OverlapCriteria, StrandRule, VcfSpan, Window};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::subtract::{SubtractOptions, subtract};
use gloci_ops::{RelateInputs, RelateOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, naive_subtract};
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

fn text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| format!("{}\n", r.bed_line())).collect()
}

/// A deterministic shuffle by `seed`.
fn shuffle(mut v: Vec<TestRec>, seed: u64) -> Vec<TestRec> {
    let key = |i: usize| {
        let mut z = seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z ^ (z >> 27)
    };
    let mut keyed: Vec<(u64, TestRec)> =
        v.drain(..).enumerate().map(|(i, r)| (key(i), r)).collect();
    keyed.sort_by_key(|k| k.0);
    keyed.into_iter().map(|k| k.1).collect()
}

fn fraction() -> impl Strategy<Value = Option<Fraction>> {
    prop::option::of(
        prop::sample::select(vec!["0.1", "0.25", "0.5", "1"])
            .prop_map(|s| Fraction::parse(s).unwrap()),
    )
}

fn criteria() -> impl Strategy<Value = OverlapCriteria> {
    (
        fraction(),
        fraction(),
        any::<bool>(),
        any::<bool>(),
        prop::option::of(0u64..8),
        prop::sample::select(vec![
            StrandRule::Any,
            StrandRule::Any,
            StrandRule::Same,
            StrandRule::Opposite,
        ]),
        any::<bool>(),
    )
        .prop_map(
            |(min_frac_a, min_frac_b, reciprocal, either, min_bases, strand, inclusive_points)| {
                OverlapCriteria {
                    min_frac_a,
                    min_frac_b,
                    reciprocal,
                    either,
                    min_bases,
                    strand,
                    inclusive_points,
                }
            },
        )
}

fn window() -> impl Strategy<Value = Option<Window>> {
    prop_oneof![
        3 => Just(None),
        1 => (0u64..60, 0u64..60, any::<bool>()).prop_map(|(left, right, strand_aware)| {
            Some(Window { left, right, strand_aware })
        }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn subtract_matches_naive(
        input in relation_inputs(),
        mut c in criteria(),
        win in window(),
        kind in 0u8..3,
        seed in any::<u64>(),
    ) {
        let (whole, sum_overlaps) = (kind == 1, kind == 2);
        if win.is_some() {
            // D37: a window takes no thresholds
            c = OverlapCriteria {
                min_frac_a: None,
                min_frac_b: None,
                reciprocal: false,
                either: false,
                min_bases: None,
                ..c
            };
        }
        if sum_overlaps {
            // D60: --sum-overlaps takes only -f
            c = OverlapCriteria {
                min_frac_b: None,
                reciprocal: false,
                either: false,
                min_bases: None,
                ..c
            };
        }
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
            let o = SubtractOptions {
                whole,
                sum_overlaps,
                relate: RelateOptions {
                    criteria: c,
                    window: win,
                    split: false,
                    split_deletions: false,
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
                    .map(|(i, f)| mem_reader(&format!("b{i}.bed"), &text(f), &ro))
                    .collect(),
            };
            let mut sink = Sink::default();
            subtract(inputs, &o, &mut sink).unwrap();
            let want = naive_subtract(&a, &b, &c, win, whole, sum_overlaps);
            prop_assert_eq!(sink.0, want, "unsorted={}", unsorted);
        }
    }
}
