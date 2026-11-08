//! `closest` against the brute-force `naive_closest`, over sorted inputs
//! and shuffled ones with `--unsorted`.

use gloci_core::{ChromOrder, OverlapCriteria, SignFrame, StrandRule, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::closest::{ClosestOptions, Ties, closest};
use gloci_ops::{RelateInputs, RelateOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, NaiveClosest, NaiveTies, naive_closest};
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

#[derive(Debug, Clone)]
struct Params {
    k: usize,
    ties: u8,
    distance: bool,
    signed: Option<u8>,
    ignore_overlaps: bool,
    ignore_upstream: bool,
    ignore_downstream: bool,
    max_distance: Option<u64>,
    each_side: bool,
    per_file: bool,
    strand: u8,
    inclusive_points: bool,
}

fn params() -> impl Strategy<Value = Params> {
    (
        (1usize..=3, 0u8..3, any::<bool>(), prop::option::of(0u8..3)),
        (any::<bool>(), any::<bool>(), any::<bool>()),
        (
            prop::option::of(prop_oneof![0u64..4, 4u64..100, 100u64..100_000]),
            any::<bool>(),
            any::<bool>(),
            prop::sample::select(vec![0u8, 0, 1, 2]),
            any::<bool>(),
        ),
    )
        .prop_map(
            |(
                (k, ties, distance, signed),
                (ignore_overlaps, ignore_upstream, ignore_downstream),
                (max_distance, each_side, per_file, strand, inclusive_points),
            )| Params {
                k,
                ties,
                distance,
                signed,
                ignore_overlaps,
                ignore_upstream,
                ignore_downstream,
                max_distance,
                each_side,
                per_file,
                strand,
                inclusive_points,
            },
        )
}

fn frame(f: u8) -> SignFrame {
    match f {
        0 => SignFrame::Ref,
        1 => SignFrame::A,
        _ => SignFrame::B,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn closest_matches_naive(
        input in relation_inputs(),
        p in params(),
        seed in any::<u64>(),
    ) {
        let strand = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite][usize::from(p.strand)];
        let labels: Vec<String> = (0..input.b.len()).map(|i| format!("L{i}")).collect();
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
            let o = ClosestOptions {
                k: p.k,
                ties: [Ties::All, Ties::First, Ties::Last][usize::from(p.ties)],
                distance: p.distance,
                signed: p.signed.map(frame),
                ignore_overlaps: p.ignore_overlaps,
                ignore_upstream: p.ignore_upstream,
                ignore_downstream: p.ignore_downstream,
                max_distance: p.max_distance,
                each_side: p.each_side,
                per_file: p.per_file,
                relate: RelateOptions {
                    criteria: OverlapCriteria {
                        min_frac_a: None,
                        min_frac_b: None,
                        reciprocal: false,
                        either: false,
                        min_bases: None,
                        strand,
                        inclusive_points: p.inclusive_points,
                    },
                    window: None,
                    split: false,
                    unsorted,
                    sort: None,
                    labels: labels.clone(),
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
            closest(inputs, &o, &mut sink).unwrap();
            let q = NaiveClosest {
                k: p.k,
                ties: [NaiveTies::All, NaiveTies::First, NaiveTies::Last][usize::from(p.ties)],
                distance: p.distance,
                signed: p.signed.map(frame),
                ignore_overlaps: p.ignore_overlaps,
                ignore_upstream: p.ignore_upstream,
                ignore_downstream: p.ignore_downstream,
                max_distance: p.max_distance,
                each_side: p.each_side,
                per_file: p.per_file,
                strand,
                inclusive_points: p.inclusive_points,
                labels: labels.clone(),
            };
            let want = naive_closest(&a, &b, &q);
            prop_assert_eq!(sink.0, want, "unsorted={} {:?}", unsorted, p);
        }
    }
}
