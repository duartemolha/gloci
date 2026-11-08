//! `intersect` against the brute-force `naive_intersect`: any, none, pairs
//! (with overlap bases and unmatched rows) and count, over sorted inputs
//! and over shuffled ones with `--unsorted`.

use gloci_core::{ChromOrder, Fraction, OverlapCriteria, StrandRule, VcfSpan, Window};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::intersect::{IntersectOptions, Match, Mode, intersect};
use gloci_ops::{RelateInputs, RelateOptions};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, NaiveMode, naive_intersect};
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

fn mode() -> impl Strategy<Value = (Mode, NaiveMode)> {
    prop_oneof![
        Just((Mode::Any, NaiveMode::Any)),
        Just((Mode::None, NaiveMode::None)),
        Just((Mode::Count, NaiveMode::Count)),
        (any::<bool>(), any::<bool>()).prop_map(|(overlap_bases, keep_unmatched)| (
            Mode::Pairs {
                overlap_bases,
                keep_unmatched
            },
            NaiveMode::Pairs {
                overlap_bases,
                keep_unmatched
            }
        )),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn intersect_matches_naive(
        input in relation_inputs(),
        mut c in criteria(),
        win in window(),
        (m, nm) in mode(),
        seed in any::<u64>(),
    ) {
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
        let rank = naive::karyotype_rank(CHROMS);
        let a_sorted = naive::sort(input.a.clone(), &rank);
        let b_sorted: Vec<Vec<TestRec>> =
            input.b.iter().map(|f| naive::sort(f.clone(), &rank)).collect();
        let labels: Vec<String> = (0..b_sorted.len()).map(|i| format!("b{i}")).collect();
        let multi = labels.len() > 1;
        for unsorted in [false, true] {
            let (a, b) = if unsorted {
                (
                    shuffle(a_sorted.clone(), seed),
                    b_sorted.iter().map(|f| shuffle(f.clone(), seed ^ 7)).collect::<Vec<_>>(),
                )
            } else {
                (a_sorted.clone(), b_sorted.clone())
            };
            let o = IntersectOptions {
                mode: m,
                matching: Match::Overlap,
                relate: RelateOptions {
                    criteria: c,
                    window: win,
                    split: false,
                    unsorted,
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
            intersect(inputs, &o, &mut sink).unwrap();
            let want = naive_intersect(&a, &b, &c, win, nm, multi.then_some(labels.as_slice()));
            prop_assert_eq!(sink.0, want, "unsorted={}", unsorted);
        }
    }
}

// ---- --match alleles ----

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[derive(Debug, Clone)]
struct V {
    pos: u64,
    reference: String,
    alts: Vec<String>,
}

impl V {
    fn line(&self) -> String {
        format!(
            "chr1\t{}\t.\t{}\t{}\t.\t.\t.",
            self.pos,
            self.reference,
            self.alts.join(",")
        )
    }
}

fn seq() -> impl Strategy<Value = String> {
    prop::collection::vec(prop::sample::select(vec!['A', 'C', 'G']), 1..=3)
        .prop_map(|v| v.into_iter().collect())
}

fn variants() -> impl Strategy<Value = Vec<V>> {
    prop::collection::vec((1u64..8, seq(), prop::collection::vec(seq(), 1..=2)), 0..10).prop_map(
        |v| {
            let mut v: Vec<V> = v
                .into_iter()
                .map(|(pos, reference, alts)| V {
                    pos,
                    reference,
                    alts,
                })
                .collect();
            v.sort_by_key(|r| (r.pos, r.reference.len()));
            v
        },
    )
}

/// (0-based position, REF, ALT) after trimming the shared prefix, then the
/// shared suffix, written out here independently of gloci-core.
fn trim(pos: u64, r: &str, a: &str) -> (u64, String, String) {
    let (r, a) = (r.as_bytes(), a.as_bytes());
    let mut p = 0;
    while p < r.len() && p < a.len() && r[p] == a[p] {
        p += 1;
    }
    let (mut re, mut ae) = (r.len(), a.len());
    while re > p && ae > p && r[re - 1] == a[ae - 1] {
        re -= 1;
        ae -= 1;
    }
    (
        pos - 1 + p as u64,
        String::from_utf8(r[p..re].to_vec()).unwrap(),
        String::from_utf8(a[p..ae].to_vec()).unwrap(),
    )
}

fn shares(a: &V, b: &V) -> bool {
    a.alts.iter().any(|x| {
        b.alts
            .iter()
            .any(|y| trim(a.pos, &a.reference, x) == trim(b.pos, &b.reference, y))
    })
}

fn vcf_text(v: &[V]) -> String {
    let mut t = VCF_HEAD.to_string();
    for r in v {
        t.push_str(&r.line());
        t.push('\n');
    }
    t
}

fn allele_options(mode: Mode, unsorted: bool) -> IntersectOptions {
    IntersectOptions {
        mode,
        matching: Match::Alleles,
        relate: RelateOptions {
            criteria: OverlapCriteria::default(),
            window: None,
            split: false,
            unsorted,
            labels: vec![],
            genome: None,
            chrom_order: ChromOrder::Karyotype,
            header: false,
        },
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// A record is kept when some trimmed ALT of it equals some trimmed ALT
    /// of some B record, under either footprint (--vcf-span) and in sorted
    /// and unsorted mode, however the footprints sit (bookended included).
    #[test]
    fn match_alleles_equals_naive(
        a in variants(),
        b in variants(),
        span_ref in any::<bool>(),
        seed in any::<u64>(),
    ) {
        let ro = ReaderOptions {
            threads: NonZeroUsize::MIN,
            format: None,
            vcf_span: if span_ref { VcfSpan::Ref } else { VcfSpan::Allele },
            invalid: InvalidPolicy::Fail,
        };
        let want_lines: Vec<String> = a
            .iter()
            .filter(|x| b.iter().any(|y| shares(x, y)))
            .map(V::line)
            .collect();
        let want_counts: Vec<String> = a
            .iter()
            .map(|x| format!("{}\t{}", x.line(), b.iter().filter(|y| shares(x, y)).count()))
            .collect();
        for unsorted in [false, true] {
            // unsorted: A and B in a scrambled order, A's order reported
            let scramble = |v: &[V]| -> Vec<V> {
                let mut k: Vec<(u64, V)> = v.iter().cloned().enumerate()
                    .map(|(i, r)| ((seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)).rotate_left(13), r))
                    .collect();
                k.sort_by_key(|x| x.0);
                k.into_iter().map(|x| x.1).collect()
            };
            let (a_in, b_in) = if unsorted { (scramble(&a), scramble(&b)) } else { (a.clone(), b.clone()) };
            let want_lines: Vec<String> = if unsorted {
                a_in.iter().filter(|x| b.iter().any(|y| shares(x, y))).map(V::line).collect()
            } else { want_lines.clone() };
            let want_counts: Vec<String> = if unsorted {
                a_in.iter().map(|x| format!("{}\t{}", x.line(), b.iter().filter(|y| shares(x, y)).count())).collect()
            } else { want_counts.clone() };
            for (mode, want) in [(Mode::Any, &want_lines), (Mode::Count, &want_counts)] {
                let inputs = RelateInputs {
                    a: mem_reader("a.vcf", &vcf_text(&a_in), &ro),
                    b: vec![mem_reader("b.vcf", &vcf_text(&b_in), &ro)],
                };
                let mut sink = Sink::default();
                intersect(inputs, &allele_options(mode, unsorted), &mut sink).unwrap();
                prop_assert_eq!(&sink.0, want, "{:?} unsorted={} ref={}", mode, unsorted, span_ref);
            }
        }
    }
}
