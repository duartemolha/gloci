//! `stats jaccard` against the brute-force `naive_jaccard`: points, bookends,
//! nesting, thresholds and strand rules, over sorted and `--sort`ed inputs.

use gloci_core::{ChromOrder, Fraction, OverlapCriteria, StrandRule, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions};
use gloci_ops::RelateInputs;
use gloci_ops::stats::{JaccardOptions, jaccard};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, naive_jaccard};
use gloci_testkit::strategies::{CHROMS, relation_inputs};
use proptest::prelude::*;
use std::num::NonZeroUsize;

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

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn jaccard_matches_naive(
        input in relation_inputs(),
        strand in 0u8..3,
        frac_a in prop::option::of(1u32..=100),
        frac_b in prop::option::of(1u32..=100),
        reciprocal in any::<bool>(),
        either in any::<bool>(),
        min_bases in prop::option::of(0u64..40),
        inclusive_points in any::<bool>(),
    ) {
        let strand = [StrandRule::Any, StrandRule::Same, StrandRule::Opposite][usize::from(strand)];
        let criteria = OverlapCriteria {
            min_frac_a: frac_a.map(|p| Fraction::parse(&format!("{p}%")).unwrap()),
            min_frac_b: frac_b.map(|p| Fraction::parse(&format!("{p}%")).unwrap()),
            reciprocal,
            either,
            min_bases,
            strand,
            inclusive_points,
        };
        let rank = naive::karyotype_rank(CHROMS);
        let a = naive::sort(input.a.clone(), &rank);
        let b = naive::sort(input.b[0].clone(), &rank);
        let want = naive_jaccard(&a, &b, &criteria);
        let o = JaccardOptions {
            criteria,
            sort: None,
            genome: None,
            chrom_order: ChromOrder::Karyotype,
        };
        let ro = reader_opts();
        let got = jaccard(
            RelateInputs {
                a: mem_reader("a.bed", &text(&a), &ro),
                b: vec![mem_reader("b.bed", &text(&b), &ro)],
            },
            &o,
        )
        .unwrap();
        prop_assert_eq!((got.intersection, got.union, got.n_intersections), want);
    }
}

/// The merged stream holds only open clusters, never the input: an early `-`
/// and `.` record must not hold back the many finished `+` clusters after them.
#[test]
fn merged_reader_memory_is_flat() {
    use gloci_core::ChromTable;
    use gloci_io::{RecordBatch, RecordReader};
    use gloci_ops::stats::MergedReader;
    let mut s = String::from("chr1\t0\t10\tm\t0\t-\nchr1\t5\t12\tu\t0\t.\n");
    for i in 0..20_000u64 {
        s += &format!("chr1\t{}\t{}\tp{i}\t0\t+\n", 100 + i * 20, 110 + i * 20);
    }
    let mut ro = reader_opts();
    ro.format = Some(Format::Bed);
    let inner = mem_reader("a.bed", &s, &ro);
    let mut r = MergedReader::new(inner, true);
    let mut chroms = ChromTable::karyotype();
    let mut batch = RecordBatch::default();
    let mut rows = 0;
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 64).unwrap();
        rows += batch.len();
        if !more {
            break;
        }
    }
    assert_eq!(rows, 20_002);
    assert!(r.peak_held() < 16, "held {}", r.peak_held());
}
