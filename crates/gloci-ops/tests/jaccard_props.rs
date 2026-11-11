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
