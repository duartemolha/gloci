//! The VCF-style records of `vcf_relation_inputs`: their spans are the
//! documented ones, and the generator really makes spans that start after
//! the key (the case BED inputs never reach).

use gloci_core::Interval;
use gloci_testkit::strategies::{
    Feature, TestFormat, VcfShape, vcf_rec, vcf_relation_inputs, vcf_with_span,
};
use proptest::strategy::{Strategy, ValueTree};
use proptest::test_runner::TestRunner;

fn iv(start: u64, end: u64) -> Interval {
    Interval { start, end }
}

#[test]
fn shapes_have_the_documented_spans() {
    // The trimmed deletion `chr1 101 . AAATTT AAA`: key 100, span [103,106).
    let del = vcf_rec("chr1", 100, VcfShape::Del { prefix: 3, len: 3 }, ".".into());
    assert_eq!(del.line(), "chr1\t101\t.\tAAATTT\tAAA\t.\t.\t.");
    assert_eq!((del.key(), del.span()), (100, iv(103, 106)));

    let ins = vcf_rec("chr1", 100, VcfShape::Ins { len: 2 }, "i".into());
    assert_eq!(ins.line(), "chr1\t101\ti\tA\tACC\t.\t.\t.");
    assert_eq!((ins.key(), ins.span()), (100, iv(101, 101)));

    let snv = vcf_rec("chr1", 100, VcfShape::Snv, "s".into());
    assert_eq!((snv.key(), snv.span()), (100, iv(100, 101)));

    let sv = vcf_rec("chr1", 100, VcfShape::Sv { len: 10 }, "d".into());
    assert_eq!(sv.line(), "chr1\t101\td\tN\t<DEL>\t.\t.\tEND=111");
    assert_eq!((sv.key(), sv.span()), (100, iv(101, 111)));

    // `vcf_with_span` hits the span it is asked for, in every shape.
    for (start, end) in [(5, 5), (5, 6), (0, 1), (5, 9), (5, 70), (1, 3)] {
        for sel in 0..12 {
            let r = vcf_with_span("chr2", start, end, sel, "x".into()).unwrap();
            assert_eq!(r.span(), iv(start, end), "{r:?}");
        }
    }
    // Spans starting at 0 other than one base have no VCF record.
    assert_eq!(vcf_with_span("chr1", 0, 0, 0, "x".into()), None);
    assert_eq!(vcf_with_span("chr1", 0, 5, 1, "x".into()), None);
}

#[test]
fn generator_makes_spans_that_start_after_the_key() {
    let mut runner = TestRunner::deterministic();
    let strategy = vcf_relation_inputs();
    let (mut a_vcf, mut vcf_recs, mut after_key, mut b_vcf) = (0, 0, 0, 0);
    for _ in 0..200 {
        let input = strategy.new_tree(&mut runner).unwrap().current();
        a_vcf += usize::from(input.a.format == TestFormat::Vcf);
        b_vcf += input
            .b
            .iter()
            .filter(|f| f.format == TestFormat::Vcf)
            .count();
        for f in std::iter::once(&input.a).chain(&input.b) {
            if f.format == TestFormat::Vcf {
                vcf_recs += f.recs.len();
                after_key += f.recs.iter().filter(|r| r.span().start > r.key()).count();
            }
        }
    }
    assert!(a_vcf > 80, "A was VCF in {a_vcf} of 200 inputs");
    assert!(b_vcf > 80, "{b_vcf} VCF B files in 200 inputs");
    assert!(
        after_key * 2 > vcf_recs,
        "{after_key} of {vcf_recs} VCF spans start after the key"
    );
}
