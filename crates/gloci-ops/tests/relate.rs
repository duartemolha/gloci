use gloci_core::{
    ChromId, ChromOrder, Fraction, Interval, OverlapCriteria, Strand, StrandRule, VcfSpan, Window,
};
use gloci_engine::{ARec, Hit, Neighbours, Visitor};
use gloci_io::{Format, InvalidPolicy, IoError, ReaderOptions, RecordMeta};
use gloci_ops::{
    RelateInputs, RelateOptions, accepted, b_placeholder, rewrite_coords, run_relation,
};
use gloci_testkit::mem::mem_reader;
use std::num::NonZeroUsize;

fn opts() -> RelateOptions {
    RelateOptions {
        criteria: OverlapCriteria::default(),
        window: None,
        split: false,
        unsorted: false,
        labels: vec![],
        genome: None,
        chrom_order: ChromOrder::Karyotype,
        header: false,
    }
}

fn a_rec<'a>(span: Interval, strand: Strand, raw: &'a [u8]) -> ARec<'a> {
    ARec {
        meta: RecordMeta {
            chrom: ChromId(0),
            key: span.start,
            span,
            strand,
            loc: 1,
            off: 0,
            len: raw.len() as u32,
        },
        raw,
    }
}

fn hit<'a>(span: Interval, strand: Strand, raw: &'a [u8]) -> Hit<'a> {
    Hit {
        file: 0,
        span,
        strand,
        raw,
        loc: 1,
    }
}

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

#[test]
fn rewrite_bed() {
    let mut out = Vec::new();
    rewrite_coords(Format::Bed, b"chr1\t5\t15\tx", iv(7, 9), &mut out);
    assert_eq!(out, b"chr1\t7\t9\tx");
}

#[test]
fn rewrite_gff_is_one_based_closed() {
    let mut out = Vec::new();
    let raw = b"chr1\tsrc\tgene\t6\t15\t.\t+\t.\tID=g";
    rewrite_coords(Format::Gff3, raw, iv(6, 9), &mut out);
    assert_eq!(out, b"chr1\tsrc\tgene\t7\t9\t.\t+\t.\tID=g");
}

#[test]
fn placeholders() {
    let p = |f, n| {
        let mut o = Vec::new();
        b_placeholder(f, n, &mut o);
        String::from_utf8(o).unwrap()
    };
    assert_eq!(p(Format::Bed, 3), ".\t-1\t-1");
    assert_eq!(p(Format::Bed, 5), ".\t-1\t-1\t.\t.");
    assert_eq!(p(Format::Vcf, 8), ".\t.\t.\t.\t.\t.\t.\t.");
    assert_eq!(p(Format::Gff3, 9).split('\t').count(), 9);
}

#[test]
fn accepted_fraction_uses_unpadded_b() {
    let mut o = opts();
    o.criteria.min_frac_a = Some(Fraction::parse("0.5").unwrap());
    let a = a_rec(iv(0, 10), Strand::Unknown, b"chr1\t0\t10");
    // 4 of 10 bases of A: refused; 5 of 10: accepted.
    assert!(!accepted(&a, &hit(iv(6, 20), Strand::Unknown, b""), &o));
    assert!(accepted(&a, &hit(iv(5, 20), Strand::Unknown, b""), &o));
}

#[test]
fn accepted_strand_rule() {
    let mut o = opts();
    o.criteria.strand = StrandRule::Same;
    let a = a_rec(iv(0, 10), Strand::Plus, b"chr1\t0\t10");
    assert!(accepted(&a, &hit(iv(5, 20), Strand::Plus, b""), &o));
    assert!(!accepted(&a, &hit(iv(5, 20), Strand::Minus, b""), &o));
}

#[test]
fn accepted_with_window_trusts_the_candidate() {
    let mut o = opts();
    o.window = Some(Window {
        left: 100,
        right: 100,
        strand_aware: false,
    });
    let a = a_rec(iv(0, 10), Strand::Unknown, b"chr1\t0\t10");
    assert!(accepted(&a, &hit(iv(50, 60), Strand::Unknown, b""), &o));
}

const A12: &[u8] = b"chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,";
const B12_MISS: &[u8] = b"chr1\t20\t80\tb\t0\t+\t20\t80\t0\t2\t10,10,\t0,50,";
const B12_HIT: &[u8] = b"chr1\t5\t100\tb\t0\t+\t5\t100\t0\t2\t10,10,\t0,85,";

#[test]
fn accepted_split_uses_blocks() {
    let mut o = opts();
    o.split = true;
    let a = a_rec(iv(0, 100), Strand::Plus, A12);
    assert!(!accepted(&a, &hit(iv(20, 80), Strand::Plus, B12_MISS), &o));
    assert!(accepted(&a, &hit(iv(5, 100), Strand::Plus, B12_HIT), &o));
    // Without --split the spans overlap.
    o.split = false;
    assert!(accepted(&a, &hit(iv(20, 80), Strand::Plus, B12_MISS), &o));
}

#[test]
fn accepted_split_thresholds_use_block_lengths() {
    let mut o = opts();
    o.split = true;
    o.criteria.min_frac_a = Some(Fraction::parse("0.5").unwrap());
    let a = a_rec(iv(0, 100), Strand::Plus, A12);
    // A's blocks total 20 bases; B_HIT shares 5 + 10 = 15 of 20.
    assert!(accepted(&a, &hit(iv(5, 100), Strand::Plus, B12_HIT), &o));
    o.criteria.min_bases = Some(16);
    assert!(!accepted(&a, &hit(iv(5, 100), Strand::Plus, B12_HIT), &o));
}

struct Collect(Vec<(u64, Vec<u64>)>);

impl Visitor for Collect {
    fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], _: &Neighbours<'_>) -> Result<(), IoError> {
        self.0.push((
            a.meta.span.start,
            hits.iter().map(|h| h.span.start).collect(),
        ));
        Ok(())
    }
}

fn ro() -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

#[test]
fn relation_from_memory() {
    for unsorted in [false, true] {
        let mut o = opts();
        o.unsorted = unsorted;
        let inputs = RelateInputs {
            a: mem_reader("a.bed", "chr1\t10\t20\nchr2\t5\t9\n", &ro()),
            b: vec![mem_reader(
                "b.bed",
                "chr1\t15\t30\nchr1\t40\t50\nchr2\t0\t6\n",
                &ro(),
            )],
        };
        let mut v = Collect(Vec::new());
        let s = run_relation(inputs, &o, None, &mut v).unwrap();
        assert_eq!(v.0, vec![(10, vec![15]), (5, vec![0])]);
        assert_eq!(s.a_records, 2);
        assert!(s.shared_chroms);
    }
}

#[test]
fn relation_strand_rule_still_runs_with_unknown_strands() {
    let mut o = opts();
    o.criteria.strand = StrandRule::Same;
    let inputs = RelateInputs {
        a: mem_reader("a.bed", "chr1\t10\t20\n", &ro()),
        b: vec![mem_reader("b.bed", "chr1\t15\t30\n", &ro())],
    };
    let mut v = Collect(Vec::new());
    run_relation(inputs, &o, None, &mut v).unwrap();
    assert_eq!(v.0.len(), 1);
}

#[test]
fn open_report_enforces_tsv_extensions() {
    use gloci_io::{OutFormat, OutputOptions, OutputSpec, WriterOptions};
    use gloci_ops::{ReportOutput, open_report};
    let a = mem_reader(
        "a.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
        &ReaderOptions {
            format: None,
            ..ro()
        },
    );
    let spec = OutputSpec {
        path: Some("out.vcf".into()),
        format: Some(OutFormat::Vcf),
        bgzf: false,
        direct: false,
    };
    let oo = OutputOptions {
        threads: NonZeroUsize::MIN,
        level: 6,
        write_index: false,
    };
    let wo = WriterOptions {
        write_header: false,
        command_line: None,
        version: "0",
    };
    let kind = ReportOutput::Tsv {
        appended: vec!["n".into()],
    };
    let e = open_report(a.as_ref(), &kind, &spec, &oo, &wo, false)
        .err()
        .unwrap();
    assert!(e.to_string().contains("TSV"), "{e}");
}
