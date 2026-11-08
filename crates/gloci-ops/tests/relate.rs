use gloci_core::{
    ChromId, ChromOrder, Fraction, Interval, OverlapCriteria, Strand, StrandRule, VcfSpan, Window,
};
use gloci_engine::{ARec, Hit, Neighbours, Visitor};
use gloci_io::{Format, InvalidPolicy, IoError, ReaderOptions, RecordMeta};
use gloci_ops::{
    Acceptor, RelateInputs, RelateOptions, b_placeholder, rewrite_coords, run_relation,
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

/// `accepted` through an Acceptor over in-memory BED (or `format`) inputs.
fn acc_as(format: Format, a: &ARec, h: Hit, o: &RelateOptions) -> Result<bool, IoError> {
    let r = ReaderOptions {
        format: Some(format),
        ..ro()
    };
    let text = if format == Format::Bed {
        "chr1\t0\t1\n"
    } else {
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\nchr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0\t0\t0\n"
    };
    let inputs = RelateInputs {
        a: mem_reader("a", text, &r),
        b: vec![mem_reader("b", text, &r)],
    };
    Acceptor::new(o, &inputs).accepted(a, &h)
}

fn acc(a: &ARec, h: Hit, o: &RelateOptions) -> bool {
    acc_as(Format::Bed, a, h, o).unwrap()
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
    assert!(!acc(&a, hit(iv(6, 20), Strand::Unknown, b""), &o));
    assert!(acc(&a, hit(iv(5, 20), Strand::Unknown, b""), &o));
}

#[test]
fn accepted_strand_rule() {
    let mut o = opts();
    o.criteria.strand = StrandRule::Same;
    let a = a_rec(iv(0, 10), Strand::Plus, b"chr1\t0\t10");
    assert!(acc(&a, hit(iv(5, 20), Strand::Plus, b""), &o));
    assert!(!acc(&a, hit(iv(5, 20), Strand::Minus, b""), &o));
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
    assert!(acc(&a, hit(iv(50, 60), Strand::Unknown, b""), &o));
}

const A12: &[u8] = b"chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,";
const B12_MISS: &[u8] = b"chr1\t20\t80\tb\t0\t+\t20\t80\t0\t2\t10,10,\t0,50,";
const B12_HIT: &[u8] = b"chr1\t5\t100\tb\t0\t+\t5\t100\t0\t2\t10,10,\t0,85,";

#[test]
fn accepted_split_uses_blocks() {
    let mut o = opts();
    o.split = true;
    let a = a_rec(iv(0, 100), Strand::Plus, A12);
    assert!(!acc(&a, hit(iv(20, 80), Strand::Plus, B12_MISS), &o));
    assert!(acc(&a, hit(iv(5, 100), Strand::Plus, B12_HIT), &o));
    // Without --split the spans overlap.
    o.split = false;
    assert!(acc(&a, hit(iv(20, 80), Strand::Plus, B12_MISS), &o));
}

#[test]
fn accepted_split_thresholds_use_block_lengths() {
    let mut o = opts();
    o.split = true;
    o.criteria.min_frac_a = Some(Fraction::parse("0.5").unwrap());
    let a = a_rec(iv(0, 100), Strand::Plus, A12);
    // A's blocks total 20 bases; B_HIT shares 5 + 10 = 15 of 20.
    assert!(acc(&a, hit(iv(5, 100), Strand::Plus, B12_HIT), &o));
    o.criteria.min_bases = Some(16);
    assert!(!acc(&a, hit(iv(5, 100), Strand::Plus, B12_HIT), &o));
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

fn spec(path: &std::path::Path, f: gloci_io::OutFormat) -> gloci_io::OutputSpec {
    gloci_io::OutputSpec {
        path: Some(path.to_path_buf()),
        format: Some(f),
        bgzf: false,
        direct: false,
    }
}

fn oo() -> gloci_io::OutputOptions {
    gloci_io::OutputOptions {
        threads: NonZeroUsize::MIN,
        level: 6,
        write_index: false,
    }
}

fn wo() -> gloci_io::WriterOptions {
    gloci_io::WriterOptions {
        write_header: false,
        command_line: Some("gloci intersect".into()),
        version: "0",
    }
}

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n";

fn vcf_reader() -> Box<dyn gloci_io::RecordReader> {
    mem_reader(
        "a.vcf",
        VCF,
        &ReaderOptions {
            format: None,
            ..ro()
        },
    )
}

#[test]
fn open_report_enforces_tsv_extensions() {
    use gloci_io::OutFormat;
    use gloci_ops::{ReportOutput, open_report};
    let kind = ReportOutput::Tsv {
        appended: vec!["n".into()],
    };
    let dir = tempfile::tempdir().unwrap();
    let spec = spec(&dir.path().join("out.vcf"), OutFormat::Vcf);
    let a = vcf_reader();
    let e = open_report(a.as_ref(), &kind, &spec, &oo(), &wo(), false)
        .err()
        .unwrap();
    assert!(e.to_string().contains("TSV"), "{e}");
}

#[test]
fn open_report_writes_bed_with_header() {
    use gloci_io::OutFormat;
    use gloci_ops::{ReportOutput, open_report};
    let dir = tempfile::tempdir().unwrap();
    let a = mem_reader("a.bed", "track name=x\nchr1\t1\t2\tn\n", &ro());
    for (header, want) in [
        (true, "track name=x\nchr1\t1\t2\tn\n"),
        (false, "chr1\t1\t2\tn\n"),
    ] {
        let path = dir.path().join(format!("out{header}.bed"));
        let mut w = open_report(
            a.as_ref(),
            &ReportOutput::SameAsA,
            &spec(&path, OutFormat::Bed),
            &oo(),
            &wo(),
            header,
        )
        .unwrap();
        w.write_record(b"chr1\t1\t2\tn").unwrap();
        w.finish().unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), want);
    }
}

#[test]
fn open_report_writes_vcf_with_header_and_command() {
    use gloci_io::OutFormat;
    use gloci_ops::{ReportOutput, open_report};
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.vcf");
    let a = vcf_reader();
    let mut w = open_report(
        a.as_ref(),
        &ReportOutput::SameAsA,
        &spec(&path, OutFormat::Vcf),
        &oo(),
        &wo(),
        false,
    )
    .unwrap();
    w.write_record(b"chr1\t5\t.\tA\tC\t.\t.\t.").unwrap();
    w.finish().unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("##glociCommand"), "{text}");
    assert!(text.contains("#CHROM\tPOS"), "{text}");
    assert!(text.ends_with("chr1\t5\t.\tA\tC\t.\t.\t.\n"), "{text}");
}

#[test]
fn open_report_writes_tsv() {
    use gloci_io::OutFormat;
    use gloci_ops::{ReportOutput, open_report};
    let dir = tempfile::tempdir().unwrap();
    let kind = ReportOutput::Tsv {
        appended: vec!["n".into()],
    };
    // BED A: no header line unless asked.
    let a = mem_reader("a.bed", "track name=x\nchr1\t1\t2\n", &ro());
    let path = dir.path().join("out.tsv");
    let mut w = open_report(
        a.as_ref(),
        &kind,
        &spec(&path, OutFormat::Tsv),
        &oo(),
        &wo(),
        false,
    )
    .unwrap();
    w.write_record(b"chr1\t1\t2\t7").unwrap();
    w.finish().unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "chr1\t1\t2\t7\n");
    // VCF A: the #CHROM line gets the appended names.
    let a = vcf_reader();
    let path = dir.path().join("v.tsv");
    let mut w = open_report(
        a.as_ref(),
        &kind,
        &spec(&path, OutFormat::Tsv),
        &oo(),
        &wo(),
        false,
    )
    .unwrap();
    w.write_record(b"chr1\t5\t.\tA\tC\t.\t.\t.\t7").unwrap();
    w.finish().unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(
        text.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tn\n"),
        "{text}"
    );
}

#[test]
fn split_ignores_vcf_columns_and_rejects_bad_bed12() {
    let mut o = opts();
    o.split = true;
    let a = a_rec(iv(0, 100), Strand::Plus, A12);
    let h = hit(iv(20, 80), Strand::Plus, B12_MISS);
    // As BED the blocks miss; as VCF the 12 columns are only a span.
    assert!(!acc_as(Format::Bed, &a, h, &o).unwrap());
    assert!(acc_as(Format::Vcf, &a, h, &o).unwrap());
    // A malformed BED12 line is a data error naming file and line.
    let bad = hit(
        iv(20, 80),
        Strand::Plus,
        b"chr1\t20\t80\tb\t0\t+\t20\t80\t0\t2\t10,10,\t0,10,",
    );
    let e = acc_as(Format::Bed, &a, bad, &o).unwrap_err();
    assert!(matches!(e, IoError::Parse { .. }), "{e}");
    let msg = e.to_string();
    assert!(msg.contains("b line 1") && msg.contains("BED12"), "{msg}");
    // Without --split the line is not parsed.
    o.split = false;
    assert!(acc_as(Format::Bed, &a, bad, &o).is_ok());
}

#[test]
fn split_with_window_pads_b_blocks() {
    let a = a_rec(iv(0, 100), Strand::Plus, A12);
    let h = hit(iv(20, 80), Strand::Plus, B12_MISS);
    let mut o = opts();
    o.split = true;
    let w = |n| {
        Some(Window {
            left: n,
            right: n,
            strand_aware: false,
        })
    };
    o.window = w(5);
    assert!(!acc(&a, h, &o));
    o.window = w(15);
    assert!(acc(&a, h, &o));
}

#[test]
fn window_with_threshold_is_a_usage_error() {
    let mut o = opts();
    o.window = Some(Window {
        left: 1,
        right: 1,
        strand_aware: false,
    });
    o.criteria.min_bases = Some(1);
    let inputs = RelateInputs {
        a: mem_reader("a.bed", "chr1\t10\t20\n", &ro()),
        b: vec![mem_reader("b.bed", "chr1\t15\t30\n", &ro())],
    };
    let e = run_relation(inputs, &o, None, &mut Collect(Vec::new())).unwrap_err();
    assert!(matches!(e, IoError::Usage(_)), "{e}");
}

struct Capture;
static LOGS: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());

impl log::Log for Capture {
    fn enabled(&self, _: &log::Metadata) -> bool {
        true
    }
    fn log(&self, r: &log::Record) {
        LOGS.lock()
            .unwrap()
            .push(format!("{}: {}", r.level(), r.args()));
    }
    fn flush(&self) {}
}

#[test]
fn unknown_strand_warns_once_per_input() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        log::set_logger(&Capture).unwrap();
        log::set_max_level(log::LevelFilter::Warn);
    });
    let mut o = opts();
    o.criteria.strand = StrandRule::Same;
    let inputs = RelateInputs {
        a: mem_reader("warn_a.bed", "chr1\t10\t20\nchr1\t12\t22\n", &ro()),
        b: vec![
            mem_reader("warn_b1.bed", "chr1\t15\t30\nchr1\t16\t30\n", &ro()),
            mem_reader("warn_b2.bed", "chr1\t15\t30\t.\t0\t+\n", &ro()),
        ],
    };
    run_relation(inputs, &o, None, &mut Collect(Vec::new())).unwrap();
    let logs = LOGS.lock().unwrap();
    let n = |name: &str| {
        logs.iter()
            .filter(|l| l.contains(name) && l.contains("unknown strand"))
            .count()
    };
    assert_eq!(
        (n("warn_a.bed"), n("warn_b1.bed"), n("warn_b2.bed")),
        (1, 1, 0),
        "{logs:?}"
    );
}
