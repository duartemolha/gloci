use gloci_core::{ChromTable, Interval, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures;
use std::io::Cursor;
use std::num::NonZeroUsize;

type Rec = (String, u64, Interval, u64, Vec<u8>);

fn opts(format: Option<Format>, invalid: InvalidPolicy) -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::new(1).unwrap(),
        format,
        vcf_span: VcfSpan::Allele,
        invalid,
    }
}
fn bed() -> ReaderOptions {
    opts(Some(Format::Bed), InvalidPolicy::Fail)
}
fn gff() -> ReaderOptions {
    opts(Some(Format::Gff3), InvalidPolicy::Fail)
}
fn skip() -> ReaderOptions {
    opts(Some(Format::Bed), InvalidPolicy::Skip)
}
fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

fn run(
    bytes: &[u8],
    o: &ReaderOptions,
    max: usize,
) -> Result<(Header, Vec<Rec>, Diagnostics), IoError> {
    let mut r = open_reader_from("x.bed", Box::new(Cursor::new(bytes.to_vec())), o)?;
    let mut chroms = ChromTable::lexicographic();
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, max)?;
        assert!(batch.len() <= max);
        for (i, m) in batch.recs.iter().enumerate() {
            out.push((
                String::from_utf8(chroms.name(m.chrom).to_vec()).unwrap(),
                m.key,
                m.span,
                m.loc,
                batch.raw(i).to_vec(),
            ));
        }
        if !more {
            break;
        }
    }
    Ok((r.header().clone(), out, r.diagnostics().clone()))
}

fn read_all(bytes: &[u8], o: ReaderOptions) -> (Header, Vec<Rec>, Diagnostics) {
    run(bytes, &o, 1000).unwrap()
}
fn read_all_max(bytes: &[u8], max: usize) -> Vec<Rec> {
    run(bytes, &bed(), max).unwrap().1
}
fn err(bytes: &[u8]) -> String {
    run(bytes, &bed(), 1000).unwrap_err().to_string()
}
fn err_gff(bytes: &[u8]) -> String {
    run(bytes, &gff(), 1000).unwrap_err().to_string()
}

#[test]
fn bed_headers_blanks_crlf_and_locs() {
    let (h, recs, _) = read_all(
        b"#c\ntrack name=x\nchr1\t5\t9\tn\r\n\nchr1\t7\t7\n#mid\nchr2\t0\t3\n",
        bed(),
    );
    assert_eq!(h, Header::Lines(b"#c\ntrack name=x\n".to_vec()));
    assert_eq!(
        recs,
        vec![
            ("chr1".into(), 5, iv(5, 9), 3, b"chr1\t5\t9\tn".to_vec()),
            (
                "chr1".into(),
                7,
                Interval::point(7),
                5,
                b"chr1\t7\t7".to_vec()
            ),
            ("chr2".into(), 0, iv(0, 3), 7, b"chr2\t0\t3".to_vec())
        ]
    );
}

#[test]
fn bed_errors() {
    assert_eq!(
        err(b"chr1\t1\t2\nchr1\t10\t5\n"),
        "x.bed line 2: start 10 is greater than end 5"
    );
    assert_eq!(
        err(b"chr1\t-1\t5\n"),
        r#"x.bed line 1: start "-1" is negative"#
    );
    assert_eq!(
        err(b"chr1\t1\tx\n"),
        r#"x.bed line 1: end "x" is not an integer"#
    );
    assert!(err(b"chr1 1 5\n").ends_with("(the line contains spaces; BED must be tab-separated)"));
    assert!(err(b"chr1 1 5\n").contains("expected at least 3 tab-separated columns, found 1"));
}

#[test]
fn skip_policy_counts_and_keeps_examples() {
    let (_, recs, d) = read_all(b"chr1\t1\t2\nbad\nchr1\t3\t1\nchr1\t4\t5\n", skip());
    assert_eq!((recs.len(), d.skipped), (2, 2));
    assert!(d.examples[0].starts_with("line 2:"));
    assert!(
        d.summary("x.bed")
            .unwrap()
            .starts_with("skipped 2 invalid records in x.bed; first 2:\n  line 2:")
    );
    let many: String = (0..25).map(|_| "bad\n").collect();
    let (_, recs, d) = read_all(many.as_bytes(), skip());
    assert_eq!((recs.len(), d.skipped, d.examples.len()), (0, 25, 10));
    assert!(Diagnostics::default().summary("x").is_none());
}

#[test]
fn batch_size_does_not_change_results() {
    let input: String = (0..50)
        .map(|i| format!("chr1\t{i}\t{}\tn{i}\n", i + 3))
        .collect();
    let a = read_all_max(input.as_bytes(), 1);
    assert_eq!(a.len(), 50);
    assert_eq!(a, read_all_max(input.as_bytes(), 1000));
    assert_eq!(a, read_all_max(input.as_bytes(), 7));
}

#[test]
fn closed_table_rejects_unknown_chrom() {
    let g = gloci_core::Genome::parse_sizes(b"chr1\t100\n").unwrap();
    let mut chroms = ChromTable::closed(&g, gloci_core::OrderSource::GenomeFile("g.genome".into()));
    let mut r = open_reader_from(
        "x.bed",
        Box::new(Cursor::new(b"chrZ\t1\t2\n".to_vec())),
        &bed(),
    )
    .unwrap();
    let mut b = RecordBatch::default();
    let e = r
        .read_batch(&mut chroms, &mut b, 10)
        .unwrap_err()
        .to_string();
    assert_eq!(
        e,
        r#"x.bed line 1: chromosome "chrZ" is not in the genome file g.genome"#
    );
}

#[test]
fn gff_reader() {
    let (h, recs, _) = read_all(
        b"##gff-version 3\nchr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g\n###\n",
        gff(),
    );
    assert_eq!(h, Header::Lines(b"##gff-version 3\n".to_vec()));
    assert_eq!((recs[0].1, recs[0].2), (99, iv(99, 200)));
    assert_eq!(recs.len(), 1);
    assert!(
        err_gff(b"##gff-version 3\n##FASTA\n>x\nACGT\n")
            .contains("##FASTA sections are not supported yet")
    );
}

#[test]
fn bed12() {
    assert_eq!(
        bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t2\t10,20,\t0,80,")
            .unwrap()
            .unwrap(),
        vec![iv(100, 110), iv(180, 200)]
    );
    assert!(bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t2\t10,20,\t0,70,").is_err());
    assert_eq!(bed12_blocks(b"c\t1\t2\tn").unwrap(), None);
    assert!(bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t3\t10,20,\t0,80,").is_err());
    assert!(bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t2\t10,20,\t5,80,").is_err());
    assert!(bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t2\t50,20,\t0,40,").is_err());
}

#[test]
fn long_lines_across_bgzf_blocks() {
    let long = (0..20)
        .map(|i| format!("chr1\t{i}\t{}\t{}\n", i + 1, "x".repeat(40_000)))
        .collect::<String>();
    let (_, recs, _) = read_all(&fixtures::bgzip(long.as_bytes()), bed());
    assert!(recs.iter().all(|r| r.4.len() > 40_000) && recs.len() == 20);
}

#[test]
fn bedgraph_and_vcf_dispatch() {
    let r = open_reader_from(
        "x",
        Box::new(Cursor::new(
            b"track type=bedGraph\nchr1\t1\t2\t0.5\n".to_vec(),
        )),
        &opts(None, InvalidPolicy::Fail),
    )
    .unwrap();
    assert_eq!(r.format(), Format::BedGraph);
    let v = open_reader_from(
        "x.vcf",
        Box::new(Cursor::new(
            b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec(),
        )),
        &opts(None, InvalidPolicy::Fail),
    )
    .unwrap();
    assert_eq!(v.format(), Format::Vcf);
    assert!(matches!(v.header(), Header::Vcf(_)));
}

#[test]
fn final_line_cr_without_lf_strips_cr() {
    let bytes = b"chr1\t10\t20\r";
    let (_, recs, _) = read_all(bytes, bed());
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].1, 10);
    assert_eq!(recs[0].2, iv(10, 20));
}

#[test]
#[should_panic(expected = "read_batch max must be > 0")]
fn read_batch_zero_panics() {
    let mut r = open_reader_from(
        "x.bed",
        Box::new(Cursor::new(b"chr1\t10\t20\n".to_vec())),
        &bed(),
    )
    .unwrap();
    let mut chroms = ChromTable::lexicographic();
    let mut batch = RecordBatch::default();
    let _ = r.read_batch(&mut chroms, &mut batch, 0);
}
