//! The BEDPE reader (D146): column rules, absent ends, strands, headers,
//! `--chrom-map` and `--skip-invalid`.

use gloci_core::{ChromTable, Interval, Strand};
use gloci_io::{BedpeBatch, BedpeReader, ChromMap, InvalidPolicy, ReaderOptions, bedpe_name};
use std::num::NonZeroUsize;
use std::sync::Arc;

fn opts(invalid: InvalidPolicy, map: Option<&str>) -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid,
        chrom_map: map.map(|m| Arc::new(ChromMap::parse("m.tsv", m.as_bytes()).unwrap())),
        accept: Default::default(),
    }
}

fn open(text: &str, o: &ReaderOptions) -> Result<BedpeReader, gloci_io::IoError> {
    BedpeReader::open_from(
        "x.bedpe",
        Box::new(std::io::Cursor::new(text.as_bytes().to_vec())),
        o,
    )
}

fn read(text: &str, o: &ReaderOptions) -> Result<(BedpeBatch, ChromTable, BedpeReader), String> {
    let mut r = open(text, o).map_err(|e| e.to_string())?;
    let mut chroms = ChromTable::karyotype();
    let mut b = BedpeBatch::default();
    // Small batches, so the loop crosses batch boundaries.
    let mut all = BedpeBatch::default();
    loop {
        b.clear();
        let more = r
            .read_batch(&mut chroms, &mut b, 2)
            .map_err(|e| e.to_string())?;
        for i in 0..b.len() {
            let mut m = b.recs[i];
            m.off = all.buf.len();
            all.buf.extend_from_slice(b.raw(i));
            all.recs.push(m);
        }
        if !more {
            break;
        }
    }
    Ok((all, chroms, r))
}

fn err(text: &str) -> String {
    match read(text, &opts(InvalidPolicy::Fail, None)) {
        Err(e) => e,
        Ok(_) => panic!("{text:?} should have failed"),
    }
}

fn iv(s: u64, e: u64) -> Interval {
    Interval { start: s, end: e }
}

#[test]
fn bedpe_reader_table() {
    // Six columns: both ends, unknown strands.
    let (b, c, _) = read(
        "chr1\t10\t20\tchr2\t30\t40\n",
        &opts(InvalidPolicy::Fail, None),
    )
    .unwrap();
    assert_eq!(b.len(), 1);
    let [e1, e2] = b.recs[0].ends;
    let (e1, e2) = (e1.unwrap(), e2.unwrap());
    assert_eq!(c.name(e1.chrom), b"chr1");
    assert_eq!(c.name(e2.chrom), b"chr2");
    assert_eq!((e1.span, e2.span), (iv(10, 20), iv(30, 40)));
    assert_eq!((e1.strand, e2.strand), (Strand::Unknown, Strand::Unknown));
    assert_eq!(b.recs[0].loc, 1);

    // Ten columns: strands set; the name is column 7.
    let line = "chr1\t10\t20\tchr2\t30\t40\tr1\t5\t+\t-\tx\ty";
    let (b, _, _) = read(&format!("{line}\n"), &opts(InvalidPolicy::Fail, None)).unwrap();
    let [e1, e2] = b.recs[0].ends;
    assert_eq!(e1.unwrap().strand, Strand::Plus);
    assert_eq!(e2.unwrap().strand, Strand::Minus);
    assert_eq!(bedpe_name(b.raw(0)), Some(&b"r1"[..]));
    assert_eq!(b.raw(0), line.as_bytes());

    // 7 and 8 columns are fine; no name below 7.
    let (b, _, _) = read(
        "chr1\t1\t2\tchr1\t3\t4\nchr1\t1\t2\tchr1\t3\t4\tn\nchr1\t1\t2\tchr1\t3\t4\tn\t9\n",
        &opts(InvalidPolicy::Fail, None),
    )
    .unwrap();
    assert_eq!(b.len(), 3);
    assert_eq!(bedpe_name(b.raw(0)), None);
    assert_eq!(bedpe_name(b.raw(2)), Some(&b"n"[..]));

    // An absent end: `.` with -1 coordinates; the other end is a point on +.
    let (b, c, _) = read(
        ".\t-1\t-1\tchr1\t5\t5\tr1\t0\t.\t+\n",
        &opts(InvalidPolicy::Fail, None),
    )
    .unwrap();
    let [e1, e2] = b.recs[0].ends;
    assert!(e1.is_none());
    let e2 = e2.unwrap();
    assert_eq!(
        (c.name(e2.chrom), e2.span, e2.strand),
        (&b"chr1"[..], iv(5, 5), Strand::Plus)
    );
    // Both ends absent is a valid record that overlaps nothing.
    let (b, _, _) = read(".\t-1\t-1\t.\t-1\t-1\n", &opts(InvalidPolicy::Fail, None)).unwrap();
    assert_eq!(b.recs[0].ends, [None, None]);

    // Errors.
    let e = err("chr1\t1\t2\tchr1\t3\t4\tn\t9\t+\n");
    assert!(e.contains("6, 7, 8 or at least 10 columns"), "{e}");
    assert!(e.contains("line 1"), "{e}");
    let e = err("chr1\t20\t10\tchr2\t30\t40\n");
    assert!(e.contains("start1"), "{e}");
    let e = err("chr1\t1\t2\tchr2\t40\t30\n");
    assert!(e.contains("start2"), "{e}");
    let e = err("chr1\t-5\t10\tchr2\t30\t40\n");
    assert!(e.contains("negative"), "{e}");
    let e = err("chr1\t1\t2\tchr2\t3\t4\tr\t0\tx\t+\n");
    assert!(e.contains("column 9"), "{e}");
    let e = err("chr1\t1\t2\tchr2\t3\t4\tr\t0\t+\t?\n");
    assert!(e.contains("column 10"), "{e}");
    let e = err("chr1\t1\t2\tchr2\t3\tx\n");
    assert!(e.contains("end2"), "{e}");
    let e = err(".\tx\t-1\tchr1\t5\t5\n");
    assert!(e.contains("start1"), "{e}");
    let e = err("chr1\t1\t2\tchr2\t3\t4\n\nchr1\t1\n");
    assert!(e.contains("line 3"), "{e}");
}

#[test]
fn bedpe_headers_crlf_and_blank_lines() {
    let text = "#c1\r\ntrack name=t\r\nbrowser position chr1\r\nchr1\t10\t20\tchr2\t30\t40\r\n\r\n#late\r\nchr1\t1\t2\tchr1\t3\t4\r\n";
    let (b, _, r) = read(text, &opts(InvalidPolicy::Fail, None)).unwrap();
    let h: Vec<&[u8]> = r.header().iter().map(Vec::as_slice).collect();
    assert_eq!(h, [&b"#c1"[..], b"track name=t", b"browser position chr1"]);
    assert_eq!(b.len(), 2);
    assert_eq!(b.raw(0), b"chr1\t10\t20\tchr2\t30\t40");
    assert_eq!(b.recs[1].loc, 7);
}

#[test]
fn bedpe_chrom_map_rewrites_columns_1_and_4() {
    let (b, c, _) = read(
        "1\t10\t20\t2\t30\t40\tr\t1\t+\t-\t1\n1\t1\t2\tX\t3\t4\n",
        &opts(InvalidPolicy::Fail, Some("1\tchr1\n2\tchr2\n")),
    )
    .unwrap();
    assert_eq!(b.raw(0), b"chr1\t10\t20\tchr2\t30\t40\tr\t1\t+\t-\t1");
    assert_eq!(b.raw(1), b"chr1\t1\t2\tX\t3\t4");
    let [e1, e2] = b.recs[0].ends;
    assert_eq!(c.name(e1.unwrap().chrom), b"chr1");
    assert_eq!(c.name(e2.unwrap().chrom), b"chr2");
    // Offsets stay valid across batches.
    assert_eq!(bedpe_name(b.raw(0)), Some(&b"r"[..]));
}

#[test]
fn bedpe_skip_invalid_counts() {
    let text = "chr1\t1\t2\tchr2\t3\t4\nchr1\t9\t2\tchr2\t3\t4\nbad\nchr1\t5\t6\tchr2\t7\t8\n";
    let (b, _, r) = read(text, &opts(InvalidPolicy::Skip, None)).unwrap();
    assert_eq!(b.len(), 2);
    assert_eq!(b.raw(1), b"chr1\t5\t6\tchr2\t7\t8");
    assert_eq!(r.diagnostics().skipped, 2);
}

#[test]
fn bedpe_refuses_other_formats() {
    let vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let e = open(vcf, &opts(InvalidPolicy::Fail, None))
        .err()
        .unwrap()
        .to_string();
    assert!(e.contains("VCF"), "{e}");
    let e = open(
        "##gff-version 3\nchr1\ts\tg\t1\t2\t.\t+\t.\tID=a\n",
        &opts(InvalidPolicy::Fail, None),
    )
    .err()
    .unwrap()
    .to_string();
    assert!(e.contains("GFF3"), "{e}");
}

#[test]
fn bedpe_empty_input_has_no_records() {
    let (b, _, _) = read("", &opts(InvalidPolicy::Fail, None)).unwrap();
    assert!(b.is_empty());
}
