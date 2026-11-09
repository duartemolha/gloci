use gloci_core::{ChromTable, Interval, VcfSpan};
use gloci_io::*;
use std::io::Cursor;
use std::num::NonZeroUsize;

type Rec = (String, u64, Interval, u64);

const H: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr2,length=500>\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

fn vcf(span: VcfSpan) -> ReaderOptions {
    ReaderOptions {
        chrom_map: None,
        threads: NonZeroUsize::new(1).unwrap(),
        format: Some(Format::Vcf),
        vcf_span: span,
        invalid: InvalidPolicy::Fail,
    }
}
fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}
fn run(bytes: &[u8], o: ReaderOptions) -> Result<(Header, Vec<Rec>, Diagnostics), IoError> {
    let mut r = open_reader_from("x.vcf", Box::new(Cursor::new(bytes.to_vec())), &o)?;
    let mut chroms = ChromTable::lexicographic();
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 1000)?;
        for m in &batch.recs {
            let name = String::from_utf8(chroms.name(m.chrom).to_vec()).unwrap();
            out.push((name, m.key, m.span, m.loc));
        }
        if !more {
            break;
        }
    }
    let d = r.diagnostics().clone();
    Ok((r.header().clone(), out, d))
}
fn read_all(bytes: impl AsRef<[u8]>, o: ReaderOptions) -> (Header, Vec<Rec>, Diagnostics) {
    run(bytes.as_ref(), o).unwrap()
}
fn err(bytes: impl AsRef<[u8]>) -> String {
    run(bytes.as_ref(), vcf(VcfSpan::Allele))
        .unwrap_err()
        .to_string()
}

#[test]
fn chromx_is_not_the_chrom_line() {
    let bad_header = "##fileformat=VCFv4.3\n#CHROMX\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let e = err(bad_header.as_bytes());
    assert!(
        e.contains("no #CHROM line"),
        "expected no #CHROM line error, got: {e}"
    );
}

#[test]
fn header_records_spans_and_provenance() {
    let same_input = format!("{H}chr1\t100\t.\tACGT\tA\t.\tPASS\t.\tGT\t0/1\t0/0\n");
    let (h, recs, _) = read_all(same_input.as_bytes(), vcf(VcfSpan::Allele));
    let Header::Vcf(h) = h else { panic!() };
    assert_eq!(
        h.contigs,
        vec![(b"chr2".to_vec(), Some(500)), (b"chr1".to_vec(), None)]
    );
    assert_eq!(h.samples, vec![b"s1".to_vec(), b"s2".to_vec()]);
    assert_eq!((recs[0].1, recs[0].2, recs[0].3), (99, iv(100, 103), 5));
    assert_eq!(read_all(&same_input, vcf(VcfSpan::Ref)).1[0].2, iv(99, 103));
    let out = String::from_utf8(h.with_command_line("gloci sort in.vcf", "0.1.0")).unwrap();
    assert!(
        out.contains("##contig=<ID=chr1>\n##glociCommand=gloci sort in.vcf; Version=0.1.0\n#CHROM")
    );
    assert_eq!(
        err(format!("{H}chr1\t12a4\t.\tA\tG\t.\t.\t.\n")),
        r#"x.vcf line 5: POS "12a4" is not an integer"#
    );
    assert!(
        err(format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\t.\n"))
            .contains("needs INFO/END or INFO/SVLEN")
    );
    assert_eq!(
        read_all(
            format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=500\n"),
            vcf(VcfSpan::Allele)
        )
        .1[0]
            .2,
        iv(100, 500)
    );
    assert!(err(format!("{H}chr1\t0\t.\tA\tG\t.\t.\t.\n")).contains("POS 0"));
    assert_eq!(
        shell_quote_args(&["gloci".into(), "view".into(), "a b.vcf".into()]),
        "gloci view 'a b.vcf'"
    );
    // header-only VCF
    let (h, recs, _) = read_all(H.as_bytes(), vcf(VcfSpan::Allele));
    assert!(recs.is_empty() && matches!(h, Header::Vcf(_)));
}

#[test]
fn more_errors_and_edge_cases() {
    assert!(
        err(format!("{H}chr1\t100\t.\tA\tG\n"))
            .contains("expected at least 8 tab-separated columns, found 5")
    );
    assert_eq!(
        err("chr1\t1\t.\tA\tG\t.\t.\t.\n"),
        "x.vcf: VCF header has no #CHROM line"
    );
    // HOMEND= must not be read as END=
    assert!(
        err(format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\tHOMEND=900\n"))
            .contains("needs INFO/END or INFO/SVLEN")
    );
    assert!(err(format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\tEND=x\n")).contains("END"));
    // SVLEN, END only scanned for symbolic ALT, CRLF normalised
    let r = read_all(
        format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\tSVLEN=-50\n").replace('\n', "\r\n"),
        vcf(VcfSpan::Allele),
    );
    assert_eq!(r.1[0].2, iv(100, 150));
    let Header::Vcf(h) = r.0 else { panic!() };
    assert_eq!(h.raw, H.as_bytes());
    let bad = "##fileformat=VCFv4.3\n##contig=<ID=c,length=abc>\n#CHROM\tPOS\n".to_string();
    let e = err(bad);
    assert!(
        e.starts_with("x.vcf: header line ") && e.contains("length=abc"),
        "{e}"
    );
}

#[test]
fn quoting() {
    let q = |a: &[&str]| shell_quote_args(&a.iter().map(|s| (*s).into()).collect::<Vec<_>>());
    assert_eq!(q(&["it's"]), r"'it'\''s'");
    assert_eq!(q(&["a\tb\nc"]), r"'a\tb\nc'");
    assert_eq!(
        q(&["plain-1.vcf", "x;y", "$HOME"]),
        "plain-1.vcf 'x;y' '$HOME'"
    );
}
