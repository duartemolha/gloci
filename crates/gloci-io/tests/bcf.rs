use gloci_core::{ChromTable, Genome, Interval, OrderSource, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures;
use noodles_vcf::header::StringMaps;
use noodles_vcf::header::string_maps::StringMap;
use noodles_vcf::variant::RecordBuf;
use std::io::{Cursor, Read};
use std::num::NonZeroUsize;

/// chrom, key, span, loc, raw bytes
type Rec = (String, u64, Interval, u64, Vec<u8>);

/// Nine records: SNV, deletion, insertion, multi-allelic, <DEL> with END,
/// <DEL> with SVLEN only, <*> gVCF block, breakend, ALT "." — on contigs
/// declared with IDX= out of order, with two samples. Descriptions hold
/// commas and a fake `IDX=` inside quotes.
const FIXTURE_VCF: &[u8] = b"##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth, not IDX=9\">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length difference, ALT minus REF\">\n\
##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n\
##ALT=<ID=DEL,Description=\"Deletion\">\n\
##contig=<ID=chr2,length=100000,IDX=1>\n\
##contig=<ID=chr1,length=100000,IDX=0>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
chr2\t100\trs1\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:5\t1/1:7\n\
chr2\t200\t.\tACGT\tA\t30\tq10\tDP=3\tGT:DP\t0/1:2\t0/0:1\n\
chr2\t300\t.\tA\tACGT\t.\t.\t.\tGT\t0/1\t./.\n\
chr1\t400\t.\tACGT\tA,AC\t.\tPASS\t.\tGT\t1/2\t0/1\n\
chr1\t500\t.\tN\t<DEL>\t.\tPASS\tDP=4;SVTYPE=DEL;END=40000\tGT\t0/1\t0/0\n\
chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50\tGT\t0/1\t0/0\n\
chr1\t2000\t.\tA\t<*>\t.\t.\tEND=2500\tGT:DP\t0/0:20\t0/0:30\n\
chr1\t3000\tbnd1\tG\tG]chr2:5000]\t.\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\n\
chr1\t4000\t.\tA\t.\t.\t.\tDP=1\tGT\t0/0\t0/0\n";

fn opts(format: Option<Format>, span: VcfSpan, invalid: InvalidPolicy) -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        format,
        vcf_span: span,
        invalid,
    }
}
fn vcf(span: VcfSpan) -> ReaderOptions {
    opts(Some(Format::Vcf), span, InvalidPolicy::Fail)
}
/// BCF is detected from its content.
fn bcf(span: VcfSpan) -> ReaderOptions {
    opts(None, span, InvalidPolicy::Fail)
}
fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

fn run_with(
    bytes: &[u8],
    o: &ReaderOptions,
    chroms: &mut ChromTable,
    max: usize,
) -> Result<(Header, Vec<Rec>, Diagnostics), IoError> {
    let mut r = open_reader_from("x.bcf", Box::new(Cursor::new(bytes.to_vec())), o)?;
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(chroms, &mut batch, max)?;
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
fn run(bytes: &[u8], o: &ReaderOptions) -> Result<(Header, Vec<Rec>, Diagnostics), IoError> {
    run_with(bytes, o, &mut ChromTable::lexicographic(), 1000)
}
fn read_all(bytes: impl AsRef<[u8]>, o: ReaderOptions) -> (Header, Vec<Rec>, Diagnostics) {
    run(bytes.as_ref(), &o).unwrap()
}
fn err_bcf(bytes: impl AsRef<[u8]>) -> String {
    run(bytes.as_ref(), &bcf(VcfSpan::Allele))
        .unwrap_err()
        .to_string()
}
fn strip_raw(recs: Vec<Rec>) -> Vec<(String, u64, Interval)> {
    recs.into_iter().map(|(c, k, s, _, _)| (c, k, s)).collect()
}
fn read_header(bytes: impl AsRef<[u8]>) -> VcfHeader {
    match read_all(bytes, bcf(VcfSpan::Allele)).0 {
        Header::Vcf(h) => h,
        other => panic!("expected a VCF header, got {other:?}"),
    }
}

fn header_text(vcf: &[u8]) -> &str {
    let s = std::str::from_utf8(vcf).unwrap();
    let end = s.find("\n#CHROM").unwrap();
    let end = end + 1 + s[end + 1..].find('\n').unwrap() + 1;
    &s[..end]
}
fn noodles_maps(vcf: &[u8]) -> StringMaps {
    header_text(vcf).parse().unwrap()
}
fn names(m: &StringMap) -> Vec<Vec<u8>> {
    (0..)
        .map_while(|i| m.get_index(i))
        .map(|s| s.as_bytes().to_vec())
        .collect()
}
fn noodles_contig_names(vcf: &[u8]) -> Vec<Vec<u8>> {
    names(noodles_maps(vcf).contigs())
}
fn noodles_string_names(vcf: &[u8]) -> Vec<Vec<u8>> {
    names(noodles_maps(vcf).strings())
}
fn noodles_records(bcf_bytes: &[u8]) -> Vec<RecordBuf> {
    let mut r = noodles_bcf::io::Reader::new(bcf_bytes);
    let header = r.read_header().unwrap();
    r.records()
        .map(|rec| RecordBuf::try_from_variant_record(&header, &rec.unwrap()).unwrap())
        .collect()
}
fn bgzf_decompress(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(bytes)
        .read_to_end(&mut out)
        .unwrap();
    out
}
/// Header and records copied through gloci's raw BCF framing.
fn copy_via_raw(bcf_bytes: &[u8]) -> Vec<u8> {
    let (h, recs, _) = read_all(bcf_bytes, bcf(VcfSpan::Allele));
    let Header::Vcf(h) = h else { panic!() };
    let mut out = Vec::new();
    write_bcf_header(&mut out, &h.raw).unwrap();
    for r in &recs {
        out.extend_from_slice(&r.4);
    }
    fixtures::bgzip(&out)
}
/// Offset of the first record in decompressed BCF.
fn first_record(plain: &[u8]) -> usize {
    9 + u32::from_le_bytes(plain[5..9].try_into().unwrap()) as usize
}

#[test]
fn spans_and_keys_match_the_vcf_reader() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    for span in [VcfSpan::Allele, VcfSpan::Ref] {
        let from_vcf = read_all(FIXTURE_VCF, vcf(span)).1;
        let from_bcf = read_all(&bcf_bytes, bcf(span)).1;
        assert_eq!(from_bcf.len(), 9);
        assert_eq!(
            from_bcf.iter().map(|r| r.3).collect::<Vec<_>>(),
            (1..=9).collect::<Vec<u64>>(),
            "loc is the 1-based record number"
        );
        assert_eq!(strip_raw(from_vcf), strip_raw(from_bcf), "{span:?}");
    }
    // spot checks so a shared mistake in both readers can't hide
    let recs = strip_raw(read_all(&bcf_bytes, bcf(VcfSpan::Allele)).1);
    assert_eq!(recs[0], ("chr2".into(), 99, iv(99, 100)));
    assert_eq!(recs[2], ("chr2".into(), 299, iv(300, 300)));
    assert_eq!(recs[4], ("chr1".into(), 499, iv(500, 40000)));
    assert_eq!(recs[5], ("chr1".into(), 999, iv(1000, 1050)));
    assert_eq!(recs[6], ("chr1".into(), 1999, iv(1999, 2500)));
    assert_eq!(recs[8], ("chr1".into(), 3999, iv(3999, 4000)));
    // explicit --format bcf, uncompressed BCF (htslib -Ou), and BCF 2.1 magic
    let forced = opts(Some(Format::Bcf), VcfSpan::Allele, InvalidPolicy::Fail);
    assert_eq!(strip_raw(read_all(&bcf_bytes, forced).1), recs);
    let mut plain = bgzf_decompress(&bcf_bytes);
    assert_eq!(strip_raw(read_all(&plain, bcf(VcfSpan::Allele)).1), recs);
    plain[4] = 1;
    assert_eq!(strip_raw(read_all(&plain, bcf(VcfSpan::Allele)).1), recs);
}

#[test]
fn reader_metadata_and_batches() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    let r = open_reader_from(
        "x.bcf",
        Box::new(Cursor::new(bcf_bytes.clone())),
        &bcf(VcfSpan::Allele),
    )
    .unwrap();
    assert_eq!(
        (r.name(), r.format(), r.compression()),
        ("x.bcf", Format::Bcf, Compression::Bgzf)
    );
    let all = read_all(&bcf_bytes, bcf(VcfSpan::Allele)).1;
    let mut threaded = bcf(VcfSpan::Allele);
    threaded.threads = NonZeroUsize::new(4).unwrap();
    assert_eq!(read_all(&bcf_bytes, threaded).1, all, "parallel BGZF input");
    for max in [1, 2, 4] {
        let o = bcf(VcfSpan::Allele);
        let got = run_with(&bcf_bytes, &o, &mut ChromTable::lexicographic(), max)
            .unwrap()
            .1;
        assert_eq!(got, all, "max {max}");
    }
}

#[test]
fn header_and_dictionaries_match_noodles() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    let h = read_header(&bcf_bytes);
    assert_eq!(h.raw, header_text(FIXTURE_VCF).as_bytes());
    assert_eq!(
        h.contigs,
        vec![
            (b"chr2".to_vec(), Some(100000)),
            (b"chr1".to_vec(), Some(100000))
        ]
    );
    assert_eq!(h.samples, vec![b"s1".to_vec(), b"s2".to_vec()]);
    let d = h.bcf.unwrap();
    assert_eq!(d.contigs, noodles_contig_names(FIXTURE_VCF));
    assert_eq!(d.contigs, vec![b"chr1".to_vec(), b"chr2".to_vec()]);
    assert_eq!(d.strings, noodles_string_names(FIXTURE_VCF));
    assert_eq!(d.strings[0], b"PASS");
    // text VCF headers carry no BCF dictionaries
    let Header::Vcf(t) = read_all(FIXTURE_VCF, vcf(VcfSpan::Allele)).0 else {
        panic!()
    };
    assert!(t.bcf.is_none());
}

#[test]
fn idx_overrides_and_dedup_in_string_dictionary() {
    // IDX= everywhere, out of header order; FORMAT/DP shares INFO/DP's slot
    let v = b"##fileformat=VCFv4.3\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\",IDX=4>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=2>\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=2>\n\
##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\",IDX=1>\n\
##FILTER=<ID=q10,Description=\"Quality below 10\",IDX=3>\n\
##contig=<ID=c1,IDX=0>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
c1\t10\t.\tA\t<*>\t.\tq10\tDP=3;END=20\tGT:DP\t0/0:3\n";
    let bcf_bytes = fixtures::vcf_to_bcf(v);
    let d = read_header(&bcf_bytes).bcf.unwrap();
    assert_eq!(d.strings, noodles_string_names(v));
    assert_eq!(
        d.strings,
        [&b"PASS"[..], b"END", b"DP", b"q10", b"GT"].map(<[u8]>::to_vec)
    );
    assert_eq!(
        strip_raw(read_all(&bcf_bytes, bcf(VcfSpan::Allele)).1),
        vec![("c1".to_string(), 9, iv(9, 20))]
    );
}

#[test]
fn raw_records_pass_through_byte_for_byte() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    let copied = copy_via_raw(&bcf_bytes);
    assert_eq!(bgzf_decompress(&copied), bgzf_decompress(&bcf_bytes));
    let records = noodles_records(&bcf_bytes);
    assert_eq!(records.len(), 9);
    assert_eq!(noodles_records(&copied), records);
    // the raw record is the whole framed record: l_shared + l_indiv + 8
    for r in read_all(&bcf_bytes, bcf(VcfSpan::Allele)).1 {
        let l_shared = u32::from_le_bytes(r.4[0..4].try_into().unwrap()) as usize;
        let l_indiv = u32::from_le_bytes(r.4[4..8].try_into().unwrap()) as usize;
        assert_eq!(r.4.len(), 8 + l_shared + l_indiv);
        assert!(l_indiv > 0, "two samples");
    }
}

#[test]
fn write_bcf_header_framing() {
    let mut out = Vec::new();
    write_bcf_header(&mut out, b"##fileformat=VCFv4.3\n#CHROM\n").unwrap();
    let mut want = b"BCF\x02\x02".to_vec();
    want.extend_from_slice(&29u32.to_le_bytes());
    want.extend_from_slice(b"##fileformat=VCFv4.3\n#CHROM\n\0");
    assert_eq!(out, want);
    // read_bcf_header reads it back with empty dictionaries (PASS only)
    let h = read_bcf_header(&mut &out[..]).unwrap();
    assert_eq!(h.raw, b"##fileformat=VCFv4.3\n#CHROM\n");
    let d = h.bcf.unwrap();
    assert_eq!((d.contigs.len(), d.strings), (0, vec![b"PASS".to_vec()]));
    // CRLF header text, NUL padding after the terminator
    let mut out = Vec::new();
    write_bcf_header(
        &mut out,
        b"##fileformat=VCFv4.3\r\n##contig=<ID=c1,length=5>\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\r\n\0\0",
    )
    .unwrap();
    let h = read_bcf_header(&mut &out[..]).unwrap();
    assert_eq!(h.contigs, vec![(b"c1".to_vec(), Some(5))]);
    assert_eq!(h.samples, vec![b"s1".to_vec()]);
    assert_eq!(h.bcf.unwrap().contigs, vec![b"c1".to_vec()]);
}

#[test]
fn truncated_record_names_its_number() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    let mut t = bgzf_decompress(&bcf_bytes);
    t.truncate(t.len() - 5);
    let e = err_bcf(fixtures::bgzip(&t));
    assert!(e.contains("record 9"), "{e}");
    assert!(e.starts_with("x.bcf record 9: record is truncated"), "{e}");
    // cut inside the 8-byte length prefix of record 9
    let last = read_all(&bcf_bytes, bcf(VcfSpan::Allele)).1[8].4.len();
    let mut t = bgzf_decompress(&bcf_bytes);
    t.truncate(t.len() - last + 3);
    assert!(err_bcf(&t).starts_with("x.bcf record 9: record is truncated"));
}

#[test]
fn header_only_and_header_errors() {
    let header_only = fixtures::vcf_to_bcf(header_text(FIXTURE_VCF).as_bytes());
    let (h, recs, _) = read_all(&header_only, bcf(VcfSpan::Allele));
    assert!(recs.is_empty());
    let Header::Vcf(h) = h else { panic!() };
    assert!(h.bcf.is_some());

    let forced = opts(Some(Format::Bcf), VcfSpan::Allele, InvalidPolicy::Fail);
    let e = run(FIXTURE_VCF, &forced).unwrap_err().to_string();
    assert!(e.starts_with("x.bcf: not a BCF file"), "{e}");

    let plain = bgzf_decompress(&header_only);
    let e = err_bcf(&plain[..40]);
    assert_eq!(e, "x.bcf: BCF header is truncated");
    let e = err_bcf(&plain[..7]);
    assert_eq!(e, "x.bcf: BCF header is truncated");

    let mut v3 = plain.clone();
    v3[3] = 3;
    let e = run(&v3, &forced).unwrap_err().to_string();
    assert!(e.starts_with("x.bcf: unsupported BCF version 3.2"), "{e}");
    let mut v23 = plain.clone();
    v23[4] = 3;
    let e = err_bcf(&v23);
    assert!(e.contains("unsupported BCF version 2.3"), "{e}");

    // the BCF header text goes through the same checks as a VCF header
    let mut no_chrom = Vec::new();
    write_bcf_header(&mut no_chrom, b"##fileformat=VCFv4.3\n").unwrap();
    assert_eq!(err_bcf(&no_chrom), "x.bcf: VCF header has no #CHROM line");
    let mut bad_idx = Vec::new();
    write_bcf_header(
        &mut bad_idx,
        b"##fileformat=VCFv4.3\n##contig=<ID=a,IDX=0>\n##contig=<ID=b,IDX=0>\n#CHROM\n",
    )
    .unwrap();
    let e = err_bcf(&bad_idx);
    assert!(e.starts_with("x.bcf: ") && e.contains("IDX=0"), "{e}");
}

#[test]
fn invalid_records_follow_the_policy() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    let mut plain = bgzf_decompress(&bcf_bytes);
    let chrom = first_record(&plain) + 8;
    plain[chrom..chrom + 4].copy_from_slice(&7i32.to_le_bytes());
    let e = err_bcf(&plain);
    assert!(
        e.starts_with("x.bcf record 1: CHROM index 7 is not in the contig dictionary"),
        "{e}"
    );
    let skip = opts(None, VcfSpan::Allele, InvalidPolicy::Skip);
    let (_, recs, d) = read_all(&plain, skip);
    assert_eq!(recs.len(), 8);
    assert_eq!(recs[0].3, 2);
    assert_eq!(d.skipped, 1);
    assert!(
        d.examples[0].starts_with("record 1: CHROM index 7"),
        "{d:?}"
    );

    // POS -1 on disk is VCF POS 0
    let mut plain = bgzf_decompress(&bcf_bytes);
    let pos = first_record(&plain) + 12;
    plain[pos..pos + 4].copy_from_slice(&(-1i32).to_le_bytes());
    assert!(err_bcf(&plain).contains("record 1: POS 0"));

    // a symbolic ALT without END or SVLEN is the same error as in text
    let v = format!(
        "{}chr1\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL\tGT\t0/1\t0/0\n",
        header_text(FIXTURE_VCF)
    );
    let e = err_bcf(fixtures::vcf_to_bcf(v.as_bytes()));
    assert_eq!(
        e,
        "x.bcf record 1: symbolic ALT <DEL> needs INFO/END or INFO/SVLEN"
    );
}

#[test]
fn closed_tables_map_contigs_by_name() {
    let bcf_bytes = fixtures::vcf_to_bcf(FIXTURE_VCF);
    let genome = Genome::new(vec![(b"chr2".to_vec(), None), (b"chr1".to_vec(), None)]).unwrap();
    let mut chroms = ChromTable::closed(&genome, OrderSource::GenomeFile("g.txt".into()));
    let o = bcf(VcfSpan::Allele);
    let recs = run_with(&bcf_bytes, &o, &mut chroms, 3).unwrap().1;
    assert_eq!(
        recs.iter().map(|r| r.0.as_str()).collect::<Vec<_>>(),
        [
            "chr2", "chr2", "chr2", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"
        ]
    );
    let genome = Genome::new(vec![(b"chr2".to_vec(), None)]).unwrap();
    let mut chroms = ChromTable::closed(&genome, OrderSource::GenomeFile("g.txt".into()));
    let e = run_with(&bcf_bytes, &o, &mut chroms, 1000)
        .unwrap_err()
        .to_string();
    assert_eq!(
        e,
        "x.bcf record 4: chromosome \"chr1\" is not in the genome file g.txt"
    );
}

/// INFO fields of other types before END and SVLEN, alleles of 15 bytes or
/// more (typed-length overflow encoding) and vectors/strings longer than 14
/// elements (overflow-length typed values): the BCF decoder must skip them
/// correctly to find END and SVLEN, and give the footprints of the VCF text path.
const OVERFLOW_VCF: &[u8] = b"##fileformat=VCFv4.3\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=XD,Number=.,Type=Integer,Description=\"Depths\">\n\
##INFO=<ID=STR,Number=1,Type=String,Description=\"A long string\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"SV length\">\n\
##ALT=<ID=DEL,Description=\"Deletion\">\n\
##contig=<ID=chr1,length=100000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tACGTACGTACGTACGTACGT\tA\t.\t.\tDB\n\
chr1\t200\t.\tA\tACGTACGTACGTACGTA\t.\t.\tAF=0.5\n\
chr1\t250\t.\tACGTACGTACGTACG\tACGTACGTACGTACT\t.\t.\t.\n\
chr1\t300\t.\tN\t<DEL>\t.\t.\tDB;AF=0.5;XD=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;STR=aaaaaaaaaaaaaaaaaaaaaaaaa;END=500\n\
chr1\t600\t.\tN\t<DEL>\t.\t.\tDB;AF=0.25;XD=100000,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16;STR=bbbbbbbbbbbbbbbb;SVLEN=-100\n\
chr1\t700\t.\tN\t<DEL>,ACGTACGTACGTACGTACGTACGT\t.\t.\tAF=0.1,0.2;XD=70000;END=900\n\
chr1\t1000\t.\tA\tG\t.\t.\tXD=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15\n";

#[test]
fn overflow_encodings_and_earlier_info_types_match_the_vcf_text_path() {
    let bcf_bytes = fixtures::vcf_to_bcf(OVERFLOW_VCF);
    // the fixture really holds the encodings under test
    let plain = bgzf_decompress(&bcf_bytes);
    assert!(
        plain.windows(2).any(|w| w[0] == 0xF7 && w[1] == 0x11),
        "an allele string with an overflow length (type 0xF7 then an int8 length)"
    );
    for span in [VcfSpan::Allele, VcfSpan::Ref] {
        let from_vcf = strip_raw(read_all(OVERFLOW_VCF, vcf(span)).1);
        let from_bcf = strip_raw(read_all(&bcf_bytes, bcf(span)).1);
        assert_eq!(from_bcf.len(), 7);
        assert_eq!(from_vcf, from_bcf, "{span:?}");
    }
    // spot checks so a shared mistake in both readers cannot hide
    let recs = strip_raw(read_all(&bcf_bytes, bcf(VcfSpan::Allele)).1);
    assert_eq!(
        recs[0],
        ("chr1".into(), 99, iv(100, 119)),
        "19-base deletion"
    );
    assert_eq!(
        recs[1],
        ("chr1".into(), 199, iv(200, 200)),
        "16-base insertion is a point"
    );
    assert_eq!(
        recs[2],
        ("chr1".into(), 249, iv(263, 264)),
        "15-base REF, SNV at its end"
    );
    assert_eq!(
        recs[3],
        ("chr1".into(), 299, iv(300, 500)),
        "END after Flag, Float, int vector, string"
    );
    assert_eq!(
        recs[4],
        ("chr1".into(), 599, iv(600, 700)),
        "SVLEN after overflow-length values"
    );
    assert_eq!(
        recs[5],
        ("chr1".into(), 699, iv(699, 900)),
        "END with a 4-byte int before it"
    );
    let ref_span = strip_raw(read_all(&bcf_bytes, bcf(VcfSpan::Ref)).1);
    assert_eq!(ref_span[0].2, iv(99, 119));
    // the records pass through byte for byte
    assert_eq!(bgzf_decompress(&copy_via_raw(&bcf_bytes)), plain);
    assert_eq!(noodles_records(&bcf_bytes).len(), 7);
}
