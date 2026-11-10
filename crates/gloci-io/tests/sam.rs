//! SAM input, `AlnSource`, the alignment writers and the `@PG` line (D104, D116, D120).

use gloci_core::{ChromTable, VcfSpan};
use gloci_io::*;
use std::io::Cursor;
use std::num::NonZeroUsize;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn ropts(ao: AlignOptions, invalid: InvalidPolicy) -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(ao),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid,
    }
}

fn source(name: &str, ao: AlignOptions) -> Box<dyn AlnSource> {
    open_aln_source(
        fixture(name).to_str().unwrap(),
        &ropts(ao, InvalidPolicy::Fail),
    )
    .unwrap()
}

/// Every record of the source, as `(record number, raw bytes)`.
fn all(src: &mut dyn AlnSource) -> Vec<(u64, Vec<u8>)> {
    let mut v = Vec::new();
    let mut buf = Vec::new();
    while let Some(n) = src.next_record(&mut buf).unwrap() {
        v.push((n, std::mem::take(&mut buf)));
    }
    v
}

fn names(src: &dyn AlnSource, recs: &[(u64, Vec<u8>)]) -> Vec<String> {
    recs.iter()
        .map(|(_, r)| {
            let a = aln_ref(src.layout(), r, src.header()).unwrap();
            String::from_utf8_lossy(a.name()).into_owned()
        })
        .collect()
}

#[test]
fn aln_source_keeps_every_record() {
    let mut expected: Option<Vec<String>> = None;
    for (file, layout) in [("flags.bam", Format::Bam), ("flags.sam", Format::Sam)] {
        let mut s = source(file, AlignOptions::default());
        assert_eq!(s.layout(), layout);
        let recs = all(s.as_mut());
        let got = names(s.as_ref(), &recs);
        assert_eq!(got.len(), 14, "{file}");
        // File order, unmapped reads included, numbered from 1.
        assert_eq!(got[10], "r_unmapped_placed");
        assert_eq!(got[13], "r_unmapped_unplaced");
        assert_eq!(
            recs.iter().map(|r| r.0).collect::<Vec<_>>(),
            (1..=14).collect::<Vec<_>>()
        );
        assert_eq!(s.diagnostics().unmapped, 0, "{file}");
        match &expected {
            None => expected = Some(got),
            Some(e) => assert_eq!(&got, e, "BAM and SAM give the same records"),
        }
        // `--exclude-flags UNMAP` drops the unmapped ones, and counts them as filtered.
        let ao = AlignOptions {
            exclude_flags: 0x4,
            ..AlignOptions::default()
        };
        let mut s = source(file, ao);
        let recs = all(s.as_mut());
        let kept = names(s.as_ref(), &recs);
        assert_eq!(kept.len(), 12, "{file}");
        assert!(kept.iter().all(|n| !n.starts_with("r_unmapped")));
        assert_eq!(s.diagnostics().filtered, 2);
        // The numbers are still the file's.
        assert_eq!(recs.last().unwrap().0, 13);
    }
}

/// `(chromosome, key, start, end, strand, record number)` of a kept record.
type Meta = (String, u64, u64, u64, i8, u64);

fn metas(name: &str, ao: AlignOptions) -> (Vec<Meta>, Diagnostics) {
    let o = ropts(ao, InvalidPolicy::Fail);
    let mut r = open_reader(fixture(name).to_str().unwrap(), &o).unwrap();
    let mut table = chrom_table_for(None, &[&*r], gloci_core::ChromOrder::Karyotype).unwrap();
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut table, &mut batch, 4).unwrap();
        for m in &batch.recs {
            out.push((
                String::from_utf8_lossy(table.name(m.chrom)).into_owned(),
                m.key,
                m.span.start,
                m.span.end,
                m.strand as i8,
                m.loc,
            ));
        }
        if !more {
            break;
        }
    }
    (out, r.diagnostics().clone())
}

#[test]
fn sam_reader_equals_bam_reader() {
    for ao in [
        AlignOptions::default(),
        AlignOptions {
            min_mapq: 20,
            exclude_flags: 0x100,
            ..AlignOptions::default()
        },
    ] {
        for set in ["flags", "cigars", "order", "tags"] {
            let (mut b, bd) = metas(&format!("{set}.bam"), ao.clone());
            let (mut s, sd) = metas(&format!("{set}.sam"), ao.clone());
            assert!(!b.is_empty(), "{set}");
            // A BAM record is placed by its record number, a SAM record by its file line
            // (D101); everything else is equal.
            let text = std::fs::read_to_string(fixture(&format!("{set}.sam"))).unwrap();
            let lines: Vec<u64> = text
                .lines()
                .enumerate()
                .filter(|(_, l)| !l.starts_with('@'))
                .map(|(i, _)| i as u64 + 1)
                .collect();
            for m in &mut s {
                assert!(lines.contains(&m.5), "{set}: {} is a record line", m.5);
                m.5 = 0;
            }
            for m in &mut b {
                m.5 = 0;
            }
            assert_eq!(b, s, "{set}: metas");
            assert_eq!(
                (bd.unmapped, bd.filtered, bd.skipped),
                (sd.unmapped, sd.filtered, sd.skipped),
                "{set}: skip counts"
            );
        }
    }
}

/// The reader's raw SAM records are the file's lines, byte for byte.
#[test]
fn sam_reader_keeps_the_line() {
    let o = ropts(AlignOptions::default(), InvalidPolicy::Fail);
    let text = std::fs::read_to_string(fixture("tags.sam")).unwrap();
    let mut r = open_reader(fixture("tags.sam").to_str().unwrap(), &o).unwrap();
    let mut table = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    r.read_batch(&mut table, &mut batch, 100).unwrap();
    let want: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('@') && l.split('\t').nth(1) != Some("4"))
        .collect();
    assert_eq!(batch.len(), want.len());
    for (i, w) in want.iter().enumerate() {
        assert_eq!(batch.raw(i), w.as_bytes());
    }
}

fn sam_bytes(refs: &str, records: &[&str]) -> Vec<u8> {
    let mut s = format!("@HD\tVN:1.6\n{refs}");
    for r in records {
        s.push_str(r);
        s.push('\n');
    }
    s.into_bytes()
}

fn src_of(bytes: Vec<u8>, invalid: InvalidPolicy) -> Result<Box<dyn AlnSource>, IoError> {
    open_aln_source_from(
        "t.sam",
        Box::new(Cursor::new(bytes)),
        &ropts(AlignOptions::default(), invalid),
    )
}

#[test]
fn sam_errors_name_the_line() {
    let refs = "@SQ\tSN:chr1\tLN:1000\n";
    let good = "ok\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII";
    let cases: [(&str, String, &str); 7] = [
        (
            "RNAME not in @SQ",
            "r\t0\tchrZ\t10\t60\t4M\t*\t0\t0\tACGT\tIIII".into(),
            "RNAME 'chrZ' is not among the @SQ lines",
        ),
        (
            "too few fields",
            "r\t0\tchr1\t10".into(),
            "4 tab-separated fields",
        ),
        (
            "bad FLAG",
            "r\tx\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII".into(),
            "FLAG 'x'",
        ),
        (
            "bad CIGAR",
            "r\t0\tchr1\t10\t60\t4Q\t*\t0\t0\tACGT\tIIII".into(),
            "CIGAR '4Q'",
        ),
        (
            "CIGAR and SEQ disagree",
            "r\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGT\tIIII".into(),
            "covers 5 bases",
        ),
        (
            "header line after records",
            "@SQ\tSN:chr2\tLN:5".into(),
            "header line after the first record",
        ),
        (
            "long name",
            format!(
                "{}\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "n".repeat(255)
            ),
            "255 characters",
        ),
    ];
    for (shape, bad, text) in &cases {
        let bytes = sam_bytes(refs, &[good, bad]);
        let mut s = src_of(bytes.clone(), InvalidPolicy::Fail).unwrap();
        let mut buf = Vec::new();
        assert_eq!(s.next_record(&mut buf).unwrap(), Some(1), "{shape}");
        let e = s.next_record(&mut buf).unwrap_err().to_string();
        assert!(e.contains("t.sam line 4:"), "{shape}: {e}");
        assert!(e.contains(text), "{shape}: {e}");
        // The interval reader names the same line.
        let o = ropts(AlignOptions::default(), InvalidPolicy::Fail);
        let mut r = open_reader_from("t.sam", Box::new(Cursor::new(bytes)), &o).unwrap();
        let mut table = ChromTable::first_appearance();
        let e = r
            .read_batch(&mut table, &mut RecordBatch::default(), 10)
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("t.sam line 4:") && e.contains(text),
            "{shape}: {e}"
        );
    }
    // Under --skip-invalid the bad record is counted with its line, and the rest is read.
    let bytes = sam_bytes(refs, &[good, "r\t0\tchrZ\t10", good]);
    let mut s = src_of(bytes, InvalidPolicy::Skip).unwrap();
    let recs = all(s.as_mut());
    assert_eq!(recs.iter().map(|r| r.0).collect::<Vec<_>>(), [1, 3]);
    assert_eq!(s.diagnostics().skipped, 1);
    assert!(s.diagnostics().examples[0].starts_with("line 4:"));
}

#[test]
fn sam_blank_lines_are_skipped_and_headerless_sam_needs_no_refs() {
    let refs = "@SQ\tSN:chr1\tLN:1000\n";
    let good = "ok\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII";
    let mut bytes = sam_bytes(refs, &[good]);
    bytes.extend_from_slice(b"\n");
    bytes.extend_from_slice(good.as_bytes());
    let mut s = src_of(bytes, InvalidPolicy::Fail).unwrap();
    assert_eq!(all(s.as_mut()).len(), 2);
    assert_eq!(s.diagnostics().blank_lines, 1);
    // Unplaced reads need no @SQ lines.
    let mut s = src_of(
        b"u\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n".to_vec(),
        InvalidPolicy::Fail,
    )
    .unwrap();
    assert_eq!(all(s.as_mut()).len(), 1);
}

#[test]
fn sam_chrom_map_is_refused_with_a_fix() {
    let mut o = ropts(AlignOptions::default(), InvalidPolicy::Fail);
    o.chrom_map = Some(std::sync::Arc::new(
        ChromMap::parse("m.tsv", b"chr1\t1\n").unwrap(),
    ));
    let e = open_aln_source(fixture("flags.sam").to_str().unwrap(), &o)
        .err()
        .unwrap()
        .to_string();
    assert!(
        e.contains("--chrom-map") && e.contains("convert it to BAM"),
        "{e}"
    );
}

fn header_of(text: &str) -> AlnHeader {
    AlnHeader::parse_sam(text.as_bytes().to_vec()).unwrap()
}

fn pg_lines(h: &AlnHeader) -> Vec<String> {
    String::from_utf8_lossy(&h.text)
        .lines()
        .filter(|l| l.starts_with("@PG"))
        .map(str::to_string)
        .collect()
}

#[test]
fn pg_line_ids() {
    // None yet: appended at the end, after @CO.
    let h =
        header_of("@HD\tVN:1.6\n@SQ\tSN:c\tLN:9\n@CO\tnote\n").with_program("gloci view x", "1.2");
    let text = String::from_utf8_lossy(&h.text).into_owned();
    assert!(
        text.ends_with("@CO\tnote\n@PG\tID:gloci\tPN:gloci\tVN:1.2\tCL:gloci view x\n"),
        "{text}"
    );
    // An existing gloci line: ID:gloci.1, PP:gloci, inserted after the last @PG, before @CO.
    let h =
        header_of("@SQ\tSN:c\tLN:9\n@PG\tID:gloci\tPN:gloci\n@CO\tnote\n").with_program("c", "v");
    assert_eq!(
        pg_lines(&h),
        [
            "@PG\tID:gloci\tPN:gloci",
            "@PG\tID:gloci.1\tPN:gloci\tPP:gloci\tVN:v\tCL:c"
        ]
    );
    assert!(String::from_utf8_lossy(&h.text).ends_with("@CO\tnote\n"));
    // Again: .2 follows .1, and PP is the end of the chain.
    let h2 = h.with_program("c2", "v");
    let last = pg_lines(&h2).pop().unwrap();
    assert!(
        last.starts_with("@PG\tID:gloci.2\tPN:gloci\tPP:gloci.1\t"),
        "{last}"
    );
    // Two chain ends: the last one in the file is the predecessor (as samtools does).
    let h = header_of("@PG\tID:x\n@PG\tID:y\n@PG\tID:z\tPP:x\n").with_program("c", "v");
    assert!(pg_lines(&h)[3].contains("\tPP:z\t"), "{:?}", pg_lines(&h));
    // A cycle has no end.
    let h = header_of("@PG\tID:x\tPP:y\n@PG\tID:y\tPP:x\n").with_program("c", "v");
    assert!(!pg_lines(&h)[2].contains("PP:"), "{:?}", pg_lines(&h));
    // The references and sort order are kept.
    let base = header_of("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:9\n");
    let w = base.with_program("c", "v");
    assert_eq!((w.refs.clone(), w.sort_order), (base.refs, base.sort_order));
}

fn out_to(path: &std::path::Path, f: OutFormat) -> Output {
    let spec = OutputSpec {
        path: Some(path.to_path_buf()),
        format: Some(f),
        bgzf: f == OutFormat::Bam,
        direct: false,
    };
    open_output(
        &spec,
        f,
        &OutputOptions {
            threads: NonZeroUsize::MIN,
            level: 6,
            write_index: false,
        },
    )
    .unwrap()
}

fn wopts(cmd: Option<&str>, write_header: bool) -> WriterOptions {
    WriterOptions {
        write_header,
        command_line: cmd.map(str::to_string),
        version: "9.9",
    }
}

/// Write every record of `file` as `to`, returning the output file's bytes.
fn write_as(file: &str, to: OutFormat, cmd: Option<&str>) -> Vec<u8> {
    let mut s = source(file, AlignOptions::default());
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out");
    let h = s.header().clone();
    let mut w =
        open_alignment_writer(out_to(&path, to), s.layout(), to, &h, &wopts(cmd, true)).unwrap();
    for (_, r) in all(s.as_mut()) {
        w.write_record(&r).unwrap();
    }
    w.finish().unwrap();
    std::fs::read(&path).unwrap()
}

fn sam_text_without_pg(b: &[u8]) -> String {
    String::from_utf8(b.to_vec())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(|l| format!("{l}\n"))
        .collect()
}

#[test]
fn sam_passthrough_bytes() {
    // SAM to SAM copies each record line exactly; only the @PG line is added.
    for set in ["flags", "cigars", "order", "tags"] {
        let out = write_as(
            &format!("{set}.sam"),
            OutFormat::Sam,
            Some("gloci view x.sam"),
        );
        let want = std::fs::read_to_string(fixture(&format!("{set}.sam"))).unwrap();
        let got = String::from_utf8(out).unwrap();
        assert_eq!(
            got.lines()
                .filter(|l| l.starts_with("@PG\tID:gloci"))
                .count(),
            1
        );
        assert_eq!(sam_text_without_pg(got.as_bytes()), want, "{set}");
        // And BAM to SAM writes the SAM the BAM was made from (samtools' own text).
        let out = String::from_utf8(write_as(&format!("{set}.bam"), OutFormat::Sam, None)).unwrap();
        if set == "tags" {
            // Floats come back in `%g` form, as samtools writes them.
            assert!(out.contains("\tXF:f:0.1\tXG:f:1e-07\t"), "{out}");
            assert!(
                out.contains("\tXJ:f:0.123457\tXK:f:1e+10\tXL:f:1.23457e+08\t"),
                "{out}"
            );
        } else {
            assert_eq!(out, want, "{set} from BAM");
        }
    }
}

/// The BAM records of an output file.
fn bam_records(bgzf: &[u8]) -> (AlnHeader, Vec<Vec<u8>>) {
    let mut s = open_aln_source_from(
        "o.bam",
        Box::new(Cursor::new(bgzf.to_vec())),
        &ropts(AlignOptions::default(), InvalidPolicy::Fail),
    )
    .unwrap();
    let h = (**s.header()).clone();
    (h, all(s.as_mut()).into_iter().map(|r| r.1).collect())
}

#[test]
fn bam_output_is_the_input_records() {
    for set in ["flags", "cigars", "order", "tags"] {
        let from_bam = write_as(&format!("{set}.bam"), OutFormat::Bam, Some("c"));
        let from_sam = write_as(&format!("{set}.sam"), OutFormat::Bam, Some("c"));
        let (hb, rb) = bam_records(&from_bam);
        let (hs, rs) = bam_records(&from_sam);
        let (_, original) = {
            let mut s = source(&format!("{set}.bam"), AlignOptions::default());
            let h = (**s.header()).clone();
            (
                h,
                all(s.as_mut()).into_iter().map(|r| r.1).collect::<Vec<_>>(),
            )
        };
        assert_eq!(rb, original, "{set}: BAM to BAM copies the bytes");
        // SAM to BAM encodes to the bytes samtools encodes.
        assert_eq!(rs, original, "{set}: SAM to BAM");
        assert_eq!((hb.refs.clone(), hb.text.clone()), (hs.refs, hs.text));
        assert!(
            String::from_utf8_lossy(&hb.text).contains("@PG\tID:gloci\tPN:gloci\tVN:9.9\tCL:c")
        );
    }
}

#[test]
fn sam_line_to_bam_matches_the_bam_records() {
    let mut s = source("tags.bam", AlignOptions::default());
    let h = (**s.header()).clone();
    let bam = all(s.as_mut());
    let text = std::fs::read_to_string(fixture("tags.sam")).unwrap();
    let enc = SamToBam::new(&h);
    for (line, (_, raw)) in text.lines().filter(|l| !l.starts_with('@')).zip(&bam) {
        let mut out = Vec::new();
        sam_line_to_bam(line.as_bytes(), &h, &mut out).unwrap();
        assert_eq!(&out, raw, "{line}");
        let mut again = Vec::new();
        enc.encode(line.as_bytes(), &mut again).unwrap();
        assert_eq!(&again, raw);
    }
    // A line that is not a record is an error, not a panic.
    assert!(sam_line_to_bam(b"junk", &h, &mut Vec::new()).is_err());
}

#[test]
fn alignment_writers_refuse_what_cannot_be_written() {
    let h = header_of("@SQ\tSN:c\tLN:9\n");
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("o.bam");
    // BAM needs its header.
    let e = open_alignment_writer(
        out_to(&p, OutFormat::Bam),
        Format::Bam,
        OutFormat::Bam,
        &h,
        &wopts(None, false),
    )
    .err()
    .unwrap()
    .to_string();
    assert!(e.contains("--no-header") && e.contains("BAM"), "{e}");
    // BED output of alignments is a usage error that says how to get BED (D104).
    let e = check_family(Format::Bam, OutFormat::Bed)
        .unwrap_err()
        .to_string();
    assert!(e.contains("gloci convert"), "{e}");
    assert!(check_family(Format::Cram, OutFormat::Sam).is_ok());
    assert!(check_family(Format::Sam, OutFormat::Bam).is_ok());
    // Text records are not alignments.
    assert!(check_family(Format::Bed, OutFormat::Bam).is_err());
    // SAM output with no header writes the records alone.
    let q = dir.path().join("o.sam");
    let mut s = source("flags.sam", AlignOptions::default());
    let mut w = open_alignment_writer(
        out_to(&q, OutFormat::Sam),
        Format::Sam,
        OutFormat::Sam,
        s.header(),
        &wopts(Some("c"), false),
    )
    .unwrap();
    for (_, r) in all(s.as_mut()) {
        w.write_record(&r).unwrap();
    }
    w.finish().unwrap();
    let text = std::fs::read_to_string(&q).unwrap();
    assert!(!text.contains('@') && text.lines().count() == 14);
}

#[test]
fn bam_as_sam_resolves_a_long_cigar_like_samtools() {
    // 70,000 operations live in the CG tag. SAM text shows them as the CIGAR, and not as a tag.
    let mut tag = b"CGBI".to_vec();
    tag.extend_from_slice(&70_000u32.to_le_bytes());
    for _ in 0..35_000 {
        tag.extend_from_slice(&(1u32 << 4).to_le_bytes());
        tag.extend_from_slice(&((1u32 << 4) | 1).to_le_bytes());
    }
    tag.extend_from_slice(b"NMC\x07");
    let mut rec = Vec::new();
    rec.extend_from_slice(&0i32.to_le_bytes());
    rec.extend_from_slice(&99i32.to_le_bytes());
    rec.extend_from_slice(&[5, 60, 0, 0]);
    rec.extend_from_slice(&2u16.to_le_bytes());
    rec.extend_from_slice(&0u16.to_le_bytes());
    rec.extend_from_slice(&70_000u32.to_le_bytes());
    rec.extend_from_slice(&(-1i32).to_le_bytes());
    rec.extend_from_slice(&(-1i32).to_le_bytes());
    rec.extend_from_slice(&0i32.to_le_bytes());
    rec.extend_from_slice(b"long\0");
    rec.extend_from_slice(&(70_000u32 << 4 | 4).to_le_bytes());
    rec.extend_from_slice(&(35_000u32 << 4 | 3).to_le_bytes());
    rec.extend(std::iter::repeat_n(0x11u8, 35_000));
    rec.extend(std::iter::repeat_n(30u8, 70_000));
    rec.extend_from_slice(&tag);
    let h = header_of("@SQ\tSN:chr1\tLN:1000000\n");
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("o.sam");
    let mut w = open_alignment_writer(
        out_to(&path, OutFormat::Sam),
        Format::Bam,
        OutFormat::Sam,
        &h,
        &wopts(None, true),
    )
    .unwrap();
    w.write_record(&rec).unwrap();
    w.finish().unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let line = text.lines().last().unwrap();
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!(f[5].len(), 35_000 * 4, "the resolved CIGAR 1M1I x 35,000");
    assert!(f[5].starts_with("1M1I1M1I"));
    assert_eq!(&f[11..], ["NM:i:7"], "the CG tag is not repeated");
}

#[test]
fn sam_checks_every_field_where_it_reads_them() {
    // The encoder's checks run at read, for SAM output as for BAM output, naming the line and
    // honouring --skip-invalid.
    let refs = "@SQ\tSN:chr1\tLN:1000\n";
    let good = "ok\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII";
    let cases: [(&str, &str, &str); 6] = [
        (
            "bad tag",
            "r\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:Q:1",
            "type 'Q'",
        ),
        (
            "short QUAL",
            "r\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tII",
            "quality string has 2",
        ),
        (
            "unknown RNEXT",
            "r\t0\tchr1\t10\t60\t4M\tchrZ\t5\t0\tACGT\tIIII",
            "RNEXT 'chrZ'",
        ),
        (
            "TLEN not a number",
            "r\t0\tchr1\t10\t60\t4M\t*\t0\tx\tACGT\tIIII",
            "TLEN 'x'",
        ),
        (
            "bad B array",
            "r\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tXB:B:x,1",
            "B field",
        ),
        (
            "empty optional field",
            "r\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i",
            "TAG:TYPE:VALUE",
        ),
    ];
    for (shape, bad, text) in cases {
        let bytes = sam_bytes(refs, &[good, bad]);
        let mut s = src_of(bytes.clone(), InvalidPolicy::Fail).unwrap();
        let mut buf = Vec::new();
        s.next_record(&mut buf).unwrap();
        let e = s.next_record(&mut buf).unwrap_err().to_string();
        assert!(
            e.contains("t.sam line 4:") && e.contains(text),
            "{shape}: {e}"
        );
        assert!(!e.contains("disk"), "{shape}: {e}");
        let o = ropts(AlignOptions::default(), InvalidPolicy::Fail);
        let mut r = open_reader_from("t.sam", Box::new(Cursor::new(bytes.clone())), &o).unwrap();
        let e = r
            .read_batch(
                &mut ChromTable::first_appearance(),
                &mut RecordBatch::default(),
                9,
            )
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("t.sam line 4:") && e.contains(text),
            "{shape}: {e}"
        );
        let mut s = src_of(bytes, InvalidPolicy::Skip).unwrap();
        assert_eq!(all(s.as_mut()).len(), 1, "{shape}");
        assert_eq!(s.diagnostics().skipped, 1, "{shape}");
    }
}

#[test]
fn sam_crlf_is_read_as_lf() {
    let lf = sam_bytes(
        "@SQ\tSN:chr1\tLN:1000\n",
        &[
            "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tXZ:Z:x",
            "b\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
        ],
    );
    let crlf: Vec<u8> = String::from_utf8(lf.clone())
        .unwrap()
        .replace('\n', "\r\n")
        .into_bytes();
    let mut a = src_of(lf, InvalidPolicy::Fail).unwrap();
    let mut b = src_of(crlf, InvalidPolicy::Fail).unwrap();
    assert_eq!(all(a.as_mut()), all(b.as_mut()));
    assert_eq!(a.header().text, b.header().text);
    assert!(!b.header().text.contains(&b'\r'));
}

#[test]
fn chrom_map_renames_sam_for_interval_readers_and_is_refused_for_raw_records() {
    let map = std::sync::Arc::new(ChromMap::parse("m.tsv", b"chr1\t1\n").unwrap());
    let bytes = sam_bytes(
        "@SQ\tSN:chr1\tLN:1000\n",
        &["a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII"],
    );
    let mut o = ropts(AlignOptions::default(), InvalidPolicy::Fail);
    o.chrom_map = Some(map);
    let mut r = open_reader_from("t.sam", Box::new(Cursor::new(bytes.clone())), &o).unwrap();
    let mut table = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    r.read_batch(&mut table, &mut batch, 9).unwrap();
    assert_eq!(table.name(batch.recs[0].chrom), b"1");
    // The record still names chr1; the renamed header finds it through the file's name.
    let Header::Alignment(h) = r.header().clone() else {
        panic!()
    };
    assert_eq!(h.refs[0].0, b"1");
    assert_eq!(h.ref_id(b"chr1"), Some(0));
    assert!(aln_ref(Format::Sam, batch.raw(0), &h).is_ok());
    // Raw records are written as they are, so the rename cannot apply.
    let e = open_aln_source_from("t.sam", Box::new(Cursor::new(bytes)), &o)
        .err()
        .unwrap();
    assert!(matches!(e, IoError::Usage(_)), "{e:?}");
    assert!(e.to_string().contains("convert it to BAM"), "{e}");
}

#[test]
fn ref_id_is_a_map_lookup() {
    let n = 200_000usize;
    let refs: Vec<(Vec<u8>, u64)> = (0..n)
        .map(|i| (format!("scaffold_{i}").into_bytes(), 1000))
        .collect();
    let h = AlnHeader::new(Vec::new(), refs, SortOrder::Unknown);
    let t = std::time::Instant::now();
    for i in (0..n).rev().take(50_000) {
        assert_eq!(h.ref_id(format!("scaffold_{i}").as_bytes()), Some(i));
    }
    assert_eq!(h.ref_id(b"nope"), None);
    // 50,000 lookups of the last names would take minutes as a scan of 200,000.
    assert!(t.elapsed().as_secs() < 10, "{:?}", t.elapsed());
    // A header made equal by text, refs and order is equal, whatever was looked up.
    let same = AlnHeader::new(h.text.clone(), h.refs.clone(), SortOrder::Unknown);
    assert_eq!(h, same);
}

#[test]
fn sam_with_many_contigs_is_not_quadratic() {
    let n = 100_000usize;
    let mut sam = String::from("@HD\tVN:1.6\n");
    for i in 0..n {
        sam.push_str(&format!("@SQ\tSN:scaffold_{i}\tLN:1000\n"));
    }
    for i in (0..n).rev() {
        sam.push_str(&format!(
            "r{i}\t0\tscaffold_{i}\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        ));
    }
    let t = std::time::Instant::now();
    let mut s = src_of(sam.into_bytes(), InvalidPolicy::Fail).unwrap();
    assert_eq!(all(s.as_mut()).len(), n);
    // A scan per record is 10^10 name comparisons; the map makes it a second or two.
    assert!(t.elapsed().as_secs() < 20, "{:?}", t.elapsed());
}

fn sam_to_bam_line(line: &str) -> Result<Vec<u8>, String> {
    let h = header_of("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:50\n");
    let mut out = Vec::new();
    sam_line_to_bam(line.as_bytes(), &h, &mut out).map(|()| out)
}

/// `(ref_id, pos, flag, mate ref_id, mate pos)` of a BAM record.
type Place = (i32, i32, u16, i32, i32);

fn fields(rec: &[u8]) -> Place {
    let i = |at: usize| i32::from_le_bytes(rec[at..at + 4].try_into().unwrap());
    (
        i(0),
        i(4),
        u16::from_le_bytes([rec[14], rec[15]]),
        i(20),
        i(24),
    )
}

#[test]
fn encoder_normalises_like_htslib() {
    // (line, ref_id, pos, flag, mate ref, mate pos)
    let cases: [(&str, Place); 7] = [
        // A mapped read with CIGAR `*` is unmapped.
        (
            "a\t0\tchr1\t10\t60\t*\t*\t0\t0\tACGT\tIIII",
            (0, 9, 4, -1, -1),
        ),
        // A reference with POS 0 is no reference; `=` follows the normalised reference.
        (
            "a\t0\tchr1\t0\t60\t4M\t=\t100\t0\tACGT\tIIII",
            (-1, -1, 4, -1, 99),
        ),
        // PNEXT 0 unsets the mate.
        (
            "a\t0\tchr1\t5\t60\t4M\t=\t0\t0\tACGT\tIIII",
            (0, 4, 0, -1, -1),
        ),
        (
            "a\t0\tchr1\t5\t60\t4M\tchr2\t7\t0\tACGT\tIIII",
            (0, 4, 0, 1, 6),
        ),
        // Mapped with no reference: unmapped.
        (
            "a\t0\t*\t5\t60\t4M\t*\t0\t0\tACGT\tIIII",
            (-1, 4, 4, -1, -1),
        ),
        // Hex and octal FLAG, as strtol base 0.
        (
            "a\t0x10\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII",
            (0, 4, 16, -1, -1),
        ),
        (
            "a\t020\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII",
            (0, 4, 16, -1, -1),
        ),
    ];
    for (line, want) in cases {
        let rec = sam_to_bam_line(line).unwrap();
        assert_eq!(fields(&rec), want, "{line}");
    }
    // POS 2^31 is the last base BAM can hold.
    assert_eq!(
        fields(&sam_to_bam_line("a\t0\tchr1\t2147483648\t60\t4M\t*\t0\t0\tACGT\tIIII").unwrap()).1,
        i32::MAX
    );
    // htslib's base table: IUPAC either case, digits 0 to 3 as A C G T, anything else N.
    let rec = sam_to_bam_line("a\t4\t*\t0\t0\t*\t*\t0\t0\tUu0123=.mY\t*").unwrap();
    let seq = &rec[34..];
    assert_eq!(&seq[..5], [0xff, 0x12, 0x48, 0x0f, 0x3a]);
}

#[test]
fn encoder_accepts_what_htslib_accepts_and_refuses_what_it_refuses() {
    let ok = [
        "a\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tXI:I:5\tXC:c:5",
        "a\t0\tchr1\t5\t60\t2M1B2M\t*\t0\t0\tACGT\tIIII",
    ];
    for l in ok {
        assert!(sam_to_bam_line(l).is_ok(), "{l}");
    }
    let bad = [
        (
            "a\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tXH:H:1AE",
            "odd number of digits",
        ),
        ("a\t4\t*\t0\t0\t4M\t*\t0\t0\tACG\t*", "CIGAR covers 4"),
        (
            "a\t0\tchr1\t5\t60\t4M\t*\t3000000000\t0\tACGT\tIIII",
            "PNEXT 3000000000",
        ),
        (
            "a\t0\tchr1\t5\t60\t4M\t*\t0\t3000000000\tACGT\tIIII",
            "TLEN",
        ),
        (
            "a\t0\tchrZ\t5\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "--skip-invalid",
        ),
        (
            "a\t0\tchr1\t5\t60\t4M\tchrZ\t5\t0\tACGT\tIIII",
            "--skip-invalid",
        ),
    ];
    for (l, text) in bad {
        let e = sam_to_bam_line(l).unwrap_err();
        assert!(e.contains(text), "{l}: {e}");
    }
    // 70,000 operations with a reference span of 2^28 or more cannot be written in BAM.
    let cigar = "8000M1D".repeat(35_000);
    let e = sam_to_bam_line(&format!("a\t0\tchr1\t5\t60\t{cigar}\t*\t0\t0\t*\t*")).unwrap_err();
    assert!(e.contains("cannot be written as BAM"), "{e}");
}

/// A BAM record with the given reference, flag, CIGAR words, `l_seq` and aux bytes, written to
/// SAM text.
fn bam_to_sam_text(
    ref_id: i32,
    flag: u16,
    cigar: &[u32],
    l_seq: u32,
    qual: u8,
    tags: &[u8],
) -> String {
    let mut rec = Vec::new();
    rec.extend_from_slice(&ref_id.to_le_bytes());
    rec.extend_from_slice(&(if ref_id < 0 { -1i32 } else { 9 }).to_le_bytes());
    rec.extend_from_slice(&[2, 60, 0, 0]);
    rec.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    rec.extend_from_slice(&flag.to_le_bytes());
    rec.extend_from_slice(&l_seq.to_le_bytes());
    rec.extend_from_slice(&(-1i32).to_le_bytes());
    rec.extend_from_slice(&(-1i32).to_le_bytes());
    rec.extend_from_slice(&0i32.to_le_bytes());
    rec.extend_from_slice(b"r\0");
    for c in cigar {
        rec.extend_from_slice(&c.to_le_bytes());
    }
    rec.extend(std::iter::repeat_n(0x11u8, (l_seq as usize).div_ceil(2)));
    rec.extend(std::iter::repeat_n(qual, l_seq as usize));
    rec.extend_from_slice(tags);
    let h = header_of("@SQ\tSN:chr1\tLN:1000\n");
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("o.sam");
    let mut w = open_alignment_writer(
        out_to(&path, OutFormat::Sam),
        Format::Bam,
        OutFormat::Sam,
        &h,
        &wopts(None, false),
    )
    .unwrap();
    w.write_record(&rec).unwrap();
    w.finish().unwrap();
    std::fs::read_to_string(&path).unwrap()
}

#[test]
fn bam_to_sam_text_matches_htslib_in_the_odd_places() {
    // A negative NaN is `-nan`, a positive one `nan`.
    let mut tags = b"XfBf".to_vec();
    tags.extend_from_slice(&2u32.to_le_bytes());
    tags.extend_from_slice(&f32::from_bits(0xffc0_0000).to_le_bytes());
    tags.extend_from_slice(&f32::from_bits(0x7fc0_0000).to_le_bytes());
    let t = bam_to_sam_text(0, 0, &[1 << 4], 1, 30, &tags);
    assert!(
        t.trim_end().ends_with("XF:B:f,-nan,nan") || t.contains("Xf:B:f,-nan,nan"),
        "{t}"
    );
    // QUAL + 33 wraps as htslib's does.
    let t = bam_to_sam_text(0, 0, &[1 << 4], 1, 0xfe, b"");
    assert!(t.contains("\t\x1f\n") || t.contains("\t\u{1f}"), "{t:?}");
    // CIGAR op B (9) is read and shown, and consumes nothing.
    let t = bam_to_sam_text(0, 0, &[1 << 4, 3 << 4 | 9, 1 << 4], 2, 30, b"");
    assert!(t.contains("\t1M3B1M\t"), "{t}");
    // An unplaced read keeps its placeholder CIGAR and its CG tag.
    let mut cg = b"CGBI".to_vec();
    cg.extend_from_slice(&1u32.to_le_bytes());
    cg.extend_from_slice(&(1u32 << 4).to_le_bytes());
    let t = bam_to_sam_text(-1, 4, &[1 << 4 | 4, 5 << 4 | 3], 1, 30, &cg);
    assert!(t.contains("\t1S5N\t") && t.contains("CG:B:I,16"), "{t}");
    // A placed read has it resolved and the tag dropped.
    let t = bam_to_sam_text(0, 0, &[1 << 4 | 4, 5 << 4 | 3], 1, 30, &cg);
    assert!(t.contains("\t1M\t") && !t.contains("CG:"), "{t}");
}

#[test]
fn bam_header_without_sq_text_shows_its_references_in_sam_only() {
    let h = AlnHeader::new(
        b"@HD\tVN:1.6\n".to_vec(),
        vec![(b"c1".to_vec(), 100), (b"c2".to_vec(), 50)],
        SortOrder::Unknown,
    );
    let shown = h.with_sq_lines();
    assert_eq!(
        String::from_utf8_lossy(&shown.text),
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\n@SQ\tSN:c2\tLN:50\n"
    );
    // Text that has @SQ lines is as it was.
    assert_eq!(shown.with_sq_lines().text, shown.text);
    let dir = tempfile::tempdir().unwrap();
    for (f, has_sq) in [(OutFormat::Sam, true), (OutFormat::Bam, false)] {
        let path = dir.path().join(format!("o.{f:?}"));
        let w = open_alignment_writer(out_to(&path, f), Format::Bam, f, &h, &wopts(None, true))
            .unwrap();
        w.finish().unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let header = if f == OutFormat::Bam {
            let mut s = open_aln_source_from(
                "o.bam",
                Box::new(Cursor::new(bytes)),
                &ropts(AlignOptions::default(), InvalidPolicy::Fail),
            )
            .unwrap();
            let h = (**s.header()).clone();
            let _ = all(s.as_mut());
            h
        } else {
            header_of(&String::from_utf8(bytes).unwrap())
        };
        assert_eq!(
            String::from_utf8_lossy(&header.text).contains("@SQ\tSN:c1\tLN:100"),
            has_sq,
            "{f:?}"
        );
    }
}

#[test]
fn pg_lines_follow_the_chain_ends_as_htslib_adds_them() {
    let pgs = |text: &str| pg_lines(&header_of(text).with_program("c", "v"));
    let new = |l: &String| {
        l.split('\t')
            .filter(|f| f.starts_with("ID:") || f.starts_with("PP:"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    // One line per chain end (a line with a PP that exists and that no line names).
    let p = pgs("@PG\tID:a\n@PG\tID:b\tPP:a\n@PG\tID:c\n@PG\tID:d\tPP:c\n");
    assert_eq!(
        p[4..].iter().map(new).collect::<Vec<_>>(),
        ["ID:gloci PP:b", "ID:gloci.1 PP:d"]
    );
    // A root and a chain: the chain end only.
    let p = pgs("@PG\tID:bwa\n@PG\tID:aln\tPP:bwa\n@PG\tID:STAR\n");
    assert_eq!(
        p[3..].iter().map(new).collect::<Vec<_>>(),
        ["ID:gloci PP:aln"]
    );
    // No chain end: the last line nobody names (a root, or a PP that is missing).
    assert_eq!(
        pgs("@PG\tID:x\n@PG\tID:y\n")[2..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:y"]
    );
    assert_eq!(
        pgs("@PG\tID:a\tPP:zz\n@PG\tID:b\n")[2..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:b"]
    );
    assert_eq!(
        pgs("@PG\tID:a\tPP:zz\n")[1..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:a"]
    );
    assert_eq!(
        pgs("@PG\tID:a\tPP:zz\n@PG\tID:x\tPP:y\n@PG\tID:y\tPP:x\n")[3..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:a"]
    );
    // A line that names itself is a root: the chain end wins over it, and the last root is
    // the predecessor when there is no chain.
    assert_eq!(
        pgs("@PG\tID:r\n@PG\tID:c\tPP:r\n@PG\tID:s\tPP:s\n")[3..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:c"]
    );
    assert_eq!(
        pgs("@PG\tID:a\tPP:a\n@PG\tID:b\n")[2..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:b"]
    );
    // A cycle with nothing else: a line with no PP. A self-reference is a chain end.
    assert_eq!(
        pgs("@PG\tID:x\tPP:y\n@PG\tID:y\tPP:x\n")[2..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci"]
    );
    assert_eq!(
        pgs("@PG\tID:a\tPP:a\n")[1..]
            .iter()
            .map(new)
            .collect::<Vec<_>>(),
        ["ID:gloci PP:a"]
    );
    // IDs skip every ID in use, among them the ones just added.
    let p = pgs(
        "@PG\tID:gloci\n@PG\tID:gloci.1\tPP:gloci\n@PG\tID:gloci.3\n@PG\tID:w\n@PG\tID:v\tPP:w\n",
    );
    assert_eq!(
        p[5..].iter().map(new).collect::<Vec<_>>(),
        ["ID:gloci.2 PP:gloci.1", "ID:gloci.4 PP:v"]
    );
}

#[test]
fn bam_output_with_a_pg_line_shows_sq_lines_and_without_it_keeps_the_text() {
    let h = AlnHeader::new(
        b"@HD\tVN:1.6\n".to_vec(),
        vec![(b"c1".to_vec(), 100)],
        SortOrder::Unknown,
    );
    let dir = tempfile::tempdir().unwrap();
    for (cmd, want_sq) in [(Some("c"), true), (None, false)] {
        let path = dir.path().join(format!("o{}.bam", cmd.is_some()));
        let w = open_alignment_writer(
            out_to(&path, OutFormat::Bam),
            Format::Bam,
            OutFormat::Bam,
            &h,
            &wopts(cmd, true),
        )
        .unwrap();
        w.finish().unwrap();
        let (hdr, _) = bam_records(&std::fs::read(&path).unwrap());
        let text = String::from_utf8_lossy(&hdr.text).into_owned();
        assert_eq!(text.contains("@SQ\tSN:c1\tLN:100"), want_sq, "{text}");
        assert_eq!(text.contains("@PG\tID:gloci"), cmd.is_some(), "{text}");
    }
}

#[test]
fn sam_source_hands_over_the_bam_bytes_it_checked() {
    let mut s = source("tags.sam", AlignOptions::default());
    let bam: Vec<Vec<u8>> = {
        let mut b = source("tags.bam", AlignOptions::default());
        all(b.as_mut()).into_iter().map(|r| r.1).collect()
    };
    let mut buf = Vec::new();
    let mut i = 0;
    while s.next_record(&mut buf).unwrap().is_some() {
        assert_eq!(s.last_bam(), Some(&bam[i][..]), "record {i}");
        buf.clear();
        i += 1;
    }
    assert_eq!(i, bam.len());
    // A BAM source has nothing to hand over: its raw record is the BAM record.
    assert!(
        source("tags.bam", AlignOptions::default())
            .last_bam()
            .is_none()
    );
}

#[test]
fn normalised_lines_are_written_as_samtools_writes_them() {
    let refs = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n";
    let lines = [
        // Kept as written: spelling samtools would normalise.
        "lower\t0\tchr1\t100\t60\t4M\t*\t0\t0\tacgt\tIIII\tXF:f:0.10\tXP:i:+7",
        // Normalised.
        "nocigar\t0\tchr1\t100\t60\t*\t*\t0\t0\tACGT\tIIII",
        "pos0\t0\tchr1\t0\t60\t4M\t=\t50\t0\tACGT\tIIII",
        "norname\t0\t*\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
        "rnext0\t1\tchr1\t100\t60\t4M\tchr1\t0\t0\tACGT\tIIII",
    ];
    let mut s = src_of(sam_bytes(refs, &lines), InvalidPolicy::Fail).unwrap();
    let got: Vec<String> = all(s.as_mut())
        .into_iter()
        .map(|r| String::from_utf8(r.1).unwrap())
        .collect();
    assert_eq!(got[0], lines[0]);
    assert_eq!(got[1], "nocigar\t4\tchr1\t100\t60\t*\t*\t0\t0\tACGT\tIIII");
    assert_eq!(got[2], "pos0\t4\t*\t0\t60\t4M\t*\t50\t0\tACGT\tIIII");
    assert_eq!(got[3], "norname\t4\t*\t100\t60\t4M\t*\t0\t0\tACGT\tIIII");
    assert_eq!(got[4], "rnext0\t1\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII");
    // The filters judge the flag htslib reads.
    let ao = AlignOptions {
        require_flags: 0x4,
        ..AlignOptions::default()
    };
    let mut s = open_aln_source_from(
        "t.sam",
        Box::new(Cursor::new(sam_bytes(refs, &lines))),
        &ropts(ao, InvalidPolicy::Fail),
    )
    .unwrap();
    assert_eq!(all(s.as_mut()).len(), 3);
}

#[test]
fn records_bam_cannot_hold_are_fine_in_sam_and_an_error_for_bam() {
    let refs = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n";
    let pnext = "p\t1\tchr1\t100\t60\t4M\tchr2\t3000000000\t0\tACGT\tIIII";
    let tlen = "t\t1\tchr1\t100\t60\t4M\tchr2\t100\t3000000000\tACGT\tIIII";
    let mut s = src_of(sam_bytes(refs, &[pnext, tlen]), InvalidPolicy::Fail).unwrap();
    let recs = all(s.as_mut());
    assert_eq!(recs.len(), 2);
    // Nothing for a BAM writer to reuse; encoding it for BAM says why not.
    assert!(s.last_bam().is_none());
    let h = (**s.header()).clone();
    for (l, why) in [(pnext, "PNEXT 3000000000"), (tlen, "TLEN 3000000000")] {
        let e = sam_line_to_bam(l.as_bytes(), &h, &mut Vec::new()).unwrap_err();
        assert!(
            e.contains("cannot be written as BAM") && e.contains(why),
            "{e}"
        );
        assert!(e.contains("-O sam"), "{e}");
    }
    // And a BAM writer reports it as a data error, naming the fix.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("o.bam");
    let mut w = open_alignment_writer(
        out_to(&path, OutFormat::Bam),
        Format::Sam,
        OutFormat::Bam,
        &h,
        &wopts(None, true),
    )
    .unwrap();
    let e = w.write_record(pnext.as_bytes()).unwrap_err();
    assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn flag_with_a_sign_is_refused() {
    for f in ["+16", "-1", "0x+10", "+0x10"] {
        let line = format!("a\t{f}\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII");
        let e = sam_to_bam_line(&line).unwrap_err();
        assert!(e.contains("FLAG"), "{f}: {e}");
    }
}
