//! The TSV results writer: column-name line, header lines, D42 format rule,
//! and the D43 BED-preset index.

use gloci_io::*;
use noodles_core::Position;
use noodles_csi::BinningIndex;
use std::num::NonZeroUsize;
use std::path::Path;

fn opts(write_index: bool) -> OutputOptions {
    OutputOptions {
        threads: NonZeroUsize::new(2).unwrap(),
        level: 6,
        write_index,
    }
}

fn open(dir: &Path, name: &str, write_index: bool) -> Result<Output, IoError> {
    let spec = parse_output_spec(Some(&dir.join(name)), None).unwrap();
    open_output(&spec, spec.format.unwrap(), &opts(write_index))
}

fn write_rows(
    dir: &Path,
    name: &str,
    a: Format,
    appended: &[&str],
    header: TsvHeader,
    rows: &[&str],
) -> Vec<u8> {
    let appended: Vec<String> = appended.iter().map(|s| s.to_string()).collect();
    let mut w =
        open_tsv_writer(open(dir, name, false).unwrap(), a, &appended, header, false).unwrap();
    for r in rows {
        w.write_record(r.as_bytes()).unwrap();
    }
    w.finish().unwrap();
    std::fs::read(dir.join(name)).unwrap()
}

fn samples(s: &[&str]) -> Vec<Vec<u8>> {
    s.iter().map(|x| x.as_bytes().to_vec()).collect()
}

#[test]
fn vcf_a_column_line() {
    let d = tempfile::tempdir().unwrap();
    let got = write_rows(
        d.path(),
        "o.tsv",
        Format::Vcf,
        &["count"],
        TsvHeader::VcfColumns {
            samples: samples(&["s1", "s2"]),
        },
        &["1\t5\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1|1\t3"],
    );
    let text = String::from_utf8(got).unwrap();
    assert_eq!(
        text.lines().next().unwrap(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\tcount"
    );
    assert_eq!(text.lines().count(), 2);
}

#[test]
fn vcf_a_without_samples() {
    let d = tempfile::tempdir().unwrap();
    let got = write_rows(
        d.path(),
        "o.tsv",
        Format::Vcf,
        &["count", "mean"],
        TsvHeader::VcfColumns { samples: vec![] },
        &[],
    );
    assert_eq!(
        got,
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tcount\tmean\n"
    );
}

#[test]
fn bed_a_without_header() {
    let d = tempfile::tempdir().unwrap();
    let got = write_rows(
        d.path(),
        "o.tsv",
        Format::Bed,
        &["count"],
        TsvHeader::None,
        &["chr1\t0\t10\t2", "chr1\t20\t30\t0"],
    );
    assert_eq!(got, b"chr1\t0\t10\t2\nchr1\t20\t30\t0\n");
}

#[test]
fn bed_a_with_header() {
    let d = tempfile::tempdir().unwrap();
    let got = write_rows(
        d.path(),
        "o.tsv",
        Format::Bed,
        &["count"],
        TsvHeader::Lines(b"track name=x\n#note\n".to_vec()),
        &["chr1\t0\t10\t2"],
    );
    assert_eq!(got, b"track name=x\n#note\nchr1\t0\t10\t2\n");
}

#[test]
fn tsv_output_ok_rules() {
    assert!(tsv_output_ok(Format::Bed, OutFormat::Tsv).is_ok());
    assert!(tsv_output_ok(Format::Bed, OutFormat::Bed).is_ok());
    assert!(tsv_output_ok(Format::Bed, OutFormat::BedGraph).is_ok());
    assert!(tsv_output_ok(Format::BedGraph, OutFormat::Bed).is_ok());
    assert!(tsv_output_ok(Format::Vcf, OutFormat::Tsv).is_ok());
    assert!(tsv_output_ok(Format::Gff3, OutFormat::Tsv).is_ok());
    for (a, to) in [
        (Format::Vcf, OutFormat::Vcf),
        (Format::Gff3, OutFormat::Bed),
        (Format::Bed, OutFormat::Vcf),
    ] {
        let msg = tsv_output_ok(a, to).unwrap_err().0;
        assert!(msg.contains("appends columns"), "{msg}");
        assert!(msg.contains("TSV"), "{msg}");
        assert!(msg.contains(".tsv") && msg.contains(".tsv.gz"), "{msg}");
        assert!(msg.contains(".bed"), "{msg}");
    }
}

#[test]
fn tsv_index_for_bed_a() {
    let d = tempfile::tempdir().unwrap();
    let appended = vec!["count".to_string()];
    let mut w = open_tsv_writer(
        open(d.path(), "o.tsv.gz", true).unwrap(),
        Format::Bed,
        &appended,
        TsvHeader::Lines(b"track name=x\n".to_vec()),
        true,
    )
    .unwrap();
    for r in ["chr1\t0\t10\t2", "chr1\t20\t30\t0", "chr2\t5\t9\t1"] {
        w.write_record(r.as_bytes()).unwrap();
    }
    w.finish().unwrap();
    let data = d.path().join("o.tsv.gz");
    let idx = d.path().join("o.tsv.gz.tbi");
    assert!(data.exists() && idx.exists());
    let index = noodles_tabix::fs::read(&idx).unwrap();
    let header = index.header().unwrap();
    let names: Vec<String> = header
        .reference_sequence_names()
        .iter()
        .map(|n| n.to_string())
        .collect();
    assert_eq!(names, ["chr1", "chr2"]);
    assert_eq!(header.reference_sequence_name_index(), 0);
    assert_eq!(header.start_position_index(), 1);
    assert_eq!(header.end_position_index(), Some(2));
    let start = Position::try_from(21).unwrap();
    let end = Position::try_from(25).unwrap();
    let chunks = index.query(0, (start..=end).into()).unwrap();
    assert!(!chunks.is_empty());
}

#[test]
fn tsv_index_refused_for_vcf_a() {
    let d = tempfile::tempdir().unwrap();
    let out = open(d.path(), "o.tsv.gz", true).unwrap();
    let err = open_tsv_writer(out, Format::Vcf, &[], TsvHeader::None, true)
        .err()
        .expect("must refuse");
    let msg = err.to_string();
    assert!(msg.contains("gloci index -s 1 -b 2 -e 2"), "{msg}");
}

#[test]
fn tsv_index_refused_for_gff_a() {
    let d = tempfile::tempdir().unwrap();
    let out = open(d.path(), "o.tsv.gz", true).unwrap();
    let err = open_tsv_writer(out, Format::Gff3, &[], TsvHeader::None, true)
        .err()
        .expect("must refuse");
    assert!(
        err.to_string().contains("gloci index -s 1 -b 4 -e 5"),
        "{err}"
    );
}
