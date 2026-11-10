//! BCF conformance matrix against bcftools 1.21 (see `bcf-conformance/README.md`):
//! every FORMAT type and Number, missing, cut short, mixed and present, at
//! several FORMAT positions, plus INFO missing values.
//!
//! * the record bytes gloci writes for `input.vcf` equal `expected.bcf`;
//! * gloci reads `expected.bcf` back as the VCF records bcftools prints;
//! * `gloci check` passes on both BCF files.
mod common;
use common::*;
use predicates::prelude::*;
use std::io::Read;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/bcf-conformance")
}

/// The framed records of a BGZF BCF (after the header).
fn records(bcf: &[u8]) -> Vec<Vec<u8>> {
    let mut all = Vec::new();
    noodles_bgzf::io::Reader::new(bcf)
        .read_to_end(&mut all)
        .unwrap();
    assert_eq!(&all[..5], b"BCF\x02\x02");
    let l_text = u32::from_le_bytes(all[5..9].try_into().unwrap()) as usize;
    let mut at = 9 + l_text;
    let mut out = Vec::new();
    while at < all.len() {
        let a = u32::from_le_bytes(all[at..at + 4].try_into().unwrap()) as usize;
        let b = u32::from_le_bytes(all[at + 4..at + 8].try_into().unwrap()) as usize;
        out.push(all[at..at + 8 + a + b].to_vec());
        at += 8 + a + b;
    }
    out
}

fn hex(b: &[u8]) -> String {
    b.iter()
        .map(|x| format!("{x:02x}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn convert(d: &Path, input: &str) -> PathBuf {
    let out = d.join(format!("{input}.bcf"));
    gloci()
        .args(["view", "-O", "bcf"])
        .arg(dir().join(format!("{input}.vcf")))
        .arg("-o")
        .arg(&out)
        .assert()
        .success();
    out
}

/// `input` is `input` (the matrix) or `gt-later` (GT after another key).
fn compare_bytes(input: &str, expected: &str) {
    let d = tempdir();
    let mine = records(&std::fs::read(convert(d.path(), input)).unwrap());
    let theirs = records(&std::fs::read(dir().join(expected)).unwrap());
    assert_eq!(mine.len(), theirs.len());
    let text = std::fs::read_to_string(dir().join(format!("{input}.vcf"))).unwrap();
    let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    let bad: Vec<String> = mine
        .iter()
        .zip(&theirs)
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, (a, b))| format!("{}\n  gloci    {}\n  bcftools {}", lines[i], hex(a), hex(b)))
        .collect();
    assert!(
        bad.is_empty(),
        "{} of {} records differ; first ones:\n{}",
        bad.len(),
        mine.len(),
        bad.iter().take(8).cloned().collect::<Vec<_>>().join("\n")
    );
}

#[test]
fn bcf_records_equal_bcftools() {
    compare_bytes("input", "expected.bcf");
}

#[test]
fn bcf_records_with_gt_after_another_key_equal_bcftools() {
    compare_bytes("gt-later", "gt-later.bcf");
}

#[test]
fn gloci_reads_bcftools_records_as_bcftools_prints_them() {
    let out = gloci()
        .args(["view", "-O", "vcf", "--no-header"])
        .arg(dir().join("expected.bcf"))
        .assert()
        .success();
    let mine = String::from_utf8(out.get_output().stdout.clone()).unwrap();
    let theirs = std::fs::read_to_string(dir().join("expected.vcf")).unwrap();
    let (a, b): (Vec<&str>, Vec<&str>) = (mine.lines().collect(), theirs.lines().collect());
    assert_eq!(a.len(), b.len());
    let bad: Vec<String> = a
        .iter()
        .zip(&b)
        .filter(|(x, y)| x != y)
        .map(|(x, y)| format!("gloci    {x}\nbcftools {y}"))
        .collect();
    assert!(
        bad.is_empty(),
        "{} lines differ:\n{}",
        bad.len(),
        bad.iter().take(8).cloned().collect::<Vec<_>>().join("\n")
    );
}

#[test]
fn check_passes_on_both_bcf_files() {
    let d = tempdir();
    for f in [
        convert(d.path(), "input"),
        dir().join("expected.bcf"),
        convert(d.path(), "gt-later"),
        dir().join("gt-later.bcf"),
    ] {
        gloci().arg("check").arg(&f).assert().success();
    }
}

const HEAD: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ad\">\n\
##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"gl\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n";

fn text(args: &[&std::ffi::OsStr]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// An integer or float vector field that is `.` in every sample and not the
/// last key used to give a count-0 descriptor and corrupt BCF (or a failure).
#[test]
fn all_missing_vector_fields_convert_to_readable_bcf() {
    let d = tempdir();
    let vcf = format!(
        "{HEAD}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
chr1\t1\t.\tA\tG\t.\t.\t.\tGT:AD:DP\t./.:.:5\n\
chr1\t2\t.\tA\tG\t.\t.\t.\tGT:GL:DP\t./.:.:5\n"
    );
    let input = plain_file(d.path(), "t.vcf", &vcf);
    let bcf = d.path().join("t.bcf");
    gloci()
        .args(["view", "-O", "bcf"])
        .arg(&input)
        .arg("-o")
        .arg(&bcf)
        .assert()
        .success();
    gloci().arg("check").arg(&bcf).assert().success();
    let back = text(&[
        "view".as_ref(),
        "-O".as_ref(),
        "vcf".as_ref(),
        "--no-header".as_ref(),
        bcf.as_os_str(),
    ]);
    assert_eq!(
        back,
        "chr1\t1\t.\tA\tG\t.\t.\t.\tGT:AD:DP\t./.:.:5\nchr1\t2\t.\tA\tG\t.\t.\t.\tGT:GL:DP\t./.:.:5\n"
    );
}

/// `vcf filter --samples a --recalc-ac -O bcf` writes `AF=.` for a site with
/// no called allele; that is a scalar-or-vector INFO value `.` in BCF.
#[test]
fn recalc_ac_writes_a_missing_af_to_bcf() {
    let d = tempdir();
    let vcf = format!(
        "{HEAD}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n\
chr1\t10\t.\tA\tG\t.\t.\tAC=1;AN=4;AF=0.25\tGT\t./.\t0/1\n\
chr1\t20\t.\tA\tG\t.\t.\tAC=1;AN=4;AF=0.25\tGT\t0/1\t0/0\n"
    );
    let input = plain_file(d.path(), "r.vcf", &vcf);
    let bcf = d.path().join("r.bcf");
    gloci()
        .args([
            "vcf",
            "filter",
            "--samples",
            "a",
            "--recalc-ac",
            "-O",
            "bcf",
        ])
        .arg(&input)
        .arg("-o")
        .arg(&bcf)
        .assert()
        .success();
    gloci().arg("check").arg(&bcf).assert().success();
    let back = text(&[
        "view".as_ref(),
        "-O".as_ref(),
        "vcf".as_ref(),
        "--no-header".as_ref(),
        bcf.as_os_str(),
    ]);
    assert!(back.contains("AC=0;AN=0;AF=.\tGT\t./."), "{back}");
    assert!(back.contains("AC=1;AN=2;AF=0.5\tGT\t0/1"), "{back}");
}

#[test]
fn a_line_without_samples_for_a_header_with_samples_is_refused() {
    let d = tempdir();
    let vcf = format!(
        "{HEAD}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
chr1\t1\t.\tA\tG\t.\t.\t.\n"
    );
    let input = plain_file(d.path(), "t.vcf", &vcf);
    gloci()
        .args(["view", "-O", "bcf"])
        .arg(&input)
        .assert()
        .failure()
        .stderr(
            predicates::str::contains("line 11")
                .and(predicates::str::contains("no FORMAT and sample columns")),
        );
}
