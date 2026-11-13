//! A BCF record is written as VCF text byte for byte as bcftools 1.21 writes it: the records
//! of `bcf-text-parity/input.bcf` (see the README beside it) through every command that
//! renders BCF records as text.
mod common;
use common::*;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/bcf-text-parity")
}

fn records(args: &[&str]) -> String {
    let o = gloci()
        .args(args)
        .arg(dir().join("input.bcf"))
        .assert()
        .success()
        .get_output()
        .clone();
    String::from_utf8(o.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect()
}

#[test]
fn view_writes_bcf_records_as_bcftools_does() {
    let want = std::fs::read_to_string(dir().join("expected.txt")).unwrap();
    assert_eq!(want.lines().count(), 19);
    // `--vcf-span ref`: one record has a symbolic allele without END, which the default
    // span refuses.
    let got = records(&["view", "-O", "vcf", "--vcf-span", "ref"]);
    for (g, w) in got.lines().zip(want.lines()) {
        assert_eq!(g, w);
    }
    assert_eq!(got, want);
}

/// `vcf filter` with a filter that keeps every record renders through the same code, at
/// one thread and several.
#[test]
fn vcf_filter_writes_bcf_records_as_bcftools_does() {
    let want = std::fs::read_to_string(dir().join("expected.txt")).unwrap();
    for t in ["1", "4"] {
        let got = records(&["vcf", "filter", "-O", "vcf", "--min-qual", "0", "-t", t]);
        // QUAL `.` fails --min-qual: the records with a QUAL are compared.
        let kept: String = want
            .lines()
            .filter(|l| l.split('\t').nth(5) != Some("."))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(got, kept, "-t {t}");
    }
}
