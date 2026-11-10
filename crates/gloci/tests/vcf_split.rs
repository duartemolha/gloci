//! `gloci vcf split`: BCF, stdin, determinism, errors and the round trip
//! with `subtract --match alleles`.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/vcf/split")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn records(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(String::from)
        .collect()
}

const SIMPLE: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";

#[test]
fn split_matches_golden_and_is_thread_independent() {
    let f = s(&fixture("multi.vcf"));
    let one = out(&["vcf", "split", &f, "-t", "1"]);
    assert_eq!(records(&one).len(), 25);
    for t in ["2", "8"] {
        assert_eq!(out(&["vcf", "split", &f, "-t", t]), one);
    }
}

#[test]
fn bcf_input_and_output_and_stdin() {
    let dir = tempdir();
    let f = s(&fixture("multi.vcf"));
    let want = records(&out(&["vcf", "split", &f]));
    let bcf = s(&dir.path().join("in.bcf"));
    out(&["view", &f, "--vcf-span", "ref", "-o", &bcf]);
    // BCF in, VCF out
    assert_eq!(records(&out(&["vcf", "split", &bcf, "-O", "vcf"])), want);
    // BCF in, BCF out
    let obcf = s(&dir.path().join("out.bcf"));
    out(&["vcf", "split", &bcf, "-o", &obcf]);
    assert_eq!(
        records(&out(&["view", &obcf, "--vcf-span", "ref", "-O", "vcf"])),
        want
    );
    // VCF in, BCF out; stdin
    let o2 = s(&dir.path().join("o2.bcf"));
    gloci()
        .args(["vcf", "split", "-", "-o", &o2])
        .write_stdin(std::fs::read(&f).unwrap())
        .assert()
        .success();
    assert_eq!(
        records(&out(&["view", &o2, "--vcf-span", "ref", "-O", "vcf"])),
        want
    );
}

#[test]
fn wrong_number_of_values_is_an_error_that_says_where() {
    gloci()
        .args(["vcf", "split", &s(&fixture("bad-count.vcf"))])
        .assert()
        .code(1)
        .stderr(contains("INFO/RC").and(contains("line 19")));
    // --skip-invalid drops the record instead
    let o = out(&[
        "vcf",
        "split",
        &s(&fixture("bad-count.vcf")),
        "--skip-invalid",
    ]);
    assert!(records(&o).is_empty());
}

use predicates::prelude::*;

#[test]
fn text_output_is_not_binary_and_header_is_the_inputs() {
    let f = fixture("star.vcf");
    let o = out(&["vcf", "split", &s(&f)]);
    let hdr = |t: &str| {
        t.lines()
            .filter(|l| l.starts_with('#'))
            .map(String::from)
            .collect::<Vec<_>>()
    };
    assert_eq!(hdr(&o), hdr(&std::fs::read_to_string(&f).unwrap()));
}

#[test]
fn split_round_trip_with_subtract() {
    // A: C -> A,G   B: C -> G; subtract by allele leaves C -> A.
    let dir = tempdir();
    let a = plain_file(
        dir.path(),
        "a.vcf",
        &format!("{SIMPLE}chr1\t100\t.\tC\tA,G\t.\t.\t.\tGT\t1/2\n"),
    );
    let b = plain_file(
        dir.path(),
        "b.vcf",
        &format!("{SIMPLE}chr1\t100\t.\tC\tG\t.\t.\t.\tGT\t0/1\n"),
    );
    let split = out(&["vcf", "split", &s(&a)]);
    let sp = plain_file(dir.path(), "a.split.vcf", &split);
    let o = gloci()
        .args(["subtract", "-a", "-", "-b", &s(&b), "--match", "alleles"])
        .write_stdin(std::fs::read(&sp).unwrap())
        .assert()
        .success();
    let text = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    let recs = records(&text);
    assert_eq!(recs.len(), 1, "{text}");
    assert!(recs[0].starts_with("chr1\t100\t.\tC\tA\t"), "{}", recs[0]);
    // without the split, the whole multi-allelic record is dropped
    let o = out(&["subtract", "-a", &s(&a), "-b", &s(&b), "--match", "alleles"]);
    assert!(records(&o).is_empty());
}
