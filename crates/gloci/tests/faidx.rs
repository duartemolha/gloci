//! `gloci faidx` (D131): the command's own rules beyond the samtools goldens.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fa(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/fa")
        .join(name)
}

/// ref.fa and its index in a fresh directory.
fn setup() -> (tempfile::TempDir, PathBuf) {
    let d = tempdir();
    let p = d.path().join("ref.fa");
    std::fs::copy(fa("ref.fa"), &p).unwrap();
    std::fs::copy(fa("ref.fa.fai"), d.path().join("ref.fa.fai")).unwrap();
    (d, p)
}

fn out(args: &[&str]) -> std::process::Output {
    gloci().args(args).output().unwrap()
}

#[test]
fn no_region_exits_2_naming_index() {
    let (_d, r) = setup();
    gloci()
        .args(["faidx"])
        .arg(&r)
        .assert()
        .code(2)
        .stdout("")
        .stderr(contains("gloci index"));
}

#[test]
fn missing_index_exits_2_and_builds_nothing() {
    let d = tempdir();
    let p = d.path().join("input.fa");
    std::fs::copy(fa("ref.fa"), &p).unwrap();
    gloci()
        .args(["faidx"])
        .arg(&p)
        .arg("chr1:1-5")
        .assert()
        .code(2)
        .stderr(contains("gloci index").and(contains("never builds")));
    assert!(!d.path().join("input.fa.fai").exists());
}

#[test]
fn unknown_sequence_fails_before_any_output() {
    let (_d, r) = setup();
    gloci()
        .args(["faidx"])
        .arg(&r)
        .args(["chr1:1-5", "nope:1-5", "chrM:1-4"])
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains("nope").and(contains("--skip-invalid")));
    gloci()
        .args(["faidx"])
        .arg(&r)
        .args(["chr1:1-5", "nope:1-5", "chrM:1-4", "--skip-invalid"])
        .assert()
        .success()
        .stdout(">chr1:1-5\nACGCC\n>chrM:1-4\nacgt\n")
        .stderr(contains("skipped 1 region"));
}

#[test]
fn bad_regions_are_usage_errors() {
    let (_d, r) = setup();
    for bad in [
        "chr1:0-5",
        "chr1:5-1",
        "chr1:abc",
        "chr1:1-2-3",
        ":5-6",
        "{chr1:1-5",
    ] {
        let o = out(&["faidx", r.to_str().unwrap(), "chr1:1-5", bad]);
        assert_eq!(o.status.code(), Some(2), "{bad}");
        assert!(o.stdout.is_empty(), "{bad}: wrote output before failing");
    }
    gloci()
        .args(["faidx"])
        .arg(&r)
        .args(["chr1:1-2", "--mark-strand", "plus"])
        .assert()
        .code(2)
        .stderr(contains("custom,POS,NEG"));
}

#[test]
fn open_start_is_one() {
    let (_d, r) = setup();
    let a = out(&["faidx", r.to_str().unwrap(), "chr1:-10"]);
    let b = out(&["faidx", r.to_str().unwrap(), "chr1:1-10"]);
    assert_eq!(
        String::from_utf8_lossy(&a.stdout),
        ">chr1:-10\nACGCCGGTAG\n"
    );
    assert_eq!(
        a.stdout[a.stdout.iter().position(|&c| c == b'\n').unwrap()..],
        b.stdout[b.stdout.iter().position(|&c| c == b'\n').unwrap()..]
    );
}

#[test]
fn edges_warn_and_continue() {
    let (_d, r) = setup();
    gloci()
        .args(["faidx"])
        .arg(&r)
        .args(["chr1:995-1100", "chr1:1001-1010", "chrM:140-150"])
        .assert()
        .success()
        .stdout(">chr1:995-1100\nGAGGCG\n>chr1:1001-1010\n>chrM:140-150\nACTTTAGTACC\n")
        .stderr(contains("ends past").and(contains("starts past")));
}

#[test]
fn regions_file_is_bed_with_a_counted_warning() {
    let (d, r) = setup();
    let bed = d.path().join("r.bed");
    std::fs::write(&bed, "chr1\t4\t20\nchr2\t10\t10\nchrM\t3\t3\nchrM\t0\t4\n").unwrap();
    gloci()
        .args(["faidx"])
        .arg(&r)
        .arg("--regions-file")
        .arg(&bed)
        .assert()
        .success()
        .stdout(">chr1:5-20\nCGGTAGGGTGTCGGAG\n>chrM:1-4\nacgt\n")
        .stderr(contains("2 zero-length regions"));
    // a samtools-style text file is not BED
    let txt = d.path().join("r.txt");
    std::fs::write(&txt, "chr1:5-70\n").unwrap();
    gloci()
        .args(["faidx"])
        .arg(&r)
        .arg("--regions-file")
        .arg(&txt)
        .assert()
        .failure();
}

#[test]
fn threads_and_bgzf_output_with_index() {
    let (d, r) = setup();
    let mut outs = Vec::new();
    for t in ["1", "2", "8"] {
        let o = d.path().join(format!("o{t}.fa.gz"));
        gloci()
            .args(["faidx"])
            .arg(&r)
            .args(["chr1:5-300", "chr2:1-500", "chrM", "--reverse-complement"])
            .args(["-t", t, "--write-index", "-o"])
            .arg(&o)
            .assert()
            .success();
        outs.push((
            std::fs::read(&o).unwrap(),
            std::fs::read(format!("{}.fai", o.display())).unwrap(),
            std::fs::read(format!("{}.gzi", o.display())).unwrap(),
        ));
    }
    assert!(outs.windows(2).all(|w| w[0] == w[1]), "differs across -t");
    // the output is a FASTA with a usable index: its last record reads back
    let plain = out(&[
        "faidx",
        r.to_str().unwrap(),
        "chr1:5-300",
        "chr2:1-500",
        "chrM",
        "--reverse-complement",
    ]);
    let text = String::from_utf8(plain.stdout).unwrap();
    let want = &text[text.find(">chrM/rc").unwrap()..];
    let o = d.path().join("o1.fa.gz");
    let back = out(&["faidx", o.to_str().unwrap(), "chrM/rc"]);
    assert_eq!(String::from_utf8(back.stdout).unwrap(), want);
}

#[test]
fn output_format_must_be_fasta() {
    let (_d, r) = setup();
    gloci()
        .args(["faidx"])
        .arg(&r)
        .args(["chr1:1-5", "-O", "bed"])
        .assert()
        .code(2)
        .stderr(contains("FASTA"));
}
