//! `gloci nuc` (D112): behaviour beyond the bedtools goldens.
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

#[test]
fn threads_do_not_change_the_output() {
    let (_d, r) = setup();
    for (bed, extra) in [
        ("nuc.bed", vec!["-s", "--seq", "--pattern", "AC"]),
        ("gf-outside.bed", vec!["--skip-invalid"]),
        ("gf.vcf", vec!["--vcf-span", "ref"]),
    ] {
        let outs: Vec<_> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = gloci()
                    .args(["nuc", "--fasta"])
                    .arg(&r)
                    .arg(fa(bed))
                    .args(["-t", t])
                    .args(&extra)
                    .output()
                    .unwrap();
                assert!(o.status.success(), "{extra:?}: {o:?}");
                o.stdout
            })
            .collect();
        assert!(!outs[0].is_empty());
        assert_eq!(outs[0], outs[1]);
        assert_eq!(outs[0], outs[2]);
    }
}

#[test]
fn points_warn_once_with_a_count() {
    let (_d, r) = setup();
    gloci()
        .args(["nuc", "--fasta"])
        .arg(&r)
        .arg(fa("gf-points.bed"))
        .assert()
        .success()
        .stderr(contains("skipped 2 zero-length intervals"));
}

#[test]
fn outside_the_fasta_fails_by_default_and_skips_on_request() {
    let (d, r) = setup();
    let out = d.path().join("out.tsv");
    gloci()
        .args(["nuc", "--fasta"])
        .arg(&r)
        .arg(fa("gf-outside.bed"))
        .arg("-o")
        .arg(&out)
        .assert()
        .code(1)
        .stderr(
            contains("gf-outside.bed line 1")
                .and(contains("chr1:990-1010 ends past the 1000 bases of chr1"))
                .and(contains("--skip-invalid")),
        );
    assert!(!out.exists(), "a failed run leaves no output file");
    gloci()
        .args(["nuc", "--skip-invalid", "--fasta"])
        .arg(&r)
        .arg(fa("gf-missing.bed"))
        .assert()
        .success()
        .stderr(contains(
            "skipped 2 intervals outside the FASTA (first: chrZ is not in",
        ));
}

#[test]
fn usage_errors() {
    let (d, r) = setup();
    let bed = fa("gf.bed");
    let lone = d.path().join("lone.fa");
    std::fs::copy(fa("ref.fa"), &lone).unwrap();
    gloci()
        .args(["nuc", "--fasta"])
        .arg(&lone)
        .arg(&bed)
        .assert()
        .code(2)
        .stderr(contains("gloci index"));
    for (args, msg) in [
        (vec!["--pattern", ""], "--pattern is empty"),
        (vec!["-C"], "--pattern"),
        (vec!["-o", "x.bam"], "not BAM"),
        (vec!["--write-index", "-o", "x.tsv.gz"], "--write-index"),
    ] {
        gloci()
            .current_dir(d.path())
            .args(["nuc", "--fasta"])
            .arg(&r)
            .arg(&bed)
            .args(&args)
            .assert()
            .code(2)
            .stderr(contains(msg));
    }
}

#[test]
fn tsv_output_file() {
    let (d, r) = setup();
    let out = d.path().join("n.tsv");
    gloci()
        .args(["nuc", "--fasta"])
        .arg(&r)
        .arg(fa("gf3.bed"))
        .arg("-o")
        .arg(&out)
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.starts_with("#1_usercol\t2_usercol\t3_usercol\t4_pct_at\t"));
    assert!(text.contains("chr1\t100\t129\t0.931035\t0.068966\t14\t1\t1\t13\t0\t0\t29\n"));
}
