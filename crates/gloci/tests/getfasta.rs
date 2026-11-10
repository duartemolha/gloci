//! `gloci getfasta` (D111): behaviour beyond the bedtools goldens.
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

fn run(args: &[&str]) -> std::process::Output {
    gloci().args(args).output().unwrap()
}

#[test]
fn threads_do_not_change_the_output() {
    let (_d, r) = setup();
    let r = r.to_str().unwrap();
    for (bed, extra) in [
        ("gf.bed", vec!["-s", "--name"]),
        ("gf12.bed", vec!["--split", "-s", "--bed-out"]),
        ("gf-outside.bed", vec!["--skip-invalid"]),
    ] {
        let bed = fa(bed);
        let outs: Vec<_> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let mut a = vec!["getfasta", "--fasta", r, bed.to_str().unwrap(), "-t", t];
                a.extend(extra.iter());
                let o = run(&a);
                assert!(o.status.success(), "{extra:?}: {:?}", o);
                o.stdout
            })
            .collect();
        assert!(!outs[0].is_empty());
        assert_eq!(outs[0], outs[1]);
        assert_eq!(outs[0], outs[2]);
    }
}

#[test]
fn points_are_skipped_with_one_counted_warning() {
    let (_d, r) = setup();
    gloci()
        .args(["getfasta", "--fasta"])
        .arg(&r)
        .arg(fa("gf-points.bed"))
        .assert()
        .success()
        .stderr(contains("skipped 2 zero-length intervals"));
    let o = run(&[
        "getfasta",
        "--fasta",
        r.to_str().unwrap(),
        fa("gf-points.bed").to_str().unwrap(),
    ]);
    let text = String::from_utf8(o.stdout).unwrap();
    assert_eq!(text.matches('>').count(), 2, "{text}");
    assert_eq!(text.lines().next(), Some(">chr1:20-25"));
    assert_eq!(text.lines().nth(1), Some("GTGTC"));
}

#[test]
fn outside_the_fasta_fails_by_default_and_skips_on_request() {
    let (d, r) = setup();
    let out = d.path().join("out.fa");
    gloci()
        .args(["getfasta", "--fasta"])
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
        .args(["getfasta", "--fasta"])
        .arg(&r)
        .arg(fa("gf-missing.bed"))
        .assert()
        .code(1)
        .stderr(contains("chrZ is not in the FASTA index").and(contains("--skip-invalid")));
    // Skipped, counted once, with the first cause.
    gloci()
        .args(["getfasta", "--skip-invalid", "--name-only", "--fasta"])
        .arg(&r)
        .arg(fa("gf-missing.bed"))
        .assert()
        .success()
        .stdout(predicate::str::starts_with(">ok\nGTGTC\n>ok2\n"))
        .stderr(contains(
            "skipped 2 intervals outside the FASTA (first: chrZ is not in",
        ));
}

#[test]
fn bgzf_output_and_chrom_map() {
    let (d, r) = setup();
    let map = d.path().join("map.tsv");
    std::fs::write(&map, "one\tchr1\ntwo\tchr2\n").unwrap();
    let bed = d.path().join("in.bed");
    std::fs::write(&bed, "one\t20\t25\tx\t0\t+\ntwo\t0\t5\ty\t0\t-\n").unwrap();
    let gz = d.path().join("out.fa.gz");
    gloci()
        .args(["getfasta", "--fasta"])
        .arg(&r)
        .arg(&bed)
        .arg("--chrom-map")
        .arg(&map)
        .args(["-s", "-o"])
        .arg(&gz)
        .assert()
        .success();
    let bytes = std::fs::read(&gz).unwrap();
    // A BGZF block: gzip magic and the BC extra field.
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 8, 4]);
    assert_eq!(&bytes[12..14], b"BC");
    // The renamed intervals read the FASTA's sequences; gloci reads it back.
    gloci()
        .args(["getfasta", "--fasta"])
        .arg(&r)
        .arg(&bed)
        .arg("--chrom-map")
        .arg(&map)
        .arg("-s")
        .assert()
        .success()
        .stdout(">chr1:20-25(+)\nGTGTC\n>chr2:0-5(-)\nGGCGT\n");
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut noodles_bgzf::io::Reader::new(std::fs::File::open(&gz).unwrap()),
        &mut text,
    )
    .unwrap();
    assert_eq!(text, ">chr1:20-25(+)\nGTGTC\n>chr2:0-5(-)\nGGCGT\n");
}

#[test]
fn usage_errors() {
    let (d, r) = setup();
    let bed = fa("gf.bed");
    // No index: exit 2 naming gloci index.
    let lone = d.path().join("lone.fa");
    std::fs::copy(fa("ref.fa"), &lone).unwrap();
    gloci()
        .args(["getfasta", "--fasta"])
        .arg(&lone)
        .arg(&bed)
        .assert()
        .code(2)
        .stderr(contains("gloci index"));
    for (args, msg) in [
        (vec!["--tab", "--bed-out"], "cannot be used with"),
        (vec!["--name", "--name-only"], "cannot be used with"),
        (vec!["--write-index", "-o", "x.fa.gz"], "--write-index"),
        (vec!["-o", "x.bam"], "not BAM"),
    ] {
        gloci()
            .current_dir(d.path())
            .args(["getfasta", "--fasta"])
            .arg(&r)
            .arg(&bed)
            .args(&args)
            .assert()
            .code(2)
            .stderr(contains(msg));
    }
}

#[test]
fn bed12_without_twelve_columns_is_read_whole_and_bad_blocks_fail() {
    let (d, r) = setup();
    let bed = d.path().join("m.bed");
    std::fs::write(&bed, "chr1\t20\t25\tx\t0\t+\n").unwrap();
    gloci()
        .args(["getfasta", "--split", "--fasta"])
        .arg(&r)
        .arg(&bed)
        .assert()
        .success()
        .stdout(">chr1:20-25\nGTGTC\n");
    std::fs::write(&bed, "chr1\t20\t30\tx\t0\t+\t20\t30\t0\t2\t5,5\t0,4\n").unwrap();
    gloci()
        .args(["getfasta", "--split", "--fasta"])
        .arg(&r)
        .arg(&bed)
        .assert()
        .code(1)
        .stderr(contains("invalid BED12"));
}
