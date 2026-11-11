//! `gloci samples` (D135): the rules beyond the samtools goldens (`samples-*-vs-samtools`).
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/samtools")
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[&str]) -> (i32, String, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code().unwrap(),
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
    )
}

fn ok(args: &[&str]) -> String {
    let (code, out, err) = run(args);
    assert_eq!(code, 0, "{err}");
    out
}

#[test]
fn samples_rows() {
    let rg = s(&inputs().join("rg.bam"));
    // Header order, each distinct sample once.
    assert_eq!(
        ok(&["samples", &rg]),
        format!("NA12878\t{rg}\nNA12891\t{rg}\nHG002\t{rg}\n")
    );
    // Another tag; the line with no LB of its own adds nothing, and the one with no SM does.
    assert_eq!(
        ok(&["samples", &rg, "--tag", "LB"]),
        format!("libA\t{rg}\nlibB\t{rg}\nlibC\t{rg}\n")
    );
    assert_eq!(
        ok(&["samples", &rg, "--column-names", "--tag", "LB"]),
        format!("#LB\tPATH\nlibA\t{rg}\nlibB\t{rg}\nlibC\t{rg}\n")
    );
    // A tag must have two characters.
    for bad in ["X", "SMS", ""] {
        let (code, out, err) = run(&["samples", &rg, "--tag", bad]);
        assert_eq!((code, out.as_str()), (2, ""), "{bad}");
        assert!(err.contains("--tag must be two characters"), "{err}");
    }
    // No @RG line, or none with the tag: one line with ".".
    let nosq = s(&inputs().join("nosq.sam"));
    assert_eq!(
        ok(&["samples", &nosq, "--tag", "ID"]),
        format!("g1\t{nosq}\n")
    );
    let nosm = s(&inputs().join("nosm.sam"));
    assert_eq!(ok(&["samples", &nosm]), format!(".\t{nosm}\n"));
    let notext = s(&inputs().join("../aln/notext.bam"));
    assert_eq!(ok(&["samples", &notext]), format!(".\t{notext}\n"));
    // Two files, in the order given; SAM works too.
    let rgs = s(&inputs().join("rg.sam"));
    assert_eq!(
        ok(&["samples", &nosm, &rgs]),
        format!(".\t{nosm}\nNA12878\t{rgs}\nNA12891\t{rgs}\nHG002\t{rgs}\n")
    );
    // The same at any thread count, and -o writes the file.
    let d = tempdir();
    let o = s(&d.path().join("o.tsv"));
    ok(&["samples", &rg, "-o", &o]);
    assert_eq!(
        std::fs::read_to_string(&o).unwrap(),
        format!("NA12878\t{rg}\nNA12891\t{rg}\nHG002\t{rg}\n")
    );
    // No --chrom-map: samples describes the file as it is.
    gloci()
        .args(["samples", &rg, "--chrom-map", "x"])
        .assert()
        .code(2);
    gloci()
        .args(["samples", "/nonexistent/x.bam"])
        .assert()
        .code(1)
        .stderr(contains("/nonexistent/x.bam"));
}

#[test]
fn samples_show_index() {
    let d = tempdir();
    let bam = inputs().join("rg.bam");
    // A fresh index, whatever the checkout's file times say.
    let f = s(&copy_with_fresh_index(
        d.path(),
        &bam,
        &inputs().join("rg.bam.bai"),
    ));
    let out = ok(&["samples", &f, "--show-index", "--column-names"]);
    assert_eq!(
        out,
        format!("#SM\tPATH\tINDEX\nNA12878\t{f}\tY\nNA12891\t{f}\tY\nHG002\t{f}\tY\n")
    );
    let (_, _, err) = run(&["samples", &f, "--show-index"]);
    assert!(!err.contains("older"), "{err}");
    // No index: N. SAM never has one.
    let lone = s(&d.path().join("lone.bam"));
    std::fs::copy(&bam, &lone).unwrap();
    assert!(
        ok(&["samples", &lone, "--show-index"])
            .lines()
            .all(|l| l.ends_with("\tN"))
    );
    let sam = s(&inputs().join("rg.sam"));
    assert!(
        ok(&["samples", &sam, "--show-index"])
            .lines()
            .all(|l| l.ends_with("\tN"))
    );
}

/// A stale index warns once, however many samples the file has, and the column is still Y.
#[test]
fn samples_show_index_stale() {
    let d = tempdir();
    let f = d.path().join("rg.bam");
    std::fs::copy(inputs().join("rg.bam"), &f).unwrap();
    let i = d.path().join("rg.bam.bai");
    std::fs::copy(inputs().join("rg.bam.bai"), &i).unwrap();
    std::fs::File::options()
        .write(true)
        .open(&i)
        .unwrap()
        .set_modified(std::time::SystemTime::now() - std::time::Duration::from_secs(3600))
        .unwrap();
    let f = s(&f);
    let (code, out, err) = run(&["samples", &f, "--show-index"]);
    assert_eq!(code, 0, "{err}");
    assert_eq!(out.lines().count(), 3);
    assert!(out.lines().all(|l| l.ends_with("\tY")), "{out}");
    assert_eq!(err.matches("is older than").count(), 1, "{err}");
}

#[test]
fn samples_match_reference() {
    // A copy beside no index, so that the INDEX column is N and nothing depends on file times.
    let d = tempdir();
    let rg = s(&d.path().join("rg.bam"));
    std::fs::copy(inputs().join("rg.bam"), &rg).unwrap();
    let r1 = s(&inputs().join("ref.fa"));
    let r2 = s(&inputs().join("ref2.fa"));
    let other = s(&inputs().join("other.fa"));
    let rows = |fa: &str| format!("NA12878\t{rg}\t{fa}\nNA12891\t{rg}\t{fa}\nHG002\t{rg}\t{fa}\n");
    assert_eq!(ok(&["samples", &rg, "--match-reference", &r1]), rows(&r1));
    // Two that match: the last given wins, as samtools tries it first.
    assert_eq!(
        ok(&[
            "samples",
            &rg,
            "--match-reference",
            &r1,
            "--match-reference",
            &r2
        ]),
        rows(&r2)
    );
    assert_eq!(
        ok(&[
            "samples",
            &rg,
            "--match-reference",
            &r2,
            "--match-reference",
            &r1
        ]),
        rows(&r1)
    );
    // One that matches and a later one that does not: the match.
    assert_eq!(
        ok(&[
            "samples",
            &rg,
            "--match-reference",
            &r1,
            "--match-reference",
            &other
        ]),
        rows(&r1)
    );
    // Nothing matches: ".".
    assert_eq!(
        ok(&["samples", &rg, "--match-reference", &other]),
        rows(".")
    );
    // With the other columns.
    assert_eq!(
        ok(&[
            "samples",
            &rg,
            "--match-reference",
            &r1,
            "--show-index",
            "--column-names",
            "--tag",
            "LB"
        ]),
        format!(
            "#LB\tPATH\tINDEX\tREFERENCE\nlibA\t{rg}\tN\t{r1}\nlibB\t{rg}\tN\t{r1}\nlibC\t{rg}\tN\t{r1}\n"
        )
    );
    // A FASTA with no .fai is a usage error naming the fix; nothing is printed.
    let bare = d.path().join("bare.fa");
    std::fs::copy(&r1, &bare).unwrap();
    let (code, out, err) = run(&["samples", &rg, "--match-reference", &s(&bare)]);
    assert_eq!((code, out.as_str()), (2, ""));
    assert!(err.contains("gloci index"), "{err}");
    assert!(!d.path().join("bare.fa.fai").exists());
}

/// 200,000 `@SQ` lines against a 200,000-entry `.fai`: names are compared once, not looked up
/// per reference.
#[test]
fn samples_huge_header_match() {
    let d = tempdir();
    let n = 200_000;
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    let mut fai = String::new();
    for i in 0..n {
        sam.push_str(&format!("@SQ\tSN:c{i}\tLN:1000\n"));
        fai.push_str(&format!("c{i}\t1000\t{}\t60\t61\n", i * 10));
    }
    sam.push_str("@RG\tID:a\tSM:s1\n@RG\tID:b\tSM:s2\n");
    let sam_path = plain_file(d.path(), "big.sam", &sam);
    let bam = d.path().join("big.bam");
    gloci()
        .args(["view", "-O", "bam", &s(&sam_path), "-o", &s(&bam)])
        .assert()
        .success();
    let fa = plain_file(d.path(), "ref.fa", ">c0\nA\n");
    plain_file(d.path(), "ref.fa.fai", &fai);
    // A second FASTA whose last entry is one base short matches nothing.
    let short = fai.replace(&format!("c{}\t1000", n - 1), &format!("c{}\t999", n - 1));
    let fa2 = plain_file(d.path(), "short.fa", ">c0\nA\n");
    plain_file(d.path(), "short.fa.fai", &short);
    let (bam, fa, fa2) = (s(&bam), s(&fa), s(&fa2));
    let t = Instant::now();
    let out = ok(&[
        "samples",
        &bam,
        "--match-reference",
        &fa,
        "--match-reference",
        &fa2,
    ]);
    let took = t.elapsed();
    // A loose bound for unoptimised builds on a loaded box; quadratic work would take minutes.
    assert!(took.as_secs() < 15, "took {took:?}");
    assert_eq!(out, format!("s1\t{bam}\t{fa}\ns2\t{bam}\t{fa}\n"));
}
