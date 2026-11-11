//! `gloci pair-pair` at the command line: what its inputs must be, its
//! warnings and exit codes, and its output formats. The semantics are pinned
//! in gloci-ops' tests and by the bedtools goldens.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

const A: &str = "chr1\t10\t20\tchr1\t100\t110\tp1\t0\t+\t-\n";
const B: &str = "chr1\t12\t14\tchr1\t102\t104\tb1\t0\t+\t-\n";

fn files(dir: &std::path::Path) -> (String, String) {
    let a = plain_file(dir, "a.bedpe", A);
    let b = plain_file(dir, "b.bedpe", B);
    (
        a.to_str().unwrap().to_string(),
        b.to_str().unwrap().to_string(),
    )
}

#[test]
fn pair_pair_default_is_both() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    gloci()
        .args(["pair-pair", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout(format!("{}\t{}", A.trim_end(), B));
    // The A-alone types write BEDPE; neither finds nothing here.
    gloci()
        .args(["pair-pair", "--type", "neither", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout("");
    gloci()
        .args(["pair-pair", "--type", "bogus", "-a", &a, "-b", &b])
        .assert()
        .code(2);
}

#[test]
fn pair_pair_window_with_threshold_is_usage() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    gloci()
        .args(["pair-pair", "-w", "10", "-f", "0.5", "-a", &a, "-b", &b])
        .assert()
        .code(2);
    gloci()
        .args(["pair-pair", "--split", "-a", &a, "-b", &b])
        .assert()
        .code(2)
        .stderr(contains("--split"));
    gloci()
        .args(["pair-pair", "--help"])
        .assert()
        .success()
        .stdout(contains("--split").not());
}

#[test]
fn pair_pair_names_need_names() {
    let dir = tempdir();
    let a = plain_file(dir.path(), "a.bedpe", "chr1\t10\t20\tchr1\t100\t110\n");
    let b = plain_file(dir.path(), "b.bedpe", B);
    gloci()
        .args(["pair-pair", "--require-different-names", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .assert()
        .code(1)
        .stderr(contains("a.bedpe line 1"));
    let (a, b) = files(dir.path());
    gloci()
        .args(["pair-pair", "--require-different-names", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout(format!("{}\t{}", A.trim_end(), B));
}

#[test]
fn pair_pair_unknown_strand_warns_once() {
    let dir = tempdir();
    let a = plain_file(
        dir.path(),
        "a.bedpe",
        "chr1\t10\t20\tchr1\t100\t110\tp\t0\t.\t.\nchr1\t10\t20\tchr1\t100\t110\tq\t0\t.\t.\n",
    );
    let b = plain_file(dir.path(), "b.bedpe", B);
    let out = gloci()
        .args(["pair-pair", "-s", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .unwrap();
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(err.matches("unknown strand").count(), 1, "{err}");
}

#[test]
fn pair_pair_b_must_be_bedpe() {
    let dir = tempdir();
    let a = plain_file(dir.path(), "a.bedpe", A);
    let b = plain_file(dir.path(), "b.bed", "chr1\t12\t14\tb1\n");
    gloci()
        .args(["pair-pair", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .assert()
        .code(1);
}

#[test]
fn pair_pair_output_formats() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    let want = format!("{}\t{}", A.trim_end(), B);
    for name in ["o.tsv", "o.bedpe"] {
        let o = dir.path().join(name);
        gloci()
            .args(["pair-pair", "-a", &a, "-b", &b, "-o"])
            .arg(&o)
            .assert()
            .success();
        assert_eq!(std::fs::read_to_string(&o).unwrap(), want);
    }
    gloci()
        .args(["pair-pair", "-a", &a, "-b", &b, "--write-index", "-o"])
        .arg(dir.path().join("o.bedpe.gz"))
        .assert()
        .code(2)
        .stderr(contains("gloci index -s 1 -b 2 -e 3 --zero-based"));
}
