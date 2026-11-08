mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

#[test]
fn tsv_index_is_written() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.tsv.gz", "chr1\t10\t20\tx\nchr1\t30\t40\ty\n");
    gloci()
        .args(["index", "-s", "1", "-b", "2", "-e", "3", "--zero-based"])
        .arg(&p)
        .assert()
        .success()
        .stdout("");
    assert!(d.path().join("a.tsv.gz.tbi").exists());
}

#[test]
fn unsorted_tsv_names_line_3() {
    let d = tempdir();
    let p = bgzf_file(
        d.path(),
        "u.tsv.gz",
        "chr1\t10\t20\nchr1\t50\t60\nchr1\t30\t40\n",
    );
    gloci()
        .args(["index", "-s", "1", "-b", "2", "-e", "3", "--zero-based"])
        .arg(&p)
        .assert()
        .code(1)
        .stderr(contains("error: ").and(contains("line 3")));
}

#[test]
fn seq_without_begin_is_a_usage_error() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.tsv.gz", "chr1\t10\t20\n");
    gloci()
        .args(["index", "-s", "1"])
        .arg(&p)
        .assert()
        .code(2)
        .stderr(contains("-s").and(contains("-b")));
}

#[test]
fn zero_based_without_columns_is_a_usage_error() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.bed.gz", "chr1\t10\t20\n");
    gloci()
        .args(["index", "--zero-based"])
        .arg(&p)
        .assert()
        .code(2)
        .stderr(contains("--zero-based"));
}

#[test]
fn skip_lines_and_comment_need_columns() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.bed.gz", "chr1\t10\t20\n");
    for flag in [["--skip-lines", "1"], ["--comment", "%"]] {
        gloci().arg("index").args(flag).arg(&p).assert().code(2);
    }
}

#[test]
fn comment_must_be_one_character() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.tsv.gz", "chr1\t10\t20\n");
    gloci()
        .args(["index", "-s", "1", "-b", "2", "--comment", "ab"])
        .arg(&p)
        .assert()
        .code(2);
}

#[test]
fn tsv_mode_non_bgzf_keeps_recompress_hint() {
    let d = tempdir();
    let p = plain_file(d.path(), "a.tsv", "chr1\t10\t20\n");
    gloci()
        .args(["index", "-s", "1", "-b", "2", "-e", "3"])
        .arg(&p)
        .assert()
        .code(2)
        .stderr(contains("recompress with: gloci view"));
}

#[test]
fn column_zero_is_a_usage_error() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.tsv.gz", "chr1\t10\t20\n");
    for args in [["-s", "0", "-b", "2"], ["-s", "1", "-b", "0"]] {
        gloci()
            .arg("index")
            .args(args)
            .arg(&p)
            .assert()
            .code(2)
            .stderr(contains("column is numbered from 1"));
    }
    gloci()
        .args(["index", "-s", "1", "-b", "2", "-e", "0"])
        .arg(&p)
        .assert()
        .code(2)
        .stderr(contains("end column is numbered from 1"));
}
