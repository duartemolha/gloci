mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

const SORTED_BED: &str = "chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t6\n";
const UNSORTED_BED: &str = "chr1\t30\t40\nchr1\t10\t20\n";

#[test]
fn index_writes_tbi() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "a.bed.gz", SORTED_BED);
    gloci().arg("index").arg(&p).assert().success().stdout("");
    assert!(d.path().join("a.bed.gz.tbi").exists());
}

#[test]
fn index_plain_file_explains() {
    let d = tempdir();
    let p = plain_file(d.path(), "a.bed", SORTED_BED);
    gloci().arg("index").arg(&p).assert().code(2).stderr(
        contains("error: ").and(contains("not BGZF-compressed; recompress with: gloci view")),
    );
}

#[test]
fn index_unsorted_names_lines() {
    let d = tempdir();
    let p = bgzf_file(d.path(), "u.bed.gz", UNSORTED_BED);
    gloci()
        .arg("index")
        .arg(&p)
        .assert()
        .code(1)
        .stderr(contains("error: ").and(contains("is not sorted: line 2")));
}
