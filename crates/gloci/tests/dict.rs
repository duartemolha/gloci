//! `gloci dict` (D132): the command's own rules beyond the samtools goldens.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn inputs(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs")
        .join(name)
}

#[test]
fn dict_round_trip() {
    let d = tempdir();
    let dict = d.path().join("ref.dict");
    gloci()
        .args(["dict"])
        .arg(inputs("fa/ref.fa"))
        .arg("-o")
        .arg(&dict)
        .assert()
        .success()
        .stdout("");
    let text = std::fs::read_to_string(&dict).unwrap();
    assert!(text.starts_with("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chrM\tLN:150\tM5:"));
    assert!(text.contains("\tUR:file://"));
    // Lines of the three chromosomes the dictionary and the index share.
    let bed = d.path().join("in.bed");
    std::fs::write(
        &bed,
        "chr10\t5\t9\nchrM\t1\t2\nchr2\t3\t4\nchr1\t7\t8\nchrX\t1\t2\n",
    )
    .unwrap();
    let sorted = |g: &Path| {
        gloci()
            .args(["sort"])
            .arg(&bed)
            .arg("-g")
            .arg(g)
            .output()
            .unwrap()
    };
    let by_dict = sorted(&dict);
    let by_fai = sorted(&inputs("fa/ref.fa.fai"));
    assert!(by_dict.status.success(), "{:?}", by_dict);
    assert_eq!(by_dict.stdout, by_fai.stdout);
    assert!(String::from_utf8_lossy(&by_dict.stdout).starts_with("chrM\t1\t2\n"));
}

#[test]
fn dict_duplicate_names_exit_1() {
    let d = tempdir();
    let out = d.path().join("o.dict");
    gloci()
        .args(["dict", "--uri", "x"])
        .arg(inputs("samtools/dup.fa"))
        .arg("-o")
        .arg(&out)
        .assert()
        .code(1)
        .stderr(
            contains("line 5")
                .and(contains("line 1"))
                .and(contains("'a'")),
        );
    assert!(!out.exists(), "no dictionary is left behind on error");
}

#[test]
fn dict_reads_stdin_without_uri() {
    gloci()
        .args(["dict", "-"])
        .write_stdin(">a desc\nAC GT\r\n")
        .assert()
        .success()
        .stdout("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:a\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\n");
}

#[test]
fn dict_missing_alt_file_names_it() {
    gloci()
        .args(["dict"])
        .arg(inputs("fa/ref.fa"))
        .args(["--alt", "/nonexistent/alt.txt"])
        .assert()
        .code(1)
        .stderr(contains("/nonexistent/alt.txt"));
}

/// Blank lines before the first header are skipped, as samtools does (M9); a non-FASTA input
/// is refused in dict's own words.
#[test]
fn dict_skips_leading_blank_lines_and_names_its_input() {
    let d = tempdir();
    let fa = plain_file(d.path(), "blank.fa", "\n\n>a\nACGT\n");
    let o = gloci()
        .args(["dict", "--no-header"])
        .arg(&fa)
        .args(["--uri", "x"])
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "@SQ\tSN:a\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\tUR:x\n"
    );
    let crlf = plain_file(d.path(), "crlf.fa", "\r\n>a\r\nACGT\r\n");
    gloci()
        .args(["dict", "--no-header"])
        .arg(&crlf)
        .args(["--uri", "x"])
        .assert()
        .success()
        .stdout("@SQ\tSN:a\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\tUR:x\n");
    let bam = inputs("samtools/depth.bam");
    gloci()
        .arg("dict")
        .arg(&bam)
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains(format!(
            "{} is BAM, not FASTA; dict makes a sequence dictionary from a FASTA file",
            bam.display()
        )))
        .stderr(contains("intervals").not());
}
