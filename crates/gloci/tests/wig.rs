//! WIG input (D63): fixedStep and variableStep read as bedGraph records.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// `tests/wig-golden/input.wig` and its expected bedGraph. UCSC
/// `wigToBigWig`/`bigWigToBedGraph` are not in the reference image, so the
/// expected file is checked by hand. WIG positions are 1-based and a data
/// line covers `span` bases, so 0-based start = position - 1 and
/// end = start + span:
/// - `fixedStep chrom=chr1 start=101 step=10 span=5`: lines at 1-based
///   101, 111, 121 give [100,105), [110,115), [120,125) with values
///   1.5, 2, 0.25.
/// - `variableStep chrom=chr2 span=3`: `11 7` gives [10,13) and `20 -1.5`
///   gives [19,22).
/// - `variableStep chrom=chr2` (span 1): `30 4` gives [29,30).
///
/// Values are copied as written; the `track` and `browser` lines are skipped.
fn golden() -> (String, String) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/wig-golden");
    (
        s(&root.join("input.wig")),
        std::fs::read_to_string(root.join("expected.bedgraph")).unwrap(),
    )
}

#[test]
fn wig_fixed_and_variable() {
    let (input, expected) = golden();
    assert_eq!(out(&["convert", &input]), expected);
}

#[test]
fn convert_wig_to_bedgraph() {
    let (input, expected) = golden();
    let d = tempdir();
    let o = s(&d.path().join("o.bedgraph"));
    gloci()
        .args(["convert", &input, "-o", &o])
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&o).unwrap(), expected);
    // Equal neighbouring values are not merged.
    let same = s(&plain_file(
        d.path(),
        "same.wig",
        "variableStep chrom=c\n1 5\n2 5\n",
    ));
    assert_eq!(out(&["convert", &same]), "c\t0\t1\t5\nc\t1\t2\t5\n");
}

#[test]
fn wig_output_is_identical_for_every_thread_count() {
    let (input, expected) = golden();
    for t in ["1", "2", "4"] {
        assert_eq!(out(&["convert", &input, "--threads", t]), expected);
    }
}

#[test]
fn wig_as_b_in_map() {
    let d = tempdir();
    let a = s(&plain_file(
        d.path(),
        "a.bed",
        "chr1\t100\t125\nchr1\t200\t300\n",
    ));
    let (b, _) = golden();
    assert_eq!(
        out(&["map", "-a", &a, "-b", &b, "-c", "4:mean"]),
        "chr1\t100\t125\t1.25\nchr1\t200\t300\t.\n"
    );
}

#[test]
fn wig_with_chrom_map() {
    let d = tempdir();
    let wig = s(&plain_file(
        d.path(),
        "x.wig",
        "fixedStep chrom=1 start=1 step=1\n3\n4\nvariableStep chrom=2\n5 6\n",
    ));
    let map = s(&plain_file(d.path(), "m.tsv", "1\tchr1\n2\tchr2\n"));
    assert_eq!(
        out(&["convert", &wig, "--chrom-map", &map]),
        "chr1\t0\t1\t3\nchr1\t1\t2\t4\nchr2\t4\t5\t6\n"
    );
}

#[test]
fn wig_format_flag_and_crlf() {
    let d = tempdir();
    let wig = s(&plain_file(
        d.path(),
        "x.txt",
        "fixedStep chrom=c start=5\r\n9\r\n",
    ));
    assert_eq!(out(&["convert", &wig, "--format", "wig"]), "c\t4\t5\t9\n");
}

#[test]
fn wig_errors() {
    let d = tempdir();
    let fail = |text: &str, want: &str| {
        let f = s(&plain_file(d.path(), "bad.wig", text));
        gloci()
            .args(["convert", &f, "--format", "wig"])
            .assert()
            .code(1)
            .stderr(contains(want));
    };
    fail(
        "track type=wiggle_0\n1.5\nfixedStep chrom=c start=1 step=1\n",
        "line 2",
    );
    fail("1.5\n", "before any fixedStep or variableStep line");
    fail("fixedStep start=1 step=1\n1\n", "no chrom=NAME");
    fail("fixedStep chrom=c step=1\n1\n", "no start=POSITION");
    fail("fixedStep chrom=c start=0 step=1\n1\n", "1-based");
    fail(
        "variableStep chrom=c span=0\n1 1\n",
        "span must be at least 1",
    );
    fail("variableStep chrom=c\n0 1\n", "1-based");
    fail("variableStep chrom=c\n5 abc\n", "not a number");
    fail("variableStep chrom=c\n5\n", "POSITION VALUE");
    fail("fixedStep chrom=c start=1 step=1\n1 2\n", "one value");
}

#[test]
fn wig_skip_invalid_keeps_positions() {
    let d = tempdir();
    let f = s(&plain_file(
        d.path(),
        "x.wig",
        "fixedStep chrom=c start=1 step=2\n1\nx\n3\n",
    ));
    assert_eq!(
        out(&["convert", &f, "--skip-invalid"]),
        "c\t0\t1\t1\nc\t4\t5\t3\n"
    );
}

#[test]
fn wig_cannot_be_indexed() {
    let (input, _) = golden();
    gloci()
        .args(["index", &input])
        .assert()
        .failure()
        .stderr(contains("WIG"));
}
