//! Regression tests for the final-review fix wave (groups A-G).
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::Path;

fn w(d: &Path, name: &str, text: &str) -> String {
    plain_file(d, name, text).to_string_lossy().into_owned()
}

fn p(d: &Path, name: &str) -> String {
    d.join(name).to_string_lossy().into_owned()
}

// ---- Group A: sort/view format handling ----

const BG: &str = "chr1\t10\t20\t0.5\nchr1\t30\t40\t1.5\n";

#[test]
fn trackless_bedgraph_to_bedgraph_output() {
    let d = tempdir();
    let i = w(d.path(), "sig.bedgraph", BG);
    for cmd in ["sort", "view"] {
        let o = p(d.path(), &format!("{cmd}.out.bedgraph"));
        gloci().args([cmd, &i, "-o", &o]).assert().success();
        assert_eq!(std::fs::read_to_string(&o).unwrap(), BG, "{cmd}");
    }
}

#[test]
fn bed_input_to_bedgraph_output_and_back() {
    let d = tempdir();
    let i = w(d.path(), "x.bed", "chr1\t1\t2\n");
    let o = p(d.path(), "x.bedgraph");
    gloci().args(["sort", &i, "-o", &o]).assert().success();
    let bg = w(d.path(), "y.bedgraph", "chr1\t1\t2\t3\n");
    let o2 = p(d.path(), "y.bed");
    gloci().args(["view", &bg, "-o", &o2]).assert().success();
}

#[test]
fn sort_accepts_format_flag() {
    let d = tempdir();
    let i = w(
        d.path(),
        "sig.txt",
        "track type=bedGraph\nchr1\t30\t40\t1\nchr1\t10\t20\t2\n",
    );
    gloci()
        .args(["sort", &i, "--format", "bedgraph", "--header"])
        .assert()
        .success()
        .stdout("track type=bedGraph\nchr1\t10\t20\t2\nchr1\t30\t40\t1\n");
    gloci()
        .args(["sort", &i, "--format", "nope"])
        .assert()
        .code(2);
}

#[test]
fn family_mismatch_message_is_identical_for_sort_and_view() {
    let d = tempdir();
    let v = w(
        d.path(),
        "a.vcf",
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    );
    let o = p(d.path(), "x.bed");
    let msg = |cmd: &str| {
        let out = gloci()
            .args([cmd, &v, "-o", &o])
            .assert()
            .code(2)
            .get_output()
            .stderr
            .clone();
        String::from_utf8(out).unwrap()
    };
    let (s, vw) = (msg("sort"), msg("view"));
    assert_eq!(s, vw);
    assert!(s.contains("the input is VCF but the output format is BED"), "{s}");
    assert!(s.contains("-O"), "{s}");
    assert!(!Path::new(&o).exists());
}
