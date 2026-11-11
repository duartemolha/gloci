//! `gloci flagstat` (D128): layouts, format inference, CRAM, thread counts; the byte-level
//! comparison with samtools is in the goldens (`flagstat-*-vs-samtools`).
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success().get_output().clone();
    String::from_utf8(o.stdout).unwrap()
}

fn fixture() -> String {
    s(&inputs().join("samtools/flagstat.bam"))
}

#[test]
fn flagstat_text_has_samtools_percentages() {
    let t = out(&["flagstat", &fixture()]);
    assert!(t.contains("306 + 10 mapped (95.63% : 83.33%)\n"), "{t}");
    assert!(t.contains("7 + 0 singletons (4.37% : N/A)\n"), "{t}");
    assert!(t.contains("40 + 1 secondary\n"), "{t}");
    assert!(
        t.contains("22 + 0 with mate mapped to a different chr (mapQ>=5)\n"),
        "{t}"
    );
}

#[test]
fn flagstat_format_inference() {
    let d = tempdir();
    let bam = fixture();
    let json = d.path().join("s.json");
    gloci()
        .args(["flagstat", &bam, "-o", &s(&json)])
        .assert()
        .success()
        .stdout("");
    let j = std::fs::read_to_string(&json).unwrap();
    assert!(j.starts_with("{\n \"QC-passed reads\": { \n"), "{j}");
    // -O wins over the extension.
    let tsv = d.path().join("t.json");
    gloci()
        .args(["flagstat", &bam, "-O", "tsv", "-o", &s(&tsv)])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&tsv)
            .unwrap()
            .starts_with("320\t12\ttotal (QC-passed reads + QC-failed reads)\n")
    );
    // .tsv infers TSV; anything else is text.
    let t = d.path().join("u.tsv");
    gloci()
        .args(["flagstat", &bam, "-o", &s(&t)])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&t)
            .unwrap()
            .contains("\tmapped %\n")
    );
    let txt = d.path().join("v.txt");
    gloci()
        .args(["flagstat", &bam, "-o", &s(&txt)])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&txt)
            .unwrap()
            .starts_with("320 + 12 in total")
    );
    // An unknown -O is a usage error that lists the names.
    gloci()
        .args(["flagstat", &bam, "-O", "csv"])
        .assert()
        .code(2)
        .stdout("")
        .stderr(contains("text, tsv, json"));
}

#[test]
fn flagstat_cram_equals_bam() {
    let cram = s(&inputs().join("aln/flags.cram"));
    let reference = s(&inputs().join("fa/ref.fa"));
    let from_cram = out(&["flagstat", &cram, "--reference", &reference]);
    let from_sam = out(&["flagstat", &s(&inputs().join("aln/flags.sam"))]);
    assert_eq!(from_cram, from_sam);
    assert!(from_cram.starts_with("13 + 1 in total"), "{from_cram}");
    gloci()
        .args(["flagstat", &cram])
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains("--reference"));
}

#[test]
fn flagstat_sam_equals_bam() {
    let sam = out(&["flagstat", &s(&inputs().join("samtools/flagstat.sam"))]);
    assert_eq!(sam, out(&["flagstat", &fixture()]));
}

#[test]
fn flagstat_is_identical_at_any_thread_count() {
    let big = s(&inputs().join("aln/cat-big.bam"));
    let cram = s(&inputs().join("aln/many.cram"));
    let reference = s(&inputs().join("fa/ref.fa"));
    let one_bam = fixture();
    for args in [
        vec!["flagstat", big.as_str()],
        vec!["flagstat", cram.as_str(), "--reference", reference.as_str()],
        vec!["flagstat", one_bam.as_str()],
    ] {
        let mut one = args.clone();
        one.extend(["-t", "1"]);
        let want = out(&one);
        for t in ["2", "8"] {
            let mut a = args.clone();
            a.extend(["-t", t]);
            assert_eq!(out(&a), want, "{args:?} -t {t}");
        }
    }
}

#[test]
fn flagstat_stdin_and_missing_file() {
    let sam = std::fs::read(inputs().join("samtools/flagstat.sam")).unwrap();
    gloci()
        .args(["flagstat", "-"])
        .write_stdin(sam)
        .assert()
        .success()
        .stdout(contains("320 + 12 in total"));
    gloci()
        .args(["flagstat", "/nonexistent/x.bam"])
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains("cannot open /nonexistent/x.bam"));
}
