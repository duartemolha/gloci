//! `gloci sort` of BAM and SAM (D121): usage errors, `--check`, the output format, the index,
//! thread counts and the `@PG` line. The orders themselves are pinned against samtools by the
//! `aln-sort-*` golden cases.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn view_sam(bam: &Path) -> String {
    let o = gloci()
        .args(["view", &s(bam), "-O", "sam"])
        .assert()
        .success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn names(sam: &str) -> Vec<&str> {
    sam.lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').next().unwrap())
        .collect()
}

#[test]
fn aln_sort_cli_errors() {
    let dir = tempdir();
    let bam = s(&fixture("sort.bam"));
    let bed = plain_file(dir.path(), "x.bed", "chr1\t1\t2\n");
    let two = |args: &[&str]| gloci().arg("sort").args(args).assert().code(2);
    two(&[&bam, &s(&fixture("flags.bam"))]).stderr(contains("sort -m").and(contains("gloci cat")));
    two(&["-n", &s(&bed)]).stderr(contains("-n and -N sort alignments"));
    two(&["-N", &s(&bed)]).stderr(contains("-n and -N sort alignments"));
    two(&[
        "-n",
        &bam,
        "-o",
        &s(&dir.path().join("o.bam")),
        "--write-index",
    ])
    .stderr(contains("--write-index needs coordinate order"));
    two(&["-g", "hg38", &bam]).stderr(contains("-g does not apply to alignments"));
    two(&["--full-key", &bam]).stderr(contains("--full-key"));
    two(&["--unique", &bam]).stderr(contains("--unique"));
    two(&[&bam, "-o", &s(&dir.path().join("x.cram"))]).stderr(contains("does not write CRAM"));
    // -n and -N are two orders, not both
    two(&["-n", "-N", &bam]);
}

use predicates::prelude::PredicateBooleanExt;

#[test]
fn aln_sort_check() {
    let sorted = tempdir();
    let out = sorted.path().join("s.bam");
    gloci()
        .args(["sort", &s(&fixture("sort.bam")), "-o", &s(&out)])
        .assert()
        .success();
    gloci().args(["sort", "--check", &s(&out)]).assert().code(0);
    gloci()
        .args(["sort", "--check", &s(&fixture("sort.bam"))])
        .assert()
        .code(1)
        .stderr(contains("record 2 (read10)").and(contains("record 1 (read20)")))
        .stderr(contains("-o sorted.bam"));
    // The hint's output follows the input's format (P8).
    gloci()
        .args(["sort", "--check", &s(&fixture("sort.sam"))])
        .assert()
        .code(1)
        .stderr(contains("-o sorted.sam"));
    // -n: the coordinate-sorted file is not in name order, and sorted by name it is
    gloci()
        .args(["sort", "--check", "-n", &s(&out)])
        .assert()
        .code(1);
    let by_name = sorted.path().join("n.bam");
    gloci()
        .args(["sort", "-n", &s(&out), "-o", &s(&by_name)])
        .assert()
        .success();
    gloci()
        .args(["sort", "--check", "-n", &s(&by_name)])
        .assert()
        .code(0);
    gloci()
        .args(["sort", "--check", &s(&by_name)])
        .assert()
        .code(1);
}

#[test]
fn aln_sort_output_formats_and_index() {
    let dir = tempdir();
    // the format follows the input (SAM stays SAM) unless -o or -O say otherwise
    let bam = dir.path().join("o.bam");
    gloci()
        .args([
            "sort",
            &s(&fixture("sort.sam")),
            "-o",
            &s(&bam),
            "--write-index",
        ])
        .assert()
        .success();
    assert!(dir.path().join("o.bam.bai").exists());
    let sam_out = dir.path().join("o.sam");
    gloci()
        .args(["sort", &s(&fixture("sort.bam")), "-o", &s(&sam_out)])
        .assert()
        .success();
    let text = std::fs::read_to_string(&sam_out).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"), "{text}");
    assert_eq!(names(&text), names(&view_sam(&bam)));
    assert!(text.contains("@PG\tID:gloci"));
    // samtools reads the index through a region
    gloci()
        .args(["view", "-O", "sam", "--no-header", &s(&bam), "chr1:99-99"])
        .assert()
        .success()
        .stdout(contains("read30"));
    let o = gloci()
        .args(["sort", &s(&fixture("sort.sam"))])
        .assert()
        .success();
    assert!(
        o.get_output()
            .stdout
            .starts_with(b"@HD\tVN:1.6\tSO:coordinate\n")
    );
    // SAM output cannot be indexed
    gloci()
        .args([
            "sort",
            &s(&fixture("sort.bam")),
            "-O",
            "sam",
            "--write-index",
        ])
        .assert()
        .failure();
    // stdin
    let o = gloci()
        .args(["sort", "-O", "sam", "-"])
        .write_stdin(std::fs::read(fixture("sort.bam")).unwrap())
        .assert()
        .success();
    assert_eq!(
        String::from_utf8(o.get_output().stdout.clone())
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with("@PG"))
            .collect::<Vec<_>>(),
        text.lines()
            .filter(|l| !l.starts_with("@PG"))
            .collect::<Vec<_>>()
    );
}

#[test]
fn aln_sort_threads_and_program_line() {
    let dir = tempdir();
    let mut outs = Vec::new();
    for t in ["1", "2", "8"] {
        let out = dir.path().join(format!("t{t}.bam"));
        gloci()
            .args([
                "sort",
                "-n",
                &s(&fixture("sort.bam")),
                "-o",
                &s(&out),
                "-t",
                t,
            ])
            .assert()
            .success();
        let sam = view_sam(&out);
        let pg = sam
            .lines()
            .find(|l| l.starts_with("@PG\tID:gloci"))
            .unwrap()
            .to_string();
        assert!(
            pg.contains("CL:") && pg.contains(" sort -n ") && !pg.contains(" -t"),
            "{pg}"
        );
        outs.push(
            sam.lines()
                .filter(|l| !l.starts_with("@PG"))
                .map(String::from)
                .collect::<Vec<_>>(),
        );
    }
    assert_eq!(outs[0], outs[1]);
    assert_eq!(outs[0], outs[2]);
}

#[test]
fn aln_sort_in_place() {
    let dir = tempdir();
    let p = dir.path().join("x.bam");
    std::fs::copy(fixture("sort.bam"), &p).unwrap();
    gloci()
        .args(["sort", &s(&p), "-o", &s(&p)])
        .assert()
        .success();
    gloci().args(["sort", "--check", &s(&p)]).assert().code(0);
}
