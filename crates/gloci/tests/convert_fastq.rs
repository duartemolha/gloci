//! `gloci convert --fastq` (D166): reads as FASTQ, in pairs with `--fastq2`. Byte parity with
//! bedtools is in the golden cases `bamtofastq-*`; here are the files, compression, filters,
//! formats and errors those cannot show.
mod common;
use common::*;
use predicates::str::contains;
use std::io::Read;
use std::path::{Path, PathBuf};

fn aln(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn golden(case: &str) -> Vec<u8> {
    std::fs::read(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/golden")
            .join(case)
            .join("expected"),
    )
    .unwrap()
}

fn stdout(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

/// The records of a FASTQ text whose names end in `/1` and in `/2`, as two texts.
fn split_ends(fastq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let text = String::from_utf8(fastq.to_vec()).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let (mut a, mut b) = (String::new(), String::new());
    for rec in lines.chunks(4) {
        let dst = if rec[0].ends_with("/1") {
            &mut a
        } else {
            &mut b
        };
        for l in rec {
            dst.push_str(l);
            dst.push('\n');
        }
    }
    (a.into_bytes(), b.into_bytes())
}

fn gunzip(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut out)
        .unwrap();
    out
}

#[test]
fn two_files_hold_the_ends_of_the_golden_interleave() {
    let d = tempdir();
    let (r1, r2) = (d.path().join("r1.fq"), d.path().join("r2.fq"));
    gloci()
        .args(["convert", &aln("fastq-pairs.bam"), "--fastq", "-o"])
        .arg(&r1)
        .arg("--fastq2")
        .arg(&r2)
        .assert()
        .success()
        .stdout("");
    let (want1, want2) = split_ends(&golden("bamtofastq-fq2-vs-bedtools"));
    assert!(!want1.is_empty());
    assert_eq!(std::fs::read(&r1).unwrap(), want1);
    assert_eq!(std::fs::read(&r2).unwrap(), want2);
}

#[test]
fn tags_files_hold_the_ends_of_the_golden_interleave() {
    let d = tempdir();
    let (r1, r2) = (d.path().join("r1.fq"), d.path().join("r2.fq"));
    gloci()
        .args(["convert", &aln("fastq-tags.bam"), "--fastq", "--tags", "-o"])
        .arg(&r1)
        .arg("--fastq2")
        .arg(&r2)
        .assert()
        .success();
    let (want1, want2) = split_ends(&golden("bamtofastq-tags-vs-bedtools"));
    assert_eq!(std::fs::read(&r1).unwrap(), want1);
    assert_eq!(std::fs::read(&r2).unwrap(), want2);
}

#[test]
fn gz_names_are_bgzf_for_both_files_and_stdout_is_text() {
    let d = tempdir();
    let (r1, r2) = (d.path().join("r1.fq.gz"), d.path().join("r2.fastq.gz"));
    gloci()
        .args(["convert", &aln("fastq-pairs.bam"), "--fastq", "-o"])
        .arg(&r1)
        .arg("--fastq2")
        .arg(&r2)
        .assert()
        .success();
    let (want1, want2) = split_ends(&golden("bamtofastq-fq2-vs-bedtools"));
    assert_eq!(gunzip(&r1), want1);
    assert_eq!(gunzip(&r2), want2);
    // stdout is text, whatever -t says
    let text = stdout(&["convert", &aln("fastq-reads.bam"), "--fastq"]);
    assert!(text.starts_with(b"@r1\nACGTACGT\n+\nIIIIHHGG\n"));
}

#[test]
fn first_ends_to_stdout_when_only_fastq2_is_a_file() {
    let d = tempdir();
    let r2 = d.path().join("r2.fq");
    let out = gloci()
        .args(["convert", &aln("fastq-pairs.bam"), "--fastq", "--fastq2"])
        .arg(&r2)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let (want1, want2) = split_ends(&golden("bamtofastq-fq2-vs-bedtools"));
    assert_eq!(out, want1);
    assert_eq!(std::fs::read(&r2).unwrap(), want2);
}

#[test]
fn the_same_file_for_both_ends_interleaves_them() {
    let d = tempdir();
    let both = d.path().join("both.fq");
    gloci()
        .args(["convert", &aln("fastq-pairs.bam"), "--fastq", "-o"])
        .arg(&both)
        .arg("--fastq2")
        .arg(&both)
        .assert()
        .success();
    assert_eq!(
        std::fs::read(&both).unwrap(),
        golden("bamtofastq-fq2-vs-bedtools")
    );
    // Stdout as a BGZF stream for one end and as text for the other is not one output.
    gloci()
        .args([
            "convert",
            &aln("fastq-pairs.bam"),
            "--fastq",
            "-O",
            "fq.gz",
            "--fastq2",
            "-",
        ])
        .assert()
        .code(2)
        .stderr(contains("different compression"));
}

#[test]
fn sam_bam_and_cram_agree() {
    let want = golden("bamtofastq-fq2-vs-bedtools");
    for f in [
        "fastq-pairs.sam",
        "fastq-pairs.bam",
        "fastq-pairs.noref.cram",
    ] {
        let got = stdout(&[
            "convert",
            &aln(f),
            "--fastq",
            "-o",
            "/dev/stdout",
            "--fastq2",
            "/dev/stdout",
        ]);
        assert_eq!(got, want, "{f}");
    }
}

#[test]
fn output_is_the_same_at_every_thread_count() {
    let reads = aln("fastq-reads.bam");
    let pairs = aln("fastq-pairs.bam");
    let tags = aln("fastq-tags.bam");
    let runs: [Vec<&str>; 3] = [
        vec!["--fastq", &reads],
        vec![
            "--fastq",
            &pairs,
            "-o",
            "/dev/stdout",
            "--fastq2",
            "/dev/stdout",
        ],
        vec![
            "--fastq",
            &tags,
            "--tags",
            "-o",
            "/dev/stdout",
            "--fastq2",
            "/dev/stdout",
        ],
    ];
    for args in runs {
        let run = |t: &str| {
            let mut v = vec!["convert", "-t", t];
            v.extend_from_slice(&args);
            stdout(&v)
        };
        let one = run("1");
        assert!(!one.is_empty());
        assert_eq!(run("2"), one, "{args:?}");
        assert_eq!(run("8"), one, "{args:?}");
    }
}

#[test]
fn exclude_flags_drop_records_before_pairing() {
    let all = stdout(&["convert", &aln("fastq-reads.bam"), "--fastq"]);
    let text = String::from_utf8(all).unwrap();
    assert!(text.contains("@r4\n") && text.contains("@r5\n") && text.contains("@r3\n"));
    // 0x900 is secondary and supplementary; unmapped (0x4), QC fail and duplicates stay.
    let kept = String::from_utf8(stdout(&[
        "convert",
        &aln("fastq-reads.bam"),
        "--fastq",
        "--exclude-flags",
        "0x900",
    ]))
    .unwrap();
    assert!(!kept.contains("@r4\n") && !kept.contains("@r5\n"));
    for name in ["@r3\n", "@r6\n", "@r7\n"] {
        assert!(kept.contains(name), "{name}");
    }
}

#[test]
fn a_third_record_is_an_error_that_the_filter_cures() {
    let d = tempdir();
    let r1 = d.path().join("r1.fq");
    let r2 = d.path().join("r2.fq");
    gloci()
        .args(["convert", &aln("fastq-three.bam"), "--fastq", "-o"])
        .arg(&r1)
        .arg("--fastq2")
        .arg(&r2)
        .assert()
        .code(1)
        .stderr(contains("the read name p2 has more than two records"))
        .stderr(contains("--exclude-flags 0x900"))
        .stderr(contains("record 5 of chr1"));
    // Nothing is left behind.
    assert!(!r1.exists() && !r2.exists());
    gloci()
        .args(["convert", &aln("fastq-three.bam"), "--fastq", "-o"])
        .arg(&r1)
        .arg("--fastq2")
        .arg(&r2)
        .args(["--exclude-flags", "0x900"])
        .assert()
        .success();
    let text = String::from_utf8(std::fs::read(&r1).unwrap()).unwrap();
    assert_eq!(text.matches("/1\n").count(), 2);
}

#[test]
fn a_coordinate_sorted_file_is_a_usage_error_with_fastq2() {
    gloci()
        .args([
            "convert",
            &aln("cigars.bam"),
            "--fastq",
            "--fastq2",
            "/dev/null",
        ])
        .assert()
        .code(2)
        .stderr(contains("sorted by coordinate"))
        .stderr(contains("gloci sort -n"));
    // Without --fastq2 the order does not matter.
    gloci()
        .args(["convert", &aln("cigars.bam"), "--fastq"])
        .assert()
        .success();
}

#[test]
fn tags_need_the_tags_and_fastq2() {
    gloci()
        .args(["convert", &aln("fastq-notags.bam"), "--fastq", "--tags"])
        .assert()
        .code(2)
        .stderr(contains("--fastq2"));
    gloci()
        .args([
            "convert",
            &aln("fastq-notags.bam"),
            "--fastq",
            "--tags",
            "-o",
            "/dev/null",
            "--fastq2",
            "/dev/null",
        ])
        .assert()
        .code(1)
        .stderr(contains("read p1 has no R2 string tag"))
        .stderr(contains("drop --tags"));
}

#[test]
fn fastq_is_for_alignments_and_alone() {
    gloci()
        .args(["convert", "--fastq", &aln("targets.bed")])
        .assert()
        .code(2)
        .stderr(contains("is BED, but this command reads alignments"));
    for other in [
        "--bedpe",
        "--bed12",
        "--blocks",
        "--cigar",
        "--split-alleles",
    ] {
        gloci()
            .args(["convert", &aln("fastq-reads.bam"), "--fastq", other])
            .assert()
            .code(2)
            .stderr(contains("cannot be used with"));
    }
    gloci()
        .args(["convert", &aln("fastq-reads.bam"), "--fastq2", "x.fq"])
        .assert()
        .code(2)
        .stderr(contains("--fastq"));
}

#[test]
fn the_output_must_be_fastq() {
    let d = tempdir();
    gloci()
        .args(["convert", &aln("fastq-reads.bam"), "--fastq", "-o"])
        .arg(d.path().join("x.bed"))
        .assert()
        .code(2)
        .stderr(contains("convert --fastq writes FASTQ, not BED"));
    gloci()
        .args([
            "convert",
            &aln("fastq-reads.bam"),
            "--fastq",
            "--write-index",
            "-o",
        ])
        .arg(d.path().join("x.fq.gz"))
        .assert()
        .code(2)
        .stderr(contains("cannot index FASTQ"));
    assert!(!d.path().join("x.bed").exists());
}

#[test]
fn an_input_that_is_not_readable_creates_no_output() {
    let d = tempdir();
    let (r1, r2) = (d.path().join("r1.fq"), d.path().join("r2.fq"));
    gloci()
        .args(["convert", "/no/such.bam", "--fastq", "-o"])
        .arg(&r1)
        .arg("--fastq2")
        .arg(&r2)
        .assert()
        .failure();
    assert!(!r1.exists() && !r2.exists());
}

#[test]
fn reads_are_written_once_in_input_order() {
    let want = String::from_utf8(golden("bamtofastq-vs-bedtools")).unwrap();
    let got = String::from_utf8(stdout(&["convert", &aln("fastq-reads.bam"), "--fastq"])).unwrap();
    assert_eq!(got, want);
    // a reverse read, with IUPAC codes and mixed qualities
    assert!(got.contains("@r2\nHVKMRYNACGTT\n+\nDCBA~JI?5#\"!\n"));
    // QUAL and SEQ missing
    assert!(got.contains("@r8\nACGTAC\n+\n!!!!!!\n"));
    assert!(got.contains("@r9\n\n+\n\n"));
}
