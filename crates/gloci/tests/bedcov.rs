//! `gloci bedcov` (D140): errors, filters, columns, threads.
mod common;
use common::*;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;
use std::path::PathBuf;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/samtools")
}

fn f(name: &str) -> String {
    dir().join(name).display().to_string()
}

/// `gloci bedcov ARGS` in `cwd`: (stdout, stderr, exit code).
fn bedcov_in(cwd: &std::path::Path, args: &[&str]) -> (String, String, i32) {
    let o = gloci()
        .current_dir(cwd)
        .arg("bedcov")
        .args(args)
        .output()
        .unwrap();
    (
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
        o.status.code().unwrap(),
    )
}

/// The fixtures copied beside fresh indexes, so the stale check never speaks.
fn copied(names: &[&str]) -> tempfile::TempDir {
    let t = tempdir();
    for n in names {
        copy_with_fresh_index(t.path(), &dir().join(n), &dir().join(format!("{n}.bai")));
    }
    t
}

#[test]
fn bedcov_needs_index() {
    let t = tempdir();
    std::fs::copy(dir().join("depth2.bam"), t.path().join("depth2.bam")).unwrap();
    let bed = plain_file(t.path(), "r.bed", "chr1\t0\t10\n");
    let (out, err, code) = bedcov_in(t.path(), &[bed.to_str().unwrap(), "depth2.bam"]);
    assert_eq!(code, 2, "{err}");
    assert!(out.is_empty());
    assert!(err.contains("gloci index depth2.bam"), "{err}");
}

#[test]
fn bedcov_bad_chrom() {
    let t = copied(&["depth.bam"]);
    let bad = f("bedcov-badchrom.bed");
    let bam = t.path().join("depth.bam").display().to_string();
    let (out, err, code) = bedcov_in(t.path(), &[&bad, &bam]);
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("line 9") && err.contains("chrZ"), "{err}");
    // The lines before it were written, and none after.
    assert_eq!(out.lines().count(), 5, "{out}");
    // --skip-invalid: one counted warning, every other line written.
    let (out, err, code) = bedcov_in(t.path(), &[&bad, &bam, "--skip-invalid"]);
    assert_eq!(code, 0, "{err}");
    assert_eq!(err.matches("skipped 1 BED line").count(), 1, "{err}");
    let (good, _, _) = bedcov_in(t.path(), &[&f("bedcov.bed"), &bam]);
    assert_eq!(out, good);
}

#[test]
fn bedcov_unmapped_never_counted() {
    let t = copied(&["depth.bam"]);
    let bam = t.path().join("depth.bam").display().to_string();
    let bed = plain_file(
        t.path(),
        "r.bed",
        "chr1\t170\t180\nchr1\t0\t500\nchr1\t0\t0\n",
    );
    let bed = bed.to_str().unwrap();
    // The placed unmapped read at chr1:175 does not count, with the default exclusion or none.
    let (out, err, code) = bedcov_in(
        t.path(),
        &[bed, &bam, "--read-count", "--exclude-flags", "0"],
    );
    assert_eq!(code, 0, "{err}");
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines[0], "chr1\t170\t180\t0\t0");
    // 15 reads pass the default exclusion; the duplicate, secondary and QC-fail ones join them.
    assert_eq!(lines[1].rsplit('\t').next(), Some("18"));
    assert_eq!(lines[2], "chr1\t0\t0\t0\t0");
    let (out, _, _) = bedcov_in(t.path(), &[bed, &bam, "--read-count"]);
    assert_eq!(out.lines().nth(1).unwrap().rsplit('\t').next(), Some("15"));
}

#[test]
fn bedcov_columns_and_options() {
    let t = copied(&["depth.bam", "depth2.bam"]);
    let a = t.path().join("depth.bam").display().to_string();
    let b = t.path().join("depth2.bam").display().to_string();
    let bed = plain_file(t.path(), "r.bed", "# c\nchr1\t0\t100\tx\n");
    let bed = bed.to_str().unwrap();
    let (out, err, code) = bedcov_in(
        t.path(),
        &[
            bed,
            &a,
            &b,
            "--column-names",
            "--bases-at-depth",
            "1",
            "--read-count",
        ],
    );
    assert_eq!(code, 0, "{err}");
    let mut lines = out.lines();
    assert_eq!(
        lines.next().unwrap(),
        format!("#chrom\tstart\tend\t{a}_cov\t{b}_cov\t{a}_depth\t{b}_depth\t{a}_count\t{b}_count")
    );
    let row = lines.next().unwrap();
    assert!(row.starts_with("chr1\t0\t100\tx\t"), "{row}");
    assert_eq!(row.split('\t').count(), 4 + 6);
    // stdin BED, and -o.
    let o = t.path().join("out.tsv");
    gloci()
        .current_dir(t.path())
        .args(["bedcov", "-", &a, "-o"])
        .arg(&o)
        .write_stdin("chr1\t0\t100\tx\n")
        .assert()
        .success();
    let written = std::fs::read_to_string(&o).unwrap();
    assert_eq!(written.lines().count(), 1);
    assert!(written.starts_with("chr1\t0\t100\tx\t"));
    // Differing @SQ lists are refused.
    let other = t.path().join("reordered.bam");
    std::fs::copy(dir().join("depth-reordered.bam"), &other).unwrap();
    gloci().arg("index").arg(&other).assert().success();
    gloci()
        .args(["bedcov", bed, other.to_str().unwrap(), &a])
        .assert()
        .failure()
        .stderr(contains("different @SQ lists"));
}

#[test]
fn bedcov_filters_need_valid_flags() {
    let t = copied(&["depth.bam"]);
    let a = t.path().join("depth.bam").display().to_string();
    let bed = plain_file(t.path(), "r.bed", "chr1\t0\t100\n");
    gloci()
        .args([
            "bedcov",
            bed.to_str().unwrap(),
            &a,
            "--exclude-flags",
            "bogus",
        ])
        .assert()
        .code(2)
        .stderr(contains("--exclude-flags").and(contains("bogus")));
}

/// One BAM record (without `block_size`): `len` M bases at `pos` (0-based) on `rid`, with
/// `l_seq` bases of sequence; `l_seq != len` makes the record invalid.
fn bam_record(rid: i32, pos: i32, name: &str, len: u32, l_seq: u32) -> Vec<u8> {
    let mut r = Vec::new();
    r.extend_from_slice(&rid.to_le_bytes());
    r.extend_from_slice(&pos.to_le_bytes());
    r.push(name.len() as u8 + 1);
    r.push(60);
    r.extend_from_slice(&0u16.to_le_bytes());
    r.extend_from_slice(&1u16.to_le_bytes());
    r.extend_from_slice(&0u16.to_le_bytes());
    r.extend_from_slice(&l_seq.to_le_bytes());
    r.extend_from_slice(&(-1i32).to_le_bytes());
    r.extend_from_slice(&(-1i32).to_le_bytes());
    r.extend_from_slice(&0i32.to_le_bytes());
    r.extend_from_slice(name.as_bytes());
    r.push(0);
    r.extend_from_slice(&(len << 4).to_le_bytes());
    r.extend(std::iter::repeat_n(0x11u8, l_seq.div_ceil(2) as usize));
    r.extend(std::iter::repeat_n(30u8, l_seq as usize));
    r
}

fn write_bam(dir: &std::path::Path, name: &str, bad: bool) -> PathBuf {
    let mut raw = b"BAM\x01".to_vec();
    raw.extend_from_slice(&0i32.to_le_bytes());
    raw.extend_from_slice(&2i32.to_le_bytes());
    for name in ["chr1", "chr2"] {
        raw.extend_from_slice(&5i32.to_le_bytes());
        raw.extend_from_slice(name.as_bytes());
        raw.push(0);
        raw.extend_from_slice(&100i32.to_le_bytes());
    }
    let len = if bad { 6 } else { 5 };
    for rec in [
        bam_record(0, 4, "a", 5, 5),
        bam_record(0, 9, "bad1", len, 5),
        bam_record(0, 19, "c", 5, 5),
        bam_record(1, 4, "d", 5, 5),
        bam_record(1, 9, "bad2", len, 5),
    ] {
        raw.extend_from_slice(&(rec.len() as i32).to_le_bytes());
        raw.extend_from_slice(&rec);
    }
    gloci_testkit::fixtures::write(dir, name, &gloci_testkit::fixtures::bgzip(&raw))
}

/// A BAM with two invalid records, indexed through a good twin of the same layout.
fn bad_bam(dir: &std::path::Path) -> PathBuf {
    let good = write_bam(dir, "good.bam", false);
    gloci().args(["index"]).arg(&good).assert().success();
    let bad = write_bam(dir, "bad.bam", true);
    let bai = |p: &PathBuf| PathBuf::from(format!("{}.bai", p.display()));
    std::fs::copy(bai(&good), bai(&bad)).unwrap();
    bad
}

#[test]
fn bedcov_errors_and_diagnostics_do_not_depend_on_threads() {
    let t = tempdir();
    let bam = bad_bam(t.path());
    let bam = bam.to_str().unwrap();
    // Enough lines for several batches, so more than one worker has work at -t 4.
    let mut text = String::new();
    for i in 0..600 {
        text.push_str(&format!("chr{}\t0\t{}\n", 1 + i % 2, 100 - i % 7));
    }
    let bed = plain_file(t.path(), "r.bed", &text);
    let bed = bed.to_str().unwrap();
    let run = |args: &[&str]| bedcov_in(t.path(), args);
    // Skipping: the count and the examples are the same at every -t.
    let (o1, e1, c1) = run(&[bed, bam, "--skip-invalid", "-t", "1"]);
    assert_eq!(c1, 0, "{e1}");
    assert!(e1.contains("invalid record"), "{e1}");
    for tt in ["2", "4"] {
        let (o, e, c) = run(&[bed, bam, "--skip-invalid", "-t", tt]);
        assert_eq!((c, &o, &e), (c1, &o1, &e1), "-t {tt}");
    }
    // Failing: the same record is named, by its place in the indexed form, and the same lines
    // were written.
    let (o1, e1, c1) = run(&[bed, bam, "-t", "1"]);
    assert_eq!(c1, 1, "{e1}");
    assert!(e1.contains("BGZF offset"), "{e1}");
    for tt in ["2", "4"] {
        let (o, e, c) = run(&[bed, bam, "-t", tt]);
        assert_eq!((c, &o, &e), (c1, &o1, &e1), "-t {tt}");
    }
}

#[test]
fn bedcov_unreadable_bed_line_leaves_the_same_stdout_at_every_thread_count() {
    let t = copied(&["depth.bam"]);
    let bam = t.path().join("depth.bam").display().to_string();
    let mut text = String::new();
    for i in 0..600 {
        text.push_str(&format!("chr1\t{}\t{}\n", i % 50, 60 + i % 40));
    }
    text.push_str("chr1\tx\t5\n");
    let bed = plain_file(t.path(), "r.bed", &text);
    let bed = bed.to_str().unwrap();
    let (o1, e1, c1) = bedcov_in(t.path(), &[bed, &bam, "-t", "1"]);
    assert_eq!(c1, 1, "{e1}");
    assert!(e1.contains("line 601"), "{e1}");
    assert!(o1.lines().count() >= 256, "{}", o1.lines().count());
    for tt in ["2", "4", "8"] {
        let (o, e, c) = bedcov_in(t.path(), &[bed, &bam, "-t", tt]);
        assert_eq!((c, &o, &e), (c1, &o1, &e1), "-t {tt}");
    }
}

/// Two references of 100 bases, each with a read of 10 bases at 96 (past the end).
fn clipped_bam(dir: &std::path::Path) -> PathBuf {
    let mut raw = b"BAM\x01".to_vec();
    raw.extend_from_slice(&0i32.to_le_bytes());
    raw.extend_from_slice(&2i32.to_le_bytes());
    for name in ["chr1", "chr2"] {
        raw.extend_from_slice(&5i32.to_le_bytes());
        raw.extend_from_slice(name.as_bytes());
        raw.push(0);
        raw.extend_from_slice(&100i32.to_le_bytes());
    }
    for rec in [
        bam_record(0, 95, "a", 10, 10),
        bam_record(1, 95, "b", 10, 10),
    ] {
        raw.extend_from_slice(&(rec.len() as i32).to_le_bytes());
        raw.extend_from_slice(&rec);
    }
    let p =
        gloci_testkit::fixtures::write(dir, "clipped.bam", &gloci_testkit::fixtures::bgzip(&raw));
    gloci().arg("index").arg(&p).assert().success();
    p
}

#[test]
fn bedcov_past_the_end_warning_is_the_first_in_bed_order_at_every_thread_count() {
    let t = tempdir();
    let bam = clipped_bam(t.path());
    let bam = bam.to_str().unwrap();
    // The first lines are on chr2, so batch order makes chr2 the reference the warning names.
    let mut text = String::new();
    for i in 0..1200 {
        let chrom = if i < 300 || i % 2 == 0 {
            "chr2"
        } else {
            "chr1"
        };
        text.push_str(&format!("{chrom}\t90\t100\n"));
    }
    let bed = plain_file(t.path(), "r.bed", &text);
    let bed = bed.to_str().unwrap();
    let (o1, e1, c1) = bedcov_in(t.path(), &[bed, bam, "-t", "1"]);
    assert_eq!(c1, 0, "{e1}");
    assert_eq!(e1.matches("past the end").count(), 1, "{e1}");
    assert!(e1.contains("chr2 has LN 100"), "{e1}");
    for tt in ["2", "8"] {
        for _ in 0..5 {
            let (o, e, c) = bedcov_in(t.path(), &[bed, bam, "-t", tt]);
            assert_eq!((c, &o, &e), (c1, &o1, &e1), "-t {tt}");
        }
    }
}
