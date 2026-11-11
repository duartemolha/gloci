//! `gloci depth` (D137, D138): rows, filters, regions, several inputs, threads, errors.
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

fn out(args: &[&str]) -> String {
    let o = gloci()
        .arg("depth")
        .args(args)
        .assert()
        .success()
        .get_output()
        .clone();
    String::from_utf8(o.stdout).unwrap()
}

#[test]
fn depth_rows_of_a_small_file() {
    let t = tempdir();
    let sam = plain_file(
        t.path(),
        "a.sam",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:20\n@SQ\tSN:c2\tLN:5\n\
r1\t0\tc1\t3\t60\t3M2D2M\t*\t0\t0\tACGTA\tIIIII\n\
r2\t0\tc1\t4\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let p = sam.display().to_string();
    assert_eq!(
        out(&[&p]),
        "c1\t3\t1\nc1\t4\t2\nc1\t5\t2\nc1\t6\t1\nc1\t7\t1\nc1\t8\t1\nc1\t9\t1\n"
    );
    assert_eq!(
        out(&[&p, "--count-deletions"]),
        "c1\t3\t1\nc1\t4\t2\nc1\t5\t2\nc1\t6\t2\nc1\t7\t2\nc1\t8\t1\nc1\t9\t1\n"
    );
    // -a: all of c1 (it has reads), none of c2; -aa: c2 too.
    let all = out(&[&p, "--all-positions"]);
    assert_eq!(all.lines().count(), 20);
    assert!(all.starts_with("c1\t1\t0\nc1\t2\t0\nc1\t3\t1\n"));
    assert!(all.ends_with("c1\t20\t0\n"));
    let aa = out(&[&p, "--all-references"]);
    assert_eq!(aa.lines().count(), 25);
    assert!(aa.ends_with("c2\t5\t0\n"));
    // --column-names
    let named = out(&[&p, "--column-names", &p]);
    assert!(named.starts_with(&format!("#CHROM\tPOS\t{p}\t{p}\nc1\t3\t1\t1\n")));
}

#[test]
fn depth_exclude_flags_default_is_samtools_and_zero_turns_it_off() {
    let f1 = f("depth.bam");
    let base = out(&[&f1]);
    // The placed unmapped read at chr1:175, a duplicate (200), a secondary (210) and a QC-fail
    // read (220) are not counted by default; the supplementary read at 230 is.
    for p in [175, 200, 210, 220] {
        assert!(!base.contains(&format!("chr1\t{p}\t")), "{p}");
    }
    assert!(base.contains("chr1\t230\t1\n"));
    let none = out(&[&f1, "--exclude-flags", "0"]);
    assert!(
        none.contains("chr1\t175\t0\n"),
        "a placed unmapped read prints 0"
    );
    for p in [200, 210, 220] {
        assert!(none.contains(&format!("chr1\t{p}\t1\n")), "{p}");
    }
    let same = out(&[&f1, "--exclude-flags", "unmap,secondary,qcfail,dup"]);
    assert_eq!(same, base);
    let sup = out(&[
        &f1,
        "--exclude-flags",
        "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
    ]);
    assert!(!sup.contains("chr1\t230\t"));
    gloci().args(["depth", "--help"]).assert().success().stdout(
        contains("default UNMAP,SECONDARY,QCFAIL,DUP").or(contains("UNMAP,SECONDARY,QCFAIL,DUP")),
    );
}

#[test]
fn depth_filters() {
    let f1 = f("depth.bam");
    // MAPQ 29 and 30 reads at 15 and 16: -Q 30 drops only the 29.
    let q = out(&[&f1, "--min-mapq", "30"]);
    assert!(q.contains("chr1\t16\t2\n") && q.contains("chr1\t24\t2\n"));
    assert!(out(&[&f1]).contains("chr1\t16\t3\n"));
    // -l 50: the 49-base read goes, the 50-base one stays, the clipped 30-base one goes.
    let l = out(&[&f1, "--min-read-length", "50"]);
    assert!(!l.contains("chr1\t300\t"));
    assert!(l.contains("chr1\t335\t1\n"));
    assert!(out(&[&f1, "--min-read-length", "49"]).contains("chr1\t340\t2\n"));
    assert!(!l.contains("chr1\t385\t"));
    // --require-flags: only the proper pair at 250 and 255.
    let r = out(&[&f1, "--require-flags", "PROPER_PAIR"]);
    assert_eq!(r.lines().next().unwrap(), "chr1\t250\t1");
    assert_eq!(r.lines().count(), 25);
}

#[test]
fn depth_regions_through_the_index_and_by_scanning() {
    let bam = f("depth.bam");
    let want = out(&[&bam, "--region", "chr1:100-200"]);
    assert!(want.lines().count() > 0);
    // The SAM is scanned, with a warning, and gives the same rows.
    gloci()
        .args(["depth", &f("depth.sam"), "--region", "chr1:100-200"])
        .assert()
        .success()
        .stderr(contains("scanning the whole file"))
        .stdout(want.clone());
    // A regions file, and the union of regions with a point in two of them once.
    let bed = f("targets.bed");
    let a = out(&[&bam, "--regions-file", &bed]);
    let b = out(&[
        &bam,
        "--region",
        "chr1:6-30",
        "--region",
        "chr1:21-45",
        "--region",
        "chr1:96-130",
        "--region",
        "chr1:171-180",
        "--region",
        "chr1:461-480",
        "--region",
        "chr3:11-20",
    ]);
    assert_eq!(a, b);
    // With -a, every position of the regions: chr3 has no read.
    let c = out(&[&bam, "--regions-file", &bed, "--all-positions"]);
    assert!(c.contains("chr3\t11\t0\n") && c.contains("chr3\t20\t0\n"));
    assert_eq!(c.lines().filter(|l| l.starts_with("chr3\t")).count(), 10);
    // A region on a reference the file lacks is dropped with a warning.
    gloci()
        .args([
            "depth",
            &bam,
            "--region",
            "chrZ:1-10",
            "--region",
            "chr1:10-12",
        ])
        .assert()
        .success()
        .stderr(contains("chrZ"))
        .stdout("chr1\t10\t1\nchr1\t11\t1\nchr1\t12\t1\n");
    // Regions need a file to read the references from.
    gloci()
        .args(["depth", "-", "--region", "chr1"])
        .write_stdin("")
        .assert()
        .code(2)
        .stderr(contains("standard input"));
}

#[test]
fn depth_past_the_end_of_a_reference_warns_once() {
    let o = gloci()
        .args(["depth", &f("depth-pastend.bam")])
        .assert()
        .success()
        .get_output()
        .clone();
    let err = String::from_utf8(o.stderr).unwrap();
    assert_eq!(
        err.matches("past the end of its reference").count(),
        1,
        "{err}"
    );
    assert!(err.contains("chr2") && err.contains("LN 300"));
    let rows = String::from_utf8(o.stdout).unwrap();
    assert!(rows.contains("chr2\t300\t2\n") && !rows.contains("chr2\t301\t"));
}

#[test]
fn depth_inputs_must_agree_and_be_sorted() {
    gloci()
        .args(["depth", &f("depth-reordered.bam"), &f("depth.bam")])
        .assert()
        .code(1)
        .stderr(contains("@SQ 1"))
        .stdout("");
    gloci()
        .args(["depth", &f("depth-unsorted.bam")])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"));
    gloci()
        .args(["depth", &f("depth-unsorted.sam")])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"));
    gloci().args(["depth"]).assert().code(2);
    gloci()
        .args(["depth", &f("depth.bam"), "--exclude-flags", "bogus"])
        .assert()
        .code(2)
        .stderr(contains("--exclude-flags"));
}

#[test]
fn depth_threads_give_identical_output() {
    let bam = f("depth.bam");
    let bed = f("targets.bed");
    for extra in [
        vec![],
        vec!["--all-references"],
        vec!["--all-positions", "--regions-file", bed.as_str()],
        vec!["--count-deletions", "--min-baseq", "20"],
    ] {
        let run = |t: &str| {
            let mut args = vec![bam.as_str(), "-t", t];
            args.extend(extra.iter().copied());
            out(&args)
        };
        let one = run("1");
        assert_eq!(run("2"), one, "{extra:?}");
        assert_eq!(run("8"), one, "{extra:?}");
    }
    // Two inputs, one of them in another format.
    let two = |t: &str| out(&[&bam, &f("depth2.bam"), "-t", t, "--all-references"]);
    assert_eq!(two("1"), two("4"));
}

#[test]
fn depth_reads_sam_cram_and_stdin_alike() {
    let bam = out(&[&f("depth.bam")]);
    assert_eq!(out(&[&f("depth.sam")]), bam);
    // The CRAM was written without a reference.
    assert_eq!(out(&[&f("depth.cram")]), bam);
    let o = gloci()
        .args(["depth", "-"])
        .write_stdin(std::fs::read(f("depth.sam")).unwrap())
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(String::from_utf8(o.stdout).unwrap(), bam);
}

#[test]
fn depth_writes_a_file_with_o() {
    let t = tempdir();
    let p = t.path().join("d.tsv");
    gloci()
        .args(["depth", &f("depth.bam"), "-o"])
        .arg(&p)
        .assert()
        .success()
        .stdout("");
    assert_eq!(
        std::fs::read_to_string(&p).unwrap(),
        out(&[&f("depth.bam")])
    );
}

#[test]
fn depth_noseq_reads_use_their_cigar_for_the_length() {
    // SEQ * : -l counts the bases the CIGAR puts in the read (n1 50M kept, n2 49M dropped).
    let o = out(&[&f("depth-noseq.bam"), "--min-read-length", "50"]);
    assert!(o.contains("chr1\t400\t1\n"));
    assert!(!o.contains("chr1\t445\t2\n"));
    assert!(out(&[&f("depth-noseq.bam"), "--min-read-length", "49"]).contains("chr1\t445\t2\n"));
    // Every base of such a read passes -q.
    assert!(out(&[&f("depth-noseq.bam"), "--min-baseq", "40"]).contains("chr1\t400\t1\n"));
}
