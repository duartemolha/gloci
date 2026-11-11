//! `gloci view --count` and `--regions-file` on alignments (D133): the count equals the records
//! `view` writes, for every format; the conflicts; the filters; the index is used for a regions
//! file; the output does not depend on `-t`. The memory bound is in
//! `gloci-ops/tests/view_count_memory.rs`.
mod common;
use common::*;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn p(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success().get_output().clone();
    String::from_utf8(o.stdout).unwrap()
}

fn count(args: &[&str]) -> u64 {
    let mut a = vec!["view", "--count"];
    a.extend_from_slice(args);
    let text = out(&a);
    assert!(
        text.ends_with('\n') && text.matches('\n').count() == 1,
        "{text:?}"
    );
    text.trim().parse().unwrap()
}

/// The number of records `view` writes: for alignments, the SAM lines without a header; for
/// the rest, the lines that are not header lines.
fn written(args: &[&str], sam: bool) -> u64 {
    let mut a = vec!["view"];
    if sam {
        a.extend(["--no-header", "-O", "sam"]);
    } else if args[0].ends_with(".bcf") {
        // BCF needs its header: write text and drop the header lines below.
        a.extend(["-O", "vcf"]);
    }
    a.extend_from_slice(args);
    out(&a).lines().filter(|l| !l.starts_with('#')).count() as u64
}

#[test]
fn view_count_formats() {
    let i = inputs();
    let flags_bam = p(&i.join("aln/flags.bam"));
    let flags_sam = p(&i.join("aln/flags.sam"));
    let cram = p(&i.join("aln/flags.cram"));
    let fa = p(&i.join("fa/ref.fa"));
    let vcf = p(&i.join("regions.vcf"));
    let bed = p(&i.join("regions.bed"));
    let dir = tempdir();
    let bcf = dir.path().join("x.bcf");
    gloci()
        .args(["view", &vcf, "-o", &p(&bcf)])
        .assert()
        .success();
    let bcf = p(&bcf);
    let targets = p(&i.join("samtools/targets.bed"));
    let regions = [
        vec![],
        vec!["chr1:1-150", "chr2:95-105"],
        vec!["chr1:101-101"],
    ];
    for r in &regions {
        for (name, file, extra, sam) in [
            ("bam", &flags_bam, vec![], true),
            ("sam", &flags_sam, vec![], true),
            ("cram", &cram, vec!["--reference", fa.as_str()], true),
            ("bed", &bed, vec![], false),
        ] {
            let mut a: Vec<&str> = extra.clone();
            a.push(file);
            a.extend(r.iter().copied());
            assert_eq!(count(&a), written(&a, sam), "{name} {r:?}");
        }
        // variants use chr1/chr2 of regions.vcf; BED uses its own names.
        for (name, file) in [("vcf", &vcf), ("bcf", &bcf)] {
            let mut a: Vec<&str> = vec![file];
            a.extend(r.iter().copied());
            assert_eq!(count(&a), written(&a, false), "{name} {r:?}");
        }
    }
    // A regions file, alignments and a BED.
    for (file, sam) in [(&flags_bam, true), (&flags_sam, true)] {
        let a = [file.as_str(), "--regions-file", targets.as_str()];
        assert_eq!(count(&a), written(&a, sam));
    }
    assert_eq!(count(&[&flags_bam]), 14);
    assert_eq!(count(&[&bed, "--regions-file", &bed]), 4);
    // The output goes to -o too.
    let o = dir.path().join("n.txt");
    gloci()
        .args(["view", "--count", &flags_bam, "-o", &p(&o)])
        .assert()
        .success()
        .stdout("");
    assert_eq!(std::fs::read_to_string(&o).unwrap(), "14\n");
    // ... whatever the extension says.
    let o = dir.path().join("n.bam");
    gloci()
        .args(["view", "--count", &flags_bam, "-o", &p(&o)])
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&o).unwrap(), "14\n");
}

#[test]
fn view_count_conflicts() {
    let bam = p(&inputs().join("aln/flags.bam"));
    for extra in [
        vec!["--header-only"],
        vec!["-O", "bam"],
        vec!["-O", "sam"],
        vec!["--write-index"],
    ] {
        let mut a = vec!["view", "--count", bam.as_str()];
        a.extend(extra.iter().copied());
        gloci().args(&a).assert().code(2);
    }
}

#[test]
fn view_count_filters() {
    // samtools 1.21 on flags.bam: `-c -F 4` and `-c -q 30` (committed numbers).
    let bam = p(&inputs().join("aln/flags.bam"));
    assert_eq!(count(&["--exclude-flags", "UNMAP", &bam]), 12);
    assert_eq!(count(&["--min-mapq", "30", &bam]), 10);
    assert_eq!(count(&["--exclude-flags", "UNMAP,SECONDARY", &bam]), 11);
    // The same filters limit the records view writes.
    assert_eq!(written(&["--min-mapq", "30", &bam], true), 10);
}

/// A BAM of `n` reads on chr1, with a .bai, and the BED regions file for its first reads.
fn big_bam(dir: &Path, n: usize) -> PathBuf {
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n");
    let seq = "ACGT".repeat(40);
    let qual = "I".repeat(seq.len());
    for i in 0..n {
        sam.push_str(&format!(
            "read{i}\t0\tchr1\t{}\t60\t160M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0\n",
            i * 100 + 1
        ));
    }
    let sam_path = plain_file(dir, "big.sam", &sam);
    let bam = dir.join("big.bam");
    gloci()
        .args(["view", &p(&sam_path), "-o", &p(&bam), "--write-index"])
        .assert()
        .success();
    assert!(dir.join("big.bam.bai").exists());
    bam
}

#[test]
fn regions_file_uses_index() {
    let dir = tempdir();
    let bam = big_bam(dir.path(), 60_000);
    let bed = plain_file(
        dir.path(),
        "t.bed",
        "chr1\t0\t300\nchr1\t5000\t5100\nchr1\t100\t120\n",
    );
    let indexed = {
        let o = gloci()
            .args([
                "view",
                "-O",
                "sam",
                "--no-header",
                &p(&bam),
                "--regions-file",
                &p(&bed),
            ])
            .assert()
            .success()
            .get_output()
            .clone();
        assert!(
            !String::from_utf8_lossy(&o.stderr).contains("scanning"),
            "{:?}",
            o.stderr
        );
        String::from_utf8(o.stdout).unwrap()
    };
    assert_eq!(indexed.lines().count(), 3 + 2 + 1 - 1);
    // Destroy the middle and the end of the file: only a read that stays inside the
    // index's chunks for the regions still works.
    let mut bytes = std::fs::read(&bam).unwrap();
    let (from, to) = (bytes.len() / 2, bytes.len() - 100);
    bytes[from..to].fill(0xA5);
    let broken = dir.path().join("broken.bam");
    std::fs::write(&broken, &bytes).unwrap();
    std::fs::copy(
        dir.path().join("big.bam.bai"),
        dir.path().join("broken.bam.bai"),
    )
    .unwrap();
    // The regions touch reads 0-2 and 49-50, in the first part of the file.
    let again = out(&[
        "view",
        "-O",
        "sam",
        "--no-header",
        &p(&broken),
        "--regions-file",
        &p(&bed),
    ]);
    assert_eq!(again, indexed, "an indexed read touches only its chunks");
    // Unindexed, the same file is scanned, with one warning, and fails on the damage ...
    std::fs::remove_file(dir.path().join("broken.bam.bai")).unwrap();
    gloci()
        .args([
            "view",
            "-O",
            "sam",
            "--no-header",
            &p(&broken),
            "--regions-file",
            &p(&bed),
        ])
        .assert()
        .code(1);
    // ... and on the intact file gives the same records, with one warning.
    let plain = dir.path().join("plain.bam");
    std::fs::copy(&bam, &plain).unwrap();
    let o = gloci()
        .args([
            "view",
            "-O",
            "sam",
            "--no-header",
            &p(&plain),
            "--regions-file",
            &p(&bed),
        ])
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(String::from_utf8(o.stdout).unwrap(), indexed);
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    assert_eq!(err.matches("scanning the whole file").count(), 1, "{err}");
}

#[test]
fn view_count_is_identical_at_every_thread_count() {
    let i = inputs();
    let bam = p(&i.join("aln/flags.bam"));
    let sam = p(&i.join("aln/flags.sam"));
    let targets = p(&i.join("samtools/targets.bed"));
    let cram = p(&i.join("aln/flags.cram"));
    let fa = p(&i.join("fa/ref.fa"));
    let lines: Vec<Vec<&str>> = vec![
        vec![&bam],
        vec![&sam, "chr1:1-150"],
        vec![&bam, "--regions-file", &targets],
        vec![&cram, "--reference", &fa, "--regions-file", &targets],
    ];
    for l in &lines {
        let mut first = None;
        for t in ["1", "2", "8"] {
            let mut a = vec!["-t", t];
            a.extend(l.iter().copied());
            let n = count(&a);
            assert_eq!(*first.get_or_insert(n), n, "{l:?} at -t {t}");
        }
    }
}
