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
    let targets = p(&i.join("samtools/view-targets.bed"));
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
        vec!["--no-header"],
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

/// `view --count` of a CRAM without regions reads the flags and MAPQ alone: for every CRAM
/// fixture, at one thread and several, it prints the number of records the full decode
/// (`view`) writes, with and without read filters.
#[test]
fn view_count_cram_equals_the_records_view_writes() {
    let fa = p(&inputs().join("fa/ref.fa"));
    let mut reads = 0;
    for (cram, _, records) in decodable_crams() {
        let cram = p(&cram);
        reads += records.len();
        let flag_and_mapq = |r: &String| {
            let mut f = r.split('\t');
            let flag: u16 = f.nth(1).unwrap().parse().unwrap();
            let mapq: u8 = f.nth(2).unwrap().parse().unwrap();
            (flag, mapq)
        };
        let mapped = records
            .iter()
            .filter(|r| flag_and_mapq(r).0 & 0x4 == 0)
            .count() as u64;
        let mapq30 = records.iter().filter(|r| flag_and_mapq(r).1 >= 30).count() as u64;
        for t in ["1", "4"] {
            let base = ["--reference", fa.as_str(), "-t", t, cram.as_str()];
            assert_eq!(count(&base), records.len() as u64, "{cram} -t {t}");
            let mut a = vec!["--exclude-flags", "UNMAP"];
            a.extend(base);
            assert_eq!(count(&a), mapped, "{cram} -t {t} --exclude-flags UNMAP");
            let mut a = vec!["--min-mapq", "30"];
            a.extend(base);
            assert_eq!(count(&a), mapq30, "{cram} -t {t} --min-mapq 30");
        }
    }
    assert!(reads > 1000, "{reads} reads compared");
}

/// From a CRAM, `view --count` without regions decodes no bases (#cram-limits), as flagstat:
/// a FASTA with the right names and lengths but other bases gives the same count, where
/// `view`, which decodes them, refuses it; overhang.cram's two reads are counted. With a
/// region the positions are read, so the records of that region are decoded in full and the
/// wrong FASTA is refused. A CRAM that needs a reference still needs --reference.
#[test]
fn view_count_cram_reads_no_bases() {
    use predicates::str::contains;
    let d = tempdir();
    let wrong = p(&wrong_bases_fasta(d.path()));
    let fa = p(&inputs().join("fa/ref.fa"));
    let mates = p(&inputs().join("aln/mates.cram"));
    let n = count(&["--reference", &fa, &mates]);
    assert!(n > 0);
    assert_eq!(count(&["--reference", &wrong, &mates]), n);
    gloci()
        .args(["view", "--reference", &wrong, &mates])
        .assert()
        .code(1)
        .stderr(contains("does not match"));
    gloci()
        .args(["view", "--count", "--reference", &wrong, &mates, "chr1"])
        .assert()
        .code(1)
        .stderr(contains("does not match"));
    let overhang = p(&inputs().join("aln/overhang.cram"));
    assert_eq!(count(&["--reference", &fa, &overhang]), 2);
    gloci()
        .args(["view", "--count", &mates])
        .assert()
        .code(1)
        .stderr(contains("--reference"));
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

/// A FIFO (process substitution) cannot be indexed: the scan warning suggests no `gloci index`
/// for it (O1), as for stdin, and still does for a regular file.
#[test]
fn view_count_region_of_a_fifo_suggests_no_index() {
    let bam = p(&inputs().join("aln/flags.bam"));
    let o = bash_gloci(&format!("view --count <(cat '{bam}') chr1"), "");
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(String::from_utf8_lossy(&o.stdout), "11\n");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("scanning the whole file /dev/fd/"), "{err}");
    assert!(!err.contains("gloci index"), "{err}");
    let dir = tempdir();
    let plain = dir.path().join("plain.bam");
    std::fs::copy(&bam, &plain).unwrap();
    let o = gloci()
        .args(["view", "--count", &p(&plain), "chr1"])
        .assert()
        .success()
        .get_output()
        .clone();
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("create an index with: gloci index"), "{err}");
}

#[test]
fn view_count_is_identical_at_every_thread_count() {
    let i = inputs();
    let bam = p(&i.join("aln/flags.bam"));
    let sam = p(&i.join("aln/flags.sam"));
    let targets = p(&i.join("samtools/view-targets.bed"));
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
