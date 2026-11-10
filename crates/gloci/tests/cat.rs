//! Tests for `gloci cat` (D123): BAM block copy, as samtools cat.
mod common;
use common::*;
use predicates::prelude::*;
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

/// `view -O sam` of a BAM, without the @PG lines.
fn sam(bam: &Path) -> String {
    let o = gloci()
        .args(["view", &s(bam), "-O", "sam"])
        .assert()
        .success();
    String::from_utf8(o.get_output().stdout.clone())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(|l| format!("{l}\n"))
        .collect()
}

/// The records of `sam` text, without its header.
fn records(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('@')).collect()
}

#[test]
fn two_files_join_with_the_extra_read_group() {
    let dir = tempdir();
    let out = dir.path().join("o.bam");
    gloci()
        .args(["cat", &s(&fixture("cat-a.bam")), &s(&fixture("cat-b.bam"))])
        .args(["-o", &s(&out)])
        .assert()
        .success();
    let text = sam(&out);
    let a = sam(&fixture("cat-a.bam"));
    let b = sam(&fixture("cat-b.bam"));
    // header: a's lines, then b's new @RG at the end
    let header: Vec<&str> = text.lines().take_while(|l| l.starts_with('@')).collect();
    assert_eq!(
        header,
        [
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:100000",
            "@SQ\tSN:chr2\tLN:100000",
            "@RG\tID:g1\tSM:s1",
            "@CO\tfirst",
            "@RG\tID:g2\tSM:s2"
        ]
    );
    let want: Vec<&str> = records(&a).into_iter().chain(records(&b)).collect();
    assert_eq!(records(&text), want);
    // gloci's own @PG line is there
    let o = gloci()
        .args(["view", &s(&out), "-O", "sam", "--header-only"])
        .assert()
        .success();
    let full = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(full.contains("@PG\tID:gloci\tPN:gloci"), "{full}");
}

#[test]
fn different_sq_lists_exit_1_and_write_nothing() {
    let dir = tempdir();
    let out = dir.path().join("o.bam");
    gloci()
        .args(["cat", &s(&fixture("cat-a.bam")), &s(&fixture("cat-c.bam"))])
        .args(["-o", &s(&out)])
        .assert()
        .code(1)
        .stderr(contains("cat-c.bam"))
        .stderr(contains("reference 1 is chr1"))
        .stderr(contains("but chr2"));
    assert!(!out.exists());
}

#[test]
fn cram_sam_and_other_input_are_refused() {
    let dir = tempdir();
    let cram = plain_file(
        dir.path(),
        "x.cram",
        &format!("CRAM\x03\x00{}", "x".repeat(100)),
    );
    gloci()
        .args(["cat", &s(&cram)])
        .assert()
        .code(2)
        .stderr(contains("CRAM"))
        .stderr(contains("samtools cat"));
    gloci()
        .args(["cat", &s(&fixture("cat-a.sam"))])
        .assert()
        .code(2)
        .stderr(contains("SAM"))
        .stderr(contains("gloci view"))
        .stderr(contains("-o x.bam"));
    let bed = plain_file(dir.path(), "x.bed", "chr1\t1\t2\n");
    gloci()
        .args(["cat", &s(&bed)])
        .assert()
        .code(1)
        .stderr(contains("not BAM"));
    let vcf = bgzf_file(
        dir.path(),
        "x.vcf.gz",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\n",
    );
    gloci()
        .args(["cat", &s(&vcf)])
        .assert()
        .code(1)
        .stderr(contains("not BAM"));
}

#[test]
fn damaged_input_is_an_error_not_a_panic() {
    let dir = tempdir();
    let good = std::fs::read(fixture("cat-a.bam")).unwrap();
    let out = dir.path().join("o.bam");
    for (name, bytes) in [
        ("empty.bam", Vec::new()),
        ("cut-header.bam", good[..30].to_vec()),
        ("cut-middle.bam", good[..good.len() - 40].to_vec()),
        ("no-eof.bam", good[..good.len() - 28].to_vec()),
        (
            "garbage.bam",
            vec![
                0x1f, 0x8b, 8, 4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
            ],
        ),
    ] {
        let p = dir.path().join(name);
        std::fs::write(&p, bytes).unwrap();
        let o = gloci()
            .args(["cat", &s(&fixture("cat-a.bam")), &s(&p), "-o", &s(&out)])
            .assert()
            .code(1);
        let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
        assert!(err.contains(name), "{name}: {err}");
        assert!(!err.contains("panicked"), "{name}: {err}");
        assert!(!out.exists(), "{name} left an output behind");
    }
}

/// `--write-index` on joined files that are not in coordinate order is an error naming what
/// to sort, and nothing is published (`write_index_on_bam_output` has the working case).
#[test]
fn write_index_of_unsorted_cat_output_is_an_error() {
    let dir = tempdir();
    let out = dir.path().join("o.bam");
    gloci()
        .args([
            "cat",
            &s(&fixture("cat-a.bam")),
            &s(&fixture("cat-b.bam")),
            "-o",
            &s(&out),
        ])
        .arg("--write-index")
        .assert()
        .code(1)
        .stderr(contains("gloci sort").and(contains("panicked").not()));
    assert!(!out.exists());
    assert!(!dir.path().join("o.bam.bai").exists());
}

#[test]
fn output_must_be_bam() {
    let dir = tempdir();
    gloci()
        .args([
            "cat",
            &s(&fixture("cat-a.bam")),
            "-o",
            &s(&dir.path().join("o.sam")),
        ])
        .assert()
        .code(2)
        .stderr(contains("cat writes BAM"));
}

#[test]
fn stdin_is_one_input_and_stdout_is_the_default() {
    let b = std::fs::read(fixture("cat-b.bam")).unwrap();
    let o = gloci()
        .args(["cat", &s(&fixture("cat-a.bam")), "-"])
        .write_stdin(b)
        .assert()
        .success();
    let dir = tempdir();
    let joined = dir.path().join("j.bam");
    std::fs::write(&joined, &o.get_output().stdout).unwrap();
    let text = sam(&joined);
    assert_eq!(records(&text).len(), 6);
    assert!(text.contains("@RG\tID:g2"));
}

/// A sorted BAM of `n` reads of noise, made by `gloci view`.
fn big_bam(dir: &Path, name: &str, n: usize, seed: u64) -> PathBuf {
    let mut x = seed;
    let mut next = || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (x >> 33) as usize
    };
    let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n");
    for i in 0..n {
        let seq: String = (0..100).map(|_| b"ACGT"[next() % 4] as char).collect();
        let qual: String = (0..100).map(|_| (33 + next() % 40) as u8 as char).collect();
        t.push_str(&format!(
            "{name}{i}\t0\tchr1\t{}\t60\t100M\t*\t0\t0\t{seq}\t{qual}\n",
            i + 1
        ));
    }
    let sam = plain_file(dir, &format!("{name}.sam"), &t);
    let bam = dir.join(format!("{name}.bam"));
    gloci()
        .args(["view", &s(&sam), "-o", &s(&bam)])
        .assert()
        .success();
    bam
}

#[test]
fn later_blocks_are_copied_byte_for_byte() {
    let dir = tempdir();
    let a = big_bam(dir.path(), "a", 3000, 1);
    let b = big_bam(dir.path(), "b", 3000, 2);
    let out = dir.path().join("o.bam");
    gloci()
        .args(["cat", &s(&a), &s(&b), "-o", &s(&out)])
        .assert()
        .success();
    let (bytes_a, bytes_b, joined) = (
        std::fs::read(&a).unwrap(),
        std::fs::read(&b).unwrap(),
        std::fs::read(&out).unwrap(),
    );
    let eof = gloci_io::BGZF_EOF.len();
    assert!(bytes_b.len() > 100_000, "the fixture spans many blocks");
    assert!(joined.ends_with(&gloci_io::BGZF_EOF));
    // b's blocks after the one holding its header end are in the output as they were, ...
    let tail = &bytes_b[bytes_b.len() - eof - 60_000..bytes_b.len() - eof];
    assert!(joined[..joined.len() - eof].ends_with(tail));
    // ... and so are a's, right after the new header.
    let head = &bytes_a[bytes_a.len() - eof - 60_000..bytes_a.len() - eof];
    assert!(joined.windows(head.len()).any(|w| w == head));
    // one EOF marker only
    assert_eq!(
        joined
            .windows(eof)
            .filter(|w| *w == gloci_io::BGZF_EOF)
            .count(),
        1
    );
    let (ta, tb, tj) = (sam(&a), sam(&b), sam(&out));
    let want: Vec<&str> = records(&ta).into_iter().chain(records(&tb)).collect();
    assert_eq!(records(&tj), want);
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let dir = tempdir();
    let a = big_bam(dir.path(), "a", 1500, 3);
    let b = big_bam(dir.path(), "b", 1500, 4);
    // Cat does no threaded compression, so this cannot fail on thread count by itself. It guards
    // that `-t` is stripped from the @PG CL line and that only one block is recompressed per
    // input. It writes to stdout so the CL line is the same too.
    let run = |t: &str| {
        let o = gloci()
            .args(["cat", &s(&a), &s(&b), "-t", t])
            .assert()
            .success();
        o.get_output().stdout.clone()
    };
    let one = run("1");
    assert_eq!(run("2"), one);
    assert_eq!(run("8"), one);
}

#[test]
fn vcf_concat_naive_still_joins() {
    // The copy loop is shared with `vcf concat --naive`; its own tests are in vcf_concat.rs.
    let dir = tempdir();
    let h = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let a = bgzf_file(
        dir.path(),
        "a.vcf.gz",
        &format!("{h}chr1\t1\t.\tA\tC\t.\t.\t.\n"),
    );
    let b = bgzf_file(
        dir.path(),
        "b.vcf.gz",
        &format!("{h}chr1\t2\t.\tA\tC\t.\t.\t.\n"),
    );
    let out = dir.path().join("o.vcf.gz");
    gloci()
        .args(["vcf", "concat", "--naive", &s(&a), &s(&b), "-o", &s(&out)])
        .assert()
        .success();
    let o = gloci()
        .args(["view", &s(&out), "-O", "vcf"])
        .assert()
        .success();
    let text = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(text.contains("chr1\t1\t") && text.contains("chr1\t2\t"));
}
