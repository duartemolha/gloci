//! `gloci maskfasta` (D113): streamed masking that keeps the FASTA's layout.
//! The bedtools comparison lives in the golden cases `maskfasta-*`.
#![allow(clippy::zombie_processes)] // the streaming test reaps its child with wait4
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Write;
use std::path::{Path, PathBuf};

fn fa(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/fa")
        .join(name)
}

fn run_ok(args: &[&str]) -> Vec<u8> {
    let out = gloci().arg("maskfasta").args(args).output().unwrap();
    assert!(
        out.status.success(),
        "gloci maskfasta {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn mask_with(dir: &Path, fasta: &[u8], bed: &str, extra: &[&str]) -> Vec<u8> {
    let f = gloci_testkit::fixtures::write(dir, "in.fa", fasta);
    let b = plain_file(dir, "m.bed", bed);
    let mut args = vec!["--fasta", f.to_str().unwrap(), "-b", b.to_str().unwrap()];
    args.extend_from_slice(extra);
    run_ok(&args)
}

#[test]
fn masks_overlaps_once_and_keeps_layout() {
    let d = tempdir();
    let input = b">a desc here\nACGTACGTAC\nACGT\n\n>b\r\nACGTAC\r\nAC\r\n>c\nACGT";
    // 0-based half-open; the two masks on a overlap, one on b crosses a CRLF line.
    let out = mask_with(
        d.path(),
        input,
        "a\t2\t6\na\t4\t8\nb\t4\t7\nzz\t0\t1\n",
        &[],
    );
    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">a desc here\nACNNNNNNAC\nACGT\n\n>b\r\nACGTNN\r\nNC\r\n>c\nACGT"
    );
    // The count is of distinct bases: 6 on a, 3 on b.
    let f = fa("ref.fa");
    let b = fa("masks.bed");
    gloci()
        .args(["-v", "maskfasta", "--fasta"])
        .arg(&f)
        .arg("-b")
        .arg(&b)
        .assert()
        .success()
        .stderr(contains("masked 188 bases").and(contains("chrZ")));
}

#[test]
fn soft_mask_keeps_non_letters() {
    let d = tempdir();
    let out = mask_with(d.path(), b">s\nACGT-*ACGT..NN\n", "s\t0\t12\n", &["--soft"]);
    assert_eq!(String::from_utf8(out).unwrap(), ">s\nacgt-*acgt..NN\n");
    // A hard mask replaces every base, letters or not.
    let out = mask_with(
        d.path(),
        b">s\nACGT-*ACGT..NN\n",
        "s\t0\t12\n",
        &["--mask-char", "x"],
    );
    assert_eq!(String::from_utf8(out).unwrap(), ">s\nxxxxxxxxxxxxNN\n");
}

#[test]
fn mask_char_validation() {
    let d = tempdir();
    let f = gloci_testkit::fixtures::write(d.path(), "in.fa", b">s\nACGT\n");
    let b = plain_file(d.path(), "m.bed", "s\t0\t2\n");
    for bad in ["NN", "", "\u{7}", "é"] {
        gloci()
            .args(["maskfasta", "--fasta"])
            .arg(&f)
            .arg("-b")
            .arg(&b)
            .arg(format!("--mask-char={bad}"))
            .assert()
            .code(2)
            .stderr(contains("mask character"));
    }
    gloci()
        .args(["maskfasta", "--fasta"])
        .arg(&f)
        .arg("-b")
        .arg(&b)
        .arg("--mask-char= ")
        .assert()
        .code(2)
        .stderr(contains("space").and(contains("index")));
    gloci()
        .args(["maskfasta", "--fasta"])
        .arg(&f)
        .arg("-b")
        .arg(&b)
        .args(["--soft", "--mask-char", "X"])
        .assert()
        .code(2);
}

/// A line of exactly 1 MiB makes the stream emit its line break as an empty
/// piece; layout and masks must survive that, whatever follows the line.
#[test]
fn empty_piece_at_the_chunk_boundary() {
    const MIB: usize = 1 << 20;
    let d = tempdir();
    let line: Vec<u8> = (0..MIB).map(|i| b"ACGT"[i % 4]).collect();
    for (tail, name) in [
        (&b"\n"[..], "lf"),
        (&b"\r\n"[..], "crlf"),
        (&b""[..], "no final break"),
    ] {
        let mut input = b">big\n".to_vec();
        input.extend_from_slice(&line);
        input.extend_from_slice(tail);
        if !tail.is_empty() {
            input.extend_from_slice(b"ACGT");
            input.extend_from_slice(tail);
        }
        let start = 5 + MIB - 3;
        let mut want = input.clone();
        want[start..start + 3].fill(b'N');
        if tail.is_empty() {
            // The whole sequence is the one line, and the mask ends with it.
            let out = mask_with(d.path(), &input, &format!("big\t{}\t{MIB}\n", MIB - 3), &[]);
            assert_eq!(out, want, "{name}");
            continue;
        }
        // The mask straddles the boundary into the second line.
        let out = mask_with(
            d.path(),
            &input,
            &format!("big\t{}\t{}\n", MIB - 3, MIB + 2),
            &[],
        );
        let second = 5 + MIB + tail.len();
        want[second..second + 2].fill(b'N');
        assert_eq!(out, want, "{name}");
    }
}

#[test]
fn past_the_end_is_an_error_and_publishes_nothing() {
    let d = tempdir();
    let out = d.path().join("out.fa");
    gloci()
        .args(["maskfasta", "--fasta"])
        .arg(fa("ref.fa"))
        .arg("-b")
        .arg(fa("masks-past-end.bed"))
        .arg("-o")
        .arg(&out)
        .assert()
        .code(1)
        .stderr(
            contains("chr1")
                .and(contains("1,010"))
                .and(contains("1,000 bases")),
        );
    assert!(!out.exists());
}

#[test]
fn masks_from_several_files_and_formats_with_chrom_map() {
    let d = tempdir();
    let f = gloci_testkit::fixtures::write(
        d.path(),
        "in.fa",
        b">chr1\nACGTACGTAC\n>chr2\nACGTACGTAC\n",
    );
    let bed = plain_file(d.path(), "a.bed", "chr2\t0\t2\n");
    let gff = plain_file(
        d.path(),
        "b.gff3",
        "##gff-version 3\n1\tsrc\tgene\t9\t10\t.\t+\t.\tID=g\n",
    );
    let map = plain_file(d.path(), "map.tsv", "1\tchr1\n");
    let out = run_ok(&[
        "--fasta",
        f.to_str().unwrap(),
        "-b",
        bed.to_str().unwrap(),
        gff.to_str().unwrap(),
        "--chrom-map",
        map.to_str().unwrap(),
    ]);
    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">chr1\nACGTACGTNN\n>chr2\nNNGTACGTAC\n"
    );
    // A bad mask line is an error unless --skip-invalid.
    let bad = plain_file(d.path(), "bad.bed", "chr1\t0\t1\nchr1\t5\tx\n");
    gloci()
        .args(["maskfasta", "--fasta"])
        .arg(&f)
        .arg("-b")
        .arg(&bad)
        .assert()
        .code(1);
    let out = run_ok(&[
        "--fasta",
        f.to_str().unwrap(),
        "-b",
        bad.to_str().unwrap(),
        "--skip-invalid",
    ]);
    assert_eq!(
        String::from_utf8(out).unwrap(),
        ">chr1\nNCGTACGTAC\n>chr2\nACGTACGTAC\n"
    );
}

#[test]
fn reads_compressed_and_piped_fasta_without_an_index() {
    let plain = run_ok(&[
        "--fasta",
        fa("ref.fa").to_str().unwrap(),
        "-b",
        fa("masks.bed").to_str().unwrap(),
    ]);
    let d = tempdir();
    let text = std::fs::read(fa("ref.fa")).unwrap();
    let gz = gloci_testkit::fixtures::write(
        d.path(),
        "ref.fa.gz",
        &gloci_testkit::fixtures::gzip(&text),
    );
    let bgz = gloci_testkit::fixtures::write(
        d.path(),
        "ref2.fa.gz",
        &gloci_testkit::fixtures::bgzip(&text),
    );
    for p in [gz, bgz] {
        let got = run_ok(&[
            "--fasta",
            p.to_str().unwrap(),
            "-b",
            fa("masks.bed").to_str().unwrap(),
        ]);
        assert_eq!(got, plain, "{}", p.display());
    }
    let out = gloci()
        .args(["maskfasta", "--fasta", "-", "-b"])
        .arg(fa("masks.bed"))
        .write_stdin(text)
        .output()
        .unwrap();
    assert_eq!(out.stdout, plain);
}

#[test]
fn write_index_equals_gloci_index_of_the_output() {
    for (name, gz) in [("out.fa.gz", true), ("out.fa", false)] {
        let d = tempdir();
        let out = d.path().join(name);
        gloci()
            .args(["maskfasta", "--fasta"])
            .arg(fa("ref.fa.gz"))
            .arg("-b")
            .arg(fa("masks.bed"))
            .arg("-o")
            .arg(&out)
            .arg("--write-index")
            .assert()
            .success();
        let fai = std::fs::read(format!("{}.fai", out.display())).unwrap();
        let gzi = std::fs::read(format!("{}.gzi", out.display())).ok();
        assert_eq!(gzi.is_some(), gz, "{name}");
        // The reference: index a copy of the output with `gloci index`.
        let copy = d.path().join(format!("copy-{name}"));
        std::fs::copy(&out, &copy).unwrap();
        gloci().arg("index").arg(&copy).assert().success();
        assert_eq!(
            fai,
            std::fs::read(format!("{}.fai", copy.display())).unwrap(),
            "{name}"
        );
        if gz {
            assert_eq!(
                gzi.unwrap(),
                std::fs::read(format!("{}.gzi", copy.display())).unwrap()
            );
        }
        assert!(!fai.is_empty());
    }
    // Without -o there is nothing to index.
    gloci()
        .args(["maskfasta", "--fasta"])
        .arg(fa("ref.fa"))
        .arg("-b")
        .arg(fa("masks.bed"))
        .arg("--write-index")
        .assert()
        .code(2)
        .stderr(contains("--write-index needs -o FILE"));
}

#[test]
fn output_is_identical_at_any_thread_count() {
    let d = tempdir();
    let mut plain: Vec<Vec<u8>> = Vec::new();
    let mut bgzf: Vec<Vec<u8>> = Vec::new();
    for t in ["1", "2", "8"] {
        let p = d.path().join(format!("o{t}.fa.gz"));
        gloci()
            .args(["maskfasta", "-t", t, "--fasta"])
            .arg(fa("ref.fa.gz"))
            .arg("-b")
            .arg(fa("masks.bed"))
            .arg("-o")
            .arg(&p)
            .assert()
            .success();
        bgzf.push(std::fs::read(&p).unwrap());
        plain.push(
            gloci()
                .args(["maskfasta", "-t", t, "--fasta"])
                .arg(fa("ref.fa.gz"))
                .arg("-b")
                .arg(fa("masks.bed"))
                .output()
                .unwrap()
                .stdout,
        );
    }
    assert!(plain.windows(2).all(|w| w[0] == w[1]));
    assert!(bgzf.windows(2).all(|w| w[0] == w[1]));
    assert!(plain[0].contains(&b'N'));
}

/// Peak RSS of a finished child in bytes, from its rusage.
fn wait_rss(child: &std::process::Child) -> (bool, u64) {
    let pid = child.id() as libc::pid_t;
    let mut status = 0;
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    let r = unsafe { libc::wait4(pid, &mut status, 0, &mut ru) };
    assert_eq!(r, pid, "wait4 failed");
    let v = ru.ru_maxrss as u64;
    (
        libc::WIFEXITED(status) && libc::WEXITSTATUS(status) == 0,
        if cfg!(target_os = "macos") {
            v
        } else {
            v * 1024
        },
    )
}

/// A 50 MiB FASTA as one sequence line (the worst case for a line reader)
/// is masked within a few MiB more than a bare process needs.
#[test]
fn maskfasta_streams() {
    const MIB: usize = 1 << 20;
    let d = tempdir();
    let fasta = d.path().join("big.fa");
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&fasta).unwrap());
        w.write_all(b">big\n").unwrap();
        let block: Vec<u8> = (0..MIB).map(|i| b"ACGT"[i % 4]).collect();
        for _ in 0..50 {
            w.write_all(&block).unwrap();
        }
        w.write_all(b"\n").unwrap();
    }
    let bed = plain_file(
        d.path(),
        "m.bed",
        &format!(
            "big\t0\t10\nbig\t{}\t{}\nbig\t{}\t{}\n",
            MIB - 5,
            MIB + 5,
            50 * MIB - 10,
            50 * MIB
        ),
    );
    let out = d.path().join("out.fa");
    let child = std::process::Command::new(assert_cmd::cargo::cargo_bin("gloci"))
        .args(["maskfasta", "-t", "1", "--fasta"])
        .arg(&fasta)
        .arg("-b")
        .arg(&bed)
        .arg("-o")
        .arg(&out)
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    assert!(ok, "gloci maskfasta failed");
    // The same process reading a 1 MiB FASTA sets the baseline.
    let baseline = {
        let small = d.path().join("small.fa");
        std::fs::write(&small, b">big\nACGTACGTAC\n").unwrap();
        let sb = plain_file(d.path(), "s.bed", "big\t0\t1\n");
        let child = std::process::Command::new(assert_cmd::cargo::cargo_bin("gloci"))
            .args(["maskfasta", "-t", "1", "--fasta"])
            .arg(&small)
            .arg("-b")
            .arg(&sb)
            .stdout(std::process::Stdio::null())
            .spawn()
            .unwrap();
        wait_rss(&child).1
    };
    assert!(
        rss <= baseline + 12 * MIB as u64,
        "peak RSS {} MiB for a 50 MiB FASTA (baseline {} MiB)",
        rss / MIB as u64,
        baseline / MIB as u64
    );
    // And the result is right: only the masked bases differ.
    let got = std::fs::read(&out).unwrap();
    assert_eq!(got.len(), 5 + 50 * MIB + 1);
    let seq = &got[5..5 + 50 * MIB];
    assert!(seq[..10].iter().all(|&b| b == b'N'));
    assert!(seq[MIB - 5..MIB + 5].iter().all(|&b| b == b'N'));
    assert!(seq[50 * MIB - 10..].iter().all(|&b| b == b'N'));
    assert_eq!(seq[10], b"ACGT"[10 % 4]);
    assert_eq!(seq[MIB + 5], b"ACGT"[(MIB + 5) % 4]);
    assert_eq!(seq.iter().filter(|&&b| b == b'N').count(), 30);
}

#[test]
fn a_replaced_fasta_output_loses_its_stale_index() {
    let d = tempdir();
    let masks = plain_file(d.path(), "m.bed", "a\t0\t2\n");
    for (name, gz) in [("o.fa", false), ("o.fa.gz", true)] {
        let out = d.path().join(name);
        let fai = d.path().join(format!("{name}.fai"));
        let gzi = d.path().join(format!("{name}.gzi"));
        let run = |fasta: &Path, extra: &[&str]| {
            gloci()
                .args(["maskfasta", "--fasta"])
                .arg(fasta)
                .arg("-b")
                .arg(&masks)
                .arg("-o")
                .arg(&out)
                .args(extra)
                .assert()
                .success();
        };
        let first = gloci_testkit::fixtures::write(d.path(), "w10.fa", b">a\nACGTACGTAC\nACGT\n");
        run(&first, &["--write-index"]);
        assert!(fai.exists(), "{name}");
        assert_eq!(gzi.exists(), gz, "{name}");
        // Rewritten with other line widths and no --write-index: the old index describes other data.
        let second = gloci_testkit::fixtures::write(d.path(), "w4.fa", b">a\nACGT\nACGT\nAC\n");
        run(&second, &[]);
        assert!(!fai.exists(), "{name}: stale .fai kept");
        assert!(!gzi.exists(), "{name}: stale .gzi kept");
        assert!(out.exists());
    }
}

#[test]
fn skip_invalid_clips_masks_past_the_end() {
    let d = tempdir();
    let input = b">a\nACGTAC\nGT\n>b\nACGT\n";
    // a is 8 bases: 6-12 is clipped to 6-8, 8-9 and 20-30 start at or past the end and mask
    // nothing; b is 4 bases and 3-5 is clipped.
    let bed = "a\t6\t12\na\t8\t9\na\t20\t30\nb\t3\t5\n";
    let f = gloci_testkit::fixtures::write(d.path(), "in.fa", input);
    let b = plain_file(d.path(), "m.bed", bed);
    let args = [
        "maskfasta",
        "--fasta",
        f.to_str().unwrap(),
        "-b",
        b.to_str().unwrap(),
    ];
    gloci()
        .args(args)
        .assert()
        .code(1)
        .stderr(contains("--skip-invalid"));
    let out = gloci().args(args).arg("--skip-invalid").output().unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        ">a\nACGTAC\nNN\n>b\nACGN\n"
    );
    let err = String::from_utf8(out.stderr).unwrap();
    assert_eq!(
        err.lines().filter(|l| l.contains("past the end")).count(),
        1,
        "{err}"
    );
    assert!(
        err.contains("3 masks") && err.contains("first: a mask on a covers 6-12"),
        "{err}"
    );
}

#[test]
fn an_empty_fasta_is_named_as_empty() {
    let d = tempdir();
    let empty = plain_file(d.path(), "empty.fa", "");
    let b = plain_file(d.path(), "m.bed", "a\t0\t2\n");
    gloci()
        .args(["maskfasta", "--fasta"])
        .arg(&empty)
        .arg("-b")
        .arg(&b)
        .assert()
        .code(1)
        .stderr(contains("empty.fa is empty: it has no records"));
    // An empty mask file stays valid and masks nothing.
    let f = gloci_testkit::fixtures::write(d.path(), "in.fa", b">a\nACGT\n");
    let eb = plain_file(d.path(), "e.bed", "");
    let out = run_ok(&["--fasta", f.to_str().unwrap(), "-b", eb.to_str().unwrap()]);
    assert_eq!(out, b">a\nACGT\n");
}
