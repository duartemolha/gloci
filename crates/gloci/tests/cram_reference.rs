//! CRAM on the command line (D109, D110): `view x.cram --reference ref.fa`, the reference
//! errors (missing, mismatched, BGZF without .gzi), truncated files, `gloci index x.cram`,
//! `-g x.cram`, and the warning for `--reference` with no CRAM input. Every failure is exit 1
//! or 2, never a crash (101, 134).
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn aln(name: &str) -> String {
    s(&inputs().join("aln").join(name))
}

fn reference() -> String {
    s(&inputs().join("fa/ref.fa"))
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn text(b: &[u8]) -> String {
    String::from_utf8(b.to_vec()).unwrap()
}

/// The record lines of SAM text.
fn records(sam: &[u8]) -> Vec<String> {
    text(sam)
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_string)
        .collect()
}

fn stdout_of(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

/// gloci fails with `code` and a message containing each of `want`.
fn fails(args: &[&str], code: i32, want: &[&str]) {
    let out = gloci().args(args).output().unwrap();
    let err = text(&out.stderr);
    assert_eq!(out.status.code(), Some(code), "{args:?}: {err}");
    for w in want {
        assert!(err.contains(w), "{args:?}: no {w:?} in {err}");
    }
}

/// ref.fa rewritten by `edit`, indexed, in `dir`.
fn edited_reference(dir: &Path, edit: impl Fn(&str) -> String) -> String {
    let fa = dir.join("ref.fa");
    let text = std::fs::read_to_string(reference()).unwrap();
    std::fs::write(&fa, edit(&text)).unwrap();
    gloci().args(["index", &s(&fa)]).assert().success();
    s(&fa)
}

#[test]
fn view_cram_with_reference() {
    let r = reference();
    for set in ["flags", "tags", "pastend"] {
        let want = records(&stdout_of(&[
            "view",
            "-O",
            "sam",
            &aln(&format!("{set}.bam")),
        ]));
        for v in ["", ".v31", ".multi", ".embed"] {
            let cram = aln(&format!("{set}{v}.cram"));
            if !Path::new(&cram).exists() {
                continue;
            }
            let one = stdout_of(&["view", "-O", "sam", &cram, "--reference", &r]);
            assert_eq!(records(&one), want, "{cram}");
            let header = text(&one);
            assert!(
                header.contains("\tM5:") || v == ".noref",
                "{cram}: {header}"
            );
            assert!(header.contains("@PG\tID:gloci"), "{cram}");
            for t in ["2", "8"] {
                let got = stdout_of(&["view", "-t", t, "-O", "sam", &cram, "--reference", &r]);
                assert_eq!(got, one, "{cram} -t {t}");
            }
        }
    }
    // BAM by default, and to a file.
    let dir = tempdir();
    let out = dir.path().join("o.bam");
    gloci()
        .args([
            "view",
            &aln("flags.cram"),
            "--reference",
            &r,
            "-o",
            &s(&out),
        ])
        .assert()
        .success();
    let back = stdout_of(&["view", "-O", "sam", &s(&out)]);
    assert_eq!(
        records(&back),
        records(&stdout_of(&["view", "-O", "sam", &aln("flags.bam")]))
    );
    // Regions by scan, with the read filters, and the -v counts of a BAM scanned likewise (a
    // copy without its .bai: through the index, fewer reads are visited and counted).
    let region = ["chr1:650-701", "--exclude-flags", "0x704"];
    let scanned = dir.path().join("flags.bam");
    std::fs::copy(aln("flags.bam"), &scanned).unwrap();
    let bam = gloci()
        .args(["-v", "view", "-O", "sam", &s(&scanned)])
        .args(region)
        .output()
        .unwrap();
    let cram = gloci()
        .args([
            "-v",
            "view",
            "-O",
            "sam",
            &aln("flags.cram"),
            "--reference",
            &r,
        ])
        .args(region)
        .output()
        .unwrap();
    assert!(cram.status.success());
    assert_eq!(records(&cram.stdout), records(&bam.stdout));
    let dropped = |o: &std::process::Output| {
        text(&o.stderr)
            .lines()
            .filter(|l| l.contains("reads dropped by"))
            .map(|l| l.split(": ").last().unwrap().to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(dropped(&cram), dropped(&bam));
    assert!(!dropped(&cram).is_empty());
}

#[test]
fn cram_no_ref_needed() {
    let want = records(&stdout_of(&["view", "-O", "sam", &aln("tags.bam")]));
    for v in ["tags.noref.cram", "tags.embed.cram"] {
        assert_eq!(
            records(&stdout_of(&["view", "-O", "sam", &aln(v)])),
            want,
            "{v}"
        );
    }
}

#[test]
fn cram_missing_reference() {
    for t in ["1", "4"] {
        fails(
            &["view", "-t", t, "-O", "sam", &aln("flags.cram")],
            1,
            &[
                "flags.cram needs the reference sequence of chr1 (M5 ",
                "pass it with --reference FASTA",
            ],
        );
    }
    // Written to a file: nothing is left behind.
    let dir = tempdir();
    let out = dir.path().join("o.bam");
    fails(&["view", &aln("tags.cram"), "-o", &s(&out)], 1, &["chrM"]);
    assert!(!out.exists());
}

#[test]
fn cram_mismatched_reference() {
    let dir = tempdir();
    let changed = edited_reference(dir.path(), |t| {
        // chr1:101, the first base of the first reads of flags.cram.
        let i = t.find(">chr1\n").unwrap() + 6 + 61 + 40;
        let mut b = t.as_bytes().to_vec();
        b[i] = if b[i] == b'A' { b'C' } else { b'A' };
        String::from_utf8(b).unwrap()
    });
    fails(
        &[
            "view",
            "-O",
            "sam",
            &aln("flags.cram"),
            "--reference",
            &changed,
        ],
        1,
        &[
            "the sequence chr1 in",
            "does not match the one",
            "(chr1:101-701)",
        ],
    );
    let dir = tempdir();
    let longer = edited_reference(dir.path(), |t| t.replace(">chrX\n", "A\n>chrX\n"));
    // chr10 is now one base longer: refused before any record is read.
    let out = gloci()
        .args([
            "view",
            "-O",
            "sam",
            &aln("tags.cram"),
            "--reference",
            &longer,
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1), "{}", text(&out.stderr));
    assert!(out.stdout.is_empty());
    let err = text(&out.stderr);
    assert!(
        err.contains("tags.cram was not made with") && err.contains("chr10 is 600 bases long"),
        "{err}"
    );
}

#[test]
fn cram_reference_bgzf() {
    let dir = tempdir();
    let gz = dir.path().join("ref.fa.gz");
    std::fs::write(
        &gz,
        gloci_testkit::fixtures::bgzip(&std::fs::read(reference()).unwrap()),
    )
    .unwrap();
    gloci().args(["index", &s(&gz)]).assert().success();
    let want = records(&stdout_of(&["view", "-O", "sam", &aln("flags.bam")]));
    let got = stdout_of(&[
        "view",
        "-O",
        "sam",
        &aln("flags.cram"),
        "--reference",
        &s(&gz),
    ]);
    assert_eq!(records(&got), want);
    std::fs::remove_file(dir.path().join("ref.fa.gz.gzi")).unwrap();
    fails(
        &[
            "view",
            "-O",
            "sam",
            &aln("flags.cram"),
            "--reference",
            &s(&gz),
        ],
        2,
        &["is BGZF-compressed and has no index", "gloci index"],
    );
}

#[test]
fn cram_truncated() {
    let bytes = std::fs::read(aln("flags.cram")).unwrap();
    let dir = tempdir();
    for cut in [100usize, 400, 1105, 1300, bytes.len() - 38, bytes.len() - 1] {
        let path = dir.path().join(format!("cut{cut}.cram"));
        std::fs::write(&path, &bytes[..cut]).unwrap();
        fails(
            &["view", "-O", "sam", &s(&path), "--reference", &reference()],
            1,
            &[&s(&path), "probably truncated"],
        );
        // From standard input, where the end cannot be checked first.
        let out = gloci()
            .args(["view", "-O", "sam", "-", "--reference", &reference()])
            .write_stdin(bytes[..cut].to_vec())
            .output()
            .unwrap();
        assert_eq!(
            out.status.code(),
            Some(1),
            "cut {cut}: {}",
            text(&out.stderr)
        );
        assert!(text(&out.stderr).contains("stdin"), "{}", text(&out.stderr));
    }
    // Garbage after the magic.
    let path = dir.path().join("junk.cram");
    let mut junk = b"CRAM\x03\x00".to_vec();
    junk.extend(std::iter::repeat_n(0xa5u8, 200));
    std::fs::write(&path, junk).unwrap();
    fails(&["view", "-O", "sam", &s(&path)], 1, &["junk.cram"]);
    // CRAM 2.1 is refused by its version, before its end is checked.
    let path = dir.path().join("v2.cram");
    let mut v2 = b"CRAM\x02\x01".to_vec();
    v2.extend(std::iter::repeat_n(0u8, 100));
    std::fs::write(&path, v2).unwrap();
    fails(
        &["view", "-O", "sam", &s(&path)],
        1,
        &[
            "v2.cram is CRAM 2.1, which gloci does not read",
            "samtools view -C",
        ],
    );
}

#[test]
fn index_cram_cli() {
    let dir = tempdir();
    let cram = dir.path().join("multi.cram");
    std::fs::copy(aln("tags.multi.cram"), &cram).unwrap();
    // Slices of several references need no reference to index (samtools needs one).
    gloci().args(["index", &s(&cram)]).assert().success();
    let crai = dir.path().join("multi.cram.crai");
    let read = |p: &Path| gloci_io::CramIndex::read(p).unwrap();
    let samtools = read(&inputs().join("aln/tags.multi.cram.crai"));
    assert_eq!(read(&crai), samtools);
    std::fs::remove_file(&crai).unwrap();
    gloci()
        .args(["index", "-t", "4", &s(&cram), "--reference", &reference()])
        .assert()
        .success();
    assert_eq!(read(&crai), samtools);
    // A reference of other lengths is refused.
    let other = edited_reference(dir.path(), |t| t.replace(">chrX\n", "A\n>chrX\n"));
    fails(
        &["index", &s(&cram), "--reference", &other],
        1,
        &["multi.cram was not made with"],
    );
    // A truncated CRAM is not indexed.
    let cut = dir.path().join("cut.cram");
    let bytes = std::fs::read(&cram).unwrap();
    std::fs::write(&cut, &bytes[..bytes.len() - 10]).unwrap();
    fails(&["index", &s(&cut)], 1, &["cut.cram", "probably truncated"]);
    assert!(!dir.path().join("cut.cram.crai").exists());
}

/// gloci index refuses a CRAM unsorted across slices (samtools index refuses it too), and a read
/// far past its reference's end is an error, not a crash.
#[test]
fn unsorted_and_overhang_are_errors() {
    let dir = tempdir();
    let cram = dir.path().join("u.cram");
    std::fs::copy(aln("unsorted.cram"), &cram).unwrap();
    fails(
        &["index", &s(&cram)],
        1,
        &["u.cram is not sorted by position", "gloci sort"],
    );
    assert!(!dir.path().join("u.cram.crai").exists());
    fails(
        &[
            "view",
            "-O",
            "sam",
            &aln("overhang.cram"),
            "--reference",
            &reference(),
        ],
        1,
        &[
            "past the end of chr1",
            "more than gloci's limit of 1,048,576",
        ],
    );
}

/// CRAM 3.1 with order-1 rANS and the name tokeniser decodes on every worker (16 MiB stacks).
#[test]
fn cram_31_codecs_at_every_thread_count() {
    let one = stdout_of(&["view", "-O", "sam", &aln("cat-big.v31.cram")]);
    assert_eq!(records(&one).len(), 1300);
    for t in ["2", "8"] {
        let got = stdout_of(&["view", "-t", t, "-O", "sam", &aln("cat-big.v31.cram")]);
        assert_eq!(got, one, "-t {t}");
    }
}

/// A crafted CRAM (a block of flags.noref.cram changed and its checksum recomputed) that
/// reaches a panic inside noodles-cram is reported as malformed data, pointing at #cram-limits,
/// not as a bug of gloci's (exit 1 still), on the main thread and on a worker.
#[test]
fn crafted_cram_panic_points_to_the_limits() {
    let mut b = std::fs::read(aln("flags.noref.cram")).unwrap();
    // The external block at byte 743: method, content type, then ITF8 id, size and raw size.
    let start = 743;
    let itf8 = |at: &mut usize, b: &[u8]| -> usize {
        let v = b[*at] as usize;
        assert!(v < 0x80, "a one-byte ITF8 in this block");
        *at += 1;
        v
    };
    let mut at = start + 2;
    itf8(&mut at, &b);
    let size = itf8(&mut at, &b);
    itf8(&mut at, &b);
    let end = at + size;
    b[748] ^= 1;
    // CRC32 (IEEE), bit by bit.
    let mut crc = !0u32;
    for &byte in &b[start..end] {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0xedb8_8320
            } else {
                crc >> 1
            };
        }
    }
    b[end..end + 4].copy_from_slice(&(!crc).to_le_bytes());
    let dir = tempdir();
    let path = dir.path().join("crafted.cram");
    std::fs::write(&path, b).unwrap();
    for t in ["1", "4"] {
        let out = gloci()
            .args(["view", "-t", t, "-O", "sam", &s(&path)])
            .output()
            .unwrap();
        let err = text(&out.stderr);
        assert_eq!(out.status.code(), Some(1), "-t {t}: {err}");
        assert!(
            err.contains("crafted.cram: the CRAM decoder failed on malformed data in container 1")
                && err.contains("see #cram-limits in docs/differences.md")
                && !err.contains("this is a bug"),
            "-t {t}: {err}"
        );
    }
}

#[test]
fn genome_from_cram_needs_no_reference() {
    let dir = tempdir();
    let bed = plain_file(
        dir.path(),
        "a.bed",
        "chr10\t1\t2\nchr1\t1\t2\nchrX\t1\t2\nchrM\t1\t2\nchr2\t1\t2\n",
    );
    let out = stdout_of(&["sort", &s(&bed), "-g", &aln("tags.cram")]);
    assert_eq!(
        text(&out),
        "chrM\t1\t2\nchr1\t1\t2\nchr2\t1\t2\nchr10\t1\t2\nchrX\t1\t2\n"
    );
}

#[test]
fn reference_without_cram_warns() {
    gloci()
        .args([
            "view",
            "-O",
            "sam",
            &aln("flags.bam"),
            "--reference",
            &reference(),
        ])
        .assert()
        .success()
        .stderr(contains("--reference").and(contains("is not used: no input is CRAM")));
    let dir = tempdir();
    let bed = bgzf_file(dir.path(), "a.bed.gz", "chr1\t1\t2\n");
    gloci()
        .args(["index", &s(&bed), "--reference", &reference()])
        .assert()
        .success()
        .stderr(contains("is not used: no input is CRAM"));
    // With a CRAM input there is no such warning.
    gloci()
        .args([
            "view",
            "-O",
            "sam",
            &aln("flags.cram"),
            "--reference",
            &reference(),
        ])
        .assert()
        .success()
        .stderr(contains("is not used").not());
}

#[test]
fn interval_commands_do_not_read_cram_yet() {
    // Only view reads CRAM so far; the interval commands keep the "not yet" answer.
    let dir = tempdir();
    let bed = plain_file(dir.path(), "a.bed", "chr1\t1\t100\n");
    let (b, c) = (s(&bed), aln("flags.cram"));
    for args in [
        vec!["intersect", "-a", &c, "-b", &b],
        vec!["intersect", "-a", &b, "-b", &c],
        vec!["coverage", "-a", &b, "-b", &c],
        vec!["genomecov", &c],
        vec!["convert", &c],
        vec!["sort", &c],
    ] {
        gloci().args(&args).assert().failure().stderr(contains(
            "is CRAM, and reading it is not available in this build yet",
        ));
    }
}
