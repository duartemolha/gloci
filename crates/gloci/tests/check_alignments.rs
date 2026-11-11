//! `gloci check` of BAM, CRAM and SAM (D129) and `check --quick`, which replaces samtools
//! quickcheck. The comparison with samtools is in the goldens (`check-quick-*`).
mod common;
use common::*;
use gloci_ops::{QuickOptions, quick_check_reader};
use predicates::str::contains;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn sam_fixture(name: &str) -> String {
    s(&inputs().join("samtools").join(name))
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[&str]) -> (i32, String, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code().unwrap(),
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
    )
}

/// The decompressed bytes of a BAM and the offset of its first record's `block_size`.
fn bam_bytes(path: &Path) -> (Vec<u8>, usize) {
    let mut raw = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut raw)
        .unwrap();
    let u32_at = |at: usize| u32::from_le_bytes(raw[at..at + 4].try_into().unwrap()) as usize;
    let l_text = u32_at(4);
    let mut at = 8 + l_text;
    let n_ref = u32_at(at);
    at += 4;
    for _ in 0..n_ref {
        at += 4 + u32_at(at) + 4;
    }
    (raw, at)
}

fn bgzf_of(bytes: &[u8]) -> Vec<u8> {
    let mut w = noodles_bgzf::io::Writer::new(Vec::new());
    w.write_all(bytes).unwrap();
    w.finish().unwrap()
}

/// Writes `xtask/diff/inputs/samtools/badfirst.bam`: `flagstat.bam` with the `block_size` of its
/// first record set to 20 (less than the 32 bytes of fixed fields), under a fine header and
/// a fine end-of-file block. Run: `cargo test -p gloci --test check_alignments make_broken_bams
/// -- --ignored`.
#[test]
#[ignore]
fn make_broken_bams() {
    let src = inputs().join("samtools/flagstat.bam");
    let (mut raw, at) = bam_bytes(&src);
    raw[at..at + 4].copy_from_slice(&20u32.to_le_bytes());
    std::fs::write(inputs().join("samtools/badfirst.bam"), bgzf_of(&raw)).unwrap();
}

#[test]
fn badfirst_fixture_is_what_it_says() {
    let p = inputs().join("samtools/badfirst.bam");
    let (raw, at) = bam_bytes(&p);
    assert_eq!(&raw[at..at + 4], &20u32.to_le_bytes());
    let bytes = std::fs::read(&p).unwrap();
    assert_eq!(&bytes[bytes.len() - 28..], &gloci_io::BGZF_EOF);
}

fn sam_file(dir: &Path, name: &str, text: &str) -> String {
    s(&plain_file(dir, name, text))
}

#[test]
fn check_alignment_full() {
    // Out of order under SO:coordinate: an error that names the record, the reference and both
    // positions.
    for f in ["coord-unsorted.bam", "coord-unsorted.sam"] {
        let (code, out, _) = run(&["check", &sam_fixture(f)]);
        assert_eq!(code, 1, "{f}: {out}");
        assert!(
            out.contains("record 7 of chr1: position 150 is before the previous position 600"),
            "{f}: {out}"
        );
        assert!(out.contains("NOT sorted"), "{out}");
        assert!(out.contains("errors:       1\n"), "{out}");
        assert!(out.contains("mapped:       10"), "{out}");
    }
    // The same reads under another SO: reported, not an error.
    let d = tempdir();
    let text = std::fs::read_to_string(sam_fixture("coord-unsorted.sam"))
        .unwrap()
        .replace("SO:coordinate", "SO:queryname");
    let f = sam_file(d.path(), "q.sam", &text);
    let (code, out, _) = run(&["check", &f]);
    assert_eq!(code, 0, "{out}");
    assert!(out.contains("sorted:       no (@HD SO:queryname"), "{out}");
    // A read past LN is a warning.
    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
a\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
b\t0\tchr1\t998\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n";
    let f = sam_file(d.path(), "past.sam", text);
    let (code, out, _) = run(&["check", &f]);
    assert_eq!(code, 0, "{out}");
    assert!(out.contains("warnings:     1\n"), "{out}");
    assert!(
        out.contains("record 2 of chr1: the alignment at position 998 runs past the end of chr1"),
        "{out}"
    );
    // Identical records at one position are a warning too.
    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
a\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
a\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
a\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTT\tIIIII\n";
    let f = sam_file(d.path(), "dup.sam", text);
    let (code, out, _) = run(&["check", &f]);
    assert_eq!(code, 0, "{out}");
    assert!(out.contains("warnings:     1\n"), "{out}");
    assert!(out.contains("record 2 of chr1: duplicate"), "{out}");
}

#[test]
fn check_alignment_summary_counts() {
    let d = tempdir();
    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n\
a\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
b\t16\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
c\t4\tchr1\t30\t0\t*\t*\t0\t0\tACGTA\tIIIII\n\
d\t0\tchr2\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
e\t4\tchr2\t10\t0\t*\t*\t0\t0\tACGTA\tIIIII\n\
f\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n\
g\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n";
    let f = sam_file(d.path(), "c.sam", text);
    let (code, out, _) = run(&["check", &f]);
    assert_eq!(code, 0, "{out}");
    for want in [
        "records:      7\n",
        "mapped:       3\n",
        "unmapped placed: 2\n",
        "unplaced:     2\n",
        "chromosomes:  2 (chr1, chr2)\n",
        "sorted:       yes (@HD SO:coordinate)\n",
    ] {
        assert!(out.contains(want), "{want:?} in {out}");
    }
    // A placed read after an unplaced one.
    let bad = text
        .replace("g\t4\t*\t0\t0", "g\t0\tchr2\t900\t60")
        .replace("g\t0\tchr2\t900\t60\t*", "g\t0\tchr2\t900\t60\t5M");
    let f = sam_file(d.path(), "u.sam", &bad);
    let (code, out, _) = run(&["check", &f]);
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("record 7 of chr2: a placed read follows an unplaced one"),
        "{out}"
    );
}

#[test]
fn check_alignment_files_and_threads() {
    // A BAM without its end-of-file block, and a truncated one, fail with one summary each.
    for f in ["noeof.bam", "truncated.bam"] {
        let (code, out, _) = run(&["check", &sam_fixture(f)]);
        assert_eq!(code, 1, "{f}: {out}");
        assert!(out.contains("cannot be checked"), "{f}: {out}");
    }
    // A CRAM without its end-of-file container fails after its records are read.
    let cram = s(&inputs().join("aln/flags-noeof.cram"));
    let fa = s(&inputs().join("fa/ref.fa"));
    let (code, out, _) = run(&["check", &cram, "--reference", &fa]);
    assert_eq!(code, 1, "{out}");
    assert!(out.contains("end-of-file container is missing"), "{out}");
    let cram = s(&inputs().join("aln/flags.cram"));
    let (code, out, _) = run(&["check", &cram, "--reference", &fa]);
    assert_eq!(code, 0, "{out}");
    assert!(out.contains("CRAM\n"), "{out}");
    // The same output at any thread count.
    let bam = sam_fixture("coord-unsorted.bam");
    let (_, one, _) = run(&["check", &bam, "-t", "1"]);
    for t in ["2", "8"] {
        let (_, many, _) = run(&["check", &bam, "-t", t]);
        assert_eq!(one, many, "-t {t}");
    }
}

fn quick_status(args: &[&str]) -> (i32, String) {
    let mut a = vec!["check", "--quick", "--list-failed"];
    a.extend_from_slice(args);
    let (code, out, _) = run(&a);
    (code, out)
}

#[test]
fn quick_check_matrix() {
    let d = tempdir();
    let good_bam = sam_fixture("flagstat.bam");
    let good_sam = sam_fixture("flagstat.sam");
    let good_cram = s(&inputs().join("aln/flags.cram"));
    let good_vcf = sam_fixture("calls.vcf.gz");
    for f in [&good_bam, &good_sam, &good_cram, &good_vcf] {
        assert_eq!(quick_status(&[f]), (0, String::new()), "{f}");
    }
    // The summary form, without --list-failed.
    let (code, out, _) = run(&["check", "--quick", &good_bam]);
    assert_eq!(code, 0);
    assert!(out.contains("quick check:  passed"), "{out}");
    assert!(out.ends_with("result: OK\n"), "{out}");
    let missing = s(&d.path().join("nope.bam"));
    for f in [
        sam_fixture("noeof.bam"),
        sam_fixture("truncated.bam"),
        sam_fixture("badfirst.bam"),
        missing,
        s(&inputs().join("aln/flags-noeof.cram")),
    ] {
        assert_eq!(quick_status(&[&f]), (1, format!("{f}\n")), "{f}");
    }
    // No @SQ: a failure, unless --unaligned.
    for f in ["nosq.bam", "nosq.sam"] {
        let f = sam_fixture(f);
        assert_eq!(quick_status(&[&f]), (1, format!("{f}\n")), "{f}");
        assert_eq!(
            quick_status(&["--unaligned", &f]),
            (0, String::new()),
            "{f}"
        );
        assert_eq!(quick_status(&["-u", &f]), (0, String::new()), "{f}");
    }
    // A CRAM that needs a reference passes without one.
    let (code, out, err) = run(&["check", "--quick", &good_cram]);
    assert_eq!(code, 0, "{out}{err}");
    // A CRAM cut inside its first container fails.
    let bytes = std::fs::read(&good_cram).unwrap();
    let cut = d.path().join("cut.cram");
    std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
    let cut = s(&cut);
    assert_eq!(quick_status(&[&cut]), (1, format!("{cut}\n")));
    // Standard input is not checked for the end-of-file block.
    let noeof = std::fs::read(sam_fixture("noeof.bam")).unwrap();
    let o = gloci()
        .args(["check", "--quick", "--list-failed", "-"])
        .write_stdin(noeof)
        .assert()
        .success();
    assert!(o.get_output().stdout.is_empty());
    // The problem is on stderr, in words.
    let (_, _, err) = run(&[
        "check",
        "--quick",
        "--list-failed",
        &sam_fixture("badfirst.bam"),
    ]);
    assert!(err.contains("badfirst.bam"), "{err}");
}

#[test]
fn list_failed_only_paths() {
    let good = sam_fixture("flagstat.bam");
    let bad1 = sam_fixture("noeof.bam");
    let bad2 = sam_fixture("nosq.bam");
    let (code, out) = quick_status(&[&bad1, &good, &bad2]);
    assert_eq!(code, 1);
    assert_eq!(out, format!("{bad1}\n{bad2}\n"));
    // All good: nothing at all.
    assert_eq!(quick_status(&[&good, &good]), (0, String::new()));
    // The full check takes --list-failed too.
    let (code, out, _) = run(&[
        "check",
        "--list-failed",
        &good,
        &sam_fixture("coord-unsorted.sam"),
    ]);
    assert_eq!(code, 1);
    assert_eq!(out, format!("{}\n", sam_fixture("coord-unsorted.sam")));
    // -u only goes with --quick.
    gloci()
        .args(["check", "-u", &good])
        .assert()
        .code(2)
        .stderr(contains("--quick"));
}

/// A reader that counts the bytes it hands out.
struct Counting<R> {
    inner: R,
    n: Arc<AtomicU64>,
}

impl<R: Read> Read for Counting<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let k = self.inner.read(buf)?;
        self.n.fetch_add(k as u64, Ordering::Relaxed);
        Ok(k)
    }
}

impl<R: Seek> Seek for Counting<R> {
    fn seek(&mut self, to: SeekFrom) -> std::io::Result<u64> {
        self.inner.seek(to)
    }
}

#[test]
fn quick_reads_little() {
    // flagstat.bam's blocks, then a long run of padding (more blocks), then the end-of-file
    // block: quick check must read the header, the first block and the last 28 bytes, not the
    // 12 MB between.
    let d = tempdir();
    let mut bytes = std::fs::read(sam_fixture("flagstat.bam")).unwrap();
    bytes.truncate(bytes.len() - 28);
    // Incompressible, so the padding really is 12 MB of file.
    let mut x = 0x9E37_79B9_7F4A_7C15u64;
    let pad: Vec<u8> = (0..12 << 20)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 24) as u8
        })
        .collect();
    let padded = bgzf_of(&pad);
    bytes.extend_from_slice(&padded[..padded.len() - 28]);
    bytes.extend_from_slice(&gloci_io::BGZF_EOF);
    let path = d.path().join("big.bam");
    std::fs::write(&path, &bytes).unwrap();
    let n = Arc::new(AtomicU64::new(0));
    let r = Counting {
        inner: std::fs::File::open(&path).unwrap(),
        n: n.clone(),
    };
    let ro = gloci_io::ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: std::num::NonZeroUsize::new(4).unwrap(),
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: gloci_io::InvalidPolicy::Fail,
    };
    let report = quick_check_reader("big.bam", r, &ro, &QuickOptions::default());
    assert!(!report.has_errors(), "{}", report.render());
    let read = n.load(Ordering::Relaxed);
    assert!(bytes.len() > 12 << 20, "{}", bytes.len());
    assert!(
        read < 200 * 1024,
        "read {read} bytes of a {} byte file",
        bytes.len()
    );
}

#[test]
fn check_chrom_map_is_not_for_alignments() {
    let d = tempdir();
    let map = plain_file(d.path(), "m.tsv", "chr1\tone\n");
    let sam = sam_fixture("coord-unsorted.sam");
    let mapf = s(&map);
    let files = [
        sam_fixture("coord-unsorted.bam"),
        sam,
        s(&inputs().join("aln/flags.cram")),
    ];
    for f in &files {
        for quick in [false, true] {
            let mut args = vec!["check", "--chrom-map", &mapf];
            if quick {
                args.push("--quick");
            }
            args.push(f);
            let (code, out, err) = run(&args);
            assert_eq!(code, 2, "{f} quick={quick}: {out}{err}");
            assert!(
                err.contains("--chrom-map does not apply to BAM, CRAM or SAM"),
                "{f}: {err}"
            );
        }
    }
    // Interval files keep it.
    let bed = plain_file(d.path(), "a.bed", "chr1\t1\t5\n");
    gloci()
        .args(["check", "--chrom-map", &s(&map), &s(&bed)])
        .assert()
        .success();
}
