//! `gloci idxstats` (D127): the index path and the scan give the same bytes; the rules beyond
//! the samtools goldens (`idxstats-*-vs-samtools`).
mod common;
use common::*;
use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions};
use gloci_ops::{IdxstatsPath, idxstats};
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn copy(dir: &Path, from: &Path, to: &str) -> PathBuf {
    let p = dir.join(to);
    std::fs::copy(from, &p).unwrap();
    p
}

/// What samtools 1.21 prints for aln/flags.bam: chr1 has a placed unmapped read, and one read
/// has no reference.
const FLAGS: &str = "chr1\t100000\t10\t1\nchr2\t100000\t2\t0\n*\t0\t0\t1\n";

const NO_INDEX: &str = "has no index; counting every record (index it for speed: gloci index ";

fn opts(reference: Option<PathBuf>) -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                reference,
                ..AlignOptions::default()
            }),
            convert_only: false,
        },
        chrom_map: None,
        threads: 2.try_into().unwrap(),
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn lib(path: &Path) -> (IdxstatsPath, String) {
    let mut out = Vec::new();
    let p = idxstats(&s(path), &opts(None), &mut out).unwrap();
    (p, String::from_utf8(out).unwrap())
}

fn run_ok(args: &[&str]) -> (String, String) {
    let o = gloci().args(args).assert().success().get_output().clone();
    (
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
    )
}

#[test]
fn idxstats_index_matches_samtools() {
    let d = tempdir();
    let bam = copy(d.path(), &inputs().join("aln/flags.bam"), "flags.bam");
    copy(
        d.path(),
        &inputs().join("aln/flags.bam.bai"),
        "flags.bam.bai",
    );
    let (path, text) = lib(&bam);
    // The golden bytes are equal through the index and by a scan, so this assertion is what
    // tells the two paths apart.
    assert_eq!(path, IdxstatsPath::Index);
    assert_eq!(text, FLAGS);
    let (out, err) = run_ok(&["idxstats", &s(&bam)]);
    assert_eq!(out, FLAGS);
    assert!(!err.contains("counting every record"), "{err}");
}

#[test]
fn idxstats_scan_equals_index() {
    let d = tempdir();
    let bam = copy(d.path(), &inputs().join("aln/flags.bam"), "flags.bam");
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Scan);
    assert_eq!(text, FLAGS);
    let (out, err) = run_ok(&["idxstats", &s(&bam)]);
    assert_eq!(out, FLAGS);
    assert!(
        err.contains(&format!("{} {NO_INDEX}{})", s(&bam), s(&bam))),
        "{err}"
    );
    assert_eq!(err.matches("counting every record").count(), 1, "{err}");
}

/// The BAI of aln/flags.bam with every metadata pseudo-bin (bin 37450) removed:
/// what an index written by an old or minimal tool holds.
fn strip_pseudo_bins(bai: &[u8]) -> Vec<u8> {
    const META_BIN: u32 = 37450;
    let u32_at = |at: usize| u32::from_le_bytes(bai[at..at + 4].try_into().unwrap());
    let n_ref = u32_at(4) as usize;
    let mut out = bai[..8].to_vec();
    let mut at = 8;
    for _ in 0..n_ref {
        let n_bin = u32_at(at) as usize;
        at += 4;
        let mut kept: Vec<&[u8]> = Vec::new();
        for _ in 0..n_bin {
            let n_chunk = u32_at(at + 4) as usize;
            let len = 8 + 16 * n_chunk;
            if u32_at(at) != META_BIN {
                kept.push(&bai[at..at + len]);
            }
            at += len;
        }
        out.extend_from_slice(&(kept.len() as u32).to_le_bytes());
        for b in kept {
            out.extend_from_slice(b);
        }
        let n_intv = u32_at(at) as usize;
        let len = 4 + 8 * n_intv;
        out.extend_from_slice(&bai[at..at + len]);
        at += len;
    }
    // `n_no_coor` stays: only the pseudo-bins are missing.
    out.extend_from_slice(&bai[at..]);
    out
}

/// Writes `xtask/diff/inputs/samtools/nometa/flags.bam.bai`, the fixture of the
/// `idxstats-index-without-counts` case. Run: `cargo test -p gloci --test idxstats
/// make_nometa_index -- --ignored`.
#[test]
#[ignore]
fn make_nometa_index() {
    let bai = std::fs::read(inputs().join("aln/flags.bam.bai")).unwrap();
    std::fs::write(
        inputs().join("samtools/nometa/flags.bam.bai"),
        strip_pseudo_bins(&bai),
    )
    .unwrap();
}

#[test]
fn nometa_fixture_is_the_stripped_index() {
    let bai = std::fs::read(inputs().join("aln/flags.bam.bai")).unwrap();
    assert_eq!(
        std::fs::read(inputs().join("samtools/nometa/flags.bam.bai")).unwrap(),
        strip_pseudo_bins(&bai),
        "rewrite it with: cargo test -p gloci --test idxstats make_nometa_index -- --ignored"
    );
    assert_eq!(
        std::fs::read(inputs().join("samtools/nometa/flags.bam")).unwrap(),
        std::fs::read(inputs().join("aln/flags.bam")).unwrap()
    );
}

#[test]
fn idxstats_index_without_counts() {
    let d = tempdir();
    let bam = copy_with_fresh_index(
        d.path(),
        &inputs().join("samtools/nometa/flags.bam"),
        &inputs().join("samtools/nometa/flags.bam.bai"),
    );
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Scan);
    assert_eq!(text, FLAGS);
    let (out, err) = run_ok(&["idxstats", &s(&bam)]);
    assert_eq!(out, FLAGS);
    assert_eq!(err.matches("warning").count(), 1, "{err}");
    assert!(err.contains("has no read counts for chr1"), "{err}");
}

#[test]
fn idxstats_unsorted_scan() {
    let sam = inputs().join("samtools/unsorted-chroms.sam");
    let bam = inputs().join("samtools/unsorted-chroms.bam");
    let want = "chr1\t1000\t3\t0\nchr2\t2000\t1\t1\nchr3\t3000\t0\t0\n*\t0\t0\t1\n";
    for f in [&sam, &bam] {
        let (out, _) = run_ok(&["idxstats", &s(f)]);
        assert_eq!(out, want, "{}", f.display());
    }
}

#[test]
fn idxstats_cram() {
    let cram = s(&inputs().join("aln/flags.cram"));
    let reference = s(&inputs().join("fa/ref.fa"));
    let (out, err) = run_ok(&["idxstats", &cram, "--reference", &reference]);
    assert_eq!(out, "chr1\t1000\t10\t1\nchr2\t1234\t2\t0\n*\t0\t0\t1\n");
    assert!(err.contains("counting every record"), "{err}");
    gloci()
        .args(["idxstats", &cram])
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains("needs the reference sequence of chr1").and(contains("--reference")));
    // The CRAM's SAM twin counts the same.
    let (sam, _) = run_ok(&["idxstats", &s(&inputs().join("aln/flags.sam"))]);
    assert_eq!(sam, FLAGS);
}

#[test]
fn idxstats_scan_is_identical_at_any_thread_count() {
    let bam = s(&inputs().join("aln/many.cram"));
    let reference = s(&inputs().join("fa/ref.fa"));
    let want = run_ok(&["idxstats", &bam, "--reference", &reference, "-t", "1"]).0;
    assert!(want.lines().count() > 2);
    for t in ["2", "8"] {
        let got = run_ok(&["idxstats", &bam, "--reference", &reference, "-t", t]).0;
        assert_eq!(got, want, "-t {t}");
    }
    let sam = s(&inputs().join("aln/cat-big.sam"));
    let want = run_ok(&["idxstats", &sam, "-t", "1"]).0;
    for t in ["2", "8"] {
        assert_eq!(run_ok(&["idxstats", &sam, "-t", t]).0, want, "-t {t}");
    }
    let bam = s(&inputs().join("aln/cat-big.bam"));
    let want = run_ok(&["idxstats", &bam, "-t", "1"]).0;
    for t in ["2", "8"] {
        assert_eq!(run_ok(&["idxstats", &bam, "-t", t]).0, want, "-t {t}");
    }
}

#[test]
fn idxstats_stdin_is_scanned() {
    let sam = std::fs::read(inputs().join("aln/flags.sam")).unwrap();
    gloci()
        .args(["idxstats", "-"])
        .write_stdin(sam)
        .assert()
        .success()
        .stdout(FLAGS);
}

#[test]
fn idxstats_missing_file_exits_1() {
    gloci()
        .args(["idxstats", "/nonexistent/x.bam"])
        .assert()
        .code(1)
        .stderr(contains("cannot open /nonexistent/x.bam"));
}

/// 200,000 references, a read on every 1,000th: the index path holds the header and one
/// counter pair per reference, and does nothing per reference that is quadratic.
#[test]
fn idxstats_huge_header() {
    let d = tempdir();
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for i in 0..200_000 {
        sam.push_str(&format!("@SQ\tSN:c{i}\tLN:1000\n"));
    }
    for i in (0..200_000).step_by(1000) {
        sam.push_str(&format!(
            "r{i}\t0\tc{i}\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n"
        ));
    }
    let sam_path = plain_file(d.path(), "big.sam", &sam);
    let bam = d.path().join("big.bam");
    gloci()
        .args(["view", "-O", "bam", &s(&sam_path), "-o", &s(&bam)])
        .assert()
        .success();
    gloci().args(["index", &s(&bam)]).assert().success();
    let t = Instant::now();
    let (path, text) = lib(&bam);
    let took = t.elapsed();
    assert_eq!(path, IdxstatsPath::Index);
    // A loose bound for unoptimised builds on a loaded box; quadratic work would take minutes.
    assert!(took.as_secs() < 15, "took {took:?}");
    assert_eq!(text.lines().count(), 200_001);
    assert_eq!(text.lines().next(), Some("c0\t1000\t1\t0"));
    assert_eq!(text.lines().nth(1), Some("c1\t1000\t0\t0"));
    assert_eq!(text.lines().last(), Some("*\t0\t0\t0"));
    // The scan agrees.
    let plain = copy(d.path(), &bam, "plain.bam");
    let (path, scanned) = lib(&plain);
    assert_eq!(path, IdxstatsPath::Scan);
    assert_eq!(scanned, text);
}

#[test]
fn idxstats_csi_matches_bai() {
    let d = tempdir();
    let bam = copy(d.path(), &inputs().join("aln/flags.bam"), "flags.bam");
    gloci()
        .args(["index", "--csi", &s(&bam)])
        .assert()
        .success();
    assert!(d.path().join("flags.bam.csi").exists());
    assert!(!d.path().join("flags.bam.bai").exists());
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Index);
    assert_eq!(text, FLAGS);
}

/// The CSI `samtools index -c` wrote (samtools/csi/flags.bam.csi) is read the same way.
#[test]
fn idxstats_samtools_csi() {
    let d = tempdir();
    let bam = copy_with_fresh_index(
        d.path(),
        &inputs().join("samtools/csi/flags.bam"),
        &inputs().join("samtools/csi/flags.bam.csi"),
    );
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Index);
    assert_eq!(text, FLAGS);
}

/// A reference over 2^29 gets a CSI whose metadata pseudo-bin sits deeper than a BAI's.
#[test]
fn idxstats_long_reference_csi() {
    let d = tempdir();
    let mut sam = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:big\tLN:600000000\n@SQ\tSN:small\tLN:1000\n",
    );
    sam.push_str("r1\t0\tbig\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    sam.push_str("r2\t0\tbig\t590000000\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    sam.push_str("r3\t4\tbig\t590000000\t0\t*\t*\t0\t0\tACGTA\tIIIII\n");
    sam.push_str("r4\t0\tsmall\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    sam.push_str("r5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n");
    let sam_path = plain_file(d.path(), "long.sam", &sam);
    let bam = d.path().join("long.bam");
    gloci()
        .args(["view", "-O", "bam", &s(&sam_path), "-o", &s(&bam)])
        .assert()
        .success();
    gloci().args(["index", &s(&bam)]).assert().success();
    assert!(d.path().join("long.bam.csi").exists());
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Index);
    assert_eq!(
        text,
        "big\t600000000\t2\t1\nsmall\t1000\t1\t0\n*\t0\t0\t1\n"
    );
}

#[test]
fn idxstats_bam_on_stdin_is_not_told_to_index() {
    let bam = std::fs::read(inputs().join("aln/flags.bam")).unwrap();
    let o = gloci()
        .args(["idxstats", "-"])
        .write_stdin(bam)
        .assert()
        .success()
        .stdout(FLAGS)
        .get_output()
        .clone();
    let err = String::from_utf8(o.stderr).unwrap();
    assert_eq!(
        err,
        "warning: stdin is not a regular file, so no index can be read; counting every record\n"
    );
    // A FIFO (process substitution) is not a regular file either: no index is suggested.
    let p = s(&inputs().join("aln/flags.bam"));
    let o = bash_gloci(&format!("idxstats <(cat '{p}')"), "");
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(String::from_utf8_lossy(&o.stdout), FLAGS);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("is not a regular file, so no index can be read; counting every record"),
        "{err}"
    );
    assert!(!err.contains("gloci index"), "{err}");
}

/// An index without the trailing `n_no_coor` cannot give the unplaced count: no silent zero.
#[test]
fn idxstats_index_without_unplaced_count() {
    let d = tempdir();
    let bam = copy(d.path(), &inputs().join("aln/flags.bam"), "flags.bam");
    let bai = std::fs::read(inputs().join("aln/flags.bam.bai")).unwrap();
    std::fs::write(d.path().join("flags.bam.bai"), &bai[..bai.len() - 8]).unwrap();
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Scan);
    assert_eq!(text, FLAGS);
    let (out, err) = run_ok(&["idxstats", &s(&bam)]);
    assert_eq!(out, FLAGS);
    assert_eq!(err.matches("warning").count(), 1, "{err}");
    assert!(err.contains("no count of unplaced reads"), "{err}");
}

#[test]
fn idxstats_stale_index_warns_once_and_is_used() {
    let d = tempdir();
    let bam = copy(d.path(), &inputs().join("aln/flags.bam"), "flags.bam");
    let bai = copy(
        d.path(),
        &inputs().join("aln/flags.bam.bai"),
        "flags.bam.bai",
    );
    let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
    std::fs::File::options()
        .write(true)
        .open(&bai)
        .unwrap()
        .set_modified(old)
        .unwrap();
    let (path, text) = lib(&bam);
    assert_eq!(path, IdxstatsPath::Index);
    assert_eq!(text, FLAGS);
    let (out, err) = run_ok(&["idxstats", &s(&bam)]);
    assert_eq!(out, FLAGS);
    assert_eq!(err.matches("older than").count(), 1, "{err}");
    assert_eq!(err.matches("warning").count(), 1, "{err}");
}

#[test]
fn idxstats_writes_to_o() {
    let d = tempdir();
    let bam = copy(d.path(), &inputs().join("aln/flags.bam"), "flags.bam");
    for name in ["stats.tsv", "stats.txt", "stats"] {
        let out = d.path().join(name);
        gloci()
            .args(["idxstats", &s(&bam), "-o", &s(&out)])
            .assert()
            .success()
            .stdout(predicate::str::is_empty());
        assert_eq!(std::fs::read_to_string(&out).unwrap(), FLAGS, "{name}");
    }
    // .gz compresses the table (BGZF).
    let gz = d.path().join("stats.tsv.gz");
    gloci()
        .args(["idxstats", &s(&bam), "--output", &s(&gz)])
        .assert()
        .success();
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut noodles_bgzf::io::Reader::new(std::fs::File::open(&gz).unwrap()),
        &mut text,
    )
    .unwrap();
    assert_eq!(text, FLAGS);
    // A name of another format is refused before anything is read.
    let bad = d.path().join("stats.bam");
    gloci()
        .args(["idxstats", &s(&bam), "-o", &s(&bad)])
        .assert()
        .code(2);
    assert!(!bad.exists());
    // A failing scan leaves no output file.
    let missing = d.path().join("missing.bam");
    let out = d.path().join("never.tsv");
    gloci()
        .args(["idxstats", &s(&missing), "-o", &s(&out)])
        .assert()
        .code(1);
    assert!(!out.exists());
}
