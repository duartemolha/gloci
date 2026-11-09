//! Large-file tests (nightly): `cargo test --release --test large -- --ignored`.
//!
//! Data is generated deterministically into a temp dir and removed afterwards.
//! `GLOCI_LARGE_BYTES` overrides the size of the generated BED payload (default
//! 2 GiB) for smaller local runs. Peak RSS of each child is read from
//! `wait4` (the same `ru_maxrss` that `getrusage(RUSAGE_CHILDREN)` reports, but
//! per child); it is KiB on Linux and bytes on macOS.
#![allow(clippy::zombie_processes)] // children are reaped by wait4 in wait_rss
mod common;
use common::*;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const MIB: u64 = 1 << 20;

fn payload_bytes() -> u64 {
    std::env::var("GLOCI_LARGE_BYTES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2 << 30)
}

/// splitmix64: tiny, fast, deterministic.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// Write BED6 lines (unsorted, 24 chromosomes) until `bytes` have been written.
fn gen_bed(w: &mut impl Write, bytes: u64, seed: u64) {
    let mut rng = Rng(seed);
    let mut line = Vec::with_capacity(96);
    let mut buf = Vec::with_capacity(1 << 20);
    let (mut written, mut i) = (0u64, 0u64);
    while written < bytes {
        let r = rng.next();
        let chrom = 1 + r % 24;
        let start = (r >> 8) % 249_000_000;
        let end = start + 50 + (r >> 40) % 10_000;
        line.clear();
        writeln!(
            line,
            "chr{chrom}\t{start}\t{end}\tr{i}\t{}\t{}",
            (r >> 20) % 1000,
            if r & 1 == 0 { '+' } else { '-' }
        )
        .unwrap();
        buf.extend_from_slice(&line);
        written += line.len() as u64;
        i += 1;
        if buf.len() >= 1 << 20 {
            w.write_all(&buf).unwrap();
            buf.clear();
        }
    }
    w.write_all(&buf).unwrap();
}

/// Write sorted BED6 lines in karyotype order until `bytes` have been written.
fn gen_sorted_bed(w: &mut impl Write, bytes: u64, seed: u64) {
    let mut rng = Rng(seed);
    let mut line = Vec::with_capacity(96);
    let mut buf = Vec::with_capacity(1 << 20);
    let (mut written, mut i) = (0u64, 0u64);
    let bytes_per_chrom = (bytes / 22).max(1024);
    let mut chrom = 1;
    let mut chrom_written = 0u64;
    let mut pos = 100u64;
    while written < bytes {
        if chrom_written >= bytes_per_chrom && chrom < 22 {
            chrom += 1;
            chrom_written = 0;
            pos = 100;
        }
        let step = 10 + (rng.next() % 200);
        pos += step;
        let end = pos + 50 + (rng.next() % 500);
        line.clear();
        writeln!(
            line,
            "chr{chrom}\t{pos}\t{end}\tr{i}\t{}\t{}",
            (rng.next() >> 20) % 1000,
            if rng.next() & 1 == 0 { '+' } else { '-' }
        )
        .unwrap();
        buf.extend_from_slice(&line);
        written += line.len() as u64;
        chrom_written += line.len() as u64;
        i += 1;
        if buf.len() >= 1 << 20 {
            w.write_all(&buf).unwrap();
            buf.clear();
        }
    }
    w.write_all(&buf).unwrap();
}

/// Write sorted VCF lines in karyotype order until `bytes` have been written.
fn gen_sorted_vcf(w: &mut impl Write, bytes: u64, seed: u64) {
    let mut rng = Rng(seed);
    let mut line = Vec::with_capacity(128);
    let mut buf = Vec::with_capacity(1 << 20);
    let header = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    buf.extend_from_slice(header.as_bytes());
    let mut written = buf.len() as u64;
    let mut i = 0u64;
    let bytes_per_chrom = (bytes / 22).max(1024);
    let mut chrom = 1;
    let mut chrom_written = 0u64;
    let mut pos = 100u64;
    while written < bytes {
        if chrom_written >= bytes_per_chrom && chrom < 22 {
            chrom += 1;
            chrom_written = 0;
            pos = 100;
        }
        pos += 10 + (rng.next() % 200);
        line.clear();
        let qual = 10 + (rng.next() % 80);
        writeln!(line, "chr{chrom}\t{pos}\trs{i}\tA\tC\t{qual}\tPASS\tDP=30").unwrap();
        buf.extend_from_slice(&line);
        written += line.len() as u64;
        chrom_written += line.len() as u64;
        i += 1;
        if buf.len() >= 1 << 20 {
            w.write_all(&buf).unwrap();
            buf.clear();
        }
    }
    w.write_all(&buf).unwrap();
}

/// Peak RSS of a finished child in bytes, from its rusage.
fn max_rss_bytes(ru: &libc::rusage) -> u64 {
    let v = ru.ru_maxrss as u64;
    if cfg!(target_os = "macos") {
        v
    } else {
        v * 1024
    }
}

/// Reap `child` (Child's drop does not wait) and return its exit success and peak RSS in bytes.
fn wait_rss(child: &std::process::Child) -> (bool, u64) {
    let pid = child.id() as libc::pid_t;
    let mut status = 0;
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    let r = unsafe { libc::wait4(pid, &mut status, 0, &mut ru) };
    assert_eq!(r, pid, "wait4 failed");
    (
        libc::WIFEXITED(status) && libc::WEXITSTATUS(status) == 0,
        max_rss_bytes(&ru),
    )
}

fn gloci_bin() -> std::path::PathBuf {
    assert_cmd::cargo::cargo_bin("gloci")
}

fn gen_plain(dir: &Path, bytes: u64) -> std::path::PathBuf {
    let p = dir.join("big.bed");
    let mut f = std::io::BufWriter::new(std::fs::File::create(&p).unwrap());
    gen_bed(&mut f, bytes, 42);
    f.flush().unwrap();
    p
}

struct SharedDir {
    dir: tempfile::TempDir,
}
static SHARED_DIR: std::sync::OnceLock<SharedDir> = std::sync::OnceLock::new();
fn shared_dir() -> &'static Path {
    SHARED_DIR
        .get_or_init(|| SharedDir { dir: tempdir() })
        .dir
        .path()
}

static BED_A: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
fn big_bed_a() -> &'static Path {
    BED_A.get_or_init(|| {
        let p = shared_dir().join("big_a.bed");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&p).unwrap());
        gen_sorted_bed(&mut f, payload_bytes(), 42);
        f.flush().unwrap();
        p
    })
}

static BED_B: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
fn big_bed_b() -> &'static Path {
    BED_B.get_or_init(|| {
        let p = shared_dir().join("big_b.bed");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&p).unwrap());
        gen_sorted_bed(&mut f, payload_bytes() / 4, 99);
        f.flush().unwrap();
        p
    })
}

static VCF_FILE: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
fn big_vcf() -> &'static Path {
    VCF_FILE.get_or_init(|| {
        let p = shared_dir().join("big.vcf");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&p).unwrap());
        gen_sorted_vcf(&mut f, payload_bytes(), 123);
        f.flush().unwrap();
        p
    })
}

static BED_GZ: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
fn big_bed_gz() -> &'static Path {
    BED_GZ.get_or_init(|| {
        let p = shared_dir().join("big.bed.gz");
        let f = std::fs::File::create(&p).unwrap();
        let mut w = noodles_bgzf::io::Writer::new(std::io::BufWriter::new(f));
        gen_sorted_bed(&mut w, payload_bytes(), 43);
        w.finish().unwrap();
        p
    })
}

#[test]
fn generator_and_sort_check_small() {
    let d = tempdir();
    let p = gen_plain(d.path(), MIB);
    gloci()
        .args(["sort", "--chrom-order", "lex", "--check"])
        .arg(&p)
        .assert()
        .code(1);
    let sorted = d.path().join("s.bed");
    gloci()
        .args(["sort", "--chrom-order", "lex", "-o"])
        .arg(&sorted)
        .arg(&p)
        .assert()
        .success();
    gloci()
        .args(["sort", "--chrom-order", "lex", "--check"])
        .arg(&sorted)
        .assert()
        .success();
}

#[test]
#[ignore = "nightly: generates a 2 GiB file"]
fn sort_2gib_bed_obeys_max_mem() {
    let d = tempdir();
    let input = gen_plain(d.path(), payload_bytes());
    let tmp = d.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    let t0 = std::time::Instant::now();
    // Sort to stdout, piped into a second gloci that checks the order, so the
    // sorted output is never stored on disk.
    let mut sort = Command::new(gloci_bin())
        .args([
            "sort",
            "--chrom-order",
            "lex",
            "--max-mem",
            "512M",
            "--tmp-dir",
        ])
        .arg(&tmp)
        .arg(&input)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let check = Command::new(gloci_bin())
        .args(["sort", "--chrom-order", "lex", "--check", "-"])
        .stdin(sort.stdout.take().unwrap())
        .spawn()
        .unwrap();
    let (sort_ok, sort_rss) = wait_rss(&sort);
    let (check_ok, _) = wait_rss(&check);
    eprintln!(
        "sort 2GiB: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        sort_rss / MIB
    );
    assert!(sort_ok, "gloci sort failed");
    assert!(check_ok, "sorted output failed gloci sort --check");
    assert!(
        sort_rss <= (512 + 150) * MIB,
        "peak RSS {} MiB exceeds 662 MiB",
        sort_rss / MIB
    );
    assert_eq!(
        std::fs::read_dir(&tmp).unwrap().count(),
        0,
        "spill files left"
    );
}

#[test]
#[ignore = "nightly: generates a 2 GiB (uncompressed) .bed.gz"]
fn view_2gib_bed_gz_streams() {
    let p = big_bed_gz();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .arg("view")
        .arg(p)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "view 2GiB bed.gz: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci view failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn intersect_2gib_bed_streams() {
    let a = big_bed_a();
    let b = big_bed_b();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["intersect", "-a"])
        .arg(a)
        .arg("-b")
        .arg(b)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "intersect large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci intersect failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn merge_2gib_bed_streams() {
    let a = big_bed_a();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .arg("merge")
        .arg(a)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "merge large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci merge failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn subtract_2gib_bed_streams() {
    let a = big_bed_a();
    let b = big_bed_b();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["subtract", "-a"])
        .arg(a)
        .arg("-b")
        .arg(b)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "subtract large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci subtract failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn closest_2gib_bed_streams() {
    let a = big_bed_a();
    let b = big_bed_b();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["closest", "-a"])
        .arg(a)
        .arg("-b")
        .arg(b)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "closest large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci closest failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn map_2gib_bed_streams() {
    let a = big_bed_a();
    let b = big_bed_b();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["map", "-a"])
        .arg(a)
        .arg("-b")
        .arg(b)
        .args(["-c", "5:sum"])
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "map large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci map failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn coverage_2gib_bed_streams() {
    let a = big_bed_a();
    let b = big_bed_b();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["coverage", "-a"])
        .arg(a)
        .arg("-b")
        .arg(b)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "coverage large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci coverage failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn complement_2gib_bed_streams() {
    let a = big_bed_a();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .arg("complement")
        .arg(a)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "complement large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci complement failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn vcf_filter_2gib_vcf_streams() {
    let vcf = big_vcf();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["vcf", "filter", "--min-qual", "30"])
        .arg(vcf)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "vcf filter large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci vcf filter failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn check_2gib_bed_rss() {
    let a = big_bed_a();
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .arg("check")
        .arg(a)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "check large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci check failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
}

#[test]
#[ignore = "nightly: generates large file"]
fn index_2gib_bed_gz_rss() {
    let d = tempdir();
    let src = big_bed_gz();
    let target = d.path().join("big.bed.gz");
    if std::fs::hard_link(src, &target).is_err() {
        std::fs::copy(src, &target).unwrap();
    }
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .arg("index")
        .arg(&target)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "index large: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci index failed");
    assert!(rss <= 200_000_000, "peak RSS {} bytes exceeds 200 MB", rss);
    assert!(
        d.path().join("big.bed.gz.tbi").exists() || d.path().join("big.bed.gz.csi").exists(),
        "index file missing"
    );
}

/// I7: per-chromosome mode on indexed inputs runs at most 8 chunks at once,
/// sharing one 64 MiB budget, so `-t 8` and `-t 16` stay under the 200 MB
/// ceiling even though `--mode pairs` writes more than it reads.
#[test]
#[ignore = "nightly: generates large file"]
fn intersect_indexed_2gib_rss() {
    let d = tempdir();
    let a = d.path().join("a.bed.gz");
    if std::fs::hard_link(big_bed_gz(), &a).is_err() {
        std::fs::copy(big_bed_gz(), &a).unwrap();
    }
    let b = d.path().join("b.bed.gz");
    for (args, what) in [
        (
            vec![
                "view".into(),
                big_bed_b().display().to_string(),
                "-o".into(),
                b.display().to_string(),
            ],
            "bgzip B",
        ),
        (vec!["index".into(), a.display().to_string()], "index A"),
        (vec!["index".into(), b.display().to_string()], "index B"),
    ] {
        let st = Command::new(gloci_bin()).args(&args).status().unwrap();
        assert!(st.success(), "{what} failed");
    }
    for threads in ["8", "16"] {
        let t0 = std::time::Instant::now();
        let child = Command::new(gloci_bin())
            .args(["intersect", "--mode", "pairs", "-t", threads, "-a"])
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .arg("--tmp-dir")
            .arg(d.path())
            .stdout(Stdio::null())
            .spawn()
            .unwrap();
        let (ok, rss) = wait_rss(&child);
        eprintln!(
            "intersect indexed -t {threads} large: {:.1}s, peak RSS {} MiB",
            t0.elapsed().as_secs_f64(),
            rss / MIB
        );
        assert!(ok, "gloci intersect (indexed, -t {threads}) failed");
        assert!(
            rss <= 200_000_000,
            "-t {threads}: peak RSS {rss} bytes exceeds 200 MB"
        );
    }
}
