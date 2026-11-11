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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
    assert!(
        rss <= 200_000_000,
        "peak RSS {} bytes exceeds the ceiling",
        rss
    );
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
            "-t {threads}: peak RSS {rss} bytes exceeds the ceiling"
        );
    }
}

/// Run `args`, discard stdout, and return the peak RSS in bytes; asserts success.
fn rss_of(what: &str, args: &[&std::ffi::OsStr]) -> u64 {
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(args)
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "{what}: {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    assert!(ok, "gloci {what} failed");
    rss
}

/// A genome file for the generated BED: chr1..chr22 at 2 Gbp each, so no
/// generated position falls off the end at any payload size.
fn big_genome() -> PathBuf {
    let p = shared_dir().join("big.genome");
    if !p.exists() {
        let mut t = String::new();
        for c in 1..=22 {
            t.push_str(&format!("chr{c}\t2000000000\n"));
        }
        std::fs::write(&p, t).unwrap();
    }
    p
}

/// Peak RSS ceiling of the streaming P2 commands: about three times the
/// largest value measured (22 MiB). It is meaningful only at the default
/// 2 GiB payload: with `GLOCI_LARGE_BYTES` set small, a command that buffered
/// its whole input would still fit under it.
const RSS_CEILING: u64 = 64 << 20;

#[test]
#[ignore = "nightly: generates large file"]
fn partition_2gib_bed_streams() {
    let (a, b) = (big_bed_a(), big_bed_b());
    let rss = rss_of(
        "partition large",
        &[
            "partition".as_ref(),
            a.as_os_str(),
            b.as_os_str(),
            "--names".as_ref(),
            "A,B".as_ref(),
        ],
    );
    assert!(rss <= RSS_CEILING, "peak RSS {rss} bytes exceeds 64 MiB");
}

#[test]
#[ignore = "nightly: generates large file"]
fn genomecov_2gib_bed_streams() {
    let a = big_bed_a();
    let g = big_genome();
    for mode in ["--bedgraph", "--bga", "--hist"] {
        let rss = rss_of(
            &format!("genomecov {mode} large"),
            &[
                "genomecov".as_ref(),
                a.as_os_str(),
                "-g".as_ref(),
                g.as_os_str(),
                mode.as_ref(),
            ],
        );
        assert!(
            rss <= RSS_CEILING,
            "genomecov {mode}: peak RSS {rss} bytes exceeds 64 MiB"
        );
    }
}

#[test]
#[ignore = "nightly: generates large file"]
fn cluster_2gib_bed_streams() {
    let a = big_bed_a();
    let rss = rss_of("cluster large", &["cluster".as_ref(), a.as_os_str()]);
    assert!(rss <= RSS_CEILING, "peak RSS {rss} bytes exceeds 64 MiB");
    let rss = rss_of("spacing large", &["spacing".as_ref(), a.as_os_str()]);
    assert!(
        rss <= RSS_CEILING,
        "spacing: peak RSS {rss} bytes exceeds 64 MiB"
    );
}

const STATS_SAMPLES: usize = 50;

/// Sorted VCF with `STATS_SAMPLES` genotype columns until `bytes` have been
/// written; returns the number of records.
fn gen_sample_vcf(w: &mut impl Write, bytes: u64, seed: u64) -> u64 {
    let mut rng = Rng(seed);
    let mut buf = Vec::with_capacity(1 << 20);
    let mut head = String::from(
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for i in 1..=STATS_SAMPLES {
        head.push_str(&format!("\ts{i}"));
    }
    head.push('\n');
    buf.extend_from_slice(head.as_bytes());
    let mut written = buf.len() as u64;
    let bytes_per_chrom = (bytes / 22).max(1024);
    let (mut chrom, mut chrom_written, mut pos, mut n) = (1, 0u64, 100u64, 0u64);
    let mut line = String::new();
    while written < bytes {
        if chrom_written >= bytes_per_chrom && chrom < 22 {
            chrom += 1;
            chrom_written = 0;
            pos = 100;
        }
        pos += 10 + rng.next() % 200;
        line.clear();
        let (r, a) = [("A", "C"), ("C", "T"), ("G", "A"), ("AT", "A")][(rng.next() % 4) as usize];
        line.push_str(&format!(
            "chr{chrom}\t{pos}\trs{n}\t{r}\t{a}\t{}\tPASS\tDP=30\tGT:DP",
            10 + rng.next() % 80
        ));
        for _ in 0..STATS_SAMPLES {
            let x = rng.next();
            let gt = ["0/0", "0/1", "1/1", "./.", "0|1"][(x % 5) as usize];
            line.push_str(&format!("\t{gt}:{}", (x >> 8) % 60));
        }
        line.push('\n');
        buf.extend_from_slice(line.as_bytes());
        written += line.len() as u64;
        chrom_written += line.len() as u64;
        n += 1;
        if buf.len() >= 1 << 20 {
            w.write_all(&buf).unwrap();
            buf.clear();
        }
    }
    w.write_all(&buf).unwrap();
    n
}

fn lines_of(path: &Path) -> u64 {
    use std::io::BufRead;
    std::io::BufReader::new(std::fs::File::open(path).unwrap())
        .split(b'\n')
        .count() as u64
}

#[test]
#[ignore = "nightly: generates large file"]
fn vcf_stats_2gib_vcf_streams() {
    let vcf = shared_dir().join("stats_in.vcf");
    let sites = {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&vcf).unwrap());
        let n = gen_sample_vcf(&mut f, payload_bytes(), 321);
        f.flush().unwrap();
        n
    };
    let prefix = shared_dir().join("stats_out");
    let t0 = std::time::Instant::now();
    let child = Command::new(gloci_bin())
        .args(["vcf", "stats"])
        .arg(&vcf)
        .arg("--out-prefix")
        .arg(&prefix)
        .args([
            "--site-quality",
            "--freq",
            "--missing-indv",
            "--het",
            "--hardy",
            "--filter-summary",
            "--indel-hist",
            "--tstv-summary",
            "--summary",
            "--snp-density",
            "100000",
        ])
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    let (ok, rss) = wait_rss(&child);
    eprintln!(
        "vcf stats large ({sites} sites, {STATS_SAMPLES} samples): {:.1}s, peak RSS {} MiB",
        t0.elapsed().as_secs_f64(),
        rss / MIB
    );
    // one line per site (plus the header) in the per-site tables, and one
    // per sample in the per-sample ones
    let ext = |e: &str| shared_dir().join(format!("stats_out.{e}"));
    let counts = (
        lines_of(&ext("lqual")),
        lines_of(&ext("frq")),
        lines_of(&ext("hwe")),
        lines_of(&ext("het")),
        lines_of(&ext("imiss")),
    );
    let _ = std::fs::remove_file(&vcf);
    for e in std::fs::read_dir(shared_dir()).unwrap() {
        let p = e.unwrap().path();
        if p.file_name()
            .is_some_and(|n| n.to_string_lossy().starts_with("stats_out."))
        {
            let _ = std::fs::remove_file(p);
        }
    }
    assert!(ok, "gloci vcf stats failed");
    assert_eq!(counts.0, sites + 1, ".lqual lines");
    assert_eq!(counts.1, sites + 1, ".frq lines");
    assert_eq!(counts.2, sites + 1, ".hwe lines");
    assert_eq!(counts.3, STATS_SAMPLES as u64 + 1, ".het lines");
    assert_eq!(counts.4, STATS_SAMPLES as u64 + 1, ".imiss lines");
    assert!(rss <= RSS_CEILING, "peak RSS {rss} bytes exceeds 64 MiB");
}

// ---------------------------------------------------------------------------
// P3a: BAM, CRAM, FASTA and alignment sort
// ---------------------------------------------------------------------------

/// SAM bytes per BAM byte for the generated reads (random bases, constant qualities): a BAM of
/// `payload_bytes()` needs about this many bytes of SAM text.
const SAM_PER_BAM: u64 = 7;
const READ_LEN: usize = 100;

fn chrom_names(n: u64) -> Vec<String> {
    (1..=n).map(|c| format!("chr{c}")).collect()
}

/// SAM header over `chroms` chromosomes of `len` bases.
fn sam_header(chroms: u64, len: u64, sorted: bool) -> String {
    let mut h = format!(
        "@HD\tVN:1.6\tSO:{}\n",
        if sorted { "coordinate" } else { "unsorted" }
    );
    for c in chrom_names(chroms) {
        h.push_str(&format!("@SQ\tSN:{c}\tLN:{len}\n"));
    }
    h
}

/// Write SAM text of random-base 100 bp reads until `bytes` have been written: over `chroms`
/// chromosomes of `len` bases, in coordinate order (samtools tie order) or at random places.
fn gen_sam(w: &mut impl Write, bytes: u64, chroms: u64, len: u64, sorted: bool, seed: u64) {
    let mut rng = Rng(seed);
    let header = sam_header(chroms, len, sorted);
    w.write_all(header.as_bytes()).unwrap();
    let mut written = header.len() as u64;
    let mut buf = Vec::with_capacity(1 << 20);
    let mut line = Vec::with_capacity(512);
    let qual = "I".repeat(READ_LEN);
    let per_chrom = (bytes / chroms).max(1024);
    let (mut chrom, mut chrom_written, mut pos) = (1u64, 0u64, 1u64);
    let mut i = 0u64;
    while written < bytes {
        let r = rng.next();
        let (c, p) = if sorted {
            if chrom_written >= per_chrom && chrom < chroms {
                chrom += 1;
                chrom_written = 0;
                pos = 1;
            }
            pos += 1 + rng.next() % 60;
            (chrom, pos.min(len - READ_LEN as u64))
        } else {
            (1 + r % chroms, 1 + (r >> 8) % (len - READ_LEN as u64))
        };
        // a strand flag; sorted input has a forward read before a reverse one at a position
        let flag = if rng.next() & 1 == 0 { 0 } else { 16 };
        let mut seq = [0u8; READ_LEN];
        let mut bits = rng.next();
        for (k, s) in seq.iter_mut().enumerate() {
            if k % 30 == 29 {
                bits = rng.next();
            }
            *s = b"ACGT"[(bits >> (2 * (k % 30))) as usize & 3];
        }
        line.clear();
        write!(
            line,
            "read{i}\t{flag}\tchr{c}\t{p}\t60\t{READ_LEN}M\t*\t0\t0\t"
        )
        .unwrap();
        line.extend_from_slice(&seq);
        line.push(b'\t');
        line.extend_from_slice(qual.as_bytes());
        line.push(b'\n');
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

/// A BAM of about `bam_bytes` made from generated SAM (removed afterwards); `index` also writes
/// the `.bai`.
fn make_bam(
    dir: &Path,
    name: &str,
    bam_bytes: u64,
    sorted: bool,
    seed: u64,
    index: bool,
) -> PathBuf {
    let sam = dir.join(format!("{name}.sam"));
    let mut f = std::io::BufWriter::new(std::fs::File::create(&sam).unwrap());
    gen_sam(
        &mut f,
        bam_bytes * SAM_PER_BAM,
        22,
        249_000_000,
        sorted,
        seed,
    );
    f.flush().unwrap();
    drop(f);
    let bam = dir.join(format!("{name}.bam"));
    let st = Command::new(gloci_bin())
        .arg("view")
        .arg(&sam)
        .arg("-o")
        .arg(&bam)
        .status()
        .unwrap();
    assert!(st.success(), "SAM to BAM failed");
    let _ = std::fs::remove_file(&sam);
    if index {
        let st = Command::new(gloci_bin())
            .arg("index")
            .arg(&bam)
            .status()
            .unwrap();
        assert!(st.success(), "index BAM failed");
    }
    bam
}

static SORTED_BAM: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
fn big_sorted_bam() -> &'static Path {
    SORTED_BAM.get_or_init(|| make_bam(shared_dir(), "sorted", payload_bytes(), true, 5, true))
}

static SHUFFLED_BAM: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
fn big_shuffled_bam() -> &'static Path {
    SHUFFLED_BAM
        .get_or_init(|| make_bam(shared_dir(), "shuffled", payload_bytes(), false, 6, false))
}

/// Ceiling of the streaming alignment commands.
const ALN_RSS: u64 = 200_000_000;

/// The ceiling at `-t N`. A per-chromosome run on an indexed BAM keeps up to 8 chromosomes in
/// flight, about 20 MB each, so `-t 8` is measured against 200 MB plus that slack (the box that
/// recorded the ceiling has 4 cores; at `-t 8` it peaked between 140 and 208 MB).
fn aln_ceiling(t: &str) -> u64 {
    if t == "8" {
        ALN_RSS + 100_000_000
    } else {
        ALN_RSS
    }
}

#[test]
fn p3a_generators_small() {
    let d = tempdir();
    let bam = make_bam(d.path(), "s", MIB, true, 1, true);
    let n = bam.metadata().unwrap().len();
    assert!(
        n > MIB / 2 && n < 4 * MIB,
        "BAM of {n} bytes for a 1 MiB target"
    );
    gloci()
        .args(["sort", "--check"])
        .arg(&bam)
        .assert()
        .success();
    let shuf = make_bam(d.path(), "u", MIB, false, 2, false);
    gloci()
        .args(["sort", "--check"])
        .arg(&shuf)
        .assert()
        .code(1);
}

#[test]
#[ignore = "nightly: generates a 2 GB BAM"]
fn genomecov_2gb_bam_streams() {
    let bam = big_sorted_bam();
    for t in ["1", "4", "8"] {
        for mode in ["--bedgraph", "--bga"] {
            let rss = rss_of(
                &format!("genomecov {mode} BAM -t {t}"),
                &[
                    "genomecov".as_ref(),
                    bam.as_os_str(),
                    mode.as_ref(),
                    "-t".as_ref(),
                    t.as_ref(),
                ],
            );
            assert!(
                rss < aln_ceiling(t),
                "genomecov {mode} -t {t}: peak RSS {rss} bytes exceeds the ceiling"
            );
        }
    }
}

#[test]
#[ignore = "nightly: generates a 2 GB BAM"]
fn convert_2gb_bam_streams() {
    let bam = big_sorted_bam();
    for t in ["1", "4", "8"] {
        let rss = rss_of(
            &format!("convert BAM -t {t}"),
            &[
                "convert".as_ref(),
                bam.as_os_str(),
                "-t".as_ref(),
                t.as_ref(),
            ],
        );
        assert!(
            rss < aln_ceiling(t),
            "convert -t {t}: peak RSS {rss} bytes exceeds the ceiling"
        );
    }
}

#[test]
#[ignore = "nightly: generates a 2 GB BAM"]
fn coverage_with_2gb_bam_as_b_streams() {
    let bam = big_sorted_bam();
    let a = big_bed_b();
    for t in ["1", "4", "8"] {
        let rss = rss_of(
            &format!("coverage -b BAM -t {t}"),
            &[
                "coverage".as_ref(),
                "-a".as_ref(),
                a.as_os_str(),
                "-b".as_ref(),
                bam.as_os_str(),
                "-t".as_ref(),
                t.as_ref(),
            ],
        );
        assert!(
            rss < aln_ceiling(t),
            "coverage -t {t}: peak RSS {rss} bytes exceeds the ceiling"
        );
    }
}

/// A FASTA of `chroms` chromosomes of `len` random bases (60 per line, some lowercase), with its
/// `.fai`; reused between tests.
fn make_fasta(name: &str, chroms: u64, len: u64, seed: u64, soft: bool) -> PathBuf {
    let p = shared_dir().join(name);
    if p.exists() {
        return p;
    }
    let partial = shared_dir().join(format!("{name}.partial"));
    let mut w = std::io::BufWriter::new(std::fs::File::create(&partial).unwrap());
    let mut rng = Rng(seed);
    for c in 1..=chroms {
        writeln!(w, ">chr{c}").unwrap();
        let mut left = len;
        let mut bits = 0u64;
        let mut lower = false;
        let mut line = [0u8; 61];
        let mut k = 0;
        while left > 0 {
            let n = left.min(60) as usize;
            for s in line.iter_mut().take(n) {
                if k % 28 == 0 {
                    bits = rng.next();
                    if soft && bits & 0xFF < 24 {
                        lower = !lower;
                    }
                }
                let b = b"ACGT"[(bits >> (2 * (k % 28) + 8)) as usize & 3];
                *s = if lower { b.to_ascii_lowercase() } else { b };
                k += 1;
            }
            line[n] = b'\n';
            w.write_all(&line[..=n]).unwrap();
            left -= n as u64;
        }
    }
    w.flush().unwrap();
    drop(w);
    std::fs::rename(&partial, &p).unwrap();
    let st = Command::new(gloci_bin())
        .arg("index")
        .arg(&p)
        .status()
        .unwrap();
    assert!(st.success(), "index FASTA failed");
    p
}

/// 4 chromosomes of 300 Mb: 1.2 GB of FASTA, whatever `GLOCI_LARGE_BYTES` says (the ceiling is
/// only meaningful when a chromosome is much larger than it).
fn big_fasta() -> PathBuf {
    make_fasta("big.fa", 4, 300_000_000, 77, true)
}

#[test]
#[ignore = "nightly: generates a 1.2 GB FASTA"]
fn getfasta_whole_300mb_chromosomes_streams() {
    let fa = big_fasta();
    let bed = shared_dir().join("whole.bed");
    std::fs::write(
        &bed,
        "chr1\t0\t300000000\nchr2\t0\t300000000\nchr4\t0\t300000000\n",
    )
    .unwrap();
    for (what, extra) in [
        ("getfasta", vec![]),
        ("getfasta -s --tab", vec!["-s", "--tab"]),
    ] {
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "getfasta".as_ref(),
            "--fasta".as_ref(),
            fa.as_os_str(),
            bed.as_os_str(),
        ];
        args.extend(extra.iter().map(|s| std::ffi::OsStr::new(*s)));
        let rss = rss_of(what, &args);
        assert!(
            rss < ALN_RSS,
            "{what}: peak RSS {rss} bytes exceeds the ceiling"
        );
    }
    let rss = rss_of(
        "nuc whole chromosomes",
        &[
            "nuc".as_ref(),
            "--fasta".as_ref(),
            fa.as_os_str(),
            bed.as_os_str(),
        ],
    );
    assert!(
        rss < ALN_RSS,
        "nuc: peak RSS {rss} bytes exceeds the ceiling"
    );
}

#[test]
#[ignore = "nightly: generates a 1.2 GB FASTA"]
fn maskfasta_1gb_streams() {
    let fa = big_fasta();
    let masks = shared_dir().join("masks.bed");
    if !masks.exists() {
        let mut rng = Rng(5);
        let mut w = std::io::BufWriter::new(std::fs::File::create(&masks).unwrap());
        // 2 million sorted masks, 500 per chromosome-megabase
        for c in 1..=4u64 {
            let mut pos = 0u64;
            for _ in 0..500_000 {
                pos += 100 + rng.next() % 500;
                writeln!(w, "chr{c}\t{pos}\t{}", pos + 20 + rng.next() % 80).unwrap();
            }
        }
        w.flush().unwrap();
    }
    for (what, extra) in [("maskfasta", vec![]), ("maskfasta --soft", vec!["--soft"])] {
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "maskfasta".as_ref(),
            "--fasta".as_ref(),
            fa.as_os_str(),
            "-b".as_ref(),
            masks.as_os_str(),
        ];
        args.extend(extra.iter().map(|s| std::ffi::OsStr::new(*s)));
        let rss = rss_of(what, &args);
        assert!(
            rss < ALN_RSS,
            "{what}: peak RSS {rss} bytes exceeds the ceiling"
        );
    }
}

/// Bytes of every file under `dir`.
fn dir_bytes(dir: &Path) -> u64 {
    let mut total = 0;
    if let Ok(rd) = std::fs::read_dir(dir) {
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                total += dir_bytes(&p);
            } else if let Ok(m) = e.metadata() {
                total += m.len();
            }
        }
    }
    total
}

/// Sort the shuffled BAM at `--max-mem 512M`, writing nothing to disk but the spills: the peak
/// RSS stays under the budget plus 200 MB, the spill directory never holds more than twice the
/// input, and it is empty afterwards.
fn sort_bam_obeys_max_mem(extra: &[&str]) {
    // 512 MiB for the default 2 GiB payload, scaled down with a smaller `GLOCI_LARGE_BYTES` so
    // the sort still spills (never below the 16 MiB minimum)
    let budget_mib = (512 * payload_bytes() / (2 << 30)).clamp(16, 512);
    let budget = format!("{budget_mib}M");
    let input = big_shuffled_bam();
    let in_bytes = input.metadata().unwrap().len();
    let tmp = shared_dir().join(format!("spill-{}", extra.join("")));
    std::fs::create_dir_all(&tmp).unwrap();
    let t0 = std::time::Instant::now();
    let mut sort = Command::new(gloci_bin());
    sort.arg("sort").args(extra);
    sort.args(["--max-mem", &budget, "--tmp-dir"])
        .arg(&tmp)
        .arg(input);
    let mut sort = sort
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    // the sorted output is checked by a second process, so it never touches the disk
    let check = Command::new(gloci_bin())
        .args(["sort", "--check"])
        .args(extra)
        .arg("-")
        .stdin(sort.stdout.take().unwrap())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let done = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
    let peak = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
    let monitor = {
        let (done, peak, tmp) = (done.clone(), peak.clone(), tmp.clone());
        std::thread::spawn(move || {
            while !done.load(std::sync::atomic::Ordering::Relaxed) {
                peak.fetch_max(dir_bytes(&tmp), std::sync::atomic::Ordering::Relaxed);
                std::thread::sleep(std::time::Duration::from_millis(100));
            }
        })
    };
    let (sort_ok, rss) = wait_rss(&sort);
    let (check_ok, _) = wait_rss(&check);
    done.store(true, std::sync::atomic::Ordering::Relaxed);
    monitor.join().unwrap();
    let spill = peak.load(std::sync::atomic::Ordering::Relaxed);
    eprintln!(
        "sort {extra:?} BAM ({} MiB): {:.1}s, peak RSS {} MiB, peak spill {} MiB",
        in_bytes / MIB,
        t0.elapsed().as_secs_f64(),
        rss / MIB,
        spill / MIB
    );
    assert!(sort_ok, "gloci sort {extra:?} failed");
    assert!(check_ok, "sorted output failed the order check");
    assert!(
        rss <= (budget_mib << 20) + ALN_RSS,
        "peak RSS {} MiB exceeds {budget_mib} MiB + 200 MB",
        rss / MIB
    );
    assert!(
        spill <= 2 * in_bytes,
        "spill {} MiB exceeds twice the input ({} MiB)",
        spill / MIB,
        in_bytes / MIB
    );
    assert_eq!(
        std::fs::read_dir(&tmp).unwrap().count(),
        0,
        "spill files left"
    );
}

#[test]
#[ignore = "nightly: generates a 2 GB BAM"]
fn sort_2gb_bam_obeys_max_mem() {
    sort_bam_obeys_max_mem(&[]);
}

#[test]
#[ignore = "nightly: generates a 2 GB BAM"]
fn sort_n_2gb_bam_obeys_max_mem() {
    sort_bam_obeys_max_mem(&["-n"]);
}

/// A CRAM of `reads` reads sampled from a 3 x 20 Mb reference (written with noodles), and the
/// reference FASTA with its index. One reference sequence is 20 MB.
const CRAM_CHROM: usize = 20_000_000;

/// The CRAM and its reference, made by a child copy of this test binary: writing a CRAM holds the
/// whole reference in memory, and a child's peak RSS (`wait4`) includes its parent's memory at the
/// time it was forked, so doing that here would inflate every measurement made after it.
fn make_cram(name: &str, reads: u64) -> (PathBuf, PathBuf) {
    let fa = make_fasta("cram-ref.fa", 3, CRAM_CHROM as u64, 31, false);
    let cram = shared_dir().join(format!("{name}.cram"));
    if !cram.exists() {
        let st = Command::new(std::env::current_exe().unwrap())
            .args([
                "--ignored",
                "--exact",
                "cram_fixture_writer",
                "--test-threads=1",
            ])
            .env(
                "GLOCI_CRAM_WRITE",
                format!("{reads}\n{}\n{}", fa.display(), cram.display()),
            )
            .stdout(Stdio::null())
            .status()
            .unwrap();
        assert!(st.success(), "writing the CRAM failed");
        let st = Command::new(gloci_bin())
            .arg("index")
            .arg(&cram)
            .status()
            .unwrap();
        assert!(st.success(), "index CRAM failed");
    }
    (cram, fa)
}

/// The child of [`make_cram`]: `reads` reads sampled from the reference, written with noodles.
/// Does nothing unless `GLOCI_CRAM_WRITE` is set.
#[test]
#[ignore = "helper of the CRAM tests"]
fn cram_fixture_writer() {
    use noodles_sam::alignment::io::Write as _;
    let Ok(spec) = std::env::var("GLOCI_CRAM_WRITE") else {
        return;
    };
    let mut it = spec.lines();
    let reads: u64 = it.next().unwrap().parse().unwrap();
    let fa = PathBuf::from(it.next().unwrap());
    let cram = PathBuf::from(it.next().unwrap());
    let mut records = Vec::new();
    {
        let mut r = noodles_fasta::io::Reader::new(std::io::BufReader::new(
            std::fs::File::open(&fa).unwrap(),
        ));
        for rec in r.records() {
            records.push(rec.unwrap());
        }
    }
    let repo = noodles_fasta::Repository::new(records);
    let header: noodles_sam::Header = sam_header(3, CRAM_CHROM as u64, true).parse().unwrap();
    let partial = cram.with_extension("partial.cram");
    let mut w = noodles_cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repo.clone())
        .build_from_path(&partial)
        .unwrap();
    w.write_header(&header).unwrap();
    let per_chrom = reads / 3;
    let step = ((CRAM_CHROM - READ_LEN - 2) as u64 / per_chrom.max(1)).max(1);
    let mut rng = Rng(9);
    let mut i = 0u64;
    let mut text = Vec::with_capacity(512);
    for c in 0..3usize {
        let seq = repo
            .get(format!("chr{}", c + 1).as_bytes())
            .unwrap()
            .unwrap();
        let seq: &[u8] = seq.as_ref().as_ref();
        let mut pos = 1u64;
        for _ in 0..per_chrom {
            pos += rng.next() % (2 * step);
            let p = pos.min((CRAM_CHROM - READ_LEN) as u64);
            let s = &seq[(p - 1) as usize..(p - 1) as usize + READ_LEN];
            text.clear();
            write!(
                text,
                "read{i}\t0\tchr{}\t{p}\t60\t{READ_LEN}M\t*\t0\t0\t",
                c + 1
            )
            .unwrap();
            text.extend_from_slice(&s.to_ascii_uppercase());
            text.push(b'\t');
            text.extend(std::iter::repeat_n(b'I', READ_LEN));
            text.push(b'\n');
            let mut rd = noodles_sam::io::Reader::new(&text[..]);
            let mut rec = noodles_sam::Record::default();
            rd.read_record(&mut rec).unwrap();
            w.write_alignment_record(&header, &rec).unwrap();
            i += 1;
        }
    }
    w.finish(&header).unwrap();
    drop(w);
    std::fs::rename(&partial, &cram).unwrap();
}

#[test]
fn p3a_cram_generator_small() {
    let (cram, fa) = {
        // a 30 k read CRAM in its own place: make_cram writes into the shared directory
        make_cram("small", 30_000)
    };
    let o = gloci()
        .args(["convert"])
        .arg(&cram)
        .arg("--reference")
        .arg(&fa)
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(o.stdout.iter().filter(|&&b| b == b'\n').count(), 30_000);
    let _ = std::fs::remove_file(&cram);
    let _ = std::fs::remove_file(format!("{}.crai", cram.display()));
}

#[test]
#[ignore = "nightly: generates a large CRAM"]
fn convert_cram_streams() {
    let (cram, fa) = make_cram("big", payload_bytes() / 100);
    for t in ["1", "4"] {
        let rss = rss_of(
            &format!("convert CRAM -t {t}"),
            &[
                "convert".as_ref(),
                cram.as_os_str(),
                "--reference".as_ref(),
                fa.as_os_str(),
                "-t".as_ref(),
                t.as_ref(),
            ],
        );
        // `-t 1` holds the ceiling of the spec: 200 MB plus one reference sequence. With a `.crai`
        // a CRAM is read one chromosome per worker, each holding its reference (and loading it
        // takes about twice its size), so `-t 4` over 3 chromosomes is bounded by 200 MB plus
        // 4 reference sequences for each of them. Measured: 109-113 MiB at `-t 1` and 263-301
        // MiB at `-t 4`, flat from 2 M to 10 M reads.
        let ceiling = if t == "1" {
            ALN_RSS + CRAM_CHROM as u64
        } else {
            ALN_RSS + 3 * 4 * CRAM_CHROM as u64
        };
        assert!(
            rss < ceiling,
            "convert CRAM -t {t}: peak RSS {rss} bytes exceeds {ceiling}"
        );
    }
}
