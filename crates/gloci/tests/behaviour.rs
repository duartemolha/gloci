//! Cross-command behaviour: thread determinism, broken pipe, truncation,
//! atomic outputs, signals and the error-message format.
mod common;
use common::*;
use gloci_testkit::fixtures::{bgzip, vcf_to_bcf, write};
use predicates::prelude::*;
use predicates::str::contains;
use std::fmt::Write as _;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

fn bin() -> PathBuf {
    assert_cmd::cargo::cargo_bin("gloci")
}

/// Deterministic unsorted BED with `n` records over 12 chromosomes.
fn bed_text(n: usize) -> String {
    let mut s = String::with_capacity(n * 40);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let start = (x >> 33) % 5_000_000;
        writeln!(
            s,
            "chr{}\t{start}\t{}\tr{i}\t{}\t+",
            1 + (x >> 20) % 12,
            start + 50 + (x >> 50) % 500,
            i % 1000
        )
        .unwrap();
    }
    s
}

fn vcf_text(n: usize, sorted: bool) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    let mut x: u64 = 7;
    for i in 0..n {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
        let (chrom, pos) = if sorted {
            (1 + i * 2 / n, 1 + i * 20)
        } else {
            (1 + (x >> 40) as usize % 2, 1 + (x >> 20) as usize % 100_000)
        };
        writeln!(s, "chr{chrom}\t{pos}\t.\tA\tC\t.\t.\t.").unwrap();
    }
    s
}

/// Run with a thread option in the given spelling (`-t 1`, `-t2`, `--threads=8`).
fn stdout_with(args: &[&str], t: &str) -> Vec<u8> {
    let opt: Vec<String> = match t {
        "1" => vec!["-t".into(), "1".into()],
        "2" => vec!["-t2".into()],
        n => vec![format!("--threads={n}")],
    };
    let out = gloci().args(args).args(&opt).output().unwrap();
    assert_eq!(
        out.status.code(),
        Some(0),
        "{args:?} -t {t}: {:?}",
        out.stderr
    );
    out.stdout
}

fn index_bytes(src: &Path, t: &str) -> Vec<u8> {
    let d = tempdir();
    let f = d.path().join("in.vcf.gz");
    std::fs::copy(src, &f).unwrap();
    gloci()
        .args(["index", f.to_str().unwrap(), "-t", t])
        .assert()
        .success();
    std::fs::read(d.path().join("in.vcf.gz.tbi")).unwrap()
}

#[test]
fn thread_count_never_changes_bytes() {
    let d = tempdir();
    let p = |f: &PathBuf| f.to_str().unwrap().to_string();
    let bed = p(&plain_file(d.path(), "a.bed", &bed_text(100_000)));
    let sorted_bed = d.path().join("sorted.bed");
    gloci()
        .args(["sort", &bed, "-o", sorted_bed.to_str().unwrap()])
        .assert()
        .success();
    let sorted_bed = p(&sorted_bed);
    let vcf_txt = vcf_text(20_000, false);
    let vcf = p(&plain_file(d.path(), "a.vcf", &vcf_txt));
    let vcf_gz = p(&write(d.path(), "u.vcf.gz", &bgzip(vcf_txt.as_bytes())));
    let bcf = p(&write(d.path(), "a.bcf", &vcf_to_bcf(vcf_txt.as_bytes())));
    let sorted_txt = vcf_text(20_000, true);
    let sorted_gz = d.path().join("s.vcf.gz");
    write(d.path(), "s.vcf.gz", &bgzip(sorted_txt.as_bytes()));
    gloci()
        .args(["index", sorted_gz.to_str().unwrap()])
        .assert()
        .success();
    let sorted_gz = p(&sorted_gz);

    let cases: Vec<Vec<&str>> = vec![
        vec!["sort", &bed],
        vec!["sort", &vcf_gz, "-O", "vcf.gz"],
        vec!["sort", &bed, "--max-mem", "16M", "-O", "bed.gz"],
        vec!["view", &bcf, "-O", "vcf"],
        vec!["view", &vcf, "-O", "bcf"],
        vec!["view", &sorted_gz, "chr1:1-50000"],
        vec!["check", &sorted_bed],
    ];
    for args in &cases {
        let r: Vec<_> = ["1", "2", "8"]
            .iter()
            .map(|t| stdout_with(args, t))
            .collect();
        assert!(r.windows(2).all(|w| w[0] == w[1]), "{args:?}");
    }
    // sort that spills with a 16M budget
    // The input (~65 MB) is far above the 16M budget, so it must spill; the
    // default (2G) budget sorts in memory and must give the same bytes.
    let big_path = plain_file(d.path(), "big.bed", &bed_text(1_600_000));
    assert!(std::fs::metadata(&big_path).unwrap().len() > 2 * (16 << 20));
    let big = p(&big_path);
    let spill = vec!["sort", &big, "--max-mem", "16M", "-O", "bed.gz"];
    let r: Vec<_> = ["1", "8"].iter().map(|t| stdout_with(&spill, t)).collect();
    assert_eq!(r[0], r[1], "spilling sort");
    let in_mem = stdout_with(&["sort", &big, "-O", "bed.gz"], "2");
    assert_eq!(r[0], in_mem, "spilled vs in-memory");
    // index bytes
    let src = PathBuf::from(&sorted_gz);
    assert_eq!(index_bytes(&src, "1"), index_bytes(&src, "8"));
}

#[test]
fn broken_pipe_is_quiet_success() {
    let d = tempdir();
    let big = plain_file(d.path(), "big.bed", &bed_text(400_000));
    let mut child = Command::new(bin())
        .args(["view", big.to_str().unwrap()])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut first = [0u8; 16];
    child
        .stdout
        .as_mut()
        .unwrap()
        .read_exact(&mut first)
        .unwrap();
    drop(child.stdout.take());
    let out = child.wait_with_output().unwrap();
    assert_eq!(out.status.code(), Some(0));
    assert!(
        out.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

#[test]
fn truncated_bgzf_errors() {
    let d = tempdir();
    let vcf = vcf_text(50_000, true);
    let payloads: Vec<(&str, Vec<u8>)> = vec![
        ("bed.gz", bgzip(bed_text(50_000).as_bytes())),
        ("vcf.gz", bgzip(vcf.as_bytes())),
        ("bcf", vcf_to_bcf(vcf.as_bytes())),
    ];
    for (ext, full) in payloads {
        let name = format!("t.{ext}");
        // 1. cut mid-block, then a valid EOF marker: only the block is wrong
        let mut cut = full[..full.len() / 2].to_vec();
        cut.extend_from_slice(&BGZF_EOF);
        let mid = write(d.path(), &format!("mid-{name}"), &cut);
        let out = gloci()
            .args(["view", mid.to_str().unwrap()])
            .output()
            .unwrap();
        let err = String::from_utf8_lossy(&out.stderr).into_owned();
        assert_eq!(out.status.code(), Some(1), "{ext}: {err}");
        assert!(err.starts_with("error: "), "{ext}: {err}");
        assert!(!err.contains("end-of-file marker missing"), "{ext}: {err}");
        // same corruption through stdin
        gloci()
            .args(["view", "-"])
            .write_stdin(cut)
            .assert()
            .code(1)
            .stderr(predicates::str::contains("error: "));
        // 2. whole file but without the 28-byte EOF marker
        let no_eof = write(d.path(), &format!("noeof-{name}"), &full[..full.len() - 28]);
        gloci()
            .args(["view", no_eof.to_str().unwrap()])
            .assert()
            .code(1)
            .stderr(contains("end-of-file marker missing").and(contains("truncated")));
    }
}

/// Names in `dir` other than `keep`.
fn leftovers(dir: &Path, keep: &[&str]) -> Vec<String> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| !keep.contains(&n.as_str()))
        .collect()
}

#[test]
fn no_partial_outputs_on_error() {
    let d = tempdir();
    let mut text = bed_text(50_000);
    let at = text.len() * 9 / 10;
    let at = text[at..].find('\n').unwrap() + at + 1;
    text.insert_str(at, "chr1\tnotanumber\t5\n");
    let input = plain_file(d.path(), "bad.bed", &text);
    let out = d.path().join("out.bed");
    gloci()
        .args(["sort", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .assert()
        .code(1);
    assert!(!out.exists());
    assert_eq!(leftovers(d.path(), &["bad.bed"]), Vec::<String>::new());
}

#[test]
fn no_partial_outputs_when_error_follows_spilling() {
    let d = tempdir();
    let mut text = bed_text(1_500_000);
    let at = text.len() * 9 / 10;
    let at = text[at..].find('\n').unwrap() + at + 1;
    let line_no = text[..at].matches('\n').count() + 1;
    text.insert_str(at, "chr1\tnotanumber\t5\n");
    assert!(text.len() > 2 * (16 << 20));
    let input = plain_file(d.path(), "bad.bed", &text);
    let spill = d.path().join("spill");
    std::fs::create_dir(&spill).unwrap();
    let out = d.path().join("out.bed");
    gloci()
        .current_dir(d.path())
        .args(["sort", "bad.bed", "--max-mem", "16M", "--tmp-dir"])
        .arg(&spill)
        .arg("-o")
        .arg(&out)
        .assert()
        .code(1)
        .stderr(contains(format!(
            "bad.bed line {}",
            gloci_core::thousands(line_no as u64)
        )));
    assert!(!out.exists());
    assert_eq!(
        leftovers(d.path(), &["bad.bed", "spill"]),
        Vec::<String>::new()
    );
    assert!(std::fs::read_dir(&spill).unwrap().next().is_none());
    drop(input);
}

fn wait_until(what: &str, mut f: impl FnMut() -> bool) {
    let t = Instant::now();
    while !f() {
        assert!(t.elapsed() < Duration::from_secs(60), "timed out: {what}");
        std::thread::sleep(Duration::from_millis(1));
    }
}

#[test]
fn killed_mid_write_leaves_nothing() {
    let d = tempdir();
    let big = plain_file(d.path(), "big.bed", &bed_text(3_000_000));
    let tmp = d.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    let out = d.path().join("out.bed");
    let mut c = Command::new(bin())
        .args([
            "sort",
            big.to_str().unwrap(),
            "--max-mem",
            "16M",
            "--tmp-dir",
        ])
        .arg(&tmp)
        .arg("-o")
        .arg(&out)
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    wait_until("spilling to start", || {
        std::fs::read_dir(&tmp).unwrap().next().is_some()
    });
    assert!(
        c.try_wait().unwrap().is_none(),
        "sort finished before the signal"
    );
    unsafe { libc::kill(c.id() as i32, libc::SIGTERM) };
    let st = c.wait().unwrap();
    assert_eq!(st.code(), Some(130));
    assert!(!out.exists());
    assert_eq!(
        leftovers(d.path(), &["big.bed", "spill"]),
        Vec::<String>::new()
    );
    assert!(std::fs::read_dir(&tmp).unwrap().next().is_none());
}

#[test]
fn missing_index_region_query_warns() {
    let d = tempdir();
    let f = write(d.path(), "u.vcf.gz", &bgzip(vcf_text(100, true).as_bytes()));
    let name = f.to_str().unwrap();
    gloci()
        .args(["view", name, "chr1:1-1000"])
        .assert()
        .code(0)
        .stderr(contains(format!("warning: no index for {name}")));
}

#[test]
fn error_format() {
    let d = tempdir();
    let text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t12a4\t.\tA\tC\t.\t.\t.\n";
    let p = plain_file(d.path(), "bad.vcf", text);
    gloci()
        .current_dir(d.path())
        .args(["sort", p.file_name().unwrap().to_str().unwrap()])
        .assert()
        .code(1)
        .stderr("error: bad.vcf line 3: POS \"12a4\" is not an integer\n");
}

#[test]
fn error_format_thousands_separator_in_line_number() {
    let d = tempdir();
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for i in 0..1231 {
        writeln!(text, "chr1\t{}\t.\tA\tC\t.\t.\t.", i + 1).unwrap();
    }
    text.push_str("chr1\t12a4\t.\tA\tC\t.\t.\t.\n");
    assert_eq!(text.matches('\n').count(), 1234);
    plain_file(d.path(), "bad.vcf", &text);
    gloci()
        .current_dir(d.path())
        .args(["sort", "bad.vcf"])
        .assert()
        .code(1)
        .stderr("error: bad.vcf line 1,234: POS \"12a4\" is not an integer\n");
}

#[test]
fn vcf_output_ignores_thread_option_spelling() {
    let d = tempdir();
    let f = plain_file(d.path(), "x.vcf", &vcf_text(500, false));
    let x = f.to_str().unwrap();
    // -q and -v stay in the recorded command, so compare within a flag set
    let groups: [&[&[&str]]; 3] = [
        &[
            &[],
            &["-t", "1"],
            &["-t=2"],
            &["-t2"],
            &["--threads", "3"],
            &["--threads=4"],
        ],
        &[&["-q"], &["-qt2"], &["-qt", "3"], &["-qt=4"]],
        &[&["-v"], &["-vt2"], &["-vt=4"]],
    ];
    for g in groups {
        let outs: Vec<Vec<u8>> = g
            .iter()
            .map(|t| {
                let o = gloci().args(["sort", x]).args(*t).output().unwrap();
                assert_eq!(o.status.code(), Some(0), "{t:?}");
                o.stdout
            })
            .collect();
        assert!(outs.windows(2).all(|w| w[0] == w[1]), "{g:?}");
    }
}
