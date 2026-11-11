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
        "{args:?} -t {t}: {}",
        String::from_utf8_lossy(&out.stderr)
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

#[test]
fn determinism_matrix() {
    let d = tempdir();
    let p = |f: &PathBuf| f.to_str().unwrap().to_string();

    let genome_txt = "chr1\t10000\nchr2\t10000\nchr3\t10000\nchr4\t10000\n";
    let genome = p(&plain_file(d.path(), "g.genome", genome_txt));

    // Multi-chromosome sorted BED files (chr1..chr4)
    let mut a_text = String::new();
    let mut b_text = String::new();
    let mut bed12_text = String::new();
    for c in 1..=4 {
        for i in 0..50 {
            let start = 100 + i * 40;
            let end = start + 30;
            writeln!(
                a_text,
                "chr{c}\t{start}\t{end}\ta_{c}_{i}\t{}\t+",
                10 + i % 100
            )
            .unwrap();
            writeln!(
                bed12_text,
                "chr{c}\t{start}\t{end}\ta_{c}_{i}\t{}\t+\t{start}\t{end}\t0\t2\t10,10,\t0,20,",
                10 + i % 100
            )
            .unwrap();
            let b_start = 110 + i * 40;
            let b_end = b_start + 35;
            writeln!(
                b_text,
                "chr{c}\t{b_start}\t{b_end}\tb_{c}_{i}\t{}\t+",
                20 + i % 100
            )
            .unwrap();
        }
    }
    let a_plain = p(&plain_file(d.path(), "a.bed", &a_text));
    let b_plain = p(&plain_file(d.path(), "b.bed", &b_text));

    let a_bgzf = p(&bgzf_file(d.path(), "a.bed.gz", &a_text));
    let b_bgzf = p(&bgzf_file(d.path(), "b.bed.gz", &b_text));

    let a_idx = d.path().join("a_idx.bed.gz");
    std::fs::copy(&a_bgzf, &a_idx).unwrap();
    gloci()
        .args(["index", a_idx.to_str().unwrap()])
        .assert()
        .success();
    let a_idx = p(&a_idx);

    let b_idx = d.path().join("b_idx.bed.gz");
    std::fs::copy(&b_bgzf, &b_idx).unwrap();
    gloci()
        .args(["index", b_idx.to_str().unwrap()])
        .assert()
        .success();
    let b_idx = p(&b_idx);

    // Multi-chromosome unsorted BED files for --unsorted relation modes
    let mut ua_text = String::new();
    let mut ub_text = String::new();
    for i in (0..50).rev() {
        for c in (1..=4).rev() {
            let start = 100 + i * 40;
            let end = start + 30;
            writeln!(
                ua_text,
                "chr{c}\t{start}\t{end}\tua_{c}_{i}\t{}\t+",
                10 + i % 100
            )
            .unwrap();
            let b_start = 110 + i * 40;
            let b_end = b_start + 35;
            writeln!(
                ub_text,
                "chr{c}\t{b_start}\t{b_end}\tub_{c}_{i}\t{}\t+",
                20 + i % 100
            )
            .unwrap();
        }
    }
    let ua_plain = p(&plain_file(d.path(), "ua.bed", &ua_text));
    let ub_plain = p(&plain_file(d.path(), "ub.bed", &ub_text));
    let ua_bgzf = p(&bgzf_file(d.path(), "ua.bed.gz", &ua_text));
    let ub_bgzf = p(&bgzf_file(d.path(), "ub.bed.gz", &ub_text));

    // BED12 plain, BGZF, and indexed BGZF
    let bed12_plain = p(&plain_file(d.path(), "a12.bed", &bed12_text));
    let bed12_bgzf = p(&bgzf_file(d.path(), "a12.bed.gz", &bed12_text));
    let bed12_idx = d.path().join("a12_idx.bed.gz");
    std::fs::copy(&bed12_bgzf, &bed12_idx).unwrap();
    gloci()
        .args(["index", bed12_idx.to_str().unwrap()])
        .assert()
        .success();
    let bed12_idx = p(&bed12_idx);

    // Multi-chromosome sorted GFF3 for convert
    let mut gff_text = String::from("##gff-version 3\n");
    for c in 1..=4 {
        for i in 0..50 {
            let start = 101 + i * 40;
            let end = start + 29;
            writeln!(
                gff_text,
                "chr{c}\tsrc\tgene\t{start}\t{end}\t{}\t+\t.\tID=g_{c}_{i};Name=Gene_{c}_{i}",
                10 + i % 100
            )
            .unwrap();
        }
    }
    let gff_plain = p(&plain_file(d.path(), "a.gff3", &gff_text));
    let gff_bgzf = p(&bgzf_file(d.path(), "a.gff3.gz", &gff_text));
    let gff_idx = d.path().join("a_idx.gff3.gz");
    std::fs::copy(&gff_bgzf, &gff_idx).unwrap();
    gloci()
        .args(["index", gff_idx.to_str().unwrap()])
        .assert()
        .success();
    let gff_idx = p(&gff_idx);

    // Multi-chromosome sorted VCF
    let mut vcf_text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n##contig=<ID=chr3>\n##contig=<ID=chr4>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for c in 1..=4 {
        for i in 0..30 {
            let pos = 100 + i * 50;
            let qual = 15 + i % 60;
            writeln!(
                vcf_text,
                "chr{c}\t{pos}\trs_{c}_{i}\tA\tC\t{qual}\tPASS\tDP=30"
            )
            .unwrap();
        }
    }
    let vcf_plain = p(&plain_file(d.path(), "a.vcf", &vcf_text));
    let vcf_bgzf = p(&bgzf_file(d.path(), "a.vcf.gz", &vcf_text));
    let vcf_idx = d.path().join("a_idx.vcf.gz");
    std::fs::copy(&vcf_bgzf, &vcf_idx).unwrap();
    gloci()
        .args(["index", vcf_idx.to_str().unwrap()])
        .assert()
        .success();
    let vcf_idx = p(&vcf_idx);

    // Command templates: each takes (a, b) and returns arguments
    enum Kind<'a> {
        BedRel(&'a [&'a str]),
        BedUnsortedRel(&'a [&'a str]),
        BedSingle(&'a [&'a str]),
        Bed12Single(&'a [&'a str]),
        GffSingle(&'a [&'a str]),
        GenomeSingle(&'a [&'a str]),
        VcfSingle(&'a [&'a str]),
    }

    let cases: &[(&str, Kind)] = &[
        // intersect modes: any/pairs/count/clip/wao/left-outer/v
        (
            "intersect any",
            Kind::BedRel(&["intersect", "-a", "{a}", "-b", "{b}", "--mode", "any"]),
        ),
        (
            "intersect pairs",
            Kind::BedRel(&["intersect", "-a", "{a}", "-b", "{b}", "--mode", "pairs"]),
        ),
        (
            "intersect count",
            Kind::BedRel(&["intersect", "-a", "{a}", "-b", "{b}", "--mode", "count"]),
        ),
        (
            "intersect clip",
            Kind::BedRel(&["intersect", "-a", "{a}", "-b", "{b}", "--mode", "clip"]),
        ),
        (
            "intersect wao",
            Kind::BedRel(&[
                "intersect",
                "-a",
                "{a}",
                "-b",
                "{b}",
                "--mode",
                "pairs",
                "--overlap-bases",
                "--keep-unmatched",
            ]),
        ),
        (
            "intersect left-outer",
            Kind::BedRel(&[
                "intersect",
                "-a",
                "{a}",
                "-b",
                "{b}",
                "--mode",
                "pairs",
                "--keep-unmatched",
            ]),
        ),
        (
            "intersect v",
            Kind::BedRel(&["intersect", "-a", "{a}", "-b", "{b}", "--mode", "none"]),
        ),
        (
            "intersect unsorted",
            Kind::BedUnsortedRel(&["intersect", "-a", "{a}", "-b", "{b}", "--unsorted"]),
        ),
        // merge
        ("merge", Kind::BedSingle(&["merge", "{a}"])),
        (
            "merge agg",
            Kind::BedSingle(&["merge", "{a}", "-c", "5:sum"]),
        ),
        // subtract: piece/whole/any/unsorted
        (
            "subtract piece",
            Kind::BedRel(&["subtract", "-a", "{a}", "-b", "{b}"]),
        ),
        (
            "subtract whole",
            Kind::BedRel(&["subtract", "-a", "{a}", "-b", "{b}", "--whole"]),
        ),
        (
            "subtract any",
            Kind::BedRel(&[
                "subtract",
                "-a",
                "{a}",
                "-b",
                "{b}",
                "--sum-overlaps",
                "-f",
                "0.5",
            ]),
        ),
        (
            "subtract unsorted",
            Kind::BedUnsortedRel(&["subtract", "-a", "{a}", "-b", "{b}", "--unsorted"]),
        ),
        // closest: standard/k/d/unsorted
        (
            "closest standard",
            Kind::BedRel(&["closest", "-a", "{a}", "-b", "{b}"]),
        ),
        (
            "closest k",
            Kind::BedRel(&["closest", "-a", "{a}", "-b", "{b}", "-k", "2"]),
        ),
        (
            "closest d",
            Kind::BedRel(&["closest", "-a", "{a}", "-b", "{b}", "--distance"]),
        ),
        (
            "closest unsorted",
            Kind::BedUnsortedRel(&["closest", "-a", "{a}", "-b", "{b}", "--unsorted"]),
        ),
        // map: standard/unsorted
        (
            "map",
            Kind::BedRel(&["map", "-a", "{a}", "-b", "{b}", "-c", "5:sum"]),
        ),
        (
            "map unsorted",
            Kind::BedUnsortedRel(&["map", "-a", "{a}", "-b", "{b}", "-c", "5:sum", "--unsorted"]),
        ),
        // coverage: summary/hist/per-base/counts/unsorted
        (
            "coverage summary",
            Kind::BedRel(&["coverage", "-a", "{a}", "-b", "{b}"]),
        ),
        (
            "coverage hist",
            Kind::BedRel(&["coverage", "-a", "{a}", "-b", "{b}", "--hist"]),
        ),
        (
            "coverage per-base",
            Kind::BedRel(&["coverage", "-a", "{a}", "-b", "{b}", "--per-base"]),
        ),
        (
            "coverage counts",
            Kind::BedRel(&["coverage", "-a", "{a}", "-b", "{b}", "--mean"]),
        ),
        (
            "coverage unsorted",
            Kind::BedUnsortedRel(&["coverage", "-a", "{a}", "-b", "{b}", "--unsorted"]),
        ),
        // complement
        ("complement", Kind::BedSingle(&["complement", "{a}"])),
        // slop, flank, windows, convert
        (
            "slop",
            Kind::BedSingle(&["slop", "{a}", "-g", &genome, "-b", "10"]),
        ),
        (
            "flank",
            Kind::BedSingle(&["flank", "{a}", "-g", &genome, "-b", "10"]),
        ),
        (
            "windows -g",
            Kind::GenomeSingle(&["windows", "-g", &genome, "-w", "100"]),
        ),
        (
            "windows -i",
            Kind::BedSingle(&["windows", "-i", "{a}", "-w", "5"]),
        ),
        (
            "convert bed12-to-bed6",
            Kind::Bed12Single(&["convert", "{a}", "--blocks"]),
        ),
        ("convert gff3-to-bed", Kind::GffSingle(&["convert", "{a}"])),
        // vcf filter
        (
            "vcf filter",
            Kind::VcfSingle(&["vcf", "filter", "{a}", "--min-qual", "30"]),
        ),
    ];

    let strip_gloci_command = |bytes: &[u8]| -> Vec<u8> {
        String::from_utf8_lossy(bytes)
            .lines()
            .filter(|l| !l.starts_with("##glociCommand="))
            .collect::<Vec<_>>()
            .join("\n")
            .into_bytes()
    };

    for (name, kind) in cases {
        let (file_sets, is_vcf): (Vec<(&str, &str)>, bool) = match kind {
            Kind::BedRel(_) => (
                vec![(&a_plain, &b_plain), (&a_bgzf, &b_bgzf), (&a_idx, &b_idx)],
                false,
            ),
            Kind::BedUnsortedRel(_) => (vec![(&ua_plain, &ub_plain), (&ua_bgzf, &ub_bgzf)], false),
            Kind::BedSingle(_) => (vec![(&a_plain, ""), (&a_bgzf, ""), (&a_idx, "")], false),
            Kind::Bed12Single(_) => (
                vec![(&bed12_plain, ""), (&bed12_bgzf, ""), (&bed12_idx, "")],
                false,
            ),
            Kind::GffSingle(_) => (
                vec![(&gff_plain, ""), (&gff_bgzf, ""), (&gff_idx, "")],
                false,
            ),
            Kind::GenomeSingle(_) => (vec![(&genome, "")], false),
            Kind::VcfSingle(_) => (
                vec![(&vcf_plain, ""), (&vcf_bgzf, ""), (&vcf_idx, "")],
                true,
            ),
        };

        let mut format_outputs = Vec::new();

        for (a_file, b_file) in file_sets {
            let template = match kind {
                Kind::BedRel(tmpl)
                | Kind::BedUnsortedRel(tmpl)
                | Kind::BedSingle(tmpl)
                | Kind::Bed12Single(tmpl)
                | Kind::GffSingle(tmpl)
                | Kind::GenomeSingle(tmpl)
                | Kind::VcfSingle(tmpl) => tmpl,
            };
            let args: Vec<String> = template
                .iter()
                .map(|arg| {
                    if *arg == "{a}" {
                        a_file.to_string()
                    } else if *arg == "{b}" {
                        b_file.to_string()
                    } else {
                        arg.to_string()
                    }
                })
                .collect();
            let str_args: Vec<&str> = args.iter().map(|s| s.as_str()).collect();

            let t_outputs: Vec<Vec<u8>> = ["1", "2", "8"]
                .iter()
                .map(|t| stdout_with(&str_args, t))
                .collect();

            assert_eq!(
                t_outputs[0], t_outputs[1],
                "{name} at {a_file}: -t 1 vs -t 2"
            );
            assert_eq!(
                t_outputs[1], t_outputs[2],
                "{name} at {a_file}: -t 2 vs -t 8"
            );

            format_outputs.push(t_outputs[0].clone());
        }

        if is_vcf {
            for w in format_outputs.windows(2) {
                assert_eq!(
                    strip_gloci_command(&w[0]),
                    strip_gloci_command(&w[1]),
                    "{name}: plain vs bgzf vs indexed bgzf"
                );
            }
        } else {
            for w in format_outputs.windows(2) {
                assert_eq!(w[0], w[1], "{name}: plain vs bgzf vs indexed bgzf");
            }
        }
    }
}

#[test]
fn broken_pipe_matrix() {
    let d = tempdir();
    let genome = plain_file(d.path(), "g.genome", "chr1\t1000000\n");
    let mut a_text = String::new();
    let mut b_text = String::new();
    let mut bed12_text = String::new();
    for i in 0..10_000 {
        let s = i * 20;
        let e = s + 15;
        writeln!(a_text, "chr1\t{s}\t{e}\ta{i}\t100\t+").unwrap();
        writeln!(b_text, "chr1\t{}\t{}\tb{i}\t50\t+", s + 5, s + 25).unwrap();
        writeln!(
            bed12_text,
            "chr1\t{s}\t{e}\ta{i}\t100\t+\t{s}\t{e}\t0\t2\t5,5,\t0,10,"
        )
        .unwrap();
    }
    let a_bed = plain_file(d.path(), "a.bed", &a_text);
    let b_bed = plain_file(d.path(), "b.bed", &b_text);
    let bed12 = plain_file(d.path(), "a12.bed", &bed12_text);

    let mut vcf_text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..10_000 {
        writeln!(
            vcf_text,
            "chr1\t{}\trs{i}\tA\tC\t50\tPASS\tDP=30",
            1 + i * 20
        )
        .unwrap();
    }
    let vcf = plain_file(d.path(), "v.vcf", &vcf_text);

    let genome_str = genome.to_str().unwrap();
    let a_str = a_bed.to_str().unwrap();
    let b_str = b_bed.to_str().unwrap();
    let bed12_str = bed12.to_str().unwrap();
    let vcf_str = vcf.to_str().unwrap();

    let cmds: &[&[&str]] = &[
        &["intersect", "-a", a_str, "-b", b_str],
        &["merge", a_str],
        &["subtract", "-a", a_str, "-b", b_str],
        &["closest", "-a", a_str, "-b", b_str],
        &["map", "-a", a_str, "-b", b_str, "-c", "5:sum"],
        &["coverage", "-a", a_str, "-b", b_str],
        &["complement", a_str],
        &["slop", a_str, "-g", genome_str, "-b", "2"],
        &["flank", a_str, "-g", genome_str, "-b", "2"],
        &["windows", "-g", genome_str, "-w", "10"],
        &["windows", "-i", a_str, "-w", "5"],
        &["convert", bed12_str, "--blocks"],
        &["convert", vcf_str],
        &["vcf", "filter", vcf_str],
    ];

    for args in cmds {
        let mut child = Command::new(bin())
            .args(*args)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut first = [0u8; 64];
        let _ = child.stdout.as_mut().unwrap().read(&mut first);
        drop(child.stdout.take());
        let out = child.wait_with_output().unwrap();
        assert_eq!(
            out.status.code(),
            Some(0),
            "command {args:?} exited with code {:?}, stderr: {}",
            out.status.code(),
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            out.stderr.is_empty(),
            "command {args:?} printed on stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
}

#[test]
fn empty_inputs_matrix() {
    let d = tempdir();
    let genome = plain_file(d.path(), "g.genome", "chr1\t1000\nchr2\t1000\n");
    let empty_bed = plain_file(d.path(), "empty.bed", "");
    let a_bed = plain_file(
        d.path(),
        "a.bed",
        "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n",
    );
    let b_bed = plain_file(
        d.path(),
        "b.bed",
        "chr1\t20\t40\tb1\t50\t+\nchr1\t80\t120\tb2\t60\t+\n",
    );
    let vcf_head =
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let empty_vcf = plain_file(d.path(), "empty.vcf", vcf_head);
    let b_diff = plain_file(d.path(), "b_diff.bed", "chr2\t20\t40\tb1\t50\t+\n");

    let genome_str = genome.to_str().unwrap();
    let empty_bed_str = empty_bed.to_str().unwrap();
    let a_str = a_bed.to_str().unwrap();
    let b_str = b_bed.to_str().unwrap();
    let empty_vcf_str = empty_vcf.to_str().unwrap();
    let b_diff_str = b_diff.to_str().unwrap();

    // 1. Empty A -> exits 0 with empty or header-only output
    let empty_a_cmds: &[&[&str]] = &[
        &["intersect", "-a", empty_bed_str, "-b", b_str],
        &["merge", empty_bed_str],
        &["subtract", "-a", empty_bed_str, "-b", b_str],
        &["closest", "-a", empty_bed_str, "-b", b_str],
        &["map", "-a", empty_bed_str, "-b", b_str, "-c", "5:sum"],
        &["coverage", "-a", empty_bed_str, "-b", b_str],
        &["complement", empty_bed_str],
        &["slop", empty_bed_str, "-g", genome_str, "-b", "10"],
        &["flank", empty_bed_str, "-g", genome_str, "-b", "10"],
        &["windows", "-i", empty_bed_str, "-w", "10"],
        &["convert", empty_bed_str, "--blocks"],
        &["convert", empty_vcf_str],
    ];
    for args in empty_a_cmds {
        let o = gloci().args(*args).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{args:?}");
        assert!(o.stdout.is_empty(), "{args:?}: expected empty stdout");
    }
    // empty VCF A through vcf filter
    let o = gloci()
        .args(["vcf", "filter", empty_vcf_str])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    let s = String::from_utf8(o.stdout).unwrap();
    assert!(s.starts_with("##fileformat=VCFv4.2\n"));
    assert!(s.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert_eq!(s.lines().filter(|l| !l.starts_with('#')).count(), 0);

    // 2. Empty B -> exits 0 with expected output
    // intersect --mode any: empty
    let o = gloci()
        .args([
            "intersect",
            "-a",
            a_str,
            "-b",
            empty_bed_str,
            "--mode",
            "any",
        ])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert!(o.stdout.is_empty());

    // intersect --mode none: all A
    let o = gloci()
        .args([
            "intersect",
            "-a",
            a_str,
            "-b",
            empty_bed_str,
            "--mode",
            "none",
        ])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n"
    );

    // intersect --mode count: A with 0
    let o = gloci()
        .args([
            "intersect",
            "-a",
            a_str,
            "-b",
            empty_bed_str,
            "--mode",
            "count",
        ])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "chr1\t10\t50\ta1\t100\t+\t0\nchr1\t60\t100\ta2\t200\t-\t0\n"
    );

    // subtract: all A
    let o = gloci()
        .args(["subtract", "-a", a_str, "-b", empty_bed_str])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n"
    );

    // closest: A with placeholder B
    let o = gloci()
        .args(["closest", "-a", a_str, "-b", empty_bed_str])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "chr1\t10\t50\ta1\t100\t+\t.\t-1\t-1\nchr1\t60\t100\ta2\t200\t-\t.\t-1\t-1\n"
    );

    // map: A with empty aggregation value (.)
    let o = gloci()
        .args(["map", "-a", a_str, "-b", empty_bed_str, "-c", "5:sum"])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "chr1\t10\t50\ta1\t100\t+\t.\nchr1\t60\t100\ta2\t200\t-\t.\n"
    );

    // coverage: A with 0 hits
    let o = gloci()
        .args(["coverage", "-a", a_str, "-b", empty_bed_str])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        "chr1\t10\t50\ta1\t100\t+\t0\t0\t40\t0\nchr1\t60\t100\ta2\t200\t-\t0\t0\t40\t0\n"
    );

    // 3. Header-only VCF A -> exits 0 with header-only output
    let vcf_passthrough_cmds: &[&[&str]] = &[
        &[
            "intersect",
            "-a",
            empty_vcf_str,
            "-b",
            b_str,
            "--mode",
            "any",
        ],
        &["subtract", "-a", empty_vcf_str, "-b", b_str],
        &["vcf", "filter", empty_vcf_str],
    ];
    for args in vcf_passthrough_cmds {
        let o = gloci().args(*args).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{args:?}");
        let s = String::from_utf8(o.stdout).unwrap();
        assert!(s.starts_with("##fileformat=VCFv4.2\n"), "{args:?}: {s}");
        assert!(
            s.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"),
            "{args:?}: {s}"
        );
        assert_eq!(
            s.lines().filter(|l| !l.starts_with('#')).count(),
            0,
            "{args:?}"
        );
    }

    let vcf_tsv_cmds: &[&[&str]] = &[
        &[
            "intersect",
            "-a",
            empty_vcf_str,
            "-b",
            b_str,
            "--mode",
            "pairs",
        ],
        &[
            "intersect",
            "-a",
            empty_vcf_str,
            "-b",
            b_str,
            "--mode",
            "count",
        ],
        &["closest", "-a", empty_vcf_str, "-b", b_str],
        &["map", "-a", empty_vcf_str, "-b", b_str, "-c", "5:sum"],
        &["coverage", "-a", empty_vcf_str, "-b", b_str],
    ];
    for args in vcf_tsv_cmds {
        let o = gloci().args(*args).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{args:?}");
        let s = String::from_utf8(o.stdout).unwrap();
        assert!(!s.starts_with("##"), "{args:?}: {s}");
        assert!(
            s.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
            "{args:?}: {s}"
        );
        assert_eq!(
            s.lines().filter(|l| !l.starts_with('#')).count(),
            0,
            "{args:?}"
        );
    }

    // 4. B with no records on A's chromosomes -> exits 0 with exact expected output
    let diff_cases: &[(&[&str], &str)] = &[
        (
            &["intersect", "-a", a_str, "-b", b_diff_str, "--mode", "any"],
            "",
        ),
        (
            &["intersect", "-a", a_str, "-b", b_diff_str, "--mode", "none"],
            "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n",
        ),
        (
            &[
                "intersect",
                "-a",
                a_str,
                "-b",
                b_diff_str,
                "--mode",
                "count",
            ],
            "chr1\t10\t50\ta1\t100\t+\t0\nchr1\t60\t100\ta2\t200\t-\t0\n",
        ),
        (
            &["subtract", "-a", a_str, "-b", b_diff_str],
            "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n",
        ),
        (
            &["closest", "-a", a_str, "-b", b_diff_str],
            "chr1\t10\t50\ta1\t100\t+\t.\t-1\t-1\t.\t-1\t.\nchr1\t60\t100\ta2\t200\t-\t.\t-1\t-1\t.\t-1\t.\n",
        ),
        (
            &["map", "-a", a_str, "-b", b_diff_str, "-c", "5:sum"],
            "chr1\t10\t50\ta1\t100\t+\t.\nchr1\t60\t100\ta2\t200\t-\t.\n",
        ),
        (
            &["coverage", "-a", a_str, "-b", b_diff_str],
            "chr1\t10\t50\ta1\t100\t+\t0\t0\t40\t0\nchr1\t60\t100\ta2\t200\t-\t0\t0\t40\t0\n",
        ),
    ];
    for (args, expected) in diff_cases {
        let o = gloci().args(*args).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{args:?}");
        assert_eq!(String::from_utf8(o.stdout).unwrap(), *expected, "{args:?}");
    }
}

#[test]
fn vcf_a_matrix() {
    let d = tempdir();
    let vcf_head =
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let vcf_txt = format!("{vcf_head}chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\n");
    let a_vcf = plain_file(d.path(), "a.vcf", &vcf_txt);
    let b_bed = plain_file(d.path(), "b.bed", "chr1\t90\t110\tb1\t50\t+\n");

    let a_str = a_vcf.to_str().unwrap();
    let b_str = b_bed.to_str().unwrap();

    let passthrough: &[&[&str]] = &[
        &["intersect", "-a", a_str, "-b", b_str, "--mode", "any"],
        &["subtract", "-a", a_str, "-b", b_str],
    ];
    for args in passthrough {
        let o = gloci().args(*args).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{args:?}");
        let s = String::from_utf8(o.stdout).unwrap();
        assert!(s.starts_with("##fileformat=VCFv4.2\n"), "{args:?}");
        assert!(s.contains("\n##glociCommand="), "{args:?}");
        assert!(
            s.contains("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"),
            "{args:?}"
        );
    }

    let tsv_modes: &[&[&str]] = &[
        &["intersect", "-a", a_str, "-b", b_str, "--mode", "pairs"],
        &["intersect", "-a", a_str, "-b", b_str, "--mode", "count"],
        &["closest", "-a", a_str, "-b", b_str],
        &["map", "-a", a_str, "-b", b_str, "-c", "5:sum"],
        &["coverage", "-a", a_str, "-b", b_str],
    ];
    for args in tsv_modes {
        let o = gloci().args(*args).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{args:?}");
        let s = String::from_utf8(o.stdout).unwrap();
        assert!(
            !s.starts_with("##"),
            "{args:?}: expected TSV without ##, got {s}"
        );
        let first_line = s.lines().next().unwrap();
        assert!(
            first_line.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"),
            "{args:?}: {first_line}"
        );
        let second_line = s.lines().nth(1).unwrap();
        assert!(
            second_line.starts_with("chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\t"),
            "{args:?}: {second_line}"
        );
    }
}

#[test]
fn failure_matrix() {
    let d = tempdir();
    let sorted_a = plain_file(
        d.path(),
        "sorted_a.bed",
        "chr1\t10\t20\tr1\t10\t+\nchr1\t200\t300\tr2\t20\t+\n",
    );
    let unsorted_a = plain_file(
        d.path(),
        "unsorted_a.bed",
        "chr1\t100\t200\tr1\t10\t+\nchr1\t50\t150\tr2\t20\t+\n",
    );
    let sorted_b = plain_file(
        d.path(),
        "sorted_b.bed",
        "chr1\t10\t20\tb1\t10\t+\nchr1\t100\t200\tb2\t20\t+\nchr1\t250\t350\tb3\t30\t+\n",
    );
    let unsorted_b = plain_file(
        d.path(),
        "unsorted_b.bed",
        "chr1\t100\t200\tb1\t10\t+\nchr1\t50\t150\tb2\t20\t+\n",
    );
    let malformed_b = plain_file(
        d.path(),
        "malformed_b.bed",
        "chr1\t100\t200\tb1\t10\t+\nchr1\tNOT_A_NUM\t250\tb2\t20\t+\n",
    );
    let full_bgzf = bgzip(b"chr1\t10\t20\tb1\t10\t+\nchr1\t100\t200\tb2\t20\t+\n");
    let trunc_bgzf = d.path().join("trunc.bed.gz");
    std::fs::write(&trunc_bgzf, &full_bgzf[..full_bgzf.len() / 2]).unwrap();

    let sorted_a_str = sorted_a.to_str().unwrap();
    let unsorted_a_str = unsorted_a.to_str().unwrap();
    let sorted_b_str = sorted_b.to_str().unwrap();
    let unsorted_b_str = unsorted_b.to_str().unwrap();
    let malformed_b_str = malformed_b.to_str().unwrap();
    let trunc_bgzf_str = trunc_bgzf.to_str().unwrap();

    let assert_no_tmp = |dir: &Path| {
        let tmps: Vec<String> = std::fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .filter(|n| n.contains(".tmp-"))
            .collect();
        assert!(tmps.is_empty(), "leftover .tmp- files in {dir:?}: {tmps:?}");
    };

    let commands: &[&[&str]] = &[
        &["intersect", "-a", "{a}", "-b", "{b}"],
        &["subtract", "-a", "{a}", "-b", "{b}"],
        &["closest", "-a", "{a}", "-b", "{b}"],
        &["map", "-a", "{a}", "-b", "{b}", "-c", "5:sum"],
        &["coverage", "-a", "{a}", "-b", "{b}"],
    ];

    let mut case_counter = 0;
    for base_cmd in commands {
        // 1. Unsorted A
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| {
                    if s == "{a}" {
                        unsorted_a_str
                    } else if s == "{b}" {
                        sorted_b_str
                    } else {
                        s
                    }
                })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("is not sorted") && err.contains("run: gloci sort"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 2. Unsorted B
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| {
                    if s == "{a}" {
                        sorted_a_str
                    } else if s == "{b}" {
                        unsorted_b_str
                    } else {
                        s
                    }
                })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("is not sorted") && err.contains("run: gloci sort"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 3. Malformed B line
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| {
                    if s == "{a}" {
                        sorted_a_str
                    } else if s == "{b}" {
                        malformed_b_str
                    } else {
                        s
                    }
                })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("is not an integer") || err.contains("format"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 4. Truncated BGZF
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| {
                    if s == "{a}" {
                        sorted_a_str
                    } else if s == "{b}" {
                        trunc_bgzf_str
                    } else {
                        s
                    }
                })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("marker missing")
                    || err.contains("truncated")
                    || err.contains("error: "),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 5. Output to a missing directory
        {
            let out_p = d.path().join("missing_dir/subdir/out.bed");
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| {
                    if s == "{a}" {
                        sorted_a_str
                    } else if s == "{b}" {
                        sorted_b_str
                    } else {
                        s
                    }
                })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("No such file or directory") || err.contains("cannot create"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }
    }

    // Single-input commands: merge and complement
    let single_cmds: &[&[&str]] = &[&["merge", "{a}"], &["complement", "{a}"]];
    for base_cmd in single_cmds {
        // 1. Unsorted input
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| if s == "{a}" { unsorted_a_str } else { s })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("is not sorted") && err.contains("gloci sort"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 2. Malformed input line
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| if s == "{a}" { malformed_b_str } else { s })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("is not an integer") || err.contains("format"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 3. Truncated BGZF
        {
            case_counter += 1;
            let out_p = d.path().join(format!("out_{case_counter}.bed"));
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| if s == "{a}" { trunc_bgzf_str } else { s })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("marker missing")
                    || err.contains("truncated")
                    || err.contains("error: "),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }

        // 4. Output to a missing directory
        {
            let out_p = d.path().join("missing_dir/subdir/out.bed");
            let args: Vec<&str> = base_cmd
                .iter()
                .map(|&s| if s == "{a}" { sorted_a_str } else { s })
                .collect();
            let o = gloci()
                .args(&args)
                .args(["-o", out_p.to_str().unwrap()])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(1), "{args:?}");
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(
                err.contains("No such file or directory") || err.contains("cannot create"),
                "{err}"
            );
            assert!(!out_p.exists(), "partial output file left at {:?}", out_p);
            assert_no_tmp(d.path());
        }
    }
}

// ---------------------------------------------------------------------------
// P2 commands: shift, cluster, spacing, groupby, partition, genomecov,
// vcf stats, vcf extract, vcf concat.
// ---------------------------------------------------------------------------

fn s_of(p: &Path) -> String {
    p.to_str().unwrap().to_string()
}

/// A text VCF without the `##glociCommand` line.
fn without_command(bytes: &[u8]) -> Vec<u8> {
    String::from_utf8_lossy(bytes)
        .lines()
        .filter(|l| !l.starts_with("##glociCommand"))
        .flat_map(|l| [l, "\n"])
        .collect::<String>()
        .into_bytes()
}

const P2_VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n\
##contig=<ID=chr2,length=10000>\n##contig=<ID=chr3,length=10000>\n##contig=<ID=chr4,length=10000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n##FILTER=<ID=q10,Description=\"low\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n";

/// Sorted multi-chromosome VCF with three samples and a mix of SNPs and indels.
fn p2_vcf(chroms: std::ops::RangeInclusive<usize>, per_chrom: usize) -> String {
    let mut t = String::from(P2_VCF_HEAD);
    for c in chroms {
        for i in 0..per_chrom {
            let (r, a) = match i % 5 {
                3 => ("AT", "A"),
                4 => ("A", "ATT"),
                k => (["A", "C", "G"][k % 3], ["C", "G", "T"][k % 3]),
            };
            let gt = ["0/1", "1/1", "0/0", "./.", "1|0"];
            writeln!(
                t,
                "chr{c}\t{}\trs{c}_{i}\t{r}\t{a}\t{}\t{}\tDP={}\tGT:DP\t{}:{}\t{}:{}\t{}:{}",
                100 + i * 37,
                20 + i % 60,
                if i % 7 == 0 { "q10" } else { "PASS" },
                i % 40,
                gt[i % 5],
                i % 11,
                gt[(i + c) % 5],
                (i + 3) % 17,
                gt[(i + 2 * c) % 5],
                (i + 5) % 23
            )
            .unwrap();
        }
    }
    t
}

/// Every P2 fixture, as plain / BGZF / indexed BGZF where the command takes it.
struct P2 {
    dir: tempfile::TempDir,
}

impl P2 {
    fn new() -> P2 {
        P2 { dir: tempdir() }
    }
    fn plain(&self, n: &str, t: &str) -> String {
        s_of(&plain_file(self.dir.path(), n, t))
    }
    fn bgzf(&self, n: &str, t: &str) -> String {
        s_of(&bgzf_file(self.dir.path(), n, t))
    }
    fn indexed(&self, n: &str, t: &str) -> String {
        let p = self.bgzf(n, t);
        gloci().args(["index", &p]).assert().success();
        p
    }
    fn bcf(&self, n: &str, t: &str) -> String {
        s_of(&write(self.dir.path(), n, &vcf_to_bcf(t.as_bytes())))
    }
}

fn bed_a(chroms: usize) -> String {
    let mut t = String::new();
    for c in 1..=chroms {
        for i in 0..60 {
            let st = 100 + i * 40 + (i % 3) * 7;
            writeln!(
                t,
                "chr{c}\t{st}\t{}\ta_{c}_{i}\t{}\t{}",
                st + 20 + i % 45,
                10 + i % 90,
                if i % 2 == 0 { '+' } else { '-' }
            )
            .unwrap();
        }
    }
    t
}

fn bed_b(chroms: usize) -> String {
    let mut t = String::new();
    for c in 1..=chroms {
        for i in 0..40 {
            let st = 130 + i * 55;
            writeln!(
                t,
                "chr{c}\t{st}\t{}\tb_{c}_{i}\t{}\t+",
                st + 30 + i % 25,
                i % 50
            )
            .unwrap();
        }
    }
    t
}

fn run_stdout(args: &[&str]) -> Vec<u8> {
    let o = gloci().args(args).output().unwrap();
    assert_eq!(o.status.code(), Some(0), "{args:?}: {:?}", o.stderr);
    o.stdout
}

/// Run `args` plus `-t T` for T in 1, 2, 8 (three spellings) and require
/// byte-identical output. This does not exercise chunk boundaries:
/// `GLOCI_TEST_CHUNK_MEM_LIMIT` only shrinks the chunks of the relation
/// commands' chunk driver (intersect, subtract, closest, map, coverage),
/// which the P2 commands do not use.
fn same_at_every_thread_count(what: &str, args: &[&str]) -> Vec<u8> {
    let outs: Vec<Vec<u8>> = ["1", "2", "8"]
        .iter()
        .map(|t| stdout_with(args, t))
        .collect();
    assert_eq!(outs[0], outs[1], "{what}: -t 1 vs -t 2");
    assert_eq!(outs[1], outs[2], "{what}: -t 2 vs -t 8");
    outs[0].clone()
}

/// The BED-shaped commands of the determinism matrices; `{a}` and `{b}` are
/// the inputs.
fn p2_bed_cmds(g: &str) -> Vec<(&str, Vec<&str>)> {
    vec![
        ("shift", vec!["shift", "{a}", "--by", "7", "-g", g]),
        (
            "shift -s --pct",
            vec!["shift", "{a}", "--by", "-0.25", "--pct", "-s"],
        ),
        ("cluster", vec!["cluster", "{a}"]),
        ("cluster -d -s", vec!["cluster", "{a}", "-d", "30", "-s"]),
        ("spacing", vec!["spacing", "{a}"]),
        (
            "groupby",
            vec!["groupby", "{a}", "-g", "1", "-c", "5:sum,5:max,4:collapse"],
        ),
        (
            "partition",
            vec!["partition", "{a}", "{b}", "--names", "A,B"],
        ),
        (
            "partition --all",
            vec![
                "partition",
                "{a}",
                "{b}",
                "--names",
                "A,B",
                "--all",
                "--merge-adjacent",
            ],
        ),
        ("genomecov", vec!["genomecov", "{a}", "-g", g]),
        (
            "genomecov --bedgraph",
            vec!["genomecov", "{a}", "-g", g, "--bedgraph"],
        ),
        (
            "genomecov --bga",
            vec!["genomecov", "{a}", "-g", g, "--bga"],
        ),
        (
            "genomecov --per-base",
            vec!["genomecov", "{a}", "-g", g, "--per-base"],
        ),
        (
            "genomecov --split",
            vec!["genomecov", "{a}", "-g", g, "--bedgraph", "--split"],
        ),
    ]
}

#[test]
fn p2_determinism_matrix_streams() {
    let fx = P2::new();
    let genome = fx.plain(
        "g.genome",
        "chr1\t10000\nchr2\t10000\nchr3\t10000\nchr4\t10000\n",
    );
    let (a, b) = (bed_a(4), bed_b(4));
    let a_plain = fx.plain("a.bed", &a);
    let a_gz = fx.bgzf("a.bed.gz", &a);
    let a_idx = fx.indexed("a_idx.bed.gz", &a);
    let b_plain = fx.plain("b.bed", &b);
    let b_gz = fx.bgzf("b.bed.gz", &b);
    // unsorted, for --sort
    let mut rev: Vec<&str> = a.lines().collect();
    rev.reverse();
    let a_unsorted = fx.plain("ua.bed", &(rev.join("\n") + "\n"));
    let bg = {
        let mut t = String::new();
        for c in 1..=4 {
            for i in 0..30 {
                writeln!(t, "chr{c}\t{}\t{}\t{}", i * 100, i * 100 + 60, 1 + i % 5).unwrap();
            }
        }
        t
    };
    let bg1 = fx.plain("one.bedgraph", &bg);
    let bg2 = fx.plain("two.bedgraph", &bg.replace("\t60\t", "\t80\t"));
    let vcf = p2_vcf(1..=4, 80);
    let vcf_plain = fx.plain("v.vcf", &vcf);
    let vcf_gz = fx.bgzf("v.vcf.gz", &vcf);
    let vcf_idx = fx.indexed("v_idx.vcf.gz", &vcf);
    let vcf_bcf = fx.bcf("v.bcf", &vcf);
    let g = genome.as_str();

    // BED-shaped input: one command per row, every format spelling identical.
    let bed_cmds = p2_bed_cmds(g);
    for (name, tpl) in &bed_cmds {
        let mut per_format = Vec::new();
        for (ai, bi) in [(&a_plain, &b_plain), (&a_gz, &b_gz), (&a_idx, &b_gz)] {
            let args: Vec<&str> = tpl
                .iter()
                .map(|x| match *x {
                    "{a}" => ai.as_str(),
                    "{b}" => bi.as_str(),
                    o => o,
                })
                .collect();
            per_format.push(same_at_every_thread_count(name, &args));
        }
        assert!(!per_format[0].is_empty(), "{name}: empty output");
        assert_eq!(per_format[0], per_format[1], "{name}: plain vs bgzf");
        assert_eq!(
            per_format[0], per_format[2],
            "{name}: plain vs indexed bgzf"
        );
    }

    // --sort on unsorted input is a different engine path.
    for (name, tpl) in [
        ("cluster --sort", vec!["cluster", &a_unsorted, "--sort"]),
        ("spacing --sort", vec!["spacing", &a_unsorted, "--sort"]),
        (
            "partition --sort",
            vec![
                "partition",
                &a_unsorted,
                &b_plain,
                "--sort",
                "--names",
                "A,B",
            ],
        ),
        (
            "genomecov --sort",
            vec!["genomecov", &a_unsorted, "-g", g, "--sort", "--bedgraph"],
        ),
    ] {
        let out = same_at_every_thread_count(name, &tpl);
        assert!(!out.is_empty(), "{name}");
    }
    // bedGraph partition
    let out = same_at_every_thread_count(
        "partition --bedgraph",
        &["partition", "--bedgraph", &bg1, &bg2],
    );
    assert!(!out.is_empty());

    // VCF input: identical across plain, bgzf, indexed and BCF.
    let vcf_cmds: Vec<(&str, Vec<&str>)> = vec![
        ("cluster", vec!["cluster", "{v}"]),
        ("spacing", vec!["spacing", "{v}"]),
        (
            "partition",
            vec!["partition", "{v}", &b_plain, "--names", "A,B"],
        ),
        ("genomecov", vec!["genomecov", "{v}", "-g", g, "--bedgraph"]),
        (
            "extract --info",
            vec!["vcf", "extract", "{v}", "--info", "DP"],
        ),
        (
            "extract --format wide",
            vec!["vcf", "extract", "{v}", "--format", "GT,DP"],
        ),
        (
            "extract --format long",
            vec![
                "vcf", "extract", "{v}", "--format", "GT,DP", "--layout", "long",
            ],
        ),
        (
            "extract filtered",
            vec![
                "vcf",
                "extract",
                "{v}",
                "--info",
                "DP",
                "--min-qual",
                "30",
                "--biallelic",
            ],
        ),
        ("concat", vec!["vcf", "concat", "{v}", "{v2}"]),
    ];
    for cmd in ["cluster", "spacing"] {
        gloci()
            .args([cmd, &vcf_bcf])
            .assert()
            .code(2)
            .stderr(contains("convert first"));
    }
    // concat needs the second file to start where the first ends
    let second = fx.plain(
        "v_tail.vcf",
        &p2_vcf(4..=4, 1).replace("chr4\t100", "chr4\t99999"),
    );
    for (name, tpl) in &vcf_cmds {
        let mut outs = Vec::new();
        // cluster and spacing append a column, which a BCF record has no place for
        let takes_bcf = !matches!(*name, "cluster" | "spacing");
        for v in [&vcf_plain, &vcf_gz, &vcf_idx, &vcf_bcf] {
            if !takes_bcf && *v == vcf_bcf {
                continue;
            }
            let args: Vec<&str> = tpl
                .iter()
                .map(|x| match *x {
                    "{v}" => v.as_str(),
                    "{v2}" => second.as_str(),
                    o => o,
                })
                .collect();
            // concat's header differs between plain and BCF in nothing but the
            // glociCommand line
            outs.push(without_command(&same_at_every_thread_count(name, &args)));
        }
        assert!(!outs[0].is_empty(), "{name}");
        // concat's default output format follows its first input (see
        // p2_concat_output_format), so only its thread determinism is
        // compared here
        if *name != "concat" {
            for w in outs.windows(2) {
                assert_eq!(w[0], w[1], "{name}: input format changed the output");
            }
        }
    }
}

/// `vcf stats`, `vcf extract --genotype-matrix` and `vcf concat` write files.
#[test]
fn p2_determinism_matrix_files() {
    let fx = P2::new();
    let vcf = p2_vcf(1..=4, 120);
    let plain = fx.plain("v.vcf", &vcf);
    let gz = fx.bgzf("v.vcf.gz", &vcf);
    let bcf = fx.bcf("v.bcf", &vcf);
    let flags = [
        "--freq",
        "--counts",
        "--depth",
        "--site-depth",
        "--site-mean-depth",
        "--geno-depth",
        "--site-quality",
        "--missing-site",
        "--missing-indv",
        "--het",
        "--hardy",
        "--singletons",
        "--filter-summary",
        "--tstv-summary",
        "--indel-hist",
        "--summary",
    ];
    let mut reference: Option<Vec<(String, Vec<u8>)>> = None;
    for input in [&plain, &gz, &bcf] {
        for t in ["1", "2", "8"] {
            let d = tempdir();
            let prefix = s_of(&d.path().join("o"));
            let mut args = vec![
                "vcf",
                "stats",
                input.as_str(),
                "--out-prefix",
                &prefix,
                "-t",
                t,
            ];
            args.extend(flags);
            args.extend(["--snp-density", "500"]);
            gloci().args(&args).assert().success();
            let mut files: Vec<(String, Vec<u8>)> = std::fs::read_dir(d.path())
                .unwrap()
                .map(|e| {
                    let e = e.unwrap();
                    (
                        e.file_name().to_string_lossy().into_owned(),
                        std::fs::read(e.path()).unwrap(),
                    )
                })
                .collect();
            files.sort();
            assert!(
                files.len() >= 16,
                "{input} -t {t}: wrote {} files",
                files.len()
            );
            assert!(files.iter().all(|(_, b)| !b.is_empty()));
            match &reference {
                None => reference = Some(files),
                Some(r) => assert_eq!(&files, r, "{input} -t {t}"),
            }
        }
    }
    // --json is a single document and is deterministic too
    let outs: Vec<Vec<u8>> = ["1", "2", "8"]
        .iter()
        .map(|t| {
            let d = tempdir();
            let prefix = s_of(&d.path().join("o"));
            gloci()
                .args([
                    "vcf",
                    "stats",
                    &plain,
                    "--out-prefix",
                    &prefix,
                    "--summary",
                    "--json",
                    "-t",
                    t,
                ])
                .assert()
                .success();
            let mut names: Vec<_> = std::fs::read_dir(d.path())
                .unwrap()
                .map(|e| e.unwrap().path())
                .collect();
            names.sort();
            names
                .iter()
                .flat_map(|p| std::fs::read(p).unwrap())
                .collect()
        })
        .collect();
    assert_eq!(outs[0], outs[1]);
    assert_eq!(outs[1], outs[2]);

    // genotype matrix: three files, identical across thread counts and inputs
    let mut reference: Option<Vec<Vec<u8>>> = None;
    for input in [&plain, &gz, &bcf] {
        for t in ["1", "2", "8"] {
            let d = tempdir();
            let prefix = s_of(&d.path().join("m"));
            gloci()
                .args([
                    "vcf",
                    "extract",
                    input,
                    "--genotype-matrix",
                    "012",
                    "--out-prefix",
                    &prefix,
                    "-t",
                    t,
                ])
                .assert()
                .success();
            let files: Vec<Vec<u8>> = ["012", "012.indv", "012.pos"]
                .iter()
                .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap())
                .collect();
            assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 3);
            assert!(files.iter().all(|f| !f.is_empty()));
            match &reference {
                None => reference = Some(files),
                Some(r) => assert_eq!(&files, r, "{input} -t {t}"),
            }
        }
    }

    // extract to a file; concat checked and --naive to files
    let part1 = fx.bgzf("p1.vcf.gz", &p2_vcf(1..=2, 120));
    let part2 = fx.bgzf("p2.vcf.gz", &p2_vcf(3..=4, 120));
    for (name, extra) in [("checked", vec![]), ("naive", vec!["--naive"])] {
        for ext in ["vcf.gz", "bcf"] {
            if name == "naive" && ext == "bcf" {
                continue;
            }
            let mut outs = Vec::new();
            for t in ["1", "2", "8"] {
                let d = tempdir();
                let out = s_of(&d.path().join(format!("o.{ext}")));
                let mut args = vec!["vcf", "concat", &part1, &part2, "-o", &out, "-t", t];
                args.extend(extra.iter().copied());
                gloci().args(&args).assert().success();
                outs.push(std::fs::read(&out).unwrap());
            }
            assert_eq!(outs[0], outs[1], "concat {name} {ext}: -t 1 vs 2");
            assert_eq!(outs[1], outs[2], "concat {name} {ext}: -t 2 vs 8");
        }
    }
    let mut outs = Vec::new();
    for t in ["1", "2", "8"] {
        let d = tempdir();
        let out = s_of(&d.path().join("o.tsv"));
        gloci()
            .args([
                "vcf", "extract", &gz, "--format", "GT,DP", "-o", &out, "-t", t,
            ])
            .assert()
            .success();
        outs.push(std::fs::read(&out).unwrap());
    }
    assert_eq!(outs[0], outs[1]);
    assert_eq!(outs[1], outs[2]);
}

#[test]
fn p2_broken_pipe_matrix() {
    let d = tempdir();
    let genome = plain_file(d.path(), "g.genome", "chr1\t1000000\n");
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..20_000 {
        let s = i * 20;
        writeln!(a, "chr1\t{s}\t{}\ta{i}\t100\t+", s + 15).unwrap();
        writeln!(b, "chr1\t{}\t{}\tb{i}\t50\t+", s + 5, s + 25).unwrap();
    }
    let mut v = String::from(P2_VCF_HEAD);
    for i in 0..20_000 {
        writeln!(
            v,
            "chr1\t{}\t.\tA\tC\t50\tPASS\tDP=30\tGT:DP\t0/1:5\t1/1:6\t0/0:7",
            1 + i * 20
        )
        .unwrap();
    }
    let a = s_of(&plain_file(d.path(), "a.bed", &a));
    let b = s_of(&plain_file(d.path(), "b.bed", &b));
    let v = s_of(&plain_file(d.path(), "v.vcf", &v));
    let tail = tail_after_big_vcf(d.path());
    let g = s_of(&genome);

    let cmds: Vec<Vec<&str>> = vec![
        vec!["shift", &a, "--by", "3"],
        vec!["cluster", &a],
        vec!["spacing", &a],
        vec!["groupby", &a, "-g", "4", "-c", "5:sum"],
        vec!["partition", &a, &b],
        vec!["partition", &a, &b, "--no-membership"],
        vec!["genomecov", &a, "-g", &g, "--per-base"],
        vec!["genomecov", &a, "-g", &g, "--bga"],
        vec!["vcf", "extract", &v, "--info", "DP"],
        vec![
            "vcf", "extract", &v, "--format", "GT,DP", "--layout", "long",
        ],
        vec!["vcf", "concat", &v, &tail],
        vec!["cluster", &v],
        vec!["spacing", &v],
    ];
    for args in &cmds {
        for t in ["1", "4"] {
            let mut child = Command::new(bin())
                .args(args)
                .args(["-t", t])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            let mut first = [0u8; 64];
            let _ = child.stdout.as_mut().unwrap().read(&mut first);
            drop(child.stdout.take());
            let out = child.wait_with_output().unwrap();
            assert_eq!(
                out.status.code(),
                Some(0),
                "{args:?} -t {t}: stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert!(
                out.stderr.is_empty(),
                "{args:?} -t {t} printed on stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
    }
}

/// (Writes a plain `.vcf`.) A one-record VCF that follows the 20 000-record one at the same
/// chromosome end (for `vcf concat`).
fn tail_after_big_vcf(dir: &Path) -> String {
    let mut t = String::from(P2_VCF_HEAD);
    t.push_str("chr1\t900000\t.\tA\tC\t50\tPASS\tDP=30\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n");
    s_of(&plain_file(dir, "tail.vcf", &t))
}

#[test]
fn p2_empty_inputs_matrix() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t1000\nchr2\t1000\n");
    let a = fx.plain(
        "a.bed",
        "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n",
    );
    let empty = fx.plain("empty.bed", "");
    let empty_gz = fx.bgzf("empty.bed.gz", "");
    let hdr_only = fx.plain("hdr.vcf", &p2_vcf(1..=1, 0));
    let hdr_only_gz = fx.bgzf("hdr.vcf.gz", &p2_vcf(1..=1, 0));
    let hdr_only_bcf = fx.bcf("hdr.bcf", &p2_vcf(1..=1, 0));
    let one = fx.plain("one.vcf", &p2_vcf(1..=1, 3));
    let (g, a) = (genome.as_str(), a.as_str());

    // empty BED (zero bytes, and an empty BGZF stream): no records out
    for e in [&empty, &empty_gz] {
        let silent: Vec<Vec<&str>> = vec![
            vec!["shift", e, "--by", "5"],
            vec!["cluster", e],
            vec!["spacing", e],
            vec!["groupby", e, "-g", "1", "-c", "2:sum"],
            vec!["partition", e, "--names", "A", "--min-files", "1"],
            vec!["genomecov", e, "-g", g, "--bedgraph"],
        ];
        for args in &silent {
            for t in ["1", "8"] {
                assert!(stdout_with(args, t).is_empty(), "{args:?} -t {t}");
            }
        }
        // genomecov: an empty input still covers every chromosome with zero
        assert_eq!(
            String::from_utf8(run_stdout(&["genomecov", e, "-g", g, "--bga"])).unwrap(),
            "chr1\t0\t1000\t0\nchr2\t0\t1000\t0\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&["genomecov", e, "-g", g, "--hist"])).unwrap(),
            "chr1\t0\t1000\t1000\t1\nchr2\t0\t1000\t1000\t1\ngenome\t0\t2000\t2000\t1\n"
        );
        assert!(run_stdout(&["genomecov", e, "-g", g, "--per-base"]).starts_with(b"chr1\t1\t0\n"));
        // one empty and one non-empty input to partition
        assert_eq!(
            String::from_utf8(run_stdout(&["partition", e, a, "--names", "E,A"])).unwrap(),
            "chr1\t10\t50\t1\tA\nchr1\t60\t100\t1\tA\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&["partition", a, e, "--names", "A,E"])).unwrap(),
            "chr1\t10\t50\t1\tA\nchr1\t60\t100\t1\tA\n"
        );
        assert!(run_stdout(&["partition", e, e, "--names", "X,Y"]).is_empty());
        // -o with empty input: the file exists and is empty
        let d = tempdir();
        let out = s_of(&d.path().join("o.bed"));
        gloci().args(["cluster", e, "-o", &out]).assert().success();
        assert_eq!(std::fs::metadata(&out).unwrap().len(), 0);
    }

    // header-only VCF: cluster/spacing keep the column header, nothing else
    for v in [&hdr_only, &hdr_only_gz] {
        for cmd in ["cluster", "spacing"] {
            let out = String::from_utf8(run_stdout(&[cmd, v])).unwrap();
            assert_eq!(out.lines().count(), 1, "{cmd}: {out}");
            assert!(out.starts_with("#CHROM\tPOS"), "{out}");
        }
        assert!(run_stdout(&["partition", v, "--names", "V"]).is_empty());
        assert_eq!(
            String::from_utf8(run_stdout(&["vcf", "extract", v, "--info", "DP"])).unwrap(),
            "CHROM\tPOS\tID\tREF\tALT\tDP\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&["vcf", "extract", v, "--format", "GT"])).unwrap(),
            "CHROM\tPOS\tID\tREF\tALT\ts1:GT\ts2:GT\ts3:GT\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&[
                "vcf", "extract", v, "--format", "GT", "--layout", "long"
            ]))
            .unwrap()
            .lines()
            .count(),
            1
        );
    }
    assert_eq!(
        String::from_utf8(run_stdout(&[
            "vcf",
            "extract",
            &hdr_only_bcf,
            "--info",
            "DP"
        ]))
        .unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\tDP\n"
    );
    assert_eq!(
        String::from_utf8(run_stdout(&[
            "vcf",
            "extract",
            &hdr_only_bcf,
            "--format",
            "GT"
        ]))
        .unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\ts1:GT\ts2:GT\ts3:GT\n"
    );

    // vcf stats on a header-only VCF: every file is written, with its header
    let d = tempdir();
    let prefix = s_of(&d.path().join("o"));
    gloci()
        .args([
            "vcf",
            "stats",
            &hdr_only,
            "--out-prefix",
            &prefix,
            "--freq",
            "--counts",
            "--depth",
            "--site-depth",
            "--site-mean-depth",
            "--geno-depth",
            "--site-quality",
            "--missing-site",
            "--missing-indv",
            "--het",
            "--hardy",
            "--singletons",
            "--filter-summary",
            "--tstv-summary",
            "--indel-hist",
            "--summary",
        ])
        .assert()
        .success();
    let written: Vec<String> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert!(written.len() >= 14, "{written:?}");
    let imiss = std::fs::read_to_string(format!("{prefix}.imiss")).unwrap();
    assert_eq!(
        imiss.lines().count(),
        4,
        "one line per sample plus header: {imiss}"
    );
    let frq = std::fs::read_to_string(format!("{prefix}.frq")).unwrap();
    assert_eq!(frq.lines().count(), 1, "{frq}");

    // genotype matrix of a header-only VCF: empty matrix, samples listed
    let prefix = s_of(&d.path().join("m"));
    gloci()
        .args([
            "vcf",
            "extract",
            &hdr_only,
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read(format!("{prefix}.012")).unwrap(),
        b"0\n1\n2\n"
    );
    assert_eq!(std::fs::read(format!("{prefix}.012.pos")).unwrap(), b"");
    assert_eq!(
        std::fs::read_to_string(format!("{prefix}.012.indv")).unwrap(),
        "s1\ns2\ns3\n"
    );

    // concat: empty files are skipped without changing the result
    let want = String::from_utf8(without_command(&run_stdout(&["vcf", "concat", &one]))).unwrap();
    for args in [
        vec!["vcf", "concat", &hdr_only, &one],
        vec!["vcf", "concat", &one, &hdr_only],
        vec!["vcf", "concat", &hdr_only, &one, &hdr_only],
    ] {
        let got = without_command(&run_stdout(&args));
        assert_eq!(String::from_utf8(got).unwrap(), want, "{args:?}");
    }
    let both = String::from_utf8(run_stdout(&["vcf", "concat", &hdr_only, &hdr_only])).unwrap();
    assert!(
        both.contains("#CHROM") && both.lines().all(|l| l.starts_with('#')),
        "{both}"
    );
    let naive_empty = fx.bgzf("hdr_n.vcf.gz", &p2_vcf(1..=1, 0));
    let d2 = tempdir();
    let out = s_of(&d2.path().join("n.vcf.gz"));
    gloci()
        .args([
            "vcf",
            "concat",
            "--naive",
            &naive_empty,
            &hdr_only_gz,
            "-o",
            &out,
        ])
        .assert()
        .success();
}

#[test]
fn p2_vcf_input_matrix() {
    let fx = P2::new();
    // an indel and a symbolic deletion: the footprints differ by --vcf-span
    let vcf = format!(
        "{P2_VCF_HEAD}chr1\t100\t.\tA\tC\t50\tPASS\tDP=3\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n\
chr1\t103\t.\tATGC\tA\t50\tPASS\tDP=3\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n\
chr1\t500\t.\tA\t<DEL>\t50\tPASS\tDP=3;END=700\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n"
    );
    let v = fx.plain("v.vcf", &vcf);
    let g = fx.plain("g.genome", "chr1\t10000\n");
    // cluster appends a TSV column, keeps the record text and the header
    let out = String::from_utf8(run_stdout(&["cluster", &v])).unwrap();
    let lines: Vec<&str> = out.lines().collect();
    assert!(lines[0].starts_with("#CHROM") && lines[0].ends_with("\tcluster"));
    assert!(lines[1].ends_with("\t1") && lines[2].ends_with("\t2") && lines[3].ends_with("\t3"));
    assert!(!out.contains("##"), "TSV output carries no ## meta lines");
    // spacing: gap between footprints
    let out = String::from_utf8(run_stdout(&["spacing", &v])).unwrap();
    let last: Vec<&str> = out
        .lines()
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    assert_eq!(last, ["spacing", ".", "3", "394"], "{out}");
    // genomecov counts footprints: the span mode decides what a record covers
    let cov = |span: &str| {
        String::from_utf8(run_stdout(&[
            "genomecov",
            &v,
            "-g",
            &g,
            "--bedgraph",
            "--vcf-span",
            span,
        ]))
        .unwrap()
    };
    assert_eq!(
        cov("allele"),
        "chr1\t99\t100\t1\nchr1\t103\t106\t1\nchr1\t500\t700\t1\n"
    );
    assert_eq!(
        cov("ref"),
        "chr1\t99\t100\t1\nchr1\t102\t106\t1\nchr1\t499\t700\t1\n"
    );
    // partition labels the VCF by file name and lists breakpoints
    let out = String::from_utf8(run_stdout(&["partition", &v, "--names", "V"])).unwrap();
    assert!(out.lines().all(|l| l.ends_with("\tV")), "{out}");
    // --sort accepts an unsorted VCF
    let mut rows: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
    rows.reverse();
    let uv = fx.plain("u.vcf", &format!("{P2_VCF_HEAD}{}\n", rows.join("\n")));
    assert_eq!(
        run_stdout(&["cluster", &uv, "--sort"]),
        run_stdout(&["cluster", &v])
    );
    // shift and groupby do not take VCF, and say what to do
    gloci()
        .args(["shift", &v, "--by", "1"])
        .assert()
        .code(2)
        .stderr(contains("gloci view -O bed"));
    gloci()
        .args(["groupby", &v, "-g", "1", "-c", "2:count"])
        .assert()
        .code(2)
        .stderr(contains("gloci vcf extract"));
    // vcf stats and extract take --samples and filters on the VCF
    let d = tempdir();
    let prefix = s_of(&d.path().join("o"));
    gloci()
        .args([
            "vcf",
            "stats",
            &v,
            "--out-prefix",
            &prefix,
            "--missing-indv",
            "--samples",
            "s2",
        ])
        .assert()
        .success();
    let imiss = std::fs::read_to_string(format!("{prefix}.imiss")).unwrap();
    assert_eq!(imiss.lines().count(), 2, "{imiss}");
    let out = String::from_utf8(run_stdout(&[
        "vcf",
        "extract",
        &v,
        "--format",
        "GT",
        "--samples",
        "s3",
        "--biallelic",
    ]))
    .unwrap();
    assert_eq!(
        out.lines().next().unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\ts3:GT"
    );
    assert_eq!(out.lines().count(), 4, "{out}");
}

/// One failing run of a P2 command: which output it writes decides what
/// "nothing left behind" means.
enum Sink {
    /// `-o DIR/out.EXT`
    File(&'static str),
    /// `--out-prefix DIR/out`
    Prefix,
}

#[test]
fn p2_failure_matrix() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t10000\nchr2\t10000\n");
    let good_bed = fx.plain("good.bed", "chr1\t10\t50\ta\t1\t+\nchr1\t60\t90\tb\t2\t-\n");
    let other_bed = fx.plain("other.bed", "chr1\t20\t70\tc\t3\t+\n");
    let good_vcf_text = p2_vcf(1..=2, 4);
    let good_vcf = fx.plain("good.vcf", &good_vcf_text);
    let good_vcf_gz = fx.bgzf("good.vcf.gz", &good_vcf_text);

    let unsorted_bed = fx.plain(
        "u.bed",
        "chr1\t100\t200\tr1\t1\t+\nchr1\t50\t150\tr2\t2\t+\n",
    );
    let malformed_bed = fx.plain(
        "m.bed",
        "chr1\t100\t200\tr1\t1\t+\nchr1\tNOPE\t250\tr2\t2\t+\n",
    );
    let bed_gz = bgzip(bed_a(1).as_bytes());
    let trunc_bed = s_of(&write(
        fx.dir.path(),
        "t.bed.gz",
        &bed_gz[..bed_gz.len() / 2],
    ));
    let nofooter_bed = s_of(&write(
        fx.dir.path(),
        "nf.bed.gz",
        &bed_gz[..bed_gz.len() - 28],
    ));

    let rows: Vec<&str> = good_vcf_text.lines().collect();
    let hdr = rows.iter().filter(|l| l.starts_with('#')).count();
    let body: Vec<&str> = rows[hdr..].to_vec();
    let swapped = {
        let mut b = body.clone();
        b.swap(1, 2);
        fx.plain("u.vcf", &format!("{P2_VCF_HEAD}{}\n", b.join("\n")))
    };
    let malformed_vcf = fx.plain(
        "m.vcf",
        &format!(
            "{P2_VCF_HEAD}{}\n{}\n",
            body[0],
            body[1].replacen("\t137\t", "\tXX\t", 1)
        ),
    );
    let good_bcf = vcf_to_bcf(good_vcf_text.as_bytes());
    let trunc_bcf = s_of(&write(
        fx.dir.path(),
        "t.bcf",
        &good_bcf[..good_bcf.len() / 2],
    ));
    let vgz = bgzip(good_vcf_text.as_bytes());
    let trunc_vcf = s_of(&write(fx.dir.path(), "t.vcf.gz", &vgz[..vgz.len() / 2]));
    let nofooter_vcf = s_of(&write(fx.dir.path(), "nf.vcf.gz", &vgz[..vgz.len() - 28]));

    let after_gz = fx.bgzf("after.vcf.gz", &p2_vcf(3..=4, 4));
    let before_gz = fx.bgzf("before.vcf.gz", &p2_vcf(1..=1, 0));
    let g = genome.as_str();
    let (gb, ob, gv, gvz) = (&good_bed, &other_bed, &good_vcf, &good_vcf_gz);
    let unsorted_msg = "is not sorted";
    let eof = "BGZF end-of-file marker missing";

    // (name, args with {IN} marking the bad input, sink, per-kind messages)
    type Case<'a> = (&'a str, Vec<&'a str>, Sink);
    let cases: Vec<Case> = vec![
        (
            "shift",
            vec!["shift", "{IN}", "--by", "3"],
            Sink::File("bed"),
        ),
        ("cluster", vec!["cluster", "{IN}"], Sink::File("tsv")),
        ("spacing", vec!["spacing", "{IN}"], Sink::File("tsv")),
        (
            "groupby",
            vec!["groupby", "{IN}", "-g", "1", "-c", "2:sum"],
            Sink::File("tsv"),
        ),
        (
            "partition",
            vec!["partition", "{IN}", ob, "--names", "A,B"],
            Sink::File("bed"),
        ),
        (
            "partition B",
            vec!["partition", gb, "{IN}", "--names", "A,B"],
            Sink::File("bed"),
        ),
        (
            "genomecov",
            vec!["genomecov", "{IN}", "-g", g],
            Sink::File("tsv"),
        ),
    ];
    // BED-shaped failures: (input, needs sorted input, message)
    let bed_bad: [(&str, &String, &str); 4] = [
        ("unsorted", &unsorted_bed, unsorted_msg),
        ("malformed", &malformed_bed, "NOPE"),
        ("truncated", &trunc_bed, eof),
        ("no BGZF footer", &nofooter_bed, eof),
    ];
    let mut ran = 0;
    for (name, tpl, sink) in &cases {
        for (kind, bad, msg) in &bed_bad {
            // shift and groupby never check the order
            // shift and groupby never check the order, and must keep working
            if *kind == "unsorted" && matches!(*name, "shift" | "groupby") {
                assert_succeeds(name, tpl, bad, sink);
                continue;
            }
            // groupby sums column 2, so the malformed row is a bad number there
            let want = if *name == "groupby" && *kind == "malformed" {
                "is not a number"
            } else {
                msg
            };
            assert_fails_clean(name, kind, tpl, bad, sink, want);
            ran += 1;
        }
        // missing output directory, on good input
        assert_missing_dir(name, tpl, gb, sink);
    }

    // VCF-shaped failures
    let vcf_cases: Vec<Case> = vec![
        ("cluster", vec!["cluster", "{IN}"], Sink::File("tsv")),
        ("spacing", vec!["spacing", "{IN}"], Sink::File("tsv")),
        (
            "partition",
            vec!["partition", "{IN}", ob, "--names", "A,B"],
            Sink::File("bed"),
        ),
        (
            "genomecov",
            vec!["genomecov", "{IN}", "-g", g],
            Sink::File("tsv"),
        ),
        (
            "extract --info",
            vec!["vcf", "extract", "{IN}", "--info", "DP"],
            Sink::File("tsv"),
        ),
        (
            "extract --format",
            vec!["vcf", "extract", "{IN}", "--format", "GT"],
            Sink::File("tsv"),
        ),
        (
            "extract --012",
            vec!["vcf", "extract", "{IN}", "--genotype-matrix", "012"],
            Sink::Prefix,
        ),
        (
            "stats --freq",
            vec!["vcf", "stats", "{IN}", "--freq"],
            Sink::Prefix,
        ),
        (
            "stats --summary",
            vec![
                "vcf",
                "stats",
                "{IN}",
                "--summary",
                "--het",
                "--missing-indv",
            ],
            Sink::Prefix,
        ),
        (
            "concat first",
            vec!["vcf", "concat", "{IN}", &after_gz],
            Sink::File("vcf.gz"),
        ),
        (
            "concat second",
            vec!["vcf", "concat", &before_gz, "{IN}"],
            Sink::File("vcf.gz"),
        ),
    ];
    // extract and stats read an unsorted VCF without complaint (extract
    // writes in input order; stats accumulates), so they must succeed on it:
    // adding an order check to either would be a visible decision here
    let vcf_bad: [(&str, &String, &str); 4] = [
        ("unsorted", &swapped, unsorted_msg),
        ("malformed", &malformed_vcf, "XX"),
        ("truncated", &trunc_vcf, eof),
        ("no BGZF footer", &nofooter_vcf, eof),
    ];
    for (name, tpl, sink) in &vcf_cases {
        for (kind, bad, msg) in &vcf_bad {
            let order_free = name.starts_with("extract") || name.starts_with("stats");
            if *kind == "unsorted" && order_free {
                assert_succeeds(name, tpl, bad, sink);
                continue;
            }
            assert_fails_clean(name, kind, tpl, bad, sink, msg);
            ran += 1;
        }
        assert_missing_dir(name, tpl, gv, sink);
    }
    // truncated BCF, for the commands that take BCF
    for (name, tpl, sink) in vcf_cases.iter().filter(|(n, ..)| {
        n.starts_with("extract") || n.starts_with("stats") || n.starts_with("concat")
    }) {
        assert_fails_clean(name, "truncated BCF", tpl, &trunc_bcf, sink, eof);
        ran += 1;
    }
    // a BED where a VCF is required says so
    gloci()
        .args([
            "vcf",
            "stats",
            gb,
            "--freq",
            "--out-prefix",
            &s_of(&fx.dir.path().join("x")),
        ])
        .assert()
        .failure()
        .stderr(contains("BED"));
    // concat: the second file starting before the first ends, and mismatched
    // contigs on --naive input
    let d = tempdir();
    let out = s_of(&d.path().join("o.vcf.gz"));
    gloci()
        .args(["vcf", "concat", gvz, gvz, "-o", &out])
        .assert()
        .code(1)
        .stderr(contains("each must start where the one before ends"));
    assert!(!Path::new(&out).exists());
    gloci()
        .args(["vcf", "concat", "--naive", gvz, gvz, "-o", &out])
        .assert()
        .failure()
        .stderr(contains("is not sorted"));
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
    gloci()
        .args(["vcf", "concat", "--naive", gv, gvz, "-o", &out])
        .assert()
        .code(2)
        .stderr(contains("compress it with: gloci view"));
    assert!(ran > 60, "ran {ran} failure cases");
}

/// Run `tpl` with `{IN}` replaced by the bad input, once per thread count:
/// it exits 1, says `msg`, and writes nothing.
fn assert_fails_clean(name: &str, kind: &str, tpl: &[&str], bad: &str, sink: &Sink, msg: &str) {
    for t in ["1", "4"] {
        let d = tempdir();
        let out_dir = d.path().join("out");
        std::fs::create_dir(&out_dir).unwrap();
        let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", bad)).collect();
        match sink {
            Sink::File(ext) => {
                args.extend(["-o".into(), s_of(&out_dir.join(format!("out.{ext}")))])
            }
            Sink::Prefix => args.extend(["--out-prefix".into(), s_of(&out_dir.join("out"))]),
        }
        args.extend(["-t".into(), t.into()]);
        let o = gloci().args(&args).output().unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert_eq!(o.status.code(), Some(1), "{name} / {kind} -t {t}: {err}");
        assert!(
            err.contains(msg),
            "{name} / {kind} -t {t}: wanted `{msg}` in: {err}"
        );
        assert!(err.starts_with("error: "), "{name} / {kind} -t {t}: {err}");
        let left: Vec<_> = std::fs::read_dir(&out_dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert!(
            left.is_empty(),
            "{name} / {kind} -t {t}: left behind {left:?}"
        );
    }
}

/// With the output directory missing the run exits 1, names the path and
/// how to fix it, and creates nothing.
fn assert_missing_dir(name: &str, tpl: &[&str], good: &str, sink: &Sink) {
    let d = tempdir();
    let missing = d.path().join("no/such/dir");
    let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", good)).collect();
    match sink {
        Sink::File(ext) => args.extend(["-o".into(), s_of(&missing.join(format!("out.{ext}")))]),
        Sink::Prefix => args.extend(["--out-prefix".into(), s_of(&missing.join("out"))]),
    }
    let o = gloci().args(&args).output().unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert_eq!(o.status.code(), Some(1), "{name}: {args:?}: {err}");
    assert!(
        err.contains("no/such/dir")
            && err.contains("No such file or directory")
            && err.contains("check that the directory exists"),
        "{name}: {err}"
    );
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0, "{name}");
}

/// `--chrom-map` is a rename applied as records are read: running with the
/// map on a file in the source names equals running without it on the file
/// already in the target names, and the inverse map takes the output back.
#[test]
fn p2_chrom_map_round_trip() {
    let fx = P2::new();
    let fwd = fx.plain("fwd.map", "1\tchr1\n2\tchr2\n");
    let inv = fx.plain("inv.map", "chr1\t1\nchr2\t2\n");
    let genome_chr = fx.plain("g_chr.genome", "chr1\t10000\nchr2\t10000\n");
    let chr_bed = bed_a(2);
    let num_bed = chr_bed.replace("chr", "");
    let bed_chr = fx.plain("chr.bed", &chr_bed);
    let bed_num = fx.plain("num.bed", &num_bed);
    let other_chr = fx.plain("o_chr.bed", &bed_b(2));
    let other_num = fx.plain("o_num.bed", &bed_b(2).replace("chr", ""));
    let vcf_chr = p2_vcf(1..=2, 30);
    let vcf_num = vcf_chr.replace("chr", "");
    let vchr = fx.plain("chr.vcf", &vcf_chr);
    let vnum = fx.plain("num.vcf", &vcf_num);
    let vnum_gz = fx.bgzf("num.vcf.gz", &vcf_num);
    let vnum_bcf = fx.bcf("num.bcf", &vcf_num);

    // commands that write input lines back: the output round-trips exactly
    for cmd in ["cluster", "spacing"] {
        for t in ["1", "8"] {
            let renamed = stdout_with(&[cmd, &bed_num, "--chrom-map", &fwd], t);
            assert_eq!(
                renamed,
                stdout_with(&[cmd, &bed_chr], t),
                "{cmd} map vs renamed"
            );
            let back = stdout_with(&[cmd, &bed_chr, "--chrom-map", &inv], t);
            assert_eq!(
                back,
                stdout_with(&[cmd, &bed_num], t),
                "{cmd}: the inverse map restores the source names"
            );
        }
    }
    // the others
    let g = genome_chr.as_str();
    let cases: Vec<Vec<&str>> = vec![
        vec!["genomecov", "{I}", "-g", g, "--bedgraph"],
        vec!["genomecov", "{I}", "-g", g],
        vec!["partition", "{I}", "{O}", "--names", "A,B"],
    ];
    for tpl in &cases {
        let build = |i: &str, o: &str, map: bool| -> Vec<String> {
            let mut v: Vec<String> = tpl
                .iter()
                .map(|a| match *a {
                    "{I}" => i.to_string(),
                    "{O}" => o.to_string(),
                    x => x.to_string(),
                })
                .collect();
            if map {
                v.extend(["--chrom-map".into(), fwd.clone()]);
            }
            v
        };
        let run = |args: Vec<String>| {
            let refs: Vec<&str> = args.iter().map(String::as_str).collect();
            same_at_every_thread_count("chrom-map", &refs)
        };
        let with_map = run(build(&bed_num, &other_num, true));
        let renamed = run(build(&bed_chr, &other_chr, false));
        assert!(!with_map.is_empty());
        assert_eq!(with_map, renamed, "{tpl:?}");
        assert!(
            !String::from_utf8_lossy(&with_map)
                .lines()
                .any(|l| l.starts_with("1\t"))
        );
    }
    // a name the map does not list passes through; mixed names work
    let mixed = fx.plain("mixed.bed", "1\t10\t20\nchr2\t5\t9\nchrX\t1\t2\n");
    let out = String::from_utf8(run_stdout(&[
        "cluster",
        &mixed,
        "--chrom-map",
        &fwd,
        "--sort",
    ]))
    .unwrap();
    assert_eq!(out, "chr1\t10\t20\t1\nchr2\t5\t9\t2\nchrX\t1\t2\t3\n");

    // VCF input (plain, BGZF and BCF carry the source names in the header):
    // extract and stats write the target names
    for input in [&vnum, &vnum_gz, &vnum_bcf] {
        let mapped = stdout_with(
            &["vcf", "extract", input, "--info", "DP", "--chrom-map", &fwd],
            "2",
        );
        let renamed = stdout_with(&["vcf", "extract", &vchr, "--info", "DP"], "2");
        assert_eq!(mapped, renamed, "extract {input}");
        assert!(mapped.starts_with(b"CHROM\tPOS\tID\tREF\tALT\tDP\nchr1\t"));

        let files = |vcf: &str, map: bool, t: &str| -> Vec<Vec<u8>> {
            let d = tempdir();
            let prefix = s_of(&d.path().join("o"));
            let mut args = vec![
                "vcf",
                "stats",
                vcf,
                "--out-prefix",
                &prefix,
                "--freq",
                "--site-quality",
                "--singletons",
                "--snp-density",
                "500",
                "-t",
                t,
            ];
            if map {
                args.extend(["--chrom-map", &fwd]);
            }
            gloci().args(&args).assert().success();
            ["frq", "lqual", "singletons", "snpden"]
                .iter()
                .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap())
                .collect()
        };
        let a = files(input, true, "1");
        assert_eq!(a, files(input, true, "8"), "stats {input}: threads");
        assert_eq!(a, files(&vchr, false, "1"), "stats {input}: map vs renamed");
        assert!(a[0].windows(5).any(|w| w == b"\nchr1"), "{input}");

        let m = |vcf: &str, map: bool| {
            let d = tempdir();
            let prefix = s_of(&d.path().join("m"));
            let mut args = vec![
                "vcf",
                "extract",
                vcf,
                "--genotype-matrix",
                "012",
                "--out-prefix",
                &prefix,
            ];
            if map {
                args.extend(["--chrom-map", &fwd]);
            }
            gloci().args(&args).assert().success();
            std::fs::read(format!("{prefix}.012.pos")).unwrap()
        };
        assert_eq!(m(input, true), m(&vchr, false), "012.pos {input}");
    }
    // cluster/spacing on a VCF with the map: the renamed record text
    let out = String::from_utf8(run_stdout(&["cluster", &vnum, "--chrom-map", &fwd])).unwrap();
    assert_eq!(
        out,
        String::from_utf8(run_stdout(&["cluster", &vchr])).unwrap()
    );
    // a map file with a bad line says where
    let bad = fx.plain("bad.map", "1\tchr1\nlonely\n");
    gloci()
        .args(["cluster", &bed_num, "--chrom-map", &bad])
        .assert()
        .code(2)
        .stderr(contains("bad.map:2:"));
}

/// `--sort` over a budget smaller than the input spills to disk and gives the
/// bytes of the in-memory sort. (The input is asserted to be more than twice
/// the budget, so the spill is not optional; the spill directory is empty
/// again afterwards.)
#[test]
fn p2_sort_spill_matches_in_memory() {
    // --sort over a budget smaller than the input spills to disk and gives
    // the bytes of the in-memory sort.
    let fx = P2::new();
    let big = fx.plain("big.bed", &bed_text(1_600_000));
    assert!(std::fs::metadata(&big).unwrap().len() > 2 * (16 << 20));
    let tmp = tempdir();
    let spill = [
        "cluster",
        &big,
        "--sort",
        "--max-mem",
        "16M",
        "--tmp-dir",
        tmp.path().to_str().unwrap(),
    ];
    let spilled: Vec<_> = ["1", "8"].iter().map(|t| stdout_with(&spill, t)).collect();
    assert_eq!(spilled[0], spilled[1], "cluster --sort spilled");
    assert_eq!(
        spilled[0],
        stdout_with(&["cluster", &big, "--sort"], "2"),
        "spilled vs in-memory"
    );
    assert_eq!(
        std::fs::read_dir(tmp.path()).unwrap().count(),
        0,
        "spill files left"
    );
}

/// Fixtures large enough to be split: 60k records is more than one 16384-record
/// batch, the BGZF copies span dozens of blocks, and chr1's 40k records make a
/// groupby run that crosses batch boundaries. Output at `-t 1` must equal
/// `-t 8`.
#[test]
fn p2_determinism_large_fixture() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t1000000\nchr2\t1000000\n");
    let mut a = String::new();
    let mut b = String::new();
    for (chrom, n) in [("chr1", 40_000usize), ("chr2", 20_000)] {
        for i in 0..n {
            let st = 100 + i * 23 + (i % 5) * 3;
            writeln!(
                a,
                "{chrom}\t{st}\t{}\ta{i}\t{}\t{}",
                st + 10 + (i % 7) * 9,
                i % 100,
                if i % 3 == 0 { '-' } else { '+' }
            )
            .unwrap();
            if i % 2 == 0 {
                let bs = 105 + i * 31;
                writeln!(b, "{chrom}\t{bs}\t{}\tb{i}\t1\t+", bs + 40).unwrap();
            }
        }
    }
    let a_plain = fx.plain("a.bed", &a);
    let a_gz = fx.bgzf("a.bed.gz", &a);
    let a_idx = fx.indexed("a_idx.bed.gz", &a);
    let b_plain = fx.plain("b.bed", &b);
    let b_gz = fx.bgzf("b.bed.gz", &b);
    assert!(
        std::fs::metadata(&a_gz).unwrap().len() > 3 * 65_536,
        "several BGZF blocks"
    );
    let run = |args: &[&str], t: &str| -> Vec<u8> {
        let mut c = gloci();
        c.args(args).args(["-t", t]);
        let o = c.output().unwrap();
        assert_eq!(
            o.status.code(),
            Some(0),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        o.stdout
    };
    for (name, tpl) in p2_bed_cmds(&genome) {
        // per-base output is one line per genome base: too long to repeat
        if name == "genomecov --per-base" {
            continue;
        }
        let sub = |ai: &str, bi: &str| -> Vec<String> {
            tpl.iter()
                .map(|x| match *x {
                    "{a}" => ai.to_string(),
                    "{b}" => bi.to_string(),
                    o => o.to_string(),
                })
                .collect()
        };
        let plain_args = sub(&a_plain, &b_plain);
        let refs: Vec<&str> = plain_args.iter().map(String::as_str).collect();
        let base = run(&refs, "1");
        // the default genomecov is a short histogram; the others are long
        let want = if name == "genomecov" { 50 } else { 100_000 };
        assert!(base.len() > want, "{name}: small output");
        assert_eq!(base, run(&refs, "8"), "{name}: -t 8");
        for (ai, bi) in [(&a_gz, &b_gz), (&a_idx, &b_gz)] {
            let args = sub(ai, bi);
            let refs: Vec<&str> = args.iter().map(String::as_str).collect();
            assert_eq!(base, run(&refs, "8"), "{name}: {ai} -t 8");
        }
    }
    // the groupby run over chr1 crosses the batch boundary and is one line
    let gb =
        String::from_utf8(run(&["groupby", &a_plain, "-g", "1", "-c", "4:count"], "8")).unwrap();
    assert_eq!(gb, "chr1\t40000\nchr2\t20000\n");

    // the same for a VCF of 40k records
    let mut v = String::from(P2_VCF_HEAD);
    for i in 0..40_000usize {
        writeln!(
            v,
            "chr1\t{}\t.\tA\tC\t{}\tPASS\tDP={}\tGT:DP\t0/1:{}\t1/1:3\t./.:4",
            100 + i * 7,
            20 + i % 50,
            i % 30,
            i % 9
        )
        .unwrap();
    }
    let v_plain = fx.plain("v.vcf", &v);
    let v_gz = fx.bgzf("v.vcf.gz", &v);
    let v_idx = fx.indexed("v_idx.vcf.gz", &v);
    for tpl in [
        vec!["cluster", "-d", "3"],
        vec!["spacing"],
        vec!["vcf", "extract", "--info", "DP"],
        vec!["vcf", "extract", "--format", "GT,DP"],
    ] {
        let with = |f: &str| -> Vec<String> {
            let mut a: Vec<String> = tpl.iter().map(|s| s.to_string()).collect();
            a.insert(if tpl[0] == "vcf" { 2 } else { 1 }, f.to_string());
            a
        };
        let args = with(&v_plain);
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let base = run(&refs, "1");
        assert!(base.len() > 100_000);
        assert_eq!(base, run(&refs, "8"), "{tpl:?}");
        for f in [&v_gz, &v_idx] {
            let args = with(f);
            let refs: Vec<&str> = args.iter().map(String::as_str).collect();
            assert_eq!(base, run(&refs, "8"), "{tpl:?} {f}");
        }
    }
    // vcf concat: a second file that starts after the first ends, so the
    // 40k-record first file spans several batches and BGZF blocks
    let mut tail = String::from(P2_VCF_HEAD);
    for i in 0..3_000usize {
        writeln!(
            tail,
            "chr1\t{}\t.\tA\tC\t{}\tPASS\tDP={}\tGT:DP\t0/1:{}\t1/1:3\t./.:4",
            1_000_000 + i * 7,
            20 + i % 50,
            i % 30,
            i % 9
        )
        .unwrap();
    }
    let tail_plain = fx.plain("tail.vcf", &tail);
    let tail_gz = fx.bgzf("tail.vcf.gz", &tail);
    for (first, second) in [
        (&v_plain, &tail_plain),
        (&v_gz, &tail_gz),
        (&v_idx, &tail_gz),
    ] {
        let args = ["vcf", "concat", first.as_str(), second.as_str()];
        let base = run(&args, "1");
        assert!(base.len() > 100_000, "concat output");
        assert_eq!(base, run(&args, "8"), "concat {first}");
    }
    // vcf stats files
    let stats = |f: &str, t: &str| -> Vec<Vec<u8>> {
        let d = tempdir();
        let prefix = s_of(&d.path().join("o"));
        gloci()
            .args([
                "vcf",
                "stats",
                f,
                "--out-prefix",
                &prefix,
                "--freq",
                "--het",
                "--site-quality",
                "-t",
                t,
            ])
            .assert()
            .success();
        ["frq", "het", "lqual"]
            .iter()
            .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap())
            .collect()
    };
    let base = stats(&v_plain, "1");
    assert_eq!(base[0].iter().filter(|&&c| c == b'\n').count(), 40_001);
    assert_eq!(base, stats(&v_plain, "8"));
    assert_eq!(base, stats(&v_gz, "8"));
    assert_eq!(base, stats(&v_idx, "8"));
}

/// `tpl` with `{IN}` replaced by `input` exits 0 and writes its output.
fn assert_succeeds(name: &str, tpl: &[&str], input: &str, sink: &Sink) {
    let d = tempdir();
    let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", input)).collect();
    match sink {
        Sink::File(ext) => args.extend(["-o".into(), s_of(&d.path().join(format!("out.{ext}")))]),
        Sink::Prefix => args.extend(["--out-prefix".into(), s_of(&d.path().join("out"))]),
    }
    let o = gloci().args(&args).output().unwrap();
    assert_eq!(
        o.status.code(),
        Some(0),
        "{name}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(
        std::fs::read_dir(d.path()).unwrap().count() > 0,
        "{name}: no output"
    );
}

/// The BCF magic after BGZF decompression.
fn is_bcf(bytes: &[u8]) -> bool {
    let mut all = Vec::new();
    noodles_bgzf::io::Reader::new(bytes)
        .read_to_end(&mut all)
        .is_ok()
        && all.starts_with(b"BCF\x02\x02")
}

/// `vcf concat` writes the format `-O` or the `-o` extension names; the first
/// input's format is only the default.
#[test]
fn p2_concat_output_format() {
    let fx = P2::new();
    let (t1, t2) = (p2_vcf(1..=1, 5), p2_vcf(2..=2, 5));
    let (v1, v2) = (fx.plain("1.vcf", &t1), fx.plain("2.vcf", &t2));
    let (b1, b2) = (fx.bcf("1.bcf", &t1), fx.bcf("2.bcf", &t2));
    let (z1, z2) = (fx.bgzf("1.vcf.gz", &t1), fx.bgzf("2.vcf.gz", &t2));
    let d = tempdir();
    let to = |name: &str, inputs: [&str; 2], extra: &[&str]| -> Vec<u8> {
        let out = s_of(&d.path().join(name));
        let mut full = vec!["vcf", "concat", inputs[0], inputs[1]];
        full.extend_from_slice(extra);
        full.extend(["-o", &out]);
        gloci().args(&full).assert().success();
        std::fs::read(&out).unwrap()
    };
    // -o extension decides, whatever the inputs are
    for inputs in [[&v1, &v2], [&z1, &z2], [&b1, &b2]] {
        let inputs = [inputs[0].as_str(), inputs[1].as_str()];
        assert!(is_bcf(&to("x.bcf", inputs, &[])), "{inputs:?} -o x.bcf");
        assert!(
            is_bcf(&to("y.out", inputs, &["-O", "bcf"])),
            "{inputs:?} -O bcf"
        );
        let gz = to("x.vcf.gz", inputs, &[]);
        assert!(!is_bcf(&gz), "{inputs:?} -o x.vcf.gz is VCF");
        let mut text = Vec::new();
        noodles_bgzf::io::Reader::new(&gz[..])
            .read_to_end(&mut text)
            .unwrap();
        assert!(
            text.starts_with(b"##fileformat=VCF"),
            "{inputs:?}: BGZF VCF text"
        );
        // -O vcf is plain text even for BCF input
        let plain = to("z.out", inputs, &["-O", "vcf"]);
        assert!(plain.starts_with(b"##fileformat=VCF"), "{inputs:?} -O vcf");
    }
    // the default (stdout) follows the first input
    let o = gloci().args(["vcf", "concat", &b1, &b2]).output().unwrap();
    assert!(o.status.success());
    assert!(is_bcf(&o.stdout), "BCF in, BCF out by default");
    let o = gloci().args(["vcf", "concat", &v1, &v2]).output().unwrap();
    assert!(o.stdout.starts_with(b"##fileformat=VCF"));
}

// ---------------------------------------------------------------------------
// P3a: alignments (BAM, CRAM, SAM) and FASTA commands
// ---------------------------------------------------------------------------

fn p3a_inputs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn p3a_aln(name: &str) -> String {
    s_of(&p3a_inputs().join("aln").join(name))
}

fn p3a_fa(name: &str) -> String {
    s_of(&p3a_inputs().join("fa").join(name))
}

/// SAM text, record lines only, with `@PG` command lines dropped (they name `-t`).
/// BAM input (gzip magic) is decoded through `gloci view` first.
fn p3a_norm(bytes: &[u8]) -> Vec<u8> {
    let text = if bytes.starts_with(&[0x1f, 0x8b]) {
        let d = tempdir();
        let f = write(d.path(), "x.bam", bytes);
        gloci()
            .args(["view", "-O", "sam", &s_of(&f)])
            .output()
            .unwrap()
            .stdout
    } else {
        bytes.to_vec()
    };
    String::from_utf8_lossy(&text)
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .flat_map(|l| [l, "\n"])
        .collect::<String>()
        .into_bytes()
}

/// The three alignment formats over the same 1210 reads, plus a FASTA and BED to go with them.
struct P3a {
    dir: tempfile::TempDir,
    cram: String,
    bam: String,
    sam: String,
    fa: String,
    fa_gz: String,
    bed: String,
    masks: String,
}

impl P3a {
    fn new() -> P3a {
        let dir = tempdir();
        let d = dir.path();
        let cram = d.join("many.cram");
        std::fs::copy(p3a_aln("many.cram"), &cram).unwrap();
        std::fs::copy(p3a_aln("many.cram.crai"), d.join("many.cram.crai")).unwrap();
        let fa = d.join("ref.fa");
        std::fs::copy(p3a_fa("ref.fa"), &fa).unwrap();
        std::fs::copy(p3a_fa("ref.fa.fai"), d.join("ref.fa.fai")).unwrap();
        let fa_gz = d.join("ref.fa.gz");
        for ext in ["", ".fai", ".gzi"] {
            std::fs::copy(
                p3a_fa(&format!("ref.fa.gz{ext}")),
                d.join(format!("ref.fa.gz{ext}")),
            )
            .unwrap();
        }
        let sam = d.join("many.sam");
        let bam = d.join("many.bam");
        for out in [&sam, &bam] {
            gloci()
                .args(["view", &s_of(&cram), "--reference", &s_of(&fa), "-o"])
                .arg(out)
                .assert()
                .success();
        }
        gloci().args(["index", &s_of(&bam)]).assert().success();
        // strictly sorted copies, for `sort -m` (many.cram is only position-sorted)
        for (src, name) in [(&bam, "sorted.bam"), (&sam, "sorted.sam")] {
            gloci()
                .args(["sort", &s_of(src), "-o", &s_of(&d.join(name))])
                .assert()
                .success();
        }
        let bed = plain_file(
            d,
            "t.bed",
            "chrM\t0\t100\nchr1\t0\t300\nchr1\t600\t900\nchr2\t100\t900\nchr10\t0\t500\nchrX\t100\t700\n",
        );
        let masks = plain_file(
            d,
            "masks.bed",
            "chrM\t5\t20\nchr1\t0\t30\nchr1\t100\t140\nchr2\t500\t700\nchr10\t10\t20\n",
        );
        P3a {
            cram: s_of(&cram),
            bam: s_of(&bam),
            sam: s_of(&sam),
            fa: s_of(&fa),
            fa_gz: s_of(&fa_gz),
            bed: s_of(&bed),
            masks: s_of(&masks),
            dir,
        }
    }

    fn path(&self, name: &str) -> String {
        s_of(&self.dir.path().join(name))
    }

    /// `(kind, input, extra args)` for each alignment format.
    fn kinds(&self) -> Vec<(&'static str, String, Vec<String>)> {
        vec![
            ("bam", self.bam.clone(), vec![]),
            ("sam", self.sam.clone(), vec![]),
            (
                "cram",
                self.cram.clone(),
                vec!["--reference".into(), self.fa.clone()],
            ),
        ]
    }

    /// Commands over one alignment input, as `(name, args)`. `-t` is added by the caller.
    fn alignment_cmds(&self, kind: &str, input: &str) -> Vec<(String, Vec<String>)> {
        let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
        let bed = self.bed.as_str();
        let mut c: Vec<(String, Vec<String>)> = vec![
            ("convert".into(), v(&["convert", input])),
            ("convert --bed12".into(), v(&["convert", input, "--bed12"])),
            ("view".into(), v(&["view", "-O", "sam", input])),
            (
                "view --min-mapq".into(),
                v(&["view", "-O", "sam", "--min-mapq", "30", input]),
            ),
            (
                "intersect A".into(),
                v(&["intersect", "-a", input, "-b", bed]),
            ),
            (
                "intersect B".into(),
                v(&["intersect", "-a", bed, "-b", input]),
            ),
            (
                "intersect counts".into(),
                v(&["intersect", "-a", bed, "-b", input, "--mode", "count"]),
            ),
            ("coverage".into(), v(&["coverage", "-a", bed, "-b", input])),
            ("genomecov".into(), v(&["genomecov", input])),
            (
                "genomecov bg".into(),
                v(&["genomecov", input, "--bedgraph"]),
            ),
            ("sort".into(), v(&["sort", "-O", "sam", input])),
            ("sort -n".into(), v(&["sort", "-n", "-O", "sam", input])),
            ("sort -N".into(), v(&["sort", "-N", "-O", "sam", input])),
            (
                "sort small budget".into(),
                v(&["sort", "--max-mem", "16M", "-O", "sam", input]),
            ),
        ];
        if kind != "cram" {
            let sorted = input.replace("many.", "sorted.");
            c.push((
                "sort -m".into(),
                v(&["sort", "-m", "-O", "sam", &sorted, &sorted]),
            ));
        }
        if kind != "sam" {
            c.push((
                "view region".into(),
                v(&["view", "-O", "sam", input, "chr1:1-500", "chr2:100-900"]),
            ));
        }
        if kind == "bam" {
            let sorted = input.replace("many.", "sorted.");
            c.push(("cat".into(), v(&["cat", &sorted, &sorted])));
        }
        c
    }

    fn fasta_cmds(&self) -> Vec<(String, Vec<String>)> {
        let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
        let bed12 = p3a_fa("gf12.bed");
        let gf = p3a_fa("gf.bed");
        let nuc = p3a_fa("nuc.bed");
        let mut c = Vec::new();
        for (tag, fa) in [("plain", &self.fa), ("bgzf", &self.fa_gz)] {
            c.push((
                format!("getfasta {tag}"),
                v(&["getfasta", "--fasta", fa, &gf, "-s", "--name"]),
            ));
            c.push((
                format!("getfasta --split {tag}"),
                v(&[
                    "getfasta",
                    "--fasta",
                    fa,
                    &bed12,
                    "--split",
                    "-s",
                    "--bed-out",
                ]),
            ));
            c.push((
                format!("getfasta --tab {tag}"),
                v(&["getfasta", "--fasta", fa, &gf, "--tab"]),
            ));
            c.push((
                format!("nuc {tag}"),
                v(&["nuc", "--fasta", fa, &nuc, "-s", "--seq", "--pattern", "AC"]),
            ));
            c.push((
                format!("maskfasta {tag}"),
                v(&["maskfasta", "--fasta", fa, "-b", &self.masks]),
            ));
            c.push((
                format!("maskfasta --soft {tag}"),
                v(&["maskfasta", "--fasta", fa, "-b", &self.masks, "--soft"]),
            ));
        }
        c
    }
}

fn p3a_run_t(args: &[String], t: &str) -> std::process::Output {
    let opt: Vec<String> = match t {
        "1" => vec!["-t".into(), "1".into()],
        "2" => vec!["-t2".into()],
        n => vec![format!("--threads={n}")],
    };
    gloci().args(args).args(&opt).output().unwrap()
}

/// Every alignment and FASTA command gives the same bytes at `-t 1`, `-t2` and `--threads=8`,
/// on BAM (indexed), SAM and CRAM (indexed); the three formats of the same reads agree for the
/// commands whose output does not name the format.
#[test]
fn p3a_determinism_matrix() {
    let fx = P3a::new();
    let mut ran = 0;
    let mut by_cmd: std::collections::BTreeMap<String, Vec<(String, Vec<u8>)>> = Default::default();
    for (kind, input, extra) in fx.kinds() {
        for (name, mut args) in fx.alignment_cmds(kind, &input) {
            args.extend(extra.iter().cloned());
            let outs: Vec<Vec<u8>> = ["1", "2", "8"]
                .iter()
                .map(|t| {
                    let o = p3a_run_t(&args, t);
                    assert_eq!(
                        o.status.code(),
                        Some(0),
                        "{kind} {name} -t {t}: {}",
                        String::from_utf8_lossy(&o.stderr)
                    );
                    p3a_norm(&o.stdout)
                })
                .collect();
            assert!(
                outs[0] == outs[1] && outs[0] == outs[2],
                "{kind} {name}: output differs between thread counts"
            );
            ran += 1;
            by_cmd
                .entry(name)
                .or_default()
                .push((kind.to_string(), outs[0].clone()));
        }
    }
    // BED output does not name the format: BAM, SAM and CRAM of the same reads agree.
    for name in [
        "convert",
        "convert --bed12",
        "coverage",
        "genomecov",
        "genomecov bg",
        "intersect B",
        "intersect counts",
    ] {
        let runs = &by_cmd[name];
        assert_eq!(runs.len(), 3, "{name}");
        assert!(!runs[0].1.is_empty(), "{name} wrote nothing");
        assert!(
            runs.iter().all(|r| r.1 == runs[0].1),
            "{name}: formats disagree"
        );
    }
    // The same reads sorted in each format are the same reads.
    for name in ["sort", "sort -n", "view"] {
        let runs = &by_cmd[name];
        let recs = |b: &[u8]| -> Vec<String> {
            String::from_utf8_lossy(b)
                .lines()
                .filter(|l| !l.starts_with('@'))
                .map(String::from)
                .collect()
        };
        assert!(
            runs.iter().all(|r| recs(&r.1) == recs(&runs[0].1)),
            "{name}: formats disagree"
        );
    }
    for (name, args) in fx.fasta_cmds() {
        let outs: Vec<Vec<u8>> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = p3a_run_t(&args, t);
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{name} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                o.stdout
            })
            .collect();
        assert!(!outs[0].is_empty(), "{name}: no output");
        assert!(
            outs[0] == outs[1] && outs[0] == outs[2],
            "{name}: output differs between thread counts"
        );
        ran += 1;
    }
    assert!(ran > 50, "ran {ran} cases");

    // Written files (-o) are the same bytes at every thread count too: BAM output compressed
    // by the parallel writer, and a BGZF FASTA.
    let files: Vec<(&str, Vec<String>)> = vec![
        ("sort", vec!["sort".into(), fx.bam.clone()]),
        (
            "sort spill",
            vec![
                "sort".into(),
                "--max-mem".into(),
                "16M".into(),
                fx.sam.clone(),
                "-O".into(),
                "bam".into(),
            ],
        ),
        (
            "convert-view",
            vec![
                "view".into(),
                fx.cram.clone(),
                "--reference".into(),
                fx.fa.clone(),
            ],
        ),
        (
            "maskfasta",
            vec![
                "maskfasta".into(),
                "--fasta".into(),
                fx.fa.clone(),
                "-b".into(),
                fx.masks.clone(),
            ],
        ),
    ];
    for (name, args) in files {
        let ext = if name == "maskfasta" { "fa.gz" } else { "bam" };
        let outs: Vec<Vec<u8>> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let out = fx.path(&format!("det-{name}-{t}.{ext}"));
                let mut a = args.clone();
                a.extend(["-o".into(), out.clone()]);
                let o = p3a_run_t(&a, t);
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{name} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                let bytes = std::fs::read(&out).unwrap();
                // the @PG line names -t; compare the records and the compressed layout of FASTA
                if ext == "bam" {
                    p3a_norm(&bytes)
                } else {
                    bytes
                }
            })
            .collect();
        assert!(
            outs[0] == outs[1] && outs[0] == outs[2],
            "{name}: written file differs between thread counts"
        );
    }
}

/// A SAM (and its BAM, indexed when sorted) of `n` 100 bp reads over chr1 and chr2 (10 Mb each),
/// in coordinate order or shuffled.
fn p3a_big_alignments(dir: &Path, n: usize, shuffled: bool) -> (String, String) {
    let mut x: u64 = 0xC0FFEE;
    let mut next = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x >> 20
    };
    let mut reads: Vec<(usize, u64, u64)> = (0..n)
        .map(|_| ((next() % 2) as usize, 1 + next() % 9_999_000, next()))
        .collect();
    if !shuffled {
        reads.sort_by_key(|r| (r.0, r.1, r.2 & 1));
    }
    let mut sam = format!(
        "@HD\tVN:1.6\tSO:{}\n@SQ\tSN:chr1\tLN:10000000\n@SQ\tSN:chr2\tLN:10000000\n",
        if shuffled { "unsorted" } else { "coordinate" }
    );
    for (i, (c, pos, r)) in reads.iter().enumerate() {
        let seq: String = (0..100)
            .map(|k| b"ACGT"[((r >> (k % 30)) as usize + k) & 3] as char)
            .collect();
        writeln!(
            sam,
            "read{i}\t{}\tchr{}\t{pos}\t60\t100M\t*\t0\t0\t{seq}\t{}",
            if r & 1 == 0 { 0 } else { 16 },
            c + 1,
            "I".repeat(100)
        )
        .unwrap();
    }
    let name = if shuffled { "shuf" } else { "big" };
    let sam_path = plain_file(dir, &format!("{name}.sam"), &sam);
    let bam_path = dir.join(format!("{name}.bam"));
    gloci()
        .args(["view", &s_of(&sam_path), "-o", &s_of(&bam_path)])
        .assert()
        .success();
    if !shuffled {
        gloci().args(["index", &s_of(&bam_path)]).assert().success();
    }
    (s_of(&sam_path), s_of(&bam_path))
}

/// A 3 Mb FASTA (`big.fa`, indexed) and BED of 6000 intervals over it, for output that outgrows a
/// pipe.
fn p3a_big_fasta(dir: &Path) -> (String, String) {
    let mut x: u64 = 99;
    let mut fa = String::from(">big\n");
    let bases = b"ACGT";
    for line in 0..50_000 {
        for _ in 0..60 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            fa.push(bases[(x >> 60) as usize & 3] as char);
        }
        fa.push('\n');
        let _ = line;
    }
    let f = plain_file(dir, "big.fa", &fa);
    gloci().args(["index", &s_of(&f)]).assert().success();
    let mut bed = String::new();
    for i in 0..6000 {
        writeln!(bed, "big\t{}\t{}\tn{i}\t0\t+", i * 400, i * 400 + 380).unwrap();
    }
    let b = plain_file(dir, "big.bed", &bed);
    (s_of(&f), s_of(&b))
}

/// Closing stdout after the first bytes is a quiet success (exit 0, nothing on stderr) for every
/// alignment and FASTA command; each command writes far more than a pipe holds, so the reader
/// really does leave early.
#[test]
fn p3a_broken_pipe_matrix() {
    let d = tempdir();
    let (big_fa, big_bed) = p3a_big_fasta(d.path());
    let (big_sam, big_bam) = p3a_big_alignments(d.path(), 12_000, false);
    let sorted_bam = big_bam.clone();
    let cram = p3a_aln("manyref.cram");
    let cram_fa = p3a_fa("manyref.fa");
    let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
    let one_per_chrom = d.path().join("chroms.bed");
    std::fs::write(&one_per_chrom, "chr1\t0\t10000000\nchr2\t0\t10000000\n").unwrap();
    let chroms = s_of(&one_per_chrom);
    let mut cmds: Vec<(String, Vec<String>)> = Vec::new();
    for (kind, input) in [("bam", &big_bam), ("sam", &big_sam)] {
        for (name, a) in [
            ("convert", v(&["convert", input])),
            ("view", v(&["view", "-O", "sam", input])),
            ("view bam", v(&["view", "-O", "bam", input])),
            ("intersect A", v(&["intersect", "-a", input, "-b", &chroms])),
            (
                "intersect B",
                v(&["intersect", "-a", &chroms, "-b", input, "--mode", "pairs"]),
            ),
            ("genomecov per-base", v(&["genomecov", input, "--per-base"])),
            ("sort", v(&["sort", input])),
            ("sort -n", v(&["sort", "-n", input])),
        ] {
            cmds.push((format!("{kind} {name}"), a));
        }
    }
    cmds.push((
        "bam sort -m".into(),
        v(&["sort", "-m", &sorted_bam, &sorted_bam]),
    ));
    cmds.push(("bam cat".into(), v(&["cat", &sorted_bam, &sorted_bam])));
    cmds.push((
        "bam view region".into(),
        v(&["view", "-O", "sam", &big_bam, "chr1:1-10000000"]),
    ));
    for (name, a) in [
        (
            "view",
            v(&["view", "-O", "sam", &cram, "--reference", &cram_fa]),
        ),
        (
            "genomecov per-base",
            v(&["genomecov", &cram, "--per-base", "--reference", &cram_fa]),
        ),
        (
            "sort",
            v(&["sort", &cram, "--reference", &cram_fa, "-O", "sam"]),
        ),
    ] {
        cmds.push((format!("cram {name}"), a));
    }
    cmds.push((
        "getfasta".into(),
        v(&["getfasta", "--fasta", &big_fa, &big_bed]),
    ));
    cmds.push((
        "getfasta --tab".into(),
        v(&["getfasta", "--fasta", &big_fa, &big_bed, "--tab"]),
    ));
    cmds.push((
        "nuc".into(),
        v(&["nuc", "--fasta", &big_fa, &big_bed, "--seq"]),
    ));
    cmds.push((
        "maskfasta".into(),
        v(&["maskfasta", "--fasta", &big_fa, "-b", &big_bed]),
    ));

    for (name, args) in &cmds {
        let full = gloci().args(args).output().unwrap();
        assert_eq!(
            full.status.code(),
            Some(0),
            "{name}: {}",
            String::from_utf8_lossy(&full.stderr)
        );
        assert!(
            full.stdout.len() > 100_000,
            "{name}: only {} bytes of output",
            full.stdout.len()
        );
        for t in ["1", "4"] {
            let mut child = Command::new(bin())
                .args(args)
                .args(["-t", t])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            let mut first = [0u8; 64];
            child
                .stdout
                .as_mut()
                .unwrap()
                .read_exact(&mut first)
                .unwrap();
            drop(child.stdout.take());
            let out = child.wait_with_output().unwrap();
            assert_eq!(
                out.status.code(),
                Some(0),
                "{name} -t {t}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert!(
                out.stderr.is_empty(),
                "{name} -t {t}: stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
    }
}

/// Records of SAM text (no header lines).
fn p3a_records(sam: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(&p3a_norm(sam))
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(String::from)
        .collect()
}

/// A header-only BAM and SAM, an empty BED and FASTA inputs without sequences: every command
/// exits 0 (or says clearly why not), writes a valid empty result and prints nothing on stderr.
#[test]
fn p3a_empty_inputs_matrix() {
    let d = tempdir();
    let dir = d.path();
    let hdr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n";
    let sam = s_of(&plain_file(dir, "h.sam", hdr));
    let bam = s_of(&dir.join("h.bam"));
    gloci().args(["view", &sam, "-o", &bam]).assert().success();
    gloci().args(["index", &bam]).assert().success();
    let empty_bed = s_of(&plain_file(dir, "e.bed", ""));
    let bed = s_of(&plain_file(dir, "t.bed", "chr1\t0\t100\nchr2\t10\t20\n"));

    for t in ["1", "4"] {
        for input in [&sam, &bam] {
            let ok = |args: &[&str]| -> std::process::Output {
                let o = gloci().args(args).args(["-t", t]).output().unwrap();
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{args:?} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                assert!(
                    o.stderr.is_empty(),
                    "{args:?} -t {t}: stderr {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                o
            };
            // BED-writing commands: no rows
            for args in [
                vec!["convert", input],
                vec!["convert", input, "--bed12"],
                vec!["convert", input, "--blocks"],
                vec!["intersect", "-a", &bed, "-b", input],
                vec!["intersect", "-a", &bed, "-b", input, "--mode", "pairs"],
                vec!["intersect", "-a", &empty_bed, "-b", input],
                vec!["coverage", "-a", &empty_bed, "-b", input],
            ] {
                assert!(ok(&args).stdout.is_empty(), "{args:?} -t {t}: wrote rows");
            }
            // commands that report on A: one zero row per A record
            let cov = ok(&["coverage", "-a", &bed, "-b", input]);
            assert_eq!(
                String::from_utf8_lossy(&cov.stdout),
                "chr1\t0\t100\t0\t0\t100\t0\nchr2\t10\t20\t0\t0\t10\t0\n",
                "coverage -t {t}"
            );
            let gc = ok(&["genomecov", input]);
            let text = String::from_utf8_lossy(&gc.stdout).into_owned();
            assert!(
                text.contains("chr1\t0\t1000\t1000\t1\n"),
                "genomecov {text}"
            );
            assert!(
                text.contains("genome\t0\t1500\t1500\t1\n"),
                "genomecov {text}"
            );
            let bg = ok(&["genomecov", input, "--bedgraph"]);
            assert!(bg.stdout.is_empty(), "genomecov --bedgraph wrote rows");
            let bga = ok(&["genomecov", input, "--bga"]);
            assert_eq!(
                String::from_utf8_lossy(&bga.stdout),
                "chr1\t0\t1000\t0\nchr2\t0\t500\t0\n"
            );
            // alignment-writing commands: the header and no reads
            let sam_out = |args: &[&str]| {
                let o = ok(args);
                let norm = p3a_norm(&o.stdout);
                let text = String::from_utf8_lossy(&norm).into_owned();
                assert!(
                    text.contains("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n"),
                    "{args:?}: {text}"
                );
                assert!(
                    p3a_records(&o.stdout).is_empty(),
                    "{args:?} -t {t}: wrote reads"
                );
            };
            sam_out(&["view", "-O", "sam", input]);
            sam_out(&["view", "-O", "bam", input]);
            if input == &bam {
                sam_out(&["view", "-O", "sam", input, "chr1:1-1000"]);
            }
            sam_out(&["sort", input, "-O", "sam"]);
            sam_out(&["sort", "-n", input, "-O", "sam"]);
            sam_out(&["sort", "-N", input, "-O", "sam"]);
            sam_out(&["sort", "-m", input, input, "-O", "sam"]);
            sam_out(&["intersect", "-a", input, "-b", &empty_bed, "-O", "sam"]);
            sam_out(&["intersect", "-a", input, "-b", &bed, "-O", "sam"]);
            if input == &bam {
                sam_out(&["cat", input, input]);
            }
            // -n and sort -m of the header-only file written to a file is a BAM `view` reads back
            let out = s_of(&dir.join(format!("sorted-{t}.bam")));
            ok(&["sort", input, "-o", &out, "--write-index"]);
            assert!(dir.join(format!("sorted-{t}.bam.bai")).exists());
            assert!(p3a_records(&std::fs::read(&out).unwrap()).is_empty());
        }
    }

    // reads on one contig only: the other contig's rows are still there
    let one = format!("{hdr}r1\t0\tchr1\t11\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    let one_sam = s_of(&plain_file(dir, "one.sam", &one));
    let o = gloci()
        .args(["genomecov", &one_sam, "--bga"])
        .output()
        .unwrap();
    assert_eq!(
        String::from_utf8_lossy(&o.stdout),
        "chr1\t0\t10\t0\nchr1\t10\t15\t1\nchr1\t15\t1000\t0\nchr2\t0\t500\t0\n"
    );

    // FASTA without sequences. An empty file is not a FASTA to gloci: it cannot be indexed, and
    // getfasta, nuc and maskfasta say so (exit 1 or 2) instead of writing nothing. index and
    // maskfasta give the one empty-input message (M1).
    let fa = s_of(&plain_file(dir, "empty.fa", ""));
    gloci()
        .args(["index", &fa])
        .assert()
        .code(1)
        .stderr(starts_with_error().and(contains("empty.fa is empty: it has no records")));
    assert!(!dir.join("empty.fa.fai").exists());
    for cmd in ["getfasta", "nuc"] {
        gloci()
            .args([cmd, "--fasta", &fa, &empty_bed])
            .assert()
            .code(2)
            .stderr(contains("has no index").and(contains("gloci index")));
    }
    for b in [&empty_bed, &bed] {
        gloci()
            .args(["maskfasta", "--fasta", &fa, "-b", b])
            .assert()
            .code(1)
            .stderr(contains("is empty: it has no records"));
    }
    // With an (empty) index beside it: nothing to read, so an empty BED is fine and any
    // interval names the chromosome the FASTA lacks.
    plain_file(dir, "empty.fa.fai", "");
    for t in ["1", "4"] {
        for cmd in ["getfasta", "nuc"] {
            let o = gloci()
                .args([cmd, "--fasta", &fa, &empty_bed, "-t", t])
                .output()
                .unwrap();
            assert_eq!(
                o.status.code(),
                Some(0),
                "{cmd} -t {t}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            assert!(o.stdout.is_empty() && o.stderr.is_empty(), "{cmd}");
            gloci()
                .args([cmd, "--fasta", &fa, &bed, "-t", t])
                .assert()
                .code(1)
                .stderr(contains("chr1 is not in the FASTA index"));
            gloci()
                .args([cmd, "--fasta", &fa, &bed, "--skip-invalid", "-t", t])
                .assert()
                .success()
                .stdout(predicate::function(|b: &[u8]| {
                    b.is_empty()
                        || b.starts_with(b"#1_usercol")
                            && b.iter().filter(|&&c| c == b'\n').count() == 1
                }))
                .stderr(contains("skipped 2 intervals outside the FASTA"));
        }
    }
    // A real FASTA with an empty BED (or empty mask list): getfasta and nuc write nothing,
    // maskfasta copies the FASTA unchanged.
    let ref_fa = p3a_fa("ref.fa");
    for t in ["1", "4"] {
        for cmd in ["getfasta", "nuc"] {
            let o = gloci()
                .args([cmd, "--fasta", &ref_fa, &empty_bed, "-t", t])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(0), "{cmd} -t {t}");
            assert!(o.stdout.is_empty() && o.stderr.is_empty(), "{cmd} -t {t}");
        }
        let o = gloci()
            .args(["maskfasta", "--fasta", &ref_fa, "-b", &empty_bed, "-t", t])
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(0));
        assert_eq!(
            o.stdout,
            std::fs::read(&ref_fa).unwrap(),
            "maskfasta with no masks"
        );
    }
    // A FASTA whose only record has no bases is refused by the index with the line number, and
    // maskfasta keeps the empty record.
    let hdr_only = s_of(&plain_file(dir, "hdr.fa", ">x\n"));
    gloci()
        .args(["index", &hdr_only])
        .assert()
        .code(1)
        .stderr(contains("hdr.fa line 1").and(contains("has no bases")));
    assert!(!dir.join("hdr.fa.fai").exists());
    gloci()
        .args(["maskfasta", "--fasta", &hdr_only, "-b", &empty_bed])
        .assert()
        .success()
        .stdout(">x\n");
}

fn starts_with_error() -> predicates::str::StartsWithPredicate {
    predicates::str::starts_with("error: ")
}

/// A BAM cut inside its last record but ending with a valid BGZF end-of-file block (so only the
/// record itself is short).
fn p3a_cut_record_bam(dir: &Path) -> String {
    let mut raw = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(p3a_aln("cat-big.bam")).unwrap())
        .read_to_end(&mut raw)
        .unwrap();
    raw.truncate(raw.len() - 7);
    s_of(&write(dir, "cut-record.bam", &bgzip(&raw)))
}

/// Each bad input fails every alignment command cleanly: exit 1, one `error:` line that says
/// what is wrong, nothing written, at `-t 1` and `-t 4`.
#[test]
fn p3a_failure_matrix_bad_inputs() {
    let fx = P3a::new();
    let dir = fx.path("bad");
    std::fs::create_dir(&dir).unwrap();
    let cut = p3a_cut_record_bam(Path::new(&dir));
    let half_cram = {
        let b = std::fs::read(&fx.cram).unwrap();
        s_of(&write(Path::new(&dir), "half.cram", &b[..b.len() / 2]))
    };
    let noeof_cram = {
        let b = std::fs::read(&fx.cram).unwrap();
        s_of(&write(Path::new(&dir), "noeof.cram", &b[..b.len() - 38]))
    };
    let bad_bams: Vec<(&str, String, &str)> = vec![
        (
            "no EOF block",
            p3a_aln("cat-noeof.bam"),
            "BGZF end-of-file marker missing",
        ),
        ("cut record", cut, "truncated"),
    ];
    let bad_crams: Vec<(&str, String, &str)> = vec![
        (
            "no EOF container",
            noeof_cram,
            "the CRAM end-of-file container is missing",
        ),
        (
            "half a file",
            half_cram,
            "the CRAM end-of-file container is missing",
        ),
        (
            "corrupt container",
            p3a_aln("many.badtail.cram"),
            "container block checksum mismatch",
        ),
    ];
    let bam_bed = s_of(&plain_file(
        Path::new(&dir),
        "chr12.bed",
        "chr1\t0\t1000\nchr2\t0\t1000\n",
    ));
    let ref_args = ["--reference", fx.fa.as_str()];
    let mut ran = 0;
    for (what, bad, msg) in bad_bams.iter().chain(bad_crams.iter()) {
        let is_cram = bad.ends_with(".cram");
        let bed = if is_cram {
            fx.bed.as_str()
        } else {
            bam_bed.as_str()
        };
        let cmds: Vec<(&str, Vec<&str>, Sink)> = vec![
            ("convert", vec!["convert", "{IN}"], Sink::File("bed")),
            ("view", vec!["view", "{IN}"], Sink::File("sam")),
            ("view bam", vec!["view", "{IN}"], Sink::File("bam")),
            (
                "intersect A",
                vec!["intersect", "-a", "{IN}", "-b", bed],
                Sink::File("bam"),
            ),
            (
                "intersect B",
                vec!["intersect", "-a", bed, "-b", "{IN}"],
                Sink::File("bed"),
            ),
            (
                "coverage",
                vec!["coverage", "-a", bed, "-b", "{IN}"],
                Sink::File("tsv"),
            ),
            ("genomecov", vec!["genomecov", "{IN}"], Sink::File("tsv")),
            ("sort", vec!["sort", "{IN}"], Sink::File("bam")),
            ("sort -n", vec!["sort", "-n", "{IN}"], Sink::File("bam")),
            (
                "sort -m",
                vec!["sort", "-m", "{IN}", "{IN}"],
                Sink::File("bam"),
            ),
            ("cat", vec!["cat", "{IN}", "{IN}"], Sink::File("bam")),
        ];
        for (name, mut tpl, sink) in cmds {
            if is_cram && matches!(name, "sort -m" | "cat") {
                continue; // BAM only
            }
            if *what == "cut record" && name == "cat" {
                continue; // cat copies compressed blocks and never decodes a record
            }
            if is_cram {
                tpl.extend(ref_args);
            }
            assert_fails_clean(name, what, &tpl, bad, &sink, msg);
            ran += 1;
        }
    }
    assert!(ran > 40, "ran {ran}");
}

/// An unsorted BAM or CRAM is refused with the `gloci sort` / `--sort` hint, and `--sort` gives
/// what the same file sorted first gives.
#[test]
fn p3a_failure_matrix_unsorted_alignments() {
    let fx = P3a::new();
    let d = fx.path("uns");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let ref_fa = p3a_fa("ref.fa");
    let bam_bed = s_of(&plain_file(d, "b.bed", "chrM\t0\t1000\nchr1\t0\t100000\n"));
    let cram_bed = s_of(&plain_file(d, "c.bed", "chr1\t0\t1000\n"));
    let sorted_bam = s_of(&d.join("sorted.bam"));
    gloci()
        .args(["sort", &p3a_aln("sort.bam"), "-o", &sorted_bam])
        .assert()
        .success();
    let sorted_cram = s_of(&d.join("sorted-cram.bam"));
    gloci()
        .args([
            "sort",
            &p3a_aln("unsorted.cram"),
            "--reference",
            &ref_fa,
            "-o",
            &sorted_cram,
        ])
        .assert()
        .success();
    let cases = [
        ("bam", p3a_aln("sort.bam"), sorted_bam, bam_bed, vec![]),
        (
            "cram",
            p3a_aln("unsorted.cram"),
            sorted_cram,
            cram_bed,
            vec!["--reference", ref_fa.as_str()],
        ),
    ];
    let mut ran = 0;
    for (kind, bad, sorted, bed, extra) in &cases {
        let tpls: Vec<(&str, Vec<&str>, Sink)> = vec![
            (
                "intersect A",
                vec!["intersect", "-a", "{IN}", "-b", bed],
                Sink::File("bam"),
            ),
            (
                "intersect B",
                vec!["intersect", "-a", bed, "-b", "{IN}"],
                Sink::File("bed"),
            ),
            (
                "coverage",
                vec!["coverage", "-a", bed, "-b", "{IN}"],
                Sink::File("tsv"),
            ),
            ("genomecov", vec!["genomecov", "{IN}"], Sink::File("tsv")),
        ];
        for (name, mut tpl, sink) in tpls {
            tpl.extend(extra.iter());
            // refused, with the hint, and nothing written
            assert_fails_clean(name, kind, &tpl, bad, &sink, "is not sorted");
            let probe: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", bad)).collect();
            let err = String::from_utf8_lossy(&gloci().args(&probe).output().unwrap().stderr)
                .into_owned();
            assert!(
                err.contains("hint: run: gloci sort ") && err.contains("or pass --sort"),
                "{kind} {name}: {err}"
            );
            // --sort equals the sorted file, at every thread count
            let want = {
                let a: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", sorted)).collect();
                let o = gloci().args(&a).output().unwrap();
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{kind} {name} sorted: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                p3a_norm(&o.stdout)
            };
            for t in ["1", "4"] {
                let mut a = probe.clone();
                a.extend(["--sort".into(), "-t".into(), t.into()]);
                let o = gloci().args(&a).output().unwrap();
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{kind} {name} --sort -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                assert!(!want.is_empty(), "{kind} {name}: nothing to compare");
                let got = p3a_norm(&o.stdout);
                if name == "intersect A" {
                    // The reads come out in (reference, start) order with ties in input order, not
                    // samtools' strand order, and the header keeps SO:unsorted: compare the reads
                    // as a set (see the report).
                    let set = |b: &[u8]| {
                        let mut v = p3a_records(b);
                        v.sort();
                        v
                    };
                    assert_eq!(set(&got), set(&want), "{kind} {name} --sort -t {t}");
                } else {
                    assert_eq!(got, want, "{kind} {name} --sort -t {t}");
                }
            }
            ran += 1;
        }
    }
    assert_eq!(ran, 8);
    // sort -m refuses an unsorted BAM, naming the record and the fix
    for t in ["1", "4"] {
        gloci()
            .args([
                "sort",
                "-m",
                &p3a_aln("sort.bam"),
                &p3a_aln("sort.bam"),
                "-t",
                t,
            ])
            .assert()
            .code(1)
            .stderr(
                contains("is not in coordinate order")
                    .and(contains("sort it first with: gloci sort")),
            );
    }
}

/// A missing index: a region query on a BAM or CRAM scans the file with one warning and gives the
/// indexed answer; getfasta and nuc refuse a FASTA without its `.fai` (or `.gzi`) and say how to
/// make it; per-chromosome runs at `-t 8` fall back to streaming with the same bytes.
#[test]
fn p3a_failure_matrix_missing_index() {
    let fx = P3a::new();
    let d = fx.path("noidx");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let bare_bam = d.join("many.bam");
    std::fs::copy(&fx.bam, &bare_bam).unwrap();
    let bare_cram = d.join("many.cram");
    std::fs::copy(&fx.cram, &bare_cram).unwrap();
    let (bare_bam, bare_cram) = (s_of(&bare_bam), s_of(&bare_cram));
    assert!(!Path::new(&format!("{bare_bam}.bai")).exists());
    for (kind, indexed, bare, extra) in [
        ("bam", &fx.bam, &bare_bam, vec![]),
        (
            "cram",
            &fx.cram,
            &bare_cram,
            vec!["--reference", fx.fa.as_str()],
        ),
    ] {
        for t in ["1", "8"] {
            for tpl in [
                vec!["view", "-O", "sam", "{IN}", "chr1:1-500", "chr2:100-900"],
                vec![
                    "intersect",
                    "-a",
                    "{IN}",
                    "-b",
                    fx.bed.as_str(),
                    "-O",
                    "sam",
                ],
                vec!["coverage", "-a", fx.bed.as_str(), "-b", "{IN}"],
                vec!["genomecov", "{IN}"],
                vec!["convert", "{IN}"],
            ] {
                let run = |input: &String| {
                    let mut a: Vec<String> = tpl.iter().map(|x| x.replace("{IN}", input)).collect();
                    a.extend(extra.iter().map(|s| s.to_string()));
                    a.extend(["-t".into(), t.into()]);
                    gloci().args(&a).output().unwrap()
                };
                let (with, without) = (run(indexed), run(bare));
                assert_eq!(
                    without.status.code(),
                    Some(0),
                    "{kind} {tpl:?} -t {t}: {}",
                    String::from_utf8_lossy(&without.stderr)
                );
                assert_eq!(
                    p3a_norm(&with.stdout),
                    p3a_norm(&without.stdout),
                    "{kind} {tpl:?} -t {t}"
                );
                let err = String::from_utf8_lossy(&without.stderr);
                if tpl[0] == "view" {
                    assert_eq!(
                        err.matches("scanning the whole file").count(),
                        1,
                        "{kind} {tpl:?} -t {t}: {err}"
                    );
                    assert!(err.contains(&format!("gloci index {bare}")), "{err}");
                }
            }
        }
    }
    // FASTA without .fai / .gzi
    let fa = d.join("ref.fa");
    std::fs::copy(&fx.fa, &fa).unwrap();
    let fagz = d.join("ref.fa.gz");
    std::fs::copy(&fx.fa_gz, &fagz).unwrap();
    std::fs::copy(format!("{}.fai", fx.fa_gz), d.join("ref.fa.gz.fai")).unwrap();
    let gf = p3a_fa("gf.bed");
    for cmd in ["getfasta", "nuc"] {
        for t in ["1", "4"] {
            let o = gloci()
                .args([cmd, "--fasta", &s_of(&fa), &gf, "-t", t])
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(2), "{cmd}: {err}");
            assert!(
                err.contains("has no index")
                    && err.contains("ref.fa.fai")
                    && err.contains("gloci index"),
                "{err}"
            );
            assert!(o.stdout.is_empty());
            let o = gloci()
                .args([cmd, "--fasta", &s_of(&fagz), &gf, "-t", t])
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(2), "{cmd}: {err}");
            assert!(
                err.contains("ref.fa.gz.gzi") && err.contains("gloci index"),
                "{err}"
            );
            assert!(o.stdout.is_empty());
        }
    }
    assert!(
        !d.join("ref.fa.fai").exists(),
        "gloci never builds an index implicitly"
    );
    // maskfasta streams and needs no index
    gloci()
        .args(["maskfasta", "--fasta", &s_of(&fa), "-b", &fx.masks])
        .assert()
        .success();
}

/// With the output directory missing every alignment and FASTA command exits 1, names the
/// directory and how to fix it, and creates nothing; `sort`'s own `--tmp-dir` is checked before
/// any input is read.
#[test]
fn p3a_failure_matrix_missing_output_dir() {
    let fx = P3a::new();
    let bed = fx.bed.as_str();
    let sorted_bam = fx.path("sorted.bam");
    let ref_args = ["--reference", fx.fa.as_str()];
    let gf = p3a_fa("gf.bed");
    let nuc = p3a_fa("nuc.bed");
    let mut cases: Vec<(&str, Vec<&str>, &str, Sink)> = vec![
        (
            "convert",
            vec!["convert", "{IN}"],
            &fx.bam,
            Sink::File("bed"),
        ),
        ("view", vec!["view", "{IN}"], &fx.bam, Sink::File("sam")),
        ("view bam", vec!["view", "{IN}"], &fx.sam, Sink::File("bam")),
        (
            "intersect A",
            vec!["intersect", "-a", "{IN}", "-b", bed],
            &fx.bam,
            Sink::File("bam"),
        ),
        (
            "intersect B",
            vec!["intersect", "-a", bed, "-b", "{IN}"],
            &fx.sam,
            Sink::File("bed"),
        ),
        (
            "coverage",
            vec!["coverage", "-a", bed, "-b", "{IN}"],
            &fx.bam,
            Sink::File("tsv"),
        ),
        (
            "genomecov",
            vec!["genomecov", "{IN}"],
            &fx.sam,
            Sink::File("tsv"),
        ),
        ("sort", vec!["sort", "{IN}"], &fx.bam, Sink::File("bam")),
        (
            "sort -n",
            vec!["sort", "-n", "{IN}"],
            &fx.sam,
            Sink::File("sam"),
        ),
        (
            "sort --write-index",
            vec!["sort", "--write-index", "{IN}"],
            &fx.bam,
            Sink::File("bam"),
        ),
        (
            "sort -m",
            vec!["sort", "-m", "{IN}", "{IN}"],
            &sorted_bam,
            Sink::File("bam"),
        ),
        (
            "cat",
            vec!["cat", "{IN}", "{IN}"],
            &sorted_bam,
            Sink::File("bam"),
        ),
        (
            "getfasta",
            vec!["getfasta", "--fasta", fx.fa.as_str(), "{IN}"],
            &gf,
            Sink::File("fa"),
        ),
        (
            "getfasta bgzf",
            vec!["getfasta", "--fasta", fx.fa.as_str(), "{IN}"],
            &gf,
            Sink::File("fa.gz"),
        ),
        (
            "nuc",
            vec!["nuc", "--fasta", fx.fa.as_str(), "{IN}"],
            &nuc,
            Sink::File("tsv"),
        ),
        (
            "maskfasta",
            vec!["maskfasta", "--fasta", fx.fa.as_str(), "-b", "{IN}"],
            &fx.masks,
            Sink::File("fa"),
        ),
        (
            "maskfasta --write-index",
            vec![
                "maskfasta",
                "--write-index",
                "--fasta",
                fx.fa.as_str(),
                "-b",
                "{IN}",
            ],
            &fx.masks,
            Sink::File("fa.gz"),
        ),
    ];
    let mut cram_cases: Vec<(&str, Vec<&str>, &str, Sink)> = vec![
        (
            "cram convert",
            vec!["convert", "{IN}"],
            &fx.cram,
            Sink::File("bed"),
        ),
        (
            "cram view",
            vec!["view", "{IN}"],
            &fx.cram,
            Sink::File("bam"),
        ),
        (
            "cram intersect",
            vec!["intersect", "-a", "{IN}", "-b", bed],
            &fx.cram,
            Sink::File("bam"),
        ),
        (
            "cram coverage",
            vec!["coverage", "-a", bed, "-b", "{IN}"],
            &fx.cram,
            Sink::File("tsv"),
        ),
        (
            "cram genomecov",
            vec!["genomecov", "{IN}"],
            &fx.cram,
            Sink::File("tsv"),
        ),
        (
            "cram sort",
            vec!["sort", "{IN}"],
            &fx.cram,
            Sink::File("bam"),
        ),
    ];
    for c in &mut cram_cases {
        c.1.extend(ref_args);
    }
    cases.extend(cram_cases);
    for (name, tpl, good, sink) in &cases {
        assert_missing_dir(name, tpl, good, sink);
    }

    // sort's --tmp-dir must exist
    let d = tempdir();
    let out = s_of(&d.path().join("o.bam"));
    for extra in [vec![], vec!["-n"], vec!["--max-mem", "16M"]] {
        let o = gloci()
            .args(["sort", &fx.bam, "--tmp-dir", "/no/such/tmp", "-o", &out])
            .args(&extra)
            .output()
            .unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert_eq!(o.status.code(), Some(1), "{err}");
        assert!(err.contains("/no/such/tmp"), "{err}");
    }
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
}

/// FASTA commands fail cleanly: intervals past a sequence or on a missing one, a truncated BGZF
/// FASTA, a mask past the end; exit 1, one `error:` line, nothing written, at every `-t`.
#[test]
fn p3a_failure_matrix_fasta() {
    let fx = P3a::new();
    let d = fx.path("fabad");
    let dir = Path::new(&d);
    std::fs::create_dir(dir).unwrap();
    let past = s_of(&plain_file(
        dir,
        "past.bed",
        "chr1\t10\t20\nchr1\t990\t1010\n",
    ));
    let absent = s_of(&plain_file(dir, "absent.bed", "chr1\t10\t20\nchrZ\t0\t5\n"));
    let gz = std::fs::read(&fx.fa_gz).unwrap();
    let cut = write(dir, "cut.fa.gz", &gz[..gz.len() / 2]);
    let cut_noeof = write(dir, "noeof.fa.gz", &gz[..gz.len() - 28]);
    let mut ran = 0;
    for (cmd, extra) in [
        ("getfasta", vec![]),
        ("getfasta", vec!["-s", "--name"]),
        ("getfasta", vec!["--split"]),
        ("nuc", vec![]),
    ] {
        for (bad, msg) in [
            (&past, "ends past the 1000 bases of chr1"),
            (&absent, "chrZ is not in the FASTA index"),
        ] {
            let mut tpl = vec![cmd, "--fasta", fx.fa.as_str(), "{IN}"];
            tpl.extend(extra.iter());
            let ext = if cmd == "nuc" { "tsv" } else { "fa" };
            assert_fails_clean(cmd, "interval", &tpl, bad, &Sink::File(ext), msg);
            ran += 1;
        }
    }
    // skipping instead: one warning with the count, and the rest are written
    for t in ["1", "4"] {
        let o = gloci()
            .args([
                "getfasta",
                "--fasta",
                &fx.fa,
                &absent,
                "--skip-invalid",
                "-t",
                t,
            ])
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(0));
        assert_eq!(
            String::from_utf8_lossy(&o.stdout).lines().count(),
            2,
            "-t {t}"
        );
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("skipped 1 interval"),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
    }
    // maskfasta
    let mask_past = s_of(&plain_file(
        dir,
        "mask-past.bed",
        "chr1\t10\t20\nchr1\t990\t1010\n",
    ));
    assert_fails_clean(
        "maskfasta",
        "mask past end",
        &["maskfasta", "--fasta", &fx.fa, "-b", "{IN}"],
        &mask_past,
        &Sink::File("fa"),
        "covers 990-1,010",
    );
    for (what, bad) in [("cut BGZF", &cut), ("BGZF without EOF block", &cut_noeof)] {
        assert_fails_clean(
            "maskfasta",
            what,
            &["maskfasta", "--fasta", "{IN}", "-b", &fx.masks],
            &s_of(bad),
            &Sink::File("fa"),
            "",
        );
        ran += 1;
    }
    assert!(ran >= 10);
}

/// Run gloci with its file size limited to `limit` bytes (writes past it fail with EFBIG, as on a
/// full disk, instead of killing the process with SIGXFSZ).
#[cfg(unix)]
fn p3a_run_size_limited(args: &[String], limit: u64) -> std::process::Output {
    use std::os::unix::process::CommandExt;
    let mut cmd = Command::new(bin());
    cmd.args(args);
    unsafe {
        cmd.pre_exec(move || {
            libc::signal(libc::SIGXFSZ, libc::SIG_IGN);
            let rlim = libc::rlimit {
                rlim_cur: limit,
                rlim_max: limit,
            };
            if libc::setrlimit(libc::RLIMIT_FSIZE, &rlim) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            Ok(())
        });
    }
    cmd.output().unwrap()
}

/// A spill file that cannot be written (the temp directory is "full": every file is capped at
/// 1 MiB) fails `sort`, `sort -n` and `intersect --sort` with exit 1 and a message naming the
/// directory, and leaves nothing behind: the spill directory is empty and no output exists.
#[cfg(unix)]
#[test]
fn p3a_full_tmp_dir_during_sort_leaves_nothing() {
    let d = tempdir();
    let (_sam, bam) = p3a_big_alignments(d.path(), 200_000, true);
    let tmp = d.path().join("spill");
    let out_dir = d.path().join("out");
    std::fs::create_dir(&tmp).unwrap();
    std::fs::create_dir(&out_dir).unwrap();
    let bed = s_of(&plain_file(d.path(), "t.bed", "chr1\t0\t5000000\n"));
    let cases: Vec<(&str, Vec<String>)> = vec![
        (
            "sort",
            vec![
                "sort".into(),
                bam.clone(),
                "-o".into(),
                s_of(&out_dir.join("o.bam")),
            ],
        ),
        (
            "sort -n",
            vec![
                "sort".into(),
                "-n".into(),
                bam.clone(),
                "-o".into(),
                s_of(&out_dir.join("o.bam")),
            ],
        ),
        (
            "sort sam",
            vec![
                "sort".into(),
                bam.clone(),
                "-O".into(),
                "sam".into(),
                "-o".into(),
                s_of(&out_dir.join("o.sam")),
            ],
        ),
        (
            "intersect --sort",
            vec![
                "intersect".into(),
                "-a".into(),
                bam.clone(),
                "-b".into(),
                bed.clone(),
                "--sort".into(),
                "-o".into(),
                s_of(&out_dir.join("o.bam")),
            ],
        ),
        (
            "coverage --sort",
            vec![
                "coverage".into(),
                "-a".into(),
                bed.clone(),
                "-b".into(),
                bam.clone(),
                "--sort".into(),
                "-o".into(),
                s_of(&out_dir.join("o.tsv")),
            ],
        ),
    ];
    let mut spilled_once = false;
    for (name, base) in cases {
        for t in if name == "sort" {
            vec!["1", "4"]
        } else {
            vec!["4"]
        } {
            let mut args = base.clone();
            args.extend(["--max-mem", "16M", "--tmp-dir", &s_of(&tmp), "-t", t].map(String::from));
            // unlimited: succeeds and removes its spill files
            let ok = gloci().args(&args).output().unwrap();
            assert_eq!(
                ok.status.code(),
                Some(0),
                "{name} -t {t}: {}",
                String::from_utf8_lossy(&ok.stderr)
            );
            assert_eq!(
                std::fs::read_dir(&tmp).unwrap().count(),
                0,
                "{name} -t {t}: spill files left after success"
            );
            for f in std::fs::read_dir(&out_dir).unwrap() {
                std::fs::remove_file(f.unwrap().path()).unwrap();
            }
            // capped: fails cleanly
            let o = p3a_run_size_limited(&args, 1 << 20);
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(1), "{name} -t {t}: {err}");
            assert!(err.starts_with("error: "), "{name} -t {t}: {err}");
            assert!(
                err.contains("File too large") || err.contains("No space left"),
                "{name} -t {t}: {err}"
            );
            assert_eq!(
                std::fs::read_dir(&tmp).unwrap().count(),
                0,
                "{name} -t {t}: spill files left: {err}"
            );
            assert_eq!(
                std::fs::read_dir(&out_dir).unwrap().count(),
                0,
                "{name} -t {t}: output left behind"
            );
            spilled_once = true;
        }
    }
    assert!(spilled_once);
}

/// `--chrom-map` renames as records are read: BAM and CRAM headers and records come out with the
/// new names, the inverse map takes them back, and getfasta, nuc and maskfasta find intervals
/// written with the other naming.
#[test]
fn p3a_chrom_map_round_trip() {
    let fx = P3a::new();
    let d = fx.path("cm");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let fwd = s_of(&plain_file(
        d,
        "fwd.tsv",
        "# to Ensembl\nchr1\t1\nchr2\t2\nchrM\tMT\n",
    ));
    let rev = s_of(&plain_file(d, "rev.tsv", "1\tchr1\n2\tchr2\nMT\tchrM\n"));
    let sam_names = |b: &[u8]| -> Vec<String> {
        String::from_utf8_lossy(&p3a_norm(b))
            .lines()
            .filter_map(|l| l.strip_prefix("@SQ\tSN:"))
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    };
    let ref_args = ["--reference", fx.fa.as_str()];
    for (kind, input, _) in fx.kinds() {
        if kind == "sam" {
            // SAM text written as it is cannot have its @SQ names rewritten: refused, not
            // half-renamed
            gloci()
                .args(["view", &input, "--chrom-map", &fwd])
                .assert()
                .code(2)
                .stderr(contains(
                    "--chrom-map cannot rename the reference names in SAM text",
                ));
            continue;
        }
        let to_bam = |src: &str, map: &str, name: &str, extra: &[&str]| -> String {
            let out = s_of(&d.join(name));
            gloci()
                .args(["view", src, "--chrom-map", map, "-o", &out])
                .args(extra)
                .assert()
                .success();
            out
        };
        let ex: Vec<&str> = if kind == "cram" {
            ref_args.to_vec()
        } else {
            vec![]
        };
        let renamed = to_bam(&input, &fwd, &format!("{kind}-ens.bam"), &ex);
        let hdr = |p: &str| {
            let o = gloci()
                .args(["view", p, "--header-only", "-O", "sam"])
                .output()
                .unwrap();
            sam_names(&o.stdout)
        };
        assert_eq!(hdr(&renamed), ["MT", "1", "2", "chr10", "chrX"], "{kind}");
        let recs = |p: &str| {
            p3a_records(
                &gloci()
                    .args(["view", p, "-O", "sam"])
                    .output()
                    .unwrap()
                    .stdout,
            )
        };
        let ens = recs(&renamed);
        assert!(
            ens.iter().any(|r| r.split('\t').nth(2) == Some("MT")),
            "{kind}"
        );
        assert!(
            ens.iter()
                .all(|r| !r.split('\t').nth(2).unwrap().starts_with("chr1\t")),
            "{kind}"
        );
        // and back: the original names, the same reads, in a BAM that view reads
        let back = to_bam(&renamed, &rev, &format!("{kind}-back.bam"), &[]);
        assert_eq!(
            hdr(&back),
            ["chrM", "chr1", "chr2", "chr10", "chrX"],
            "{kind}"
        );
        let orig = p3a_records(
            &gloci()
                .args(["view", &input, "-O", "sam"])
                .args(&ex)
                .output()
                .unwrap()
                .stdout,
        );
        assert_eq!(recs(&back), orig, "{kind}: round trip changed the reads");
        // other commands see the new names: coverage against a BED in the new names, genomecov
        // rows in the new names, equal to the original run under the inverse naming
        let ens_bed = s_of(&plain_file(
            d,
            "ens.bed",
            "MT\t0\t100\n1\t0\t300\n1\t600\t900\n2\t100\t900\n",
        ));
        let want = {
            let bed = s_of(&plain_file(
                d,
                "orig.bed",
                "chrM\t0\t100\nchr1\t0\t300\nchr1\t600\t900\nchr2\t100\t900\n",
            ));
            let o = gloci()
                .args(["coverage", "-a", &bed, "-b", &input])
                .args(&ex)
                .output()
                .unwrap();
            assert!(o.status.success());
            String::from_utf8_lossy(&o.stdout)
                .replace("chrM\t", "MT\t")
                .replace("chr1\t", "1\t")
                .replace("chr2\t", "2\t")
        };
        for t in ["1", "4"] {
            let o = gloci()
                .args([
                    "coverage",
                    "-a",
                    &ens_bed,
                    "-b",
                    &input,
                    "--chrom-map",
                    &fwd,
                    "-t",
                    t,
                ])
                .args(&ex)
                .output()
                .unwrap();
            // the map renames B's reads; A is already in the new names
            assert_eq!(
                o.status.code(),
                Some(0),
                "{kind}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            assert_eq!(String::from_utf8_lossy(&o.stdout), want, "{kind} -t {t}");
            let g = gloci()
                .args(["genomecov", &input, "--chrom-map", &fwd, "-t", t])
                .args(&ex)
                .output()
                .unwrap();
            let text = String::from_utf8_lossy(&g.stdout).into_owned();
            assert!(
                text.contains("MT\t0\t") && text.contains("\n1\t0\t"),
                "{kind}: {text}"
            );
            assert!(
                !text.contains("chr1\t") && !text.contains("chrM\t"),
                "{kind}: {text}"
            );
        }
    }
    // a map that would merge two contigs is refused
    let bad = s_of(&plain_file(d, "bad.tsv", "chr1\tchr2\n"));
    gloci()
        .args(["view", &fx.bam, "--chrom-map", &bad])
        .assert()
        .code(2)
        .stderr(contains("which would merge them"));

    // FASTA: intervals in the other naming find their sequences through the map, and the
    // output uses the FASTA's names
    let plain = {
        let bed = s_of(&plain_file(
            d,
            "plain.bed",
            "chr1\t10\t40\tx\t0\t+\nchr2\t100\t130\ty\t0\t-\nchrM\t5\t25\tz\t0\t+\n",
        ));
        let o = gloci()
            .args(["getfasta", "--fasta", &fx.fa, &bed, "-s", "--name"])
            .output()
            .unwrap();
        assert!(o.status.success());
        o.stdout
    };
    let ens_iv = s_of(&plain_file(
        d,
        "ens-iv.bed",
        "1\t10\t40\tx\t0\t+\n2\t100\t130\ty\t0\t-\nMT\t5\t25\tz\t0\t+\n",
    ));
    let to_fa_names = s_of(&plain_file(d, "to-fa.tsv", "1\tchr1\n2\tchr2\nMT\tchrM\n"));
    for t in ["1", "4"] {
        let o = gloci()
            .args([
                "getfasta",
                "--fasta",
                &fx.fa,
                &ens_iv,
                "-s",
                "--name",
                "--chrom-map",
                &to_fa_names,
                "-t",
                t,
            ])
            .output()
            .unwrap();
        assert_eq!(
            o.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert_eq!(o.stdout, plain, "getfasta -t {t}");
        // without the map the names are not in the FASTA
        gloci()
            .args(["getfasta", "--fasta", &fx.fa, &ens_iv, "-t", t])
            .assert()
            .code(1)
            .stderr(contains("is not in the FASTA index"));
    }
    let nuc_plain = {
        let bed = s_of(&plain_file(
            d,
            "nplain.bed",
            "chr1\t10\t40\nchr2\t100\t130\n",
        ));
        gloci()
            .args(["nuc", "--fasta", &fx.fa, &bed])
            .output()
            .unwrap()
            .stdout
    };
    let nuc_map = {
        let bed = s_of(&plain_file(d, "nens.bed", "1\t10\t40\n2\t100\t130\n"));
        gloci()
            .args(["nuc", "--fasta", &fx.fa, &bed, "--chrom-map", &to_fa_names])
            .output()
            .unwrap()
            .stdout
    };
    assert_eq!(nuc_plain, nuc_map, "nuc through a map");
    let masks_ens = s_of(&plain_file(
        d,
        "masks-ens.bed",
        "MT\t5\t20\n1\t0\t30\n1\t100\t140\n2\t500\t700\n10\t10\t20\n",
    ));
    let to_fa_names10 = s_of(&plain_file(
        d,
        "to-fa10.tsv",
        "1\tchr1\n2\tchr2\nMT\tchrM\n10\tchr10\n",
    ));
    let want = gloci()
        .args(["maskfasta", "--fasta", &fx.fa, "-b", &fx.masks])
        .output()
        .unwrap()
        .stdout;
    let got = gloci()
        .args([
            "maskfasta",
            "--fasta",
            &fx.fa,
            "-b",
            &masks_ens,
            "--chrom-map",
            &to_fa_names10,
        ])
        .output()
        .unwrap();
    assert_eq!(
        got.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&got.stderr)
    );
    assert_eq!(got.stdout, want, "maskfasta through a map");
}
