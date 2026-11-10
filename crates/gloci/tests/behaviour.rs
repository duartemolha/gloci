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
/// byte-identical output.
fn same_at_every_thread_count(what: &str, args: &[&str]) -> Vec<u8> {
    let outs: Vec<Vec<u8>> = ["1", "2", "8"]
        .iter()
        .map(|t| stdout_with(args, t))
        .collect();
    assert_eq!(outs[0], outs[1], "{what}: -t 1 vs -t 2");
    assert_eq!(outs[1], outs[2], "{what}: -t 2 vs -t 8");
    outs[0].clone()
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
    let bed_cmds: Vec<(&str, Vec<&str>)> = vec![
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
    ];
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
    // --sort over a budget smaller than the input spills to disk and gives
    // the bytes of the in-memory sort.
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
    for (name, tpl) in &vcf_cmds {
        let mut outs = Vec::new();
        // cluster and spacing append a column, which a BCF record has no place for
        let takes_bcf = !matches!(*name, "cluster" | "spacing");
        for v in [&vcf_plain, &vcf_gz, &vcf_idx, &vcf_bcf] {
            if !takes_bcf && *v == vcf_bcf {
                continue;
            }
            // concat needs the second file to start where the first ends
            let second = fx.plain(
                "v_tail.vcf",
                &p2_vcf(4..=4, 1).replace("chr4\t100", "chr4\t99999"),
            );
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
        // concat writes in the format of its first input, so only its
        // thread determinism is compared (above)
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
    let tail = v_gz_tail(d.path());
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

/// A one-record VCF that follows the 20 000-record one at the same
/// chromosome end (for `vcf concat`).
fn v_gz_tail(dir: &Path) -> String {
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
    gloci()
        .args(["vcf", "extract", &hdr_only_bcf, "--info", "DP"])
        .assert()
        .success();

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
            if *kind == "unsorted" && matches!(*name, "shift" | "groupby") {
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
    // extract and stats read an unsorted VCF without complaint (their output
    // does not depend on the order), so only the sorted-input commands fail
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
                continue;
            }
            assert_fails_clean(name, kind, tpl, bad, sink, msg);
            ran += 1;
        }
        assert_missing_dir(name, tpl, gv, sink);
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
