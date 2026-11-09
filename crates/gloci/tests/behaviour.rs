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
