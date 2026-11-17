//! `gloci depth --summary` (D173): samtools coverage 1.21's table. The expected rows are
//! samtools'; the live comparison is in `xtask/diff/cases.toml`.
mod common;
use common::*;
use predicates::str::contains;
use std::path::PathBuf;

const H: &str =
    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq";

fn f(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/samtools")
        .join(name)
        .display()
        .to_string()
}

fn summary(args: &[&str]) -> String {
    let o = gloci()
        .arg("depth")
        .args(args)
        .arg("--summary")
        .assert()
        .success()
        .get_output()
        .clone();
    String::from_utf8(o.stdout).unwrap()
}

/// `H` and then `rows`, each ended by a newline.
fn table(rows: &[&str]) -> String {
    let mut t = format!("{H}\n");
    for r in rows {
        t.push_str(r);
        t.push('\n');
    }
    t
}

const DEPTH_ROWS: [&str; 3] = [
    "chr1\t1\t500\t15\t190\t38\t0.546\t47.4\t55.9",
    "chr2\t1\t300\t3\t60\t20\t0.2\t40\t60",
    "chr3\t1\t200\t0\t0\t0\t0\t0\t0",
];

#[test]
fn summary_of_depth_bam() {
    assert_eq!(summary(&[&f("depth.bam")]), table(&DEPTH_ROWS));
}

#[test]
fn summary_sums_several_inputs() {
    assert_eq!(
        summary(&[&f("depth.bam"), &f("depth2.bam")]),
        table(&[
            "chr1\t1\t500\t19\t247\t49.4\t0.736\t45.5\t56.8",
            "chr2\t1\t300\t4\t90\t30\t0.3\t40\t60",
            "chr3\t1\t200\t1\t10\t5\t0.05\t40\t60",
        ])
    );
}

#[test]
fn summary_lists_uncovered_references_last() {
    assert_eq!(
        summary(&[&f("mosdepth.bam")]),
        table(&[
            "chr1\t1\t1000\t19\t530\t53\t0.76\t0\t53.9",
            "chr4\t1\t20000\t322\t31\t0.155\t0.018\t0\t60",
            "chr2\t1\t300\t0\t0\t0\t0\t0\t0",
            "chr3\t1\t200\t0\t0\t0\t0\t0\t0",
        ])
    );
    assert_eq!(
        summary(&[&f("mosdepth.bam"), "--exclude-flags", "0"]),
        table(&[
            "chr1\t1\t1000\t23\t549\t54.9\t0.82\t0\t52.4",
            "chr2\t1\t300\t2\t20\t6.66667\t0.0666667\t0\t30",
            "chr4\t1\t20000\t322\t31\t0.155\t0.018\t0\t60",
            "chr3\t1\t200\t0\t0\t0\t0\t0\t0",
        ])
    );
}

#[test]
fn summary_region_window_as_given() {
    let b = f("depth.bam");
    for (region, row) in [
        (
            "chr1:100-200",
            "chr1\t100\t200\t1\t3\t2.9703\t0.029703\t40\t60",
        ),
        ("chr1:600-700", "chr1\t600\t700\t0\t0\t0\t0\t0\t0"),
        ("chr1:450", "chr1\t450\t500\t0\t0\t0\t0\t0\t0"),
        ("chr3", "chr3\t1\t200\t0\t0\t0\t0\t0\t0"),
    ] {
        assert_eq!(
            summary(&[&b, "--region", region]),
            table(&[row]),
            "{region}"
        );
    }
}

#[test]
fn summary_min_read_length_counts_soft_clips() {
    assert_eq!(
        summary(&[&f("depth.bam"), "--min-read-length", "50"]),
        table(&[
            "chr1\t1\t500\t2\t80\t16\t0.16\t40\t60",
            "chr2\t1\t300\t0\t0\t0\t0\t0\t0",
            "chr3\t1\t200\t0\t0\t0\t0\t0\t0",
        ])
    );
}

#[test]
fn summary_seq_star_reads_count_depth_not_quality() {
    let b = f("depth-noseq.bam");
    let chr1 = |args: &[&str]| {
        summary(&[&[b.as_str()], args].concat())
            .lines()
            .nth(1)
            .unwrap()
            .to_string()
    };
    assert_eq!(chr1(&[]), "chr1\t1\t500\t3\t89\t17.8\t0.238\t40\t60");
    assert_eq!(
        chr1(&["--min-read-length", "50"]),
        "chr1\t1\t500\t1\t50\t10\t0.1\t0\t60"
    );
}

#[test]
fn summary_qual_star_counts_255() {
    let t = tempdir();
    let sam = std::fs::read_to_string(f("depth.sam")).unwrap();
    let mut text = String::new();
    for l in sam.lines() {
        if l.starts_with('@') || l.starts_with("p1\t") || l.starts_with("nq\t") {
            text.push_str(l);
            text.push('\n');
        }
    }
    let p = plain_file(t.path(), "nq.sam", &text);
    let got = summary(&[&p.display().to_string()]);
    assert_eq!(
        got.lines().nth(1).unwrap(),
        "chr1\t1\t500\t3\t35\t7\t0.1\t83\t60"
    );
}

#[test]
fn summary_filters() {
    let row = |args: &[&str], n: usize| summary(args).lines().nth(n).unwrap().to_string();
    let b = f("depth.bam");
    assert_eq!(
        row(&[&b, "--min-mapq", "30", "--min-baseq", "30"], 1),
        "chr1\t1\t500\t14\t189\t37.8\t0.524\t47.9\t57.9"
    );
    assert_eq!(
        row(&[&b, "--exclude-flags", "0"], 1),
        "chr1\t1\t500\t20\t220\t44\t0.606\t46.7\t51.2"
    );
    let p = f("depth-pastend.bam");
    assert_eq!(row(&[&p], 1), "chr1\t1\t500\t1\t20\t4\t0.04\t40\t60");
    assert_eq!(
        row(&[&p], 2),
        "chr2\t1\t300\t3\t31\t10.3333\t0.106667\t40\t60"
    );
}

#[test]
fn summary_cram_reads_qualities() {
    // samtools coverage prints 255 for CRAM base qualities unless -Q is given.
    assert_eq!(summary(&[&f("depth.cram")]), table(&DEPTH_ROWS));
}

#[test]
fn summary_usage_errors() {
    let b = f("depth.bam");
    gloci()
        .args([
            "depth",
            &b,
            "--summary",
            "--region",
            "chr1:1-5",
            "--region",
            "chr2",
        ])
        .assert()
        .code(2)
        .stderr(contains("one --region"));
    gloci()
        .args(["depth", &b, "--summary", "--region", "chrZ:1-10"])
        .assert()
        .code(2)
        .stderr(contains("chrZ"));
    for extra in [
        &["--all-positions"][..],
        &["--all-references"],
        &["--count-deletions"],
        &["--skip-mate-overlap"],
        &["--column-names"],
        &["--regions-file", "x.bed"],
    ] {
        gloci()
            .args(["depth", &b, "--summary"])
            .args(extra)
            .assert()
            .code(2);
    }
}

#[test]
fn summary_threads_identical() {
    let run = |t: &str| summary(&[&f("depth.bam"), &f("depth2.bam"), "-t", t]);
    assert_eq!(run("1"), run("4"));
}

#[test]
fn summary_writes_o_file() {
    let t = tempdir();
    let b = f("depth.bam");
    let want = summary(&[&b]);
    let p = t.path().join("t.tsv");
    gloci()
        .args(["depth", &b, "--summary", "-o"])
        .arg(&p)
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&p).unwrap(), want);
    let gz = t.path().join("t.tsv.gz");
    gloci()
        .args(["depth", &b, "--summary", "-o"])
        .arg(&gz)
        .assert()
        .success();
    let bytes = std::fs::read(&gz).unwrap();
    // BGZF: a gzip member with the BC extra field.
    assert_eq!(&bytes[..4], [0x1f, 0x8b, 8, 4]);
    assert_eq!(&bytes[12..14], b"BC");
}

#[test]
fn summary_huge_header() {
    let t = tempdir();
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for i in 0..200_000 {
        text.push_str(&format!("@SQ\tSN:c{i}\tLN:10\n"));
    }
    text.push_str("r\t0\tc199999\t1\t60\t5M\t*\t0\t0\t*\t*\n");
    let p = plain_file(t.path(), "huge.sam", &text);
    let start = std::time::Instant::now();
    let got = summary(&[&p.display().to_string()]);
    assert!(start.elapsed().as_secs() < 20, "took {:?}", start.elapsed());
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines.len(), 200_001);
    assert_eq!(lines[1], "c199999\t1\t10\t1\t5\t50\t0.5\t0\t60");
    assert_eq!(lines[2], "c0\t1\t10\t0\t0\t0\t0\t0\t0");
    assert_eq!(lines[200_000], "c199998\t1\t10\t0\t0\t0\t0\t0\t0");
}
