//! `gloci stats jaccard`: the merged-set table, strands, fractions, empty
//! inputs, points, VCF and alignment inputs, order errors, `--sort` and
//! determinism across thread counts.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

const HEADER: &str = "intersection\tunion\tjaccard\tn_intersections\n";

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

struct Fx {
    dir: tempfile::TempDir,
}

impl Fx {
    fn new() -> Fx {
        Fx { dir: tempdir() }
    }
    fn file(&self, name: &str, text: &str) -> String {
        s(&plain_file(self.dir.path(), name, text))
    }
    /// `stats jaccard -a A -b B extra...` over BED texts: stdout.
    fn run(&self, a: &str, b: &str, extra: &[&str]) -> String {
        let af = self.file("a.bed", a);
        let bf = self.file("b.bed", b);
        let o = gloci()
            .args(["stats", "jaccard", "-a", &af, "-b", &bf])
            .args(extra)
            .assert()
            .success();
        String::from_utf8(o.get_output().stdout.clone()).unwrap()
    }
}

fn row(r: &str) -> String {
    format!("{HEADER}{r}\n")
}

#[test]
fn jaccard_table() {
    let fx = Fx::new();
    // A merges to 10-30 and 50-60 (bookended 10-20 and 20-30 join); unmerged A would count 3
    // intersections.
    let got = fx.run(
        "chr1\t10\t20\nchr1\t20\t30\nchr1\t50\t60\n",
        "chr1\t15\t25\nchr1\t55\t70\n",
        &[],
    );
    assert_eq!(got, row("15\t40\t0.375\t2"));
}

#[test]
fn jaccard_strand() {
    let fx = Fx::new();
    let a = "chr1\t10\t20\t.\t0\t+\nchr1\t30\t40\t.\t0\t-\n";
    let b = "chr1\t15\t35\t.\t0\t+\n";
    assert_eq!(fx.run(a, b, &[]), row("10\t30\t0.333333\t2"));
    assert_eq!(fx.run(a, b, &["-s"]), row("5\t35\t0.142857\t1"));
    assert_eq!(fx.run(a, b, &["-S"]), row("5\t35\t0.142857\t1"));
    // + and - clusters stay apart under -s even when they overlap each other.
    let a2 = "chr1\t10\t20\t.\t0\t+\nchr1\t15\t30\t.\t0\t-\n";
    let b2 = "chr1\t0\t100\t.\t0\t+\n";
    assert_eq!(fx.run(a2, b2, &["-s"]), row("10\t115\t0.0869565\t1"));
    assert_eq!(fx.run(a2, b2, &[]), row("20\t100\t0.2\t1"));
}

#[test]
fn jaccard_unknown_strand_warns_once() {
    let fx = Fx::new();
    let af = fx.file("a.bed", "chr1\t10\t20\nchr1\t30\t40\n");
    let bf = fx.file("b.bed", "chr1\t10\t40\t.\t0\t+\n");
    let o = gloci()
        .args(["stats", "jaccard", "-s", "-a", &af, "-b", &bf])
        .assert()
        .success()
        .stdout(row("0\t50\t0\t0"));
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    assert_eq!(err.matches("unknown strand").count(), 1, "{err}");
}

#[test]
fn jaccard_fraction() {
    let fx = Fx::new();
    // Merged A is 0-20 (two overlapping records); the 5-base overlap with B is a quarter of it.
    let a = "chr1\t0\t12\nchr1\t8\t20\n";
    let b = "chr1\t15\t30\n";
    assert_eq!(fx.run(a, b, &[]), row("5\t30\t0.166667\t1"));
    assert_eq!(fx.run(a, b, &["-f", "0.5"]), row("0\t35\t0\t0"));
    assert_eq!(fx.run(a, b, &["-f", "0.25"]), row("5\t30\t0.166667\t1"));
    assert_eq!(fx.run(a, b, &["-f", "0.3", "-r"]), row("0\t35\t0\t0"));
    assert_eq!(fx.run(a, b, &["-F", "0.3"]), row("5\t30\t0.166667\t1"));
    assert_eq!(fx.run(a, b, &["--min-bases", "6"]), row("0\t35\t0\t0"));
}

#[test]
fn jaccard_empty_and_nan() {
    let fx = Fx::new();
    assert_eq!(fx.run("", "", &[]), row("0\t0\tnan\t0"));
    assert_eq!(fx.run("chr1\t10\t20\n", "", &[]), row("0\t10\t0\t0"));
    assert_eq!(
        fx.run("chr1\t10\t20\n", "chr1\t10\t20\n", &[]),
        row("10\t10\t1\t1")
    );
}

#[test]
fn jaccard_single_precision() {
    let fx = Fx::new();
    // 3 / 16777217 divides in f32: 16777217 is not representable.
    assert_eq!(
        fx.run("chr1\t0\t16777217\n", "chr1\t0\t3\n", &[]),
        row("3\t16777217\t1.78814e-07\t1")
    );
}

#[test]
fn jaccard_points_count_as_intersections() {
    let fx = Fx::new();
    // A point covers no bases (bedtools widens it to two); inside B it is one intersection.
    assert_eq!(
        fx.run("chr1\t25\t25\n", "chr1\t10\t30\n", &[]),
        row("0\t20\t0\t1")
    );
    assert_eq!(
        fx.run("chr1\t25\t25\n", "chr1\t25\t25\n", &[]),
        row("0\t0\tnan\t1")
    );
    // On the edge only with --inclusive-points.
    assert_eq!(
        fx.run("chr1\t30\t30\n", "chr1\t10\t30\n", &[]),
        row("0\t20\t0\t0")
    );
    assert_eq!(
        fx.run("chr1\t30\t30\n", "chr1\t10\t30\n", &["--inclusive-points"]),
        row("0\t20\t0\t1")
    );
}

#[test]
fn jaccard_vcf_uses_footprints() {
    let fx = Fx::new();
    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t11\t.\tACG\tA\t.\t.\t.\n";
    let a = fx.file("a.vcf", vcf);
    let b = fx.file("b.bed", "chr1\t12\t20\n");
    // REF span is 10-13: one shared base... two (12-13 is one base), union 3 + 8 - 1.
    gloci()
        .args(["stats", "jaccard", "-a", &a, "-b", &b, "--vcf-span", "ref"])
        .assert()
        .success()
        .stdout(row("1\t10\t0.1\t1"));
}

#[test]
fn jaccard_alignments() {
    let fx = Fx::new();
    let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
r1\t0\tchr1\t11\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n";
    let a = fx.file("a.sam", sam);
    let b = fx.file("b.bed", "chr1\t15\t40\n");
    gloci()
        .args(["stats", "jaccard", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout(row("5\t30\t0.166667\t1"));
}

#[test]
fn jaccard_unsorted_input_is_refused_and_sort_fixes_it() {
    let fx = Fx::new();
    let a = "chr1\t50\t60\nchr1\t10\t20\n";
    let b = "chr1\t15\t25\nchr1\t55\t70\n";
    let af = fx.file("a.bed", a);
    let bf = fx.file("b.bed", b);
    gloci()
        .args(["stats", "jaccard", "-a", &af, "-b", &bf])
        .assert()
        .code(1)
        .stderr(contains("not sorted"))
        .stderr(contains("--sort"));
    let sorted = fx.run("chr1\t10\t20\nchr1\t50\t60\n", b, &[]);
    assert_eq!(fx.run(a, b, &["--sort"]), sorted);
    assert_eq!(sorted, row("10\t35\t0.285714\t2"));
    // B out of order is refused too.
    let af2 = fx.file("a2.bed", "chr1\t10\t20\n");
    let bf2 = fx.file("b2.bed", "chr1\t55\t70\nchr1\t15\t25\n");
    gloci()
        .args(["stats", "jaccard", "-a", &af2, "-b", &bf2])
        .assert()
        .code(1)
        .stderr(contains("not sorted"));
}

#[test]
fn jaccard_usage_errors() {
    let fx = Fx::new();
    let af = fx.file("a.bed", "chr1\t10\t20\n");
    gloci()
        .args(["stats", "jaccard", "-a", &af])
        .assert()
        .code(2);
    gloci()
        .args(["stats", "jaccard", "-a", &af, "-b", &af, "--split"])
        .assert()
        .code(2);
    gloci()
        .args(["stats", "jaccard", "-a", &af, "-b", "/nonexistent/b.bed"])
        .assert()
        .failure()
        .stderr(contains("/nonexistent/b.bed"));
}

#[test]
fn jaccard_is_identical_at_every_thread_count() {
    let fx = Fx::new();
    let mut a = String::new();
    let mut b = String::new();
    for c in ["chr1", "chr2", "chr10", "chrX"] {
        for i in 0..400u64 {
            a += &format!(
                "{c}\t{}\t{}\t.\t0\t{}\n",
                i * 37,
                i * 37 + 50,
                ['+', '-'][(i % 2) as usize]
            );
            b += &format!(
                "{c}\t{}\t{}\t.\t0\t{}\n",
                i * 41,
                i * 41 + 30,
                ['+', '-'][((i / 3) % 2) as usize]
            );
        }
    }
    let af = fx.file("a.bed", &a);
    let bf = fx.file("b.bed", &b);
    for extra in [&[][..], &["-s"][..], &["-f", "0.3", "-r"][..]] {
        let outs: Vec<String> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = gloci()
                    .args(["stats", "jaccard", "-a", &af, "-b", &bf, "-t", t])
                    .args(extra)
                    .assert()
                    .success();
                String::from_utf8(o.get_output().stdout.clone()).unwrap()
            })
            .collect();
        assert_eq!(outs[0], outs[1]);
        assert_eq!(outs[0], outs[2]);
        assert!(outs[0].starts_with(HEADER));
    }
}

#[test]
fn jaccard_writes_to_a_file() {
    let fx = Fx::new();
    let af = fx.file("a.bed", "chr1\t10\t20\n");
    let out = s(&fx.dir.path().join("j.tsv"));
    gloci()
        .args(["stats", "jaccard", "-a", &af, "-b", &af, "-o", &out])
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&out).unwrap(), row("10\t10\t1\t1"));
}
