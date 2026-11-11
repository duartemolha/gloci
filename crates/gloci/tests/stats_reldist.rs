//! `gloci stats reldist`: the bin table, `--detail`, the records bedtools skips, unsorted and
//! `--sort` runs, order errors, midpoints past the `int` range and determinism across thread
//! counts. The bedtools parity is in the golden cases (`stats-reldist-*-vs-bedtools`).
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

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
    /// `stats reldist -a A -b B extra...` over BED texts: stdout.
    fn run(&self, a: &str, b: &str, extra: &[&str]) -> String {
        let af = self.file("a.bed", a);
        let bf = self.file("b.bed", b);
        let o = gloci()
            .args(["stats", "reldist", "-a", &af, "-b", &bf])
            .args(extra)
            .assert()
            .success();
        String::from_utf8(o.get_output().stdout.clone()).unwrap()
    }
}

/// B midpoints 100, 200, 300.
const B: &str = "chr1\t90\t110\nchr1\t190\t210\nchr1\t290\t310\n";

/// One A interval of midpoint `m`.
fn at(chrom: &str, m: u64) -> String {
    format!("{chrom}\t{}\t{}\n", m - 1, m + 1)
}

/// A midpoints 150, 129, 100, 300, 50, 400, in this order.
fn a_unsorted() -> String {
    [150, 129, 100, 300, 50, 400]
        .iter()
        .map(|&m| at("chr1", m))
        .collect()
}

fn a_sorted() -> String {
    [50, 100, 129, 150, 300, 400]
        .iter()
        .map(|&m| at("chr1", m))
        .collect()
}

#[test]
fn reldist_table() {
    let got = Fx::new().run(&a_unsorted(), B, &["--unsorted"]);
    // 129 is 29/100 = 0.29, but 0.29 * 100 is 28.999..., so it lands in bin 0.28 (as in
    // bedtools). 50 and 400 lie outside the B midpoints and are skipped.
    let want = "\
reldist\tcount\ttotal\tfraction
0.00\t2\t4\t0.500
0.28\t1\t4\t0.250
0.50\t1\t4\t0.250
";
    assert_eq!(got, want);
}

#[test]
fn reldist_detail() {
    let got = Fx::new().run(&a_unsorted(), B, &["--unsorted", "--detail"]);
    // In A's order; the 129 record prints its value, 0.29, to three places.
    let want = "\
chr1\t149\t151\t0.500
chr1\t128\t130\t0.290
chr1\t99\t101\t0.000
chr1\t299\t301\t0.000
";
    assert_eq!(got, want);
}

#[test]
fn reldist_unsorted_equals_sorted() {
    let fx = Fx::new();
    let sorted = fx.run(&a_sorted(), B, &[]);
    assert_eq!(sorted, fx.run(&a_unsorted(), B, &["--unsorted"]));
    assert_eq!(sorted, fx.run(&a_unsorted(), B, &["--sort"]));
    assert_eq!(sorted, fx.run(&a_sorted(), B, &["--unsorted"]));
}

#[test]
fn reldist_detail_with_sort_is_in_sorted_order() {
    let got = Fx::new().run(&a_unsorted(), B, &["--sort", "--detail"]);
    assert_eq!(got.lines().next(), Some("chr1\t99\t101\t0.000"));
    assert_eq!(got.lines().count(), 4);
}

#[test]
fn reldist_skips_a_without_b_and_prints_only_the_header_when_nothing_qualifies() {
    let fx = Fx::new();
    let a = format!("{}{}", at("chr1", 150), at("chr2", 150));
    let r = fx.run(&a, B, &[]);
    assert_eq!(r, "reldist\tcount\ttotal\tfraction\n0.50\t1\t1\t1.000\n");
    let r = fx.run(&at("chr2", 150), B, &[]);
    assert_eq!(r, "reldist\tcount\ttotal\tfraction\n");
    assert_eq!(fx.run(&at("chr2", 150), B, &["--detail"]), "");
    // An empty A, an empty B.
    assert_eq!(fx.run("", B, &[]), "reldist\tcount\ttotal\tfraction\n");
    assert_eq!(
        fx.run(&at("chr1", 150), "", &[]),
        "reldist\tcount\ttotal\tfraction\n"
    );
}

#[test]
fn reldist_a_chromosome_with_one_b_midpoint_is_skipped() {
    // bedtools leaves out every A on a chromosome with fewer than two B midpoints, even one
    // that sits on the midpoint.
    let r = Fx::new().run(&at("chr1", 100), "chr1\t99\t101\n", &[]);
    assert_eq!(r, "reldist\tcount\ttotal\tfraction\n");
}

#[test]
fn reldist_chromosomes_in_step() {
    // chr1 has no B, chr2 and chr3 do; A has chr1, chr2, chr3 and B skips chr1.
    let fx = Fx::new();
    let a = format!(
        "{}{}{}{}",
        at("chr1", 5),
        at("chr2", 150),
        at("chr3", 100),
        at("chr3", 200)
    );
    let b =
        "chr2\t90\t110\nchr2\t190\t210\nchr3\t90\t110\nchr3\t190\t210\nchr4\t1\t3\nchr4\t5\t7\n";
    let r = fx.run(&a, b, &[]);
    assert_eq!(
        r,
        "reldist\tcount\ttotal\tfraction\n0.00\t2\t3\t0.667\n0.50\t1\t3\t0.333\n"
    );
    assert_eq!(r, fx.run(&a, b, &["--unsorted"]));
}

#[test]
fn reldist_midpoint_past_int_range() {
    // bedtools adds B's ends in a 32-bit int, which wraps past 2^31; gloci does not
    // (differences.md#reldist-midpoint-overflow). B midpoints 150, 1500000050, 1500000250.
    let b = "chr1\t100\t200\nchr1\t1500000000\t1500000100\nchr1\t1500000200\t1500000300\n";
    let a = "chr1\t1500000000\t1500000100\nchr1\t1500000050\t1500000150\n";
    let got = Fx::new().run(a, b, &["--detail"]);
    assert_eq!(
        got,
        "chr1\t1500000000\t1500000100\t0.000\nchr1\t1500000050\t1500000150\t0.250\n"
    );
}

#[test]
fn reldist_unsorted_input_names_the_flags() {
    let fx = Fx::new();
    let a = fx.file("a.bed", &a_unsorted());
    let b = fx.file("b.bed", B);
    gloci()
        .args(["stats", "reldist", "-a", &a, "-b", &b])
        .assert()
        .failure()
        .stderr(contains("not sorted"))
        .stderr(contains("--sort"))
        .stderr(contains("--unsorted"));
    // B out of order is found too.
    let a = fx.file("a2.bed", &a_sorted());
    let b = fx.file("b2.bed", "chr1\t290\t310\nchr1\t90\t110\n");
    gloci()
        .args(["stats", "reldist", "-a", &a, "-b", &b])
        .assert()
        .failure()
        .stderr(contains("b2.bed is not sorted"));
}

#[test]
fn reldist_sort_and_unsorted_conflict() {
    let fx = Fx::new();
    let a = fx.file("a.bed", &a_sorted());
    let b = fx.file("b.bed", B);
    gloci()
        .args([
            "stats",
            "reldist",
            "-a",
            &a,
            "-b",
            &b,
            "--sort",
            "--unsorted",
        ])
        .assert()
        .failure();
}

#[test]
fn reldist_is_identical_at_every_thread_count() {
    let fx = Fx::new();
    let a = a_sorted();
    let want = fx.run(&a, B, &["--detail"]);
    for t in ["1", "2", "8"] {
        assert_eq!(fx.run(&a, B, &["-t", t, "--detail"]), want, "-t {t}");
        assert_eq!(
            fx.run(&a_unsorted(), B, &["-t", t, "--sort"]),
            fx.run(&a, B, &[]),
            "-t {t}"
        );
    }
}

#[test]
fn reldist_detail_of_a_vcf_keeps_its_header_columns() {
    let fx = Fx::new();
    let vcf = fx.file(
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t150\t.\tA\tC\t.\t.\t.\n",
    );
    let b = fx.file("b.bed", B);
    let o = gloci()
        .args(["stats", "reldist", "--detail", "-a", &vcf, "-b", &b])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(
        out.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\treldist\n"),
        "{out}"
    );
    assert!(
        out.ends_with("chr1\t150\t.\tA\tC\t.\t.\t.\t0.490\n"),
        "{out}"
    );
}

#[test]
fn reldist_midpoints_above_2_pow_63_do_not_wrap() {
    // B midpoints 1.0e19 and 1.8e19 (past i64::MAX); A at 1.2e19 is 0.25 of the way.
    let b = "chr1\t9999999999999999990\t10000000000000000010\nchr1\t17999999999999999990\t18000000000000000010\n";
    let a = "chr1\t11999999999999999990\t12000000000000000010\n";
    let got = Fx::new().run(a, b, &["--detail"]);
    assert_eq!(
        got,
        "chr1\t11999999999999999990\t12000000000000000010\t0.250\n"
    );
}
