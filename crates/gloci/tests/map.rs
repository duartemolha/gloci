//! `gloci map`: aggregations over the overlapping B, windows, fractions,
//! VCF on either side, skip-unmatched, usage and data errors, determinism.
mod common;
use common::*;
use predicates::prelude::*;
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
    /// `map -a a -b b extra...` over two BED texts, as output lines.
    fn run(&self, a: &str, b: &str, extra: &[&str]) -> Vec<String> {
        let (af, bf) = (self.file("a.bed", a), self.file("b.bed", b));
        let mut v = vec!["map", "-a", &af, "-b", &bf];
        v.extend_from_slice(extra);
        let o = gloci().args(&v).assert().success();
        String::from_utf8(o.get_output().stdout.clone())
            .unwrap()
            .lines()
            .map(String::from)
            .collect()
    }
}

fn code(args: &[&str]) -> i32 {
    gloci().args(args).output().unwrap().status.code().unwrap()
}

const A: &str = "chr1\t100\t200\ta1\t0\t+\nchr1\t300\t400\ta2\t0\t-\nchr2\t5\t15\ta3\t0\t+\n";
const B: &str = "chr1\t150\t160\tb1\t2.5\t+\nchr1\t180\t260\tb2\t4\t-\nchr1\t395\t500\tb3\t.\t+\n";

#[test]
fn sum_count() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(A, B, &["-c", "5:sum", "-c", "0:count"]),
        [
            "chr1\t100\t200\ta1\t0\t+\t6.5\t2",
            "chr1\t300\t400\ta2\t0\t-\t.\t1",
            "chr2\t5\t15\ta3\t0\t+\t.\t0",
        ]
    );
    // commas separate too; integers print as integers
    let ints = "chr1\t150\t160\tb1\t2\t+\nchr1\t180\t260\tb2\t4\t-\n";
    assert_eq!(
        fx.run(
            A,
            ints,
            &["-c", "5:sum,5:min,5:max,5:mean,5:median,5:collapse"]
        )[0],
        "chr1\t100\t200\ta1\t0\t+\t6\t2\t4\t3\t3\t2,4"
    );
}

#[test]
fn no_overlap() {
    let fx = Fx::new();
    let rows = fx.run(
        "chr1\t0\t5\ta\t0\t+\n",
        B,
        &[
            "-c",
            "5:sum",
            "-c",
            "0:count",
            "-c",
            "5:count_distinct",
            "-c",
            "5:stdev",
            "-c",
            "0:bases",
            "-c",
            "0:indicator",
            "-c",
            "4:collapse",
        ],
    );
    assert_eq!(rows, ["chr1\t0\t5\ta\t0\t+\t.\t0\t0\t.\t0\t0\t."]);
    let rows = fx.run(
        "chr1\t0\t5\ta\t0\t+\n",
        B,
        &["-c", "5:sum", "--empty", "NA"],
    );
    assert_eq!(rows, ["chr1\t0\t5\ta\t0\t+\tNA"]);
}

#[test]
fn window_range() {
    let fx = Fx::new();
    // a2 (300-400) reaches b2 (180-260) only through the window, which adds
    // no bases to `bases`
    assert_eq!(
        fx.run(A, B, &["-w", "100", "-c", "0:count", "-c", "0:bases"]),
        [
            "chr1\t100\t200\ta1\t0\t+\t2\t30",
            "chr1\t300\t400\ta2\t0\t-\t2\t5",
            "chr2\t5\t15\ta3\t0\t+\t0\t0",
        ]
    );
    // --left/--right pad one side
    assert_eq!(
        fx.run(
            A,
            "chr1\t210\t220\tx\t1\t+\n",
            &["--left", "11", "-c", "0:count"]
        )[0],
        "chr1\t100\t200\ta1\t0\t+\t1"
    );
    assert_eq!(
        fx.run(
            A,
            "chr1\t210\t220\tx\t1\t+\n",
            &["--right", "10", "-c", "0:count"]
        )[0],
        "chr1\t100\t200\ta1\t0\t+\t0"
    );
}

#[test]
fn fractions() {
    let fx = Fx::new();
    // b1 covers 10 of 100 bases of a1; b2 covers 20
    assert_eq!(
        fx.run(A, B, &["-f", "0.15", "-c", "0:count", "-c", "4:collapse"])[0],
        "chr1\t100\t200\ta1\t0\t+\t1\tb2"
    );
    assert_eq!(
        fx.run(A, B, &["-f", "0.5", "-c", "0:count"])[0],
        "chr1\t100\t200\ta1\t0\t+\t0"
    );
    assert_eq!(
        fx.run(A, B, &["-s", "-c", "4:collapse"])[0],
        "chr1\t100\t200\ta1\t0\t+\tb1"
    );
}

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn vcf_b_info() {
    let fx = Fx::new();
    let v = fx.file(
        "b.vcf",
        &format!(
            "{VCF_HEAD}chr1\t155\t.\tA\tC\t10\t.\tDP=10\nchr1\t170\t.\tA\tC\t30\t.\tDP=20\nchr1\t450\t.\tA\tC\t.\t.\tNODP\n"
        ),
    );
    let a = fx.file(
        "a.bed",
        "chr1\t100\t200\ta1\nchr1\t440\t460\ta2\nchr1\t600\t700\ta3\n",
    );
    let o = gloci()
        .args([
            "map",
            "-a",
            &a,
            "-b",
            &v,
            "-c",
            "INFO/DP:mean",
            "-c",
            "QUAL:sum",
            "-c",
            "0:count",
        ])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out,
        "chr1\t100\t200\ta1\t15\t40\t2\nchr1\t440\t460\ta2\t.\t.\t1\nchr1\t600\t700\ta3\t.\t.\t0\n"
    );
}

#[test]
fn vcf_a_header_line() {
    let fx = Fx::new();
    let a = fx.file(
        "a.vcf",
        &format!("{VCF_HEAD}chr1\t150\t.\tA\tC\t.\t.\t.\nchr1\t900\t.\tA\tC\t.\t.\t.\n"),
    );
    let b = fx.file("b.bed", "chr1\t100\t200\tx\t3\nchr1\t120\t130\ty\t4\n");
    let o = gloci()
        .args(["map", "-a", &a, "-b", &b, "-c", "5:sum", "-c", "0:count"])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t5:sum\t0:count\n\
         chr1\t150\t.\tA\tC\t.\t.\t.\t3\t1\n\
         chr1\t900\t.\tA\tC\t.\t.\t.\t.\t0\n"
    );
}

#[test]
fn bases_and_indicator() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(A, B, &["-c", "0:bases", "-c", "0:indicator"]),
        [
            "chr1\t100\t200\ta1\t0\t+\t30\t1",
            "chr1\t300\t400\ta2\t0\t-\t5\t1",
            "chr2\t5\t15\ta3\t0\t+\t0\t0",
        ]
    );
}

#[test]
fn skip_unmatched() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(A, B, &["-c", "0:count", "--skip-unmatched"]),
        ["chr1\t100\t200\ta1\t0\t+\t2", "chr1\t300\t400\ta2\t0\t-\t1"]
    );
}

#[test]
fn several_b_files_in_order() {
    let fx = Fx::new();
    let a = fx.file("a.bed", A);
    let b1 = fx.file("b1.bed", "chr1\t150\t160\tx\t1\n");
    let b2 = fx.file("b2.bed", "chr1\t120\t130\ty\t5\n");
    let o = gloci()
        .args([
            "map",
            "-a",
            &a,
            "-b",
            &b1,
            &b2,
            "-c",
            "4:collapse",
            "-c",
            "5:sum",
        ])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(
        out.starts_with("chr1\t100\t200\ta1\t0\t+\tx,y\t6\n"),
        "{out}"
    );
}

#[test]
fn errors() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", B));
    let v = fx.file("b.vcf", &format!("{VCF_HEAD}chr1\t155\t.\tA\tC\t.\t.\t.\n"));
    // usage
    assert_eq!(code(&["map", "-a", &a, "-b", &b]), 2);
    assert_eq!(code(&["map", "-a", &a, "-b", &b, "-c", "5"]), 2);
    assert_eq!(code(&["map", "-a", &a, "-b", &b, "-c", "5:quantile"]), 2);
    assert_eq!(code(&["map", "-a", &a, "-b", &b, "-c", "0:sum"]), 2);
    assert_eq!(code(&["map", "-a", &a, "-b", &b, "-c", "5:tmean:0.5"]), 2);
    assert_eq!(
        code(&[
            "map", "-a", &a, "-b", &b, "-c", "5:sum", "-w", "5", "-f", "0.5"
        ]),
        2
    );
    gloci()
        .args(["map", "-a", &a, "-b", &b, &v, "-c", "0:count"])
        .assert()
        .code(2)
        .stderr(
            contains("format family")
                .and(contains("b.bed is BED"))
                .and(contains("b.vcf is VCF")),
        );
    // and the other way round
    gloci()
        .args(["map", "-a", &a, "-b", &v, &b, "-c", "0:count"])
        .assert()
        .code(2)
        .stderr(
            contains("b.vcf is VCF")
                .and(contains("b.bed is BED"))
                .and(contains("Vcf").not()),
        );
    // a value that is not a number is a data error
    gloci()
        .args(["map", "-a", &a, "-b", &b, "-c", "4:sum"])
        .assert()
        .code(1)
        .stderr(contains("b.bed line 1").and(contains("4:sum")));
    // unsorted A
    let bad = fx.file("bad.bed", "chr1\t300\t400\tx\nchr1\t100\t200\ty\n");
    gloci()
        .args(["map", "-a", &bad, "-b", &b, "-c", "0:count"])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"));
    // --unsorted accepts it, in A's order
    let o = gloci()
        .args(["map", "-a", &bad, "-b", &b, "-c", "0:count", "--unsorted"])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(out, "chr1\t300\t400\tx\t1\nchr1\t100\t200\ty\t2\n");
}

#[test]
fn help_has_examples() {
    let o = gloci().args(["map", "--help"]).assert().success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(out.contains("Examples:"), "{out}");
    assert!(
        out.contains("--skip-unmatched") && out.contains("--unsorted"),
        "{out}"
    );
}

#[test]
fn determinism() {
    let fx = Fx::new();
    let a = fx.file("a.bed", A);
    let b = fx.file("b.bed", B);
    for extra in [
        &["-c", "5:sum,5:mean,0:count,4:collapse"][..],
        &["-c", "0:bases", "-w", "50"],
        &["-c", "5:max", "--unsorted", "--skip-unmatched"],
    ] {
        let run = |t: &str| {
            let mut v = vec!["map", "-a", &a, "-b", &b, "-t", t];
            v.extend_from_slice(extra);
            gloci()
                .args(&v)
                .assert()
                .success()
                .get_output()
                .stdout
                .clone()
        };
        let one = run("1");
        assert!(!one.is_empty());
        assert_eq!(one, run("2"), "{extra:?}");
        assert_eq!(one, run("8"), "{extra:?}");
    }
}

#[test]
fn two_input_ops_end_to_end() {
    let fx = Fx::new();
    // a1 [100,200): b1 [150,160) shares 10, b2 [180,260) shares 20
    // a3 [5,15): no B
    let got = fx.run(
        A,
        B,
        &[
            "-c",
            "5:wmean,5:bases_uniq,5:frac_a,5:frac_b,0:records,5:max_record",
        ],
    );
    assert_eq!(
        got,
        [
            "chr1\t100\t200\ta1\t0\t+\t3.5\t30\t0.3\t0.3333333333333333\t\
chr1|150|160|b1|2.5|+,chr1|180|260|b2|4|-\tchr1|180|260|b2|4|-",
            "chr1\t300\t400\ta2\t0\t-\t.\t5\t0.05\t0.047619047619047616\tchr1|395|500|b3|.|+\t.",
            "chr2\t5\t15\ta3\t0\t+\t.\t0\t.\t.\t.\t.",
        ]
    );
}

#[test]
fn spread_ops_end_to_end() {
    let fx = Fx::new();
    let got = fx.run(
        A,
        B,
        &[
            "-c",
            "5:mode,5:antimode,5:mad,5:quantile:0.5,5:tmean:0-1,5:cv",
        ],
    );
    assert_eq!(
        got[0],
        "chr1\t100\t200\ta1\t0\t+\t2.5\t2.5\t0.75\t3.25\t3.25\t0.32635697593225266"
    );
}

#[test]
fn new_ops_are_identical_for_every_thread_count() {
    let fx = Fx::new();
    let spec =
        "5:cv,5:mad,5:quantile:0.3,5:tmean:0.1-0.9,5:wmean,5:mode,5:bases_uniq,5:frac_a,0:records";
    let one = fx.run(A, B, &["--threads", "1", "-c", spec]);
    for t in ["2", "4"] {
        assert_eq!(
            fx.run(A, B, &["--threads", t, "-c", spec]),
            one,
            "--threads {t}"
        );
    }
}

#[test]
fn usage_errors_for_op_arguments() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", B));
    for (spec, want) in [
        ("5:quantile:2", "must be a number from 0 to 1"),
        ("5:tmean:0.9-0.1", "LO below HI"),
    ] {
        gloci()
            .args(["map", "-a", &a, "-b", &b, "-c", spec])
            .assert()
            .code(2)
            .stderr(contains(want));
    }
}

#[test]
fn split_bases_uniq_and_frac_a_use_the_blocks_of_a() {
    let fx = Fx::new();
    // identical A coordinates, different blocks: 10+10 at the ends, 50+50 contiguous, one block
    let a = "c\t0\t100\tA1\t0\t+\t0\t100\t0\t2\t10,10\t0,90\n\
             c\t0\t100\tA2\t0\t+\t0\t100\t0\t2\t50,50\t0,50\n\
             c\t0\t100\tA3\t0\t+\t0\t100\t0\t1\t100\t0\n";
    let b = "c\t5\t95\tb1\t1\t+\nc\t40\t60\tb2\t1\t+\n";
    let got = fx.run(a, b, &["--split", "-c", "0:bases_uniq", "-c", "0:frac_a"]);
    let tails: Vec<&str> = got
        .iter()
        .map(|l| l.splitn(13, '\t').last().unwrap()) // after the 12 BED columns
        .collect();
    assert_eq!(tails, ["10\t0.5", "90\t0.9", "90\t0.9"], "{got:?}");
}
