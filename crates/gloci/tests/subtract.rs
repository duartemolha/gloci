//! `gloci subtract`: pieces, criteria, `--whole`, `--sum-overlaps`, GFF and
//! VCF A, `--split`, windows, usage errors and determinism.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn code(args: &[&str]) -> i32 {
    gloci().args(args).output().unwrap().status.code().unwrap()
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
    /// `subtract -a a -b b extra...` over two BED texts.
    fn run(&self, a: &str, b: &str, extra: &[&str]) -> String {
        let (af, bf) = (self.file("a.bed", a), self.file("b.bed", b));
        let mut v = vec!["subtract", "-a", &af, "-b", &bf];
        v.extend_from_slice(extra);
        out(&v)
    }
}

#[test]
fn pieces() {
    let fx = Fx::new();
    let a = "chr1\t0\t100\ta1\t0\t+\n";
    assert_eq!(
        fx.run(a, "chr1\t10\t20\nchr1\t30\t40\n", &[]),
        "chr1\t0\t10\ta1\t0\t+\nchr1\t20\t30\ta1\t0\t+\nchr1\t40\t100\ta1\t0\t+\n"
    );
    // a B covering A removes it
    assert_eq!(fx.run(a, "chr1\t0\t200\n", &[]), "");
    // a bookended B leaves A unchanged
    assert_eq!(fx.run(a, "chr1\t100\t110\n", &[]), a);
    // a B point inside A leaves A unchanged
    assert_eq!(fx.run(a, "chr1\t50\t50\n", &[]), a);
    // an A point inside B is removed; outside B it stays
    let pt = "chr1\t50\t50\tp\n";
    assert_eq!(fx.run(pt, "chr1\t40\t60\n", &[]), "");
    assert_eq!(fx.run(pt, "chr1\t60\t70\n", &[]), pt);
    // overlapping Bs remove their union; other chromosomes are untouched
    assert_eq!(
        fx.run(
            "chr1\t0\t100\ta\nchr2\t0\t10\tb\n",
            "chr1\t10\t50\nchr1\t40\t90\n",
            &[]
        ),
        "chr1\t0\t10\ta\nchr1\t90\t100\ta\nchr2\t0\t10\tb\n"
    );
}

#[test]
fn criteria() {
    let fx = Fx::new();
    let a = "chr1\t0\t100\ta1\t0\t+\n";
    // B1 covers 30% of A (and all of B1), B2 covers 50%
    let b = "chr1\t0\t30\tb1\t0\t+\nchr1\t50\t100\tb2\t0\t-\n";
    assert_eq!(fx.run(a, b, &[]), "chr1\t30\t50\ta1\t0\t+\n");
    assert_eq!(
        fx.run(a, b, &["-f", "0.5"]),
        "chr1\t0\t50\ta1\t0\t+\n",
        "only B2 meets -f 0.5 of A"
    );
    assert_eq!(fx.run(a, b, &["-s"]), "chr1\t30\t100\ta1\t0\t+\n");
    assert_eq!(fx.run(a, b, &["-S"]), "chr1\t0\t50\ta1\t0\t+\n");
}

#[test]
fn whole_and_sum() {
    let fx = Fx::new();
    let a = "chr1\t0\t100\ta1\t0\t+\n";
    let b = "chr1\t0\t30\nchr1\t50\t80\n";
    assert_eq!(fx.run(a, b, &["--sum-overlaps", "-f", "0.5"]), "");
    assert_eq!(fx.run(a, b, &["-f", "0.5"]), a, "no single B reaches 50%");
    assert_eq!(fx.run(a, b, &["--whole"]), "");
    // the sum falls short of the fraction: A stays, unchanged
    assert_eq!(fx.run(a, b, &["--sum-overlaps", "-f", "0.7"]), a);
    // no -f: like --whole
    assert_eq!(fx.run(a, b, &["--sum-overlaps"]), "");
    assert_eq!(fx.run(a, "chr1\t200\t300\n", &["--sum-overlaps"]), a);
    // the strand rule still applies
    let bs = "chr1\t0\t60\tb\t0\t-\n";
    assert_eq!(fx.run(a, bs, &["--sum-overlaps", "-f", "0.5", "-s"]), a);
    assert_eq!(fx.run(a, bs, &["--sum-overlaps", "-f", "0.5", "-S"]), "");
    // --whole with a threshold uses it
    assert_eq!(fx.run(a, b, &["--whole", "-f", "0.5"]), a);
}

#[test]
fn sum_overlaps_option_errors() {
    let fx = Fx::new();
    let (a, b) = (
        fx.file("a.bed", "chr1\t0\t10\n"),
        fx.file("b.bed", "chr1\t0\t5\n"),
    );
    for extra in [
        &["--sum-overlaps", "-F", "0.5"][..],
        &["--sum-overlaps", "-r", "-f", "0.5"],
        &["--sum-overlaps", "-e", "-f", "0.5"],
        &["--sum-overlaps", "--min-bases", "2"],
        &["--sum-overlaps", "--whole"],
        &["--sum-overlaps", "-f", "0.5", "-w", "5"],
    ] {
        let mut v = vec!["subtract", "-a", &a, "-b", &b];
        v.extend_from_slice(extra);
        assert_eq!(code(&v), 2, "{extra:?}");
    }
}

#[test]
fn windows() {
    let fx = Fx::new();
    let a = "chr1\t0\t100\ta1\t0\t+\n";
    assert_eq!(
        fx.run(a, "chr1\t10\t20\n", &["-w", "5"]),
        "chr1\t0\t5\ta1\t0\t+\nchr1\t25\t100\ta1\t0\t+\n"
    );
    // a B that only touches A once padded still counts
    assert_eq!(
        fx.run(a, "chr1\t102\t110\n", &["--left", "5"]),
        "chr1\t0\t97\ta1\t0\t+\n"
    );
    // a padded point B removes bases
    assert_eq!(
        fx.run(a, "chr1\t50\t50\n", &["-w", "3"]),
        "chr1\t0\t47\ta1\t0\t+\nchr1\t53\t100\ta1\t0\t+\n"
    );
}

#[test]
fn gff_pieces() {
    let fx = Fx::new();
    // 1-based closed 1..100 is [0, 100)
    let a = fx.file("a.gff", "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=x\n");
    let b = fx.file("b.bed", "chr1\t9\t20\n");
    assert_eq!(
        out(&["subtract", "-a", &a, "-b", &b]),
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=x\nchr1\tsrc\tgene\t21\t100\t.\t+\t.\tID=x\n"
    );
}

const VHEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn vcf_a_whole_records() {
    let fx = Fx::new();
    let a = fx.file(
        "a.vcf",
        &format!(
            "{VHEAD}chr1\t12\t.\tA\tC\t.\t.\t.\nchr1\t15\t.\tACGT\tA\t.\t.\t.\nchr1\t100\t.\tG\tT\t.\t.\t.\n"
        ),
    );
    let b = fx.file("b.bed", "chr1\t10\t20\tb\n");
    let sub = out(&["subtract", "-a", &a, "-b", &b]);
    let lines: Vec<&str> = sub.lines().collect();
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert!(lines.iter().any(|l| l.starts_with("##glociCommand=")));
    assert_eq!(lines.last().unwrap(), &"chr1\t100\t.\tG\tT\t.\t.\t.");
    assert_eq!(lines.len(), 4, "{sub}");
    // exactly what intersect --mode none keeps
    let strip = |t: &str| {
        t.lines()
            .filter(|l| !l.starts_with("##glociCommand"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    for extra in [&[][..], &["-w", "3"], &["--vcf-span", "ref"]] {
        let mut x = vec!["intersect", "--mode", "none", "-a", &a, "-b", &b];
        x.extend_from_slice(extra);
        let mut y = vec!["subtract", "-a", &a, "-b", &b];
        y.extend_from_slice(extra);
        assert_eq!(strip(&out(&y)), strip(&out(&x)), "{extra:?}");
    }
}

#[test]
fn split_blocks() {
    let fx = Fx::new();
    let a = "chr1\t0\t100\ta1\t0\t+\n";
    let b12 = "chr1\t0\t100\tb\t0\t+\t0\t100\t0\t2\t10,10\t0,90\n";
    assert_eq!(fx.run(a, b12, &[]), "");
    assert_eq!(fx.run(a, b12, &["--split"]), "chr1\t10\t90\ta1\t0\t+\n");
    // a B below BED12 is one block
    assert_eq!(
        fx.run(a, "chr1\t10\t20\n", &["--split"]),
        "chr1\t0\t10\ta1\t0\t+\nchr1\t20\t100\ta1\t0\t+\n"
    );
}

#[test]
fn errors_and_unsorted() {
    let fx = Fx::new();
    let a = "chr1\t50\t60\nchr1\t0\t10\n";
    let b = "chr1\t5\t8\n";
    let (af, bf) = (fx.file("a.bed", a), fx.file("b.bed", b));
    gloci()
        .args(["subtract", "-a", &af, "-b", &bf])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"));
    assert_eq!(
        out(&["subtract", "-a", &af, "-b", &bf, "--unsorted"]),
        "chr1\t50\t60\nchr1\t0\t5\nchr1\t8\t10\n"
    );
}

#[test]
fn determinism() {
    let fx = Fx::new();
    let a = fx.file(
        "a.bed",
        "chr1\t0\t100\ta1\t0\t+\nchr1\t120\t200\ta2\t0\t-\nchr2\t5\t15\ta3\t0\t+\nchr3\t0\t9\ta4\t0\t+\n",
    );
    let b = fx.file("b.bed", "chr1\t10\t20\nchr1\t30\t40\nchr2\t0\t8\n");
    let b2 = fx.file("b2.bed", "chr1\t150\t160\n");
    for extra in [
        &[][..],
        &["--whole"],
        &["--sum-overlaps", "-f", "0.5"],
        &["-w", "4"],
        &["--unsorted"],
    ] {
        let run = |t: &str| {
            let mut v = vec!["subtract", "-a", &a, "-b", &b, &b2, "-t", t];
            v.extend_from_slice(extra);
            out(&v)
        };
        let one = run("1");
        assert!(!one.is_empty());
        assert_eq!(one, run("2"), "{extra:?}");
        assert_eq!(one, run("8"), "{extra:?}");
    }
}

#[test]
fn help_has_examples() {
    let h = out(&["subtract", "--help"]);
    assert!(h.contains("Examples:"));
    assert!(h.contains("--sum-overlaps"));
    assert!(h.contains("--whole"));
}
