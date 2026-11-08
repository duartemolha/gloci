//! `gloci complement`: with and without a genome, points, bookends, built-in
//! genomes, features past the chromosome end, usage errors and determinism.
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

#[test]
fn with_genome() {
    let dir = tempdir();
    let g = s(&plain_file(
        dir.path(),
        "g.genome",
        "chr2\t100\nchr1\t50\nchr3\t30\n",
    ));
    // Genome order, an empty chromosome (chr2) as one full row, leading and
    // trailing gaps, and a feature that covers a whole chromosome (chr3).
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t10\t20\nchr1\t15\t30\nchr3\t0\t30\n",
    ));
    assert_eq!(
        out(&["complement", &f, "-g", &g]),
        "chr2\t0\t100\nchr1\t0\t10\nchr1\t30\t50\n"
    );
}

#[test]
fn without_genome() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t9\nchr3\t1\t2\nchr3\t7\t8\n",
    ));
    assert_eq!(
        out(&["complement", &f]),
        "chr1\t20\t30\nchr3\t2\t7\n",
        "no leading or trailing gap, nothing for a chromosome with one feature"
    );
}

#[test]
fn points_ignored() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "p.bed",
        "chr1\t5\t5\nchr1\t10\t20\nchr1\t25\t25\nchr1\t30\t40\nchr2\t7\t7\n",
    ));
    // The points neither split the gap nor open a chromosome.
    assert_eq!(out(&["complement", &f]), "chr1\t20\t30\n");
    let g = s(&plain_file(dir.path(), "g.genome", "chr1\t50\nchr2\t9\n"));
    assert_eq!(
        out(&["complement", &f, "-g", &g]),
        "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\nchr2\t0\t9\n"
    );
}

#[test]
fn bookended_no_gap() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "b.bed",
        "chr1\t0\t10\nchr1\t10\t20\nchr1\t20\t50\n",
    ));
    let g = s(&plain_file(dir.path(), "g.genome", "chr1\t50\n"));
    assert_eq!(out(&["complement", &f]), "");
    assert_eq!(out(&["complement", &f, "-g", &g]), "");
}

#[test]
fn builtin_genome() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.bed", "chrI\t100\t200\n"));
    let o = out(&["complement", &f, "-g", "sacCer3"]);
    let rows: Vec<&str> = o.lines().collect();
    // chrI gives two rows, the other 16 chromosomes one each.
    assert_eq!(rows.len(), 18);
    assert_eq!(rows[0], "chrI\t0\t100");
    assert_eq!(rows[1], "chrI\t200\t230218");
    assert_eq!(rows[2], "chrII\t0\t813184");
    assert_eq!(
        rows.iter()
            .map(|r| r.split('\t').next().unwrap())
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        17
    );
}

#[test]
fn past_end_error() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t0\t10\nchr1\t40\t60\n",
    ));
    let g = s(&plain_file(dir.path(), "g.genome", "chr1\t50\n"));
    gloci()
        .args(["complement", &f, "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("a.bed line 2"))
        .stderr(contains("extends past the end of chr1 (length 50)"));
}

#[test]
fn genome_without_lengths_is_usage_error() {
    let dir = tempdir();
    let vcf = s(&plain_file(
        dir.path(),
        "h.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ));
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t0\t10\n"));
    gloci()
        .args(["complement", &f, "-g", &vcf])
        .assert()
        .code(2)
        .stderr(contains("no length for chr1"));
}

#[test]
fn unsorted_hint_names_sort_flag() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "u.bed",
        "chr1\t10\t20\nchr1\t1\t5\n",
    ));
    gloci()
        .args(["complement", &f])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"))
        .stderr(contains("--sort"));
    gloci()
        .args(["merge", &f])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"))
        .stderr(contains("--sort"));
}

#[test]
fn determinism() {
    let dir = tempdir();
    let a = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\t8\t9\nchr2\t1\t4\n",
    ));
    let b = s(&plain_file(dir.path(), "b.bed", "chr1\t4\t7\nchr3\t0\t3\n"));
    let run = |t: &str| out(&["complement", &a, &b, "-g", "hg19", "-t", t]);
    let base = run("1");
    assert!(base.starts_with("chr1\t0\t1\nchr1\t7\t8\n"));
    assert_eq!(run("2"), base);
    assert_eq!(run("8"), base);
    let run = |t: &str| out(&["complement", &a, &b, "--sort", "-t", t]);
    let base = run("1");
    assert_eq!(run("2"), base);
    assert_eq!(run("8"), base);
}
