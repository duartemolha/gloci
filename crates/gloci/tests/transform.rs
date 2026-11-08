//! `gloci slop` and `gloci flank`: tables, strand, --pct, GFF rewriting,
//! usage errors, the order warning and determinism.
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

/// chr1 is 1000 long, chr2 300.
fn genome(dir: &Path) -> String {
    s(&plain_file(dir, "g.genome", "chr1\t1000\nchr2\t300\n"))
}

#[test]
fn slop_table() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t100\t200\tx\n"));
    assert_eq!(
        out(&["slop", &f, "-b", "10", "-g", &g]),
        "chr1\t90\t210\tx\n"
    );
    let f = s(&plain_file(
        dir.path(),
        "c.bed",
        "chr1\t5\t20\nchr1\t950\t990\n",
    ));
    assert_eq!(
        out(&["slop", &f, "-b", "10", "-g", &g]),
        "chr1\t0\t30\nchr1\t940\t1000\n",
        "clipped at 0 and at the chromosome end"
    );
    assert_eq!(
        out(&["slop", &f, "-b", "10"]),
        "chr1\t0\t30\nchr1\t940\t1000\n",
        "without -g only 0 clips"
    );
    let f = s(&plain_file(dir.path(), "far.bed", "chr1\t1990\t2000\n"));
    assert_eq!(out(&["slop", &f, "-b", "10"]), "chr1\t1980\t2010\n");
    let f = s(&plain_file(
        dir.path(),
        "m.bed",
        "chr1\t100\t200\tn\t0\t-\n",
    ));
    assert_eq!(
        out(&["slop", &f, "--left", "5", "--right", "20", "-s", "-g", &g]),
        "chr1\t80\t205\tn\t0\t-\n"
    );
    assert_eq!(
        out(&["slop", &f, "--left", "5", "--right", "20", "-g", &g]),
        "chr1\t95\t220\tn\t0\t-\n",
        "without -s the strand is ignored"
    );
    let f = s(&plain_file(dir.path(), "p.bed", "chr1\t1000\t2000\n"));
    assert_eq!(
        out(&["slop", &f, "--pct", "-b", "0.5"]),
        "chr1\t500\t2500\n"
    );
    assert_eq!(
        out(&["slop", &f, "--pct", "-b", "50%"]),
        "chr1\t500\t2500\n"
    );
    assert_eq!(out(&["slop", &f, "--pct", "-b", "2"]), "chr1\t0\t4000\n");
    let f = s(&plain_file(dir.path(), "q.bed", "chr1\t0\t7\n"));
    assert_eq!(
        out(&["slop", &f, "--pct", "--right", "0.5"]),
        "chr1\t0\t10\n",
        "3.5 is floored to 3"
    );
}

#[test]
fn flank_table() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t100\t200\tx\n"));
    assert_eq!(
        out(&["flank", &f, "-b", "10", "-g", &g]),
        "chr1\t90\t100\tx\nchr1\t200\t210\tx\n"
    );
    assert_eq!(
        out(&["flank", &f, "--left", "10", "-g", &g]),
        "chr1\t90\t100\tx\n"
    );
    assert_eq!(out(&["flank", &f, "--right", "10"]), "chr1\t200\t210\tx\n");
    let f = s(&plain_file(
        dir.path(),
        "e.bed",
        "chr1\t0\t20\nchr1\t980\t1000\n",
    ));
    assert_eq!(
        out(&["flank", &f, "-b", "10", "-g", &g]),
        "chr1\t20\t30\nchr1\t970\t980\n",
        "a strip clipped to nothing is not written"
    );
    let f = s(&plain_file(
        dir.path(),
        "s.bed",
        "chr1\t100\t200\tn\t0\t-\n",
    ));
    assert_eq!(
        out(&["flank", &f, "--left", "5", "-s"]),
        "chr1\t200\t205\tn\t0\t-\n",
        "-s turns --left into the strip after a minus-strand record"
    );
}

#[test]
fn gff_rewrite() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(
        dir.path(),
        "a.gff",
        "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n",
    ));
    assert_eq!(
        out(&["slop", &f, "-b", "10", "-g", &g]),
        "chr1\tsrc\tgene\t91\t210\t.\t+\t.\tID=g1\n"
    );
    assert_eq!(
        out(&["flank", &f, "-b", "10", "-g", &g]),
        "chr1\tsrc\tgene\t91\t100\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\t201\t210\t.\t+\t.\tID=g1\n"
    );
}

#[test]
fn usage_errors() {
    let dir = tempdir();
    let g = genome(dir.path());
    let v = s(&plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    ));
    for cmd in ["slop", "flank"] {
        gloci()
            .args([cmd, &v, "-b", "5"])
            .assert()
            .code(2)
            .stderr(contains("BED and GFF/GTF"));
    }
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t2\n"));
    gloci().args(["slop", &f]).assert().code(2);
    gloci()
        .args(["slop", &f, "-b", "1", "--left", "2"])
        .assert()
        .code(2);
    gloci()
        .args(["flank", &f, "-b", "1.5"])
        .assert()
        .code(2)
        .stderr(contains("--pct"));
    gloci()
        .args(["slop", &f, "--pct", "-b", "abc"])
        .assert()
        .code(2);
    let f = s(&plain_file(dir.path(), "x.bed", "chrZ\t1\t2\n"));
    gloci()
        .args(["slop", &f, "-b", "1", "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("chrZ"));
}

#[test]
fn beyond_chromosome_end() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t990\t1500\nchr1\t1100\t1200\n",
    ));
    for cmd in ["slop", "flank"] {
        gloci()
            .args([cmd, &f, "-b", "10", "-g", &g])
            .assert()
            .code(1)
            .stderr(contains("a.bed line 2"))
            .stderr(contains("beyond the end of chr1 (length 1000)"));
    }
}

#[test]
fn order_warning_once() {
    let dir = tempdir();
    let g = genome(dir.path());
    // Plain slop keeps the order: no warning.
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t100\t200\nchr1\t150\t160\nchr1\t300\t400\n",
    ));
    gloci()
        .args(["slop", &f, "-b", "10", "-g", &g])
        .assert()
        .success()
        .stderr(predicates::str::is_empty());
    // Flank of a long feature followed by a short one breaks the order twice;
    // the warning comes once and names the first offending output line.
    let f = s(&plain_file(
        dir.path(),
        "b.bed",
        "chr1\t100\t500\nchr1\t120\t130\nchr1\t140\t150\n",
    ));
    let o = gloci()
        .args(["flank", &f, "-b", "50", "-g", &g])
        .assert()
        .success()
        .stderr(contains("output line 3 is out of order"))
        .stderr(contains("gloci sort"));
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    assert_eq!(err.matches("out of order").count(), 1, "{err}");
    // -s on mixed strands in slop breaks the order too.
    let f = s(&plain_file(
        dir.path(),
        "c.bed",
        "chr1\t100\t200\tn\t0\t-\nchr1\t110\t120\tm\t0\t+\n",
    ));
    gloci()
        .args(["slop", &f, "--left", "50", "-s"])
        .assert()
        .success()
        .stderr(contains("output line 2 is out of order"));
}

#[test]
fn determinism() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t100\t200\tn\t0\t-\nchr1\t300\t310\tm\t0\t+\nchr2\t5\t50\n",
    ));
    for cmd in ["slop", "flank"] {
        let run = |t: &str| out(&[cmd, &f, "-b", "20", "-s", "-g", &g, "-t", t]);
        let base = run("1");
        assert!(!base.is_empty());
        assert_eq!(run("2"), base);
        assert_eq!(run("8"), base);
    }
}
