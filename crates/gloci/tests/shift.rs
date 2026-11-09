//! `gloci shift`: tables, strand, --pct, GFF rewriting, usage errors, the
//! order warning, dropped records and determinism.
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
fn shift_table() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t100\t200\tx\n"));
    assert_eq!(
        out(&["shift", &f, "--by", "10", "-g", &g]),
        "chr1\t110\t210\tx\n"
    );
    assert_eq!(
        out(&["shift", &f, "--by", "-10", "-g", &g]),
        "chr1\t90\t190\tx\n"
    );
    let f = s(&plain_file(
        dir.path(),
        "c.bed",
        "chr1\t5\t20\nchr1\t950\t990\n",
    ));
    assert_eq!(
        out(&["shift", &f, "--by", "-10", "-g", &g]),
        "chr1\t0\t10\nchr1\t940\t980\n",
        "clipped at 0"
    );
    assert_eq!(
        out(&["shift", &f, "--by", "30", "-g", &g]),
        "chr1\t35\t50\nchr1\t980\t1000\n",
        "clipped at the chromosome end"
    );
    assert_eq!(
        out(&["shift", &f, "--by", "30"]),
        "chr1\t35\t50\nchr1\t980\t1020\n",
        "without -g only 0 clips"
    );
    let f = s(&plain_file(
        dir.path(),
        "m.bed",
        "chr1\t100\t200\tn\t0\t-\nchr1\t300\t310\tp\t0\t+\nchr1\t400\t410\tu\t0\t.\n",
    ));
    assert_eq!(
        out(&["shift", &f, "--by", "20", "-s", "-g", &g]),
        "chr1\t80\t180\tn\t0\t-\nchr1\t320\t330\tp\t0\t+\nchr1\t420\t430\tu\t0\t.\n",
        "-s negates on - and treats . as +"
    );
    assert_eq!(
        out(&["shift", &f, "--by", "20", "-g", &g]),
        "chr1\t120\t220\tn\t0\t-\nchr1\t320\t330\tp\t0\t+\nchr1\t420\t430\tu\t0\t.\n",
        "without -s the strand is ignored"
    );
}

#[test]
fn shift_pct() {
    let dir = tempdir();
    let g = genome(dir.path());
    // Lengths 100 and 5 and 91 (odd): the product is floored, also when
    // negative (bedtools 2.31.1: -0.25 * 5 shifts by -2).
    let f = s(&plain_file(
        dir.path(),
        "p.bed",
        "chr1\t100\t200\nchr1\t500\t505\nchr1\t600\t691\n",
    ));
    assert_eq!(
        out(&["shift", &f, "--pct", "--by", "0.5", "-g", &g]),
        "chr1\t150\t250\nchr1\t502\t507\nchr1\t645\t736\n"
    );
    assert_eq!(
        out(&["shift", &f, "--pct", "--by", "-0.25", "-g", &g]),
        "chr1\t75\t175\nchr1\t498\t503\nchr1\t577\t668\n"
    );
    assert_eq!(
        out(&["shift", &f, "--pct", "--by", "25%", "-g", &g]),
        "chr1\t125\t225\nchr1\t501\t506\nchr1\t622\t713\n"
    );
    let m = s(&plain_file(
        dir.path(),
        "m.bed",
        "chr1\t100\t200\tn\t0\t-\n",
    ));
    assert_eq!(
        out(&["shift", &m, "--pct", "--by", "0.5", "-s", "-g", &g]),
        "chr1\t50\t150\tn\t0\t-\n"
    );
}

#[test]
fn gff_rewrite() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(
        dir.path(),
        "a.gff3",
        "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n",
    ));
    assert_eq!(
        out(&["shift", &f, "--by", "10", "-g", &g]),
        "chr1\tsrc\tgene\t111\t210\t.\t+\t.\tID=g1\n"
    );
}

#[test]
fn usage_errors() {
    let dir = tempdir();
    let v = s(&plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    ));
    gloci()
        .args(["shift", &v, "--by", "5"])
        .assert()
        .code(2)
        .stderr(contains("BED and GFF/GTF"));
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t2\n"));
    gloci().args(["shift", &f]).assert().code(2);
    gloci()
        .args(["shift", &f, "--by", "1.5"])
        .assert()
        .code(2)
        .stderr(contains("--pct"));
    gloci()
        .args(["shift", &f, "--pct", "--by", "abc"])
        .assert()
        .code(2);
    let g = genome(dir.path());
    let x = s(&plain_file(dir.path(), "x.bed", "chrZ\t1\t2\n"));
    gloci()
        .args(["shift", &x, "--by", "1", "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("chrZ"));
    let b = s(&plain_file(dir.path(), "b.bed", "chr1\t1100\t1200\n"));
    gloci()
        .args(["shift", &b, "--by", "1", "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("b.bed line 1"))
        .stderr(contains("beyond the end of chr1 (length 1000)"));
}

#[test]
fn wholly_past_an_edge() {
    let dir = tempdir();
    let g = genome(dir.path());
    // Non-empty features pushed wholly over an edge keep one base there
    // (as bedtools does); a point stays a point at the edge.
    let f = s(&plain_file(
        dir.path(),
        "d.bed",
        "chr1\t0\t5\tq\nchr1\t995\t1000\tp\nchr1\t10\t10\tz\nchr1\t500\t510\tk\n",
    ));
    gloci()
        .args(["shift", &f, "--by", "-10", "-g", &g])
        .assert()
        .success()
        .stderr(predicates::str::is_empty())
        .stdout("chr1\t0\t1\tq\nchr1\t985\t990\tp\nchr1\t0\t0\tz\nchr1\t490\t500\tk\n");
    gloci()
        .args(["shift", &f, "--by", "10", "-g", &g])
        .assert()
        .success()
        .stderr(predicates::str::is_empty())
        .stdout("chr1\t10\t15\tq\nchr1\t999\t1000\tp\nchr1\t20\t20\tz\nchr1\t510\t520\tk\n");
    let pt = s(&plain_file(dir.path(), "pt.bed", "chr1\t1000\t1000\ty\n"));
    assert_eq!(
        out(&["shift", &pt, "--by", "5", "-g", &g]),
        "chr1\t1000\t1000\ty\n"
    );
}

#[test]
fn gff_edge_interval() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(
        dir.path(),
        "e.gff3",
        "chr1\tsrc\tgene\t996\t1000\t.\t+\t.\tID=a\nchr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=b\n",
    ));
    assert_eq!(
        out(&["shift", &f, "--by", "20", "-g", &g]),
        "chr1\tsrc\tgene\t1000\t1000\t.\t+\t.\tID=a\nchr1\tsrc\tgene\t21\t25\t.\t+\t.\tID=b\n"
    );
    assert_eq!(
        out(&["shift", &f, "--by", "-20", "-g", &g]),
        "chr1\tsrc\tgene\t976\t980\t.\t+\t.\tID=a\nchr1\tsrc\tgene\t1\t1\t.\t+\t.\tID=b\n"
    );
}

#[test]
fn order_warning_only_when_shift_breaks_order() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t100\t200\nchr1\t150\t160\nchr1\t300\t400\n",
    ));
    gloci()
        .args(["shift", &f, "--by", "10", "-g", &g])
        .assert()
        .success()
        .stderr(predicates::str::is_empty());
    // -s on mixed strands breaks the order.
    let c = s(&plain_file(
        dir.path(),
        "c.bed",
        "chr1\t100\t200\tn\t0\t+\nchr1\t110\t120\tm\t0\t-\n",
    ));
    gloci()
        .args(["shift", &c, "--by", "50", "-s"])
        .assert()
        .success()
        .stderr(contains("output record 2 is out of order"))
        .stderr(contains("gloci sort"));
    // --pct moves a long feature further than a later short one.
    let p = s(&plain_file(
        dir.path(),
        "p.bed",
        "chr1\t100\t500\nchr1\t120\t130\n",
    ));
    gloci()
        .args(["shift", &p, "--pct", "--by", "0.5"])
        .assert()
        .success()
        .stderr(contains("output record 2 is out of order"));
    // An input that is already unsorted is not the shift's fault.
    let u = s(&plain_file(
        dir.path(),
        "u.bed",
        "chr1\t500\t600\nchr1\t100\t200\n",
    ));
    gloci()
        .args(["shift", &u, "--by", "3"])
        .assert()
        .success()
        .stderr(predicates::str::is_empty());
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
    let run = |t: &str| out(&["shift", &f, "--by", "20", "-s", "-g", &g, "-t", t]);
    let base = run("1");
    assert!(!base.is_empty());
    assert_eq!(run("2"), base);
    assert_eq!(run("8"), base);
}
