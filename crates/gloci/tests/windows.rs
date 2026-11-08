//! `gloci windows`: size/step, count, ids, intervals input, points,
//! usage errors and determinism.
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

fn genome(dir: &Path) -> String {
    s(&plain_file(dir, "g.genome", "c1\t250\nc2\t10\n"))
}

#[test]
fn size_step() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", "c1\t250\n"));
    assert_eq!(
        out(&["windows", "-g", &g, "-w", "100", "--step", "50"]),
        "c1\t0\t100\nc1\t50\t150\nc1\t100\t200\nc1\t150\t250\nc1\t200\t250\n"
    );
}

#[test]
fn step_larger_than_size_leaves_gaps() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", "c1\t250\n"));
    assert_eq!(
        out(&["windows", "-g", &g, "-w", "50", "--step", "100"]),
        "c1\t0\t50\nc1\t100\t150\nc1\t200\t250\n"
    );
}

#[test]
fn short_chromosome_gets_one_window_in_genome_order() {
    let dir = tempdir();
    let g = genome(dir.path());
    assert_eq!(
        out(&["windows", "-g", &g, "-w", "100"]),
        "c1\t0\t100\nc1\t100\t200\nc1\t200\t250\nc2\t0\t10\n"
    );
}

#[test]
fn count() {
    let dir = tempdir();
    let g = genome(dir.path());
    // The last window takes the remainder (bedtools).
    assert_eq!(
        out(&["windows", "-g", &g, "-n", "3"]),
        "c1\t0\t83\nc1\t83\t166\nc1\t166\t250\nc2\t0\t3\nc2\t3\t6\nc2\t6\t10\n"
    );
}

#[test]
fn count_larger_than_length_skips_the_interval() {
    let dir = tempdir();
    let g = genome(dir.path());
    let o = gloci()
        .args(["windows", "-g", &g, "-n", "20"])
        .assert()
        .success();
    let text = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(text.lines().count(), 20);
    assert!(text.lines().all(|l| l.starts_with("c1\t")));
    o.stderr(contains("c2:0-10"));
}

#[test]
fn id_and_reverse() {
    let dir = tempdir();
    let g = genome(dir.path());
    assert_eq!(
        out(&["windows", "-g", &g, "-n", "3", "--id-column"]),
        "c1\t0\t83\t1\nc1\t83\t166\t2\nc1\t166\t250\t3\nc2\t0\t3\t1\nc2\t3\t6\t2\nc2\t6\t10\t3\n"
    );
    assert_eq!(
        out(&["windows", "-g", &g, "-w", "100", "--reverse", "--id-column"]),
        "c1\t0\t100\t3\nc1\t100\t200\t2\nc1\t200\t250\t1\nc2\t0\t10\t1\n"
    );
}

#[test]
fn intervals_input() {
    let dir = tempdir();
    // Unsorted, with extra columns: windows come in input order, from the
    // footprint only.
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "c2\t0\t3\tx\ty\nc1\t10\t20\nc1\t30\t37\n",
    ));
    assert_eq!(
        out(&["windows", "-i", &f, "-w", "4", "--id-column"]),
        "c2\t0\t3\t1\nc1\t10\t14\t1\nc1\t14\t18\t2\nc1\t18\t20\t3\nc1\t30\t34\t1\nc1\t34\t37\t2\n"
    );
    assert_eq!(
        out(&["windows", "-i", &f, "-n", "4", "--id-column", "--reverse"]),
        "c1\t10\t12\t4\nc1\t12\t14\t3\nc1\t14\t16\t2\nc1\t16\t20\t1\n\
c1\t30\t31\t4\nc1\t31\t32\t3\nc1\t32\t33\t2\nc1\t33\t37\t1\n"
    );
}

#[test]
fn intervals_from_stdin() {
    let o = gloci()
        .args(["windows", "-i", "-", "-w", "5"])
        .write_stdin("c1\t0\t7\n")
        .assert()
        .success();
    assert_eq!(
        String::from_utf8(o.get_output().stdout.clone()).unwrap(),
        "c1\t0\t5\nc1\t5\t7\n"
    );
}

#[test]
fn points_skipped() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "p.bed",
        "c1\t5\t5\nc1\t10\t12\nc1\t20\t20\n",
    ));
    assert_eq!(out(&["windows", "-i", &f, "-w", "10"]), "c1\t10\t12\n");
}

#[test]
fn usage_errors() {
    let dir = tempdir();
    let g = genome(dir.path());
    let f = s(&plain_file(dir.path(), "a.bed", "c1\t0\t5\n"));
    let fail = |args: &[&str]| gloci().args(args).assert().code(2);
    fail(&["windows", "-w", "5"]);
    fail(&["windows", "-g", &g]);
    fail(&["windows", "-g", &g, "-i", &f, "-w", "5"]);
    fail(&["windows", "-g", &g, "-w", "5", "-n", "2"]);
    fail(&["windows", "-g", &g, "-n", "2", "--step", "1"]);
    fail(&["windows", "-g", &g, "-w", "0"]).stderr(contains("at least 1"));
    fail(&["windows", "-g", &g, "-n", "0"]).stderr(contains("at least 1"));
    fail(&["windows", "-g", &g, "-w", "5", "--step", "0"]).stderr(contains("at least 1"));
    // A genome without lengths cannot be cut.
    let nolen = s(&plain_file(
        dir.path(),
        "h.vcf",
        "##fileformat=VCFv4.3\n##contig=<ID=c1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ));
    fail(&["windows", "-g", &nolen, "-w", "5"]).stderr(contains("no length"));
}

#[test]
fn determinism() {
    let dir = tempdir();
    let g = genome(dir.path());
    let run = |t: &str| out(&["windows", "-g", &g, "-w", "7", "--id-column", "-t", t]);
    let one = run("1");
    assert_eq!(one, run("2"));
    assert_eq!(one, run("8"));
}
