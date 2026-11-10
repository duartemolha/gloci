//! `gloci partition`: membership columns, count filters, merge-adjacent,
//! coordinates only, bedGraph value columns and genome fill, points,
//! `--sort`, errors and determinism.
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

const A: &str = "chr1\t0\t10\nchr1\t20\t30\nchr2\t5\t15\n";
const B: &str = "chr1\t5\t15\nchr1\t10\t20\nchr1\t25\t26\nchr3\t0\t5\n";
const C: &str = "chr1\t0\t50\nchr2\t5\t15\n";

struct Three {
    _dir: tempfile::TempDir,
    a: String,
    b: String,
    c: String,
}

fn three() -> Three {
    let dir = tempdir();
    let a = s(&plain_file(dir.path(), "a.bed", A));
    let b = s(&plain_file(dir.path(), "b.bed", B));
    let c = s(&plain_file(dir.path(), "c.bed", C));
    Three { _dir: dir, a, b, c }
}

#[test]
fn membership_columns_for_three_inputs() {
    let t = three();
    assert_eq!(
        out(&["partition", &t.a, &t.b, &t.c]),
        "chr1\t0\t5\t2\ta.bed,c.bed\n\
chr1\t5\t10\t3\ta.bed,b.bed,c.bed\n\
chr1\t10\t15\t2\tb.bed,c.bed\n\
chr1\t15\t20\t2\tb.bed,c.bed\n\
chr1\t20\t25\t2\ta.bed,c.bed\n\
chr1\t25\t26\t3\ta.bed,b.bed,c.bed\n\
chr1\t26\t30\t2\ta.bed,c.bed\n\
chr1\t30\t50\t1\tc.bed\n\
chr2\t5\t15\t2\ta.bed,c.bed\n\
chr3\t0\t5\t1\tb.bed\n"
    );
}

#[test]
fn names_replace_the_file_names() {
    let t = three();
    let o = out(&["partition", &t.a, &t.b, "--names", "tumour,normal"]);
    assert!(o.starts_with("chr1\t0\t5\t1\ttumour\n"), "{o}");
    assert!(o.contains("chr1\t5\t10\t2\ttumour,normal\n"), "{o}");
    gloci()
        .args(["partition", &t.a, &t.b, "--names", "x"])
        .assert()
        .code(2)
        .stderr(contains("one label per file"));
    gloci()
        .args(["partition", &t.a, &t.b, "--names", "x,x"])
        .assert()
        .code(2)
        .stderr(contains("--names"));
}

#[test]
fn count_filters_and_all() {
    let t = three();
    let min3 = out(&[
        "partition",
        &t.a,
        &t.b,
        &t.c,
        "--min-files",
        "3",
        "--no-membership",
    ]);
    assert_eq!(min3, "chr1\t5\t10\nchr1\t25\t26\n");
    let max1 = out(&[
        "partition",
        &t.a,
        &t.b,
        &t.c,
        "--max-files",
        "1",
        "--merge-adjacent",
    ]);
    assert_eq!(max1, "chr1\t30\t50\t1\tc.bed\nchr3\t0\t5\t1\tb.bed\n");
    let all = out(&[
        "partition",
        &t.a,
        &t.c,
        "--all",
        "--merge-adjacent",
        "--no-membership",
    ]);
    assert_eq!(all, "chr1\t0\t10\nchr1\t20\t30\nchr2\t5\t15\n");
    gloci()
        .args(["partition", &t.a, &t.c, "--all", "--min-files", "1"])
        .assert()
        .code(2);
    gloci()
        .args([
            "partition",
            &t.a,
            &t.c,
            "--min-files",
            "3",
            "--max-files",
            "1",
        ])
        .assert()
        .code(2)
        .stderr(contains("--min-files 3 is above --max-files 1"));
}

#[test]
fn min_files_zero_needs_empty() {
    let t = three();
    gloci()
        .args(["partition", &t.a, &t.c, "--min-files", "0"])
        .assert()
        .code(2)
        .stderr(contains("--min-files 0"))
        .stderr(contains("add --empty"));
}

#[test]
fn merge_adjacent_joins_equal_columns_only() {
    let dir = tempdir();
    let x = s(&plain_file(
        dir.path(),
        "x.bed",
        "chr1\t0\t10\nchr1\t10\t20\n",
    ));
    let y = s(&plain_file(
        dir.path(),
        "y.bed",
        "chr1\t0\t20\nchr1\t30\t35\n",
    ));
    assert_eq!(
        out(&["partition", &x, &y]),
        "chr1\t0\t10\t2\tx.bed,y.bed\nchr1\t10\t20\t2\tx.bed,y.bed\nchr1\t30\t35\t1\ty.bed\n"
    );
    assert_eq!(
        out(&["partition", &x, &y, "--merge-adjacent"]),
        "chr1\t0\t20\t2\tx.bed,y.bed\nchr1\t30\t35\t1\ty.bed\n"
    );
    // Different membership, touching: kept apart with labels, joined without.
    let p = s(&plain_file(dir.path(), "p.bed", "chr1\t0\t10\n"));
    let q = s(&plain_file(dir.path(), "q.bed", "chr1\t10\t20\n"));
    assert_eq!(
        out(&["partition", &p, &q, "--merge-adjacent"]),
        "chr1\t0\t10\t1\tp.bed\nchr1\t10\t20\t1\tq.bed\n"
    );
    assert_eq!(
        out(&["partition", &p, &q, "--merge-adjacent", "--no-membership"]),
        "chr1\t0\t20\n"
    );
}

#[test]
fn nested_features_count_a_file_once() {
    let dir = tempdir();
    let n = s(&plain_file(
        dir.path(),
        "n.bed",
        "chr1\t0\t30\nchr1\t5\t10\nchr1\t8\t40\n",
    ));
    let u = s(&plain_file(dir.path(), "u.bed", "chr1\t20\t25\n"));
    assert_eq!(
        out(&["partition", &n, &u, "--merge-adjacent"]),
        "chr1\t0\t20\t1\tn.bed\nchr1\t20\t25\t2\tn.bed,u.bed\nchr1\t25\t40\t1\tn.bed\n"
    );
    assert_eq!(
        out(&["partition", &n, "--no-membership"]),
        "chr1\t0\t5\nchr1\t5\t8\nchr1\t8\t10\nchr1\t10\t30\nchr1\t30\t40\n"
    );
}

#[test]
fn points_contribute_nothing() {
    let dir = tempdir();
    let x = s(&plain_file(dir.path(), "x.bed", "chr1\t0\t50\n"));
    let p = s(&plain_file(
        dir.path(),
        "p.bed",
        "chr1\t20\t20\nchr4\t7\t7\n",
    ));
    assert_eq!(
        out(&["partition", &x, &p]),
        "chr1\t0\t50\t1\tx.bed\n",
        "no boundary at the point and no row for its chromosome"
    );
}

#[test]
fn bedgraph_value_columns_and_empty_fill() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "f.bg",
        "chr1\t0\t10\t1\nchr1\t10\t20\t1\nchr1\t30\t40\t2\n",
    ));
    let h = s(&plain_file(
        dir.path(),
        "h.bg",
        "chr1\t5\t15\t3.5\nchr2\t0\t5\t4\n",
    ));
    assert_eq!(
        out(&["partition", &f, &h, "--bedgraph"]),
        "chr1\t0\t5\t1\t0\nchr1\t5\t10\t1\t3.5\nchr1\t10\t15\t1\t3.5\nchr1\t15\t20\t1\t0\n\
chr1\t30\t40\t2\t0\nchr2\t0\t5\t0\t4\n"
    );
    assert_eq!(
        out(&[
            "partition",
            &f,
            &h,
            "--bedgraph",
            "--merge-adjacent",
            "--empty",
            "NA"
        ]),
        "chr1\t0\t5\t1\tNA\nchr1\t5\t15\t1\t3.5\nchr1\t15\t20\t1\tNA\nchr1\t20\t30\tNA\tNA\n\
chr1\t30\t40\t2\tNA\nchr2\t0\t5\tNA\t4\n"
    );
    // -g with --empty covers the whole genome, empty chromosomes included.
    let g = s(&plain_file(
        dir.path(),
        "g.txt",
        "chr1\t100\nchr2\t50\nchr3\t10\nchr4\t7\n",
    ));
    assert_eq!(
        out(&["partition", &f, &h, "--bedgraph", "--empty", "0", "-g", &g]),
        "chr1\t0\t5\t1\t0\nchr1\t5\t10\t1\t3.5\nchr1\t10\t15\t1\t3.5\nchr1\t15\t20\t1\t0\n\
chr1\t20\t30\t0\t0\nchr1\t30\t40\t2\t0\nchr1\t40\t100\t0\t0\nchr2\t0\t5\t0\t4\n\
chr2\t5\t50\t0\t0\nchr3\t0\t10\t0\t0\nchr4\t0\t7\t0\t0\n"
    );
    // -g alone only fixes the chromosome order; nothing is filled.
    assert_eq!(
        out(&["partition", &f, &h, "--bedgraph", "-g", &g]),
        out(&["partition", &f, &h, "--bedgraph"])
    );
    gloci()
        .args(["partition", &f, &h, "--empty", "0"])
        .assert()
        .code(2);
}

#[test]
fn bedgraph_errors_say_where() {
    let dir = tempdir();
    let ov = s(&plain_file(
        dir.path(),
        "ov.bg",
        "chr1\t0\t10\t1\nchr1\t5\t15\t2\n",
    ));
    let h = s(&plain_file(dir.path(), "h.bg", "chr1\t5\t15\t3\n"));
    gloci()
        .args(["partition", &ov, &h, "--bedgraph"])
        .assert()
        .code(1)
        .stderr(contains("ov.bg"))
        .stderr(contains("overlaps an earlier interval"))
        .stderr(contains(
            "split the file so that no interval overlaps another",
        ));
    // Under --sort the location is a record count, and the message says so.
    gloci()
        .args(["partition", &ov, &h, "--bedgraph", "--sort"])
        .assert()
        .code(1)
        .stderr(contains("record 2"))
        .stderr(contains("counts records after --sort"));
    let nov = s(&plain_file(dir.path(), "nov.bed", "chr1\t0\t10\n"));
    gloci()
        .args(["partition", &nov, &h, "--bedgraph"])
        .assert()
        .code(1)
        .stderr(contains("nov.bed"))
        .stderr(contains("no value in column 4"));
    let long = s(&plain_file(dir.path(), "l.bg", "chr1\t0\t200\t1\n"));
    let g = s(&plain_file(dir.path(), "g.txt", "chr1\t100\n"));
    gloci()
        .args([
            "partition",
            &long,
            &h,
            "--bedgraph",
            "--empty",
            "0",
            "-g",
            &g,
        ])
        .assert()
        .code(1)
        .stderr(contains("extends past the end of chr1"));
}

#[test]
fn unsorted_input_and_sort_flag() {
    let dir = tempdir();
    let u = s(&plain_file(
        dir.path(),
        "u.bed",
        "chr2\t5\t15\nchr1\t20\t30\nchr1\t0\t10\n",
    ));
    let b = s(&plain_file(dir.path(), "b.bed", "chr1\t5\t25\n"));
    gloci()
        .args(["partition", &u, &b])
        .assert()
        .code(1)
        .stderr(contains("--sort"));
    let sorted = s(&plain_file(
        dir.path(),
        "s.bed",
        "chr1\t0\t10\nchr1\t20\t30\nchr2\t5\t15\n",
    ));
    // Each input is sorted on its own: labels stay with their file.
    let via_sort = out(&["partition", &u, &b, "--sort", "--names", "u,b"]);
    let direct = out(&["partition", &sorted, &b, "--names", "u,b"]);
    assert_eq!(via_sort, direct);
    assert!(via_sort.contains("chr1\t5\t10\t2\tu,b\n"), "{via_sort}");
}

#[test]
fn determinism_across_threads() {
    let t = three();
    let run = |n: &str| out(&["partition", &t.a, &t.b, &t.c, "--merge-adjacent", "-t", n]);
    let base = run("1");
    assert_eq!(run("2"), base);
    assert_eq!(run("8"), base);
}
