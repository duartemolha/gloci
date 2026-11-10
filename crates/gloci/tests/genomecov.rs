//! `gloci genomecov`: every mode with exact output, scale, strand, split,
//! 5'/3' ends, max depth, empty chromosomes, errors and determinism.
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

const A: &str = "chr1\t0\t30\nchr1\t10\t20\nchr1\t10\t20\nchr1\t40\t50\nchr2\t5\t8\nchr4\t0\t3\n";
const G: &str = "chr1\t100\nchr2\t20\nchr3\t10\nchr4\t12\n";

struct Fx {
    _dir: tempfile::TempDir,
    a: String,
    g: String,
}

fn fx() -> Fx {
    let dir = tempdir();
    let a = s(&plain_file(dir.path(), "a.bed", A));
    let g = s(&plain_file(dir.path(), "g.genome", G));
    Fx { _dir: dir, a, g }
}

#[test]
fn hist_has_chromosome_and_genome_rows() {
    let f = fx();
    assert_eq!(
        out(&["genomecov", &f.a, "-g", &f.g]),
        "chr1\t0\t60\t100\t0.6\nchr1\t1\t30\t100\t0.3\nchr1\t3\t10\t100\t0.1\n\
chr2\t0\t17\t20\t0.85\nchr2\t1\t3\t20\t0.15\n\
chr3\t0\t10\t10\t1\n\
chr4\t0\t9\t12\t0.75\nchr4\t1\t3\t12\t0.25\n\
genome\t0\t96\t142\t0.676056\ngenome\t1\t36\t142\t0.253521\ngenome\t3\t10\t142\t0.0704225\n"
    );
}

#[test]
fn bedgraph_and_bga() {
    let f = fx();
    assert_eq!(
        out(&["genomecov", &f.a, "-g", &f.g, "--bedgraph"]),
        "chr1\t0\t10\t1\nchr1\t10\t20\t3\nchr1\t20\t30\t1\nchr1\t40\t50\t1\n\
chr2\t5\t8\t1\nchr4\t0\t3\t1\n"
    );
    assert_eq!(
        out(&["genomecov", &f.a, "-g", &f.g, "--bga"]),
        "chr1\t0\t10\t1\nchr1\t10\t20\t3\nchr1\t20\t30\t1\nchr1\t30\t40\t0\n\
chr1\t40\t50\t1\nchr1\t50\t100\t0\nchr2\t0\t5\t0\nchr2\t5\t8\t1\nchr2\t8\t20\t0\n\
chr3\t0\t10\t0\nchr4\t0\t3\t1\nchr4\t3\t12\t0\n"
    );
}

#[test]
fn per_base_is_one_based_with_zero_depth_bases() {
    let f = fx();
    let o = out(&["genomecov", &f.a, "-g", &f.g, "--per-base"]);
    let lines: Vec<&str> = o.lines().collect();
    assert_eq!(lines.len(), 142);
    assert_eq!(lines[0], "chr1\t1\t1");
    assert_eq!(lines[10], "chr1\t11\t3");
    assert_eq!(lines[30], "chr1\t31\t0");
    assert_eq!(lines[99], "chr1\t100\t0");
    assert_eq!(lines[100], "chr2\t1\t0");
    assert_eq!(lines[141], "chr4\t12\t0");
}

#[test]
fn scale_formats_like_bedtools() {
    let f = fx();
    let o = out(&[
        "genomecov",
        &f.a,
        "-g",
        &f.g,
        "--bedgraph",
        "--scale",
        "2.5",
    ]);
    assert!(
        o.starts_with("chr1\t0\t10\t2.5\nchr1\t10\t20\t7.5\n"),
        "{o}"
    );
    let o = out(&[
        "genomecov",
        &f.a,
        "-g",
        &f.g,
        "--bedgraph",
        "--scale",
        "1e-7",
    ]);
    assert!(
        o.starts_with("chr1\t0\t10\t1e-07\nchr1\t10\t20\t3e-07\n"),
        "{o}"
    );
    gloci()
        .args(["genomecov", &f.a, "-g", &f.g, "--scale", "2"])
        .assert()
        .code(2)
        .stderr(contains("--scale needs --bedgraph"));
}

#[test]
fn max_depth_caps_into_the_last_bin() {
    let f = fx();
    let o = out(&["genomecov", &f.a, "-g", &f.g, "--max-depth", "2"]);
    assert!(
        o.starts_with("chr1\t0\t60\t100\t0.6\nchr1\t1\t30\t100\t0.3\nchr1\t2\t10\t100\t0.1\n"),
        "{o}"
    );
    assert!(o.contains("genome\t2\t10\t142\t0.0704225\n"), "{o}");
    gloci()
        .args(["genomecov", &f.a, "-g", &f.g, "--bga", "--max-depth", "2"])
        .assert()
        .code(2)
        .stderr(contains("--max-depth only applies to the histogram"));
}

#[test]
fn strand_five_and_three_prime() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", G));
    let b = s(&plain_file(
        dir.path(),
        "s.bed",
        "chr1\t0\t30\tn\t0\t+\nchr1\t10\t20\tn\t0\t-\nchr1\t12\t14\tn\t0\t.\nchr1\t40\t50\tn\t0\t+\n",
    ));
    assert_eq!(
        out(&["genomecov", &b, "-g", &g, "--bedgraph", "--strand", "-"]),
        "chr1\t10\t20\t1\n"
    );
    assert_eq!(
        out(&["genomecov", &b, "-g", &g, "--bedgraph", "--strand", "+"]),
        "chr1\t0\t30\t1\nchr1\t40\t50\t1\n"
    );
    // An unknown strand counts as plus for the ends.
    assert_eq!(
        out(&["genomecov", &b, "-g", &g, "--bedgraph", "--5p"]),
        "chr1\t0\t1\t1\nchr1\t12\t13\t1\nchr1\t19\t20\t1\nchr1\t40\t41\t1\n"
    );
    assert_eq!(
        out(&["genomecov", &b, "-g", &g, "--bedgraph", "--3p"]),
        "chr1\t10\t11\t1\nchr1\t13\t14\t1\nchr1\t29\t30\t1\nchr1\t49\t50\t1\n"
    );
    gloci()
        .args(["genomecov", &b, "-g", &g, "--5p", "--3p"])
        .assert()
        .code(2);
}

#[test]
fn strand_needs_column_six() {
    let f = fx();
    gloci()
        .args(["genomecov", &f.a, "-g", &f.g, "--strand", "+"])
        .assert()
        .code(1)
        .stderr(contains("a.bed line 1"))
        .stderr(contains("column 6"));
}

#[test]
fn split_counts_bed12_blocks() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", G));
    let b = s(&plain_file(
        dir.path(),
        "b12.bed",
        "chr1\t0\t30\tn1\t0\t+\t0\t30\t0\t2\t5,5\t0,25\nchr1\t10\t40\tn2\t0\t-\t10\t40\t0\t2\t10,10\t0,20\n",
    ));
    assert_eq!(
        out(&["genomecov", &b, "-g", &g, "--bedgraph", "--split"]),
        "chr1\t0\t5\t1\nchr1\t10\t20\t1\nchr1\t25\t40\t1\n"
    );
    assert_eq!(
        out(&["genomecov", &b, "-g", &g, "--bedgraph"]),
        "chr1\t0\t10\t1\nchr1\t10\t30\t2\nchr1\t30\t40\t1\n"
    );
}

#[test]
fn empty_input_gives_every_chromosome_at_depth_zero() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", G));
    let e = s(&plain_file(dir.path(), "e.bed", ""));
    assert_eq!(
        out(&["genomecov", &e, "-g", &g]),
        "chr1\t0\t100\t100\t1\nchr2\t0\t20\t20\t1\nchr3\t0\t10\t10\t1\nchr4\t0\t12\t12\t1\n\
genome\t0\t142\t142\t1\n"
    );
}

#[test]
fn errors_say_what_and_where() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", G));
    let x = s(&plain_file(
        dir.path(),
        "x.bed",
        "chr1\t0\t10\nchrX\t5\t8\n",
    ));
    gloci()
        .args(["genomecov", &x, "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("chrX"))
        .stderr(contains("x.bed line 2"));
    let p = s(&plain_file(dir.path(), "p.bed", "chr2\t10\t30\n"));
    gloci()
        .args(["genomecov", &p, "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("p.bed line 1"))
        .stderr(contains("extends past the end of chr2 (length 20)"));
    gloci().args(["genomecov", &p]).assert().code(2);
}

#[test]
fn unsorted_input_and_sort_flag() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", G));
    let u = s(&plain_file(
        dir.path(),
        "u.bed",
        "chr1\t20\t30\nchr1\t0\t10\n",
    ));
    gloci()
        .args(["genomecov", &u, "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("--sort"));
    let o = out(&["genomecov", &u, "-g", &g, "--bedgraph", "--sort"]);
    assert_eq!(o, "chr1\t0\t10\t1\nchr1\t20\t30\t1\n");
}

#[test]
fn points_are_ignored() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", G));
    let p = s(&plain_file(
        dir.path(),
        "pt.bed",
        "chr1\t5\t5\nchr1\t7\t9\n",
    ));
    assert_eq!(
        out(&["genomecov", &p, "-g", &g, "--bedgraph"]),
        "chr1\t7\t9\t1\n"
    );
}

#[test]
fn determinism_across_threads() {
    let f = fx();
    for mode in ["--hist", "--bga", "--per-base"] {
        let one = out(&["genomecov", &f.a, "-g", &f.g, mode, "--threads", "1"]);
        for t in ["2", "4", "8"] {
            assert_eq!(
                out(&["genomecov", &f.a, "-g", &f.g, mode, "--threads", t]),
                one,
                "{mode} --threads {t}"
            );
        }
    }
}
