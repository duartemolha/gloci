//! `--skip-invalid` reports what it skipped on every command: per input,
//! the count and the first examples.
mod common;
use common::*;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// Run `args` and return stderr, asserting success.
fn stderr_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stderr.clone()).unwrap()
}

fn assert_reports(cmd: &str, stderr: &str, files: &[&str]) {
    for f in files {
        assert!(
            stderr.contains(&format!("skipped 1 invalid records in {f}")),
            "{cmd}: no summary for {f}:\n{stderr}"
        );
    }
}

#[test]
fn every_command_reports_skipped_records() {
    let d = tempdir();
    let a = s(&plain_file(
        d.path(),
        "a.bed",
        "chr1\t10\t50\ta1\t0\t+\nchr1\tBAD\t60\tbad\t0\t+\nchr1\t100\t200\ta2\t0\t-\n",
    ));
    let b = s(&plain_file(
        d.path(),
        "b.bed",
        "chr1\t20\t30\tb1\t5\t+\nchr1\t40\tNOPE\tbad\t0\t+\nchr1\t150\t160\tb2\t7\t-\n",
    ));
    let g = s(&plain_file(d.path(), "g.genome", "chr1\t1000\n"));
    let two = |cmd: &str, extra: &[&str]| {
        let mut args = vec![cmd, "-a", a.as_str(), "-b", b.as_str(), "--skip-invalid"];
        args.extend_from_slice(extra);
        let e = stderr_of(&args);
        assert_reports(cmd, &e, &[&a, &b]);
    };
    two("intersect", &[]);
    two("intersect", &["--unsorted"]);
    two("intersect", &["--mode", "pairs", "--keep-unmatched"]);
    two("subtract", &[]);
    two("closest", &[]);
    two("map", &["-c", "5:sum"]);
    two("coverage", &[]);
    for cmd in ["merge", "complement"] {
        let e = stderr_of(&[cmd, &a, &b, "--skip-invalid"]);
        assert_reports(cmd, &e, &[&a, &b]);
    }
    let e = stderr_of(&["slop", &a, "-b", "5", "--skip-invalid"]);
    assert_reports("slop", &e, &[&a]);
    let e = stderr_of(&["flank", &a, "-b", "5", "-g", &g, "--skip-invalid"]);
    assert_reports("flank", &e, &[&a]);
    let e = stderr_of(&["windows", "-i", &a, "-w", "20", "--skip-invalid"]);
    assert_reports("windows", &e, &[&a]);
}

/// The summary names the bad line and is printed once.
#[test]
fn summary_has_examples_once() {
    let d = tempdir();
    let a = s(&plain_file(
        d.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\tx\t9\nchr1\t7\t9\n",
    ));
    let b = s(&plain_file(d.path(), "b.bed", "chr1\t2\t3\n"));
    let e = stderr_of(&["intersect", "-a", &a, "-b", &b, "--skip-invalid"]);
    assert_eq!(e.matches("skipped 1 invalid records").count(), 1, "{e}");
    assert!(e.contains("line 2:"), "{e}");
}
