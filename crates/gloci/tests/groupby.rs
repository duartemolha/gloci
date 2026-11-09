//! `gloci groupby`: consecutive grouping, key columns, header names, every
//! scalar op, usage and data errors, determinism.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

struct Fx {
    dir: tempfile::TempDir,
}

impl Fx {
    fn new() -> Fx {
        Fx { dir: tempdir() }
    }
    fn path(&self, text: &str) -> String {
        plain_file(self.dir.path(), "in.tsv", text)
            .to_string_lossy()
            .into_owned()
    }
    /// `groupby FILE args...` over `text`, as the output text.
    fn run(&self, text: &str, args: &[&str]) -> String {
        let f = self.path(text);
        let o = gloci().arg("groupby").arg(&f).args(args).assert().success();
        String::from_utf8(o.get_output().stdout.clone()).unwrap()
    }
}

const UNSORTED: &str = "a\t1\tx\t5\na\t1\ty\t7\nb\t2\tz\t2.5\na\t1\tw\t1\n";

#[test]
fn consecutive_runs_not_whole_keys() {
    let fx = Fx::new();
    // `a` comes back after `b`, so it is a second group, as in bedtools groupby
    assert_eq!(
        fx.run(UNSORTED, &["-g", "1", "-c", "4:sum"]),
        "a\t12\nb\t2.5\na\t1\n"
    );
}

#[test]
fn several_key_columns_and_ops() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(
            UNSORTED,
            &["-g", "1,2", "-c", "4:sum,4:mean", "-c", "3:collapse"]
        ),
        "a\t1\t12\t6\tx,y\nb\t2\t2.5\t2.5\tz\na\t1\t1\t1\tw\n"
    );
    // the key columns come out in -g order
    assert_eq!(
        fx.run(UNSORTED, &["-g", "2,1", "-c", "0:count"]),
        "1\ta\t2\n2\tb\t1\n1\ta\t1\n"
    );
}

#[test]
fn keys_keep_the_bytes_of_the_first_line() {
    let fx = Fx::new();
    // equal bytes only: `1.0` and `1` are different keys
    assert_eq!(
        fx.run("1.0\t2\n1.0\t3\n1\t4\n", &["-g", "1", "-c", "2:sum"]),
        "1.0\t5\n1\t4\n"
    );
}

#[test]
fn every_scalar_op() {
    let fx = Fx::new();
    let text = "k\t3\nk\t1\nk\t2\nk\t2\n";
    for (op, want) in [
        ("count", "4"),
        ("count_distinct", "3"),
        ("sum", "8"),
        ("mean", "2"),
        ("median", "2"),
        ("min", "1"),
        ("max", "3"),
        ("var", "0.6666666666666666"),
        ("stdev", "0.816496580927726"),
        ("var_pop", "0.5"),
        ("stdev_pop", "0.7071067811865476"),
        ("collapse", "3,1,2,2"),
        ("distinct", "3,1,2"),
        ("distinct_sorted", "1,2,3"),
        ("first", "3"),
        ("last", "2"),
        ("indicator", "1"),
        ("mode", "2"),
        ("antimode", "1"),
        ("mad", "0.5"),
        ("quantile:0.5", "2"),
        ("tmean:0-1", "2"),
        ("cv", "0.408248290463863"),
    ] {
        assert_eq!(
            fx.run(text, &["-g", "1", "-c", &format!("2:{op}")]),
            format!("k\t{want}\n"),
            "{op}"
        );
    }
    // an empty aggregation prints --empty
    assert_eq!(
        fx.run("k\t.\n", &["-g", "1", "-c", "2:sum", "--empty", "NA"]),
        "k\tNA\n"
    );
    assert_eq!(
        fx.run(text, &["-g", "1", "-c", "2:collapse", "--delim", ";"]),
        "k\t3;1;2;2\n"
    );
}

#[test]
fn header_names_keys_and_columns() {
    let fx = Fx::new();
    let text = "#gene\tsample\texpr\ng1\ts1\t2\ng1\ts2\t4\ng2\ts1\t10\n";
    assert_eq!(
        fx.run(text, &["--header", "-g", "gene", "-c", "expr:sum"]),
        "#gene\texpr:sum\ng1\t6\ng2\t10\n"
    );
    // numbers work with --header too; a header without `#` stays without
    assert_eq!(
        fx.run(
            "gene\tn\ng1\t1\ng1\t2\n",
            &["--header", "-g", "1", "-c", "n:sum", "-c", "2:max"]
        ),
        "gene\tn:sum\t2:max\ng1\t3\t2\n"
    );
    // without --header the first line is data and `#` lines are skipped
    assert_eq!(
        fx.run(text, &["-g", "2", "-c", "3:sum"]),
        "s1\t2\ns2\t4\ns1\t10\n"
    );
}

#[test]
fn names_need_the_header_and_must_exist() {
    let fx = Fx::new();
    let f = fx.path("gene\tn\ng1\t1\n");
    gloci()
        .args(["groupby", &f, "-g", "gene", "-c", "2:sum"])
        .assert()
        .code(2)
        .stderr(contains("needs --header"));
    gloci()
        .args(["groupby", &f, "--header", "-g", "nope", "-c", "2:sum"])
        .assert()
        .code(2)
        .stderr(contains("no column named `nope`").and(contains("available names: gene, n")));
    gloci()
        .args(["groupby", &f, "--header", "-g", "gene", "-c", "nope:sum"])
        .assert()
        .code(2)
        .stderr(contains("no column named `nope`"));
}

#[test]
fn skips_blank_comment_and_track_lines_and_strips_cr() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(
            "#c\tn\nk\t1\r\n\nbrowser x\ntrack name=t\nk\t2\r\n",
            &["-g", "1", "-c", "2:collapse"]
        ),
        "k\t1,2\n"
    );
    assert_eq!(fx.run("", &["-g", "1", "-c", "2:sum"]), "");
}

#[test]
fn works_on_stdin_and_bed() {
    let o = gloci()
        .args(["groupby", "-", "-g", "1", "-c", "5:max"])
        .write_stdin("chr1\t1\t2\tn\t5\t+\nchr1\t3\t4\tm\t9\t+\n")
        .assert()
        .success();
    assert_eq!(String::from_utf8_lossy(&o.get_output().stdout), "chr1\t9\n");
}

#[test]
fn a_missing_key_column_is_an_error_naming_the_line() {
    let fx = Fx::new();
    let f = fx.path("a\t1\nb\n");
    gloci()
        .args(["groupby", &f, "-g", "2", "-c", "1:count"])
        .assert()
        .code(1)
        .stderr(contains("line 2").and(contains("-g needs column 2")));
}

#[test]
fn two_input_ops_are_refused() {
    let fx = Fx::new();
    let f = fx.path("a\t1\n");
    for op in [
        "wmean",
        "records",
        "min_record",
        "max_record",
        "frac_a",
        "frac_b",
        "bases_uniq",
        "bases",
    ] {
        gloci()
            .args(["groupby", &f, "-g", "1", "-c", &format!("2:{op}")])
            .assert()
            .code(2)
            .stderr(contains("groupby has none").and(contains("use map")));
    }
}

#[test]
fn vcf_input_is_a_usage_error() {
    let fx = Fx::new();
    let f = fx.path("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    gloci()
        .args(["groupby", &f, "-g", "1", "-c", "2:sum"])
        .assert()
        .code(2)
        .stderr(contains("VCF").and(contains("gloci vcf extract")));
    let g = bgzf_file(
        fx.dir.path(),
        "v.vcf.gz",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    gloci()
        .args(["groupby", &g.to_string_lossy(), "-g", "1", "-c", "2:sum"])
        .assert()
        .code(2)
        .stderr(contains("gloci vcf extract"));
}

#[test]
fn bgzf_input_and_output_file() {
    let fx = Fx::new();
    let f = bgzf_file(fx.dir.path(), "x.tsv.gz", UNSORTED);
    let out = fx.dir.path().join("o.tsv");
    gloci()
        .args([
            "groupby",
            &f.to_string_lossy(),
            "-g",
            "1",
            "-c",
            "4:sum",
            "-o",
        ])
        .arg(&out)
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(out).unwrap(),
        "a\t12\nb\t2.5\na\t1\n"
    );
}

#[test]
fn identical_for_every_thread_count() {
    let fx = Fx::new();
    let mut text = String::new();
    for i in 0..5000 {
        text.push_str(&format!("k{}\t{}\t{}\n", i / 7 % 13, i, (i * 31) % 17));
    }
    let args = |t: &'static str| {
        [
            "--threads",
            t,
            "-g",
            "1",
            "-c",
            "2:sum,3:median,3:mode,2:collapse",
        ]
    };
    let one = fx.run(&text, &args("1"));
    for t in ["2", "4"] {
        assert_eq!(fx.run(&text, &args(t)), one, "--threads {t}");
    }
}

#[test]
fn whitespace_only_lines_are_skipped() {
    let fx = Fx::new();
    assert_eq!(
        fx.run("k\t1\n\t\n  \nk\t2\n", &["-g", "1", "-c", "2:sum"]),
        "k\t3\n"
    );
}

#[test]
fn only_line_one_is_the_header() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(
            "#gene\tn\n#a comment\ng1\t1\n#another\ng1\t2\n",
            &["--header", "-g", "gene", "-c", "n:sum"]
        ),
        "#gene\tn:sum\ng1\t3\n"
    );
}

#[test]
fn an_absent_c_column_is_an_error_but_empty_and_dot_are_missing() {
    let fx = Fx::new();
    let f = fx.path("k\t1\nk\t2\nk\n");
    gloci()
        .args(["groupby", &f, "-g", "1", "-c", "2:sum"])
        .assert()
        .code(1)
        .stderr(
            contains("line 3")
                .and(contains("column 2"))
                .and(contains("line has 1 tab-separated column")),
        );
    let g = fx.path("k\t1\nk\t\nk\t.\n");
    gloci()
        .args(["groupby", &g, "-g", "1", "-c", "2:sum,2:count"])
        .assert()
        .success()
        .stdout("k\t1\t3\n");
}

#[test]
fn a_non_numeric_value_is_an_error_with_a_header_hint_on_the_first_line() {
    let fx = Fx::new();
    let f = fx.path("k\tn\nk\t2\n");
    gloci()
        .args(["groupby", &f, "-g", "1", "-c", "2:sum"])
        .assert()
        .code(1)
        .stderr(contains("line 1").and(contains("if the first line is a header, pass --header")));
    gloci()
        .args(["groupby", &f, "--header", "-g", "k", "-c", "n:sum"])
        .assert()
        .success()
        .stdout("k\tn:sum\nk\t2\n");
    let g = fx.path("k\t1\nk\tx\n");
    gloci()
        .args(["groupby", &g, "-g", "1", "-c", "2:sum"])
        .assert()
        .code(1)
        .stderr(contains("line 2").and(contains("--header").not()));
}
