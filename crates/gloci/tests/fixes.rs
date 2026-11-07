//! Regression tests for the final-review fix wave (groups A-G).
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::Path;

fn w(d: &Path, name: &str, text: &str) -> String {
    plain_file(d, name, text).to_string_lossy().into_owned()
}

fn p(d: &Path, name: &str) -> String {
    d.join(name).to_string_lossy().into_owned()
}

// ---- Group A: sort/view format handling ----

const BG: &str = "chr1\t10\t20\t0.5\nchr1\t30\t40\t1.5\n";

#[test]
fn trackless_bedgraph_to_bedgraph_output() {
    let d = tempdir();
    let i = w(d.path(), "sig.bedgraph", BG);
    for cmd in ["sort", "view"] {
        let o = p(d.path(), &format!("{cmd}.out.bedgraph"));
        gloci().args([cmd, &i, "-o", &o]).assert().success();
        assert_eq!(std::fs::read_to_string(&o).unwrap(), BG, "{cmd}");
    }
}

#[test]
fn bed_input_to_bedgraph_output_and_back() {
    let d = tempdir();
    let i = w(d.path(), "x.bed", "chr1\t1\t2\n");
    let o = p(d.path(), "x.bedgraph");
    gloci().args(["sort", &i, "-o", &o]).assert().success();
    let bg = w(d.path(), "y.bedgraph", "chr1\t1\t2\t3\n");
    let o2 = p(d.path(), "y.bed");
    gloci().args(["view", &bg, "-o", &o2]).assert().success();
}

#[test]
fn sort_accepts_format_flag() {
    let d = tempdir();
    let i = w(
        d.path(),
        "sig.txt",
        "track type=bedGraph\nchr1\t30\t40\t1\nchr1\t10\t20\t2\n",
    );
    gloci()
        .args(["sort", &i, "--format", "bedgraph", "--header"])
        .assert()
        .success()
        .stdout("track type=bedGraph\nchr1\t10\t20\t2\nchr1\t30\t40\t1\n");
    gloci()
        .args(["sort", &i, "--format", "nope"])
        .assert()
        .code(2);
}

#[test]
fn family_mismatch_message_is_identical_for_sort_and_view() {
    let d = tempdir();
    let v = w(
        d.path(),
        "a.vcf",
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    );
    let o = p(d.path(), "x.bed");
    let msg = |cmd: &str| {
        let out = gloci()
            .args([cmd, &v, "-o", &o])
            .assert()
            .code(2)
            .get_output()
            .stderr
            .clone();
        String::from_utf8(out).unwrap()
    };
    let (s, vw) = (msg("sort"), msg("view"));
    assert_eq!(s, vw);
    assert!(s.contains("the input is VCF but the output format is BED"), "{s}");
    assert!(s.contains("-O"), "{s}");
    assert!(!Path::new(&o).exists());
}

// ---- Group B: undeclared VCF contigs ----

const VCF_COLS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
/// The header declares chr1 only; the records use chr2 as well.
fn undeclared_vcf() -> String {
    format!(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n{VCF_COLS}\
chr1\t5\t.\tA\tC\t.\t.\t.\nchr2\t7\t.\tA\tC\t.\t.\t.\nchr2\t9\t.\tG\tT\t.\t.\t.\n"
    )
}
const WARN2: &str = "chromosome \"chr2\" is not in the ##contig lines; add ##contig=<ID=chr2> to the header";
const SORT_HINT: &str =
    "add ##contig=<ID=chr2> to the header, or pass -g with a genome that lists it";

fn records(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

#[test]
fn view_plain_and_recompress_accept_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.vcf.gz");
    gloci()
        .args(["view", &x, "-o", &o, "--write-index"])
        .assert()
        .success()
        .stderr(contains(WARN2).and(contains("x.vcf line 5")))
        // one warning per contig, not per record
        .stderr(predicate::function(|s: &str| s.matches("warning:").count() == 1));
    assert!(d.path().join("x.vcf.gz.tbi").exists());
    let out = gloci().args(["view", &o]).assert().success().get_output().clone();
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(records(&text).len(), 3);
}

#[test]
fn view_quiet_silences_the_undeclared_warning() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    gloci()
        .args(["-q", "view", &x])
        .assert()
        .success()
        .stderr("");
}

#[test]
fn view_region_on_undeclared_contig_via_index_and_scan() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.vcf.gz");
    gloci().args(["view", &x, "-o", &o]).assert().success();
    gloci().args(["index", &o]).assert().success();
    for file in [&o, &x] {
        let out = gloci()
            .args(["view", file, "chr2:8-9"])
            .assert()
            .success()
            .stderr(contains(WARN2))
            .get_output()
            .clone();
        let text = String::from_utf8(out.stdout).unwrap();
        let recs = records(&text);
        assert_eq!(recs.len(), 1, "{file}: {text}");
        assert!(recs[0].starts_with("chr2\t9\t"), "{recs:?}");
    }
    // the whole contig, and a declared one alongside it
    let out = gloci()
        .args(["view", &o, "chr1", "chr2"])
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(records(&String::from_utf8(out.stdout).unwrap()).len(), 3);
}

#[test]
fn view_assume_sorted_accepts_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let out = gloci()
        .args(["view", &x, "--assume-sorted", "chr2"])
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(records(&String::from_utf8(out.stdout).unwrap()).len(), 2);
}

#[test]
fn index_accepts_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.vcf.gz");
    gloci().args(["view", &x, "-o", &o]).assert().success();
    gloci().args(["index", &o]).assert().success();
    assert!(d.path().join("x.vcf.gz.tbi").exists());
}

#[test]
fn sort_check_and_merge_refuse_undeclared_contigs_with_a_hint() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    for args in [
        vec!["sort", x.as_str()],
        vec!["sort", "--check", x.as_str()],
        vec!["sort", "-m", x.as_str()],
    ] {
        gloci()
            .args(&args)
            .assert()
            .code(1)
            .stderr(contains(SORT_HINT).and(contains("x.vcf line 5")));
    }
    let o = p(d.path(), "s.vcf");
    gloci().args(["sort", &x, "-o", &o]).assert().code(1);
    assert!(!Path::new(&o).exists());
}

#[test]
fn sort_accepts_an_undeclared_contig_that_g_lists() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let g = w(d.path(), "g.genome", "chr2\t1000\nchr1\t1000\n");
    let out = gloci()
        .args(["sort", &x, "-g", &g])
        .assert()
        .success()
        .get_output()
        .clone();
    let text = String::from_utf8(out.stdout).unwrap();
    let recs = records(&text);
    assert!(recs[0].starts_with("chr2\t7") && recs[2].starts_with("chr1\t5"), "{recs:?}");
}

#[test]
fn bcf_output_refuses_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.bcf");
    gloci()
        .args(["view", &x, "-o", &o])
        .assert()
        .code(1)
        .stderr(contains("BCF needs every contig").and(contains("add ##contig=<ID=chr2>")));
    assert!(!Path::new(&o).exists());
    // a declared file still converts
    let ok = w(
        d.path(),
        "ok.vcf",
        &format!("##fileformat=VCFv4.3\n##contig=<ID=chr1>\n{VCF_COLS}chr1\t5\t.\tA\tC\t.\t.\t.\n"),
    );
    gloci().args(["view", &ok, "-o", &o]).assert().success();
}

// ---- Group C: multi-VCF sort header definitions ----

fn vcf_with(defs: &str, recs: &[&str]) -> String {
    let mut s = format!("##fileformat=VCFv4.3\n{defs}##contig=<ID=chr1,length=1000>\n{VCF_COLS}");
    for r in recs {
        s.push_str(r);
        s.push('\n');
    }
    s
}

const DP: &str = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
const REC1: &str = "chr1\t5\t.\tA\tC\t.\t.\tDP=3";
const REC2: &str = "chr1\t2\t.\tA\tC\t.\t.\t.";

#[test]
fn multi_vcf_sort_refuses_a_definition_only_one_file_has() {
    let d = tempdir();
    let a = w(d.path(), "a.vcf", &vcf_with("", &[REC2]));
    let b = w(d.path(), "b.vcf", &vcf_with(DP, &[REC1]));
    let (an, bn) = (a.as_str(), b.as_str());
    gloci()
        .args(["sort", an, bn])
        .assert()
        .code(1)
        .stderr(contains(format!(
            "{b} defines ##INFO=<ID=DP> but {a} does not; give the files the same header definitions, or sort them separately"
        )));
    // the other way round names the file that has it too
    gloci()
        .args(["sort", bn, an])
        .assert()
        .code(1)
        .stderr(contains(format!("{b} defines ##INFO=<ID=DP> but {a} does not")));
    // -m is a multi-file path too, and --check does not merge anything
    gloci().args(["sort", "-m", an, bn]).assert().code(1);
    gloci().args(["sort", "--check", an, bn]).assert().success();
}

#[test]
fn multi_vcf_sort_refuses_conflicting_definitions() {
    let d = tempdir();
    let dp_float = "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth\">\n";
    let a = w(d.path(), "a.vcf", &vcf_with(DP, &[REC1]));
    let b = w(d.path(), "b.vcf", &vcf_with(dp_float, &[REC2]));
    gloci()
        .args(["sort", &a, &b])
        .assert()
        .code(1)
        .stderr(contains("define ##INFO=<ID=DP> differently").and(contains("Type=Float")));
}

#[test]
fn multi_vcf_sort_treats_pass_as_implicit_and_ignores_other_meta_lines() {
    let d = tempdir();
    let pass = "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    let a = w(d.path(), "a.vcf", &vcf_with(&format!("{pass}{DP}"), &[REC1]));
    let b = w(
        d.path(),
        "b.vcf",
        &vcf_with(&format!("##source=other\n{DP}"), &[REC2]),
    );
    let out = gloci()
        .args(["sort", &a, &b])
        .assert()
        .success()
        .get_output()
        .clone();
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(records(&text), vec![REC2, REC1]);
    // the first file's header is kept
    assert!(text.contains("##FILTER=<ID=PASS") && !text.contains("##source=other"));
}

#[test]
fn multi_vcf_sort_accepts_identical_headers() {
    let d = tempdir();
    let a = w(d.path(), "a.vcf", &vcf_with(DP, &[REC1]));
    let b = w(d.path(), "b.vcf", &vcf_with(DP, &[REC2]));
    gloci().args(["sort", &a, &b]).assert().success();
    gloci().args(["sort", "-m", &a, &b]).assert().success();
}
