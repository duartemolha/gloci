//! Tests for `gloci vcf filter --include/--exclude` (decision D71).
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::str::contains;
use std::path::Path;

const HEAD: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
";

/// Rows 100..500. MAF over all four samples: 100 -> 0.5, 200 -> 0.125 (7 ref, 1 alt),
/// 300 -> 0.5, 400 -> 0.25 (6 ref, 2 alt), 500 -> 0.5. On s3,s4 alone: 100 -> 0,
/// 200 -> 0.25, 300 -> 0.25, 400 -> 0.25, 500 -> 0.
const ROWS: &str = "\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=40;AF=0.5\tGT:GQ\t0/1:30\t0/1:30\t0/0:30\t0/0:30\n\
chr1\t200\trs2\tG\tT\t20\tPASS\tDP=16;AF=0.125\tGT:GQ\t0/0:20\t0/0:20\t0/0:20\t0/1:20\n\
chr1\t300\trs3\tC\tG\t40\tPASS\tDP=30;AF=0.5\tGT:GQ\t0/1:25\t0/0:25\t0/1:25\t0/0:25\n\
chr1\t400\trs4\tA\tT\t60\tPASS\tDP=80;AF=0.25\tGT:GQ\t0/0:40\t0/0:40\t0/1:5\t1/1:40\n\
chr1\t500\trs5\tT\tA\t55\tPASS\tDP=20;AF=0.5\tGT:GQ\t0/1:35\t0/1:35\t0/0:35\t0/0:35\n";

fn vcf(d: &Path) -> String {
    plain_file(d, "in.vcf", &format!("{HEAD}{ROWS}"))
        .to_string_lossy()
        .into_owned()
}

fn positions(args: &[&str]) -> Vec<u64> {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect()
}

fn filter(file: &str, extra: &[&str]) -> Vec<u64> {
    let mut a = vec!["vcf", "filter", file, "-O", "vcf"];
    a.extend_from_slice(extra);
    positions(&a)
}

#[test]
fn include_selects_and_exclude_negates() {
    let d = tempdir();
    let f = vcf(d.path());
    let e = "QUAL>=50 && INFO/DP>20";
    assert_eq!(filter(&f, &["--include", e]), vec![100, 400]);
    // --exclude keeps exactly the records --include drops.
    assert_eq!(filter(&f, &["--exclude", e]), vec![200, 300, 500]);
    assert_eq!(
        filter(&f, &["--exclude", "!(QUAL>=50 && INFO/DP>20)"]),
        vec![100, 400]
    );
}

#[test]
fn expression_is_anded_with_the_other_filters_before_thin() {
    let d = tempdir();
    let f = vcf(d.path());
    // QUAL>=40 holds at 100,300,400,500; --min-qual 50 drops 300.
    assert_eq!(
        filter(&f, &["--include", "QUAL>=40", "--min-qual", "50"]),
        vec![100, 400, 500]
    );
    // Thinning runs on the sites the expression passed: 200 is dropped by the
    // expression, so it cannot occupy the window of 300.
    assert_eq!(
        filter(&f, &["--include", "QUAL>=40", "--thin", "150"]),
        vec![100, 300, 500]
    );
}

#[test]
fn maf_sees_the_selected_samples() {
    let d = tempdir();
    let f = vcf(d.path());
    // All samples: MAF=0.125 only at 200.
    assert_eq!(filter(&f, &["--include", "MAF<0.2"]), vec![200]);
    // s3,s4 only: 100 and 500 have MAF 0, the others 0.25.
    assert_eq!(
        filter(&f, &["--samples", "s3,s4", "--include", "MAF<0.2"]),
        vec![100, 500]
    );
    // N_PASS counts only the selected samples.
    assert_eq!(
        filter(
            &f,
            &["--samples", "s3,s4", "--include", "N_PASS(GT=\"het\")=2"]
        ),
        Vec::<u64>::new()
    );
    assert_eq!(
        filter(
            &f,
            &["--samples", "s1,s2", "--include", "N_PASS(GT=\"het\")=2"]
        ),
        vec![100, 500]
    );
}

#[test]
fn expression_sees_masked_genotypes() {
    let d = tempdir();
    let f = vcf(d.path());
    // --min-gq 30 masks every genotype with GQ<30; 100 keeps s1,s2 het, s3,s4 hom.
    // 400 keeps s1,s2 (GQ 40), s4 (GQ 40); s3 (GQ 5) is masked.
    assert_eq!(
        filter(&f, &["--min-gq", "30", "--include", "N_MISSING>=3"]),
        vec![200, 300]
    );
    assert_eq!(
        filter(
            &f,
            &["--min-gq", "30", "--include", "GT=\"mis\" && N_MISSING=1"]
        ),
        vec![400]
    );
}

#[test]
fn include_and_exclude_together_are_a_usage_error() {
    let d = tempdir();
    let f = vcf(d.path());
    gloci()
        .args([
            "vcf",
            "filter",
            &f,
            "--include",
            "QUAL>1",
            "--exclude",
            "QUAL>2",
        ])
        .assert()
        .code(2);
}

#[test]
fn bad_expression_exits_2_before_reading_records() {
    let d = tempdir();
    // The only record is malformed: reading it would be a data error (exit 1).
    let head = HEAD.replace(
        "##FORMAT=<ID=GT",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"d\">\n##FORMAT=<ID=GT",
    );
    let f = plain_file(
        d.path(),
        "bad.vcf",
        &format!("{head}chr1\tnot-a-position\n"),
    )
    .to_string_lossy()
    .into_owned();
    let out = d.path().join("out.vcf");
    let out = out.to_string_lossy().into_owned();
    gloci()
        .args(["vcf", "filter", &f, "--include", "QUAL>>1", "-o", &out])
        .assert()
        .code(2)
        .stderr(contains("column"));
    assert!(
        !Path::new(&out).exists(),
        "no output before the expression is checked"
    );
    gloci()
        .args(["vcf", "filter", &f, "--exclude", "FMT/AD>1"])
        .assert()
        .code(2)
        .stderr(contains("FMT/AD is not supported"));
    // Sanity: without an expression the same input is a data error.
    gloci().args(["vcf", "filter", &f]).assert().code(1);
}

#[test]
fn bare_key_defined_as_info_and_format_is_ambiguous() {
    let d = tempdir();
    let head = HEAD.replace(
        "##FORMAT=<ID=GT",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n##FORMAT=<ID=GT",
    );
    let f = plain_file(d.path(), "amb.vcf", &format!("{head}{ROWS}"))
        .to_string_lossy()
        .into_owned();
    gloci()
        .args(["vcf", "filter", &f, "--include", "DP>10"])
        .assert()
        .code(2)
        .stderr(contains("ambiguous"))
        .stderr(contains("write INFO/DP or FMT/DP"));
}

#[test]
fn bcf_input_matches_text_input() {
    let d = tempdir();
    let text = format!("{HEAD}{ROWS}");
    let bcf = d.path().join("in.bcf");
    std::fs::write(&bcf, fixtures::vcf_to_bcf(text.as_bytes())).unwrap();
    let b = bcf.to_string_lossy().into_owned();
    let f = vcf(d.path());
    for e in [
        "INFO/DP>20",
        "INFO/AF[0]=0.5 && FMT/GQ>=30",
        "MAF<0.3",
        "N_PASS(GT=\"het\")>=2",
        "TYPE=\"snp\" && QUAL>45",
    ] {
        for flag in ["--include", "--exclude"] {
            assert_eq!(filter(&b, &[flag, e]), filter(&f, &[flag, e]), "{flag} {e}");
        }
    }
    assert_eq!(
        filter(&b, &["--include", "INFO/DP>20"]),
        vec![100, 300, 400]
    );
}

#[test]
fn deterministic_across_thread_counts() {
    let d = tempdir();
    let f = vcf(d.path());
    let run = |t: &str| {
        let o = gloci()
            .args([
                "vcf",
                "filter",
                &f,
                "-t",
                t,
                "--include",
                "MAF>0.1 && N_PASS(FMT/GQ>=30)>=1",
            ])
            .assert()
            .success();
        o.get_output().stdout.clone()
    };
    let one = run("1");
    assert_eq!(one, run("2"));
    assert_eq!(one, run("4"));
}

#[test]
fn deeply_nested_expression_exits_2_instead_of_crashing() {
    let d = tempdir();
    let f = plain_file(d.path(), "in.vcf", &format!("{HEAD}{ROWS}"))
        .to_string_lossy()
        .into_owned();
    let out = d.path().join("out.vcf");
    let out = out.to_string_lossy().into_owned();
    let parens = format!("{}QUAL>1{}", "(".repeat(5000), ")".repeat(5000));
    let nots = format!("{}QUAL>1", "!".repeat(20000));
    for expr in [parens, nots] {
        gloci()
            .args(["vcf", "filter", &f, "--include", &expr, "-o", &out])
            .assert()
            .code(2)
            .stderr(contains("expression nested too deeply"))
            .stderr(contains("column"));
    }
    assert!(!Path::new(&out).exists());
}
