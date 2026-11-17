//! Tests for `gloci vcf query` (decision D172): bcftools 1.21 query's format language over the
//! shared selection stage. The parity with bcftools is the goldens' (`vcf-query-*`).
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn stdout_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn expr_format() -> String {
    s(&diff_inputs().join("expr-format.vcf"))
}

#[test]
fn query_prints_the_format() {
    let out = stdout_of(&[
        "vcf",
        "query",
        &expr_format(),
        "-f",
        "%CHROM:%POS[\\t%SAMPLE=%GT]\\n",
    ]);
    assert_eq!(
        out,
        "1:10\ts1=0/1\ts2=0/0\ts3=1/1\n1:20\ts1=./.\ts2=0|1\ts3=1\n\
         1:30\ts1=0/0\ts2=0/0\ts3=./.\n1:40\ts1=0/1\ts2=1/1\ts3=./1\n"
    );
}

#[test]
fn query_vcf_and_bcf_agree() {
    let dir = tempdir();
    let fixed =
        "%CHROM\\t%POS\\t%END\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%TYPE[\\t%SAMPLE=%GT:%TGT]";
    let cases = [
        (
            "vcf/split/multi.vcf",
            "multi",
            "\\t%DP\\t%AC\\t%AF{1}[\\t%AD:%PL:%GL:%AQ:%DP]\\n",
        ),
        (
            "vcf/query/numbers.vcf",
            "numbers",
            "\\t%I\\t%F[\\t%GL:%FI:%FF]\\n",
        ),
        // Integers at and just outside the 32-bit range (D172).
        (
            "vcf/query/int-range.vcf",
            "int-range",
            "\\t%I\\t%L[\\t%FI]\\n",
        ),
    ];
    for (rel, name, tail) in cases {
        let vcf = s(&diff_inputs().join(rel));
        let bcf = s(&dir.path().join(format!("{name}.bcf")));
        gloci()
            .args(["view", &vcf, "--vcf-span", "ref", "-o", &bcf])
            .assert()
            .success();
        let f = format!("{fixed}{tail}");
        let from_vcf = stdout_of(&["vcf", "query", &vcf, "-f", &f]);
        let from_bcf = stdout_of(&["vcf", "query", &bcf, "-f", &f]);
        assert!(!from_vcf.is_empty());
        assert_eq!(from_vcf, from_bcf, "{name}");
    }
}

#[test]
fn query_samples_subset_with_a_sample_filter() {
    let out = stdout_of(&[
        "vcf",
        "query",
        &expr_format(),
        "--samples",
        "s2,s3",
        "--include",
        "GT=\"het\"",
        "-f",
        "%POS[\\t%SAMPLE]\\n",
    ]);
    assert_eq!(out, "20\ts2\n");
}

#[test]
fn query_masked_samples_fail_a_sample_filter() {
    let out = stdout_of(&[
        "vcf",
        "query",
        &expr_format(),
        "--min-dp",
        "20",
        "--include",
        "GT=\"het\"",
        "-f",
        "%POS[\\t%SAMPLE=%GT:%DP]\\n",
    ]);
    assert_eq!(out, "10\ts1=0/1:25\n40\ts1=0/1:.\n");
    let out = stdout_of(&[
        "vcf",
        "query",
        &expr_format(),
        "--min-dp",
        "20",
        "-f",
        "[%GT ]\\n",
    ]);
    assert_eq!(out.lines().next().unwrap(), "0/1 ./. 1/1 ");
}

#[test]
fn query_bad_number_is_a_data_error() {
    let dir = tempdir();
    let vcf = plain_file(
        dir.path(),
        "bad.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
         ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\t.\tA\tC\t.\t.\tDP=abc\n",
    );
    gloci()
        .args(["vcf", "query", &s(&vcf), "-f", "%DP\\n"])
        .assert()
        .code(1)
        .stderr(contains("record 1"))
        .stderr(contains("INFO/DP"));
}

#[test]
fn query_sites_only_input() {
    let dir = tempdir();
    let vcf = plain_file(
        dir.path(),
        "sites.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t100\t.\tA\tC\t.\t.\t.\n",
    );
    let v = s(&vcf);
    assert_eq!(
        stdout_of(&["vcf", "query", &v, "-f", "%POS[\\t%GT]\\n"]),
        "100\n"
    );
    assert_eq!(
        stdout_of(&["vcf", "query", &v, "-H", "-f", "%POS[\\t%GT]\\n"]),
        "#[1]POS\n100\n"
    );
}

#[test]
fn query_output_file_gz_and_refused_names() {
    let dir = tempdir();
    let input = expr_format();
    let f = "%POS[\\t%GT]\\n";
    let plain = stdout_of(&["vcf", "query", &input, "-f", f]);

    let gz = dir.path().join("out.txt.gz");
    gloci()
        .args(["vcf", "query", &input, "-f", f, "-o", &s(&gz)])
        .assert()
        .success();
    let raw = std::fs::read(&gz).unwrap();
    assert_eq!(&raw[..2], &[0x1f, 0x8b]);
    // The BGZF extra field: FEXTRA set, subfield id `BC`.
    assert_eq!(raw[3] & 4, 4);
    assert_eq!(&raw[12..14], b"BC");
    let mut text = String::new();
    std::io::Read::read_to_string(&mut noodles_bgzf::io::Reader::new(&raw[..]), &mut text).unwrap();
    assert_eq!(text, plain);

    let txt = dir.path().join("out.txt");
    gloci()
        .args(["vcf", "query", &input, "-f", f, "-o", &s(&txt)])
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&txt).unwrap(), plain);

    gloci()
        .args([
            "vcf",
            "query",
            &input,
            "-f",
            f,
            "-o",
            &s(&dir.path().join("out.bcf")),
        ])
        .assert()
        .code(2)
        .stderr(contains("writes text"));
}

#[test]
fn query_format_errors_exit_2() {
    gloci()
        .args(["vcf", "query", "no-such-file.vcf", "-f", "%POS %"])
        .assert()
        .code(2)
        .stderr(contains("column 6"));
    let input = expr_format();
    for (args, needle) in [
        (vec!["-f", "%SAMPLE\\n"], "inside [ ]"),
        (vec!["-H", "-f", "%LINE"], "%LINE"),
        (
            vec!["--include", "GT=\"het\" & FMT/DP>5", "-f", "%POS\\n"],
            "&",
        ),
    ] {
        gloci()
            .args(["vcf", "query", &input])
            .args(&args)
            .assert()
            .code(2)
            .stderr(contains(needle));
    }
}

#[test]
fn query_help_lists_the_fields() {
    let out = stdout_of(&["vcf", "query", "--help"]);
    // Help wraps lines; compare the words.
    let words = out.split_whitespace().collect::<Vec<_>>().join(" ");
    for needle in [
        "bcftools query",
        "%FIRST_ALT",
        "%TGT",
        "%LINE",
        "Experimental",
        "They print in the order given here",
    ] {
        assert!(words.contains(needle), "missing {needle:?}");
    }
    assert!(
        !words.contains("header's order"),
        "the shared text must not appear"
    );
    assert!(out.contains("-f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'"));
}

#[test]
fn query_threads_do_not_change_the_output() {
    let input = s(&diff_inputs().join("filter-cohort.vcf"));
    let run = |t: &str| {
        stdout_of(&[
            "vcf",
            "query",
            &input,
            "-t",
            t,
            "--exclude",
            "GT=\"alt\"",
            "-f",
            "%POS[\\t%SAMPLE]\\n",
        ])
    };
    let one = run("1");
    assert!(!one.is_empty());
    assert_eq!(one, run("2"));
    assert_eq!(one, run("8"));
}

#[test]
fn query_out_of_range_integers_print_missing_and_warn_once_per_tag() {
    let input = s(&diff_inputs().join("vcf/query/int-range.vcf"));
    let o = gloci()
        .args(["vcf", "query", &input, "-f", "%POS\\t%I\\t%L[\\t%FI]\\n"])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out,
        "5\t.\t1,.,3\t.\n6\t2147483647\t2147483647\t2147483647\n7\t.\t.,4\t.\n\
         8\t.\t.,.\t.\n9\t.\t7\t.\n10\t.\t.,.,5\t.\n"
    );
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    for tag in [
        "INFO/I value 2147483648",
        "INFO/L value 2147483648",
        "FORMAT/FI value 2147483648",
    ] {
        assert_eq!(err.matches(tag).count(), 1, "{tag} in {err}");
    }
    assert!(err.contains("outside the 32-bit range of VCF integers; printed as missing"));
}
