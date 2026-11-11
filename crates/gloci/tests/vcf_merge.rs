//! `gloci vcf merge`: the rules that need no bcftools to check (the goldens
//! compare with it), determinism, BCF, errors, the file-handle limit and
//! the memory bound.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/vcf/merge")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn records(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(String::from)
        .collect()
}

fn merge3(extra: &[&str]) -> String {
    let (a, b, c) = (
        s(&fixture("a.vcf")),
        s(&fixture("b.vcf")),
        s(&fixture("c.vcf")),
    );
    let mut args = vec!["vcf", "merge", &a, &b, &c];
    args.extend_from_slice(extra);
    out(&args)
}

fn cols(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

const HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##contig=<ID=2,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"n\">\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"c\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n";

fn vcf(sample: &str, rows: &[&str]) -> String {
    let mut t = String::from(HEAD);
    t.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ));
    for r in rows {
        t.push_str(r);
        t.push('\n');
    }
    t
}

#[test]
fn header_is_the_union_and_samples_follow_the_files() {
    let text = merge3(&[]);
    let hdr: Vec<&str> = text.lines().filter(|l| l.starts_with('#')).collect();
    assert_eq!(hdr[0], "##fileformat=VCFv4.2");
    assert_eq!(
        hdr[1],
        "##FILTER=<ID=PASS,Description=\"All filters passed\">"
    );
    assert!(hdr.contains(&"##FILTER=<ID=q10,Description=\"Quality below 10\">"));
    let chrom = hdr.last().unwrap();
    assert!(chrom.ends_with("FORMAT\tA1\tA2\tB1\tC1\tC2"), "{chrom}");
    assert!(
        !text.contains("glociCommand"),
        "no command line, as --no-version"
    );
}

#[test]
fn alleles_are_united_in_file_order_and_values_reindexed() {
    let recs = records(&merge3(&[]));
    // 1:150: a has T,A, c has T,C, b has A: a and c pair (they share T).
    let r = recs
        .iter()
        .map(|l| cols(l).into_iter().map(String::from).collect::<Vec<_>>())
        .find(|c| c[1] == "150" && c[4] == "T,A,C")
        .expect("1:150 T,A,C");
    assert_eq!(r[3], "G");
    // A1 1/2 and A2 0/1 stay, B1 1/1 (A) -> 2/2, C1 0/2 -> 0/3, C2 1/2 -> 1/3
    assert_eq!(
        &r[9..14],
        [
            "1/2:0,7,8,.:90,50,40,30,0,20,.,.,.,.:15:0.7,0.8,.",
            "0/1:6,6,0,.:0,30,60,40,70,90,.,.,.,.:12:0.3,.,.",
            "2/2:0,.,9,.:99,.,.,40,.,0,.,.,.,.:9:.,1,.",
            "0/3:6,0,.,7:0,40,90,.,.,.,10,50,.,80:13:0.1,.,0.2",
            "1/3:0,4,.,5:70,20,10,.,.,.,60,30,.,0:9:.,.,0.6"
        ]
    );
    // AC and AN are recomputed from the merged genotypes.
    assert!(r[7].ends_with("AN=10;AC=3,3,2"), "{}", r[7]);
    // With -m none only records sharing an allele pair: a and c, b stays apart.
    let none = records(&merge3(&["-m", "none"]));
    assert!(
        none.iter()
            .any(|l| l.starts_with("1\t150\t.\tG\tA\t70\tPASS"))
    );
    assert!(
        none.iter()
            .any(|l| cols(l)[1] == "150" && cols(l)[4] == "T,A,C")
    );
}

#[test]
fn site_columns_follow_bcftools() {
    let recs = records(&merge3(&[]));
    let first = cols(&recs[0]);
    // 1:100: IDs joined without repeats, QUAL the largest, FILTER united, DP summed.
    let r = recs
        .iter()
        .map(|l| cols(l))
        .find(|c| c[1] == "100" && c[4] == "C,G")
        .unwrap();
    assert_eq!(r[2], "rs100;rs101");
    assert_eq!(r[5], "45");
    assert_eq!(r[6], "PASS");
    assert!(
        r[7].starts_with("NAME=first;DB;DP=25;DP4=5,6,6,6;"),
        "{}",
        r[7]
    );
    let _ = first;
    // 1:150: PASS is dropped next to q10; the first record's FILTER stays.
    let r = recs
        .iter()
        .map(|l| cols(l))
        .find(|c| c[1] == "150" && c[4] == "T,A,C")
        .unwrap();
    assert_eq!(r[6], "q10");
    assert_eq!(r[5], "70");
}

#[test]
fn missing_samples_take_the_ploidy_of_the_record() {
    let dir = tempdir();
    let h = |n: &str, rows: &[&str]| s(&plain_file(dir.path(), n, &vcf(n, rows)));
    let a = h(
        "A",
        &[
            "1\t10\t.\tA\tC\t.\t.\t.\tGT\t1",
            "1\t20\t.\tA\tC\t.\t.\t.\tGT\t1",
        ],
    );
    let b = h(
        "B",
        &[
            "1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1",
            "1\t30\t.\tA\tC\t.\t.\t.\tGT\t0/1",
        ],
    );
    let text = out(&["vcf", "merge", &a, &b]);
    let r: Vec<Vec<&str>> = records(&text)
        .iter()
        .map(|l| l.split('\t').map(|_| "").collect())
        .collect();
    let _ = r;
    let recs = records(&text);
    // haploid next to a diploid: the absent one is diploid ./.; alone it is haploid .
    assert!(recs[0].ends_with("GT\t1\t0/1"));
    assert!(recs[1].ends_with("GT\t1\t."), "{}", recs[1]);
    assert!(recs[2].ends_with("GT\t./.\t0/1"), "{}", recs[2]);
    // -0 fills 0/0 or 0 in the same ploidy
    let z = records(&out(&["vcf", "merge", &a, &b, "-0"]));
    assert!(z[1].ends_with("GT\t1\t0"), "{}", z[1]);
    assert!(z[2].ends_with("GT\t0/0\t0/1"), "{}", z[2]);
}

#[test]
fn modes_pair_or_split_records() {
    let at = |m: &str| {
        let t = merge3(&["-m", m]);
        records(&t)
            .iter()
            .filter(|l| l.starts_with("1\t100\t"))
            .count()
    };
    assert_eq!(at("none"), 2);
    assert_eq!(at("snps"), 1);
    assert_eq!(at("both"), 1);
    assert_eq!(at("all"), 1);
    assert_eq!(at("indels"), 2);
    // -m id pairs on equal IDs only: rs100 with rs100, rs101 apart.
    assert_eq!(at("id"), 2);
    // The SNP and the indel at 1:300 merge only with -m all (REF AT).
    let both = merge3(&["-m", "both"]);
    assert_eq!(
        records(&both)
            .iter()
            .filter(|l| l.starts_with("1\t300\t"))
            .count(),
        2
    );
    let all = merge3(&["-m", "all"]);
    assert!(
        records(&all)
            .iter()
            .any(|l| l.starts_with("1\t300\t") && cols(l)[3] == "AT")
    );
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let one = merge3(&["-t", "1"]);
    for t in ["2", "8"] {
        assert_eq!(merge3(&["-t", t]), one);
    }
    // and in BCF
    let dir = tempdir();
    let (a, b, c) = (
        s(&fixture("a.vcf")),
        s(&fixture("b.vcf")),
        s(&fixture("c.vcf")),
    );
    let mut bytes = Vec::new();
    for t in ["1", "2", "8"] {
        let o = s(&dir.path().join(format!("o{t}.bcf")));
        out(&["vcf", "merge", &a, &b, &c, "-o", &o, "-t", t]);
        bytes.push(std::fs::read(&o).unwrap());
    }
    assert_eq!(bytes[0], bytes[1]);
    assert_eq!(bytes[0], bytes[2]);
}

#[test]
fn bcf_inputs_and_outputs_give_the_same_records() {
    let dir = tempdir();
    let want = records(&merge3(&[]));
    let mut bcfs = Vec::new();
    for n in ["a", "b", "c"] {
        let o = s(&dir.path().join(format!("{n}.bcf")));
        out(&[
            "view",
            &s(&fixture(&format!("{n}.vcf"))),
            "--vcf-span",
            "ref",
            "-O",
            "bcf",
            "-o",
            &o,
        ]);
        bcfs.push(o);
    }
    let mut args = vec!["vcf", "merge"];
    args.extend(bcfs.iter().map(String::as_str));
    // BCF in keeps BCF out; text with -O vcf
    let mut text = args.clone();
    text.extend(["-O", "vcf"]);
    assert_eq!(records(&out(&text)), want);
    let ob = s(&dir.path().join("m.bcf"));
    let mut b = args.clone();
    b.extend(["-o", &ob]);
    out(&b);
    assert_eq!(
        records(&out(&["view", &ob, "--vcf-span", "ref", "-O", "vcf"])),
        want
    );
}

#[test]
fn duplicate_sample_names_name_both_files_and_force_renames() {
    let (a, d) = (s(&fixture("a.vcf")), s(&fixture("dup.vcf")));
    gloci()
        .args(["vcf", "merge", &a, &d])
        .assert()
        .code(1)
        .stderr(
            contains("A1")
                .and(contains("a.vcf"))
                .and(contains("dup.vcf"))
                .and(contains("--force-samples")),
        );
    let text = out(&["vcf", "merge", &a, &d, "--force-samples"]);
    let chrom = text.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(chrom.ends_with("FORMAT\tA1\tA2\t2:A1\tD9"), "{chrom}");
}

#[test]
fn out_of_order_input_is_an_error_with_the_fix() {
    let dir = tempdir();
    let ok = s(&plain_file(
        dir.path(),
        "ok.vcf",
        &vcf("OK", &["1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1"]),
    ));
    let bad = s(&plain_file(
        dir.path(),
        "bad.vcf",
        &vcf(
            "BAD",
            &[
                "1\t30\t.\tA\tC\t.\t.\t.\tGT\t0/1",
                "1\t20\t.\tA\tC\t.\t.\t.\tGT\t0/1",
            ],
        ),
    ));
    gloci()
        .args(["vcf", "merge", &ok, &bad])
        .assert()
        .code(1)
        .stderr(contains("bad.vcf").and(contains("gloci sort")));
}

#[test]
fn unsupported_bcftools_options_are_usage_errors_naming_the_option() {
    let (a, b) = (s(&fixture("a.vcf")), s(&fixture("b.vcf")));
    for (flag, val, name) in [
        ("--gvcf", Some("ref.fa"), "-g/--gvcf"),
        ("--local-alleles", Some("2"), "-L/--local-alleles"),
        ("--print-header", None, "--print-header"),
        ("--use-header", Some("h.txt"), "--use-header"),
        ("--regions", Some("1:1-10"), "-r/--regions"),
        ("-R", Some("r.bed"), "-R/--regions-file"),
        ("--missing-rules", Some("PL:max"), "-M/--missing-rules"),
        ("--regions-overlap", Some("1"), "--regions-overlap"),
    ] {
        let mut args = vec!["vcf", "merge", &a, &b, flag];
        if let Some(v) = val {
            args.push(v);
        }
        gloci()
            .args(&args)
            .assert()
            .code(2)
            .stderr(contains(name).and(contains("bcftools merge")));
    }
}

#[test]
fn input_errors_say_what_and_how_to_fix() {
    let a = s(&fixture("a.vcf"));
    gloci()
        .args(["vcf", "merge", &a])
        .assert()
        .code(2)
        .stderr(contains("at least two").and(contains("--file-list")));
    gloci()
        .args(["vcf", "merge", &a, "-"])
        .assert()
        .code(2)
        .stderr(contains("stdin"));
    gloci()
        .args(["vcf", "merge", &a, &s(&fixture("b.vcf")), "-i", "NOPE:sum"])
        .assert()
        .code(2)
        .stderr(contains("NOPE").and(contains("not defined")));
    gloci()
        .args(["vcf", "merge", &a, &s(&fixture("b.vcf")), "-i", "DP:foo"])
        .assert()
        .code(2)
        .stderr(contains("foo"));
    gloci()
        .args(["vcf", "merge", &a, &s(&fixture("b.vcf")), "-f", ""])
        .assert()
        .success();
}

#[test]
fn file_list_reads_the_inputs_from_a_file() {
    let dir = tempdir();
    let list = dir.path().join("list.txt");
    std::fs::write(
        &list,
        format!(
            "{}\n{}\n{}\n",
            s(&fixture("a.vcf")),
            s(&fixture("b.vcf")),
            s(&fixture("c.vcf"))
        ),
    )
    .unwrap();
    let l = s(&list);
    assert_eq!(
        records(&out(&["vcf", "merge", "--file-list", &l])),
        records(&merge3(&[]))
    );
}

#[test]
fn two_hundred_files_merge_and_a_low_file_limit_is_explained() {
    let dir = tempdir();
    let mut files = Vec::new();
    for i in 0..200 {
        let name = format!("s{i:03}.vcf");
        let row = format!("1\t{}\t.\tA\tC\t.\t.\t.\tGT\t0/1", 10 + i % 5);
        files.push(s(&plain_file(
            dir.path(),
            &name,
            &vcf(&format!("S{i:03}"), &[&row]),
        )));
    }
    let mut args = vec!["vcf", "merge"];
    args.extend(files.iter().map(String::as_str));
    let text = out(&args);
    let chrom = text.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert_eq!(chrom.split('\t').count(), 9 + 200);
    assert_eq!(records(&text).len(), 5);
    assert!(records(&text)[0].matches("0/1").count() == 40);
    // Under a limit of 100 handles the run stops with the way out.
    let bin = assert_cmd::cargo::cargo_bin("gloci");
    let o = std::process::Command::new("sh")
        .arg("-c")
        .arg(format!(
            "ulimit -n 100 && exec '{}' vcf merge {}",
            bin.display(),
            files.join(" ")
        ))
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(1));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("too many open files") && err.contains("ulimit -n"),
        "{err}"
    );
}

#[test]
fn only_the_records_of_one_position_are_held() {
    // 50 positions, 3 files: no more than the records of one position (3) are held.
    let dir = tempdir();
    let rows: Vec<String> = (0..50)
        .map(|i| format!("1\t{}\t.\tA\tC\t.\t.\t.\tGT\t0/1", 10 + i))
        .collect();
    let rr: Vec<&str> = rows.iter().map(String::as_str).collect();
    let mk = |n: &str| s(&plain_file(dir.path(), &format!("{n}.vcf"), &vcf(n, &rr)));
    let (a, b, c) = (mk("A"), mk("B"), mk("C"));
    let o = gloci()
        .args(["-v", "vcf", "merge", &a, &b, &c])
        .assert()
        .success();
    let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
    let held: u64 = err
        .split("at most ")
        .nth(1)
        .and_then(|t| t.split(' ').next())
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| panic!("no peak in: {err}"));
    assert_eq!(held, 3);
}

#[test]
fn info_rules_and_filters() {
    // -i - keeps DP of the first record
    let t = merge3(&["-i", "-"]);
    let r = records(&t)
        .into_iter()
        .find(|l| l.starts_with("1\t100\t") && cols(l)[4] == "C,G")
        .unwrap();
    assert!(cols(&r)[7].contains("DP=10;"), "{}", cols(&r)[7]);
    // user rules replace the defaults
    let t = merge3(&["-i", "DP:max"]);
    let r = records(&t)
        .into_iter()
        .find(|l| l.starts_with("1\t100\t") && cols(l)[4] == "C,G")
        .unwrap();
    assert!(cols(&r)[7].contains("DP=12;"), "{}", cols(&r)[7]);
    // -f PASS drops 1:150 a (q10)
    let t = merge3(&["-f", "PASS"]);
    assert!(records(&t).iter().all(|l| cols(l)[6] != "q10"));
    // -F x: any input without a filter makes the record PASS
    let t = merge3(&["-F", "x"]);
    assert!(records(&t).iter().any(|l| cols(l)[6] == "PASS"));
}
