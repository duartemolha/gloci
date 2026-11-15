//! Without `-g` or a header that lists contigs, the chromosome order is
//! lexicographic (D32, amended); `--chrom-order karyotype` selects the other.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

const LEX: &str = "chr1\t1\t5\nchr10\t1\t5\nchr2\t1\t5\nchrX\t1\t5\n";
const KARYO: &str = "chr1\t1\t5\nchr2\t1\t5\nchr10\t1\t5\nchrX\t1\t5\n";
const VCF_HEAD: &str = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn p(d: &std::path::Path, name: &str, text: &str) -> String {
    plain_file(d, name, text).to_string_lossy().into_owned()
}

fn vcf_of(chroms: &[&str]) -> String {
    let mut s = VCF_HEAD.to_string();
    for c in chroms {
        s.push_str(&format!("{c}\t10\t.\tA\tG\t.\t.\t.\n"));
    }
    s
}

#[test]
fn relation_commands_expect_lex_by_default() {
    let d = tempdir();
    let lex = p(d.path(), "lex.bed", LEX);
    let kar = p(d.path(), "kar.bed", KARYO);
    gloci()
        .args(["intersect", "-a", &lex, "-b", &lex])
        .assert()
        .success()
        .stdout(LEX);
    gloci()
        .args(["intersect", "-a", &kar, "-b", &kar])
        .assert()
        .code(1)
        .stderr(contains("looks sorted in karyotype order"))
        .stderr(contains("--chrom-order karyotype"));
    gloci()
        .args([
            "intersect",
            "-a",
            &kar,
            "-b",
            &kar,
            "--chrom-order",
            "karyotype",
        ])
        .assert()
        .success()
        .stdout(KARYO);
    // The reverse direction keeps its hint.
    gloci()
        .args([
            "intersect",
            "-a",
            &lex,
            "-b",
            &lex,
            "--chrom-order",
            "karyotype",
        ])
        .assert()
        .code(1)
        .stderr(contains(
            "looks sorted by sort -k1,1, sort-bed or bedtools sort",
        ));
}

#[test]
fn merge_and_cluster_expect_lex_by_default() {
    let d = tempdir();
    let lex = p(d.path(), "lex.bed", LEX);
    let kar = p(d.path(), "kar.bed", KARYO);
    for cmd in ["merge", "cluster"] {
        gloci().args([cmd, &lex]).assert().success();
        gloci()
            .args([cmd, &kar])
            .assert()
            .code(1)
            .stderr(contains("looks sorted in karyotype order"));
        gloci()
            .args([cmd, &kar, "--chrom-order", "karyotype"])
            .assert()
            .success();
    }
    gloci().args(["merge", &lex]).assert().success().stdout(LEX);
}

#[test]
fn vcf_without_contigs_is_lex_by_default() {
    let d = tempdir();
    let lex = p(d.path(), "lex.vcf", &vcf_of(&["chr1", "chr10", "chr2"]));
    let kar = p(d.path(), "kar.vcf", &vcf_of(&["chr1", "chr2", "chr10"]));
    gloci().args(["check", &lex]).assert().success();
    gloci()
        .args(["check", &kar])
        .assert()
        .code(1)
        .stdout(contains("looks sorted in karyotype order"));
    gloci()
        .args(["check", &kar, "--chrom-order", "karyotype"])
        .assert()
        .success();
    gloci()
        .args(["sort", &kar])
        .assert()
        .success()
        .stdout(contains("chr1\t10").and(contains("chr10\t10")));
    let out = gloci().args(["sort", &kar]).assert().success();
    let text = String::from_utf8(out.get_output().stdout.clone()).unwrap();
    let order: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(order, ["chr1", "chr10", "chr2"]);
}

/// `vcf filter --thin` needs each chromosome in one block, positions ascending, in any
/// chromosome order: a VCF in reference order without `##contig` lines is accepted.
#[test]
fn vcf_filter_thin_takes_any_chromosome_order() {
    let d = tempdir();
    let text = vcf_of(&["chr1", "chr2", "chr10"]);
    let kar = p(d.path(), "kar.vcf", &text);
    let out = gloci()
        .args(["vcf", "filter", &kar, "--thin", "5"])
        .assert()
        .success();
    let got = String::from_utf8(out.get_output().stdout.clone()).unwrap();
    let records: Vec<&str> = got.lines().filter(|l| !l.starts_with('#')).collect();
    let want: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(records, want);
    // A chromosome that comes back is still out of order.
    let back = p(d.path(), "back.vcf", &vcf_of(&["chr1", "chr2", "chr1"]));
    gloci()
        .args(["vcf", "filter", &back, "--thin", "5"])
        .assert()
        .code(1)
        .stderr(contains("back.vcf is not sorted"));
}

/// The note of an order error names only options the command has: `vcf concat`, `merge` and
/// `annotate` take `--chrom-order` but no `-g`; `filter`, `extract`, `stats` and `popgen` take
/// neither.
#[test]
fn vcf_order_notes_name_only_the_commands_options() {
    let d = tempdir();
    let kar = p(d.path(), "kar.vcf", &vcf_of(&["chr1", "chr2", "chr10"]));
    // Out of order in every chromosome order.
    let back = p(d.path(), "back.vcf", &vcf_of(&["chr1", "chr2", "chr1"]));
    let back_dp = p(
        d.path(),
        "back_dp.vcf",
        &vcf_of(&["chr1", "chr2", "chr1"]).replace(
            "#CHROM",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n#CHROM",
        ),
    );
    let mut with_sample = String::from(
        "##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n",
    );
    for c in ["chr1", "chr2", "chr1"] {
        with_sample.push_str(&format!("{c}\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"));
    }
    let back_gt = p(d.path(), "back_gt.vcf", &with_sample);
    let prefix = d.path().join("o").to_string_lossy().into_owned();
    let cases: Vec<Vec<&str>> = vec![
        vec!["concat", &kar],
        vec!["merge", &kar, &kar],
        vec!["annotate", &kar, "-a", &kar, "-c", "ID"],
        vec!["filter", &back, "--thin", "5"],
        vec!["extract", &back_dp, "--info", "DP", "--thin", "5"],
        vec![
            "stats",
            &back,
            "--snp-density",
            "1000",
            "--out-prefix",
            &prefix,
        ],
        vec![
            "popgen",
            &back_gt,
            "--window-pi",
            "1000",
            "--out-prefix",
            &prefix,
        ],
    ];
    for args in cases {
        let out = gloci().arg("vcf").args(&args).output().unwrap();
        let err = String::from_utf8_lossy(&out.stderr).into_owned();
        assert_eq!(out.status.code(), Some(1), "vcf {args:?}: {err}");
        assert!(err.contains("gloci sort"), "vcf {args:?}: {err}");
        let help = gloci()
            .args(["vcf", args[0], "--help"])
            .output()
            .unwrap()
            .stdout;
        let help = String::from_utf8(help).unwrap();
        for (named, option) in [
            (" -g ", "-g, --genome"),
            ("--chrom-order", "--chrom-order <"),
        ] {
            assert!(
                !err.contains(named) || help.contains(option),
                "vcf {}: the message names {named:?}, which the command does not have: {err}",
                args[0]
            );
        }
    }
}

/// The migration guide's recipe for BAM input to `merge`: the BED `convert` writes keeps the
/// BAM header's order (chr2 before chr10), which `-g x.bam` names; without it the lexicographic
/// default refuses the stream.
#[test]
fn a_bed_converted_from_a_bam_merges_with_g_the_bam() {
    let d = tempdir();
    let sam = p(
        d.path(),
        "x.sam",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n\
         @SQ\tSN:chr10\tLN:1000\n\
         r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
         r2\t0\tchr2\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
         r3\t0\tchr10\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n",
    );
    let bam = d.path().join("x.bam").to_string_lossy().into_owned();
    gloci().args(["view", &sam, "-o", &bam]).assert().success();
    let bed = gloci()
        .args(["convert", &bam])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    gloci()
        .args(["merge", "-"])
        .write_stdin(bed.clone())
        .assert()
        .code(1)
        .stderr(contains("looks sorted in karyotype order"));
    gloci()
        .args(["merge", "-", "-g", &bam])
        .write_stdin(bed)
        .assert()
        .success()
        .stdout("chr1\t9\t14\nchr2\t9\t14\nchr10\t9\t14\n");
}
