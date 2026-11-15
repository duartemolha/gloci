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
