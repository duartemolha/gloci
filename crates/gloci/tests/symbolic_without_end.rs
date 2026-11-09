//! A symbolic ALT with neither INFO/END nor INFO/SVLEN: an error under the
//! default allele span (honouring `--skip-invalid`), but never under
//! `--vcf-span ref`, where the footprint is the REF span (as bcftools and tabix).
mod common;
use common::*;
use predicates::str::contains;

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10\ta\tA\tC\t.\t.\t.
chr1\t20\tdel\tN\t<DEL>\t.\t.\tSVTYPE=DEL
chr1\t30\tins\tN\t<INS>\t.\t.\t.
chr1\t40\tgvcf\tN\t<*>\t.\t.\t.
chr1\t50\tnonref\tN\t<NON_REF>\t.\t.\t.
";

fn ids(out: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(out)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(2).unwrap().to_string())
        .collect()
}

fn input() -> (tempfile::TempDir, String) {
    let d = tempdir();
    let f = plain_file(d.path(), "in.vcf", VCF)
        .to_string_lossy()
        .into_owned();
    (d, f)
}

const FIX: &str = "add INFO/END or INFO/SVLEN, or use --vcf-span ref, or --skip-invalid";

#[test]
fn default_span_errors_with_the_fix() {
    let (_d, f) = input();
    for cmd in [["view", &f[..]], ["vcf", "filter"]] {
        let mut c = gloci();
        c.args(cmd);
        if cmd[0] == "vcf" {
            c.arg(&f);
        }
        c.assert()
            .code(1)
            .stderr(contains(
                "line 5: symbolic ALT <DEL> needs INFO/END or INFO/SVLEN",
            ))
            .stderr(contains(FIX));
    }
}

#[test]
fn skip_invalid_drops_only_the_unspanned_record() {
    let (_d, f) = input();
    for args in [vec!["view", &f[..]], vec!["vcf", "filter", &f[..]]] {
        let o = gloci()
            .args(&args)
            .arg("--skip-invalid")
            .assert()
            .success()
            .stderr(contains("skipped 1 invalid records"));
        assert_eq!(
            ids(&o.get_output().stdout),
            ["a", "ins", "gvcf", "nonref"],
            "{args:?}"
        );
    }
}

#[test]
fn ref_span_keeps_every_record() {
    let (_d, f) = input();
    for args in [vec!["view", &f[..]], vec!["vcf", "filter", &f[..]]] {
        let o = gloci()
            .args(&args)
            .args(["--vcf-span", "ref"])
            .assert()
            .success();
        assert_eq!(
            ids(&o.get_output().stdout),
            ["a", "del", "ins", "gvcf", "nonref"],
            "{args:?}"
        );
    }
}

#[test]
fn sort_orders_by_ref_end_and_keeps_every_record() {
    let (_d, f) = input();
    let o = gloci().args(["sort", &f]).assert().success();
    assert_eq!(
        ids(&o.get_output().stdout),
        ["a", "del", "ins", "gvcf", "nonref"]
    );
}

#[test]
fn other_symbolic_alts_never_need_end() {
    let d = tempdir();
    let text = VCF.replace("chr1\t20\tdel\tN\t<DEL>\t.\t.\tSVTYPE=DEL\n", "");
    let f = plain_file(d.path(), "in2.vcf", &text)
        .to_string_lossy()
        .into_owned();
    let o = gloci().args(["vcf", "filter", &f]).assert().success();
    assert_eq!(ids(&o.get_output().stdout), ["a", "ins", "gvcf", "nonref"]);
}
