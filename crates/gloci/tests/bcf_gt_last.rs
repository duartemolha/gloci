//! GT need not be the first FORMAT field: a BCF whose FORMAT is `DP:PL:GT`
//! reads like the same VCF in `view`, `vcf filter` and `vcf merge`.
mod common;
use common::*;

const HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
##FILTER=<ID=q10,Description=\"q\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">\n";

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

#[test]
fn gt_last_in_a_bcf_reads_like_the_vcf() {
    let dir = tempdir();
    let a = vcf(
        "A",
        &[
            "1\t10\t.\tA\tC\t10\tPASS\t.\tDP:PL:GT\t5:1,2,3:0/1",
            "1\t20\t.\tA\tG\t10\tPASS\t.\tDP:GT\t7:1|1",
            "1\t30\t.\tA\tG\t10\tPASS\t.\tDP:PL:GT\t.:.:./.",
        ],
    );
    let b = vcf(
        "B",
        &[
            "1\t10\t.\tA\tG\t20\tPASS\t.\tPL:GT:DP\t4,5,6:1/1:9",
            "1\t20\t.\tA\tG\t10\tPASS\t.\tGT:DP\t0/1:3",
        ],
    );
    let (av, bv) = (
        plain_file(dir.path(), "a.vcf", &a),
        plain_file(dir.path(), "b.vcf", &b),
    );
    let (ab, bb) = (dir.path().join("a.bcf"), dir.path().join("b.bcf"));
    let s = |p: &std::path::Path| p.to_string_lossy().into_owned();
    let (av, bv, ab, bb) = (s(&av), s(&bv), s(&ab), s(&bb));
    out(&["view", &av, "-O", "bcf", "-o", &ab]);
    out(&["view", &bv, "-O", "bcf", "-o", &bb]);
    // view
    assert_eq!(
        records(&out(&["view", &ab, "-O", "vcf"])),
        records(&out(&["view", &av, "-O", "vcf"]))
    );
    assert_eq!(records(&out(&["view", &ab, "-O", "vcf"])), records(&a));
    // vcf filter
    assert_eq!(
        records(&out(&["vcf", "filter", &ab, "--min-dp", "1", "-O", "vcf"])),
        records(&out(&["vcf", "filter", &av, "--min-dp", "1", "-O", "vcf"]))
    );
    // vcf merge
    assert_eq!(
        records(&out(&["vcf", "merge", &ab, &bb, "-O", "vcf"])),
        records(&out(&["vcf", "merge", &av, &bv, "-O", "vcf"]))
    );
}
