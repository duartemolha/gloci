//! `convert --skip-invalid` says how many lines it dropped, whichever format it reads (D171
//! "no silent wrong answers").
mod common;
use common::*;
use predicates::str::contains;

fn skipped(name: &str, text: &str, extra: &[&str]) -> String {
    let d = tempdir();
    let f = plain_file(d.path(), name, text);
    let o = gloci()
        .arg("convert")
        .arg("--skip-invalid")
        .arg(&f)
        .args(extra)
        .assert()
        .success()
        .stderr(contains("skipped 1 invalid"));
    String::from_utf8_lossy(&o.get_output().stdout).into_owned()
}

#[test]
fn skip_invalid_warns_for_gff_vcf_and_wig() {
    let out = skipped(
        "a.gff3",
        "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\tx\t300\t.\t+\t.\tID=g2\n",
        &[],
    );
    assert_eq!(out.lines().count(), 1, "{out}");
    let out = skipped(
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\nchr1\tx\t.\tA\tC\t.\t.\t.\n",
        &[],
    );
    assert_eq!(out.lines().count(), 1, "{out}");
    let out = skipped("a.wig", "variableStep chrom=chr1\n5 1.5\nx 2\n9 3\n", &[]);
    assert_eq!(out.lines().count(), 2, "{out}");
}

#[test]
fn skip_invalid_warns_for_bed_to_bam() {
    let d = tempdir();
    let bed = plain_file(
        d.path(),
        "a.bed",
        "chr1\t0\t10\tr1\nchr1\tx\t30\tr2\nchr1\t40\t50\tr3\n",
    );
    let g = plain_file(d.path(), "g.genome", "chr1\t1000\n");
    let bam = d.path().join("o.bam");
    gloci()
        .args(["convert", "--skip-invalid", "--bam", "-g"])
        .arg(&g)
        .arg(&bed)
        .arg("-o")
        .arg(&bam)
        .assert()
        .success()
        .stderr(contains("skipped 1 invalid"));
    let v = gloci()
        .arg("view")
        .arg(&bam)
        .args(["-O", "sam", "--no-header"])
        .assert()
        .success();
    let text = String::from_utf8_lossy(&v.get_output().stdout).into_owned();
    let names: Vec<&str> = text
        .lines()
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(names, ["r1", "r3"]);
}

#[test]
fn skip_invalid_warns_for_bedpe_to_bam() {
    let d = tempdir();
    let pe = plain_file(
        d.path(),
        "a.bedpe",
        "chr1\t0\t10\tchr1\t50\t60\tp1\nchr1\tx\t10\tchr1\t50\t60\tp2\nchr1\t100\t110\tchr1\t150\t160\tp3\n",
    );
    let g = plain_file(d.path(), "g.genome", "chr1\t1000\n");
    let bam = d.path().join("o.bam");
    gloci()
        .args([
            "convert",
            "--skip-invalid",
            "--bam",
            "--format",
            "bedpe",
            "-g",
        ])
        .arg(&g)
        .arg(&pe)
        .arg("-o")
        .arg(&bam)
        .assert()
        .success()
        .stderr(contains("skipped 1 invalid"));
    let v = gloci()
        .arg("view")
        .arg(&bam)
        .args(["-O", "sam", "--no-header"])
        .assert()
        .success();
    assert_eq!(
        String::from_utf8_lossy(&v.get_output().stdout)
            .lines()
            .count(),
        4
    );
}
