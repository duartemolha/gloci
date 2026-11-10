//! BEDPE is never detected, is read only by the paired commands, and is
//! written like TSV (D146).
mod common;
use common::*;
use predicates::str::contains;

const PAIR: &str = "chr1\t10\t20\tchr2\t30\t40\n";

#[test]
fn bedpe_shaped_file_detects_as_bed() {
    let dir = tempdir();
    let p = plain_file(dir.path(), "x.bedpe", PAIR);
    gloci()
        .args(["merge", p.to_str().unwrap()])
        .assert()
        .success()
        .stdout("chr1\t10\t20\n");
}

#[test]
fn bedpe_format_is_refused_by_other_commands() {
    let dir = tempdir();
    let p = plain_file(dir.path(), "x.bedpe", PAIR);
    gloci()
        .args(["merge", "--format", "bedpe", p.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("pair-intersect"))
        .stderr(contains("pair-pair"))
        .stderr(contains("shuffle --format bedpe"));
    // A VCF forced to BEDPE is refused the same way, with the format named
    // by the paired commands' reader (gloci-io's bedpe tests).
    let vcf = plain_file(
        dir.path(),
        "x.vcf",
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    gloci()
        .args(["merge", "--format", "bedpe", vcf.to_str().unwrap()])
        .assert()
        .code(2);
}

#[test]
fn bedpe_output_is_text_and_cannot_be_indexed_by_format() {
    let dir = tempdir();
    let a = plain_file(dir.path(), "a.bed", "chr1\t1\t5\nchr1\t7\t9\n");
    let b = plain_file(dir.path(), "b.bed", "chr1\t2\t3\n");
    let (a, b) = (a.to_str().unwrap(), b.to_str().unwrap());
    let want = "chr1\t1\t5\t1\nchr1\t7\t9\t0\n";
    let args = ["intersect", "-a", a, "-b", b, "--mode", "count"];
    let out = dir.path().join("o.bedpe");
    gloci()
        .args(args)
        .args(["-o", out.to_str().unwrap()])
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&out).unwrap(), want);
    let gz = dir.path().join("o.bedpe.gz");
    gloci()
        .args(args)
        .args(["-o", gz.to_str().unwrap()])
        .assert()
        .success();
    let raw = std::fs::read(&gz).unwrap();
    assert_eq!(
        gloci_io::detect_compression(&raw),
        gloci_io::Compression::Bgzf
    );
    gloci()
        .args(args)
        .args(["-o", gz.to_str().unwrap(), "--write-index"])
        .assert()
        .code(2)
        .stderr(contains("gloci index -s 1 -b 2 -e 3 --zero-based"));
    gloci()
        .args(args)
        .args(["-O", "bedpe"])
        .assert()
        .success()
        .stdout(want);
}
