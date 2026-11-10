//! The alignment/FASTA/PSL input gate (D100) and the unwritable CRAM output.
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::str::contains;

#[test]
fn input_gate() {
    let dir = tempdir();
    let bam = fixtures::write(
        dir.path(),
        "x.bam",
        &fixtures::bgzip(b"BAM\x01\0\0\0\0\0\0\0\0"),
    );
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n");
    let fa = plain_file(dir.path(), "ref.fa", ">chr1 d\nACGT\n");
    let psl = plain_file(
        dir.path(),
        "x.psl",
        "50\t1\t0\t0\t0\t0\t0\t0\t+\tq1\t51\t0\t51\tchr1\t1000\t100\t151\t1\t51,\t0,\t100,\n",
    );
    let bed = plain_file(dir.path(), "a.bed", "chr1\t1\t2\n");

    let x = bam.to_str().unwrap();
    gloci()
        .args(["merge", x])
        .assert()
        .code(2)
        .stderr(contains("does not read alignments"))
        .stderr(contains(format!("gloci convert {x} -o x.bed")))
        .stderr(contains("gloci sort -m"));
    gloci()
        .args(["sort", sam.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("does not read alignments"));
    gloci()
        .args(["sort", fa.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("getfasta, maskfasta or nuc"));
    let p = psl.to_str().unwrap();
    gloci()
        .args(["merge", p])
        .assert()
        .code(2)
        .stderr(contains(format!("gloci convert {p}")));
    gloci()
        .args(["sort", bed.to_str().unwrap(), "-o"])
        .arg(dir.path().join("o.cram"))
        .assert()
        .code(2)
        .stderr(contains("samtools view -C"));
}

#[test]
fn convert_accepts_alignments_past_the_gate() {
    // Until the readers land, `convert` is past the gate and says so (exit 1, not 2).
    let dir = tempdir();
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n");
    gloci()
        .args(["convert", sam.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(contains("not available in this build yet"));
}
