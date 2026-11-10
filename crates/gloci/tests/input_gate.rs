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
        .args(["merge", sam.to_str().unwrap()])
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
fn convert_reads_alignments() {
    // `convert` opts in to alignments: a header-only SAM has no reads and writes nothing.
    let dir = tempdir();
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n");
    gloci()
        .args(["convert", sam.to_str().unwrap()])
        .assert()
        .success()
        .stdout("");
}

#[test]
fn genomecov_reads_alignments() {
    // `genomecov` reads alignments: a SAM with no reads covers its header's chromosomes with zero.
    let dir = tempdir();
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:5\n");
    gloci()
        .args(["genomecov", sam.to_str().unwrap(), "--bga"])
        .assert()
        .success()
        .stdout("chr1\t0\t5\t0\n");
}

#[test]
fn planned_commands_say_not_yet_and_others_say_never() {
    let dir = tempdir();
    let bam = fixtures::write(
        dir.path(),
        "x.bam",
        &fixtures::bgzip(b"BAM\x01\0\0\0\0\0\0\0\0"),
    );
    let bed = plain_file(dir.path(), "a.bed", "chr1\t1\t2\n");
    let x = bam.to_str().unwrap();
    let b = bed.to_str().unwrap();
    let not_yet = |args: &[&str]| {
        let o = gloci()
            .args(args)
            .assert()
            .code(2)
            .stderr(contains("not available in this build yet"));
        let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
        assert!(!err.contains("does not read alignments"), "{args:?}: {err}");
    };
    not_yet(&["sort", x]);
    gloci()
        .args(["subtract", "-a", b, "-b", x])
        .assert()
        .code(2)
        .stderr(contains("does not read alignments"));
}
