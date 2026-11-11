//! `gloci head` (D134): the rules beyond the samtools goldens (`head-*-vs-samtools`).
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Read;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn rg(ext: &str) -> String {
    s(&inputs().join("samtools").join(format!("rg.{ext}")))
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[&str]) -> (i32, String, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code().unwrap(),
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
    )
}

const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n\
@RG\tID:g1\tSM:NA12878\tLB:libA\n@RG\tID:g2\tSM:NA12891\n@RG\tID:g3\tSM:NA12878\tLB:libB\n\
@RG\tID:g4\tSM:HG002\n@RG\tID:g5\tLB:libC\n";
const R1: &str = "r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g1\n";
const R2: &str = "r2\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g2\n";
const R3: &str = "r3\t0\tchr2\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g4\n";

#[test]
fn head_defaults_and_limits() {
    for ext in ["bam", "sam"] {
        let f = rg(ext);
        // The default: the whole header, no record, and no @PG line.
        let (code, out, _) = run(&["head", &f]);
        assert_eq!((code, out.as_str()), (0, HEADER), "{ext}");
        // Records only.
        let (code, out, _) = run(&["head", &f, "--headers", "0", "--records", "3"]);
        assert_eq!((code, out), (0, format!("{R1}{R2}{R3}")), "{ext}");
        // A cut header, then records.
        let (_, out, _) = run(&["head", &f, "--headers", "2", "--records", "1"]);
        assert_eq!(
            out,
            format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n{R1}"),
            "{ext}"
        );
        // More records than the file has: all of them, exit 0.
        let (code, out, _) = run(&["head", &f, "--headers", "0", "--records", "100"]);
        assert_eq!((code, out), (0, format!("{R1}{R2}{R3}")), "{ext}");
        // More header lines than there are.
        let (_, out, _) = run(&["head", &f, "--headers", "100"]);
        assert_eq!(out, HEADER, "{ext}");
        assert!(!out.contains("@PG"));
    }
}

#[test]
fn head_is_the_same_at_any_thread_count() {
    let f = rg("bam");
    let one = run(&["head", &f, "--records", "2", "-t", "1"]);
    for t in ["2", "8"] {
        assert_eq!(run(&["head", &f, "--records", "2", "-t", t]), one);
    }
}

#[test]
fn head_notext_bam_writes_the_sq_lines() {
    let f = s(&inputs().join("aln/notext.bam"));
    let (code, out, _) = run(&["head", &f, "--records", "1"]);
    assert_eq!(code, 0);
    assert_eq!(
        out,
        "@SQ\tSN:c1\tLN:100\nr1\t0\tc1\t11\t60\t4M\t*\t0\t0\tACGT\t????\n"
    );
}

#[test]
fn head_truncated_bam_prints_what_it_read_and_exits_1() {
    // The decompressed BAM, cut inside its third record, under a good end-of-file block.
    let mut raw = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(rg("bam")).unwrap())
        .read_to_end(&mut raw)
        .unwrap();
    raw.truncate(raw.len() - 10);
    let d = tempdir();
    let f = s(&fixtures::write(
        d.path(),
        "cut.bam",
        &fixtures::bgzip(&raw),
    ));
    let o = gloci()
        .args(["head", &f, "--records", "3"])
        .assert()
        .code(1)
        .get_output()
        .clone();
    assert_eq!(
        String::from_utf8(o.stdout).unwrap(),
        format!("{HEADER}{R1}{R2}")
    );
    assert!(!o.stderr.is_empty());
}

#[test]
fn head_vcf_and_bed() {
    let d = tempdir();
    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t5\t.\tA\tT\t.\t.\t.\nchr1\t6\t.\tA\tT\t.\t.\t.\nchr1\t7\t.\tA\tT\t.\t.\t.\n";
    let gz = s(&bgzf_file(d.path(), "a.vcf.gz", vcf));
    let (code, out, _) = run(&["head", &gz, "--records", "2"]);
    let want: String = vcf.lines().take(5).map(|l| format!("{l}\n")).collect();
    assert_eq!((code, out), (0, want.clone()));
    assert!(!want.contains("glociCommand"));
    // BCF prints VCF text, with no ##glociCommand line either.
    let bcf = s(&fixtures::write(
        d.path(),
        "a.bcf",
        &fixtures::vcf_to_bcf(vcf.as_bytes()),
    ));
    let (code, out, _) = run(&["head", &bcf, "--records", "2"]);
    assert_eq!(code, 0);
    assert!(!out.contains("glociCommand"), "{out}");
    assert!(out.starts_with("##fileformat=VCFv4.2\n"), "{out}");
    assert!(
        out.ends_with("chr1\t5\t.\tA\tT\t.\t.\t.\nchr1\t6\t.\tA\tT\t.\t.\t.\n"),
        "{out}"
    );
    assert_eq!(out.lines().count(), 5, "{out}");
    // No records by default.
    let (_, out, _) = run(&["head", &gz]);
    assert_eq!(out.lines().count(), 3);
    // BED: `#` and track lines, then the data lines as they are.
    let bed = "#note\ntrack name=x\nchr1\t10\t20\tfoo\t5\t+\nchr1\t30\t40\n";
    let b = s(&plain_file(d.path(), "a.bed", bed));
    let (_, out, _) = run(&["head", &b, "--records", "1"]);
    assert_eq!(out, "#note\ntrack name=x\nchr1\t10\t20\tfoo\t5\t+\n");
    let (_, out, _) = run(&["head", &b, "--headers", "1", "--records", "5"]);
    assert_eq!(out, format!("#note\n{}", &bed[bed.find("chr1").unwrap()..]));
    let (_, out, _) = run(&["head", &b, "--headers", "0"]);
    assert_eq!(out, "");
}

#[test]
fn head_refuses_what_it_does_not_take() {
    // head describes the whole file: no --chrom-map and no read filters.
    let f = rg("bam");
    gloci()
        .args(["head", &f, "--chrom-map", "x.tsv"])
        .assert()
        .code(2);
    gloci()
        .args(["head", &f, "--min-mapq", "5"])
        .assert()
        .code(2);
    gloci()
        .args(["head", "/nonexistent/x.bam"])
        .assert()
        .code(1)
        .stderr(contains("/nonexistent/x.bam"));
    // A FASTA is not an input of head.
    let d = tempdir();
    let fa = s(&plain_file(d.path(), "a.fa", ">a\nACGT\n"));
    gloci()
        .args(["head", &fa])
        .assert()
        .failure()
        .stdout(predicate::str::is_empty());
}
