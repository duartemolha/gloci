//! `gloci pair-intersect` at the command line: what its inputs must be, its
//! warnings and exit codes, and its output formats. The semantics are pinned
//! in gloci-ops' tests and by the bedtools goldens.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

const A: &str = "chr1\t10\t20\tchr1\t100\t110\tp1\t0\t+\t-\n";
const B: &str = "chr1\t12\t14\tb1\t0\t+\n";

fn files(dir: &std::path::Path) -> (String, String) {
    let a = plain_file(dir, "a.bedpe", A);
    let b = plain_file(dir, "b.bed", B);
    (
        a.to_str().unwrap().to_string(),
        b.to_str().unwrap().to_string(),
    )
}

#[test]
fn pair_intersect_default_is_either() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    gloci()
        .args(["pair-intersect", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout(format!("{}\t{}", A.trim_end(), B));
}

#[test]
fn pair_unknown_strand() {
    let dir = tempdir();
    // `.` strands on both ends: nothing matches -s, and gloci warns once about A.
    let a = plain_file(
        dir.path(),
        "a.bedpe",
        "chr1\t10\t20\tchr1\t100\t110\tp\t0\t.\t.\nchr1\t10\t20\tchr1\t100\t110\tq\t0\t.\t.\n",
    );
    let b = plain_file(dir.path(), "b.bed", B);
    let out = gloci()
        .args(["pair-intersect", "-s", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .unwrap();
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(err.matches("unknown strand").count(), 1, "{err}");
    assert!(err.contains("never match under -s or -S"), "{err}");
    // Without a strand rule there is nothing to warn about.
    gloci()
        .args(["pair-intersect", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .assert()
        .success()
        .stderr("");
}

#[test]
fn pair_span_types_refuse_strands() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    for flag in ["-s", "-S"] {
        gloci()
            .args([
                "pair-intersect",
                flag,
                "--type",
                "ispan",
                "-a",
                &a,
                "-b",
                &b,
            ])
            .assert()
            .code(2)
            .stderr(contains("ispan"));
    }
}

#[test]
fn pair_span_skips_are_one_warning() {
    let dir = tempdir();
    let a = plain_file(
        dir.path(),
        "a.bedpe",
        "chr1\t1\t2\tchr2\t1\t2\tp\nchr1\t1\t2\t.\t-1\t-1\tq\n",
    );
    let b = plain_file(dir.path(), "b.bed", B);
    let out = gloci()
        .args(["pair-intersect", "--type", "ospan", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .unwrap();
    assert!(out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(err.matches("left out of --type ospan").count(), 1, "{err}");
    assert!(err.contains("2 pairs"), "{err}");
}

#[test]
fn pair_a_must_be_bedpe_and_names_what_it_is() {
    let dir = tempdir();
    let b = plain_file(dir.path(), "b.bed", B);
    let vcf = plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    gloci()
        .args(["pair-intersect", "-a"])
        .arg(&vcf)
        .arg("-b")
        .arg(&b)
        .assert()
        .code(2)
        .stderr(contains("is VCF"))
        .stderr(contains("BEDPE"));
    // A SAM or BAM A points at `convert --bedpe`.
    let sam = plain_file(
        dir.path(),
        "reads.sam",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    );
    let bam = gloci_testkit::fixtures::write(
        dir.path(),
        "reads.bam",
        &gloci_testkit::fixtures::sam_to_bam(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
        ),
    );
    for a in [sam, bam] {
        let name = a.to_str().unwrap();
        gloci()
            .args(["pair-intersect", "-a", name, "-b"])
            .arg(&b)
            .assert()
            .code(2)
            .stderr(contains(format!("gloci convert {name} --bedpe")));
    }
}

#[test]
fn pair_a_is_bedpe_even_when_it_looks_like_bed() {
    let dir = tempdir();
    let a = plain_file(dir.path(), "pairs.txt", "chr1\t10\t20\tchr1\t100\t110\n");
    let b = plain_file(dir.path(), "b.bed", B);
    gloci()
        .args(["pair-intersect", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .assert()
        .success()
        .stdout("chr1\t10\t20\tchr1\t100\t110\tchr1\t12\t14\tb1\t0\t+\n");
}

#[test]
fn pair_bad_bedpe_lines_are_errors_unless_skipped() {
    let dir = tempdir();
    let b = plain_file(dir.path(), "b.bed", B);
    // Nine columns, a named chromosome with -1 coordinates, start > end, strand x.
    for bad in [
        "chr1\t10\t20\tchr1\t100\t110\tn\t0\t+\n",
        "chr1\t10\t20\tchr1\t-1\t-1\n",
        "chr1\t20\t10\tchr1\t100\t110\n",
        "chr1\t10\t20\tchr1\t100\t110\tn\t0\tx\ty\n",
    ] {
        let a = plain_file(dir.path(), "bad.bedpe", bad);
        gloci()
            .args(["pair-intersect", "-a"])
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .assert()
            .code(1)
            .stderr(contains("bad.bedpe"))
            .stderr(contains("line 1"));
        gloci()
            .args(["pair-intersect", "--skip-invalid", "-a"])
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .assert()
            .success()
            .stdout("");
    }
}

#[test]
fn pair_output_formats() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    let want = format!("{}\t{}", A.trim_end(), B);
    // TSV and BEDPE names work; BGZF for .gz; a BED name does not.
    for name in ["o.tsv", "o.bedpe"] {
        let o = dir.path().join(name);
        gloci()
            .args(["pair-intersect", "-a", &a, "-b", &b, "-o"])
            .arg(&o)
            .assert()
            .success();
        assert_eq!(std::fs::read_to_string(&o).unwrap(), want);
    }
    let gz = dir.path().join("o.bedpe.gz");
    gloci()
        .args(["pair-intersect", "-a", &a, "-b", &b, "-o"])
        .arg(&gz)
        .assert()
        .success();
    // BGZF: a gzip member with the BC extra field.
    let bytes = std::fs::read(&gz).unwrap();
    assert_eq!(&bytes[..2], &[0x1f, 0x8b]);
    assert_eq!(&bytes[12..14], b"BC");
    gloci()
        .args(["pair-intersect", "-a", &a, "-b", &b, "-o"])
        .arg(dir.path().join("o.vcf"))
        .assert()
        .code(2);
    gloci()
        .args(["pair-intersect", "-a", &a, "-b", &b, "--write-index", "-o"])
        .arg(&gz)
        .assert()
        .code(2)
        .stderr(contains("gloci index -s 1 -b 2 -e 3 --zero-based"));
}

#[test]
fn pair_split_is_refused() {
    let dir = tempdir();
    let (a, b) = files(dir.path());
    gloci()
        .args(["pair-intersect", "--split", "-a", &a, "-b", &b])
        .assert()
        .code(2)
        .stderr(contains("--split"));
    gloci()
        .args(["pair-intersect", "--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("--split").not());
}

#[test]
fn pair_threads_agree_on_the_command_line() {
    let dir = tempdir();
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..30_000u64 {
        a.push_str(&format!(
            "chr1\t{}\t{}\tchr2\t{}\t{}\tp{i}\n",
            i * 3,
            i * 3 + 20,
            i * 5,
            i * 5 + 20
        ));
        if i % 7 == 0 {
            b.push_str(&format!(
                "chr{}\t{}\t{}\tb{i}\n",
                i % 2 + 1,
                i * 4,
                i * 4 + 30
            ));
        }
    }
    let a = plain_file(dir.path(), "a.bedpe", &a);
    let b = plain_file(dir.path(), "b.bed", &b);
    for t in ["either", "notboth", "xor"] {
        let run = |threads: &str| {
            gloci()
                .args(["pair-intersect", "--type", t, "-t", threads, "-a"])
                .arg(&a)
                .arg("-b")
                .arg(&b)
                .output()
                .unwrap()
                .stdout
        };
        let one = run("1");
        assert!(!one.is_empty());
        assert!(one == run("4"), "--type {t}");
    }
}
