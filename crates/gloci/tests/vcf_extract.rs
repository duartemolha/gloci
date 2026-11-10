//! Tests for `gloci vcf extract` (decision D73): INFO and FORMAT tables and
//! the 012 genotype matrix.
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn stdout_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

const HEAD: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, as a list\">
##INFO=<ID=FLAG,Number=0,Type=Flag,Description=\"A flag\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3
";

const ROWS: &[&str] = &[
    "chr1\t100\trs1\tA\tC\t50\tPASS\tDP=40;FLAG\tGT:DP:AD:FT\t0/0:10:10,0:PASS\t0/1:.:5,5:PASS\t1/1:12:0,12:lowq",
    "chr1\t200\t.\tG\tT,C\t20\tPASS\tAF=0.5,0.2\tGT:DP\t1/2:4\t2/2:4\t0|1",
    "chr1\t300\trs3\tC\tCAA\t40\tq10\t.\tGT\t0/1\t0/0\t./.",
    "chr2\t500\trs4\tAT\tA\t60\tPASS\tDP=5\tGT:DP:AD\t1|1:20:0,20\t0/1:21\t./.:.",
];

fn vcf_text() -> String {
    format!("{HEAD}{}\n", ROWS.join("\n"))
}

fn inputs(d: &Path) -> (String, String) {
    let vcf = s(&plain_file(d, "in.vcf", &vcf_text()));
    let bcf = d.join("in.bcf");
    std::fs::write(&bcf, fixtures::vcf_to_bcf(vcf_text().as_bytes())).unwrap();
    (vcf, s(&bcf))
}

#[test]
fn info_table_missing_values_and_flags() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    assert_eq!(
        stdout_of(&["vcf", "extract", &vcf, "--info", "DP,AF,FLAG"]),
        "CHROM\tPOS\tID\tREF\tALT\tDP\tAF\tFLAG\n\
chr1\t100\trs1\tA\tC\t40\t.\t1\n\
chr1\t200\t.\tG\tT,C\t.\t0.5,0.2\t0\n\
chr1\t300\trs3\tC\tCAA\t.\t.\t0\n\
chr2\t500\trs4\tAT\tA\t5\t.\t0\n"
    );
}

#[test]
fn info_table_warns_about_an_undeclared_key() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "NOPE"])
        .assert()
        .success()
        .stderr(contains("does not declare INFO/NOPE"));
}

#[test]
fn format_wide_any_key_as_text() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    assert_eq!(
        stdout_of(&["vcf", "extract", &vcf, "--format", "GT,AD,FT"]),
        "CHROM\tPOS\tID\tREF\tALT\ts1:GT\ts1:AD\ts1:FT\ts2:GT\ts2:AD\ts2:FT\ts3:GT\ts3:AD\ts3:FT\n\
chr1\t100\trs1\tA\tC\t0/0\t10,0\tPASS\t0/1\t5,5\tPASS\t1/1\t0,12\tlowq\n\
chr1\t200\t.\tG\tT,C\t1/2\t.\t.\t2/2\t.\t.\t0|1\t.\t.\n\
chr1\t300\trs3\tC\tCAA\t0/1\t.\t.\t0/0\t.\t.\t./.\t.\t.\n\
chr2\t500\trs4\tAT\tA\t1|1\t0,20\t.\t0/1\t.\t.\t./.\t.\t.\n"
    );
}

#[test]
fn format_long_is_a_row_per_site_and_sample() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    let out = stdout_of(&[
        "vcf", "extract", &vcf, "--format", "GT,DP", "--layout", "long",
    ]);
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines[0], "CHROM\tPOS\tID\tREF\tALT\tSAMPLE\tGT\tDP");
    assert_eq!(lines.len(), 1 + 4 * 3);
    assert_eq!(lines[1], "chr1\t100\trs1\tA\tC\ts1\t0/0\t10");
    assert_eq!(lines[2], "chr1\t100\trs1\tA\tC\ts2\t0/1\t.");
    assert_eq!(lines[6], "chr1\t200\t.\tG\tT,C\ts3\t0|1\t.");
    assert_eq!(lines[12], "chr2\t500\trs4\tAT\tA\ts3\t./.\t.");
}

#[test]
fn selection_flags_apply() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    // --samples keeps the selected samples (in header order) only.
    let out = stdout_of(&[
        "vcf",
        "extract",
        &vcf,
        "--format",
        "DP",
        "--samples",
        "s3,s1",
    ]);
    assert_eq!(
        out.lines().next().unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\ts1:DP\ts3:DP"
    );
    // Site filters and the expression drop sites.
    let out = stdout_of(&["vcf", "extract", &vcf, "--info", "DP", "--pass-only"]);
    assert_eq!(out.lines().count(), 4);
    let out = stdout_of(&[
        "vcf",
        "extract",
        &vcf,
        "--info",
        "DP",
        "--include",
        "INFO/DP>10",
    ]);
    assert_eq!(out.lines().count(), 2);
    let out = stdout_of(&["vcf", "extract", &vcf, "--info", "DP", "--chrom", "chr2"]);
    assert_eq!(out.lines().nth(1).unwrap(), "chr2\t500\trs4\tAT\tA\t5");
    let out = stdout_of(&[
        "vcf",
        "extract",
        &vcf,
        "--info",
        "DP",
        "--region",
        "chr1:150-250",
    ]);
    assert_eq!(out.lines().count(), 2);
    // A masked genotype reads as missing (a sample without DP is not masked).
    let out = stdout_of(&["vcf", "extract", &vcf, "--format", "GT", "--min-dp", "11"]);
    assert_eq!(
        out.lines().nth(1).unwrap(),
        "chr1\t100\trs1\tA\tC\t./.\t0/1\t1/1"
    );
    let out = stdout_of(&["vcf", "extract", &vcf, "--info", "DP", "--thin", "250"]);
    assert_eq!(out.lines().count(), 3);
}

#[test]
fn bcf_and_vcf_give_identical_output() {
    let d = tempdir();
    let (vcf, bcf) = inputs(d.path());
    for extra in [
        vec!["--info", "DP,AF,FLAG"],
        vec!["--format", "GT,DP,AD,FT"],
        vec![
            "--format",
            "GT,AD",
            "--layout",
            "long",
            "--samples",
            "s1,s3",
        ],
        vec!["--info", "DP", "--include", "INFO/DP>3 && FMT/DP>5"],
        vec!["--format", "GT", "--min-dp", "11"],
    ] {
        let mut a = vec!["vcf", "extract", &vcf];
        a.extend(&extra);
        let mut b = vec!["vcf", "extract", &bcf];
        b.extend(&extra);
        assert_eq!(stdout_of(&a), stdout_of(&b), "{extra:?}");
    }
    for (name, input) in [("v", &vcf), ("b", &bcf)] {
        let prefix = s(&d.path().join(name));
        gloci()
            .args([
                "vcf",
                "extract",
                input,
                "--genotype-matrix",
                "012",
                "--out-prefix",
                &prefix,
            ])
            .assert()
            .success();
    }
    for ext in ["012", "012.indv", "012.pos"] {
        assert_eq!(
            std::fs::read(d.path().join(format!("v.{ext}"))).unwrap(),
            std::fs::read(d.path().join(format!("b.{ext}"))).unwrap(),
            "{ext}"
        );
    }
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let d = tempdir();
    let mut text = HEAD.to_string();
    for i in 0..3000u32 {
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tC\t50\tPASS\tDP={}\tGT:DP\t0/1:{}\t1/1:{}\t0/0:{}\n",
            100 + i,
            i % 97,
            i % 13,
            i % 7,
            i % 5
        ));
    }
    let vcf = s(&plain_file(d.path(), "big.vcf", &text));
    let bgz = d.path().join("big.vcf.gz");
    gloci()
        .args(["view", &vcf, "-o", &s(&bgz)])
        .assert()
        .success();
    let bgz = s(&bgz);
    let one = stdout_of(&[
        "vcf",
        "extract",
        &bgz,
        "--format",
        "GT,DP",
        "--threads",
        "1",
    ]);
    let four = stdout_of(&[
        "vcf",
        "extract",
        &bgz,
        "--format",
        "GT,DP",
        "--threads",
        "4",
    ]);
    assert_eq!(one, four);
    assert_eq!(one.lines().count(), 3001);
    for t in ["1", "4"] {
        let p = s(&d.path().join(format!("m{t}")));
        gloci()
            .args([
                "vcf",
                "extract",
                &bgz,
                "--genotype-matrix",
                "012",
                "--out-prefix",
                &p,
                "--threads",
                t,
            ])
            .assert()
            .success();
    }
    for ext in ["012", "012.indv", "012.pos"] {
        assert_eq!(
            std::fs::read(d.path().join(format!("m1.{ext}"))).unwrap(),
            std::fs::read(d.path().join(format!("m4.{ext}"))).unwrap()
        );
    }
}

fn read(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

fn fixture_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vcf_extract_012")
}

#[test]
fn genotype_matrix_matches_vcftools_files() {
    // expected.* were written by `vcftools --vcf input.vcf --012`.
    let d = tempdir();
    let prefix = s(&d.path().join("out"));
    gloci()
        .args([
            "vcf",
            "extract",
            &s(&fixture_dir().join("input.vcf")),
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
        ])
        .assert()
        .success()
        .stderr(contains("1 multi-allelic site"));
    for ext in ["012", "012.indv", "012.pos"] {
        assert_eq!(
            read(&d.path().join(format!("out.{ext}"))),
            read(&fixture_dir().join(format!("expected.{ext}"))),
            "{ext}"
        );
    }
}

#[test]
fn genotype_matrix_haploid_and_half_missing_differ_from_vcftools() {
    // vcftools-haploid.012 is `vcftools --012` on haploid.vcf: it writes 1 for
    // a haploid `0` and for `0/.` (docs/differences.md#vcf-extract-012).
    let d = tempdir();
    let prefix = s(&d.path().join("h"));
    gloci()
        .args([
            "vcf",
            "extract",
            &s(&fixture_dir().join("haploid.vcf")),
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
        ])
        .assert()
        .success();
    let mine = read(&d.path().join("h.012"));
    assert_eq!(mine, "0\t0\t0\t-1\n1\t1\t0\t1\n2\t1\t0\t-1\n3\t2\t1\t0\n");
    let theirs = read(&fixture_dir().join("vcftools-haploid.012"));
    assert_eq!(theirs, "0\t1\t1\t-1\n1\t1\t1\t1\n2\t1\t1\t-1\n3\t2\t1\t0\n");
    assert_ne!(mine, theirs);
    assert_eq!(read(&d.path().join("h.012.indv")), "m1\nm2\nf1\nf2\n");
    assert_eq!(
        read(&d.path().join("h.012.pos")),
        "chrX\t100\nchrX\t200\nchrX\t300\n"
    );
}

#[test]
fn genotype_matrix_follows_the_selection() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    let prefix = s(&d.path().join("sel"));
    gloci()
        .args([
            "vcf",
            "extract",
            &vcf,
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
            "--samples",
            "s1,s3",
            "--chrom",
            "chr1",
        ])
        .assert()
        .success();
    assert_eq!(read(&d.path().join("sel.012.indv")), "s1\ns3\n");
    // chr1:200 is multi-allelic and left out.
    assert_eq!(
        read(&d.path().join("sel.012.pos")),
        "chr1\t100\nchr1\t300\n"
    );
    assert_eq!(read(&d.path().join("sel.012")), "0\t0\t1\n1\t2\t-1\n");
}

#[test]
fn genotype_matrix_errors_name_the_record_and_leave_no_files() {
    let d = tempdir();
    let text = format!("{HEAD}chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/1/1\t0/0\n");
    let vcf = s(&plain_file(d.path(), "tri.vcf", &text));
    let prefix = s(&d.path().join("o"));
    gloci()
        .args([
            "vcf",
            "extract",
            &vcf,
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
        ])
        .assert()
        .code(1)
        .stderr(contains("sample s2 at chr1:100"))
        .stderr(contains("ploidy above 2"));
    assert!(!d.path().join("o.012").exists());
    assert!(!d.path().join("o.012.pos").exists());

    let text = format!("{HEAD}chr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\t2/2\t0/0\n");
    let vcf = s(&plain_file(d.path(), "bad.vcf", &text));
    gloci()
        .args([
            "vcf",
            "extract",
            &vcf,
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
        ])
        .assert()
        .code(1)
        .stderr(contains("names allele 2"));
}

#[test]
fn usage_errors() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    gloci().args(["vcf", "extract", &vcf]).assert().code(2);
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "DP", "--format", "GT"])
        .assert()
        .code(2);
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "DP", "--layout", "long"])
        .assert()
        .code(2)
        .stderr(contains("--layout applies to --format only"));
    gloci()
        .args(["vcf", "extract", &vcf, "--genotype-matrix", "012"])
        .assert()
        .code(2)
        .stderr(contains("--out-prefix"));
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "DP", "--out-prefix", "x"])
        .assert()
        .code(2)
        .stderr(contains("--genotype-matrix"));
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "DP,DP"])
        .assert()
        .code(2)
        .stderr(contains("names DP twice"));
    gloci()
        .args([
            "vcf",
            "extract",
            &vcf,
            "--info",
            "DP",
            "-o",
            &s(&d.path().join("x.vcf")),
        ])
        .assert()
        .code(2)
        .stderr(contains("writes a TSV"));
    gloci()
        .args([
            "vcf",
            "extract",
            &vcf,
            "--info",
            "DP",
            "--samples",
            "nobody",
        ])
        .assert()
        .code(2)
        .stderr(contains("unknown sample"));
    let bed = s(&plain_file(d.path(), "x.bed", "chr1\t1\t2\n"));
    gloci()
        .args(["vcf", "extract", &bed, "--info", "DP"])
        .assert()
        .code(2)
        .stderr(contains("works on VCF and BCF input"));
}

#[test]
fn table_goes_to_a_file_atomically() {
    let d = tempdir();
    let (vcf, _) = inputs(d.path());
    let out = d.path().join("info.tsv");
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "DP", "-o", &s(&out)])
        .assert()
        .success();
    assert_eq!(read(&out).lines().count(), 5);
    let gz = d.path().join("info.tsv.gz");
    gloci()
        .args(["vcf", "extract", &vcf, "--info", "DP", "-o", &s(&gz)])
        .assert()
        .success();
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut noodles_bgzf::io::Reader::new(std::fs::File::open(&gz).unwrap()),
        &mut text,
    )
    .unwrap();
    assert_eq!(text, read(&out));
}
