//! Tests for `gloci vcf filter`.
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn stdout_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn records_of(s: &str) -> Vec<&str> {
    s.lines().filter(|l| !l.starts_with('#')).collect()
}

const SAMPLE_VCF_HEAD: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=FLAG,Number=0,Type=Flag,Description=\"A flag\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3
";

#[test]
fn end_to_end_categories() {
    let d = tempdir();
    let rows = [
        "chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30;FLAG\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/1:10:30:PASS\t1/1:10:30:PASS",
        "chr1\t200\trs2\tG\tT\t20\tPASS\tDP=15\tGT:DP:GQ:FT\t0/1:5:20:PASS\t0/1:5:20:PASS\t0/0:5:20:PASS",
        "chr1\t300\trs3\tC\tCAA\t40\tq10\tDP=25\tGT:DP:GQ:FT\t0/1:8:15:q10\t0/0:8:15:PASS\t./.:9:20:PASS",
        "chr2\t500\trs4\tAT\tA\t60\tPASS\tDP=60\tGT:DP:GQ:FT\t1|1:20:40:PASS\t1|1:20:40:PASS\t1|1:20:40:PASS",
    ];
    let full_text = format!("{SAMPLE_VCF_HEAD}{}\n", rows.join("\n"));
    let vcf = s(&plain_file(d.path(), "sample.vcf", &full_text));

    // Sample selection
    let out = stdout_of(&["vcf", "filter", &vcf, "--samples", "s1,s2"]);
    let head_chrom = out.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert_eq!(
        head_chrom,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"
    );

    // Genotype masking (--min-dp 8: rs2 samples all have DP 5, so masked to ./., missingness = 1.0)
    let out = stdout_of(&[
        "vcf",
        "filter",
        &vcf,
        "--min-dp",
        "8",
        "--max-missing",
        "0.5",
    ]);
    let recs = records_of(&out);
    assert!(!recs.iter().any(|r| r.contains("rs2")));
    assert!(recs.iter().any(|r| r.contains("rs1")));

    // Type (--types indel)
    let out = stdout_of(&["vcf", "filter", &vcf, "--types", "indel"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 2);
    assert!(recs[0].contains("rs3"));
    assert!(recs[1].contains("rs4"));

    // Biallelic
    let out = stdout_of(&["vcf", "filter", &vcf, "--biallelic"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 4);

    // Quality and FILTER (--min-qual 30, --pass-only)
    let out = stdout_of(&["vcf", "filter", &vcf, "--min-qual", "30", "--pass-only"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 2);
    assert!(recs[0].contains("rs1"));
    assert!(recs[1].contains("rs4"));

    // Keep filter
    let out = stdout_of(&["vcf", "filter", &vcf, "--keep-filter", "q10"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs3"));

    // INFO flags
    let out = stdout_of(&["vcf", "filter", &vcf, "--require-info", "FLAG"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs1"));

    // Frequency: rs1 has 3 called alleles of ALT out of 6 -> MAF = 0.5
    let out = stdout_of(&[
        "vcf",
        "filter",
        &vcf,
        "--min-maf",
        "0.45",
        "--max-maf",
        "0.55",
    ]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs1"));

    // Missingness: rs3 has 1 missing out of 3 samples (missing fraction 1/3 ~ 0.333)
    let out = stdout_of(&["vcf", "filter", &vcf, "--max-missing", "0.2"]);
    let recs = records_of(&out);
    assert!(!recs.iter().any(|r| r.contains("rs3")));

    // Depth: min-mean-dp
    let out = stdout_of(&["vcf", "filter", &vcf, "--min-mean-dp", "15.0"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs4"));

    // Phased only: only rs4 is all phased
    let out = stdout_of(&["vcf", "filter", &vcf, "--phased-only"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs4"));

    // Thinning: rs1 at 100, rs2 at 200, rs3 at 300 on chr1. With --thin 150: rs1 (100) kept, rs2 (200) dropped (diff 100 < 150), rs3 (300) kept (diff 200 >= 150).
    let out = stdout_of(&["vcf", "filter", &vcf, "--thin", "150"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 3);
    assert!(recs[0].contains("rs1"));
    assert!(recs[1].contains("rs3"));
    assert!(recs[2].contains("rs4"));

    // Chrom / Exclude chrom
    let out = stdout_of(&["vcf", "filter", &vcf, "--chrom", "chr2"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs4"));

    let out = stdout_of(&["vcf", "filter", &vcf, "--exclude-chrom", "chr1"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs4"));

    // IDs
    let id_file = plain_file(d.path(), "ids.txt", "rs1\nrs4\n");
    let out = stdout_of(&["vcf", "filter", &vcf, "--ids", &s(&id_file)]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 2);
    assert!(recs[0].contains("rs1"));
    assert!(recs[1].contains("rs4"));

    // Positions
    let pos_file = plain_file(d.path(), "pos.txt", "chr1\t100\nchr2\t500\n");
    let out = stdout_of(&["vcf", "filter", &vcf, "--positions", &s(&pos_file)]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 2);
    assert!(recs[0].contains("rs1"));
    assert!(recs[1].contains("rs4"));
}

#[test]
fn pipeline_order() {
    let d = tempdir();
    // A site whose ALT carriers all have DP < 10, run with `--min-dp 10 --min-mac 1`, is dropped.
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=23\tGT:DP:GQ:FT\t0/0:15:30:PASS\t0/1:5:20:PASS\t1/1:3:20:PASS\n"
    );
    let vcf = s(&plain_file(d.path(), "site.vcf", &text));
    let out = stdout_of(&["vcf", "filter", &vcf, "--min-dp", "10", "--min-mac", "1"]);
    assert!(records_of(&out).is_empty());
}

#[test]
fn subset_then_maf() {
    let d = tempdir();
    // 3 samples: s1 is 0/0, s2 is 0/0, s3 is 1/1.
    // Overall: 4 '0', 2 '1' -> MAF = 2/6 = 0.333333.
    // On s2,s3: 2 '0', 2 '1' -> MAF = 2/4 = 0.5.
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/0:10:30:PASS\t1/1:10:30:PASS\n"
    );
    let vcf = s(&plain_file(d.path(), "subset.vcf", &text));

    // Without subset: fails min-maf 0.4
    let out = stdout_of(&["vcf", "filter", &vcf, "--min-maf", "0.4"]);
    assert!(records_of(&out).is_empty());

    // With subset s2,s3: passes min-maf 0.4
    let out = stdout_of(&[
        "vcf",
        "filter",
        &vcf,
        "--samples",
        "s2,s3",
        "--min-maf",
        "0.4",
    ]);
    assert_eq!(records_of(&out).len(), 1);
}

#[test]
fn region_with_index_and_scan() {
    let d = tempdir();
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/1:10:30:PASS\t1/1:10:30:PASS\n\
chr1\t500\trs2\tG\tT\t20\tPASS\tDP=15\tGT:DP:GQ:FT\t0/1:5:20:PASS\t0/1:5:20:PASS\t0/0:5:20:PASS\n\
chr2\t200\trs3\tC\tT\t40\tPASS\tDP=25\tGT:DP:GQ:FT\t0/1:8:15:PASS\t0/1:8:15:PASS\t./.:9:20:PASS\n"
    );

    // 1. Scan plain VCF
    let vcf = s(&plain_file(d.path(), "scan.vcf", &text));
    let out = stdout_of(&["vcf", "filter", &vcf, "--region", "chr1:400-600"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs2"));

    // 2. Indexed BGZF VCF
    let gz = s(&bgzf_file(d.path(), "indexed.vcf.gz", &text));
    gloci().args(["index", &gz]).assert().success();
    let out = stdout_of(&["vcf", "filter", &gz, "--region", "chr1:400-600"]);
    let recs = records_of(&out);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].contains("rs2"));
}

fn raw_bcf_payload(path: &Path) -> Vec<u8> {
    let mut f = std::fs::File::open(path).unwrap();
    let mut gz = noodles_bgzf::io::Reader::new(&mut f);
    let mut decompressed = Vec::new();
    gz.read_to_end(&mut decompressed).unwrap();
    // BCF header is: 5 bytes magic + 4 bytes l_text + l_text bytes
    let l_text = u32::from_le_bytes(decompressed[5..9].try_into().unwrap()) as usize;
    decompressed[9 + l_text..].to_vec()
}

#[test]
fn bcf_passthrough() {
    let d = tempdir();
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/1:10:30:PASS\t1/1:10:30:PASS\n\
chr1\t200\trs2\tG\tT\t20\tPASS\tDP=15\tGT:DP:GQ:FT\t0/1:5:20:PASS\t0/1:5:20:PASS\t0/0:5:20:PASS\n"
    );
    let bcf_in = d.path().join("in.bcf");
    std::fs::write(&bcf_in, fixtures::vcf_to_bcf(text.as_bytes())).unwrap();

    let bcf_out = d.path().join("out.bcf");
    // Filter with passthrough plan: keeps rs1, drops rs2
    gloci()
        .args([
            "vcf",
            "filter",
            &s(&bcf_in),
            "--min-qual",
            "30",
            "-o",
            &s(&bcf_out),
        ])
        .assert()
        .success();

    let in_payload = raw_bcf_payload(&bcf_in);
    let out_payload = raw_bcf_payload(&bcf_out);

    // The first record in in_payload (rs1) must match out_payload byte for byte!
    assert!(!out_payload.is_empty());
    assert_eq!(in_payload[..out_payload.len()], out_payload[..]);
}

#[test]
fn stdin_input() {
    let _d = tempdir();
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/1:10:30:PASS\t1/1:10:30:PASS\n"
    );
    let o = gloci()
        .args(["vcf", "filter", "-", "--pass-only"])
        .write_stdin(text)
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(records_of(&out).len(), 1);
}

#[test]
fn determinism() {
    let d = tempdir();
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/1:10:30:PASS\t1/1:10:30:PASS\n\
chr1\t200\trs2\tG\tT\t20\tPASS\tDP=15\tGT:DP:GQ:FT\t0/1:5:20:PASS\t0/1:5:20:PASS\t0/0:5:20:PASS\n\
chr2\t500\trs4\tAT\tA\t60\tPASS\tDP=60\tGT:DP:GQ:FT\t1|1:20:40:PASS\t1|1:20:40:PASS\t1|1:20:40:PASS\n"
    );
    let vcf = s(&plain_file(d.path(), "det.vcf", &text));

    let o1 = stdout_of(&["vcf", "filter", &vcf, "-t", "1"]);
    let o2 = stdout_of(&["vcf", "filter", &vcf, "-t", "2"]);
    let o4 = stdout_of(&["vcf", "filter", &vcf, "-t", "4"]);

    assert_eq!(o1, o2);
    assert_eq!(o1, o4);
}

#[test]
fn usage_errors() {
    let d = tempdir();
    let text = format!(
        "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP:GQ:FT\t0/0:10:30:PASS\t0/1:10:30:PASS\t1/1:10:30:PASS\n"
    );
    let vcf = s(&plain_file(d.path(), "err.vcf", &text));

    // --biallelic with --min-alleles
    gloci()
        .args(["vcf", "filter", &vcf, "--biallelic", "--min-alleles", "2"])
        .assert()
        .code(2);

    // --biallelic with --max-alleles
    gloci()
        .args(["vcf", "filter", &vcf, "--biallelic", "--max-alleles", "2"])
        .assert()
        .code(2);

    // MAF outside [0, 1]
    gloci()
        .args(["vcf", "filter", &vcf, "--min-maf", "1.5"])
        .assert()
        .code(2);

    gloci()
        .args(["vcf", "filter", &vcf, "--max-maf", "-0.1"])
        .assert()
        .code(2);

    // Missing fraction outside [0, 1]
    gloci()
        .args(["vcf", "filter", &vcf, "--max-missing", "1.1"])
        .assert()
        .code(2);

    // Samples conflict
    let samples_f = plain_file(d.path(), "samples.txt", "s1\n");
    gloci()
        .args([
            "vcf",
            "filter",
            &vcf,
            "--samples",
            "s1",
            "--samples-file",
            &s(&samples_f),
        ])
        .assert()
        .code(2);

    // Unknown sample
    gloci()
        .args(["vcf", "filter", &vcf, "--samples", "unknown_sample"])
        .assert()
        .code(2)
        .stderr(contains("unknown sample"));
}

#[test]
fn library_api_regions_and_regions_file() {
    use gloci_core::ChromTable;
    use gloci_io::{Format, ReaderOptions, RecordWriter, open_reader};
    use gloci_ops::{VcfFilterOptions, vcf_filter};

    struct TestWriter(Vec<u8>);
    impl RecordWriter for TestWriter {
        fn write_record(&mut self, record: &[u8]) -> Result<(), std::io::Error> {
            self.0.extend_from_slice(record);
            self.0.push(b'\n');
            Ok(())
        }
        fn finish(self: Box<Self>) -> Result<(), gloci_io::IoError> {
            Ok(())
        }
    }

    let d = tempfile::tempdir().unwrap();
    let vcf = plain_file(
        d.path(),
        "in.vcf",
        format!(
            "{SAMPLE_VCF_HEAD}\
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30\tGT:DP\t0/1:10\t0/1:10\t0/1:10\n\
chr1\t200\trs2\tG\tT\t50\tPASS\tDP=30\tGT:DP\t0/1:10\t0/1:10\t0/1:10\n\
chr2\t300\trs3\tT\tA\t50\tPASS\tDP=30\tGT:DP\t0/1:10\t0/1:10\t0/1:10\n"
        )
        .as_str(),
    );
    let bed = plain_file(d.path(), "regions.bed", "chr1\t150\t250\n");

    let ro = ReaderOptions {
        format: Some(Format::Vcf),
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: gloci_io::InvalidPolicy::Fail,
        threads: std::num::NonZero::new(1).unwrap(),
    };

    // Test with o.regions
    {
        let reader = open_reader(&s(&vcf), &ro).unwrap();
        let mut chroms = ChromTable::first_appearance();
        let opts = VcfFilterOptions {
            regions: vec!["chr1:150-250".to_string()],
            ..Default::default()
        };
        let mut w = TestWriter(Vec::new());
        let sum = vcf_filter(reader, &mut chroms, &opts, &mut w).unwrap();
        // Regions are selected before the filters, which see one record.
        assert_eq!(sum.records_in, 1);
        assert_eq!(sum.records_out, 1);
        let out = String::from_utf8(w.0).unwrap();
        assert!(out.contains("rs2"));
        assert!(!out.contains("rs1"));
    }

    // Test with o.regions_file
    {
        let reader = open_reader(&s(&vcf), &ro).unwrap();
        let mut chroms = ChromTable::first_appearance();
        let opts = VcfFilterOptions {
            regions_file: Some(s(&bed)),
            ..Default::default()
        };
        let mut w = TestWriter(Vec::new());
        let sum = vcf_filter(reader, &mut chroms, &opts, &mut w).unwrap();
        // Regions are selected before the filters, which see one record.
        assert_eq!(sum.records_in, 1);
        assert_eq!(sum.records_out, 1);
        let out = String::from_utf8(w.0).unwrap();
        assert!(out.contains("rs2"));
        assert!(!out.contains("rs1"));
    }
}

/// I9: BCF output of a missing haploid GT next to diploid ones is a clear
/// error naming the file and record, with VCF output suggested.
#[test]
fn bcf_mixed_ploidy_missing_is_a_clear_error() {
    let d = tempdir();
    let input = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t.\n";
    let p = plain_file(d.path(), "mixed.vcf", input);
    let p = p.to_string_lossy().into_owned();
    for args in [
        vec!["view", p.as_str(), "-o"],
        vec!["vcf", "filter", p.as_str(), "-o"],
    ] {
        let out = d.path().join("out.bcf");
        let o = gloci()
            .args(&args)
            .arg(&out)
            .assert()
            .code(1)
            .get_output()
            .stderr
            .clone();
        let e = String::from_utf8(o).unwrap();
        assert!(e.contains(&format!("{p} ")), "{e}");
        assert!(
            e.contains("chr1:10: BCF output of a record that mixes a missing haploid"),
            "{e}"
        );
        assert!(e.contains("-o out.vcf.gz"), "{e}");
        assert!(!out.exists(), "partial output left");
    }
    // VCF output of the same record works.
    gloci()
        .args(["vcf", "filter", &p, "-o"])
        .arg(d.path().join("out.vcf.gz"))
        .assert()
        .success();
}
