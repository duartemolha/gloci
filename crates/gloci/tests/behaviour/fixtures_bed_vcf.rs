//! Fixtures of the BED and VCF command matrices.

use crate::common::*;
use crate::helpers::*;
use gloci_testkit::fixtures::{vcf_to_bcf, write};
use std::fmt::Write as _;
use std::path::Path;

pub(crate) const P2_VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n\
##contig=<ID=chr2,length=10000>\n##contig=<ID=chr3,length=10000>\n##contig=<ID=chr4,length=10000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n##FILTER=<ID=q10,Description=\"low\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n";

/// Sorted multi-chromosome VCF with three samples and a mix of SNPs and indels.
pub(crate) fn p2_vcf(chroms: std::ops::RangeInclusive<usize>, per_chrom: usize) -> String {
    let mut t = String::from(P2_VCF_HEAD);
    for c in chroms {
        for i in 0..per_chrom {
            let (r, a) = match i % 5 {
                3 => ("AT", "A"),
                4 => ("A", "ATT"),
                k => (["A", "C", "G"][k % 3], ["C", "G", "T"][k % 3]),
            };
            let gt = ["0/1", "1/1", "0/0", "./.", "1|0"];
            writeln!(
                t,
                "chr{c}\t{}\trs{c}_{i}\t{r}\t{a}\t{}\t{}\tDP={}\tGT:DP\t{}:{}\t{}:{}\t{}:{}",
                100 + i * 37,
                20 + i % 60,
                if i % 7 == 0 { "q10" } else { "PASS" },
                i % 40,
                gt[i % 5],
                i % 11,
                gt[(i + c) % 5],
                (i + 3) % 17,
                gt[(i + 2 * c) % 5],
                (i + 5) % 23
            )
            .unwrap();
        }
    }
    t
}

/// Every P2 fixture, as plain / BGZF / indexed BGZF where the command takes it.
pub(crate) struct P2 {
    pub(crate) dir: tempfile::TempDir,
}

impl P2 {
    pub(crate) fn new() -> P2 {
        P2 { dir: tempdir() }
    }
    pub(crate) fn plain(&self, n: &str, t: &str) -> String {
        s_of(&plain_file(self.dir.path(), n, t))
    }
    pub(crate) fn bgzf(&self, n: &str, t: &str) -> String {
        s_of(&bgzf_file(self.dir.path(), n, t))
    }
    pub(crate) fn indexed(&self, n: &str, t: &str) -> String {
        let p = self.bgzf(n, t);
        gloci().args(["index", &p]).assert().success();
        p
    }
    pub(crate) fn bcf(&self, n: &str, t: &str) -> String {
        s_of(&write(self.dir.path(), n, &vcf_to_bcf(t.as_bytes())))
    }
}

pub(crate) fn bed_a(chroms: usize) -> String {
    let mut t = String::new();
    for c in 1..=chroms {
        for i in 0..60 {
            let st = 100 + i * 40 + (i % 3) * 7;
            writeln!(
                t,
                "chr{c}\t{st}\t{}\ta_{c}_{i}\t{}\t{}",
                st + 20 + i % 45,
                10 + i % 90,
                if i % 2 == 0 { '+' } else { '-' }
            )
            .unwrap();
        }
    }
    t
}

pub(crate) fn bed_b(chroms: usize) -> String {
    let mut t = String::new();
    for c in 1..=chroms {
        for i in 0..40 {
            let st = 130 + i * 55;
            writeln!(
                t,
                "chr{c}\t{st}\t{}\tb_{c}_{i}\t{}\t+",
                st + 30 + i % 25,
                i % 50
            )
            .unwrap();
        }
    }
    t
}

/// The BED-shaped commands of the determinism matrices; `{a}` and `{b}` are
/// the inputs.
pub(crate) fn p2_bed_cmds(g: &str) -> Vec<(&str, Vec<&str>)> {
    vec![
        ("shift", vec!["shift", "{a}", "--by", "7", "-g", g]),
        (
            "shift -s --pct",
            vec!["shift", "{a}", "--by", "-0.25", "--pct", "-s"],
        ),
        ("cluster", vec!["cluster", "{a}"]),
        ("cluster -d -s", vec!["cluster", "{a}", "-d", "30", "-s"]),
        ("spacing", vec!["spacing", "{a}"]),
        (
            "groupby",
            vec!["groupby", "{a}", "-g", "1", "-c", "5:sum,5:max,4:collapse"],
        ),
        (
            "partition",
            vec!["partition", "{a}", "{b}", "--names", "A,B"],
        ),
        (
            "partition --all",
            vec![
                "partition",
                "{a}",
                "{b}",
                "--names",
                "A,B",
                "--all",
                "--merge-adjacent",
            ],
        ),
        ("genomecov", vec!["genomecov", "{a}", "-g", g]),
        (
            "genomecov --bedgraph",
            vec!["genomecov", "{a}", "-g", g, "--bedgraph"],
        ),
        (
            "genomecov --bga",
            vec!["genomecov", "{a}", "-g", g, "--bga"],
        ),
        (
            "genomecov --per-base",
            vec!["genomecov", "{a}", "-g", g, "--per-base"],
        ),
        (
            "genomecov --split",
            vec!["genomecov", "{a}", "-g", g, "--bedgraph", "--split"],
        ),
    ]
}

/// (Writes a plain `.vcf`.) A one-record VCF that follows the 20 000-record one at the same
/// chromosome end (for `vcf concat`).
pub(crate) fn tail_after_big_vcf(dir: &Path) -> String {
    let mut t = String::from(P2_VCF_HEAD);
    t.push_str("chr1\t900000\t.\tA\tC\t50\tPASS\tDP=30\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n");
    s_of(&plain_file(dir, "tail.vcf", &t))
}
