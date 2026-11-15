//! Fixtures of the samtools-style command matrices.

use crate::common::*;
use crate::fixtures_alignments::*;
use crate::helpers::*;
use std::fmt::Write as _;
use std::path::Path;

/// One command of the P3c matrices: its name, its arguments, whether it takes `-t`, and whether
/// its output is the same for BAM, SAM and CRAM copies of the same reads.
pub(crate) struct P3cCmd {
    pub(crate) name: String,
    pub(crate) args: Vec<String>,
    pub(crate) threads: bool,
    pub(crate) format_free: bool,
}

pub(crate) fn p3c_run(c: &P3cCmd, t: &str) -> std::process::Output {
    if c.threads {
        p3a_run_t(&c.args, t)
    } else {
        gloci().args(&c.args).output().unwrap()
    }
}

impl P3a {
    /// The P3c commands over one alignment input.
    pub(crate) fn p3c_cmds(&self, kind: &str, input: &str) -> Vec<P3cCmd> {
        let s = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
        let with_ref = |mut a: Vec<String>| -> Vec<String> {
            if kind == "cram" {
                a.extend(s(&["--reference", &self.fa]));
            }
            a
        };
        let bed = self.bed.as_str();
        let mut c: Vec<P3cCmd> = Vec::new();
        let mut add = |name: &str, args: Vec<String>, threads: bool, format_free: bool| {
            c.push(P3cCmd {
                name: format!("{kind} {name}"),
                threads,
                // a column name is the path
                format_free: format_free && !args.iter().any(|a| a == "--column-names"),
                args,
            })
        };
        add("idxstats", with_ref(s(&["idxstats", input])), true, true);
        for fmt in ["text", "tsv", "json"] {
            add(
                &format!("flagstat {fmt}"),
                with_ref(s(&["flagstat", input, "-O", fmt])),
                true,
                true,
            );
        }
        add(
            "view --count",
            with_ref(s(&["view", "--count", input])),
            true,
            true,
        );
        add(
            "view --count filters",
            with_ref(s(&[
                "view",
                "--count",
                input,
                "--min-mapq",
                "30",
                "--exclude-flags",
                "0x400",
            ])),
            true,
            true,
        );
        for (name, extra) in [
            ("depth", vec![]),
            ("depth -a", vec!["--all-positions"]),
            ("depth -aa", vec!["--all-references", "--column-names"]),
            ("depth -s", vec!["--skip-mate-overlap"]),
            (
                "depth filters",
                vec!["--min-mapq", "10", "--min-baseq", "20", "--count-deletions"],
            ),
            (
                "depth regions",
                vec!["--region", "chr1:100-400", "--region", "chr2:1-50"],
            ),
            ("depth regions-file", vec!["--regions-file", bed]),
        ] {
            let mut a = s(&["depth", input]);
            a.extend(extra.iter().map(|x| x.to_string()));
            add(name, with_ref(a), true, true);
        }
        if kind != "sam" {
            add("bedcov", with_ref(s(&["bedcov", bed, input])), true, true);
            add(
                "bedcov --read-count",
                with_ref(s(&[
                    "bedcov",
                    bed,
                    input,
                    "--read-count",
                    "--bases-at-depth",
                    "2",
                    "--split",
                ])),
                true,
                true,
            );
            add(
                "bedcov --min-mapq",
                with_ref(s(&["bedcov", bed, input, "--min-mapq", "20"])),
                true,
                true,
            );
        }
        add(
            "head",
            with_ref(s(&["head", input, "--records", "4"])),
            true,
            false,
        );
        add(
            "head --headers",
            with_ref(s(&["head", input, "--headers", "2"])),
            true,
            false,
        );
        add(
            "check --quick",
            s(&["check", "--quick", input]),
            true,
            false,
        );
        add("samples", s(&["samples", input]), false, false);
        add(
            "samples --show-index",
            s(&["samples", input, "--show-index", "--column-names"]),
            false,
            false,
        );
        c
    }

    /// The P3c commands over the FASTA files, which have no format variants.
    pub(crate) fn p3c_fasta_cmds(&self) -> Vec<P3cCmd> {
        let s = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
        let mk = |name: &str, args: Vec<String>, threads: bool| P3cCmd {
            name: name.to_string(),
            args,
            threads,
            format_free: false,
        };
        vec![
            mk(
                "faidx",
                s(&[
                    "faidx",
                    &self.fa,
                    "chr1:1-200",
                    "chr2",
                    "--reverse-complement",
                ]),
                true,
            ),
            mk(
                "faidx --regions-file",
                s(&["faidx", &self.fa, "--regions-file", &self.masks]),
                true,
            ),
            mk(
                "faidx bgzf",
                s(&["faidx", &self.fa_gz, "chr1:50-300", "chrX:10-"]),
                true,
            ),
            mk("dict", s(&["dict", &self.fa, "--alias"]), true),
            mk("dict bgzf", s(&["dict", &self.fa_gz]), true),
            mk(
                "flags",
                s(&["flags", "0x704", "1796", "secondary,dup"]),
                false,
            ),
            mk("flags table", s(&["flags"]), false),
        ]
    }
}

/// SAM text of an alignment file with the `@PG` lines (which name `-t` and the paths) dropped.
pub(crate) fn p3c_sam_of(path: &str, extra: &[&str]) -> Vec<u8> {
    let o = gloci()
        .args(["head", path, "--records", "1000000"])
        .args(extra)
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0), "head {path}");
    p3a_norm(&o.stdout)
}

/// A BED of `n` regions of 3 kb over the two chromosomes of `p3a_big_alignments`.
pub(crate) fn p3c_big_bed(dir: &Path, n: usize) -> String {
    let mut bed = String::new();
    for i in 0..n {
        let start = (i / 2) * 4_900;
        writeln!(bed, "chr{}\t{start}\t{}\tr{i}", 1 + i % 2, start + 3_000).unwrap();
    }
    s_of(&plain_file(dir, "p3c-big.bed", &bed))
}

/// A header-only SAM with `n` references, so that `idxstats` and `depth --all-references` write
/// far more than a pipe holds.
pub(crate) fn p3c_many_references(dir: &Path, n: usize) -> String {
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for i in 0..n {
        writeln!(sam, "@SQ\tSN:contig{i:06}\tLN:1000").unwrap();
    }
    s_of(&plain_file(dir, "many-refs.sam", &sam))
}

/// A FASTA of `n` sequences of 120 bp, indexed.
pub(crate) fn p3c_many_sequences(dir: &Path, n: usize) -> String {
    let mut fa = String::new();
    for i in 0..n {
        writeln!(fa, ">seq{i:05} sequence number {i}").unwrap();
        let line: String = (0..60).map(|k| b"ACGT"[(i + k * k) & 3] as char).collect();
        writeln!(fa, "{line}\n{line}").unwrap();
    }
    let f = plain_file(dir, "many-seqs.fa", &fa);
    gloci().args(["index", &s_of(&f)]).assert().success();
    s_of(&f)
}
