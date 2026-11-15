//! Matrices of the BED and VCF commands.

use crate::common::*;
use crate::fixtures_bed_vcf::*;
use crate::helpers::*;
use gloci_testkit::fixtures::{bgzip, vcf_to_bcf, write};
use predicates::str::contains;
use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use std::process::Stdio;

#[test]
fn p2_determinism_matrix_streams() {
    let fx = P2::new();
    let genome = fx.plain(
        "g.genome",
        "chr1\t10000\nchr2\t10000\nchr3\t10000\nchr4\t10000\n",
    );
    let (a, b) = (bed_a(4), bed_b(4));
    let a_plain = fx.plain("a.bed", &a);
    let a_gz = fx.bgzf("a.bed.gz", &a);
    let a_idx = fx.indexed("a_idx.bed.gz", &a);
    let b_plain = fx.plain("b.bed", &b);
    let b_gz = fx.bgzf("b.bed.gz", &b);
    // unsorted, for --sort
    let mut rev: Vec<&str> = a.lines().collect();
    rev.reverse();
    let a_unsorted = fx.plain("ua.bed", &(rev.join("\n") + "\n"));
    let bg = {
        let mut t = String::new();
        for c in 1..=4 {
            for i in 0..30 {
                writeln!(t, "chr{c}\t{}\t{}\t{}", i * 100, i * 100 + 60, 1 + i % 5).unwrap();
            }
        }
        t
    };
    let bg1 = fx.plain("one.bedgraph", &bg);
    let bg2 = fx.plain("two.bedgraph", &bg.replace("\t60\t", "\t80\t"));
    let vcf = p2_vcf(1..=4, 80);
    let vcf_plain = fx.plain("v.vcf", &vcf);
    let vcf_gz = fx.bgzf("v.vcf.gz", &vcf);
    let vcf_idx = fx.indexed("v_idx.vcf.gz", &vcf);
    let vcf_bcf = fx.bcf("v.bcf", &vcf);
    let g = genome.as_str();

    // BED-shaped input: one command per row, every format spelling identical.
    let bed_cmds = p2_bed_cmds(g);
    for (name, tpl) in &bed_cmds {
        let mut per_format = Vec::new();
        for (ai, bi) in [(&a_plain, &b_plain), (&a_gz, &b_gz), (&a_idx, &b_gz)] {
            let args: Vec<&str> = tpl
                .iter()
                .map(|x| match *x {
                    "{a}" => ai.as_str(),
                    "{b}" => bi.as_str(),
                    o => o,
                })
                .collect();
            per_format.push(same_at_every_thread_count(name, &args));
        }
        assert!(!per_format[0].is_empty(), "{name}: empty output");
        assert_eq!(per_format[0], per_format[1], "{name}: plain vs bgzf");
        assert_eq!(
            per_format[0], per_format[2],
            "{name}: plain vs indexed bgzf"
        );
    }

    // --sort on unsorted input is a different engine path.
    for (name, tpl) in [
        ("cluster --sort", vec!["cluster", &a_unsorted, "--sort"]),
        ("spacing --sort", vec!["spacing", &a_unsorted, "--sort"]),
        (
            "partition --sort",
            vec![
                "partition",
                &a_unsorted,
                &b_plain,
                "--sort",
                "--names",
                "A,B",
            ],
        ),
        (
            "genomecov --sort",
            vec!["genomecov", &a_unsorted, "-g", g, "--sort", "--bedgraph"],
        ),
    ] {
        let out = same_at_every_thread_count(name, &tpl);
        assert!(!out.is_empty(), "{name}");
    }
    // bedGraph partition
    let out = same_at_every_thread_count(
        "partition --bedgraph",
        &["partition", "--bedgraph", &bg1, &bg2],
    );
    assert!(!out.is_empty());

    // VCF input: identical across plain, bgzf, indexed and BCF.
    let vcf_cmds: Vec<(&str, Vec<&str>)> = vec![
        ("cluster", vec!["cluster", "{v}"]),
        ("spacing", vec!["spacing", "{v}"]),
        (
            "partition",
            vec!["partition", "{v}", &b_plain, "--names", "A,B"],
        ),
        ("genomecov", vec!["genomecov", "{v}", "-g", g, "--bedgraph"]),
        (
            "extract --info",
            vec!["vcf", "extract", "{v}", "--info", "DP"],
        ),
        (
            "extract --format wide",
            vec!["vcf", "extract", "{v}", "--format", "GT,DP"],
        ),
        (
            "extract --format long",
            vec![
                "vcf", "extract", "{v}", "--format", "GT,DP", "--layout", "long",
            ],
        ),
        (
            "extract filtered",
            vec![
                "vcf",
                "extract",
                "{v}",
                "--info",
                "DP",
                "--min-qual",
                "30",
                "--biallelic",
            ],
        ),
        ("concat", vec!["vcf", "concat", "{v}", "{v2}"]),
    ];
    for cmd in ["cluster", "spacing"] {
        gloci()
            .args([cmd, &vcf_bcf])
            .assert()
            .code(2)
            .stderr(contains("convert first"));
    }
    // concat needs the second file to start where the first ends
    let second = fx.plain(
        "v_tail.vcf",
        &p2_vcf(4..=4, 1).replace("chr4\t100", "chr4\t99999"),
    );
    for (name, tpl) in &vcf_cmds {
        let mut outs = Vec::new();
        // cluster and spacing append a column, which a BCF record has no place for
        let takes_bcf = !matches!(*name, "cluster" | "spacing");
        for v in [&vcf_plain, &vcf_gz, &vcf_idx, &vcf_bcf] {
            if !takes_bcf && *v == vcf_bcf {
                continue;
            }
            let args: Vec<&str> = tpl
                .iter()
                .map(|x| match *x {
                    "{v}" => v.as_str(),
                    "{v2}" => second.as_str(),
                    o => o,
                })
                .collect();
            // concat's header differs between plain and BCF in nothing but the
            // glociCommand line
            outs.push(without_command(&same_at_every_thread_count(name, &args)));
        }
        assert!(!outs[0].is_empty(), "{name}");
        // concat's default output format follows its first input (see
        // p2_concat_output_format), so only its thread determinism is
        // compared here
        if *name != "concat" {
            for w in outs.windows(2) {
                assert_eq!(w[0], w[1], "{name}: input format changed the output");
            }
        }
    }
}

/// `vcf stats`, `vcf extract --genotype-matrix` and `vcf concat` write files.
#[test]
fn p2_determinism_matrix_files() {
    let fx = P2::new();
    let vcf = p2_vcf(1..=4, 120);
    let plain = fx.plain("v.vcf", &vcf);
    let gz = fx.bgzf("v.vcf.gz", &vcf);
    let bcf = fx.bcf("v.bcf", &vcf);
    let flags = [
        "--freq",
        "--counts",
        "--depth",
        "--site-depth",
        "--site-mean-depth",
        "--geno-depth",
        "--site-quality",
        "--missing-site",
        "--missing-indv",
        "--het",
        "--hardy",
        "--singletons",
        "--filter-summary",
        "--tstv-summary",
        "--indel-hist",
        "--summary",
    ];
    let mut reference: Option<Vec<(String, Vec<u8>)>> = None;
    for input in [&plain, &gz, &bcf] {
        for t in ["1", "2", "8"] {
            let d = tempdir();
            let prefix = s_of(&d.path().join("o"));
            let mut args = vec![
                "vcf",
                "stats",
                input.as_str(),
                "--out-prefix",
                &prefix,
                "-t",
                t,
            ];
            args.extend(flags);
            args.extend(["--snp-density", "500"]);
            gloci().args(&args).assert().success();
            let mut files: Vec<(String, Vec<u8>)> = std::fs::read_dir(d.path())
                .unwrap()
                .map(|e| {
                    let e = e.unwrap();
                    (
                        e.file_name().to_string_lossy().into_owned(),
                        std::fs::read(e.path()).unwrap(),
                    )
                })
                .collect();
            files.sort();
            assert!(
                files.len() >= 16,
                "{input} -t {t}: wrote {} files",
                files.len()
            );
            assert!(files.iter().all(|(_, b)| !b.is_empty()));
            match &reference {
                None => reference = Some(files),
                Some(r) => assert_eq!(&files, r, "{input} -t {t}"),
            }
        }
    }
    // --json is a single document and is deterministic too
    let outs: Vec<Vec<u8>> = ["1", "2", "8"]
        .iter()
        .map(|t| {
            let d = tempdir();
            let prefix = s_of(&d.path().join("o"));
            gloci()
                .args([
                    "vcf",
                    "stats",
                    &plain,
                    "--out-prefix",
                    &prefix,
                    "--summary",
                    "--json",
                    "-t",
                    t,
                ])
                .assert()
                .success();
            let mut names: Vec<_> = std::fs::read_dir(d.path())
                .unwrap()
                .map(|e| e.unwrap().path())
                .collect();
            names.sort();
            names
                .iter()
                .flat_map(|p| std::fs::read(p).unwrap())
                .collect()
        })
        .collect();
    assert_eq!(outs[0], outs[1]);
    assert_eq!(outs[1], outs[2]);

    // genotype matrix: three files, identical across thread counts and inputs
    let mut reference: Option<Vec<Vec<u8>>> = None;
    for input in [&plain, &gz, &bcf] {
        for t in ["1", "2", "8"] {
            let d = tempdir();
            let prefix = s_of(&d.path().join("m"));
            gloci()
                .args([
                    "vcf",
                    "extract",
                    input,
                    "--genotype-matrix",
                    "012",
                    "--out-prefix",
                    &prefix,
                    "-t",
                    t,
                ])
                .assert()
                .success();
            let files: Vec<Vec<u8>> = ["012", "012.indv", "012.pos"]
                .iter()
                .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap())
                .collect();
            assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 3);
            assert!(files.iter().all(|f| !f.is_empty()));
            match &reference {
                None => reference = Some(files),
                Some(r) => assert_eq!(&files, r, "{input} -t {t}"),
            }
        }
    }

    // extract to a file; concat checked and --naive to files
    let part1 = fx.bgzf("p1.vcf.gz", &p2_vcf(1..=2, 120));
    let part2 = fx.bgzf("p2.vcf.gz", &p2_vcf(3..=4, 120));
    for (name, extra) in [("checked", vec![]), ("naive", vec!["--naive"])] {
        for ext in ["vcf.gz", "bcf"] {
            if name == "naive" && ext == "bcf" {
                continue;
            }
            let mut outs = Vec::new();
            for t in ["1", "2", "8"] {
                let d = tempdir();
                let out = s_of(&d.path().join(format!("o.{ext}")));
                let mut args = vec!["vcf", "concat", &part1, &part2, "-o", &out, "-t", t];
                args.extend(extra.iter().copied());
                gloci().args(&args).assert().success();
                outs.push(std::fs::read(&out).unwrap());
            }
            assert_eq!(outs[0], outs[1], "concat {name} {ext}: -t 1 vs 2");
            assert_eq!(outs[1], outs[2], "concat {name} {ext}: -t 2 vs 8");
        }
    }
    let mut outs = Vec::new();
    for t in ["1", "2", "8"] {
        let d = tempdir();
        let out = s_of(&d.path().join("o.tsv"));
        gloci()
            .args([
                "vcf", "extract", &gz, "--format", "GT,DP", "-o", &out, "-t", t,
            ])
            .assert()
            .success();
        outs.push(std::fs::read(&out).unwrap());
    }
    assert_eq!(outs[0], outs[1]);
    assert_eq!(outs[1], outs[2]);
}

#[test]
fn p2_broken_pipe_matrix() {
    let d = tempdir();
    let genome = plain_file(d.path(), "g.genome", "chr1\t1000000\n");
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..20_000 {
        let s = i * 20;
        writeln!(a, "chr1\t{s}\t{}\ta{i}\t100\t+", s + 15).unwrap();
        writeln!(b, "chr1\t{}\t{}\tb{i}\t50\t+", s + 5, s + 25).unwrap();
    }
    let mut v = String::from(P2_VCF_HEAD);
    for i in 0..20_000 {
        writeln!(
            v,
            "chr1\t{}\t.\tA\tC\t50\tPASS\tDP=30\tGT:DP\t0/1:5\t1/1:6\t0/0:7",
            1 + i * 20
        )
        .unwrap();
    }
    let a = s_of(&plain_file(d.path(), "a.bed", &a));
    let b = s_of(&plain_file(d.path(), "b.bed", &b));
    let v = s_of(&plain_file(d.path(), "v.vcf", &v));
    let tail = tail_after_big_vcf(d.path());
    let g = s_of(&genome);

    let cmds: Vec<Vec<&str>> = vec![
        vec!["shift", &a, "--by", "3"],
        vec!["cluster", &a],
        vec!["spacing", &a],
        vec!["groupby", &a, "-g", "4", "-c", "5:sum"],
        vec!["partition", &a, &b],
        vec!["partition", &a, &b, "--no-membership"],
        vec!["genomecov", &a, "-g", &g, "--per-base"],
        vec!["genomecov", &a, "-g", &g, "--bga"],
        vec!["vcf", "extract", &v, "--info", "DP"],
        vec![
            "vcf", "extract", &v, "--format", "GT,DP", "--layout", "long",
        ],
        vec!["vcf", "concat", &v, &tail],
        vec!["cluster", &v],
        vec!["spacing", &v],
    ];
    for args in &cmds {
        for t in ["1", "4"] {
            let mut child = Command::new(bin())
                .args(args)
                .args(["-t", t])
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            let mut first = [0u8; 64];
            let _ = child.stdout.as_mut().unwrap().read(&mut first);
            drop(child.stdout.take());
            let out = child.wait_with_output().unwrap();
            assert_eq!(
                out.status.code(),
                Some(0),
                "{args:?} -t {t}: stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert!(
                out.stderr.is_empty(),
                "{args:?} -t {t} printed on stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
    }
}

#[test]
fn p2_empty_inputs_matrix() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t1000\nchr2\t1000\n");
    let a = fx.plain(
        "a.bed",
        "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n",
    );
    let empty = fx.plain("empty.bed", "");
    let empty_gz = fx.bgzf("empty.bed.gz", "");
    let hdr_only = fx.plain("hdr.vcf", &p2_vcf(1..=1, 0));
    let hdr_only_gz = fx.bgzf("hdr.vcf.gz", &p2_vcf(1..=1, 0));
    let hdr_only_bcf = fx.bcf("hdr.bcf", &p2_vcf(1..=1, 0));
    let one = fx.plain("one.vcf", &p2_vcf(1..=1, 3));
    let (g, a) = (genome.as_str(), a.as_str());

    // empty BED (zero bytes, and an empty BGZF stream): no records out
    for e in [&empty, &empty_gz] {
        let silent: Vec<Vec<&str>> = vec![
            vec!["shift", e, "--by", "5"],
            vec!["cluster", e],
            vec!["spacing", e],
            vec!["groupby", e, "-g", "1", "-c", "2:sum"],
            vec!["partition", e, "--names", "A", "--min-files", "1"],
            vec!["genomecov", e, "-g", g, "--bedgraph"],
        ];
        for args in &silent {
            for t in ["1", "8"] {
                assert!(stdout_with(args, t).is_empty(), "{args:?} -t {t}");
            }
        }
        // genomecov: an empty input still covers every chromosome with zero
        assert_eq!(
            String::from_utf8(run_stdout(&["genomecov", e, "-g", g, "--bga"])).unwrap(),
            "chr1\t0\t1000\t0\nchr2\t0\t1000\t0\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&["genomecov", e, "-g", g, "--hist"])).unwrap(),
            "chr1\t0\t1000\t1000\t1\nchr2\t0\t1000\t1000\t1\ngenome\t0\t2000\t2000\t1\n"
        );
        assert!(run_stdout(&["genomecov", e, "-g", g, "--per-base"]).starts_with(b"chr1\t1\t0\n"));
        // one empty and one non-empty input to partition
        assert_eq!(
            String::from_utf8(run_stdout(&["partition", e, a, "--names", "E,A"])).unwrap(),
            "chr1\t10\t50\t1\tA\nchr1\t60\t100\t1\tA\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&["partition", a, e, "--names", "A,E"])).unwrap(),
            "chr1\t10\t50\t1\tA\nchr1\t60\t100\t1\tA\n"
        );
        assert!(run_stdout(&["partition", e, e, "--names", "X,Y"]).is_empty());
        // -o with empty input: the file exists and is empty
        let d = tempdir();
        let out = s_of(&d.path().join("o.bed"));
        gloci().args(["cluster", e, "-o", &out]).assert().success();
        assert_eq!(std::fs::metadata(&out).unwrap().len(), 0);
    }

    // header-only VCF: cluster/spacing keep the column header, nothing else
    for v in [&hdr_only, &hdr_only_gz] {
        for cmd in ["cluster", "spacing"] {
            let out = String::from_utf8(run_stdout(&[cmd, v])).unwrap();
            assert_eq!(out.lines().count(), 1, "{cmd}: {out}");
            assert!(out.starts_with("#CHROM\tPOS"), "{out}");
        }
        assert!(run_stdout(&["partition", v, "--names", "V"]).is_empty());
        assert_eq!(
            String::from_utf8(run_stdout(&["vcf", "extract", v, "--info", "DP"])).unwrap(),
            "CHROM\tPOS\tID\tREF\tALT\tDP\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&["vcf", "extract", v, "--format", "GT"])).unwrap(),
            "CHROM\tPOS\tID\tREF\tALT\ts1:GT\ts2:GT\ts3:GT\n"
        );
        assert_eq!(
            String::from_utf8(run_stdout(&[
                "vcf", "extract", v, "--format", "GT", "--layout", "long"
            ]))
            .unwrap()
            .lines()
            .count(),
            1
        );
    }
    assert_eq!(
        String::from_utf8(run_stdout(&[
            "vcf",
            "extract",
            &hdr_only_bcf,
            "--info",
            "DP"
        ]))
        .unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\tDP\n"
    );
    assert_eq!(
        String::from_utf8(run_stdout(&[
            "vcf",
            "extract",
            &hdr_only_bcf,
            "--format",
            "GT"
        ]))
        .unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\ts1:GT\ts2:GT\ts3:GT\n"
    );

    // vcf stats on a header-only VCF: every file is written, with its header
    let d = tempdir();
    let prefix = s_of(&d.path().join("o"));
    gloci()
        .args([
            "vcf",
            "stats",
            &hdr_only,
            "--out-prefix",
            &prefix,
            "--freq",
            "--counts",
            "--depth",
            "--site-depth",
            "--site-mean-depth",
            "--geno-depth",
            "--site-quality",
            "--missing-site",
            "--missing-indv",
            "--het",
            "--hardy",
            "--singletons",
            "--filter-summary",
            "--tstv-summary",
            "--indel-hist",
            "--summary",
        ])
        .assert()
        .success();
    let written: Vec<String> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert!(written.len() >= 14, "{written:?}");
    let imiss = std::fs::read_to_string(format!("{prefix}.imiss")).unwrap();
    assert_eq!(
        imiss.lines().count(),
        4,
        "one line per sample plus header: {imiss}"
    );
    let frq = std::fs::read_to_string(format!("{prefix}.frq")).unwrap();
    assert_eq!(frq.lines().count(), 1, "{frq}");

    // genotype matrix of a header-only VCF: empty matrix, samples listed
    let prefix = s_of(&d.path().join("m"));
    gloci()
        .args([
            "vcf",
            "extract",
            &hdr_only,
            "--genotype-matrix",
            "012",
            "--out-prefix",
            &prefix,
        ])
        .assert()
        .success();
    assert_eq!(
        std::fs::read(format!("{prefix}.012")).unwrap(),
        b"0\n1\n2\n"
    );
    assert_eq!(std::fs::read(format!("{prefix}.012.pos")).unwrap(), b"");
    assert_eq!(
        std::fs::read_to_string(format!("{prefix}.012.indv")).unwrap(),
        "s1\ns2\ns3\n"
    );

    // concat: empty files are skipped without changing the result
    let want = String::from_utf8(without_command(&run_stdout(&["vcf", "concat", &one]))).unwrap();
    for args in [
        vec!["vcf", "concat", &hdr_only, &one],
        vec!["vcf", "concat", &one, &hdr_only],
        vec!["vcf", "concat", &hdr_only, &one, &hdr_only],
    ] {
        let got = without_command(&run_stdout(&args));
        assert_eq!(String::from_utf8(got).unwrap(), want, "{args:?}");
    }
    let both = String::from_utf8(run_stdout(&["vcf", "concat", &hdr_only, &hdr_only])).unwrap();
    assert!(
        both.contains("#CHROM") && both.lines().all(|l| l.starts_with('#')),
        "{both}"
    );
    let naive_empty = fx.bgzf("hdr_n.vcf.gz", &p2_vcf(1..=1, 0));
    let d2 = tempdir();
    let out = s_of(&d2.path().join("n.vcf.gz"));
    gloci()
        .args([
            "vcf",
            "concat",
            "--naive",
            &naive_empty,
            &hdr_only_gz,
            "-o",
            &out,
        ])
        .assert()
        .success();
}

#[test]
fn p2_vcf_input_matrix() {
    let fx = P2::new();
    // an indel and a symbolic deletion: the footprints differ by --vcf-span
    let vcf = format!(
        "{P2_VCF_HEAD}chr1\t100\t.\tA\tC\t50\tPASS\tDP=3\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n\
chr1\t103\t.\tATGC\tA\t50\tPASS\tDP=3\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n\
chr1\t500\t.\tA\t<DEL>\t50\tPASS\tDP=3;END=700\tGT:DP\t0/1:5\t1/1:6\t0/0:7\n"
    );
    let v = fx.plain("v.vcf", &vcf);
    let g = fx.plain("g.genome", "chr1\t10000\n");
    // cluster appends a TSV column, keeps the record text and the header
    let out = String::from_utf8(run_stdout(&["cluster", &v])).unwrap();
    let lines: Vec<&str> = out.lines().collect();
    assert!(lines[0].starts_with("#CHROM") && lines[0].ends_with("\tcluster"));
    assert!(lines[1].ends_with("\t1") && lines[2].ends_with("\t2") && lines[3].ends_with("\t3"));
    assert!(!out.contains("##"), "TSV output carries no ## meta lines");
    // spacing: gap between footprints
    let out = String::from_utf8(run_stdout(&["spacing", &v])).unwrap();
    let last: Vec<&str> = out
        .lines()
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    assert_eq!(last, ["spacing", ".", "3", "394"], "{out}");
    // genomecov counts footprints: the span mode decides what a record covers
    let cov = |span: &str| {
        String::from_utf8(run_stdout(&[
            "genomecov",
            &v,
            "-g",
            &g,
            "--bedgraph",
            "--vcf-span",
            span,
        ]))
        .unwrap()
    };
    assert_eq!(
        cov("allele"),
        "chr1\t99\t100\t1\nchr1\t103\t106\t1\nchr1\t500\t700\t1\n"
    );
    assert_eq!(
        cov("ref"),
        "chr1\t99\t100\t1\nchr1\t102\t106\t1\nchr1\t499\t700\t1\n"
    );
    // partition labels the VCF by file name and lists breakpoints
    let out = String::from_utf8(run_stdout(&["partition", &v, "--names", "V"])).unwrap();
    assert!(out.lines().all(|l| l.ends_with("\tV")), "{out}");
    // --sort accepts an unsorted VCF
    let mut rows: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
    rows.reverse();
    let uv = fx.plain("u.vcf", &format!("{P2_VCF_HEAD}{}\n", rows.join("\n")));
    assert_eq!(
        run_stdout(&["cluster", &uv, "--sort"]),
        run_stdout(&["cluster", &v])
    );
    // shift and groupby do not take VCF, and say what to do
    gloci()
        .args(["shift", &v, "--by", "1"])
        .assert()
        .code(2)
        .stderr(contains("gloci view -O bed"));
    gloci()
        .args(["groupby", &v, "-g", "1", "-c", "2:count"])
        .assert()
        .code(2)
        .stderr(contains("gloci vcf extract"));
    // vcf stats and extract take --samples and filters on the VCF
    let d = tempdir();
    let prefix = s_of(&d.path().join("o"));
    gloci()
        .args([
            "vcf",
            "stats",
            &v,
            "--out-prefix",
            &prefix,
            "--missing-indv",
            "--samples",
            "s2",
        ])
        .assert()
        .success();
    let imiss = std::fs::read_to_string(format!("{prefix}.imiss")).unwrap();
    assert_eq!(imiss.lines().count(), 2, "{imiss}");
    let out = String::from_utf8(run_stdout(&[
        "vcf",
        "extract",
        &v,
        "--format",
        "GT",
        "--samples",
        "s3",
        "--biallelic",
    ]))
    .unwrap();
    assert_eq!(
        out.lines().next().unwrap(),
        "CHROM\tPOS\tID\tREF\tALT\ts3:GT"
    );
    assert_eq!(out.lines().count(), 4, "{out}");
}

#[test]
fn p2_failure_matrix() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t10000\nchr2\t10000\n");
    let good_bed = fx.plain("good.bed", "chr1\t10\t50\ta\t1\t+\nchr1\t60\t90\tb\t2\t-\n");
    let other_bed = fx.plain("other.bed", "chr1\t20\t70\tc\t3\t+\n");
    let good_vcf_text = p2_vcf(1..=2, 4);
    let good_vcf = fx.plain("good.vcf", &good_vcf_text);
    let good_vcf_gz = fx.bgzf("good.vcf.gz", &good_vcf_text);

    let unsorted_bed = fx.plain(
        "u.bed",
        "chr1\t100\t200\tr1\t1\t+\nchr1\t50\t150\tr2\t2\t+\n",
    );
    let malformed_bed = fx.plain(
        "m.bed",
        "chr1\t100\t200\tr1\t1\t+\nchr1\tNOPE\t250\tr2\t2\t+\n",
    );
    let bed_gz = bgzip(bed_a(1).as_bytes());
    let trunc_bed = s_of(&write(
        fx.dir.path(),
        "t.bed.gz",
        &bed_gz[..bed_gz.len() / 2],
    ));
    let nofooter_bed = s_of(&write(
        fx.dir.path(),
        "nf.bed.gz",
        &bed_gz[..bed_gz.len() - 28],
    ));

    let rows: Vec<&str> = good_vcf_text.lines().collect();
    let hdr = rows.iter().filter(|l| l.starts_with('#')).count();
    let body: Vec<&str> = rows[hdr..].to_vec();
    let swapped = {
        let mut b = body.clone();
        b.swap(1, 2);
        fx.plain("u.vcf", &format!("{P2_VCF_HEAD}{}\n", b.join("\n")))
    };
    let malformed_vcf = fx.plain(
        "m.vcf",
        &format!(
            "{P2_VCF_HEAD}{}\n{}\n",
            body[0],
            body[1].replacen("\t137\t", "\tXX\t", 1)
        ),
    );
    let good_bcf = vcf_to_bcf(good_vcf_text.as_bytes());
    let trunc_bcf = s_of(&write(
        fx.dir.path(),
        "t.bcf",
        &good_bcf[..good_bcf.len() / 2],
    ));
    let vgz = bgzip(good_vcf_text.as_bytes());
    let trunc_vcf = s_of(&write(fx.dir.path(), "t.vcf.gz", &vgz[..vgz.len() / 2]));
    let nofooter_vcf = s_of(&write(fx.dir.path(), "nf.vcf.gz", &vgz[..vgz.len() - 28]));

    let after_gz = fx.bgzf("after.vcf.gz", &p2_vcf(3..=4, 4));
    let before_gz = fx.bgzf("before.vcf.gz", &p2_vcf(1..=1, 0));
    let g = genome.as_str();
    let (gb, ob, gv, gvz) = (&good_bed, &other_bed, &good_vcf, &good_vcf_gz);
    let unsorted_msg = "is not sorted";
    let eof = "BGZF end-of-file marker missing";

    // (name, args with {IN} marking the bad input, sink, per-kind messages)
    type Case<'a> = (&'a str, Vec<&'a str>, Sink);
    let cases: Vec<Case> = vec![
        (
            "shift",
            vec!["shift", "{IN}", "--by", "3"],
            Sink::File("bed"),
        ),
        ("cluster", vec!["cluster", "{IN}"], Sink::File("tsv")),
        ("spacing", vec!["spacing", "{IN}"], Sink::File("tsv")),
        (
            "groupby",
            vec!["groupby", "{IN}", "-g", "1", "-c", "2:sum"],
            Sink::File("tsv"),
        ),
        (
            "partition",
            vec!["partition", "{IN}", ob, "--names", "A,B"],
            Sink::File("bed"),
        ),
        (
            "partition B",
            vec!["partition", gb, "{IN}", "--names", "A,B"],
            Sink::File("bed"),
        ),
        (
            "genomecov",
            vec!["genomecov", "{IN}", "-g", g],
            Sink::File("tsv"),
        ),
    ];
    // BED-shaped failures: (input, needs sorted input, message)
    let bed_bad: [(&str, &String, &str); 4] = [
        ("unsorted", &unsorted_bed, unsorted_msg),
        ("malformed", &malformed_bed, "NOPE"),
        ("truncated", &trunc_bed, eof),
        ("no BGZF footer", &nofooter_bed, eof),
    ];
    let mut ran = 0;
    for (name, tpl, sink) in &cases {
        for (kind, bad, msg) in &bed_bad {
            // shift and groupby never check the order
            // shift and groupby never check the order, and must keep working
            if *kind == "unsorted" && matches!(*name, "shift" | "groupby") {
                assert_succeeds(name, tpl, bad, sink);
                continue;
            }
            // groupby sums column 2, so the malformed row is a bad number there
            let want = if *name == "groupby" && *kind == "malformed" {
                "is not a number"
            } else {
                msg
            };
            assert_fails_clean(name, kind, tpl, bad, sink, want);
            ran += 1;
        }
        // missing output directory, on good input
        assert_missing_dir(name, tpl, gb, sink);
    }

    // VCF-shaped failures
    let vcf_cases: Vec<Case> = vec![
        ("cluster", vec!["cluster", "{IN}"], Sink::File("tsv")),
        ("spacing", vec!["spacing", "{IN}"], Sink::File("tsv")),
        (
            "partition",
            vec!["partition", "{IN}", ob, "--names", "A,B"],
            Sink::File("bed"),
        ),
        (
            "genomecov",
            vec!["genomecov", "{IN}", "-g", g],
            Sink::File("tsv"),
        ),
        (
            "extract --info",
            vec!["vcf", "extract", "{IN}", "--info", "DP"],
            Sink::File("tsv"),
        ),
        (
            "extract --format",
            vec!["vcf", "extract", "{IN}", "--format", "GT"],
            Sink::File("tsv"),
        ),
        (
            "extract --012",
            vec!["vcf", "extract", "{IN}", "--genotype-matrix", "012"],
            Sink::Prefix,
        ),
        (
            "stats --freq",
            vec!["vcf", "stats", "{IN}", "--freq"],
            Sink::Prefix,
        ),
        (
            "stats --summary",
            vec![
                "vcf",
                "stats",
                "{IN}",
                "--summary",
                "--het",
                "--missing-indv",
            ],
            Sink::Prefix,
        ),
        (
            "concat first",
            vec!["vcf", "concat", "{IN}", &after_gz],
            Sink::File("vcf.gz"),
        ),
        (
            "concat second",
            vec!["vcf", "concat", &before_gz, "{IN}"],
            Sink::File("vcf.gz"),
        ),
    ];
    // extract and stats read an unsorted VCF without complaint (extract
    // writes in input order; stats accumulates), so they must succeed on it:
    // adding an order check to either would be a visible decision here
    let vcf_bad: [(&str, &String, &str); 4] = [
        ("unsorted", &swapped, unsorted_msg),
        ("malformed", &malformed_vcf, "XX"),
        ("truncated", &trunc_vcf, eof),
        ("no BGZF footer", &nofooter_vcf, eof),
    ];
    for (name, tpl, sink) in &vcf_cases {
        for (kind, bad, msg) in &vcf_bad {
            let order_free = name.starts_with("extract") || name.starts_with("stats");
            if *kind == "unsorted" && order_free {
                assert_succeeds(name, tpl, bad, sink);
                continue;
            }
            assert_fails_clean(name, kind, tpl, bad, sink, msg);
            ran += 1;
        }
        assert_missing_dir(name, tpl, gv, sink);
    }
    // truncated BCF, for the commands that take BCF
    for (name, tpl, sink) in vcf_cases.iter().filter(|(n, ..)| {
        n.starts_with("extract") || n.starts_with("stats") || n.starts_with("concat")
    }) {
        assert_fails_clean(name, "truncated BCF", tpl, &trunc_bcf, sink, eof);
        ran += 1;
    }
    // a BED where a VCF is required says so
    gloci()
        .args([
            "vcf",
            "stats",
            gb,
            "--freq",
            "--out-prefix",
            &s_of(&fx.dir.path().join("x")),
        ])
        .assert()
        .failure()
        .stderr(contains("BED"));
    // concat: the second file starting before the first ends, and mismatched
    // contigs on --naive input
    let d = tempdir();
    let out = s_of(&d.path().join("o.vcf.gz"));
    gloci()
        .args(["vcf", "concat", gvz, gvz, "-o", &out])
        .assert()
        .code(1)
        .stderr(contains("each must start where the one before ends"));
    assert!(!Path::new(&out).exists());
    gloci()
        .args(["vcf", "concat", "--naive", gvz, gvz, "-o", &out])
        .assert()
        .failure()
        .stderr(contains("is not sorted"));
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
    gloci()
        .args(["vcf", "concat", "--naive", gv, gvz, "-o", &out])
        .assert()
        .code(2)
        .stderr(contains("compress it with: gloci view"));
    assert!(ran > 60, "ran {ran} failure cases");
}
