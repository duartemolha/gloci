//! Cross-command behaviour: thread determinism, broken pipe, truncation,
//! atomic outputs, signals and the error-message format.
mod alignments_matrices;
mod bed_vcf_matrices;
mod bed_vcf_outputs;
#[path = "../common/mod.rs"]
mod common;
mod fixtures_alignments;
mod fixtures_bed_vcf;
mod fixtures_samtools;
mod fixtures_stats;
mod helpers;
mod matrices;
mod process;
use crate::fixtures_alignments::*;
use crate::fixtures_bed_vcf::*;
use crate::fixtures_samtools::*;
use crate::fixtures_stats::*;
use crate::helpers::*;
use common::*;
use gloci_testkit::fixtures::{bgzip, write};
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;
use std::process::{Command, Stdio};

/// `--chrom-map` is a rename applied as records are read: running with the
/// map on a file in the source names equals running without it on the file
/// already in the target names, and the inverse map takes the output back.
#[test]
fn p2_chrom_map_round_trip() {
    let fx = P2::new();
    let fwd = fx.plain("fwd.map", "1\tchr1\n2\tchr2\n");
    let inv = fx.plain("inv.map", "chr1\t1\nchr2\t2\n");
    let genome_chr = fx.plain("g_chr.genome", "chr1\t10000\nchr2\t10000\n");
    let chr_bed = bed_a(2);
    let num_bed = chr_bed.replace("chr", "");
    let bed_chr = fx.plain("chr.bed", &chr_bed);
    let bed_num = fx.plain("num.bed", &num_bed);
    let other_chr = fx.plain("o_chr.bed", &bed_b(2));
    let other_num = fx.plain("o_num.bed", &bed_b(2).replace("chr", ""));
    let vcf_chr = p2_vcf(1..=2, 30);
    let vcf_num = vcf_chr.replace("chr", "");
    let vchr = fx.plain("chr.vcf", &vcf_chr);
    let vnum = fx.plain("num.vcf", &vcf_num);
    let vnum_gz = fx.bgzf("num.vcf.gz", &vcf_num);
    let vnum_bcf = fx.bcf("num.bcf", &vcf_num);

    // commands that write input lines back: the output round-trips exactly
    for cmd in ["cluster", "spacing"] {
        for t in ["1", "8"] {
            let renamed = stdout_with(&[cmd, &bed_num, "--chrom-map", &fwd], t);
            assert_eq!(
                renamed,
                stdout_with(&[cmd, &bed_chr], t),
                "{cmd} map vs renamed"
            );
            let back = stdout_with(&[cmd, &bed_chr, "--chrom-map", &inv], t);
            assert_eq!(
                back,
                stdout_with(&[cmd, &bed_num], t),
                "{cmd}: the inverse map restores the source names"
            );
        }
    }
    // the others
    let g = genome_chr.as_str();
    let cases: Vec<Vec<&str>> = vec![
        vec!["genomecov", "{I}", "-g", g, "--bedgraph"],
        vec!["genomecov", "{I}", "-g", g],
        vec!["partition", "{I}", "{O}", "--names", "A,B"],
    ];
    for tpl in &cases {
        let build = |i: &str, o: &str, map: bool| -> Vec<String> {
            let mut v: Vec<String> = tpl
                .iter()
                .map(|a| match *a {
                    "{I}" => i.to_string(),
                    "{O}" => o.to_string(),
                    x => x.to_string(),
                })
                .collect();
            if map {
                v.extend(["--chrom-map".into(), fwd.clone()]);
            }
            v
        };
        let run = |args: Vec<String>| {
            let refs: Vec<&str> = args.iter().map(String::as_str).collect();
            same_at_every_thread_count("chrom-map", &refs)
        };
        let with_map = run(build(&bed_num, &other_num, true));
        let renamed = run(build(&bed_chr, &other_chr, false));
        assert!(!with_map.is_empty());
        assert_eq!(with_map, renamed, "{tpl:?}");
        assert!(
            !String::from_utf8_lossy(&with_map)
                .lines()
                .any(|l| l.starts_with("1\t"))
        );
    }
    // a name the map does not list passes through; mixed names work
    let mixed = fx.plain("mixed.bed", "1\t10\t20\nchr2\t5\t9\nchrX\t1\t2\n");
    let out = String::from_utf8(run_stdout(&[
        "cluster",
        &mixed,
        "--chrom-map",
        &fwd,
        "--sort",
    ]))
    .unwrap();
    assert_eq!(out, "chr1\t10\t20\t1\nchr2\t5\t9\t2\nchrX\t1\t2\t3\n");

    // VCF input (plain, BGZF and BCF carry the source names in the header):
    // extract and stats write the target names
    for input in [&vnum, &vnum_gz, &vnum_bcf] {
        let mapped = stdout_with(
            &["vcf", "extract", input, "--info", "DP", "--chrom-map", &fwd],
            "2",
        );
        let renamed = stdout_with(&["vcf", "extract", &vchr, "--info", "DP"], "2");
        assert_eq!(mapped, renamed, "extract {input}");
        assert!(mapped.starts_with(b"CHROM\tPOS\tID\tREF\tALT\tDP\nchr1\t"));

        let files = |vcf: &str, map: bool, t: &str| -> Vec<Vec<u8>> {
            let d = tempdir();
            let prefix = s_of(&d.path().join("o"));
            let mut args = vec![
                "vcf",
                "stats",
                vcf,
                "--out-prefix",
                &prefix,
                "--freq",
                "--site-quality",
                "--singletons",
                "--snp-density",
                "500",
                "-t",
                t,
            ];
            if map {
                args.extend(["--chrom-map", &fwd]);
            }
            gloci().args(&args).assert().success();
            ["frq", "lqual", "singletons", "snpden"]
                .iter()
                .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap())
                .collect()
        };
        let a = files(input, true, "1");
        assert_eq!(a, files(input, true, "8"), "stats {input}: threads");
        assert_eq!(a, files(&vchr, false, "1"), "stats {input}: map vs renamed");
        assert!(a[0].windows(5).any(|w| w == b"\nchr1"), "{input}");

        let m = |vcf: &str, map: bool| {
            let d = tempdir();
            let prefix = s_of(&d.path().join("m"));
            let mut args = vec![
                "vcf",
                "extract",
                vcf,
                "--genotype-matrix",
                "012",
                "--out-prefix",
                &prefix,
            ];
            if map {
                args.extend(["--chrom-map", &fwd]);
            }
            gloci().args(&args).assert().success();
            std::fs::read(format!("{prefix}.012.pos")).unwrap()
        };
        assert_eq!(m(input, true), m(&vchr, false), "012.pos {input}");
    }
    // cluster/spacing on a VCF with the map: the renamed record text
    let out = String::from_utf8(run_stdout(&["cluster", &vnum, "--chrom-map", &fwd])).unwrap();
    assert_eq!(
        out,
        String::from_utf8(run_stdout(&["cluster", &vchr])).unwrap()
    );
    // a map file with a bad line says where
    let bad = fx.plain("bad.map", "1\tchr1\nlonely\n");
    gloci()
        .args(["cluster", &bed_num, "--chrom-map", &bad])
        .assert()
        .code(2)
        .stderr(contains("bad.map:2:"));
}

/// `--chrom-map` renames as records are read: BAM and CRAM headers and records come out with the
/// new names, the inverse map takes them back, and getfasta, nuc and maskfasta find intervals
/// written with the other naming.
#[test]
fn p3a_chrom_map_round_trip() {
    let fx = P3a::new();
    let d = fx.path("cm");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let fwd = s_of(&plain_file(
        d,
        "fwd.tsv",
        "# to Ensembl\nchr1\t1\nchr2\t2\nchrM\tMT\n",
    ));
    let rev = s_of(&plain_file(d, "rev.tsv", "1\tchr1\n2\tchr2\nMT\tchrM\n"));
    let sam_names = |b: &[u8]| -> Vec<String> {
        String::from_utf8_lossy(&p3a_norm(b))
            .lines()
            .filter_map(|l| l.strip_prefix("@SQ\tSN:"))
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    };
    let ref_args = ["--reference", fx.fa.as_str()];
    for (kind, input, _) in fx.kinds() {
        if kind == "sam" {
            // SAM text written as it is cannot have its @SQ names rewritten: refused, not
            // half-renamed
            gloci()
                .args(["view", &input, "--chrom-map", &fwd])
                .assert()
                .code(2)
                .stderr(contains(
                    "--chrom-map cannot rename the reference names in SAM text",
                ));
            continue;
        }
        let to_bam = |src: &str, map: &str, name: &str, extra: &[&str]| -> String {
            let out = s_of(&d.join(name));
            gloci()
                .args(["view", src, "--chrom-map", map, "-o", &out])
                .args(extra)
                .assert()
                .success();
            out
        };
        let ex: Vec<&str> = if kind == "cram" {
            ref_args.to_vec()
        } else {
            vec![]
        };
        let renamed = to_bam(&input, &fwd, &format!("{kind}-ens.bam"), &ex);
        let hdr = |p: &str| {
            let o = gloci()
                .args(["view", p, "--header-only", "-O", "sam"])
                .output()
                .unwrap();
            sam_names(&o.stdout)
        };
        assert_eq!(hdr(&renamed), ["MT", "1", "2", "chr10", "chrX"], "{kind}");
        let recs = |p: &str| {
            p3a_records(
                &gloci()
                    .args(["view", p, "-O", "sam"])
                    .output()
                    .unwrap()
                    .stdout,
            )
        };
        let ens = recs(&renamed);
        assert!(
            ens.iter().any(|r| r.split('\t').nth(2) == Some("MT")),
            "{kind}"
        );
        assert!(
            ens.iter()
                .all(|r| !r.split('\t').nth(2).unwrap().starts_with("chr1\t")),
            "{kind}"
        );
        // and back: the original names, the same reads, in a BAM that view reads
        let back = to_bam(&renamed, &rev, &format!("{kind}-back.bam"), &[]);
        assert_eq!(
            hdr(&back),
            ["chrM", "chr1", "chr2", "chr10", "chrX"],
            "{kind}"
        );
        let orig = p3a_records(
            &gloci()
                .args(["view", &input, "-O", "sam"])
                .args(&ex)
                .output()
                .unwrap()
                .stdout,
        );
        assert_eq!(recs(&back), orig, "{kind}: round trip changed the reads");
        // other commands see the new names: coverage against a BED in the new names, genomecov
        // rows in the new names, equal to the original run under the inverse naming
        let ens_bed = s_of(&plain_file(
            d,
            "ens.bed",
            "MT\t0\t100\n1\t0\t300\n1\t600\t900\n2\t100\t900\n",
        ));
        let want = {
            let bed = s_of(&plain_file(
                d,
                "orig.bed",
                "chrM\t0\t100\nchr1\t0\t300\nchr1\t600\t900\nchr2\t100\t900\n",
            ));
            let o = gloci()
                .args(["coverage", "-a", &bed, "-b", &input])
                .args(&ex)
                .output()
                .unwrap();
            assert!(o.status.success());
            String::from_utf8_lossy(&o.stdout)
                .replace("chrM\t", "MT\t")
                .replace("chr1\t", "1\t")
                .replace("chr2\t", "2\t")
        };
        for t in ["1", "4"] {
            let o = gloci()
                .args([
                    "coverage",
                    "-a",
                    &ens_bed,
                    "-b",
                    &input,
                    "--chrom-map",
                    &fwd,
                    "-t",
                    t,
                ])
                .args(&ex)
                .output()
                .unwrap();
            // the map renames B's reads; A is already in the new names
            assert_eq!(
                o.status.code(),
                Some(0),
                "{kind}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            assert_eq!(String::from_utf8_lossy(&o.stdout), want, "{kind} -t {t}");
            let g = gloci()
                .args(["genomecov", &input, "--chrom-map", &fwd, "-t", t])
                .args(&ex)
                .output()
                .unwrap();
            let text = String::from_utf8_lossy(&g.stdout).into_owned();
            assert!(
                text.contains("MT\t0\t") && text.contains("\n1\t0\t"),
                "{kind}: {text}"
            );
            assert!(
                !text.contains("chr1\t") && !text.contains("chrM\t"),
                "{kind}: {text}"
            );
        }
    }
    // a map that would merge two contigs is refused
    let bad = s_of(&plain_file(d, "bad.tsv", "chr1\tchr2\n"));
    gloci()
        .args(["view", &fx.bam, "--chrom-map", &bad])
        .assert()
        .code(2)
        .stderr(contains("which would merge them"));

    // FASTA: intervals in the other naming find their sequences through the map, and the
    // output uses the FASTA's names
    let plain = {
        let bed = s_of(&plain_file(
            d,
            "plain.bed",
            "chr1\t10\t40\tx\t0\t+\nchr2\t100\t130\ty\t0\t-\nchrM\t5\t25\tz\t0\t+\n",
        ));
        let o = gloci()
            .args(["getfasta", "--fasta", &fx.fa, &bed, "-s", "--name"])
            .output()
            .unwrap();
        assert!(o.status.success());
        o.stdout
    };
    let ens_iv = s_of(&plain_file(
        d,
        "ens-iv.bed",
        "1\t10\t40\tx\t0\t+\n2\t100\t130\ty\t0\t-\nMT\t5\t25\tz\t0\t+\n",
    ));
    let to_fa_names = s_of(&plain_file(d, "to-fa.tsv", "1\tchr1\n2\tchr2\nMT\tchrM\n"));
    for t in ["1", "4"] {
        let o = gloci()
            .args([
                "getfasta",
                "--fasta",
                &fx.fa,
                &ens_iv,
                "-s",
                "--name",
                "--chrom-map",
                &to_fa_names,
                "-t",
                t,
            ])
            .output()
            .unwrap();
        assert_eq!(
            o.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert_eq!(o.stdout, plain, "getfasta -t {t}");
        // without the map the names are not in the FASTA
        gloci()
            .args(["getfasta", "--fasta", &fx.fa, &ens_iv, "-t", t])
            .assert()
            .code(1)
            .stderr(contains("is not in the FASTA index"));
    }
    let nuc_plain = {
        let bed = s_of(&plain_file(
            d,
            "nplain.bed",
            "chr1\t10\t40\nchr2\t100\t130\n",
        ));
        gloci()
            .args(["nuc", "--fasta", &fx.fa, &bed])
            .output()
            .unwrap()
            .stdout
    };
    let nuc_map = {
        let bed = s_of(&plain_file(d, "nens.bed", "1\t10\t40\n2\t100\t130\n"));
        gloci()
            .args(["nuc", "--fasta", &fx.fa, &bed, "--chrom-map", &to_fa_names])
            .output()
            .unwrap()
            .stdout
    };
    assert_eq!(nuc_plain, nuc_map, "nuc through a map");
    let masks_ens = s_of(&plain_file(
        d,
        "masks-ens.bed",
        "MT\t5\t20\n1\t0\t30\n1\t100\t140\n2\t500\t700\n10\t10\t20\n",
    ));
    let to_fa_names10 = s_of(&plain_file(
        d,
        "to-fa10.tsv",
        "1\tchr1\n2\tchr2\nMT\tchrM\n10\tchr10\n",
    ));
    let want = gloci()
        .args(["maskfasta", "--fasta", &fx.fa, "-b", &fx.masks])
        .output()
        .unwrap()
        .stdout;
    let got = gloci()
        .args([
            "maskfasta",
            "--fasta",
            &fx.fa,
            "-b",
            &masks_ens,
            "--chrom-map",
            &to_fa_names10,
        ])
        .output()
        .unwrap();
    assert_eq!(
        got.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&got.stderr)
    );
    assert_eq!(got.stdout, want, "maskfasta through a map");
}

/// Every P3b command gives the same bytes at `-t 1`, `-t2` and `--threads=8` (the seeded ones
/// with a `--seed`), and the file-writing commands the same files.
#[test]
fn p3b_determinism_matrix() {
    let fx = P3b::new();
    let cmds = fx.stdout_cmds();
    assert!(cmds.len() > 50, "{} commands", cmds.len());
    for (name, args) in &cmds {
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let out = same_at_every_thread_count(name, &refs);
        assert!(!out.is_empty(), "{name}: empty output");
    }
    for (name, args, exts) in p3b_prefix_cmds(&fx) {
        let one = p3b_prefix_files(&args, &exts, "1");
        for t in ["2", "8"] {
            assert_eq!(one, p3b_prefix_files(&args, &exts, t), "{name} -t {t}");
        }
        assert!(
            one.iter().any(|f| f.len() > 100),
            "{name}: nothing was written"
        );
    }
}

/// The output of the seeded commands is a function of `--seed` alone: the same seed gives the
/// same bytes, another seed other bytes, and without a seed the warning names the one chosen,
/// which reproduces the run.
#[test]
fn p3b_seed_matrix() {
    let fx = P3b::new();
    let seeded: Vec<(&str, Vec<String>)> = vec![
        ("random", sv(&["random", "-g", &fx.genome, "-n", "2000"])),
        ("sample", sv(&["sample", "-n", "200", &fx.a])),
        ("sample vcf", sv(&["sample", "-n", "200", &fx.vcf])),
        ("shuffle", sv(&["shuffle", "-g", &fx.genome, &fx.a])),
        (
            "shuffle bedpe",
            sv(&["shuffle", "-g", &fx.genome, "--format", "bedpe", &fx.pe_a]),
        ),
    ];
    for (name, base) in &seeded {
        let run = |seed: &str, t: &str| {
            let mut a = base.clone();
            a.extend(["--seed".into(), seed.into(), "-t".into(), t.into()]);
            let o = gloci().args(&a).output().unwrap();
            assert_eq!(o.status.code(), Some(0), "{name}: {:?}", o.stderr);
            assert!(o.stderr.is_empty(), "{name}: a seed given, no warning");
            o.stdout
        };
        let s1 = run("1", "1");
        assert_eq!(s1, run("1", "8"), "{name}: same seed, other -t");
        assert_eq!(s1, run("1", "1"), "{name}: same seed again");
        assert_ne!(s1, run("2", "1"), "{name}: another seed");
        assert_ne!(run("0", "1"), run("1", "1"), "{name}: seed 0");

        let o = gloci().args(base).args(["-t", "2"]).output().unwrap();
        assert_eq!(o.status.code(), Some(0), "{name}");
        let err = String::from_utf8_lossy(&o.stderr);
        let seed = err
            .strip_prefix("warning: no --seed given; using --seed ")
            .and_then(|r| r.split_whitespace().next())
            .unwrap_or_else(|| panic!("{name}: the warning names a seed: {err}"))
            .to_string();
        assert_eq!(err.matches('\n').count(), 1, "{name}: {err}");
        assert_eq!(
            without_command(&o.stdout),
            without_command(&run(&seed, "1")),
            "{name}: the printed seed reproduces"
        );
        let o = gloci()
            .args(base)
            .args(["-t", "2", "--quiet"])
            .output()
            .unwrap();
        assert!(o.stderr.is_empty(), "{name}: --quiet silences the warning");
    }
}

#[test]
fn p3b_broken_pipe_matrix() {
    let fx = P3b::new();
    let mut ran = 0;
    for (name, args) in fx.stdout_cmds() {
        // the statistics write a few lines; the fixtures of the small cases write little
        let skip = name.starts_with("jaccard")
            || name.starts_with("fisher")
            || name.starts_with("reldist")
            || name.contains("fixture")
            || name.ends_with(".gz")
            || name.starts_with("sample")
            || name == "vcf split --types"
            || name.starts_with("pair-intersect")
            || name.starts_with("pair-pair")
            || name.starts_with("shuffle");
        if skip {
            continue;
        }
        p3b_assert_broken_pipe(&name, &args, 100_000);
        ran += 1;
    }
    // the ones that need more output than the shared fixtures give
    let more: Vec<(&str, Vec<String>)> = vec![
        (
            "sample",
            sv(&["sample", "-n", "6000", "--seed", "1", &fx.a]),
        ),
        (
            "shuffle",
            sv(&["shuffle", "-g", &fx.genome, "--seed", "1", &fx.a]),
        ),
        (
            "shuffle bedpe",
            sv(&[
                "shuffle", "-g", &fx.genome, "--seed", "1", "--format", "bedpe", &fx.pe_a,
            ]),
        ),
        (
            "pair-intersect",
            sv(&["pair-intersect", "-a", &fx.pe_a, "-b", &fx.b]),
        ),
        (
            "pair-intersect neither",
            sv(&[
                "pair-intersect",
                "-a",
                &fx.pe_a,
                "-b",
                &fx.b,
                "--type",
                "neither",
            ]),
        ),
        (
            "pair-pair notboth",
            sv(&[
                "pair-pair",
                "-a",
                &fx.pe_a,
                "-b",
                &fx.pe_b,
                "--type",
                "notboth",
            ]),
        ),
        (
            "pair-pair neither",
            sv(&[
                "pair-pair",
                "-a",
                &fx.pe_a,
                "-b",
                &fx.pe_b,
                "--type",
                "neither",
            ]),
        ),
    ];
    for (name, args) in more {
        p3b_assert_broken_pipe(name, &args, 100_000);
        ran += 1;
    }
    assert!(ran >= 25, "ran {ran} broken-pipe cases");
}

/// Empty inputs: nothing to report is success, with the header or table a downstream reader
/// expects; asking for more records than there are is an error that writes nothing.
#[test]
fn p3b_empty_inputs_matrix() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t1000\nchr2\t1000\n");
    let a = fx.plain(
        "a.bed",
        "chr1\t10\t50\ta1\t100\t+\nchr1\t60\t100\ta2\t200\t-\n",
    );
    let pe = fx.plain("a.bedpe", "chr1\t10\t50\tchr1\t60\t100\tp\t1\t+\t-\n");
    let empty = fx.plain("empty.bed", "");
    let empty_gz = fx.bgzf("empty.bed.gz", "");
    let empty_pe = fx.plain("empty.bedpe", "");
    let empty_pe_gz = fx.bgzf("empty.bedpe.gz", "");
    let hdr = p2_vcf(1..=1, 0);
    let hv = fx.plain("h.vcf", &hdr);
    let hv_gz = fx.bgzf("h.vcf.gz", &hdr);
    let hv_bcf = fx.bcf("h.bcf", &hdr);
    let hv2 = fx.plain("h2.vcf", &p3b_rename_samples(&hdr));
    let one = fx.plain("one.vcf", &p2_vcf(1..=1, 3));
    let sam = fx.plain(
        "h.sam",
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n",
    );
    let fasta = p3b_fasta(fx.dir.path());
    let g = genome.as_str();
    let jaccard_empty = "intersection\tunion\tjaccard\tn_intersections\n0\t0\tnan\t0\n";
    let jaccard_a_only = "intersection\tunion\tjaccard\tn_intersections\n0\t80\t0\t0\n";
    let reldist_head = "reldist\tcount\ttotal\tfraction\n";

    for (kind, e, epe) in [
        ("plain", &empty, &empty_pe),
        ("bgzf", &empty_gz, &empty_pe_gz),
    ] {
        for t in ["1", "8"] {
            let run = |args: &[&str]| stdout_with(args, t);
            // statistics
            for (x, y, want) in [
                (e, e, jaccard_empty),
                (&a, e, jaccard_a_only),
                (e, &a, jaccard_a_only),
            ] {
                assert_eq!(
                    String::from_utf8(run(&["stats", "jaccard", "-a", x, "-b", y])).unwrap(),
                    want,
                    "{kind} jaccard"
                );
            }
            for (x, y) in [(e, e), (e, &a), (&a, e)] {
                let f = String::from_utf8(run(&["stats", "fisher", "-g", g, "-a", x, "-b", y]))
                    .unwrap();
                assert!(f.contains("left\tright\ttwo-tail\tratio\n1\t1\t1\t"), "{f}");
                assert!(f.contains("# Number of overlaps: 0\n"), "{f}");
                let r = String::from_utf8(run(&["stats", "reldist", "-a", x, "-b", y])).unwrap();
                assert_eq!(r, reldist_head, "{kind} reldist");
                assert!(
                    run(&["stats", "reldist", "-a", x, "-b", y, "--detail"]).is_empty(),
                    "{kind} reldist --detail"
                );
                assert_eq!(
                    String::from_utf8(run(&["stats", "reldist", "-a", x, "-b", y, "--sort"]))
                        .unwrap(),
                    reldist_head
                );
            }
            // randomness: nothing in, nothing out; nothing asked for, nothing out
            assert!(run(&["random", "-g", g, "-n", "0", "--seed", "1"]).is_empty());
            assert!(run(&["sample", "-n", "0", e, "--seed", "1"]).is_empty());
            assert!(run(&["shuffle", "-g", g, e, "--seed", "1"]).is_empty());
            assert!(run(&["shuffle", "-g", g, "--format", "bedpe", epe, "--seed", "1"]).is_empty());
            // pairs
            for ty in [
                "either", "both", "xor", "neither", "notboth", "ispan", "notospan",
            ] {
                assert!(
                    run(&["pair-intersect", "-a", epe, "-b", &a, "--type", ty]).is_empty(),
                    "{kind} pair-intersect {ty}: empty A"
                );
            }
            for ty in ["both", "either", "neither", "notboth"] {
                assert!(
                    run(&["pair-pair", "-a", epe, "-b", &pe, "--type", ty]).is_empty(),
                    "{kind} pair-pair {ty}: empty A"
                );
            }
        }
    }
    // an empty B: pairs with nothing to hit are the "neither" rows, and the rest write nothing
    let pe_row = "chr1\t10\t50\tchr1\t60\t100\tp\t1\t+\t-\n";
    for e in [&empty, &empty_gz] {
        assert_eq!(
            String::from_utf8(run_stdout(&[
                "pair-intersect",
                "-a",
                &pe,
                "-b",
                e,
                "--type",
                "neither"
            ]))
            .unwrap(),
            pe_row
        );
        assert!(run_stdout(&["pair-intersect", "-a", &pe, "-b", e]).is_empty());
        assert!(run_stdout(&["pair-intersect", "-a", &pe, "-b", e, "--type", "both"]).is_empty());
    }
    for e in [&empty_pe, &empty_pe_gz] {
        assert_eq!(
            String::from_utf8(run_stdout(&[
                "pair-pair",
                "-a",
                &pe,
                "-b",
                e,
                "--type",
                "neither"
            ]))
            .unwrap(),
            pe_row
        );
        assert!(run_stdout(&["pair-pair", "-a", &pe, "-b", e]).is_empty());
    }
    // -o with empty input: the file exists and is empty (or holds only the table)
    for (args, ext, len) in [
        (
            vec!["pair-intersect", "-a", &empty_pe, "-b", &a],
            "bedpe",
            0,
        ),
        (vec!["pair-pair", "-a", &empty_pe, "-b", &pe], "bedpe", 0),
        (vec!["sample", "-n", "0", &empty, "--seed", "1"], "bed", 0),
        (vec!["shuffle", "-g", g, &empty, "--seed", "1"], "bed", 0),
        (vec!["random", "-g", g, "-n", "0", "--seed", "1"], "bed", 0),
        (
            vec!["stats", "jaccard", "-a", &empty, "-b", &empty],
            "tsv",
            jaccard_empty.len() as u64,
        ),
    ] {
        let d = tempdir();
        let out = s_of(&d.path().join(format!("o.{ext}")));
        let o = gloci().args(&args).args(["-o", &out]).output().unwrap();
        assert_eq!(
            o.status.code(),
            Some(0),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert_eq!(std::fs::metadata(&out).unwrap().len(), len, "{args:?}");
    }
    // sampling from an input with no records: no records out, one warning, exit 0 (M1)
    for (name, input) in [("bed", &empty), ("bed.gz", &empty_gz), ("vcf", &hv)] {
        let d = tempdir();
        let out = s_of(&d.path().join(if name == "vcf" { "o.vcf" } else { "o.bed" }));
        for t in ["1", "4"] {
            let o = gloci()
                .args([
                    "sample", "-n", "1", input, "--seed", "1", "-t", t, "-o", &out,
                ])
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(0), "sample {name}: {err}");
            assert!(
                err.contains("is empty: it has no records") && err.lines().count() == 1,
                "sample {name}: {err}"
            );
            let text = std::fs::read_to_string(&out).unwrap();
            assert!(
                text.lines().all(|l| l.starts_with('#')),
                "sample {name}: {text}"
            );
        }
    }
    // a record-less alignment file gives no pairs
    for t in ["1", "8"] {
        assert!(stdout_with(&["convert", "--bedpe", &sam], t).is_empty());
    }

    // header-only VCF in every flavour: the commands keep the header and write no records
    let want_header: Vec<String> = hdr.lines().map(String::from).collect();
    for v in [&hv, &hv_gz, &hv_bcf] {
        for t in ["1", "8"] {
            let run = |args: &[&str]| stdout_with(args, t);
            let keeps: Vec<(&str, Vec<&str>)> = vec![
                ("split", vec!["vcf", "split", v, "-O", "vcf"]),
                (
                    "split --types",
                    vec!["vcf", "split", v, "--types", "snps", "-O", "vcf"],
                ),
                ("norm", vec!["vcf", "norm", v, "-f", &fasta, "-O", "vcf"]),
                (
                    "norm --split",
                    vec!["vcf", "norm", v, "-f", &fasta, "--split", "-O", "vcf"],
                ),
                (
                    "annotate",
                    vec!["vcf", "annotate", v, "-a", &one, "-c", "ID", "-O", "vcf"],
                ),
                (
                    "annotate (empty source)",
                    vec!["vcf", "annotate", &one, "-a", v, "-c", "ID", "-O", "vcf"],
                ),
                (
                    "filter",
                    vec!["vcf", "filter", v, "--min-qual", "3", "-O", "vcf"],
                ),
                (
                    "sample --header",
                    vec![
                        "sample", "-n", "0", "--header", v, "--seed", "1", "-O", "vcf",
                    ],
                ),
            ];
            for (name, args) in &keeps {
                if name.starts_with("sample") && v.ends_with(".bcf") {
                    continue; // sample writes records in the input's own format
                }
                let out = run(args);
                let got = p3b_header_lines(&out);
                let records = p3b_text(&out)
                    .lines()
                    .filter(|l| !l.starts_with('#'))
                    .count();
                assert_eq!(
                    records,
                    if name.starts_with("annotate (") { 3 } else { 0 },
                    "{name} {v}"
                );
                let want: &[String] = if name.starts_with("annotate (") {
                    &p2_vcf(1..=1, 0)
                        .lines()
                        .map(String::from)
                        .collect::<Vec<_>>()
                } else {
                    &want_header
                };
                // the output may add header lines (annotate, norm, split) but keeps these
                for l in want {
                    assert!(got.contains(l), "{name} -t {t}: lost `{l}`");
                }
            }
            // merge of two header-only files: both samples, no records
            let out = run(&["vcf", "merge", v, &hv2, "-O", "vcf"]);
            let head = p3b_header_lines(&out);
            assert!(
                head.last()
                    .unwrap()
                    .ends_with("FORMAT\ts1\ts2\ts3\tt1\tt2\tt3"),
                "{head:?}"
            );
            assert_eq!(
                p3b_text(&out)
                    .lines()
                    .filter(|l| !l.starts_with('#'))
                    .count(),
                0
            );
        }
        // extract writes its header row and no data row
        assert_eq!(
            String::from_utf8(run_stdout(&["vcf", "extract", v, "--info", "DP"])).unwrap(),
            "CHROM\tPOS\tID\tREF\tALT\tDP\n"
        );
        // the files of stats, popgen and compare hold their headers, no rows
        let d = tempdir();
        let p = s_of(&d.path().join("s"));
        gloci()
            .args([
                "vcf",
                "stats",
                v,
                "--out-prefix",
                &p,
                "--freq",
                "--site-quality",
                "--het",
                "--missing-indv",
            ])
            .assert()
            .success();
        for (ext, rows) in [("frq", 1), ("lqual", 1), ("het", 1), ("imiss", 4)] {
            let t = std::fs::read_to_string(format!("{p}.{ext}")).unwrap();
            assert_eq!(t.lines().count(), rows, ".{ext}: {t}");
        }
        gloci()
            .args(["vcf", "popgen", v, "--out-prefix", &p, "--site-pi"])
            .args(["--window-pi", "10", "--tajima-d", "10"])
            .assert()
            .success();
        for ext in ["sites.pi", "windowed.pi", "Tajima.D"] {
            let t = std::fs::read_to_string(format!("{p}.{ext}")).unwrap();
            assert_eq!(t.lines().count(), 1, ".{ext}: {t}");
        }
        for (x, y) in [(v, v), (v, &one), (&one, v)] {
            let d = tempdir();
            let p = s_of(&d.path().join("c"));
            for t in ["1", "8"] {
                gloci()
                    .args(["vcf", "compare", x, y, "--out-prefix", &p, "-t", t])
                    .assert()
                    .success();
                let sites = std::fs::read_to_string(format!("{p}.diff.sites")).unwrap();
                assert!(sites.starts_with("CHROM\tPOS\tFILES\t"), "{sites}");
                // records on one side only are not discordant: no site has two files
                assert!(
                    sites
                        .lines()
                        .skip(1)
                        .all(|l| l.split('\t').nth(2) != Some("B") || x != y),
                    "{sites}"
                );
                for ext in ["diff.indv", "diff.discordance_matrix"] {
                    assert!(
                        std::path::Path::new(&format!("{p}.{ext}")).exists(),
                        "{ext}"
                    );
                }
            }
        }
    }
    // an empty source annotates nothing and an empty VCF annotates to the header
    let out = run_stdout(&["vcf", "annotate", &one, "-a", &hv, "-c", "ID"]);
    assert_eq!(
        p3b_text(&out)
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>(),
        p2_vcf(1..=1, 3)
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
    );
}

/// The failure matrix of the P3b commands: unsorted input where order matters, malformed
/// lines, truncated BGZF/BCF, a missing output directory; every run exits 1 with the
/// `error:` format and leaves nothing behind.
#[test]
fn p3b_failure_matrix() {
    let fx = P2::new();
    let genome = fx.plain(
        "g.genome",
        "chr1\t1000000\nchr2\t1000000\nchr3\t1000000\nchr4\t1000000\n",
    );
    let good_bed = fx.plain("good.bed", &p3b_bed(2, 40, 100, 50, 30));
    let other_bed = fx.plain("other.bed", &p3b_bed(2, 30, 130, 80, 40));
    let pe_text = p3b_bedpe(40, 0);
    let good_pe = fx.plain("good.bedpe", &pe_text);
    let other_pe = fx.plain("other.bedpe", &p3b_bedpe(20, 0));
    let g = genome.as_str();

    let bed_text = p3b_bed(2, 40, 100, 50, 30);
    let bed_lines: Vec<&str> = bed_text.lines().collect();
    let unsorted_bed = fx.plain(
        "u.bed",
        &format!("{}\n{}\n{}\n", bed_lines[3], bed_lines[2], bed_lines[4]),
    );
    let malformed_bed = {
        let mut l: Vec<String> = bed_lines.iter().map(|s| s.to_string()).collect();
        let f: Vec<&str> = bed_lines[5].split('\t').collect();
        l[5] = format!("{}\tNOPE\t{}", f[0], f[2]);
        fx.plain("m.bed", &(l.join("\n") + "\n"))
    };
    let bed_gz = bgzip(bed_text.as_bytes());
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

    let pe_lines: Vec<&str> = pe_text.lines().collect();
    let malformed_pe = {
        let mut l: Vec<String> = pe_lines.iter().map(|s| s.to_string()).collect();
        let mut f: Vec<&str> = pe_lines[5].split('\t').collect();
        f[4] = "NOPE";
        l[5] = f.join("\t");
        fx.plain("m.bedpe", &(l.join("\n") + "\n"))
    };
    let pe_gz = bgzip(pe_text.as_bytes());
    let trunc_pe = s_of(&write(
        fx.dir.path(),
        "t.bedpe.gz",
        &pe_gz[..pe_gz.len() / 2],
    ));
    let nofooter_pe = s_of(&write(
        fx.dir.path(),
        "nf.bedpe.gz",
        &pe_gz[..pe_gz.len() - 28],
    ));

    let vs = p3b_bad_vcf(&fx, false);
    let vr = p3b_bad_vcf(&fx, true);
    let fasta = p3b_fasta(fx.dir.path());
    let eof = "BGZF end-of-file marker missing";
    let unsorted_msg = "is not sorted";
    let not_int = "is not an integer";

    // BED-shaped inputs: (kind, path, needs sorted input, message)
    let bed_bad: [(&str, &String, bool, &str); 4] = [
        ("unsorted", &unsorted_bed, true, unsorted_msg),
        ("malformed", &malformed_bed, false, not_int),
        ("truncated", &trunc_bed, false, eof),
        ("no BGZF footer", &nofooter_bed, false, eof),
    ];
    let pe_bad: [(&str, &String, &str); 3] = [
        ("malformed", &malformed_pe, not_int),
        ("truncated", &trunc_pe, eof),
        ("no BGZF footer", &nofooter_pe, eof),
    ];
    let (gb, ob, gp, op) = (&good_bed, &other_bed, &good_pe, &other_pe);

    // (name, template with {IN}, sink, needs sorted input)
    type Case<'a> = (&'a str, Vec<&'a str>, Sink, bool);
    let bed_cases: Vec<Case> = vec![
        (
            "jaccard A",
            vec!["stats", "jaccard", "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
            true,
        ),
        (
            "jaccard B",
            vec!["stats", "jaccard", "-a", gb, "-b", "{IN}"],
            Sink::File("tsv"),
            true,
        ),
        (
            "fisher A",
            vec!["stats", "fisher", "-g", g, "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
            true,
        ),
        (
            "fisher B",
            vec!["stats", "fisher", "-g", g, "-a", gb, "-b", "{IN}"],
            Sink::File("tsv"),
            true,
        ),
        (
            "fisher --exclude",
            vec![
                "stats",
                "fisher",
                "-g",
                g,
                "-a",
                gb,
                "-b",
                ob,
                "--exclude",
                "{IN}",
            ],
            Sink::File("tsv"),
            true,
        ),
        (
            "reldist A",
            vec!["stats", "reldist", "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
            true,
        ),
        (
            "reldist B",
            vec!["stats", "reldist", "-a", gb, "-b", "{IN}"],
            Sink::File("tsv"),
            true,
        ),
        (
            "reldist --detail",
            vec!["stats", "reldist", "--detail", "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
            true,
        ),
        (
            "sample",
            vec!["sample", "-n", "3", "--seed", "1", "{IN}"],
            Sink::File("bed"),
            false,
        ),
        (
            "shuffle",
            vec!["shuffle", "-g", g, "--seed", "1", "{IN}"],
            Sink::File("bed"),
            false,
        ),
        (
            "shuffle --exclude",
            vec!["shuffle", "-g", g, "--seed", "1", "--exclude", "{IN}", gb],
            Sink::File("bed"),
            false,
        ),
        (
            "pair-intersect B",
            vec!["pair-intersect", "-a", gp, "-b", "{IN}"],
            Sink::File("bedpe"),
            false,
        ),
    ];
    let mut ran = 0;
    for (name, tpl, sink, sorted) in &bed_cases {
        for (kind, bad, needs_sorted, msg) in &bed_bad {
            if *needs_sorted && !*sorted {
                continue;
            }
            // shuffle/sample/pair B take unsorted input without complaint
            if *kind == "unsorted" && !*sorted {
                continue;
            }
            assert_fails_clean(name, kind, tpl, bad, sink, msg);
            ran += 1;
        }
    }
    let pe_cases: Vec<Case> = vec![
        (
            "pair-intersect A",
            vec!["pair-intersect", "-a", "{IN}", "-b", ob],
            Sink::File("bedpe"),
            false,
        ),
        (
            "pair-pair A",
            vec!["pair-pair", "-a", "{IN}", "-b", op],
            Sink::File("bedpe"),
            false,
        ),
        (
            "pair-pair B",
            vec!["pair-pair", "-a", gp, "-b", "{IN}"],
            Sink::File("bedpe"),
            false,
        ),
        (
            "shuffle bedpe",
            vec![
                "shuffle", "-g", g, "--seed", "1", "--format", "bedpe", "{IN}",
            ],
            Sink::File("bedpe"),
            false,
        ),
    ];
    for (name, tpl, sink, _) in &pe_cases {
        for (kind, bad, msg) in &pe_bad {
            assert_fails_clean(name, kind, tpl, bad, sink, msg);
            ran += 1;
        }
    }

    // VCF inputs: (kind, bad input, bad input with the other sample names, needs sorted, message)
    let vcf_bad: [(&str, &String, &String, bool, &str); 5] = [
        ("unsorted", &vs.unsorted, &vr.unsorted, true, unsorted_msg),
        ("malformed", &vs.malformed, &vr.malformed, false, not_int),
        ("truncated", &vs.trunc_gz, &vr.trunc_gz, false, eof),
        (
            "no BGZF footer",
            &vs.nofooter_gz,
            &vr.nofooter_gz,
            false,
            eof,
        ),
        ("truncated BCF", &vs.trunc_bcf, &vr.trunc_bcf, false, ""),
    ];
    let (sg, rg) = (&vs.good, &vr.good);
    // (name, template, sink, takes the sample-renamed input, order matters)
    type VCase<'a> = (&'a str, Vec<&'a str>, Sink, bool, bool);
    let vcf_cases: Vec<VCase> = vec![
        (
            "vcf split",
            vec!["vcf", "split", "{IN}"],
            Sink::File("vcf"),
            false,
            false,
        ),
        (
            "vcf norm",
            vec!["vcf", "norm", "-f", &fasta, "-c", "s", "{IN}"],
            Sink::File("vcf"),
            false,
            false,
        ),
        (
            "vcf annotate",
            vec!["vcf", "annotate", "-a", sg, "-c", "ID", "{IN}"],
            Sink::File("vcf"),
            false,
            true,
        ),
        (
            "vcf annotate source",
            vec!["vcf", "annotate", "-a", "{IN}", "-c", "ID", sg],
            Sink::File("vcf"),
            false,
            true,
        ),
        (
            "vcf merge first",
            vec!["vcf", "merge", "{IN}", rg],
            Sink::File("vcf"),
            false,
            true,
        ),
        (
            "vcf merge second",
            vec!["vcf", "merge", sg, "{IN}"],
            Sink::File("vcf"),
            true,
            true,
        ),
        (
            "vcf compare A",
            vec!["vcf", "compare", "{IN}", sg],
            Sink::Prefix,
            false,
            true,
        ),
        (
            "vcf compare B",
            vec!["vcf", "compare", sg, "{IN}"],
            Sink::Prefix,
            false,
            true,
        ),
        (
            "vcf filter",
            vec!["vcf", "filter", "--min-qual", "1", "{IN}"],
            Sink::File("vcf"),
            false,
            false,
        ),
        (
            "vcf extract",
            vec!["vcf", "extract", "--info", "DP", "{IN}"],
            Sink::File("tsv"),
            false,
            false,
        ),
        (
            "vcf stats",
            vec!["vcf", "stats", "--freq", "{IN}"],
            Sink::Prefix,
            false,
            false,
        ),
        (
            "vcf popgen",
            vec!["vcf", "popgen", "--site-pi", "{IN}"],
            Sink::Prefix,
            false,
            false,
        ),
        (
            "sample vcf",
            vec!["sample", "-n", "3", "--seed", "1", "{IN}"],
            Sink::File("vcf"),
            false,
            false,
        ),
    ];
    for (name, tpl, sink, renamed, ordered) in &vcf_cases {
        for (kind, bad_s, bad_r, needs_sorted, msg) in &vcf_bad {
            if *needs_sorted && !*ordered {
                continue;
            }
            let bad = if *renamed { bad_r } else { bad_s };
            let msg = if msg.is_empty() { "" } else { msg };
            assert_fails_clean(name, kind, tpl, bad, sink, msg);
            ran += 1;
        }
    }
    // a VCF where the reference disagrees
    assert_fails_clean(
        "vcf norm",
        "REF mismatch",
        &["vcf", "norm", "-f", &fasta, "{IN}"],
        &vs.good,
        &Sink::File("vcf"),
        "REF does not match the reference",
    );
    ran += 1;
    // alignments to pairs: a bad SAM line, a truncated BAM
    let bad_sam = fx.plain(
        "bad.sam",
        "@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\nr2\tNOPE\tchr1\t10\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    );
    assert_fails_clean(
        "convert --bedpe",
        "malformed",
        &["convert", "--bedpe", "{IN}"],
        &bad_sam,
        &Sink::File("bedpe"),
        "is not a number",
    );
    let sam_text = p3b_pairs_sam(200);
    let big_sam = fx.plain("pairs.sam", &sam_text);
    let big_bam = s_of(&fx.dir.path().join("pairs.bam"));
    gloci()
        .args(["view", &big_sam, "-o", &big_bam])
        .assert()
        .success();
    let bam = std::fs::read(&big_bam).unwrap();
    let trunc_bam = s_of(&write(fx.dir.path(), "t.bam", &bam[..bam.len() / 2]));
    assert_fails_clean(
        "convert --bedpe",
        "truncated BAM",
        &["convert", "--bedpe", "{IN}"],
        &trunc_bam,
        &Sink::File("bedpe"),
        "",
    );
    ran += 2;

    // a missing output directory: every P3b command that writes a file
    let all: Vec<(&str, Vec<&str>, Sink)> = vec![
        (
            "jaccard",
            vec!["stats", "jaccard", "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
        ),
        (
            "fisher",
            vec!["stats", "fisher", "-g", g, "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
        ),
        (
            "reldist",
            vec!["stats", "reldist", "-a", "{IN}", "-b", ob],
            Sink::File("tsv"),
        ),
        (
            "random",
            vec!["random", "-g", g, "-n", "5", "--seed", "1"],
            Sink::File("bed"),
        ),
        (
            "sample",
            vec!["sample", "-n", "3", "--seed", "1", "{IN}"],
            Sink::File("bed"),
        ),
        (
            "shuffle",
            vec!["shuffle", "-g", g, "--seed", "1", "{IN}"],
            Sink::File("bed"),
        ),
        (
            "pair-intersect",
            vec!["pair-intersect", "-a", gp, "-b", "{IN}"],
            Sink::File("bedpe"),
        ),
        (
            "pair-pair",
            vec!["pair-pair", "-a", "{IN}", "-b", op],
            Sink::File("bedpe"),
        ),
        (
            "convert --bedpe",
            vec!["convert", "--bedpe", &big_sam],
            Sink::File("bedpe"),
        ),
        ("vcf split", vec!["vcf", "split", "{IN}"], Sink::File("vcf")),
        (
            "vcf norm",
            vec!["vcf", "norm", "-f", &fasta, "-c", "s", "{IN}"],
            Sink::File("vcf"),
        ),
        (
            "vcf annotate",
            vec!["vcf", "annotate", "-a", sg, "-c", "ID", "{IN}"],
            Sink::File("vcf"),
        ),
        (
            "vcf merge",
            vec!["vcf", "merge", "{IN}", rg],
            Sink::File("vcf"),
        ),
        (
            "vcf compare",
            vec!["vcf", "compare", "{IN}", sg],
            Sink::Prefix,
        ),
        (
            "vcf filter",
            vec!["vcf", "filter", "--min-qual", "1", "{IN}"],
            Sink::File("vcf"),
        ),
        (
            "vcf extract",
            vec!["vcf", "extract", "--info", "DP", "{IN}"],
            Sink::File("tsv"),
        ),
        (
            "vcf stats",
            vec!["vcf", "stats", "--freq", "{IN}"],
            Sink::Prefix,
        ),
        (
            "vcf popgen",
            vec!["vcf", "popgen", "--site-pi", "{IN}"],
            Sink::Prefix,
        ),
    ];
    for (name, tpl, sink) in &all {
        let good = match name.split(' ').next().unwrap() {
            "vcf" => sg.as_str(),
            "pair-intersect" => ob.as_str(),
            "pair-pair" => gp.as_str(),
            "jaccard" | "fisher" | "reldist" | "sample" | "shuffle" => gb.as_str(),
            _ => "",
        };
        assert_missing_dir(name, tpl, good, sink);
        ran += 1;
    }
    assert!(ran > 120, "ran {ran} failure cases");
}

/// `--chrom-map` renames as records are read: a run with the map on files in the source names
/// equals the run without it on files in the target names, and the inverse map takes the output
/// of the commands that write records back to the source names. (Seeded runs use one seed.)
#[test]
fn p3b_chrom_map_round_trip() {
    let fx = P2::new();
    let fwd = fx.plain("fwd.map", "1\tchr1\n2\tchr2\n3\tchr3\n4\tchr4\n");
    let inv = fx.plain("inv.map", "chr1\t1\nchr2\t2\nchr3\t3\nchr4\t4\n");
    let genome = fx.plain(
        "g.genome",
        "chr1\t1000000\nchr2\t1000000\nchr3\t1000000\nchr4\t1000000\n",
    );
    let strip = |t: &str| t.replace("chr", "");
    let a_txt = p3b_bed(4, 200, 300, 100, 60);
    let b_txt = p3b_bed(4, 150, 410, 130, 90);
    let pe_txt = p3b_bedpe(300, 0);
    let (a_chr, a_num) = (
        fx.plain("a_chr.bed", &a_txt),
        fx.plain("a_num.bed", &strip(&a_txt)),
    );
    let (b_chr, b_num) = (
        fx.plain("b_chr.bed", &b_txt),
        fx.plain("b_num.bed", &strip(&b_txt)),
    );
    let (pe_chr, pe_num) = (
        fx.plain("pe_chr.bedpe", &pe_txt),
        fx.plain("pe_num.bedpe", &strip(&pe_txt)),
    );
    let pe2_txt = p3b_bedpe(100, 0);
    let (pe2_chr, pe2_num) = (
        fx.plain("pe2_chr.bedpe", &pe2_txt),
        fx.plain("pe2_num.bedpe", &strip(&pe2_txt)),
    );
    let vcf_chr_txt = p3b_vcf(1..=4, 60);
    let vcf2_chr_txt = p3b_vcf(2..=4, 60);
    let v_chr = fx.plain("v_chr.vcf", &vcf_chr_txt);
    let v_num = fx.plain("v_num.vcf", &strip(&vcf_chr_txt));
    let v2_chr = fx.plain("v2_chr.vcf", &vcf2_chr_txt);
    let v2_num = fx.plain("v2_num.vcf", &strip(&vcf2_chr_txt));
    let v_num_gz = fx.bgzf("v_num.vcf.gz", &strip(&vcf_chr_txt));
    let g = genome.as_str();

    // (name, template with {A} {B} {PA} {PB}, inputs are chr or numeric)
    let cmds: Vec<(&str, Vec<&str>)> = vec![
        (
            "jaccard",
            vec!["stats", "jaccard", "-a", "{A}", "-b", "{B}"],
        ),
        (
            "fisher",
            vec!["stats", "fisher", "-g", g, "-a", "{A}", "-b", "{B}"],
        ),
        (
            "reldist",
            vec!["stats", "reldist", "-a", "{A}", "-b", "{B}", "--detail"],
        ),
        ("sample", vec!["sample", "-n", "40", "--seed", "3", "{A}"]),
        (
            "shuffle",
            vec!["shuffle", "-g", g, "--seed", "3", "--exclude", "{B}", "{A}"],
        ),
        (
            "shuffle same-chrom",
            vec!["shuffle", "-g", g, "--seed", "3", "--same-chrom", "{A}"],
        ),
        (
            "shuffle bedpe",
            vec![
                "shuffle", "-g", g, "--seed", "3", "--format", "bedpe", "{PA}",
            ],
        ),
        (
            "pair-intersect",
            vec!["pair-intersect", "-a", "{PA}", "-b", "{B}"],
        ),
        (
            "pair-intersect neither",
            vec![
                "pair-intersect",
                "-a",
                "{PA}",
                "-b",
                "{B}",
                "--type",
                "neither",
            ],
        ),
        ("pair-pair", vec!["pair-pair", "-a", "{PA}", "-b", "{PB}"]),
        (
            "pair-pair neither",
            vec!["pair-pair", "-a", "{PA}", "-b", "{PB}", "--type", "neither"],
        ),
    ];
    let fill = |tpl: &[&str], num: bool| -> Vec<String> {
        tpl.iter()
            .map(|x| match *x {
                "{A}" => (if num { &a_num } else { &a_chr }).clone(),
                "{B}" => (if num { &b_num } else { &b_chr }).clone(),
                "{PA}" => (if num { &pe_num } else { &pe_chr }).clone(),
                "{PB}" => (if num { &pe2_num } else { &pe2_chr }).clone(),
                x => x.to_string(),
            })
            .collect()
    };
    for (name, tpl) in &cmds {
        for t in ["1", "8"] {
            let run = |args: Vec<String>| -> Vec<u8> {
                let refs: Vec<&str> = args.iter().map(String::as_str).collect();
                stdout_with(&refs, t)
            };
            let mut mapped = fill(tpl, true);
            mapped.extend(["--chrom-map".into(), fwd.clone()]);
            let with_map = run(mapped);
            let renamed = run(fill(tpl, false));
            assert!(!with_map.is_empty(), "{name}: nothing written");
            assert_eq!(with_map, renamed, "{name} -t {t}: map vs renamed input");
            // the output names the target chromosomes
            if !matches!(*name, "jaccard" | "fisher") {
                assert!(
                    !String::from_utf8_lossy(&with_map)
                        .lines()
                        .any(|l| l.starts_with("1\t")),
                    "{name}: a source name in the output"
                );
            }
        }
    }
    // the inverse map takes the records the commands write back to the source names
    for (name, tpl) in cmds.iter().filter(|(n, _)| {
        matches!(
            *n,
            "sample" | "pair-intersect neither" | "pair-pair neither"
        )
    }) {
        let mut back = fill(tpl, false);
        back.extend(["--chrom-map".into(), inv.clone()]);
        let refs: Vec<&str> = back.iter().map(String::as_str).collect();
        let out = stdout_with(&refs, "2");
        let plain = fill(tpl, true);
        let refs: Vec<&str> = plain.iter().map(String::as_str).collect();
        assert_eq!(
            out,
            stdout_with(&refs, "2"),
            "{name}: the inverse map restores the names"
        );
    }

    // vcf compare: the files carry the target names, whichever side the map is for
    let files = |a: &str, b: &str, map: bool, t: &str| -> Vec<Vec<u8>> {
        let mut args = sv(&[
            "vcf",
            "compare",
            a,
            b,
            "--site-discordance",
            "--sample-discordance",
            "--discordance-matrix",
        ]);
        if map {
            args.extend(["--chrom-map".into(), fwd.clone()]);
        }
        p3b_prefix_files(
            &args,
            &["diff.sites", "diff.indv", "diff.discordance_matrix"],
            t,
        )
    };
    let want = files(&v_chr, &v2_chr, false, "1");
    assert!(want[0].windows(5).any(|w| w == b"\nchr2"));
    for input in [&v_num, &v_num_gz] {
        for t in ["1", "8"] {
            assert_eq!(
                files(input, &v2_num, true, t),
                want,
                "compare {input} -t {t}"
            );
        }
    }
    // a map file with a bad line says where
    let bad = fx.plain("bad.map", "1\tchr1\nlonely\n");
    for args in [
        sv(&["stats", "jaccard", "-a", &a_num, "-b", &b_num]),
        sv(&["sample", "-n", "3", "--seed", "1", &a_num]),
        sv(&["pair-intersect", "-a", &pe_num, "-b", &b_num]),
        sv(&["pair-pair", "-a", &pe_num, "-b", &pe2_num]),
    ] {
        gloci()
            .args(&args)
            .args(["--chrom-map", &bad])
            .assert()
            .code(2)
            .stderr(contains("bad.map:2:"));
    }
}

/// Every P3c command gives the same bytes at `-t 1`, `-t2` and `--threads=8`, on BAM (indexed),
/// SAM and CRAM (indexed); the three formats of the same reads agree for the commands whose
/// output does not name the format.
#[test]
fn p3c_determinism_matrix() {
    let fx = P3a::new();
    let mut ran = 0;
    let mut by_cmd: std::collections::BTreeMap<String, Vec<(String, Vec<u8>)>> = Default::default();
    let mut all: Vec<(String, P3cCmd)> = Vec::new();
    for (kind, input, _) in fx.kinds() {
        for c in fx.p3c_cmds(kind, &input) {
            all.push((kind.to_string(), c));
        }
    }
    for c in fx.p3c_fasta_cmds() {
        all.push(("fasta".into(), c));
    }
    for (kind, c) in &all {
        let outs: Vec<Vec<u8>> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = p3c_run(c, t);
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{} -t {t}: {}",
                    c.name,
                    String::from_utf8_lossy(&o.stderr)
                );
                o.stdout
            })
            .collect();
        assert!(
            outs.windows(2).all(|w| w[0] == w[1]),
            "{}: -t changes bytes",
            c.name
        );
        assert!(!outs[0].is_empty(), "{}: wrote nothing", c.name);
        ran += 1;
        if c.format_free {
            let key = c
                .name
                .strip_prefix(&format!("{kind} "))
                .unwrap()
                .to_string();
            by_cmd
                .entry(key)
                .or_default()
                .push((kind.clone(), outs[0].clone()));
        }
    }
    for (name, outs) in &by_cmd {
        for w in outs.windows(2) {
            assert!(w[0].1 == w[1].1, "{name}: {} and {} differ", w[0].0, w[1].0);
        }
    }
    assert!(
        by_cmd.len() >= 14 && ran > 60,
        "{} formats-free, {ran} runs",
        by_cmd.len()
    );

    // reheader: the result of each thread count reads back the same, for BAM and CRAM.
    let hdr = {
        let mut h = p3c_sam_of(&fx.bam, &[])
            .split(|&b| b == b'\n')
            .filter(|l| !l.is_empty())
            .take_while(|l| l.starts_with(b"@"))
            .flat_map(|l| l.iter().chain(b"\n".iter()).copied())
            .collect::<Vec<u8>>();
        h.extend_from_slice(b"@CO\tp3c determinism\n");
        s_of(&write(Path::new(&fx.path("")), "new-header.sam", &h))
    };
    for (kind, input, extra) in fx.kinds() {
        if kind == "sam" {
            continue;
        }
        let extra: Vec<&str> = extra.iter().map(|s| s.as_str()).collect();
        let mut seen: Vec<Vec<u8>> = Vec::new();
        for t in ["1", "2", "8"] {
            let out = fx.path(&format!("reheader-{kind}-{t}.out"));
            let o = p3a_run_t(
                &[
                    "reheader".into(),
                    hdr.clone(),
                    input.clone(),
                    "-o".into(),
                    out.clone(),
                ],
                t,
            );
            assert_eq!(
                o.status.code(),
                Some(0),
                "reheader {kind} -t {t}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            let text = p3c_sam_of(&out, &extra);
            assert!(
                String::from_utf8_lossy(&text).contains("@CO\tp3c determinism\n"),
                "reheader {kind}"
            );
            seen.push(text);
        }
        assert!(seen.windows(2).all(|w| w[0] == w[1]), "reheader {kind}");
        // and the reads are the input's
        assert_eq!(
            p3a_records(&seen[0]),
            p3a_records(&p3c_sam_of(&input, &extra)),
            "reheader {kind} changed the reads"
        );
    }
}

/// depth and bedcov split their work across references and batches of regions at `-t` above 1:
/// on an indexed BAM with enough reads and regions for several batches, every thread count gives
/// the bytes of `-t 1`, for each filter that changes the walk.
#[test]
fn p3c_determinism_large_fixture() {
    let d = tempdir();
    let (big_sam, big_bam) = p3a_big_alignments(d.path(), 12_000, false);
    let bed = p3c_big_bed(d.path(), 2_000);
    let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
    let cmds = vec![
        ("depth", v(&["depth", &big_bam])),
        ("depth sam", v(&["depth", &big_sam])),
        ("depth -s", v(&["depth", &big_bam, "--skip-mate-overlap"])),
        (
            "depth regions-file",
            v(&["depth", &big_bam, "--regions-file", &bed, "--all-positions"]),
        ),
        (
            "depth two inputs",
            v(&["depth", &big_bam, &big_bam, "--min-mapq", "10"]),
        ),
        ("bedcov", v(&["bedcov", &bed, &big_bam])),
        (
            "bedcov all columns",
            v(&[
                "bedcov",
                &bed,
                &big_bam,
                &big_bam,
                "--read-count",
                "--bases-at-depth",
                "3",
            ]),
        ),
        ("idxstats", v(&["idxstats", &big_bam])),
        ("idxstats scan", v(&["idxstats", &big_sam])),
        ("flagstat", v(&["flagstat", &big_bam])),
        ("view --count", v(&["view", "--count", &big_bam])),
    ];
    let mut last: Option<Vec<u8>> = None;
    for (name, args) in &cmds {
        let outs: Vec<Vec<u8>> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = p3a_run_t(args, t);
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{name} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                o.stdout
            })
            .collect();
        assert!(outs[0].len() >= 5, "{name}: {} bytes", outs[0].len());
        assert!(
            outs.windows(2).all(|w| w[0] == w[1]),
            "{name}: -t changes bytes"
        );
        if name.starts_with("depth") && name.len() == 5 {
            last = Some(outs[0].clone());
        }
        if *name == "depth sam" {
            assert!(Some(&outs[0]) == last.as_ref(), "SAM and BAM depth differ");
        }
    }
}

/// Closing stdout after the first bytes is a quiet success (exit 0, nothing on stderr) for every
/// P3c command whose output outgrows a pipe: `depth --all-references | head`, `faidx | head`
/// and the rest. The commands that write a few lines (`flagstat`, `view --count`, `flags` with
/// an argument, `check --quick`, `reheader` to a file) have nothing to cut short.
#[test]
fn p3c_broken_pipe_matrix() {
    let d = tempdir();
    let (big_fa, big_bed) = p3a_big_fasta(d.path());
    let (big_sam, big_bam) = p3a_big_alignments(d.path(), 12_000, false);
    let many_refs = p3c_many_references(d.path(), 120_000);
    let many_seqs = p3c_many_sequences(d.path(), 4_000);
    let cov_bed = p3c_big_bed(d.path(), 20_000);
    let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
    // Commands that write tens of millions of rows are not run to the end.
    let huge = [
        "depth --all-references",
        "depth --all-positions",
        "depth regions",
        "depth many references",
    ];
    let cmds: Vec<(&str, Vec<String>)> = vec![
        (
            "depth --all-references",
            v(&["depth", &big_bam, "--all-references"]),
        ),
        (
            "depth --all-positions",
            v(&["depth", &big_bam, "--all-positions"]),
        ),
        ("depth", v(&["depth", &big_bam])),
        ("depth sam", v(&["depth", &big_sam])),
        ("depth -s", v(&["depth", &big_bam, "--skip-mate-overlap"])),
        (
            "depth regions",
            v(&[
                "depth",
                &big_bam,
                "--region",
                "chr1:1-9000000",
                "--all-positions",
            ]),
        ),
        (
            "depth many references",
            v(&["depth", &many_refs, "--all-references"]),
        ),
        ("bedcov", v(&["bedcov", &cov_bed, &big_bam])),
        (
            "bedcov all columns",
            v(&[
                "bedcov",
                &cov_bed,
                &big_bam,
                "--read-count",
                "--bases-at-depth",
                "1",
            ]),
        ),
        ("idxstats scan", v(&["idxstats", &many_refs])),
        ("faidx whole", v(&["faidx", &big_fa, "big"])),
        (
            "faidx --reverse-complement",
            v(&["faidx", &big_fa, "big", "--reverse-complement"]),
        ),
        (
            "faidx --regions-file",
            v(&["faidx", &big_fa, "--regions-file", &big_bed]),
        ),
        (
            "faidx line width",
            v(&["faidx", &big_fa, "big", "--line-width", "10"]),
        ),
        ("dict", v(&["dict", &many_seqs])),
        ("head", v(&["head", &big_sam, "--records", "12000"])),
        ("head bam", v(&["head", &big_bam, "--records", "12000"])),
        ("head many references", v(&["head", &many_refs])),
        ("flags table x many", {
            // `flags` prints its twelve-line table or one line per argument
            let mut a = v(&["flags"]);
            a.extend((0..20_000).map(|i| (i % 4096).to_string()));
            a
        }),
        ("samples", {
            let mut a = v(&["samples"]);
            a.extend((0..4_000).map(|_| big_sam.clone()));
            a
        }),
    ];
    for (name, args) in &cmds {
        if !huge.contains(name) {
            let full = gloci().args(args).arg("-q").output().unwrap();
            assert_eq!(
                full.status.code(),
                Some(0),
                "{name}: {}",
                String::from_utf8_lossy(&full.stderr)
            );
            assert!(
                full.stdout.len() > 100_000,
                "{name}: only {} bytes of output",
                full.stdout.len()
            );
        }
        for t in ["1", "4"] {
            let mut cmd = Command::new(bin());
            cmd.args(args);
            if !matches!(args[0].as_str(), "flags" | "samples") {
                cmd.args(["-t", t]);
            }
            cmd.arg("-q");
            let mut child = cmd
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            let mut first = [0u8; 64];
            child
                .stdout
                .as_mut()
                .unwrap()
                .read_exact(&mut first)
                .unwrap();
            drop(child.stdout.take());
            let out = child.wait_with_output().unwrap();
            assert_eq!(
                out.status.code(),
                Some(0),
                "{name} -t {t}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert!(
                out.stderr.is_empty(),
                "{name} -t {t}: stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
    }
}

/// A header-only BAM and SAM, a BED with no lines, and FASTA files with no sequences: every P3c
/// command exits 0 (or says clearly why not), writes the valid empty result and prints nothing on
/// stderr, at `-t 1` and `-t 4`.
#[test]
fn p3c_empty_inputs_matrix() {
    let d = tempdir();
    let dir = d.path();
    let hdr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n";
    let sam = s_of(&plain_file(dir, "h.sam", hdr));
    let bam = s_of(&dir.join("h.bam"));
    gloci().args(["view", &sam, "-o", &bam]).assert().success();
    gloci().args(["index", &bam]).assert().success();
    let bed = s_of(&plain_file(dir, "t.bed", "chr1\t0\t100\nchr2\t10\t20\n"));
    let empty_bed = s_of(&plain_file(dir, "e.bed", ""));
    let zero_flagstat = [
        "0 + 0 in total (QC-passed reads + QC-failed reads)",
        "0 + 0 primary",
        "0 + 0 secondary",
        "0 + 0 supplementary",
        "0 + 0 duplicates",
        "0 + 0 primary duplicates",
        "0 + 0 mapped (N/A : N/A)",
        "0 + 0 primary mapped (N/A : N/A)",
        "0 + 0 paired in sequencing",
        "0 + 0 read1",
        "0 + 0 read2",
        "0 + 0 properly paired (N/A : N/A)",
        "0 + 0 with itself and mate mapped",
        "0 + 0 singletons (N/A : N/A)",
        "0 + 0 with mate mapped to a different chr",
        "0 + 0 with mate mapped to a different chr (mapQ>=5)",
    ]
    .join("\n")
        + "\n";
    let ran = std::cell::Cell::new(0);
    for t in ["1", "4"] {
        let run = |args: &[&str], code: i32| -> std::process::Output {
            let mut cmd = gloci();
            cmd.args(args).arg("-q");
            if args[0] != "samples" {
                cmd.args(["-t", t]);
            }
            let o = cmd.output().unwrap();
            assert_eq!(
                o.status.code(),
                Some(code),
                "{args:?} -t {t}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            if code == 0 {
                assert!(
                    o.stderr.is_empty(),
                    "{args:?} -t {t}: stderr {}",
                    String::from_utf8_lossy(&o.stderr)
                );
            }
            ran.set(ran.get() + 1);
            o
        };
        let text = |o: &std::process::Output| String::from_utf8_lossy(&o.stdout).into_owned();
        for input in [&bam, &sam] {
            assert_eq!(
                text(&run(&["flagstat", input], 0)),
                zero_flagstat,
                "{input}"
            );
            let j = text(&run(&["flagstat", input, "-O", "json"], 0));
            assert!(j.contains("\"total\": 0") && j.contains("null"), "{j}");
            let tsv = text(&run(&["flagstat", input, "-O", "tsv"], 0));
            assert!(
                tsv.starts_with("0\t0\ttotal (QC-passed reads + QC-failed reads)\n")
                    && tsv.contains("N/A\tN/A\tmapped %\n"),
                "{tsv}"
            );
            assert_eq!(
                text(&run(&["idxstats", input], 0)),
                "chr1\t1000\t0\t0\nchr2\t500\t0\t0\n*\t0\t0\t0\n",
                "{input}"
            );
            assert_eq!(text(&run(&["view", "--count", input], 0)), "0\n");
            assert_eq!(
                text(&run(&["view", "--count", input, "--min-mapq", "20"], 0)),
                "0\n"
            );
            // depth: no rows; every position of every reference as zeros with --all-references
            assert_eq!(text(&run(&["depth", input], 0)), "");
            assert_eq!(text(&run(&["depth", input, "--all-positions"], 0)), "");
            let aa = text(&run(&["depth", input, "--all-references"], 0));
            assert_eq!(aa.lines().count(), 1500, "{input}");
            assert!(aa.starts_with("chr1\t1\t0\n") && aa.ends_with("chr2\t500\t0\n"));
            assert!(aa.lines().all(|l| l.ends_with("\t0")));
            assert_eq!(
                text(&run(&["depth", input, "--column-names"], 0)),
                format!("#CHROM\tPOS\t{input}\n")
            );
            assert_eq!(
                text(&run(
                    &["depth", input, "--region", "chr2:1-3", "--all-positions"],
                    0
                )),
                "chr2\t1\t0\nchr2\t2\t0\nchr2\t3\t0\n"
            );
            assert_eq!(text(&run(&["depth", input, "--region", "chr2:1-3"], 0)), "");
            // head: the header, whole or cut, and no record
            let h = text(&run(&["head", input], 0));
            assert!(h.starts_with(hdr), "{h}");
            assert_eq!(h.lines().filter(|l| !l.starts_with('@')).count(), 0);
            assert_eq!(text(&run(&["head", input, "--records", "10"], 0)), h);
            assert_eq!(
                text(&run(
                    &["head", input, "--headers", "1", "--records", "10"],
                    0
                )),
                "@HD\tVN:1.6\tSO:coordinate\n"
            );
        }
        // samples: no read groups, so the sample is "."
        assert_eq!(text(&run(&["samples", &bam], 0)), format!(".\t{bam}\n"));
        assert_eq!(text(&run(&["samples", &sam], 0)), format!(".\t{sam}\n"));
        // bedcov: the BED lines with zero columns; an empty BED gives no rows
        assert_eq!(
            text(&run(&["bedcov", &bed, &bam], 0)),
            "chr1\t0\t100\t0\nchr2\t10\t20\t0\n"
        );
        assert_eq!(
            text(&run(
                &[
                    "bedcov",
                    &bed,
                    &bam,
                    "--read-count",
                    "--bases-at-depth",
                    "1"
                ],
                0
            )),
            "chr1\t0\t100\t0\t0\t0\nchr2\t10\t20\t0\t0\t0\n"
        );
        assert_eq!(text(&run(&["bedcov", &empty_bed, &bam], 0)), "");
        assert_eq!(
            text(&run(&["bedcov", &empty_bed, &bam, "--column-names"], 0)),
            format!("#chrom\tstart\tend\t{bam}_cov\n")
        );
        // check --quick: a header-only file passes (it has an @SQ line and no record to read)
        for input in [&bam, &sam] {
            let o = run(&["check", "--quick", input], 0);
            assert!(text(&o).contains("result: OK"), "{}", text(&o));
        }
        // reheader of a header-only BAM
        let out = s_of(&dir.join(format!("re-{t}.bam")));
        run(&["reheader", &sam, &bam, "-o", &out], 0);
        let o = gloci().args(["head", &out]).output().unwrap();
        assert!(String::from_utf8_lossy(&o.stdout).starts_with(hdr));
        // view --count of a region of a header-only BAM is still zero
        assert_eq!(
            text(&run(&["view", "--count", &bam, "chr1:1-1000"], 0)),
            "0\n"
        );
    }

    // A zero-byte file named .bam or .sam, or a zero-byte stdin, given to a command that reads
    // alignments: it is empty, not "BED" (content detection calls a file with no bytes an
    // empty BED, which is not the problem here).
    let zbam = s_of(&plain_file(dir, "z.bam", ""));
    let zsam = s_of(&plain_file(dir, "z.sam", ""));
    let stdin = "-".to_string();
    for t in ["1", "4"] {
        for (input, name) in [
            (&zbam, zbam.as_str()),
            (&zsam, zsam.as_str()),
            (&stdin, "stdin"),
        ] {
            let cases: Vec<Vec<&str>> = vec![
                vec!["idxstats", input],
                vec!["flagstat", input],
                vec!["samples", input],
                vec!["depth", input],
                vec!["bedcov", &bed, input],
            ];
            for args in cases {
                let mut cmd = gloci();
                cmd.args(&args).write_stdin(Vec::new());
                if args[0] != "samples" {
                    cmd.args(["-t", t]);
                }
                let o = cmd.output().unwrap();
                let err = String::from_utf8_lossy(&o.stderr);
                assert_eq!(o.status.code(), Some(1), "{args:?} -t {t}: {err}");
                assert_eq!(
                    err,
                    format!("error: {name} is empty: it has no records\n"),
                    "{args:?} -t {t}"
                );
                assert!(o.stdout.is_empty(), "{args:?} -t {t}");
                ran.set(ran.get() + 1);
            }
        }
    }

    // A FASTA with no sequences. An empty file is not a FASTA to gloci: dict says so (exit 1),
    // and faidx has nothing to look in. A header with no bases cannot be indexed, and dict
    // writes it as a sequence of length 0 with the MD5 of nothing.
    let fa = s_of(&plain_file(dir, "empty.fa", ""));
    plain_file(dir, "empty.fa.fai", "");
    let hx = s_of(&plain_file(dir, "hx.fa", ">x\n"));
    for t in ["1", "4"] {
        gloci()
            .args(["dict", &fa, "-t", t])
            .assert()
            .code(1)
            .stderr(starts_with_error().and(contains("empty.fa is empty: it has no records")))
            .stdout("");
        let o = gloci().args(["dict", &hx, "-t", t]).output().unwrap();
        assert_eq!(o.status.code(), Some(0));
        let text = String::from_utf8_lossy(&o.stdout);
        assert!(
            text.starts_with("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:x\tLN:0\tM5:d41d8cd98f00b204e9800998ecf8427e\tUR:file://"),
            "{text}"
        );
        // faidx: no region is a usage error; a region names a sequence the index lacks
        gloci()
            .args(["faidx", &fa, "-t", t])
            .assert()
            .code(2)
            .stderr(contains("no region given").and(contains("gloci index")))
            .stdout("");
        gloci()
            .args(["faidx", &fa, "chr1", "-t", t])
            .assert()
            .code(1)
            .stderr(contains("chr1").and(contains("is not in the index")))
            .stdout("");
        gloci()
            .args(["faidx", &fa, "chr1", "--skip-invalid", "-t", t])
            .assert()
            .code(0)
            .stdout("")
            .stderr(contains("skipped 1 region on sequences not in"));
        gloci()
            .args(["faidx", &fa, "--regions-file", &empty_bed, "-t", t])
            .assert()
            .code(0)
            .stdout("")
            .stderr("");
        ran.set(ran.get() + 5);
    }
    gloci()
        .args(["index", &hx])
        .assert()
        .code(1)
        .stderr(contains("has no bases"));
    assert!(!dir.join("hx.fa.fai").exists());
    assert!(ran.get() > 80, "ran {}", ran.get());
}

/// Each damaged alignment file fails every P3c command that reads it: exit 1, one `error:` line
/// that says what is wrong (`check --quick` reports it in its summary), nothing left behind, at
/// `-t 1` and `-t 4`. A BAM without its BGZF end-of-file block and a CRAM cut short or without its
/// end-of-file container are caught before any result is written.
#[test]
fn p3c_failure_matrix_bad_inputs() {
    let fx = P3a::new();
    let dir = fx.path("bad");
    std::fs::create_dir(&dir).unwrap();
    let dir = Path::new(&dir);
    let cram_bytes = std::fs::read(&fx.cram).unwrap();
    let half_cram = s_of(&write(
        dir,
        "half.cram",
        &cram_bytes[..cram_bytes.len() / 2],
    ));
    let noeof_cram = s_of(&write(
        dir,
        "noeof.cram",
        &cram_bytes[..cram_bytes.len() - 38],
    ));
    // a BAM with its index beside it, so that bedcov gets as far as the end-of-file check
    let noeof_bam = dir.join("noeof.bam");
    std::fs::copy(p3a_aln("cat-noeof.bam"), &noeof_bam).unwrap();
    std::fs::copy(format!("{}.bai", fx.bam), dir.join("noeof.bam.bai")).unwrap();
    let noeof_bam = s_of(&noeof_bam);
    let noeof_bam_msg = "BGZF end-of-file marker missing";
    let cram_msg = "the CRAM end-of-file container is missing";
    let hdr = s_of(&plain_file(
        dir,
        "new.sam",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1\n",
    ));
    let bed = fx.bed.as_str();
    let ref_args = ["--reference", fx.fa.as_str()];
    let cases: Vec<(&str, &String, &str)> = vec![
        ("BAM without its EOF block", &noeof_bam, noeof_bam_msg),
        ("truncated CRAM", &half_cram, cram_msg),
        ("CRAM without its EOF container", &noeof_cram, cram_msg),
    ];
    let mut ran = 0;
    for (what, bad, msg) in cases {
        let is_cram = bad.ends_with(".cram");
        // (name, template, whether `-o FILE` is accepted, whether the file is an alignment reader)
        let cmds: Vec<(&str, Vec<&str>, bool)> = vec![
            ("idxstats", vec!["idxstats", "{IN}"], false),
            ("flagstat", vec!["flagstat", "{IN}"], true),
            ("head", vec!["head", "{IN}", "--records", "1"], false),
            ("samples", vec!["samples", "{IN}"], true),
            ("view --count", vec!["view", "--count", "{IN}"], true),
            ("depth", vec!["depth", "{IN}"], true),
            ("bedcov", vec!["bedcov", bed, "{IN}"], true),
            ("reheader", vec!["reheader", &hdr, "{IN}"], true),
            ("check --quick", vec!["check", "--quick", "{IN}"], false),
        ];
        for (name, mut tpl, takes_o) in cmds {
            if is_cram && !matches!(name, "samples" | "reheader" | "check --quick") {
                tpl.extend(ref_args);
            }
            for t in ["1", "4"] {
                let out = tempdir();
                let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", bad)).collect();
                if takes_o {
                    args.extend(["-o".into(), s_of(&out.path().join("out.txt"))]);
                }
                args.extend(["-t".into(), t.into()]);
                if name == "samples" {
                    args.truncate(args.len() - 2); // no -t
                }
                let o = gloci().args(&args).output().unwrap();
                let (err, so) = (
                    String::from_utf8_lossy(&o.stderr),
                    String::from_utf8_lossy(&o.stdout),
                );
                assert_eq!(o.status.code(), Some(1), "{name} / {what} -t {t}: {err}");
                if name == "check --quick" {
                    assert!(
                        so.contains("failed the quick check") && so.contains("FAILED"),
                        "{name} / {what}: {so}"
                    );
                } else {
                    assert!(
                        err.starts_with("error: ") && err.contains(msg),
                        "{name} / {what} -t {t}: wanted `{msg}` in: {err}"
                    );
                    assert!(so.is_empty(), "{name} / {what} -t {t}: stdout {so}");
                }
                assert_eq!(
                    std::fs::read_dir(out.path()).unwrap().count(),
                    0,
                    "{name} / {what} -t {t}: left a file behind"
                );
                ran += 1;
            }
        }
    }
    // reheader's -o may name IN: the damaged file is left as it was
    let before = std::fs::read(&noeof_bam).unwrap();
    gloci()
        .args(["reheader", &hdr, &noeof_bam, "-o", &noeof_bam])
        .assert()
        .code(1)
        .stderr(contains(noeof_bam_msg));
    assert_eq!(std::fs::read(&noeof_bam).unwrap(), before);
    let keep = [
        "half.cram",
        "noeof.bam",
        "noeof.bam.bai",
        "new.sam",
        "noeof.cram",
    ];
    assert_eq!(leftovers(dir, &keep), Vec::<String>::new());
    assert_eq!(ran, 3 * 9 * 2);
}

/// An unsorted BAM or CRAM is refused by `depth` with the `gloci sort` hint, naming the records,
/// at every thread count and for several inputs; nothing is written. A damaged FASTA is refused
/// by `faidx` and `dict` the same way.
#[test]
fn p3c_failure_matrix_unsorted_and_fasta() {
    let fx = P3a::new();
    let d = fx.path("uns");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let ref_fa = p3a_fa("ref.fa");
    for (kind, bad, extra) in [
        ("bam", p3a_aln("sort.bam"), vec![]),
        (
            "cram",
            p3a_aln("unsorted.cram"),
            vec!["--reference", ref_fa.as_str()],
        ),
    ] {
        for inputs in [vec![&bad], vec![&bad, &bad]] {
            for t in ["1", "4"] {
                let out = d.join(format!("{kind}-{}-{t}.out", inputs.len()));
                let mut args: Vec<String> = vec!["depth".into()];
                args.extend(inputs.iter().map(|s| s.to_string()));
                args.extend(extra.iter().map(|s| s.to_string()));
                args.extend(["-o".into(), s_of(&out), "-t".into(), t.into()]);
                let o = gloci().args(&args).output().unwrap();
                let err = String::from_utf8_lossy(&o.stderr);
                assert_eq!(o.status.code(), Some(1), "{kind} {t}: {err}");
                assert!(
                    err.starts_with("error: ")
                        && err.contains("is not sorted: record")
                        && err.contains("comes after record")
                        && err.contains("hint: run: gloci sort "),
                    "{kind} {t}: {err}"
                );
                assert!(!out.exists(), "{kind} {t}: partial output");
            }
        }
    }
    assert_eq!(leftovers(d, &[]).len(), 0, "{:?}", leftovers(d, &[]));

    // FASTA files faidx and dict cannot use
    let nl = plain_file(d, "nolines.fa", "ACGT\n");
    let dup = p3a_fa("dup-name.fa");
    for t in ["1", "4"] {
        gloci()
            .args(["dict", &s_of(&nl), "-t", t])
            .assert()
            .code(1)
            .stderr(starts_with_error())
            .stdout("");
        gloci()
            .args(["dict", &dup, "-t", t])
            .assert()
            .code(1)
            .stderr(starts_with_error().and(contains("was already used on line 1")));
        // to a file, nothing is left behind
        let out = d.join(format!("dup-{t}.dict"));
        gloci()
            .args(["dict", &dup, "-t", t, "-o", &s_of(&out)])
            .assert()
            .code(1);
        assert!(!out.exists(), "partial dict");
    }
}

/// A missing index: `bedcov` (BAM or CRAM) exits 2 naming `gloci index`, `faidx` without its
/// `.fai` or `.gzi` exits 2, and `idxstats`, `depth --region` and `view --count` with a region
/// scan with one warning and give the indexed answer, at every thread count.
#[test]
fn p3c_failure_matrix_missing_index() {
    let fx = P3a::new();
    let d = fx.path("noidx");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let bare_bam = s_of(&d.join("many.bam"));
    std::fs::copy(&fx.bam, &bare_bam).unwrap();
    let bare_cram = s_of(&d.join("many.cram"));
    std::fs::copy(&fx.cram, &bare_cram).unwrap();
    assert!(!Path::new(&format!("{bare_bam}.bai")).exists());
    assert!(!Path::new(&format!("{bare_cram}.crai")).exists());
    for t in ["1", "4"] {
        for (kind, bare, ix) in [("bam", &bare_bam, "bai"), ("cram", &bare_cram, "crai")] {
            let mut args = vec!["bedcov", fx.bed.as_str(), bare.as_str(), "-t", t];
            if kind == "cram" {
                args.extend(["--reference", fx.fa.as_str()]);
            }
            let o = gloci().args(&args).output().unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(2), "bedcov {kind}: {err}");
            assert!(
                err.starts_with("error: ")
                    && err.contains("reads each region through an index")
                    && err.contains(&format!("gloci index {bare}")),
                "bedcov {kind} ({ix}): {err}"
            );
            assert!(o.stdout.is_empty());
            // a second file with an index does not rescue the first
            let mut args = vec!["bedcov", fx.bed.as_str(), bare.as_str()];
            args.push(if kind == "bam" { &fx.bam } else { &fx.cram });
            args.extend(["-t", t]);
            if kind == "cram" {
                args.extend(["--reference", fx.fa.as_str()]);
            }
            gloci().args(&args).assert().code(2);
        }
        // the scan answers give what the index gives
        let want_idx = gloci().args(["idxstats", &fx.bam]).output().unwrap().stdout;
        let o = gloci()
            .args(["idxstats", &bare_bam, "-t", t])
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(0));
        assert_eq!(o.stdout, want_idx, "idxstats scan");
        assert_eq!(
            String::from_utf8_lossy(&o.stderr)
                .matches("warning:")
                .count(),
            1,
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
        for tpl in [
            vec![
                "depth",
                "{IN}",
                "--region",
                "chr1:1-500",
                "--region",
                "chr2:100-900",
            ],
            vec!["depth", "{IN}", "--regions-file", fx.bed.as_str()],
            vec!["view", "--count", "{IN}", "chr1:1-500", "chr2:100-900"],
        ] {
            let run = |input: &String| {
                let a: Vec<String> = tpl.iter().map(|x| x.replace("{IN}", input)).collect();
                gloci().args(&a).args(["-t", t]).output().unwrap()
            };
            let (with, without) = (run(&fx.bam), run(&bare_bam));
            assert_eq!(without.status.code(), Some(0), "{tpl:?}");
            assert!(with.stdout == without.stdout, "{tpl:?} -t {t}");
            assert!(!with.stdout.is_empty());
            let err = String::from_utf8_lossy(&without.stderr);
            assert_eq!(err.matches("scanning the whole file").count(), 1, "{err}");
            assert!(err.contains(&format!("gloci index {bare_bam}")), "{err}");
            assert!(
                String::from_utf8_lossy(&with.stderr).is_empty(),
                "{}",
                String::from_utf8_lossy(&with.stderr)
            );
        }
    }
    // faidx without its .fai or .gzi
    let fa = d.join("ref.fa");
    std::fs::copy(&fx.fa, &fa).unwrap();
    let fagz = d.join("ref.fa.gz");
    std::fs::copy(&fx.fa_gz, &fagz).unwrap();
    std::fs::copy(format!("{}.fai", fx.fa_gz), d.join("ref.fa.gz.fai")).unwrap();
    for t in ["1", "4"] {
        for (path, msg) in [(&fa, "ref.fa.fai"), (&fagz, "ref.fa.gz.gzi")] {
            for extra in [vec!["chr1:1-10"], vec!["--regions-file", fx.masks.as_str()]] {
                let o = gloci()
                    .args(["faidx", &s_of(path)])
                    .args(&extra)
                    .args(["-t", t])
                    .output()
                    .unwrap();
                let err = String::from_utf8_lossy(&o.stderr);
                assert_eq!(o.status.code(), Some(2), "faidx {path:?}: {err}");
                assert!(
                    err.starts_with("error: ") && err.contains(msg) && err.contains("gloci index"),
                    "{err}"
                );
                assert!(o.stdout.is_empty());
            }
        }
    }
    assert!(
        !d.join("ref.fa.fai").exists() && !d.join("ref.fa.gz.gzi").exists(),
        "faidx never builds an index implicitly"
    );
    // dict needs no index
    gloci().args(["dict", &s_of(&fa)]).assert().success();
}

/// With the output directory missing every P3c command that writes a file exits 1, names the
/// directory and how to fix it, and creates nothing.
#[test]
fn p3c_failure_matrix_missing_output_dir() {
    let fx = P3a::new();
    let bed = fx.bed.as_str();
    let hdr = s_of(&plain_file(
        Path::new(&fx.path("")),
        "new.sam",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrM\tLN:150\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1234\n@SQ\tSN:chr10\tLN:600\n@SQ\tSN:chrX\tLN:777\n",
    ));
    let cases: Vec<(&str, Vec<&str>, &str, Sink)> = vec![
        (
            "flagstat",
            vec!["flagstat", "{IN}"],
            &fx.bam,
            Sink::File("txt"),
        ),
        (
            "flagstat json",
            vec!["flagstat", "{IN}"],
            &fx.sam,
            Sink::File("json"),
        ),
        (
            "view --count",
            vec!["view", "--count", "{IN}"],
            &fx.bam,
            Sink::File("txt"),
        ),
        ("depth", vec!["depth", "{IN}"], &fx.bam, Sink::File("txt")),
        (
            "depth -aa",
            vec!["depth", "{IN}", "--all-references"],
            &fx.sam,
            Sink::File("txt"),
        ),
        (
            "bedcov",
            vec!["bedcov", bed, "{IN}"],
            &fx.bam,
            Sink::File("txt"),
        ),
        (
            "reheader",
            vec!["reheader", &hdr, "{IN}"],
            &fx.bam,
            Sink::File("bam"),
        ),
        (
            "samples",
            vec!["samples", "{IN}"],
            &fx.bam,
            Sink::File("txt"),
        ),
        ("dict", vec!["dict", "{IN}"], &fx.fa, Sink::File("dict")),
        (
            "faidx",
            vec!["faidx", "{IN}", "chr1:1-100"],
            &fx.fa,
            Sink::File("fa"),
        ),
        (
            "faidx bgzf",
            vec!["faidx", "{IN}", "chr1:1-100"],
            &fx.fa_gz,
            Sink::File("fa.gz"),
        ),
    ];
    for (name, tpl, good, sink) in &cases {
        assert_missing_dir(name, tpl, good, sink);
    }
    // a CRAM input too
    assert_missing_dir(
        "reheader cram",
        &["reheader", &hdr, "{IN}"],
        &fx.cram,
        &Sink::File("cram"),
    );
    assert_missing_dir(
        "depth cram",
        &["depth", "{IN}", "--reference", &fx.fa],
        &fx.cram,
        &Sink::File("txt"),
    );
}

/// `--chrom-map` on `depth` and `bedcov` renames as records are read: the output rows (and the
/// BED names `bedcov` echoes) are the original run's with the names mapped, the inverse map takes
/// them back, a BED already in the new names gives the same rows, and a SAM is refused.
#[test]
fn p3c_chrom_map_round_trip() {
    let fx = P3a::new();
    let d = fx.path("cm");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let fwd = s_of(&plain_file(
        d,
        "fwd.tsv",
        "# to Ensembl\nchr1\t1\nchr2\t2\nchrM\tMT\n",
    ));
    let rev = s_of(&plain_file(d, "rev.tsv", "1\tchr1\n2\tchr2\nMT\tchrM\n"));
    let rename = |text: &str| -> String {
        // names are whole first fields
        text.lines()
            .map(|l| {
                let (name, rest) = l.split_once('\t').unwrap_or((l, ""));
                let new = match name {
                    "chrM" => "MT",
                    "chr1" => "1",
                    "chr2" => "2",
                    n => n,
                };
                if rest.is_empty() && !l.contains('\t') {
                    format!("{new}\n")
                } else {
                    format!("{new}\t{rest}\n")
                }
            })
            .collect()
    };
    let orig_bed = s_of(&plain_file(
        d,
        "orig.bed",
        "chrM\t0\t100\nchr1\t0\t300\nchr1\t600\t900\nchr2\t100\t900\nchr10\t0\t500\nchrX\t100\t700\n",
    ));
    let ens_bed = s_of(&plain_file(
        d,
        "ens.bed",
        &rename(&std::fs::read_to_string(&orig_bed).unwrap()),
    ));
    assert!(
        std::fs::read_to_string(&ens_bed)
            .unwrap()
            .starts_with("MT\t0\t100\n1\t0\t300\n")
    );
    let ref_args = ["--reference", fx.fa.as_str()];
    for (kind, input, _) in fx.kinds() {
        let ex: Vec<&str> = if kind == "cram" {
            ref_args.to_vec()
        } else {
            vec![]
        };
        if kind == "sam" {
            // a SAM is not rewritten: refused for both commands' alignment inputs, not renamed
            // by halves
            gloci()
                .args(["depth", &input, "--chrom-map", &fwd])
                .assert()
                .code(2)
                .stderr(contains(
                    "--chrom-map cannot rename the reference names in SAM text",
                ));
            continue;
        }
        let run = |args: Vec<&str>, t: &str| -> String {
            let o = gloci().args(&args).args(["-t", t]).output().unwrap();
            assert_eq!(
                o.status.code(),
                Some(0),
                "{kind} {args:?} -t {t}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            String::from_utf8(o.stdout).unwrap()
        };
        // depth: rows in the new names, and the new-name regions find them
        let plain = run([vec!["depth", &input], ex.clone()].concat(), "1");
        let want = rename(&plain);
        assert!(want.starts_with("MT\t") && want.contains("\n1\t") && want.contains("\nchr10\t"));
        for t in ["1", "4"] {
            let got = run(
                [vec!["depth", &input, "--chrom-map", &fwd], ex.clone()].concat(),
                t,
            );
            assert_eq!(got, want, "{kind} depth -t {t}");
        }
        let region = run(
            [
                vec![
                    "depth",
                    &input,
                    "--chrom-map",
                    &fwd,
                    "--region",
                    "2:100-120",
                ],
                ex.clone(),
            ]
            .concat(),
            "4",
        );
        assert!(
            !region.is_empty() && region.lines().all(|l| l.starts_with("2\t")),
            "{region}"
        );
        let orig_region = run(
            [
                vec!["depth", &input, "--region", "chr2:100-120"],
                ex.clone(),
            ]
            .concat(),
            "4",
        );
        assert_eq!(region, rename(&orig_region), "{kind}: region rows");
        // the old name is no longer a region of the renamed file
        let o = gloci()
            .args([
                "depth",
                &input,
                "--chrom-map",
                &fwd,
                "--region",
                "chr2:100-120",
            ])
            .args(&ex)
            .output()
            .unwrap();
        assert!(o.stdout.is_empty(), "{kind}");
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("skipped"),
            "{kind}"
        );
        // map a file already in the new names with the inverse map: the same rows as the original
        let bam_ens = s_of(&d.join(format!("{kind}-ens.bam")));
        gloci()
            .args(["view", &input, "--chrom-map", &fwd, "-o", &bam_ens])
            .args(&ex)
            .assert()
            .success();
        gloci().args(["index", &bam_ens]).assert().success();
        assert_eq!(
            run(vec!["depth", &bam_ens, "--chrom-map", &rev], "4"),
            plain,
            "{kind}: depth through the inverse map"
        );
        assert_eq!(
            run(vec!["depth", &bam_ens], "4"),
            want,
            "{kind}: depth of the renamed file"
        );

        // bedcov: the BED lines are renamed too, so either naming of the BED gives the new rows
        let plain = run(
            [vec!["bedcov", &orig_bed, &input], ex.clone()].concat(),
            "1",
        );
        let want = rename(&plain);
        assert!(want.starts_with("MT\t0\t100\t") && want.contains("\n1\t0\t300\t"));
        for t in ["1", "4"] {
            for bed in [&orig_bed, &ens_bed] {
                let got = run(
                    [vec!["bedcov", bed, &input, "--chrom-map", &fwd], ex.clone()].concat(),
                    t,
                );
                assert_eq!(got, want, "{kind} bedcov {bed} -t {t}");
            }
            // the renamed BAM with the inverse map: the original rows
            let got = run(vec!["bedcov", &orig_bed, &bam_ens, "--chrom-map", &rev], t);
            assert_eq!(got, plain, "{kind} bedcov inverse map -t {t}");
            // and without a map, the renamed file wants BED names in its own naming
            let got = run(vec!["bedcov", &ens_bed, &bam_ens], t);
            assert_eq!(got, want, "{kind} bedcov of the renamed file -t {t}");
        }
        // a map that gives two references one name is refused before any row is written
        let clash = s_of(&plain_file(d, "clash.tsv", "chr1\tchr2\n"));
        for cmd in [vec!["depth", &input], vec!["bedcov", &orig_bed, &input]] {
            gloci()
                .args(&cmd)
                .args(["--chrom-map", &clash])
                .args(&ex)
                .assert()
                .code(2)
                .stderr(contains("which would merge them"))
                .stdout("");
        }
    }
}
