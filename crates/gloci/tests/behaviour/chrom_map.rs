//! The chromosome-map round trip of each command family.

use crate::common::*;
use crate::fixtures_alignments::*;
use crate::fixtures_bed_vcf::*;
use crate::fixtures_stats::*;
use crate::helpers::*;
use predicates::str::contains;
use std::path::Path;

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
