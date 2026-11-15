//! Matrices of the statistics, randomness, pairs and VCF commands.

use crate::common::*;
use crate::fixtures_bed_vcf::*;
use crate::fixtures_stats::*;
use crate::helpers::*;
use gloci_testkit::fixtures::{bgzip, write};

/// Every statistics command gives the same bytes at `-t 1`, `-t2` and `--threads=8` (the seeded ones
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
            || name.starts_with("shuffle")
            // compressed BAM of the shared BED, under 100 kB; matrices.rs breaks its pipe
            || name == "convert --bam";
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
    let fst_pop1 = fx.plain("fst-pop1.txt", "s1\n");
    let fst_pop2 = fx.plain("fst-pop2.txt", "s2\ns3\n");
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
        assert!(stdout_with(&["convert", "--fastq", &sam], t).is_empty());
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
        // Fst: the column headers, and a summary of no site
        gloci()
            .args(["vcf", "popgen", v, "--out-prefix", &p])
            .args(["--fst-pop", &fst_pop1, "--fst-pop", &fst_pop2])
            .args(["--fst-window-size", "10"])
            .assert()
            .success();
        for ext in ["weir.fst", "windowed.weir.fst"] {
            let t = std::fs::read_to_string(format!("{p}.{ext}")).unwrap();
            assert_eq!(t.lines().count(), 1, ".{ext}: {t}");
        }
        assert_eq!(
            std::fs::read_to_string(format!("{p}.weir.fst.summary")).unwrap(),
            "N_SITES\tMEAN_FST\tWEIGHTED_FST\n0\tnan\tnan\n"
        );
        // LD: the column headers (vcftools stops: docs/differences.md#ld-rules)
        gloci()
            .args(["vcf", "popgen", v, "--out-prefix", &p])
            .args(["--geno-r2", "--hap-r2", "--interchrom-geno-r2"])
            .args(["--interchrom-hap-r2"])
            .assert()
            .success();
        for ext in [
            "geno.ld",
            "hap.ld",
            "interchrom.geno.ld",
            "interchrom.hap.ld",
        ] {
            let t = std::fs::read_to_string(format!("{p}.{ext}")).unwrap();
            assert_eq!(t.lines().count(), 1, ".{ext}: {t}");
        }
        // relatedness: every pair of the three samples, nan (as vcftools)
        gloci()
            .args(["vcf", "popgen", v, "--out-prefix", &p])
            .args(["--relatedness", "--relatedness2"])
            .assert()
            .success();
        for (ext, rows) in [("relatedness", 1 + 6), ("relatedness2", 1 + 9)] {
            let t = std::fs::read_to_string(format!("{p}.{ext}")).unwrap();
            assert_eq!(t.lines().count(), rows, ".{ext}: {t}");
            assert!(
                t.lines().skip(1).all(|l| l.ends_with("\tnan")),
                ".{ext}: {t}"
            );
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

/// The failure matrix of the statistics commands: unsorted input where order matters, malformed
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
    let fst_pop1 = fx.plain("fst-pop1.txt", "s1\n");
    let fst_pop2 = fx.plain("fst-pop2.txt", "s2\ns3\n");
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
            "vcf popgen fst",
            vec![
                "vcf",
                "popgen",
                "--fst-pop",
                &fst_pop1,
                "--fst-pop",
                &fst_pop2,
                "--fst-window-size",
                "1000",
                "{IN}",
            ],
            Sink::Prefix,
            false,
            // its own message for unsorted input: tests/vcf_popgen.rs
            false,
        ),
        (
            "vcf popgen ld",
            vec!["vcf", "popgen", "--geno-r2", "--ld-window", "5", "{IN}"],
            Sink::Prefix,
            false,
            // its own message for unsorted input: tests/vcf_popgen.rs
            false,
        ),
        (
            "vcf popgen relatedness",
            vec!["vcf", "popgen", "--relatedness", "--relatedness2", "{IN}"],
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

    // a missing output directory: every statistics command that writes a file
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
        (
            "convert --fastq",
            vec!["convert", "--fastq", &big_sam],
            Sink::File("fq"),
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
        (
            "vcf popgen fst",
            vec![
                "vcf",
                "popgen",
                "--fst-pop",
                &fst_pop1,
                "--fst-pop",
                &fst_pop2,
                "{IN}",
            ],
            Sink::Prefix,
        ),
        (
            "vcf popgen ld",
            vec!["vcf", "popgen", "--hap-r2", "--interchrom-geno-r2", "{IN}"],
            Sink::Prefix,
        ),
        (
            "vcf popgen relatedness",
            vec!["vcf", "popgen", "--relatedness", "--relatedness2", "{IN}"],
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
