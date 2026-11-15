//! Matrices of the alignment and FASTA commands.

use crate::common::*;
use crate::fixtures_alignments::*;
use crate::helpers::*;
use gloci_testkit::fixtures::write;
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use std::process::Stdio;

/// Every alignment and FASTA command gives the same bytes at `-t 1`, `-t2` and `--threads=8`,
/// on BAM (indexed), SAM and CRAM (indexed); the three formats of the same reads agree for the
/// commands whose output does not name the format.
#[test]
fn p3a_determinism_matrix() {
    let fx = P3a::new();
    let mut ran = 0;
    let mut by_cmd: std::collections::BTreeMap<String, Vec<(String, Vec<u8>)>> = Default::default();
    for (kind, input, extra) in fx.kinds() {
        for (name, mut args) in fx.alignment_cmds(kind, &input) {
            args.extend(extra.iter().cloned());
            let outs: Vec<Vec<u8>> = ["1", "2", "8"]
                .iter()
                .map(|t| {
                    let o = p3a_run_t(&args, t);
                    assert_eq!(
                        o.status.code(),
                        Some(0),
                        "{kind} {name} -t {t}: {}",
                        String::from_utf8_lossy(&o.stderr)
                    );
                    p3a_norm(&o.stdout)
                })
                .collect();
            assert!(
                outs[0] == outs[1] && outs[0] == outs[2],
                "{kind} {name}: output differs between thread counts"
            );
            ran += 1;
            by_cmd
                .entry(name)
                .or_default()
                .push((kind.to_string(), outs[0].clone()));
        }
    }
    // BED output does not name the format: BAM, SAM and CRAM of the same reads agree.
    for name in [
        "convert",
        "convert --bed12",
        "convert --fastq",
        "coverage",
        "genomecov",
        "genomecov bg",
        "intersect B",
        "intersect counts",
    ] {
        let runs = &by_cmd[name];
        assert_eq!(runs.len(), 3, "{name}");
        assert!(!runs[0].1.is_empty(), "{name} wrote nothing");
        assert!(
            runs.iter().all(|r| r.1 == runs[0].1),
            "{name}: formats disagree"
        );
    }
    // The same reads sorted in each format are the same reads.
    for name in ["sort", "sort -n", "view"] {
        let runs = &by_cmd[name];
        let recs = |b: &[u8]| -> Vec<String> {
            String::from_utf8_lossy(b)
                .lines()
                .filter(|l| !l.starts_with('@'))
                .map(String::from)
                .collect()
        };
        assert!(
            runs.iter().all(|r| recs(&r.1) == recs(&runs[0].1)),
            "{name}: formats disagree"
        );
    }
    for (name, args) in fx.fasta_cmds() {
        let outs: Vec<Vec<u8>> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = p3a_run_t(&args, t);
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{name} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                o.stdout
            })
            .collect();
        assert!(!outs[0].is_empty(), "{name}: no output");
        assert!(
            outs[0] == outs[1] && outs[0] == outs[2],
            "{name}: output differs between thread counts"
        );
        ran += 1;
    }
    assert!(ran > 50, "ran {ran} cases");

    // Written files (-o) are the same bytes at every thread count too: BAM output compressed
    // by the parallel writer, and a BGZF FASTA.
    let files: Vec<(&str, Vec<String>)> = vec![
        ("sort", vec!["sort".into(), fx.bam.clone()]),
        (
            "sort spill",
            vec![
                "sort".into(),
                "--max-mem".into(),
                "16M".into(),
                fx.sam.clone(),
                "-O".into(),
                "bam".into(),
            ],
        ),
        (
            "convert-view",
            vec![
                "view".into(),
                fx.cram.clone(),
                "--reference".into(),
                fx.fa.clone(),
            ],
        ),
        (
            "maskfasta",
            vec![
                "maskfasta".into(),
                "--fasta".into(),
                fx.fa.clone(),
                "-b".into(),
                fx.masks.clone(),
            ],
        ),
    ];
    for (name, args) in files {
        let ext = if name == "maskfasta" { "fa.gz" } else { "bam" };
        let outs: Vec<Vec<u8>> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let out = fx.path(&format!("det-{name}-{t}.{ext}"));
                let mut a = args.clone();
                a.extend(["-o".into(), out.clone()]);
                let o = p3a_run_t(&a, t);
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{name} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                let bytes = std::fs::read(&out).unwrap();
                // the @PG line names -t; compare the records and the compressed layout of FASTA
                if ext == "bam" {
                    p3a_norm(&bytes)
                } else {
                    bytes
                }
            })
            .collect();
        assert!(
            outs[0] == outs[1] && outs[0] == outs[2],
            "{name}: written file differs between thread counts"
        );
    }
}

/// Closing stdout after the first bytes is a quiet success (exit 0, nothing on stderr) for every
/// alignment and FASTA command; each command writes far more than a pipe holds, so the reader
/// really does leave early.
#[test]
fn p3a_broken_pipe_matrix() {
    let d = tempdir();
    let (big_fa, big_bed) = p3a_big_fasta(d.path());
    let (big_sam, big_bam) = p3a_big_alignments(d.path(), 12_000, false);
    let sorted_bam = big_bam.clone();
    let cram = p3a_aln("manyref.cram");
    let cram_fa = p3a_fa("manyref.fa");
    let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
    let one_per_chrom = d.path().join("chroms.bed");
    std::fs::write(&one_per_chrom, "chr1\t0\t10000000\nchr2\t0\t10000000\n").unwrap();
    let chroms = s_of(&one_per_chrom);
    let mut cmds: Vec<(String, Vec<String>)> = Vec::new();
    for (kind, input) in [("bam", &big_bam), ("sam", &big_sam)] {
        for (name, a) in [
            ("convert", v(&["convert", input])),
            ("convert --fastq", v(&["convert", input, "--fastq"])),
            ("view", v(&["view", "-O", "sam", input])),
            ("view bam", v(&["view", "-O", "bam", input])),
            ("intersect A", v(&["intersect", "-a", input, "-b", &chroms])),
            (
                "intersect B",
                v(&["intersect", "-a", &chroms, "-b", input, "--mode", "pairs"]),
            ),
            ("genomecov per-base", v(&["genomecov", input, "--per-base"])),
            ("sort", v(&["sort", input])),
            ("sort -n", v(&["sort", "-n", input])),
        ] {
            cmds.push((format!("{kind} {name}"), a));
        }
    }
    cmds.push((
        "bam sort -m".into(),
        v(&["sort", "-m", &sorted_bam, &sorted_bam]),
    ));
    cmds.push(("bam cat".into(), v(&["cat", &sorted_bam, &sorted_bam])));
    cmds.push((
        "bam view region".into(),
        v(&["view", "-O", "sam", &big_bam, "chr1:1-10000000"]),
    ));
    for (name, a) in [
        (
            "view",
            v(&["view", "-O", "sam", &cram, "--reference", &cram_fa]),
        ),
        (
            "genomecov per-base",
            v(&["genomecov", &cram, "--per-base", "--reference", &cram_fa]),
        ),
        (
            "sort",
            v(&["sort", &cram, "--reference", &cram_fa, "-O", "sam"]),
        ),
    ] {
        cmds.push((format!("cram {name}"), a));
    }
    cmds.push((
        "getfasta".into(),
        v(&["getfasta", "--fasta", &big_fa, &big_bed]),
    ));
    cmds.push((
        "getfasta --tab".into(),
        v(&["getfasta", "--fasta", &big_fa, &big_bed, "--tab"]),
    ));
    cmds.push((
        "nuc".into(),
        v(&["nuc", "--fasta", &big_fa, &big_bed, "--seq"]),
    ));
    cmds.push((
        "maskfasta".into(),
        v(&["maskfasta", "--fasta", &big_fa, "-b", &big_bed]),
    ));

    for (name, args) in &cmds {
        let full = gloci().args(args).output().unwrap();
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
        for t in ["1", "4"] {
            let mut child = Command::new(bin())
                .args(args)
                .args(["-t", t])
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

/// A header-only BAM and SAM, an empty BED and FASTA inputs without sequences: every command
/// exits 0 (or says clearly why not), writes a valid empty result and prints nothing on stderr.
#[test]
fn p3a_empty_inputs_matrix() {
    let d = tempdir();
    let dir = d.path();
    let hdr = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n";
    let sam = s_of(&plain_file(dir, "h.sam", hdr));
    let bam = s_of(&dir.join("h.bam"));
    gloci().args(["view", &sam, "-o", &bam]).assert().success();
    gloci().args(["index", &bam]).assert().success();
    let empty_bed = s_of(&plain_file(dir, "e.bed", ""));
    let bed = s_of(&plain_file(dir, "t.bed", "chr1\t0\t100\nchr2\t10\t20\n"));

    for t in ["1", "4"] {
        for input in [&sam, &bam] {
            let ok = |args: &[&str]| -> std::process::Output {
                let o = gloci().args(args).args(["-t", t]).output().unwrap();
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{args:?} -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                assert!(
                    o.stderr.is_empty(),
                    "{args:?} -t {t}: stderr {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                o
            };
            // BED-writing commands: no rows
            for args in [
                vec!["convert", input],
                vec!["convert", input, "--bed12"],
                vec!["convert", input, "--blocks"],
                vec!["convert", input, "--fastq"],
                vec!["intersect", "-a", &bed, "-b", input],
                vec!["intersect", "-a", &bed, "-b", input, "--mode", "pairs"],
                vec!["intersect", "-a", &empty_bed, "-b", input],
                vec!["coverage", "-a", &empty_bed, "-b", input],
            ] {
                assert!(ok(&args).stdout.is_empty(), "{args:?} -t {t}: wrote rows");
            }
            // commands that report on A: one zero row per A record
            let cov = ok(&["coverage", "-a", &bed, "-b", input]);
            assert_eq!(
                String::from_utf8_lossy(&cov.stdout),
                "chr1\t0\t100\t0\t0\t100\t0.0000000\nchr2\t10\t20\t0\t0\t10\t0.0000000\n",
                "coverage -t {t}"
            );
            let gc = ok(&["genomecov", input]);
            let text = String::from_utf8_lossy(&gc.stdout).into_owned();
            assert!(
                text.contains("chr1\t0\t1000\t1000\t1\n"),
                "genomecov {text}"
            );
            assert!(
                text.contains("genome\t0\t1500\t1500\t1\n"),
                "genomecov {text}"
            );
            let bg = ok(&["genomecov", input, "--bedgraph"]);
            assert!(bg.stdout.is_empty(), "genomecov --bedgraph wrote rows");
            let bga = ok(&["genomecov", input, "--bga"]);
            assert_eq!(
                String::from_utf8_lossy(&bga.stdout),
                "chr1\t0\t1000\t0\nchr2\t0\t500\t0\n"
            );
            // alignment-writing commands: the header and no reads
            let sam_out = |args: &[&str]| {
                let o = ok(args);
                let norm = p3a_norm(&o.stdout);
                let text = String::from_utf8_lossy(&norm).into_owned();
                assert!(
                    text.contains("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n"),
                    "{args:?}: {text}"
                );
                assert!(
                    p3a_records(&o.stdout).is_empty(),
                    "{args:?} -t {t}: wrote reads"
                );
            };
            sam_out(&["view", "-O", "sam", input]);
            sam_out(&["view", "-O", "bam", input]);
            if input == &bam {
                sam_out(&["view", "-O", "sam", input, "chr1:1-1000"]);
            }
            sam_out(&["sort", input, "-O", "sam"]);
            sam_out(&["sort", "-n", input, "-O", "sam"]);
            sam_out(&["sort", "-N", input, "-O", "sam"]);
            sam_out(&["sort", "-m", input, input, "-O", "sam"]);
            sam_out(&["intersect", "-a", input, "-b", &empty_bed, "-O", "sam"]);
            sam_out(&["intersect", "-a", input, "-b", &bed, "-O", "sam"]);
            if input == &bam {
                sam_out(&["cat", input, input]);
            }
            // -n and sort -m of the header-only file written to a file is a BAM `view` reads back
            let out = s_of(&dir.join(format!("sorted-{t}.bam")));
            ok(&["sort", input, "-o", &out, "--write-index"]);
            assert!(dir.join(format!("sorted-{t}.bam.bai")).exists());
            assert!(p3a_records(&std::fs::read(&out).unwrap()).is_empty());
        }
    }

    // reads on one contig only: the other contig's rows are still there
    let one = format!("{hdr}r1\t0\tchr1\t11\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    let one_sam = s_of(&plain_file(dir, "one.sam", &one));
    let o = gloci()
        .args(["genomecov", &one_sam, "--bga"])
        .output()
        .unwrap();
    assert_eq!(
        String::from_utf8_lossy(&o.stdout),
        "chr1\t0\t10\t0\nchr1\t10\t15\t1\nchr1\t15\t1000\t0\nchr2\t0\t500\t0\n"
    );

    // FASTA without sequences. An empty file is not a FASTA to gloci: it cannot be indexed, and
    // getfasta, nuc and maskfasta say so (exit 1 or 2) instead of writing nothing. index and
    // maskfasta give the one empty-input message (M1).
    let fa = s_of(&plain_file(dir, "empty.fa", ""));
    gloci()
        .args(["index", &fa])
        .assert()
        .code(1)
        .stderr(starts_with_error().and(contains("empty.fa is empty: it has no records")));
    assert!(!dir.join("empty.fa.fai").exists());
    for cmd in ["getfasta", "nuc"] {
        gloci()
            .args([cmd, "--fasta", &fa, &empty_bed])
            .assert()
            .code(2)
            .stderr(contains("has no index").and(contains("gloci index")));
    }
    for b in [&empty_bed, &bed] {
        gloci()
            .args(["maskfasta", "--fasta", &fa, "-b", b])
            .assert()
            .code(1)
            .stderr(contains("is empty: it has no records"));
    }
    // With an (empty) index beside it: nothing to read, so an empty BED is fine and any
    // interval names the chromosome the FASTA lacks.
    plain_file(dir, "empty.fa.fai", "");
    for t in ["1", "4"] {
        for cmd in ["getfasta", "nuc"] {
            let o = gloci()
                .args([cmd, "--fasta", &fa, &empty_bed, "-t", t])
                .output()
                .unwrap();
            assert_eq!(
                o.status.code(),
                Some(0),
                "{cmd} -t {t}: {}",
                String::from_utf8_lossy(&o.stderr)
            );
            assert!(o.stdout.is_empty() && o.stderr.is_empty(), "{cmd}");
            gloci()
                .args([cmd, "--fasta", &fa, &bed, "-t", t])
                .assert()
                .code(1)
                .stderr(contains("chr1 is not in the FASTA index"));
            gloci()
                .args([cmd, "--fasta", &fa, &bed, "--skip-invalid", "-t", t])
                .assert()
                .success()
                .stdout(predicate::function(|b: &[u8]| {
                    b.is_empty()
                        || b.starts_with(b"#1_usercol")
                            && b.iter().filter(|&&c| c == b'\n').count() == 1
                }))
                .stderr(contains("skipped 2 intervals outside the FASTA"));
        }
    }
    // A real FASTA with an empty BED (or empty mask list): getfasta and nuc write nothing,
    // maskfasta copies the FASTA unchanged.
    let ref_fa = p3a_fa("ref.fa");
    for t in ["1", "4"] {
        for cmd in ["getfasta", "nuc"] {
            let o = gloci()
                .args([cmd, "--fasta", &ref_fa, &empty_bed, "-t", t])
                .output()
                .unwrap();
            assert_eq!(o.status.code(), Some(0), "{cmd} -t {t}");
            assert!(o.stdout.is_empty() && o.stderr.is_empty(), "{cmd} -t {t}");
        }
        let o = gloci()
            .args(["maskfasta", "--fasta", &ref_fa, "-b", &empty_bed, "-t", t])
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(0));
        assert_eq!(
            o.stdout,
            std::fs::read(&ref_fa).unwrap(),
            "maskfasta with no masks"
        );
    }
    // A FASTA whose only record has no bases is refused by the index with the line number, and
    // maskfasta keeps the empty record.
    let hdr_only = s_of(&plain_file(dir, "hdr.fa", ">x\n"));
    gloci()
        .args(["index", &hdr_only])
        .assert()
        .code(1)
        .stderr(contains("hdr.fa line 1").and(contains("has no bases")));
    assert!(!dir.join("hdr.fa.fai").exists());
    gloci()
        .args(["maskfasta", "--fasta", &hdr_only, "-b", &empty_bed])
        .assert()
        .success()
        .stdout(">x\n");
}

/// Each bad input fails every alignment command cleanly: exit 1, one `error:` line that says
/// what is wrong, nothing written, at `-t 1` and `-t 4`.
#[test]
fn p3a_failure_matrix_bad_inputs() {
    let fx = P3a::new();
    let dir = fx.path("bad");
    std::fs::create_dir(&dir).unwrap();
    let cut = p3a_cut_record_bam(Path::new(&dir));
    let half_cram = {
        let b = std::fs::read(&fx.cram).unwrap();
        s_of(&write(Path::new(&dir), "half.cram", &b[..b.len() / 2]))
    };
    let noeof_cram = {
        let b = std::fs::read(&fx.cram).unwrap();
        s_of(&write(Path::new(&dir), "noeof.cram", &b[..b.len() - 38]))
    };
    let bad_bams: Vec<(&str, String, &str)> = vec![
        (
            "no EOF block",
            p3a_aln("cat-noeof.bam"),
            "BGZF end-of-file marker missing",
        ),
        ("cut record", cut, "truncated"),
    ];
    let bad_crams: Vec<(&str, String, &str)> = vec![
        (
            "no EOF container",
            noeof_cram,
            "the CRAM end-of-file container is missing",
        ),
        (
            "half a file",
            half_cram,
            "the CRAM end-of-file container is missing",
        ),
        (
            "corrupt container",
            p3a_aln("many.badtail.cram"),
            "container block checksum mismatch",
        ),
    ];
    let bam_bed = s_of(&plain_file(
        Path::new(&dir),
        "chr12.bed",
        "chr1\t0\t1000\nchr2\t0\t1000\n",
    ));
    let ref_args = ["--reference", fx.fa.as_str()];
    let mut ran = 0;
    for (what, bad, msg) in bad_bams.iter().chain(bad_crams.iter()) {
        let is_cram = bad.ends_with(".cram");
        let bed = if is_cram {
            fx.bed.as_str()
        } else {
            bam_bed.as_str()
        };
        let cmds: Vec<(&str, Vec<&str>, Sink)> = vec![
            ("convert", vec!["convert", "{IN}"], Sink::File("bed")),
            (
                "convert --fastq",
                vec!["convert", "{IN}", "--fastq"],
                Sink::File("fq"),
            ),
            ("view", vec!["view", "{IN}"], Sink::File("sam")),
            ("view bam", vec!["view", "{IN}"], Sink::File("bam")),
            (
                "intersect A",
                vec!["intersect", "-a", "{IN}", "-b", bed],
                Sink::File("bam"),
            ),
            (
                "intersect B",
                vec!["intersect", "-a", bed, "-b", "{IN}"],
                Sink::File("bed"),
            ),
            (
                "coverage",
                vec!["coverage", "-a", bed, "-b", "{IN}"],
                Sink::File("tsv"),
            ),
            ("genomecov", vec!["genomecov", "{IN}"], Sink::File("tsv")),
            ("sort", vec!["sort", "{IN}"], Sink::File("bam")),
            ("sort -n", vec!["sort", "-n", "{IN}"], Sink::File("bam")),
            (
                "sort -m",
                vec!["sort", "-m", "{IN}", "{IN}"],
                Sink::File("bam"),
            ),
            ("cat", vec!["cat", "{IN}", "{IN}"], Sink::File("bam")),
        ];
        for (name, mut tpl, sink) in cmds {
            if is_cram && matches!(name, "sort -m" | "cat") {
                continue; // BAM only
            }
            if *what == "cut record" && name == "cat" {
                continue; // cat copies compressed blocks and never decodes a record
            }
            if is_cram {
                tpl.extend(ref_args);
            }
            assert_fails_clean(name, what, &tpl, bad, &sink, msg);
            ran += 1;
        }
    }
    assert!(ran > 40, "ran {ran}");
}

/// An unsorted BAM or CRAM is refused with the `gloci sort` / `--sort` hint, and `--sort` gives
/// what the same file sorted first gives.
#[test]
fn p3a_failure_matrix_unsorted_alignments() {
    let fx = P3a::new();
    let d = fx.path("uns");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let ref_fa = p3a_fa("ref.fa");
    let bam_bed = s_of(&plain_file(d, "b.bed", "chrM\t0\t1000\nchr1\t0\t100000\n"));
    let cram_bed = s_of(&plain_file(d, "c.bed", "chr1\t0\t1000\n"));
    let sorted_bam = s_of(&d.join("sorted.bam"));
    gloci()
        .args(["sort", &p3a_aln("sort.bam"), "-o", &sorted_bam])
        .assert()
        .success();
    let sorted_cram = s_of(&d.join("sorted-cram.bam"));
    gloci()
        .args([
            "sort",
            &p3a_aln("unsorted.cram"),
            "--reference",
            &ref_fa,
            "-o",
            &sorted_cram,
        ])
        .assert()
        .success();
    let cases = [
        ("bam", p3a_aln("sort.bam"), sorted_bam, bam_bed, vec![]),
        (
            "cram",
            p3a_aln("unsorted.cram"),
            sorted_cram,
            cram_bed,
            vec!["--reference", ref_fa.as_str()],
        ),
    ];
    let mut ran = 0;
    for (kind, bad, sorted, bed, extra) in &cases {
        let tpls: Vec<(&str, Vec<&str>, Sink)> = vec![
            (
                "intersect A",
                vec!["intersect", "-a", "{IN}", "-b", bed],
                Sink::File("bam"),
            ),
            (
                "intersect B",
                vec!["intersect", "-a", bed, "-b", "{IN}"],
                Sink::File("bed"),
            ),
            (
                "coverage",
                vec!["coverage", "-a", bed, "-b", "{IN}"],
                Sink::File("tsv"),
            ),
            ("genomecov", vec!["genomecov", "{IN}"], Sink::File("tsv")),
        ];
        for (name, mut tpl, sink) in tpls {
            tpl.extend(extra.iter());
            // refused, with the hint, and nothing written
            assert_fails_clean(name, kind, &tpl, bad, &sink, "is not sorted");
            let probe: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", bad)).collect();
            let err = String::from_utf8_lossy(&gloci().args(&probe).output().unwrap().stderr)
                .into_owned();
            assert!(
                err.contains("hint: run: gloci sort ") && err.contains("or pass --sort"),
                "{kind} {name}: {err}"
            );
            // --sort equals the sorted file, at every thread count
            let want = {
                let a: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", sorted)).collect();
                let o = gloci().args(&a).output().unwrap();
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{kind} {name} sorted: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                p3a_norm(&o.stdout)
            };
            for t in ["1", "4"] {
                let mut a = probe.clone();
                a.extend(["--sort".into(), "-t".into(), t.into()]);
                let o = gloci().args(&a).output().unwrap();
                assert_eq!(
                    o.status.code(),
                    Some(0),
                    "{kind} {name} --sort -t {t}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                assert!(!want.is_empty(), "{kind} {name}: nothing to compare");
                let got = p3a_norm(&o.stdout);
                if name == "intersect A" {
                    // The reads come out in (reference, start) order with ties in input order, not
                    // samtools' strand order, and the header keeps SO:unsorted: compare the reads
                    // as a set (see the report).
                    let set = |b: &[u8]| {
                        let mut v = p3a_records(b);
                        v.sort();
                        v
                    };
                    assert_eq!(set(&got), set(&want), "{kind} {name} --sort -t {t}");
                } else {
                    assert_eq!(got, want, "{kind} {name} --sort -t {t}");
                }
            }
            ran += 1;
        }
    }
    assert_eq!(ran, 8);
    // sort -m refuses an unsorted BAM, naming the record and the fix
    for t in ["1", "4"] {
        gloci()
            .args([
                "sort",
                "-m",
                &p3a_aln("sort.bam"),
                &p3a_aln("sort.bam"),
                "-t",
                t,
            ])
            .assert()
            .code(1)
            .stderr(
                contains("is not in coordinate order")
                    .and(contains("sort it first with: gloci sort")),
            );
    }
}

/// A missing index: a region query on a BAM or CRAM scans the file with one warning and gives the
/// indexed answer; getfasta and nuc refuse a FASTA without its `.fai` (or `.gzi`) and say how to
/// make it; per-chromosome runs at `-t 8` fall back to streaming with the same bytes.
#[test]
fn p3a_failure_matrix_missing_index() {
    let fx = P3a::new();
    let d = fx.path("noidx");
    let d = Path::new(&d);
    std::fs::create_dir(d).unwrap();
    let bare_bam = d.join("many.bam");
    std::fs::copy(&fx.bam, &bare_bam).unwrap();
    let bare_cram = d.join("many.cram");
    std::fs::copy(&fx.cram, &bare_cram).unwrap();
    let (bare_bam, bare_cram) = (s_of(&bare_bam), s_of(&bare_cram));
    assert!(!Path::new(&format!("{bare_bam}.bai")).exists());
    for (kind, indexed, bare, extra) in [
        ("bam", &fx.bam, &bare_bam, vec![]),
        (
            "cram",
            &fx.cram,
            &bare_cram,
            vec!["--reference", fx.fa.as_str()],
        ),
    ] {
        for t in ["1", "8"] {
            for tpl in [
                vec!["view", "-O", "sam", "{IN}", "chr1:1-500", "chr2:100-900"],
                vec![
                    "intersect",
                    "-a",
                    "{IN}",
                    "-b",
                    fx.bed.as_str(),
                    "-O",
                    "sam",
                ],
                vec!["coverage", "-a", fx.bed.as_str(), "-b", "{IN}"],
                vec!["genomecov", "{IN}"],
                vec!["convert", "{IN}"],
                vec!["convert", "{IN}", "--fastq"],
            ] {
                let run = |input: &String| {
                    let mut a: Vec<String> = tpl.iter().map(|x| x.replace("{IN}", input)).collect();
                    a.extend(extra.iter().map(|s| s.to_string()));
                    a.extend(["-t".into(), t.into()]);
                    gloci().args(&a).output().unwrap()
                };
                let (with, without) = (run(indexed), run(bare));
                assert_eq!(
                    without.status.code(),
                    Some(0),
                    "{kind} {tpl:?} -t {t}: {}",
                    String::from_utf8_lossy(&without.stderr)
                );
                assert_eq!(
                    p3a_norm(&with.stdout),
                    p3a_norm(&without.stdout),
                    "{kind} {tpl:?} -t {t}"
                );
                let err = String::from_utf8_lossy(&without.stderr);
                if tpl[0] == "view" {
                    assert_eq!(
                        err.matches("scanning the whole file").count(),
                        1,
                        "{kind} {tpl:?} -t {t}: {err}"
                    );
                    assert!(err.contains(&format!("gloci index {bare}")), "{err}");
                }
            }
        }
    }
    // FASTA without .fai / .gzi
    let fa = d.join("ref.fa");
    std::fs::copy(&fx.fa, &fa).unwrap();
    let fagz = d.join("ref.fa.gz");
    std::fs::copy(&fx.fa_gz, &fagz).unwrap();
    std::fs::copy(format!("{}.fai", fx.fa_gz), d.join("ref.fa.gz.fai")).unwrap();
    let gf = p3a_fa("gf.bed");
    for cmd in ["getfasta", "nuc"] {
        for t in ["1", "4"] {
            let o = gloci()
                .args([cmd, "--fasta", &s_of(&fa), &gf, "-t", t])
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(2), "{cmd}: {err}");
            assert!(
                err.contains("has no index")
                    && err.contains("ref.fa.fai")
                    && err.contains("gloci index"),
                "{err}"
            );
            assert!(o.stdout.is_empty());
            let o = gloci()
                .args([cmd, "--fasta", &s_of(&fagz), &gf, "-t", t])
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(2), "{cmd}: {err}");
            assert!(
                err.contains("ref.fa.gz.gzi") && err.contains("gloci index"),
                "{err}"
            );
            assert!(o.stdout.is_empty());
        }
    }
    assert!(
        !d.join("ref.fa.fai").exists(),
        "gloci never builds an index implicitly"
    );
    // maskfasta streams and needs no index
    gloci()
        .args(["maskfasta", "--fasta", &s_of(&fa), "-b", &fx.masks])
        .assert()
        .success();
}

/// With the output directory missing every alignment and FASTA command exits 1, names the
/// directory and how to fix it, and creates nothing; `sort`'s own `--tmp-dir` is checked before
/// any input is read.
#[test]
fn p3a_failure_matrix_missing_output_dir() {
    let fx = P3a::new();
    let bed = fx.bed.as_str();
    let sorted_bam = fx.path("sorted.bam");
    let ref_args = ["--reference", fx.fa.as_str()];
    let gf = p3a_fa("gf.bed");
    let nuc = p3a_fa("nuc.bed");
    let mut cases: Vec<(&str, Vec<&str>, &str, Sink)> = vec![
        (
            "convert",
            vec!["convert", "{IN}"],
            &fx.bam,
            Sink::File("bed"),
        ),
        (
            "convert --fastq",
            vec!["convert", "{IN}", "--fastq"],
            &fx.bam,
            Sink::File("fq"),
        ),
        ("view", vec!["view", "{IN}"], &fx.bam, Sink::File("sam")),
        ("view bam", vec!["view", "{IN}"], &fx.sam, Sink::File("bam")),
        (
            "intersect A",
            vec!["intersect", "-a", "{IN}", "-b", bed],
            &fx.bam,
            Sink::File("bam"),
        ),
        (
            "intersect B",
            vec!["intersect", "-a", bed, "-b", "{IN}"],
            &fx.sam,
            Sink::File("bed"),
        ),
        (
            "coverage",
            vec!["coverage", "-a", bed, "-b", "{IN}"],
            &fx.bam,
            Sink::File("tsv"),
        ),
        (
            "genomecov",
            vec!["genomecov", "{IN}"],
            &fx.sam,
            Sink::File("tsv"),
        ),
        ("sort", vec!["sort", "{IN}"], &fx.bam, Sink::File("bam")),
        (
            "sort -n",
            vec!["sort", "-n", "{IN}"],
            &fx.sam,
            Sink::File("sam"),
        ),
        (
            "sort --write-index",
            vec!["sort", "--write-index", "{IN}"],
            &fx.bam,
            Sink::File("bam"),
        ),
        (
            "sort -m",
            vec!["sort", "-m", "{IN}", "{IN}"],
            &sorted_bam,
            Sink::File("bam"),
        ),
        (
            "cat",
            vec!["cat", "{IN}", "{IN}"],
            &sorted_bam,
            Sink::File("bam"),
        ),
        (
            "getfasta",
            vec!["getfasta", "--fasta", fx.fa.as_str(), "{IN}"],
            &gf,
            Sink::File("fa"),
        ),
        (
            "getfasta bgzf",
            vec!["getfasta", "--fasta", fx.fa.as_str(), "{IN}"],
            &gf,
            Sink::File("fa.gz"),
        ),
        (
            "nuc",
            vec!["nuc", "--fasta", fx.fa.as_str(), "{IN}"],
            &nuc,
            Sink::File("tsv"),
        ),
        (
            "maskfasta",
            vec!["maskfasta", "--fasta", fx.fa.as_str(), "-b", "{IN}"],
            &fx.masks,
            Sink::File("fa"),
        ),
        (
            "maskfasta --write-index",
            vec![
                "maskfasta",
                "--write-index",
                "--fasta",
                fx.fa.as_str(),
                "-b",
                "{IN}",
            ],
            &fx.masks,
            Sink::File("fa.gz"),
        ),
    ];
    let mut cram_cases: Vec<(&str, Vec<&str>, &str, Sink)> = vec![
        (
            "cram convert",
            vec!["convert", "{IN}"],
            &fx.cram,
            Sink::File("bed"),
        ),
        (
            "cram convert --fastq",
            vec!["convert", "{IN}", "--fastq"],
            &fx.cram,
            Sink::File("fq"),
        ),
        (
            "cram view",
            vec!["view", "{IN}"],
            &fx.cram,
            Sink::File("bam"),
        ),
        (
            "cram intersect",
            vec!["intersect", "-a", "{IN}", "-b", bed],
            &fx.cram,
            Sink::File("bam"),
        ),
        (
            "cram coverage",
            vec!["coverage", "-a", bed, "-b", "{IN}"],
            &fx.cram,
            Sink::File("tsv"),
        ),
        (
            "cram genomecov",
            vec!["genomecov", "{IN}"],
            &fx.cram,
            Sink::File("tsv"),
        ),
        (
            "cram sort",
            vec!["sort", "{IN}"],
            &fx.cram,
            Sink::File("bam"),
        ),
    ];
    for c in &mut cram_cases {
        c.1.extend(ref_args);
    }
    cases.extend(cram_cases);
    for (name, tpl, good, sink) in &cases {
        assert_missing_dir(name, tpl, good, sink);
    }

    // sort's --tmp-dir must exist
    let d = tempdir();
    let out = s_of(&d.path().join("o.bam"));
    for extra in [vec![], vec!["-n"], vec!["--max-mem", "16M"]] {
        let o = gloci()
            .args(["sort", &fx.bam, "--tmp-dir", "/no/such/tmp", "-o", &out])
            .args(&extra)
            .output()
            .unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert_eq!(o.status.code(), Some(1), "{err}");
        assert!(err.contains("/no/such/tmp"), "{err}");
    }
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
}

/// FASTA commands fail cleanly: intervals past a sequence or on a missing one, a truncated BGZF
/// FASTA, a mask past the end; exit 1, one `error:` line, nothing written, at every `-t`.
#[test]
fn p3a_failure_matrix_fasta() {
    let fx = P3a::new();
    let d = fx.path("fabad");
    let dir = Path::new(&d);
    std::fs::create_dir(dir).unwrap();
    let past = s_of(&plain_file(
        dir,
        "past.bed",
        "chr1\t10\t20\nchr1\t990\t1010\n",
    ));
    let absent = s_of(&plain_file(dir, "absent.bed", "chr1\t10\t20\nchrZ\t0\t5\n"));
    let gz = std::fs::read(&fx.fa_gz).unwrap();
    let cut = write(dir, "cut.fa.gz", &gz[..gz.len() / 2]);
    let cut_noeof = write(dir, "noeof.fa.gz", &gz[..gz.len() - 28]);
    let mut ran = 0;
    for (cmd, extra) in [
        ("getfasta", vec![]),
        ("getfasta", vec!["-s", "--name"]),
        ("getfasta", vec!["--split"]),
        ("nuc", vec![]),
    ] {
        for (bad, msg) in [
            (&past, "ends past the 1000 bases of chr1"),
            (&absent, "chrZ is not in the FASTA index"),
        ] {
            let mut tpl = vec![cmd, "--fasta", fx.fa.as_str(), "{IN}"];
            tpl.extend(extra.iter());
            let ext = if cmd == "nuc" { "tsv" } else { "fa" };
            assert_fails_clean(cmd, "interval", &tpl, bad, &Sink::File(ext), msg);
            ran += 1;
        }
    }
    // skipping instead: one warning with the count, and the rest are written
    for t in ["1", "4"] {
        let o = gloci()
            .args([
                "getfasta",
                "--fasta",
                &fx.fa,
                &absent,
                "--skip-invalid",
                "-t",
                t,
            ])
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(0));
        assert_eq!(
            String::from_utf8_lossy(&o.stdout).lines().count(),
            2,
            "-t {t}"
        );
        assert!(
            String::from_utf8_lossy(&o.stderr).contains("skipped 1 interval"),
            "{}",
            String::from_utf8_lossy(&o.stderr)
        );
    }
    // maskfasta
    let mask_past = s_of(&plain_file(
        dir,
        "mask-past.bed",
        "chr1\t10\t20\nchr1\t990\t1010\n",
    ));
    assert_fails_clean(
        "maskfasta",
        "mask past end",
        &["maskfasta", "--fasta", &fx.fa, "-b", "{IN}"],
        &mask_past,
        &Sink::File("fa"),
        "covers 990-1,010",
    );
    for (what, bad) in [("cut BGZF", &cut), ("BGZF without EOF block", &cut_noeof)] {
        assert_fails_clean(
            "maskfasta",
            what,
            &["maskfasta", "--fasta", "{IN}", "-b", &fx.masks],
            &s_of(bad),
            &Sink::File("fa"),
            "",
        );
        ran += 1;
    }
    assert!(ran >= 10);
}

/// A spill file that cannot be written (the temp directory is "full": every file is capped at
/// 1 MiB) fails `sort`, `sort -n` and `intersect --sort` with exit 1 and a message naming the
/// directory, and leaves nothing behind: the spill directory is empty and no output exists.
#[cfg(unix)]
#[test]
fn p3a_full_tmp_dir_during_sort_leaves_nothing() {
    let d = tempdir();
    let (_sam, bam) = p3a_big_alignments(d.path(), 200_000, true);
    let tmp = d.path().join("spill");
    let out_dir = d.path().join("out");
    std::fs::create_dir(&tmp).unwrap();
    std::fs::create_dir(&out_dir).unwrap();
    let bed = s_of(&plain_file(d.path(), "t.bed", "chr1\t0\t5000000\n"));
    let cases: Vec<(&str, Vec<String>)> = vec![
        (
            "sort",
            vec![
                "sort".into(),
                bam.clone(),
                "-o".into(),
                s_of(&out_dir.join("o.bam")),
            ],
        ),
        (
            "sort -n",
            vec![
                "sort".into(),
                "-n".into(),
                bam.clone(),
                "-o".into(),
                s_of(&out_dir.join("o.bam")),
            ],
        ),
        (
            "sort sam",
            vec![
                "sort".into(),
                bam.clone(),
                "-O".into(),
                "sam".into(),
                "-o".into(),
                s_of(&out_dir.join("o.sam")),
            ],
        ),
        (
            "intersect --sort",
            vec![
                "intersect".into(),
                "-a".into(),
                bam.clone(),
                "-b".into(),
                bed.clone(),
                "--sort".into(),
                "-o".into(),
                s_of(&out_dir.join("o.bam")),
            ],
        ),
        (
            "coverage --sort",
            vec![
                "coverage".into(),
                "-a".into(),
                bed.clone(),
                "-b".into(),
                bam.clone(),
                "--sort".into(),
                "-o".into(),
                s_of(&out_dir.join("o.tsv")),
            ],
        ),
    ];
    let mut spilled_once = false;
    for (name, base) in cases {
        for t in if name == "sort" {
            vec!["1", "4"]
        } else {
            vec!["4"]
        } {
            let mut args = base.clone();
            args.extend(["--max-mem", "16M", "--tmp-dir", &s_of(&tmp), "-t", t].map(String::from));
            // unlimited: succeeds and removes its spill files
            let ok = gloci().args(&args).output().unwrap();
            assert_eq!(
                ok.status.code(),
                Some(0),
                "{name} -t {t}: {}",
                String::from_utf8_lossy(&ok.stderr)
            );
            assert_eq!(
                std::fs::read_dir(&tmp).unwrap().count(),
                0,
                "{name} -t {t}: spill files left after success"
            );
            for f in std::fs::read_dir(&out_dir).unwrap() {
                std::fs::remove_file(f.unwrap().path()).unwrap();
            }
            // capped: fails cleanly
            let o = p3a_run_size_limited(&args, 1 << 20);
            let err = String::from_utf8_lossy(&o.stderr);
            assert_eq!(o.status.code(), Some(1), "{name} -t {t}: {err}");
            assert!(err.starts_with("error: "), "{name} -t {t}: {err}");
            assert!(
                err.contains("File too large") || err.contains("No space left"),
                "{name} -t {t}: {err}"
            );
            assert_eq!(
                std::fs::read_dir(&tmp).unwrap().count(),
                0,
                "{name} -t {t}: spill files left: {err}"
            );
            assert_eq!(
                std::fs::read_dir(&out_dir).unwrap().count(),
                0,
                "{name} -t {t}: output left behind"
            );
            spilled_once = true;
        }
    }
    assert!(spilled_once);
}
