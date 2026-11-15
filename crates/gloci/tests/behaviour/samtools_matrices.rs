//! Matrices of the samtools-style commands.

use crate::common::*;
use crate::fixtures_alignments::*;
use crate::fixtures_samtools::*;
use crate::helpers::*;
use gloci_testkit::fixtures::write;
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use std::process::Stdio;

/// Every samtools-essentials command gives the same bytes at `-t 1`, `-t2` and `--threads=8`, on BAM (indexed),
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
/// A samtools-essentials command whose output outgrows a pipe: `depth --all-references | head`, `faidx | head`
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

/// Each damaged alignment file fails every samtools-essentials command that reads it: exit 1, one `error:` line
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

/// With the output directory missing every samtools-essentials command that writes a file exits 1, names the
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
