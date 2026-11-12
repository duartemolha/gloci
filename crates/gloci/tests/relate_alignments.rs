//! BAM and SAM in `intersect` and `coverage` (D101, D104, D105, D115, D124): determinism, usage
//! errors, `--sort`, slim records and hand-checked blocks. Byte parity with bedtools is in the
//! golden cases `intersect-bam-*`, `coverage-bam-*` and friends.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn aln(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn run(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// Output without the @PG line, whose command line names the thread count.
fn no_pg(s: &str) -> String {
    s.lines()
        .filter(|l| !l.starts_with("@PG"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn write(dir: &Path, name: &str, text: &str) -> String {
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_string_lossy().into_owned()
}

/// A BAM made from SAM text.
fn bam_of(dir: &Path, name: &str, sam: &str) -> String {
    let s = write(dir, &format!("{name}.sam"), sam);
    let out = dir.join(name).to_string_lossy().into_owned();
    gloci()
        .args(["view", &s, "-O", "bam", "-o", &out])
        .assert()
        .success();
    out
}

const TARGETS: &str =
    "chrM\t50\t120\nchr1\t0\t300\nchr1\t400\t700\nchr2\t100\t101\nchr10\t0\t300\nchrX\t500\t520\n";

#[test]
fn alignment_relations_deterministic() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", TARGETS);
    let bam = aln("order.bam");
    let mut cases: Vec<Vec<String>> = Vec::new();
    for mode in ["any", "none", "count", "pairs"] {
        let mut a_first = vec!["intersect", "-a", &bam, "-b", &bed, "--mode", mode];
        if matches!(mode, "any" | "none") {
            a_first.extend(["-O", "sam"]);
        }
        cases.push(a_first.into_iter().map(String::from).collect());
        let as_b = ["intersect", "-a", &bed, "-b", &bam, "--mode", mode];
        cases.push(as_b.map(String::from).to_vec());
    }
    for extra in [&[][..], &["--hist"], &["--per-base"], &["--split"]] {
        for (a, b) in [(&bed, &bam), (&bam, &bed)] {
            let mut v: Vec<String> = ["coverage", "-a", a, "-b", b].map(String::from).to_vec();
            v.extend(extra.iter().map(|s| s.to_string()));
            cases.push(v);
        }
    }
    for case in cases {
        let at = |t: &str| {
            let mut v = case.clone();
            v.extend(["-t".to_string(), t.to_string()]);
            let refs: Vec<&str> = v.iter().map(String::as_str).collect();
            no_pg(&run(&refs))
        };
        let one = at("1");
        assert!(
            !one.is_empty() || case.contains(&"none".to_string()),
            "{case:?}"
        );
        assert_eq!(one, at("2"), "{case:?}");
        assert_eq!(one, at("8"), "{case:?}");
    }
}

#[test]
fn clip_refused_for_alignments() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", TARGETS);
    gloci()
        .args([
            "intersect",
            "-a",
            &aln("order.bam"),
            "-b",
            &bed,
            "--mode",
            "clip",
        ])
        .assert()
        .code(2)
        .stderr(contains("cannot be clipped"));
}

#[test]
fn split_deletions_needs_alignment() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", "chr1\t0\t10\n");
    for cmd in ["intersect", "coverage"] {
        gloci()
            .args([cmd, "-a", &bed, "-b", &bed, "--split-deletions"])
            .assert()
            .code(2)
            .stderr(contains("--split-deletions"))
            .stderr(contains("no input is BAM, CRAM or SAM"));
    }
}

#[test]
fn read_filters_without_alignment_input_warn_once() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", "chr1\t0\t10\n");
    let g = write(dir.path(), "g.genome", "chr1\t100\n");
    let bam = aln("order.bam");
    let runs: [Vec<&str>; 4] = [
        vec!["intersect", "-a", &bed, "-b", &bed],
        vec!["coverage", "-a", &bed, "-b", &bed],
        vec!["genomecov", &bed, "-g", &g],
        vec!["view", &bed],
    ];
    let one = "--min-mapq is not used: no input is BAM, CRAM or SAM";
    let two = "--min-mapq and --exclude-flags are not used: no input is BAM, CRAM or SAM";
    for args in &runs {
        let stderr = |extra: &[&str]| {
            let o = gloci().args(args).args(extra).assert().success();
            String::from_utf8(o.get_output().stderr.clone()).unwrap()
        };
        let err = stderr(&["--min-mapq", "20"]);
        assert_eq!(err.matches(one).count(), 1, "{args:?}: {err}");
        let err = stderr(&["--min-mapq", "20", "--exclude-flags", "dup"]);
        assert_eq!(err.matches(two).count(), 1, "{args:?}: {err}");
        assert!(!stderr(&[]).contains("not used"), "{args:?}");
    }
    // An alignment input uses them: no warning.
    let o = gloci()
        .args(["intersect", "-a", &bam, "-b", &bed, "--mode", "count"])
        .args(["--min-mapq", "20"])
        .assert()
        .success();
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    assert!(!err.contains("not used"), "{err}");
}

#[test]
fn bed_output_for_alignment_a_refused() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", TARGETS);
    let out = dir.path().join("x.bed");
    gloci()
        .args(["intersect", "-a", &aln("order.bam"), "-b", &bed, "-o"])
        .arg(&out)
        .assert()
        .code(2)
        .stderr(contains("gloci convert -"));
    assert!(!out.exists());
}

#[test]
fn sam_a_with_chrom_map_refused_for_alignment_output() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", TARGETS);
    let map = write(dir.path(), "map.tsv", "chr1\tX1\n");
    gloci()
        .args([
            "intersect",
            "-a",
            &aln("order.sam"),
            "-b",
            &bed,
            "--chrom-map",
            &map,
        ])
        .assert()
        .code(2)
        .stderr(contains("--chrom-map"));
    // Text rows and BAM output are fine.
    gloci()
        .args([
            "intersect",
            "-a",
            &aln("order.sam"),
            "-b",
            &bed,
            "--chrom-map",
            &map,
        ])
        .args(["--mode", "pairs"])
        .assert()
        .success();
    gloci()
        .args([
            "intersect",
            "-a",
            &aln("order.bam"),
            "-b",
            &bed,
            "--chrom-map",
            &map,
        ])
        .args(["-O", "sam"])
        .assert()
        .success();
}

/// order.sam with the chromosomes' blocks in the opposite order (reads within a chromosome
/// keep their order).
fn shuffled_bam(dir: &Path) -> String {
    let text = std::fs::read_to_string(aln("order.sam")).unwrap();
    let header: Vec<&str> = text.lines().filter(|l| l.starts_with('@')).collect();
    let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('@')).collect();
    let mut out = header.join("\n");
    out.push('\n');
    for chunk in body.chunks(3).rev() {
        out.push_str(&chunk.join("\n"));
        out.push('\n');
    }
    bam_of(dir, "shuffled.bam", &out)
}

fn dir_is_empty(p: &Path) -> bool {
    std::fs::read_dir(p).unwrap().next().is_none()
}

#[test]
fn sort_flag_for_alignments() {
    let dir = tempdir();
    let tmp = dir.path().join("tmp");
    std::fs::create_dir(&tmp).unwrap();
    let tmp = tmp.to_string_lossy().into_owned();
    let bed = write(dir.path(), "t.bed", TARGETS);
    let shuffled = shuffled_bam(dir.path());
    let sorted = aln("order.bam");
    let any = |a: &str, extra: &[&str]| {
        let mut v = vec!["intersect", "-a", a, "-b", &bed, "-O", "sam"];
        v.extend(extra);
        no_pg(&run(&v))
    };
    let want = any(&sorted, &[]);
    assert_eq!(
        want,
        any(
            &shuffled,
            &["--sort", "--max-mem", "16M", "--tmp-dir", &tmp]
        )
    );
    assert!(dir_is_empty(Path::new(&tmp)), "temp files left behind");
    let cov = |a: &str, b: &str, extra: &[&str]| {
        let mut v = vec!["coverage", "-a", a, "-b", b, "--hist"];
        v.extend(extra);
        run(&v)
    };
    assert_eq!(
        cov(&bed, &sorted, &[]),
        cov(
            &bed,
            &shuffled,
            &["--sort", "--max-mem", "16M", "--tmp-dir", &tmp]
        )
    );
    assert!(dir_is_empty(Path::new(&tmp)));
    // Without --sort the order error says how to fix it.
    gloci()
        .args(["intersect", "-a", &shuffled, "-b", &bed])
        .args(["-O", "sam"])
        .assert()
        .failure()
        .stderr(contains("gloci sort"))
        .stderr(contains("--sort"));
}

#[test]
fn sort_flag_keeps_sam_a_as_sam() {
    let dir = tempdir();
    let tmp = dir.path().to_string_lossy().into_owned();
    let bed = write(dir.path(), "t.bed", TARGETS);
    let text = std::fs::read_to_string(aln("order.sam")).unwrap();
    let mut lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('@')).collect();
    let header: Vec<&str> = text.lines().filter(|l| l.starts_with('@')).collect();
    lines.reverse();
    let shuffled = write(
        dir.path(),
        "s.sam",
        &format!("{}\n{}\n", header.join("\n"), lines.join("\n")),
    );
    let got = run(&[
        "intersect",
        "-a",
        &shuffled,
        "-b",
        &bed,
        "--sort",
        "--tmp-dir",
        &tmp,
    ]);
    let want = run(&["intersect", "-a", &aln("order.sam"), "-b", &bed]);
    // Ties between a read and its reverse mate-name are broken by input order, which the
    // reversal flipped: compare as sets of lines.
    let set = |s: &str| {
        let mut v: Vec<String> = no_pg(s).lines().map(String::from).collect();
        v.sort();
        v
    };
    assert_eq!(set(&got), set(&want));
    assert!(got.starts_with("@HD"));
}

/// Reads whose interval order (start, end) is not `gloci sort`'s (reference, position, strand),
/// in a header whose @SQ order is not lexicographic.
const UNSORTED_TIES: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n\
a1\t16\tchr1\t101\t60\t10M\t*\t0\t0\t*\t*\n\
a2\t0\tchr1\t101\t60\t20M\t*\t0\t0\t*\t*\n\
a3\t0\tchr1\t101\t60\t5M\t*\t0\t0\t*\t*\n\
a4\t0\tchr2\t50\t60\t5M\t*\t0\t0\t*\t*\n\
a5\t16\tchr1\t101\t60\t5M\t*\t0\t0\t*\t*\n";

/// SAM text without its @PG lines.
fn records_and_hd(s: &str) -> Vec<String> {
    s.lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(String::from)
        .collect()
}

#[test]
fn sort_flag_equals_gloci_sort_then_command() {
    // The help's claim: --sort equals gloci sort on each input followed by the command, for
    // the records and for @HD (SO:coordinate), with A or B an alignment, BAM or SAM.
    let dir = tempdir();
    let tmp = dir.path().to_string_lossy().into_owned();
    let bed = write(dir.path(), "x.bed", "chr2\t0\t1000\nchr1\t0\t1000\n");
    let sam = write(dir.path(), "u.sam", UNSORTED_TIES);
    let bam = bam_of(dir.path(), "u.bam", UNSORTED_TIES);
    for (input, out) in [(&sam, "s.sam"), (&bam, "s.bam")] {
        let sorted = dir.path().join(out).to_string_lossy().into_owned();
        gloci()
            .args(["sort", input, "-o", &sorted])
            .assert()
            .success();
        let a = |file: &str, extra: &[&str]| {
            let mut v = vec!["intersect", "-a", file, "-b", &bed, "-O", "sam"];
            v.extend(extra);
            records_and_hd(&run(&v))
        };
        let got = a(input, &["--sort", "--tmp-dir", &tmp]);
        assert_eq!(got, a(&sorted, &[]), "{input} as A");
        assert!(got[0].contains("SO:coordinate"), "{got:?}");
        let names: Vec<&str> = got
            .iter()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(names, ["a4", "a2", "a3", "a1", "a5"]);
        let b = |file: &str, extra: &[&str]| {
            let mut v = vec!["intersect", "-a", &bed, "-b", file, "--mode", "pairs"];
            v.extend(extra);
            run(&v)
        };
        assert_eq!(
            b(input, &["--sort", "--tmp-dir", &tmp]),
            b(&sorted, &[]),
            "{input} as B"
        );
        // A -g table that orders the references unlike the header: the records follow it,
        // which SO:coordinate would not describe.
        let g = write(dir.path(), "g.genome", "chr1\t1000\nchr2\t1000\n");
        let bed_g = write(dir.path(), "g.bed", "chr1\t0\t1000\nchr2\t0\t1000\n");
        let got = records_and_hd(&run(&[
            "intersect",
            "-a",
            input,
            "-b",
            &bed_g,
            "-g",
            &g,
            "-O",
            "sam",
            "--sort",
            "--tmp-dir",
            &tmp,
        ]));
        assert!(got[0].contains("SO:unsorted"), "{got:?}");
        let names: Vec<&str> = got
            .iter()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(names, ["a2", "a3", "a1", "a5", "a4"]);
    }
}

#[test]
fn sort_flag_with_another_genome_order_never_claims_coordinate() {
    // A coordinate-sorted BAM (its header says SO:coordinate) sorted by --sort under a -g
    // table that reverses the references: the records follow the table, so the output says
    // SO:unsorted, and gloci check passes it.
    let dir = tempdir();
    let tmp = dir.path().to_string_lossy().into_owned();
    let u = bam_of(dir.path(), "u.bam", UNSORTED_TIES);
    let sorted = dir.path().join("s.bam").to_string_lossy().into_owned();
    gloci().args(["sort", &u, "-o", &sorted]).assert().success();
    let hd = |f: &str| {
        run(&["head", f, "--headers", "1", "--records", "0"])
            .lines()
            .next()
            .unwrap()
            .to_string()
    };
    assert!(hd(&sorted).contains("SO:coordinate"));
    let g = write(dir.path(), "rev.genome", "chr1\t1000\nchr2\t1000\n");
    let bed = write(dir.path(), "g.bed", "chr1\t0\t1000\nchr2\t0\t1000\n");
    for out in ["i.bam", "i.sam"] {
        let out = dir.path().join(out).to_string_lossy().into_owned();
        gloci()
            .args(["intersect", "-a", &sorted, "-b", &bed, "-g", &g, "--sort"])
            .args(["--tmp-dir", &tmp, "-o", &out])
            .assert()
            .success();
        let h = hd(&out);
        assert!(
            h.contains("SO:unsorted") && !h.contains("coordinate"),
            "{h}"
        );
        let o = gloci().args(["check", &out]).output().unwrap();
        assert_eq!(
            o.status.code(),
            Some(0),
            "{}{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        );
    }
}

#[test]
fn slim_b_records() {
    use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions, RecordBatch, open_reader};
    let dir = tempdir();
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n");
    let seq = "ACGT".repeat(38);
    let seq = &seq[..150];
    let qual = "I".repeat(150);
    for i in 0..10_000 {
        sam.push_str(&format!(
            "r{i}\t0\tchr1\t1001\t60\t150M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0\tRG:Z:group1\n"
        ));
    }
    let bam = bam_of(dir.path(), "dense.bam", &sam);
    let sam_path = write(dir.path(), "dense.sam", &sam);
    let bytes = |file: &str, full: bool| {
        let ro = ReaderOptions {
            accept: Accept {
                alignments: Some(AlignOptions {
                    full_records: full,
                    ..Default::default()
                }),
                ..Default::default()
            },
            chrom_map: None,
            threads: std::num::NonZeroUsize::new(1).unwrap(),
            format: None,
            vcf_span: gloci_core::VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
        };
        let mut r = open_reader(file, &ro).unwrap();
        assert_eq!(r.whole_records(), full, "{file}");
        let mut chroms = gloci_core::ChromTable::first_appearance();
        let mut batch = RecordBatch::default();
        while r.read_batch(&mut chroms, &mut batch, 4096).unwrap() {}
        assert_eq!(batch.recs.len(), 10_000);
        (batch.buf.len(), batch.raw(0).to_vec())
    };
    let ((slim, _), (full, _)) = (bytes(&bam, false), bytes(&bam, true));
    assert!(slim * 3 <= full, "slim {slim} bytes against full {full}");
    // SAM lines are slimmed too (P9): SEQ and QUAL become *, and the other tags go.
    let ((slim, line), (full, _)) = (bytes(&sam_path, false), bytes(&sam_path, true));
    assert!(
        slim * 3 <= full,
        "SAM: slim {slim} bytes against full {full}"
    );
    assert_eq!(line, b"r0\t0\tchr1\t1001\t60\t150M\t*\t0\t0\t*\t*");
    // and give the same results as the BAM
    let bed = write(dir.path(), "t.bed", "chr1\t1000\t1200\nchr1\t1100\t1101\n");
    for args in [
        vec!["intersect", "-a", &bed, "-b", "@", "--mode", "pairs"],
        vec![
            "intersect",
            "-a",
            &bed,
            "-b",
            "@",
            "--mode",
            "count",
            "--split",
        ],
        vec!["coverage", "-a", &bed, "-b", "@"],
    ] {
        let with = |f: &str| {
            let v: Vec<&str> = args.iter().map(|a| if *a == "@" { f } else { a }).collect();
            run(&v)
        };
        assert_eq!(with(&sam_path), with(&bam), "{args:?}");
    }
}

#[test]
fn unsorted_alignment_b() {
    let dir = tempdir();
    let bed = write(dir.path(), "a.bed", TARGETS);
    let text = std::fs::read_to_string(aln("order.sam")).unwrap();
    let header: Vec<&str> = text.lines().filter(|l| l.starts_with('@')).collect();
    let mut body: Vec<&str> = text.lines().filter(|l| !l.starts_with('@')).collect();
    body.reverse();
    let shuffled = write(
        dir.path(),
        "b.sam",
        &format!("{}\n{}\n", header.join("\n"), body.join("\n")),
    );
    for mode in ["pairs", "count"] {
        let sorted = run(&[
            "intersect",
            "-a",
            &bed,
            "-b",
            &aln("order.sam"),
            "--mode",
            mode,
        ]);
        let got = run(&[
            "intersect",
            "-a",
            &bed,
            "-b",
            &shuffled,
            "--mode",
            mode,
            "--unsorted",
        ]);
        // A's order is kept; within one A, B arrives in B's (reversed) order.
        let set = |s: &str| {
            let mut v: Vec<String> = s.lines().map(String::from).collect();
            v.sort();
            v
        };
        assert_eq!(set(&sorted), set(&got), "{mode}");
        let a_cols = |s: &str| -> Vec<String> {
            s.lines()
                .map(|l| l.split('\t').take(3).collect::<Vec<_>>().join("\t"))
                .collect()
        };
        let mut dedup = a_cols(&got);
        dedup.dedup();
        let mut want = a_cols(&sorted);
        want.dedup();
        assert_eq!(dedup, want, "{mode}: A's order");
    }
}

const SPLIT_SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
skip\t0\tchr1\t101\t60\t5M100N5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
del\t0\tchr1\t301\t60\t5M2D5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";

#[test]
fn split_blocks_hand_checked() {
    let dir = tempdir();
    let sam = write(dir.path(), "r.sam", SPLIT_SAM);
    // skip: blocks [100,105) and [205,210); del: [300,312), or [300,305) and [307,312).
    let cov = |t: &str, extra: &[&str]| {
        let bed = write(dir.path(), "t.bed", t);
        let mut v = vec!["coverage", "-a", &bed, "-b", &sam];
        v.extend(extra);
        run(&v)
    };
    // Whole spans: the intron target is covered.
    assert_eq!(
        cov("chr1\t120\t180\n", &[]),
        "chr1\t120\t180\t1\t60\t60\t1\n"
    );
    // Blocks: inside the intron, nothing.
    assert_eq!(
        cov("chr1\t120\t180\n", &["--split"]),
        "chr1\t120\t180\t0\t0\t60\t0\n"
    );
    // Bookended to a block on both sides: no overlap.
    assert_eq!(
        cov("chr1\t105\t205\n", &["--split"]),
        "chr1\t105\t205\t0\t0\t100\t0\n"
    );
    // Over only the deletion: covered by the span and by --split, not by --split-deletions.
    assert_eq!(
        cov("chr1\t305\t307\n", &["--split"]),
        "chr1\t305\t307\t1\t2\t2\t1\n"
    );
    assert_eq!(
        cov("chr1\t305\t307\n", &["--split-deletions"]),
        "chr1\t305\t307\t0\t0\t2\t0\n"
    );
    // A spanning the deletion: 10 of 12 bases.
    assert_eq!(
        cov("chr1\t300\t312\n", &["--split-deletions"]),
        "chr1\t300\t312\t1\t10\t12\t0.8333333333333334\n"
    );
    // A point inside the deletion (a point never has bases).
    assert_eq!(cov("chr1\t306\t306\n", &[]), "chr1\t306\t306\t1\t0\t0\t.\n");
    // --split-deletions with A as the alignments.
    let bed = write(dir.path(), "d.bed", "chr1\t305\t307\n");
    let n = |extra: &[&str]| {
        let mut v = vec!["intersect", "-a", &sam, "-b", &bed, "--mode", "count"];
        v.extend(extra);
        run(&v)
    };
    assert!(n(&[]).ends_with("\t1\n") && n(&[]).contains("del\t"));
    assert!(n(&["--split-deletions"]).contains("chr1\t100\t210\tskip\t60\t+\t0\n"));
    assert!(n(&["--split-deletions"]).contains("chr1\t300\t312\tdel\t60\t+\t0\n"));
}

#[test]
fn split_uses_the_cut_at_the_end_of_the_reference() {
    // pastend.bam: `over` covers 994..1004 of a 1000-base chr1 and is cut at 1000.
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", "chr1\t990\t1000\nchr1\t1000\t1010\n");
    let out = run(&[
        "intersect",
        "-a",
        &bed,
        "-b",
        &aln("pastend.bam"),
        "--mode",
        "pairs",
        "--split",
    ]);
    assert_eq!(out, "chr1\t990\t1000\tchr1\t994\t1000\tover\t60\t+\n");
}

#[test]
fn region_restricts_an_alignment_a() {
    let dir = tempdir();
    let bed = write(dir.path(), "t.bed", TARGETS);
    let all = run(&[
        "intersect",
        "-a",
        &aln("order.bam"),
        "-b",
        &bed,
        "-O",
        "sam",
    ]);
    for a in [aln("order.bam"), aln("order.sam")] {
        let got = run(&[
            "intersect",
            "-a",
            &a,
            "-b",
            &bed,
            "--region",
            "chr1",
            "-O",
            "sam",
        ]);
        let want: Vec<&str> = all
            .lines()
            .filter(|l| l.split('\t').nth(2) == Some("chr1"))
            .collect();
        let body: Vec<&str> = got.lines().filter(|l| !l.starts_with('@')).collect();
        assert_eq!(body, want, "{a}");
        assert_eq!(body.len(), 3);
    }
    // Coverage with an alignment B: A restricted to a region.
    let out = run(&[
        "coverage",
        "-a",
        &bed,
        "-b",
        &aln("order.bam"),
        "--region",
        "chr2",
    ]);
    assert!(out.lines().all(|l| l.starts_with("chr2\t")) && !out.is_empty());
}
