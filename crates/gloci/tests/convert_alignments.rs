//! `gloci convert` of BAM and SAM to BED (D106): usage errors, tags, filters, the index and
//! determinism. Byte parity with bedtools is in the golden cases `bamtobed-*`.
mod common;
use common::*;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;
use std::path::PathBuf;

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

#[test]
fn color_needs_bed12() {
    gloci()
        .args(["convert", &aln("cigars.bam"), "--color", "1,2,3"])
        .assert()
        .code(2)
        .stderr(contains("--color"))
        .stderr(contains("--bed12"));
}

#[test]
fn color_must_be_rgb() {
    gloci()
        .args([
            "convert",
            &aln("cigars.bam"),
            "--bed12",
            "--color",
            "1,2,300",
        ])
        .assert()
        .code(2)
        .stderr(contains("R,G,B"));
}

#[test]
fn cigar_is_bed6_only() {
    for other in ["--bed12", "--blocks", "--split-deletions"] {
        gloci()
            .args(["convert", &aln("cigars.bam"), "--cigar", other])
            .assert()
            .code(2)
            .stderr(contains("--cigar"));
    }
}

#[test]
fn alignment_flags_need_alignment_input() {
    let dir = tempdir();
    let f = plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc1\t5\t.\tA\tG\t.\t.\t.\n",
    );
    for flag in [vec!["--bed12"], vec!["--score-tag", "NM"], vec!["--cigar"]] {
        gloci()
            .arg("convert")
            .arg(&f)
            .args(&flag)
            .assert()
            .code(2)
            .stderr(contains("BAM/CRAM/SAM"));
    }
}

#[test]
fn bad_tag_name_is_a_usage_error() {
    gloci()
        .args(["convert", &aln("cigars.bam"), "--score-tag", "N"])
        .assert()
        .code(2)
        .stderr(contains("two characters"));
}

#[test]
fn non_integer_tag_is_a_data_error() {
    // t_floats (reverse) has XF:f:0.1.
    gloci()
        .args([
            "convert",
            &aln("tags.bam"),
            "--score-tag",
            "XF",
            "--require-flags",
            "REVERSE",
        ])
        .assert()
        .code(1)
        .stderr(contains("type f, not an integer"))
        .stderr(contains("t_floats"));
}

#[test]
fn absent_tag_names_the_read_and_keeps_earlier_rows() {
    // AS is on the first read of tags.bam only: its row (negative score) comes out first.
    let o = gloci()
        .args(["convert", &aln("tags.bam"), "--score-tag", "AS"])
        .assert()
        .code(1)
        .stderr(contains("read t_floats has no AS tag"));
    assert_eq!(
        String::from_utf8(o.get_output().stdout.clone()).unwrap(),
        "chrM\t100\t110\tt_ints\t-5\t+\n"
    );
}

#[test]
fn sam_error_names_the_line() {
    gloci()
        .args(["convert", &aln("cigars.sam"), "--score-tag", "AS"])
        .assert()
        .code(1)
        .stderr(contains("cigars.sam line "));
}

#[test]
fn mates_get_suffixes() {
    let dir = tempdir();
    let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
both\t193\tchr1\t11\t5\t4M\t*\t0\t0\tACGT\tIIII\n\
none\t0\tchr1\t21\t5\t4M\t*\t0\t0\tACGT\tIIII\n\
m1\t65\tchr1\t31\t5\t4M\t*\t0\t0\tACGT\tIIII\n";
    let f = plain_file(dir.path(), "m.sam", sam);
    let out = run(&["convert", f.to_str().unwrap()]);
    assert_eq!(
        out,
        "chr1\t10\t14\tboth/1/2\t5\t+\nchr1\t20\t24\tnone\t5\t+\nchr1\t30\t34\tm1/1\t5\t+\n"
    );
}

#[test]
fn counts_unmapped_and_filtered_at_verbose() {
    gloci()
        .args(["-v", "convert", &aln("flags.bam"), "--min-mapq", "30"])
        .assert()
        .success()
        .stderr(contains("unmapped reads left out"))
        .stderr(contains("dropped by --min-mapq"));
}

#[test]
fn write_index_on_a_sorted_bam() {
    let dir = tempdir();
    let o = dir.path().join("o.bed.gz");
    gloci()
        .args([
            "convert",
            &aln("flags.bam"),
            "-o",
            o.to_str().unwrap(),
            "--write-index",
        ])
        .assert()
        .success();
    assert!(dir.path().join("o.bed.gz.tbi").exists());
}

#[test]
fn output_is_the_same_at_every_thread_count() {
    for file in ["cigars.bam", "order.bam", "flags.sam"] {
        for extra in [
            vec![],
            vec!["--bed12"],
            vec!["--blocks", "--split-deletions"],
        ] {
            let f = aln(file);
            let at = |t: &str| {
                let mut a = vec!["convert", &f, "-t", t];
                a.extend(&extra);
                run(&a)
            };
            let one = at("1");
            assert!(!one.is_empty());
            for t in ["2", "8"] {
                assert_eq!(at(t), one, "{file} {extra:?} -t {t}");
            }
        }
    }
}

#[test]
fn split_rows_of_a_sorted_bam_do_not_warn_about_order() {
    // c_skip's second block (705) comes after the next read's start (700): not an input problem.
    gloci()
        .args(["convert", &aln("cigars.bam"), "--blocks"])
        .assert()
        .success()
        .stderr(predicates::str::contains("out of order").not());
}

#[test]
fn bed12_closes_at_the_last_kept_block() {
    let out = run(&["convert", &aln("splice-end.bam"), "--bed12"]);
    assert!(
        out.contains("chr1\t990\t995\tspliced\t60\t+\t990\t995\t255,0,0\t1\t5\t0\n"),
        "{out}"
    );
}

#[test]
fn bed12_opens_at_the_first_kept_block() {
    let out = run(&["convert", &aln("cigars-ends.bam"), "--bed12"]);
    let want = [
        "chr1\t105\t115\tlead\t60\t+\t105\t115\t255,0,0\t1\t10\t0",
        "chr1\t200\t213\tleadD\t60\t+\t200\t213\t255,0,0\t1\t13\t0",
        "chr1\t305\t313\tsc\t60\t+\t305\t313\t255,0,0\t1\t8\t0",
        "chr1\t400\t410\ttailN\t60\t+\t400\t410\t255,0,0\t1\t10\t0",
        "chr1\t500\t513\ttailD\t60\t-\t500\t513\t255,0,0\t1\t13\t0",
    ];
    assert_eq!(out.lines().collect::<Vec<_>>(), want);
    let out = run(&[
        "convert",
        &aln("cigars-ends.bam"),
        "--bed12",
        "--split-deletions",
    ]);
    assert!(
        out.contains("chr1\t203\t213\tleadD\t60\t+\t203\t213\t255,0,0\t1\t10\t0\n")
            && out.contains("chr1\t500\t510\ttailD\t60\t-\t500\t510\t255,0,0\t1\t10\t0\n"),
        "{out}"
    );
}

#[test]
fn bed12_round_trips_through_the_split_commands() {
    // gloci reads what convert --bed12 writes: its blocks are convert --blocks' blocks, and
    // intersect --split takes the BED12 row (it refused a first block past offset 0).
    let dir = tempdir();
    let targets = dir.path().join("t.bed");
    std::fs::write(&targets, "chr1\t0\t100000\n").unwrap();
    let targets = targets.to_string_lossy().into_owned();
    let b12 = dir.path().join("x.bed");
    let b12 = b12.to_string_lossy().into_owned();
    for file in ["cigars-ends.bam", "cigars-bases.bam", "splice-end.bam"] {
        let path = aln(file);
        for extra in [&[][..], &["--split-deletions"][..]] {
            let mut args = vec!["convert", path.as_str()];
            args.extend(extra);
            let mut bed12 = args.clone();
            bed12.extend(["--bed12", "-o", &b12]);
            run(&bed12);
            let mut blocks = args.clone();
            blocks.push("--blocks");
            assert_eq!(
                run(&blocks),
                run(&["convert", &b12, "--blocks"]),
                "{file} {extra:?}"
            );
            let counts = run(&[
                "intersect",
                "-a",
                &b12,
                "-b",
                &targets,
                "--split",
                "--mode",
                "count",
            ]);
            assert!(
                counts.lines().all(|l| l.ends_with("\t1")),
                "{file} {extra:?}: {counts}"
            );
        }
    }
}

#[test]
fn negative_scores_stay_signed() {
    let out = run(&["convert", &aln("negtag.bam"), "--score-tag", "AS"]);
    assert_eq!(
        out,
        "chr1\t10\t20\tneg_one\t-1\t+\nchr1\t30\t40\tneg_five\t-5\t+\n"
    );
}
