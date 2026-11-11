//! `gloci convert --bedpe` (D156): pairs from name-sorted alignments. Byte parity with bedtools
//! is in the golden cases `bamtobed-bedpe-*`.
mod common;
use common::*;
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

const PAIRS: &str = "\
chr1\t100\t150\tchr1\t300\t350\tr1\t60\t+\t-
.\t-1\t-1\tchr1\t100\t150\tr2\t0\t.\t+
.\t-1\t-1\t.\t-1\t-1\tr3\t0\t.\t.
chr1\t1000\t1050\tchr1\t2000\t2050\tr4\t30\t+\t-
chr10\t200\t250\tchr2\t500\t550\tr5\t40\t-\t+
chr1\t900\t950\tchr1\t1000\t1050\tr6\t60\t+\t-
";

#[test]
fn bedpe_pairs_table() {
    assert_eq!(run(&["convert", &aln("pairs.sam"), "--bedpe"]), PAIRS);
}

#[test]
fn sam_bam_and_cram_agree() {
    for f in ["pairs.sam", "pairs.bam", "pairs.noref.cram"] {
        assert_eq!(run(&["convert", &aln(f), "--bedpe"]), PAIRS, "{f}");
    }
}

#[test]
fn mate1_first_keeps_read1_first() {
    let out = run(&["convert", &aln("pairs.bam"), "--bedpe", "--mate1-first"]);
    let rows: Vec<&str> = out.lines().collect();
    assert_eq!(rows[0], "chr1\t100\t150\tchr1\t300\t350\tr1\t60\t+\t-");
    // READ2 is the unmapped one: it stays second.
    assert_eq!(rows[1], "chr1\t100\t150\t.\t-1\t-1\tr2\t0\t+\t.");
    assert_eq!(rows[4], "chr2\t500\t550\tchr10\t200\t250\tr5\t40\t+\t-");
    assert_eq!(rows[5], "chr1\t1000\t1050\tchr1\t900\t950\tr6\t60\t-\t+");
}

#[test]
fn score_tag_sums_the_mapped_ends() {
    let out = run(&["convert", &aln("pairs.bam"), "--bedpe", "--score-tag", "NM"]);
    let scores: Vec<&str> = out.lines().map(|l| l.split('\t').nth(7).unwrap()).collect();
    assert_eq!(scores, ["3", "1", "0", "3", "1", "6"]);
}

#[test]
fn score_tag_missing_on_a_mapped_end_is_an_error() {
    gloci()
        .args(["convert", &aln("pairs.bam"), "--bedpe", "--score-tag", "AS"])
        .assert()
        .code(1)
        .stderr(contains("read r1 has no AS tag"))
        .stderr(contains("record 1 of chr1 (record 1 of the file)"));
    gloci()
        .args(["convert", &aln("pairs.sam"), "--bedpe", "--score-tag", "AS"])
        .assert()
        .code(1)
        .stderr(contains("read r1 has no AS tag"))
        .stderr(contains("line 5"));
}

#[test]
fn supplementary_reads_are_not_mates() {
    assert_eq!(
        run(&["convert", &aln("supp.sam"), "--bedpe"]),
        "chr1\t100\t150\tchr1\t300\t350\ts1\t60\t+\t-\n"
    );
    assert_eq!(
        run(&["convert", &aln("supp.bam"), "--bedpe"]),
        "chr1\t100\t150\tchr1\t300\t350\ts1\t60\t+\t-\n"
    );
}

#[test]
fn orphan_mate_warns_with_a_count() {
    let dir = tempdir();
    let sam = plain_file(
        dir.path(),
        "orphan.sam",
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n\
a\t65\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
b\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    );
    gloci()
        .args(["convert"])
        .arg(&sam)
        .arg("--bedpe")
        .assert()
        .success()
        .stdout("")
        .stderr(contains("1 read name has one mate"))
        .stderr(contains("1 read is neither READ1 nor READ2"));
}

#[test]
fn coordinate_sorted_input_is_refused() {
    gloci()
        .args(["convert", &aln("flags.bam"), "--bedpe"])
        .assert()
        .code(2)
        .stderr(contains("SO:coordinate"))
        .stderr(contains("gloci sort -n"));
}

#[test]
fn conflicts_and_usage() {
    for other in ["--bed12", "--blocks", "--cigar", "--split-deletions"] {
        gloci()
            .args(["convert", &aln("pairs.bam"), "--bedpe", other])
            .assert()
            .code(2);
    }
    gloci()
        .args(["convert", &aln("pairs.bam"), "--mate1-first"])
        .assert()
        .code(2)
        .stderr(contains("--bedpe"));
    gloci()
        .args(["convert", &aln("pairs.bam"), "--bedpe", "-O", "bed"])
        .assert()
        .code(2)
        .stderr(contains(".bedpe"));
    let dir = tempdir();
    let bed = plain_file(dir.path(), "x.bed", "chr1\t1\t2\n");
    gloci()
        .args(["convert"])
        .arg(&bed)
        .arg("--bedpe")
        .assert()
        .failure();
}

#[test]
fn write_index_is_refused_with_the_gloci_index_hint() {
    let dir = tempdir();
    let out = dir.path().join("p.bedpe.gz");
    gloci()
        .args([
            "convert",
            &aln("pairs.bam"),
            "--bedpe",
            "--write-index",
            "-o",
        ])
        .arg(&out)
        .assert()
        .code(2)
        .stderr(contains("gloci index -s 1 -b 2 -e 3 --zero-based"));
}

#[test]
fn threads_do_not_change_the_output() {
    for t in ["1", "2", "8"] {
        assert_eq!(
            run(&["convert", &aln("pairs.bam"), "--bedpe", "-t", t]),
            PAIRS
        );
    }
}

#[test]
fn pair_intersect_reads_the_bedpe_it_writes() {
    let dir = tempdir();
    let pairs = dir.path().join("pairs.bedpe.gz");
    gloci()
        .args(["convert", &aln("pairs.bam"), "--bedpe", "-o"])
        .arg(&pairs)
        .assert()
        .success();
    let b = plain_file(
        dir.path(),
        "b.bed",
        "chr1\t90\t110\tgeneA\nchr10\t0\t1000\tgeneB\n",
    );
    let out = run(&[
        "pair-intersect",
        "-a",
        pairs.to_str().unwrap(),
        "-b",
        b.to_str().unwrap(),
    ]);
    // r1 end 1 and r2 (its mapped end 2) hit geneA; r5's chr10 end hits geneB; pairs with an
    // absent end are read without trouble.
    assert!(out.contains("r1\t60\t+\t-\tchr1\t90\t110\tgeneA"), "{out}");
    assert!(out.contains("r2\t0\t.\t+\tchr1\t90\t110\tgeneA"), "{out}");
    assert!(out.contains("r5\t40\t-\t+\tchr10\t0\t1000\tgeneB"), "{out}");
    assert!(!out.contains("r3"), "{out}");
}

#[test]
fn past_end_reads_follow_convert() {
    for f in ["pairs-pastend.sam", "pairs-pastend.bam"] {
        // e1's READ1 straddles the end of chr1 (1,000 bases) and is cut; e2's READ1 starts at
        // the end, so it is an absent end: not in the score, no NM needed.
        let o = gloci()
            .args(["convert", &aln(f), "--bedpe"])
            .assert()
            .success()
            .stderr(contains(
                "read e1 on chr1 extends past the end of the reference (1000 bases)",
            ))
            .stderr(contains(
                "read e2 starts at 1001 on chr1, which is 1000 bases long",
            ))
            .stderr(contains("1 reads left out for starting at or past the end"));
        assert_eq!(
            String::from_utf8(o.get_output().stdout.clone()).unwrap(),
            "chr1\t100\t150\tchr1\t990\t1000\te1\t60\t-\t+\n.\t-1\t-1\tchr1\t300\t350\te2\t0\t.\t-\n",
            "{f}"
        );
        let o = gloci()
            .args(["convert", &aln(f), "--bedpe", "--score-tag", "NM"])
            .assert()
            .success();
        let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
        let scores: Vec<&str> = out.lines().map(|l| l.split('\t').nth(7).unwrap()).collect();
        assert_eq!(scores, ["3", "3"], "{f}");
    }
}

#[test]
fn duplicate_primary_mates_warn_and_use_the_first() {
    let dir = tempdir();
    let sam = plain_file(
        dir.path(),
        "dup.sam",
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n\
a\t65\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
a\t65\tchr1\t500\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
a\t129\tchr1\t50\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    );
    gloci()
        .args(["convert"])
        .arg(&sam)
        .arg("--bedpe")
        .assert()
        .success()
        .stdout("chr1\t4\t14\tchr1\t49\t59\ta\t60\t+\t+\n")
        .stderr(contains(
            "1 read name has more than one primary READ1 or READ2; the first is used",
        ));
}
