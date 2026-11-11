//! The alignment/FASTA/PSL input gate (D100) and the unwritable CRAM output.
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::str::contains;

#[test]
fn input_gate() {
    let dir = tempdir();
    let bam = fixtures::write(
        dir.path(),
        "x.bam",
        &fixtures::bgzip(b"BAM\x01\0\0\0\0\0\0\0\0"),
    );
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n");
    let fa = plain_file(dir.path(), "ref.fa", ">chr1 d\nACGT\n");
    let psl = plain_file(
        dir.path(),
        "x.psl",
        "50\t1\t0\t0\t0\t0\t0\t0\t+\tq1\t51\t0\t51\tchr1\t1000\t100\t151\t1\t51,\t0,\t100,\n",
    );
    let bed = plain_file(dir.path(), "a.bed", "chr1\t1\t2\n");

    let x = bam.to_str().unwrap();
    gloci()
        .args(["merge", x])
        .assert()
        .code(2)
        .stderr(contains("does not read alignments"))
        .stderr(contains(format!("gloci convert {x} -o x.bed")))
        .stderr(contains("gloci sort -m"));
    gloci()
        .args(["merge", sam.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("does not read alignments"));
    gloci()
        .args(["sort", fa.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("getfasta, maskfasta or nuc"));
    let p = psl.to_str().unwrap();
    gloci()
        .args(["merge", p])
        .assert()
        .code(2)
        .stderr(contains(format!("gloci convert {p}")));
    gloci()
        .args(["sort", bed.to_str().unwrap(), "-o"])
        .arg(dir.path().join("o.cram"))
        .assert()
        .code(2)
        .stderr(contains("samtools view -C"));
}

#[test]
fn convert_reads_alignments() {
    // `convert` opts in to alignments: a header-only SAM has no reads and writes nothing.
    let dir = tempdir();
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n");
    gloci()
        .args(["convert", sam.to_str().unwrap()])
        .assert()
        .success()
        .stdout("");
}

#[test]
fn genomecov_reads_alignments() {
    // `genomecov` reads alignments: a SAM with no reads covers its header's chromosomes with zero.
    let dir = tempdir();
    let sam = plain_file(dir.path(), "x.sam", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:5\n");
    gloci()
        .args(["genomecov", sam.to_str().unwrap(), "--bga"])
        .assert()
        .success()
        .stdout("chr1\t0\t5\t0\n");
}

#[test]
fn commands_that_do_not_read_alignments_say_never() {
    let dir = tempdir();
    let bam = fixtures::write(
        dir.path(),
        "x.bam",
        &fixtures::bgzip(b"BAM\x01\0\0\0\0\0\0\0\0"),
    );
    let bed = plain_file(dir.path(), "a.bed", "chr1\t1\t2\n");
    let x = bam.to_str().unwrap();
    let b = bed.to_str().unwrap();
    gloci()
        .args(["subtract", "-a", b, "-b", x])
        .assert()
        .code(2)
        .stderr(contains("does not read alignments"));
}

#[test]
fn an_empty_input_is_empty_where_records_are_needed() {
    // One message for an input with no bytes after decompression (M1).
    let dir = tempdir();
    let empty_bgzf = fixtures::bgzip(b"");
    let paths: Vec<String> = [
        ("x.sam", Vec::new()),
        ("x.bam", Vec::new()),
        ("e.bam", empty_bgzf.clone()),
        ("x.fa", Vec::new()),
        ("e.fa.gz", empty_bgzf),
    ]
    .into_iter()
    .map(|(n, b)| {
        fixtures::write(dir.path(), n, &b)
            .to_string_lossy()
            .into_owned()
    })
    .collect();
    let bed = plain_file(dir.path(), "m.bed", "a\t0\t2\n");
    let bed = bed.to_str().unwrap();
    let out = dir.path().join("o.bam");
    let out = out.to_str().unwrap();
    let good = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln/cat-a.bam");
    let good = good.to_str().unwrap();
    let mut runs: Vec<(String, Vec<&str>)> = Vec::new();
    for (i, p) in paths.iter().enumerate() {
        runs.push((p.clone(), vec!["convert", p]));
        // index refuses only a file of zero bytes; an EOF-only BGZF file is an empty interval
        // file and gets an empty index (see an_empty_bgzf_file_gets_an_empty_index).
        if ![2, 4].contains(&i) {
            runs.push((p.clone(), vec!["index", p]));
        }
        runs.push((p.clone(), vec!["maskfasta", "--fasta", p, "-b", bed]));
    }
    for p in &paths[1..3] {
        runs.push((p.clone(), vec!["cat", good, p, "-o", out]));
    }
    for (p, args) in &runs {
        gloci()
            .args(args)
            .assert()
            .code(1)
            .stderr(contains(format!("{p} is empty: it has no records")));
    }
    // Commands whose result may be empty read it as no records; so does convert --blocks, for
    // which an empty file is an empty BED.
    for args in [
        vec!["convert", paths[0].as_str(), "--blocks"],
        vec!["view", paths[0].as_str()],
        vec!["sort", paths[0].as_str()],
        vec!["intersect", "-a", paths[0].as_str(), "-b", bed],
    ] {
        gloci().args(&args).assert().success();
    }
}

#[test]
fn an_empty_bgzf_file_gets_an_empty_index() {
    // An empty result written as .bed.gz is the BGZF EOF block alone. gloci index treats it as an
    // empty interval file, as tabix -p bed does and as --write-index on the same output does.
    let dir = tempdir();
    let gz = fixtures::write(dir.path(), "e.bed.gz", &fixtures::bgzip(b""));
    let gz = gz.to_str().unwrap();
    gloci().args(["index", gz]).assert().success();
    assert!(dir.path().join("e.bed.gz.tbi").exists());
    gloci()
        .args(["view", gz, "chr1"])
        .assert()
        .success()
        .stdout("");
}
