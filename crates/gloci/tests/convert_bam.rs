//! `gloci convert --bam` (D167): features and BEDPE pairs as BAM records. Byte parity with
//! bedtools is in the golden cases `bedtobam-*` and `bedpetobam-*`; here are the files, genomes,
//! thread counts, errors and round trips those cannot show.
mod common;
use common::*;
use gloci_testkit::differential::compare;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn tobam(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/tobam")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn golden(case: &str) -> Vec<u8> {
    std::fs::read(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/golden")
            .join(case)
            .join("expected"),
    )
    .unwrap()
}

fn stdout(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

/// A BAM as `gloci view -O sam` shows it, without the `@PG` lines.
fn sam(bam: &[u8]) -> Vec<String> {
    let out = gloci()
        .args(["view", "-", "-O", "sam"])
        .write_stdin(bam.to_vec())
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(String::from)
        .collect()
}

fn records(bam: &[u8]) -> Vec<String> {
    sam(bam)
        .into_iter()
        .filter(|l| !l.starts_with('@'))
        .collect()
}

fn no_temp_files(dir: &Path) {
    let left: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.contains(".tmp-"))
        .collect();
    assert!(left.is_empty(), "temporary files left: {left:?}");
}

#[test]
fn a_file_holds_the_golden_bam() {
    let d = tempdir();
    let out = d.path().join("out.bam");
    gloci()
        .args(["convert", &tobam("features.bed"), "--bam", "-g"])
        .arg(tobam("genome.txt"))
        .arg("-o")
        .arg(&out)
        .assert()
        .success()
        .stdout("");
    // The golden's AS: names `genome.txt`; this run names the full path.
    let want = golden("bedtobam-vs-bedtools");
    let got = std::fs::read(&out).unwrap();
    assert_eq!(records(&got), records(&want));
    let out = d.path().join("pairs.bam");
    gloci()
        .args([
            "convert",
            &tobam("pairs.bedpe"),
            "--format",
            "bedpe",
            "--bam",
        ])
        .args(["-g", &tobam("genome-lex.txt"), "-o"])
        .arg(&out)
        .assert()
        .success();
    let got = std::fs::read(&out).unwrap();
    assert_eq!(records(&got), records(&golden("bedpetobam-vs-bedtools")));
    no_temp_files(d.path());
}

#[test]
fn the_header_lists_the_genome_in_its_order_and_reads_back() {
    let g = tobam("genome.txt");
    let bam = stdout(&["convert", &tobam("features.bed"), "--bam", "-g", &g]);
    let lines = sam(&bam);
    let header: Vec<&str> = lines
        .iter()
        .filter(|l| l.starts_with('@'))
        .map(String::as_str)
        .collect();
    let sq = |c: &str, n: u32| format!("@SQ\tSN:{c}\tAS:{g}\tLN:{n}");
    assert_eq!(
        header,
        [
            "@HD\tVN:1.0\tSO:unsorted".to_string(),
            sq("chr2", 50000),
            sq("chr1", 250000),
            sq("chr10", 3000),
            sq("chrX", 2000),
            sq("chrM", 16571),
        ]
    );
    let recs = records(&bam);
    assert_eq!(recs.len(), 10);
    assert_eq!(recs[0], "r1\t0\tchr1\t101\t255\t100M\t*\t0\t0\t*\t*");
    assert_eq!(recs[1], "r2\t16\tchr2\t11\t255\t10M\t*\t0\t0\t*\t*");
    // A zero-length record is one base wider on each side.
    assert_eq!(recs[4], "zero\t16\tchr1\t100\t255\t2M\t*\t0\t0\t*\t*");
}

#[test]
fn every_thread_count_writes_the_same_bytes() {
    let g = tobam("genome.txt");
    let gl = tobam("genome-lex.txt");
    let runs: [Vec<String>; 3] = [
        vec![
            "convert".into(),
            tobam("features.bed"),
            "--bam".into(),
            "-g".into(),
            g.clone(),
        ],
        vec![
            "convert".into(),
            tobam("bed12.bed"),
            "--bam".into(),
            "--bed12".into(),
            "-g".into(),
            g.clone(),
        ],
        vec![
            "convert".into(),
            tobam("pairs.bedpe"),
            "--format".into(),
            "bedpe".into(),
            "--bam".into(),
            "-g".into(),
            gl,
        ],
    ];
    let d = tempdir();
    for args in &runs {
        let mut outs = Vec::new();
        let mut files = Vec::new();
        for t in ["-t1", "-t2", "--threads=8"] {
            let o = gloci().args(args).arg(t).output().unwrap();
            assert_eq!(o.status.code(), Some(0), "{args:?} {t}");
            outs.push(o.stdout);
            // One path, so that the command lines in the @PG lines are the same.
            let f = d.path().join("o.bam");
            gloci()
                .args(args)
                .arg(t)
                .arg("-o")
                .arg(&f)
                .assert()
                .success();
            files.push(std::fs::read(&f).unwrap());
        }
        assert!(outs.iter().all(|o| *o == outs[0]), "{args:?}: stdout");
        assert!(files.iter().all(|o| *o == files[0]), "{args:?}: files");
        // The file's @PG line names -o; the rest is the same.
        assert!(
            compare("bam", &outs[0], &files[0]),
            "{args:?}: stdout and file"
        );
    }
}

#[test]
fn a_built_in_genome_or_a_fai_gives_the_references() {
    let bam = stdout(&["convert", &tobam("bed4.bed"), "--bam", "-g", "hg38"]);
    let lines = sam(&bam);
    assert_eq!(lines[1], "@SQ\tSN:chr1\tAS:hg38\tLN:248956422");
    assert!(lines.contains(&"@SQ\tSN:chrX\tAS:hg38\tLN:156040895".to_string()));
    assert_eq!(records(&bam)[0], "a\t0\tchr2\t1\t255\t100M\t*\t0\t0\t*\t*");
    let d = tempdir();
    let fai = plain_file(
        d.path(),
        "ref.fa.fai",
        "chr2\t5000\t6\t60\t61\nchr1\t300000\t5096\t60\t61\nchrX\t2000\t310000\t60\t61\n",
    );
    let fai = fai.to_str().unwrap();
    let bam = stdout(&["convert", &tobam("bed4.bed"), "--bam", "-g", fai]);
    let lines = sam(&bam);
    assert_eq!(lines[1], format!("@SQ\tSN:chr2\tAS:{fai}\tLN:5000"));
    assert_eq!(lines[3], format!("@SQ\tSN:chrX\tAS:{fai}\tLN:2000"));
    assert_eq!(records(&bam).len(), 4);
}

#[test]
fn chrom_map_renames_before_the_genome_is_asked() {
    let d = tempdir();
    let bed = plain_file(d.path(), "a.bed", "1\t10\t20\tr\t0\t-\n");
    let map = plain_file(d.path(), "map.tsv", "1\tchr1\n");
    gloci()
        .args([
            "convert",
            bed.to_str().unwrap(),
            "--bam",
            "-g",
            &tobam("genome.txt"),
        ])
        .assert()
        .failure()
        .code(1)
        .stderr(contains("chromosome \"1\" is not in the genome file"));
    let bam = stdout(&[
        "convert",
        bed.to_str().unwrap(),
        "--bam",
        "-g",
        &tobam("genome.txt"),
        "--chrom-map",
        map.to_str().unwrap(),
    ]);
    assert_eq!(records(&bam), ["r\t16\tchr1\t11\t255\t10M\t*\t0\t0\t*\t*"]);
}

#[test]
fn a_chromosome_the_genome_lacks_is_exit_1_and_writes_no_file() {
    let d = tempdir();
    for (input, extra) in [
        (tobam("missing-chrom.bed"), vec![]),
        (
            tobam("pairs-missing-chrom.bedpe"),
            vec!["--format", "bedpe"],
        ),
    ] {
        let out = d.path().join("out.bam");
        gloci()
            .args(["convert", &input, "--bam", "-g", &tobam("genome.txt")])
            .args(&extra)
            .arg("-o")
            .arg(&out)
            .assert()
            .failure()
            .code(1)
            .stderr(contains(
                "line 2: chromosome \"chrZ\" is not in the genome file",
            ))
            .stderr(contains("add chrZ to the genome"))
            .stderr(contains("--chrom-map"));
        assert!(!out.exists(), "{input}: a partial output was left");
        no_temp_files(d.path());
    }
}

#[test]
fn absent_ends_are_unmapped_mates_that_come_back_as_absent() {
    let bam = stdout(&[
        "convert",
        &tobam("pairs-absent.bedpe"),
        "--format",
        "bedpe",
        "--bam",
        "-g",
        &tobam("genome-lex.txt"),
    ]);
    assert_eq!(
        records(&bam),
        [
            "a2\t73\tchr1\t101\t255\t100M\t=\t101\t0\t*\t*",
            "a2\t133\tchr1\t101\t0\t*\t=\t101\t0\t*\t*",
            "a1\t101\tchr2\t11\t0\t*\t=\t11\t0\t*\t*",
            "a1\t153\tchr2\t11\t255\t10M\t=\t11\t0\t*\t*",
            "both\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*",
            "both\t141\t*\t0\t0\t*\t*\t0\t0\t*\t*",
            // The absent end's own `-` strand: 0x10 on it, 0x20 on its mate.
            "a2m\t105\tchr1\t301\t255\t100M\t=\t301\t0\t*\t*",
            "a2m\t149\tchr1\t301\t0\t*\t=\t301\t0\t*\t*",
        ]
    );
    // `convert --bedpe` reads an unmapped mate as an absent end again.
    let d = tempdir();
    let f = d.path().join("pairs.bam");
    std::fs::write(&f, &bam).unwrap();
    let pairs = String::from_utf8(stdout(&["convert", f.to_str().unwrap(), "--bedpe"])).unwrap();
    let ends: Vec<(&str, usize)> = pairs
        .lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let absent = [&c[0..3], &c[3..6]]
                .iter()
                .filter(|e| **e == [".", "-1", "-1"])
                .count();
            (c[6], absent)
        })
        .collect();
    assert_eq!(ends, [("a2", 1), ("a1", 1), ("both", 2), ("a2m", 1)]);
}

#[test]
fn bed12_reads_blocks_only_from_bed12() {
    gloci()
        .args(["convert", &tobam("features.gff"), "--bam", "--bed12"])
        .args(["-g", &tobam("genome.txt"), "-o", "/dev/null"])
        .assert()
        .failure()
        .code(1)
        .stderr(contains(
            "--bed12 needs BED12 records with 12 columns; this input is GFF3",
        ));
    gloci()
        .args(["convert", &tobam("bed12-unsorted.bed"), "--bam", "--bed12"])
        .args(["-g", &tobam("genome.txt"), "-o", "/dev/null"])
        .assert()
        .failure()
        .code(1)
        .stderr(contains(
            "line 1: invalid BED12: the first block must start at offset 0",
        ));
    let d = tempdir();
    let overlap = plain_file(
        d.path(),
        "o.bed",
        "chr1\t100\t300\tovl\t0\t+\t100\t300\t0\t3\t100,100,50,\t0,50,150,\n",
    );
    gloci()
        .args(["convert", overlap.to_str().unwrap(), "--bam", "--bed12"])
        .args(["-g", &tobam("genome.txt"), "-o", "/dev/null"])
        .assert()
        .failure()
        .code(1)
        .stderr(contains(
            "line 1: invalid BED12: blocks are unsorted or overlap",
        ));
    gloci()
        .args(["convert", &tobam("long-name.bed"), "--bam"])
        .args(["-g", &tobam("genome.txt"), "-o", "/dev/null"])
        .assert()
        .failure()
        .code(1)
        .stderr(contains(
            "cannot be written as BAM: the read name is 255 characters long",
        ));
}

#[test]
fn usage_errors_are_exit_2() {
    let g = tobam("genome.txt");
    let bed = tobam("features.bed");
    let d = tempdir();
    let bedpe_named = d.path().join("pairs.bedpe");
    std::fs::copy(tobam("pairs.bedpe"), &bedpe_named).unwrap();
    let bedpe_named = bedpe_named.to_str().unwrap();
    let bam_in = d.path().join("in.bam");
    std::fs::write(&bam_in, golden("bedtobam-vs-bedtools")).unwrap();
    let bam_in = bam_in.to_str().unwrap();
    let cases: Vec<(Vec<&str>, &str)> = vec![
        (
            vec!["convert", &bed, "--bam", "-g", &g, "-O", "sam"],
            "convert --bam writes BAM, not SAM",
        ),
        (
            vec!["convert", &bed, "--bam", "-g", &g, "-o", "x.bed"],
            "convert --bam writes BAM, not BED",
        ),
        (
            vec!["convert", bedpe_named, "--bam", "-g", &g],
            "add --format bedpe",
        ),
        (
            vec![
                "convert",
                bedpe_named,
                "--format",
                "bedpe",
                "--bam",
                "--bed12",
                "-g",
                &g,
            ],
            "BEDPE has none",
        ),
        (vec!["convert", &bed, "--bam"], "--genome <NAME|FILE>"),
        (vec!["convert", &bed, "-g", &g], "--bam"),
        (vec!["convert", &bed, "--mapq", "7"], "--bam"),
        (
            vec!["convert", &bed, "--bam", "-g", &g, "--mapq", "256"],
            "256",
        ),
        (
            vec!["convert", &bed, "--bam", "-g", &g, "--fastq"],
            "cannot be used with",
        ),
        (
            vec!["convert", &bed, "--bam", "-g", &g, "--bedpe"],
            "cannot be used with",
        ),
        (vec!["convert", bam_in, "--bam", "-g", &g], "BAM"),
        (
            vec!["convert", &bed, "--bam", "-g", "no-such-genome-file.txt"],
            "-g",
        ),
    ];
    for (args, msg) in cases {
        gloci()
            .current_dir(d.path())
            .args(&args)
            .assert()
            .failure()
            .code(2)
            .stderr(contains(msg));
    }
    // A genome without lengths cannot make a BAM header.
    let dict = plain_file(d.path(), "n.dict", "@SQ\tSN:chr1\n");
    gloci()
        .args(["convert", &bed, "--bam", "-g", dict.to_str().unwrap()])
        .assert()
        .failure()
        .code(2)
        .stderr(contains("the genome gives no length for chr1"));
    assert!(!d.path().join("x.bed").exists());
}

#[test]
fn write_index_needs_input_in_the_genome_order() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.txt", "chr1\t1000\nchr2\t1000\n");
    let g = g.to_str().unwrap();
    let sorted = plain_file(
        d.path(),
        "s.bed",
        "chr1\t1\t2\tr\t0\t+\nchr1\t50\t60\ts\t0\t-\nchr2\t5\t9\tt\t0\t+\n",
    );
    let out = d.path().join("s.bam");
    gloci()
        .args([
            "convert",
            sorted.to_str().unwrap(),
            "--bam",
            "-g",
            g,
            "--write-index",
            "-o",
        ])
        .arg(&out)
        .assert()
        .success();
    assert!(d.path().join("s.bam.bai").exists());
    let hit = stdout(&["view", out.to_str().unwrap(), "chr1:40-70", "-O", "sam"]);
    let hit: Vec<&str> = std::str::from_utf8(&hit)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .collect();
    assert_eq!(hit, ["s\t16\tchr1\t51\t255\t10M\t*\t0\t0\t*\t*"]);
    let unsorted = plain_file(
        d.path(),
        "u.bed",
        "chr2\t5\t9\tt\t0\t+\nchr1\t1\t2\tr\t0\t+\n",
    );
    let out = d.path().join("u.bam");
    gloci()
        .args([
            "convert",
            unsorted.to_str().unwrap(),
            "--bam",
            "-g",
            g,
            "--write-index",
            "-o",
        ])
        .arg(&out)
        .assert()
        .failure()
        .code(1)
        .stderr(contains("gloci sort"));
    assert!(!out.exists());
    no_temp_files(d.path());
}

#[test]
fn level_0_stores_the_blocks_and_reads_back_the_same() {
    let args = [
        "convert",
        &tobam("bed12.bed"),
        "--bam",
        "--bed12",
        "-g",
        &tobam("genome.txt"),
    ];
    let level6 = stdout(&args);
    let level0 = gloci()
        .args(args)
        .args(["--compression-level", "0"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    // The header text is in the file as it is.
    assert!(level0.windows(17).any(|w| w == b"@HD\tVN:1.0\tSO:uns"));
    assert!(!level6.windows(17).any(|w| w == b"@HD\tVN:1.0\tSO:uns"));
    assert!(compare("bam", &level0, &level6));
    assert_eq!(sam(&level0), sam(&level6));
}

#[test]
fn stdin_and_empty_input() {
    let g = tobam("genome.txt");
    let bam = gloci()
        .args(["convert", "-", "--bam", "-g", &g])
        .write_stdin("chr1\t100\t200\tr1\t5\t+\n")
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(
        records(&bam),
        ["r1\t0\tchr1\t101\t255\t100M\t*\t0\t0\t*\t*"]
    );
    // An empty input is a BAM with the header and no records, as bedtools writes it.
    for extra in [&[][..], &["--format", "bedpe"][..]] {
        let bam = stdout(
            &[
                &["convert", &tobam("empty.bed"), "--bam", "-g", &g][..],
                extra,
            ]
            .concat(),
        );
        assert!(records(&bam).is_empty());
        assert_eq!(sam(&bam).len(), 6);
    }
}
