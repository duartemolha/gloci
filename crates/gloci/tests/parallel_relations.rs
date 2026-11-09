mod common;
use common::*;
use predicates::str::contains;
use std::fs;
use std::path::{Path, PathBuf};

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> Vec<u8> {
    let o = gloci().args(args).assert().success();
    o.get_output().stdout.clone()
}

fn out_env(args: &[&str], key: &str, val: &str) -> Vec<u8> {
    let o = gloci().args(args).env(key, val).assert().success();
    o.get_output().stdout.clone()
}

struct Fixture {
    dir: tempfile::TempDir,
    a_indexed: PathBuf,
    b_indexed: PathBuf,
    a_unindexed: PathBuf,
    b_unindexed: PathBuf,
    bad_a_indexed: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let dir = tempdir();

        let mut a_text = String::new();
        let mut b_text = String::new();
        let mut valid_a_text = String::new();
        let mut bad_a_text = String::new();

        for i in 1..=6 {
            let chrom = format!("chr{i}");
            a_text.push_str(&format!("{chrom}\t10\t50\ta_{i}_1\t100\t+\n"));
            a_text.push_str(&format!("{chrom}\t60\t100\ta_{i}_2\t200\t-\n"));
            a_text.push_str(&format!("{chrom}\t120\t200\ta_{i}_3\t300\t+\n"));

            b_text.push_str(&format!("{chrom}\t20\t40\tb_{i}_1\t50\t+\n"));
            b_text.push_str(&format!("{chrom}\t45\t70\tb_{i}_2\t60\t+\n"));
            b_text.push_str(&format!("{chrom}\t150\t180\tb_{i}_3\t70\t-\n"));

            valid_a_text.push_str(&format!("{chrom}\t10\t50\ta_{i}_1\t100\t+\n"));
            if i == 2 {
                bad_a_text.push_str(&format!("{chrom}\tNOT_A_NUMBER\t50\ta_{i}_1\t100\t+\n"));
            } else {
                bad_a_text.push_str(&format!("{chrom}\t10\t50\ta_{i}_1\t100\t+\n"));
            }
        }

        let a_indexed = bgzf_file(dir.path(), "a.bed.gz", &a_text);
        let b_indexed = bgzf_file(dir.path(), "b.bed.gz", &b_text);
        gloci().args(["index", &s(&a_indexed)]).assert().success();
        gloci().args(["index", &s(&b_indexed)]).assert().success();

        let unindexed_dir = dir.path().join("no_index");
        fs::create_dir(&unindexed_dir).unwrap();
        let a_unindexed = bgzf_file(&unindexed_dir, "a.bed.gz", &a_text);
        let b_unindexed = bgzf_file(&unindexed_dir, "b.bed.gz", &b_text);

        // Build valid index on valid file, then replace with corrupted chr2 file
        let bad_dir = dir.path().join("bad");
        fs::create_dir(&bad_dir).unwrap();
        let bad_a_indexed = bgzf_file(&bad_dir, "bad_a.bed.gz", &valid_a_text);
        gloci()
            .args(["index", &s(&bad_a_indexed)])
            .assert()
            .success();
        // Overwrite bad_a.bed.gz with corrupted text, leaving bad_a.bed.gz.tbi intact
        let _ = bgzf_file(&bad_dir, "bad_a.bed.gz", &bad_a_text);

        Self {
            dir,
            a_indexed,
            b_indexed,
            a_unindexed,
            b_unindexed,
            bad_a_indexed,
        }
    }
}

#[test]
fn parallel_equals_streaming() {
    let fx = Fixture::new();

    let a_idx = s(&fx.a_indexed);
    let b_idx = s(&fx.b_indexed);
    let a_un = s(&fx.a_unindexed);
    let b_un = s(&fx.b_unindexed);

    // 5 commands:
    // 1. intersect
    // 2. subtract
    // 3. closest
    // 4. map
    // 5. coverage
    let test_cases: Vec<(&'static str, Vec<&str>)> = vec![
        ("intersect", vec!["intersect", "-a", "{A}", "-b", "{B}"]),
        (
            "intersect_pairs",
            vec!["intersect", "-a", "{A}", "-b", "{B}", "--mode", "pairs"],
        ),
        (
            "intersect_count",
            vec!["intersect", "-a", "{A}", "-b", "{B}", "--mode", "count"],
        ),
        ("subtract", vec!["subtract", "-a", "{A}", "-b", "{B}"]),
        ("closest", vec!["closest", "-a", "{A}", "-b", "{B}"]),
        (
            "map",
            vec!["map", "-a", "{A}", "-b", "{B}", "-c", "4:collapse"],
        ),
        ("coverage", vec!["coverage", "-a", "{A}", "-b", "{B}"]),
        (
            "coverage_hist",
            vec!["coverage", "-a", "{A}", "-b", "{B}", "--hist"],
        ),
    ];

    for (name, template) in test_cases {
        let stream_args: Vec<String> = template
            .iter()
            .map(|tok| match *tok {
                "{A}" => a_un.clone(),
                "{B}" => b_un.clone(),
                other => other.to_string(),
            })
            .collect();
        let mut stream_cmd: Vec<&str> = stream_args.iter().map(|s| s.as_str()).collect();
        stream_cmd.extend_from_slice(&["-t", "1"]);
        let baseline = out(&stream_cmd);

        for threads in ["1", "2", "8"] {
            let indexed_args: Vec<String> = template
                .iter()
                .map(|tok| match *tok {
                    "{A}" => a_idx.clone(),
                    "{B}" => b_idx.clone(),
                    other => other.to_string(),
                })
                .collect();
            let mut indexed_cmd: Vec<&str> = indexed_args.iter().map(|s| s.as_str()).collect();
            indexed_cmd.extend_from_slice(&["-t", threads]);
            let got = out(&indexed_cmd);
            assert_eq!(
                got, baseline,
                "command {name} with -t {threads} does not match baseline streaming output",
            );
        }
    }
}

#[test]
fn chunk_spill_path() {
    let fx = Fixture::new();
    let a_idx = s(&fx.a_indexed);
    let b_idx = s(&fx.b_indexed);
    let a_un = s(&fx.a_unindexed);
    let b_un = s(&fx.b_unindexed);

    let baseline = out(&["intersect", "-a", &a_un, "-b", &b_un, "-t", "1"]);

    let got = out_env(
        &["intersect", "-a", &a_idx, "-b", &b_idx, "-t", "2"],
        "GLOCI_CHUNK_MEM_LIMIT",
        "1024",
    );
    assert_eq!(got, baseline, "chunk spill output does not match baseline");
}

#[test]
fn chunk_error_cleans_up() {
    let fx = Fixture::new();
    let bad_a = s(&fx.bad_a_indexed);
    let b_idx = s(&fx.b_indexed);

    let test_dir = fx.dir.path().join("error_test");
    fs::create_dir(&test_dir).unwrap();
    let out_file = test_dir.join("out.bed");
    let tmp_dir = test_dir.join("tmp");
    fs::create_dir(&tmp_dir).unwrap();

    let cmd = gloci()
        .args([
            "intersect",
            "-a",
            &bad_a,
            "-b",
            &b_idx,
            "-t",
            "2",
            "-o",
            &s(&out_file),
        ])
        .env("GLOCI_TMPDIR", s(&tmp_dir))
        .output()
        .unwrap();

    assert!(!cmd.status.success(), "expected failure on corrupted chr2");
    assert!(!out_file.exists(), "output file should not exist on error");

    let tmp_entries: Vec<_> = fs::read_dir(&tmp_dir).unwrap().collect();
    assert!(
        tmp_entries.is_empty(),
        "tmp_dir not empty after error: {:?}",
        tmp_entries
    );
}

#[test]
fn missing_index_falls_back() {
    let fx = Fixture::new();
    let a_idx = s(&fx.a_indexed);
    let b_un = s(&fx.b_unindexed);
    let a_un = s(&fx.a_unindexed);

    let baseline = out(&["intersect", "-a", &a_un, "-b", &b_un, "-t", "1"]);
    // a has index, b does not -> falls back to streaming silently
    let got = out(&["intersect", "-a", &a_idx, "-b", &b_un, "-t", "2"]);
    assert_eq!(got, baseline);
}

#[test]
fn stale_index_warns_in_relation() {
    let fx = Fixture::new();
    let a_idx = s(&fx.a_indexed);
    let b_idx = s(&fx.b_indexed);

    let tbi = format!("{a_idx}.tbi");
    let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
    std::fs::File::options()
        .write(true)
        .open(&tbi)
        .unwrap()
        .set_modified(old)
        .unwrap();

    gloci()
        .args(["intersect", "-a", &a_idx, "-b", &b_idx, "-t", "2"])
        .assert()
        .success()
        .stderr(contains(format!("the index {tbi} is older than {a_idx}")));
}

#[test]
fn cli_tmp_dir_on_indexed_relations() {
    let fx = Fixture::new();
    let a_idx = s(&fx.a_indexed);
    let b_idx = s(&fx.b_indexed);
    let tmp = fx.dir.path().join("cli_tmp");
    fs::create_dir(&tmp).unwrap();

    gloci()
        .args([
            "intersect",
            "-a",
            &a_idx,
            "-b",
            &b_idx,
            "--tmp-dir",
            &s(&tmp),
            "-t",
            "2",
        ])
        .env("GLOCI_CHUNK_MEM_LIMIT", "1024")
        .assert()
        .success();
}

#[test]
fn closed_table_unknown_b_contig_fails() {
    let dir = tempdir();
    let a_text = "chr1\t10\t20\nchr2\t10\t20\n";
    let b_text = "chr1\t10\t20\nchrUnknown\t10\t20\n";
    let a_file = bgzf_file(dir.path(), "a.bed.gz", a_text);
    let b_file = bgzf_file(dir.path(), "b.bed.gz", b_text);
    gloci().args(["index", &s(&a_file)]).assert().success();
    gloci().args(["index", &s(&b_file)]).assert().success();

    let cmd = gloci()
        .args([
            "intersect",
            "-a",
            &s(&a_file),
            "-b",
            &s(&b_file),
            "-g",
            "hg38",
            "-t",
            "2",
        ])
        .output()
        .unwrap();

    assert!(!cmd.status.success());
    let stderr = String::from_utf8_lossy(&cmd.stderr);
    assert!(
        stderr.contains("b.bed.gz"),
        "stderr should name the B file: {stderr}"
    );
    assert!(
        stderr.contains("chrUnknown"),
        "stderr should name the unknown chromosome: {stderr}"
    );
}
