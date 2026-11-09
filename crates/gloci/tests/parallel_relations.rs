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
        "GLOCI_TEST_CHUNK_MEM_LIMIT",
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
        .env("GLOCI_TEST_CHUNK_MEM_LIMIT", "1024")
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

/// Indexed and plain copies of `text` in `dir` (`NAME.bed.gz` with a
/// `.tbi`, and `plain/NAME.bed.gz` without one).
fn indexed_pair(dir: &Path, name: &str, text: &str) -> (String, String) {
    let plain_dir = dir.join("plain");
    fs::create_dir_all(&plain_dir).unwrap();
    let idx = bgzf_file(dir, &format!("{name}.bed.gz"), text);
    gloci().args(["index", &s(&idx)]).assert().success();
    let plain = bgzf_file(&plain_dir, &format!("{name}.bed.gz"), text);
    (s(&idx), s(&plain))
}

/// Run `args` (with `{A}`/`{B}` placeholders) streamed at `-t 1` on the
/// plain copies, and per chromosome at `-t 2` and `-t 8` on the indexed
/// ones; every output must be byte-identical.
fn assert_indexed_equals_streaming(
    name: &str,
    template: &[&str],
    a: &(String, String),
    b: &(String, String),
) {
    let fill = |idx: bool| -> Vec<String> {
        template
            .iter()
            .map(|t| match *t {
                "{A}" => {
                    if idx {
                        a.0.clone()
                    } else {
                        a.1.clone()
                    }
                }
                "{B}" => {
                    if idx {
                        b.0.clone()
                    } else {
                        b.1.clone()
                    }
                }
                o => o.to_string(),
            })
            .collect()
    };
    let mut base = fill(false);
    base.extend(["-t".into(), "1".into()]);
    let base_refs: Vec<&str> = base.iter().map(String::as_str).collect();
    let want = out(&base_refs);
    for t in ["2", "8"] {
        let mut args = fill(true);
        args.extend(["-t".into(), t.into()]);
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let got = out(&refs);
        assert_eq!(
            String::from_utf8_lossy(&got),
            String::from_utf8_lossy(&want),
            "{name} -t {t}: indexed output differs from streaming"
        );
    }
}

/// C3: a zero-length record at position 0 is read in indexed mode (whole
/// chromosome reads apply no footprint filter).
#[test]
fn point_at_zero_indexed_equals_streaming() {
    let dir = tempdir();
    let a = indexed_pair(
        dir.path(),
        "a",
        "chr1\t0\t0\ta0\nchr1\t5\t10\ta1\nchr2\t3\t8\ta2\n",
    );
    let b = indexed_pair(
        dir.path(),
        "b",
        "chr1\t0\t0\tb0\nchr1\t6\t7\tb1\nchr2\t4\t5\tb2\n",
    );
    // B without the point: A's point survives subtract.
    let b2 = indexed_pair(dir.path(), "b2", "chr1\t6\t7\tb1\nchr2\t4\t5\tb2\n");
    let cases: &[(&str, &[&str])] = &[
        ("subtract", &["subtract", "-a", "{A}", "-b", "{B}"]),
        ("coverage", &["coverage", "-a", "{A}", "-b", "{B}"]),
        (
            "map",
            &["map", "-a", "{A}", "-b", "{B}", "-c", "4:collapse"],
        ),
        (
            "closest",
            &["closest", "-a", "{A}", "-b", "{B}", "--distance"],
        ),
        (
            "pairs",
            &["intersect", "-a", "{A}", "-b", "{B}", "--mode", "pairs"],
        ),
    ];
    for (name, template) in cases {
        assert_indexed_equals_streaming(name, template, &a, &b);
        assert_indexed_equals_streaming(name, template, &a, &b2);
    }
    // The point really is in the output.
    let got = out(&[
        "intersect",
        "-a",
        &a.0,
        "-b",
        &b.0,
        "--mode",
        "pairs",
        "-t",
        "2",
    ]);
    assert!(
        String::from_utf8_lossy(&got).starts_with("chr1\t0\t0\ta0\tchr1\t0\t0\tb0\n"),
        "{}",
        String::from_utf8_lossy(&got)
    );
}

/// Indexed and plain copies of a VCF.
fn indexed_vcf_pair(dir: &Path, name: &str, text: &str) -> (String, String) {
    let plain_dir = dir.join("plain");
    fs::create_dir_all(&plain_dir).unwrap();
    let idx = bgzf_file(dir, &format!("{name}.vcf.gz"), text);
    gloci().args(["index", &s(&idx)]).assert().success();
    let plain = bgzf_file(&plain_dir, &format!("{name}.vcf.gz"), text);
    (s(&idx), s(&plain))
}

/// C2: chromosomes with no B get B file 1's real width in their
/// placeholder rows, in indexed mode as in streaming mode.
#[test]
fn missing_b_chromosomes_indexed_equals_streaming() {
    let dir = tempdir();
    let mut a_text = String::new();
    for c in ["chr1", "chr2", "chr3", "chr4"] {
        a_text.push_str(&format!("{c}\t100\t200\ta_{c}\t1\t+\n"));
        a_text.push_str(&format!("{c}\t500\t600\tb_{c}\t2\t-\n"));
    }
    let a = indexed_pair(dir.path(), "a", &a_text);
    // B only on chr1 and chr3, six columns.
    let b = indexed_pair(
        dir.path(),
        "b",
        "chr1\t150\t160\tb1\t9\t+\nchr3\t700\t800\tb3\t9\t-\n",
    );
    let b2 = indexed_pair(dir.path(), "b2", "chr3\t120\t130\tc3\t8\t+\n");
    let cases: &[(&str, &[&str])] = &[
        (
            "closest --distance",
            &["closest", "-a", "{A}", "-b", "{B}", "--distance"],
        ),
        (
            "closest --each-side",
            &["closest", "-a", "{A}", "-b", "{B}", "--each-side"],
        ),
        (
            "pairs --keep-unmatched",
            &[
                "intersect",
                "-a",
                "{A}",
                "-b",
                "{B}",
                "--mode",
                "pairs",
                "--keep-unmatched",
            ],
        ),
    ];
    for (name, template) in cases {
        assert_indexed_equals_streaming(name, template, &a, &b);
    }
    // --per-file over two B files (B file 1 sets the width).
    let fill = |idx: bool, t: &str| -> Vec<String> {
        let pick = |p: &(String, String)| if idx { p.0.clone() } else { p.1.clone() };
        vec![
            "closest".into(),
            "-a".into(),
            pick(&a),
            "-b".into(),
            pick(&b),
            pick(&b2),
            "--per-file".into(),
            "--distance".into(),
            "-t".into(),
            t.into(),
        ]
    };
    let base = fill(false, "1");
    let want = out(&base.iter().map(String::as_str).collect::<Vec<_>>());
    for t in ["2", "8"] {
        let args = fill(true, t);
        let got = out(&args.iter().map(String::as_str).collect::<Vec<_>>());
        assert_eq!(
            String::from_utf8_lossy(&got),
            String::from_utf8_lossy(&want),
            "closest --per-file -t {t}"
        );
    }
    // Every placeholder has B's six columns.
    let got = out(&[
        "intersect",
        "-a",
        &a.0,
        "-b",
        &b.0,
        "--mode",
        "pairs",
        "--keep-unmatched",
        "-t",
        "2",
    ]);
    let got = String::from_utf8(got).unwrap();
    assert!(
        got.contains("chr2\t100\t200\ta_chr2\t1\t+\t.\t-1\t-1\t.\t-1\t.\n"),
        "{got}"
    );
    for line in got.lines() {
        assert_eq!(line.split('\t').count(), 12, "ragged row: {line}");
    }

    // A VCF A: the #CHROM line names six B columns and every row has them.
    let vcf = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=10000>\n\
##contig=<ID=chr2,length=10000>\n##contig=<ID=chr3,length=10000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t155\tv1\tA\tG\t.\tPASS\t.\nchr2\t155\tv2\tA\tG\t.\tPASS\t.\n\
chr3\t755\tv3\tA\tG\t.\tPASS\t.\n";
    let va = indexed_vcf_pair(dir.path(), "a", vcf);
    let vcf_cases: &[(&str, &[&str])] = &[
        (
            "vcf closest",
            &["closest", "-a", "{A}", "-b", "{B}", "--distance"],
        ),
        (
            "vcf pairs --keep-unmatched",
            &[
                "intersect",
                "-a",
                "{A}",
                "-b",
                "{B}",
                "--mode",
                "pairs",
                "--keep-unmatched",
            ],
        ),
    ];
    for (name, template) in vcf_cases {
        assert_indexed_equals_streaming(name, template, &va, &b);
    }
    let got = out(&[
        "intersect",
        "-a",
        &va.0,
        "-b",
        &b.0,
        "--mode",
        "pairs",
        "--keep-unmatched",
        "-t",
        "2",
    ]);
    let got = String::from_utf8(got).unwrap();
    let header = got.lines().next().unwrap();
    assert!(header.ends_with("\tb_1\tb_2\tb_3\tb_4\tb_5\tb_6"), "{got}");
    let width = header.split('\t').count();
    for line in got.lines() {
        assert_eq!(line.split('\t').count(), width, "ragged row: {line}");
    }
}

/// Equally near B are reported in B's input order (D38) in indexed mode
/// too: region readers number records in file order.
#[test]
fn ties_keep_b_input_order_indexed() {
    let dir = tempdir();
    let (mut a_text, mut b_text) = (String::new(), String::new());
    for c in 1..=40 {
        let chrom = format!("chr{c}");
        a_text.push_str(&format!(
            "{chrom}\t2051\t2052\ta{c}\n{chrom}\t5000\t5100\tx{c}\n"
        ));
        // Two B ending at the same base: the same distance from A.
        b_text.push_str(&format!(
            "{chrom}\t1803\t2016\tfirst{c}\n{chrom}\t1857\t2016\tsecond{c}\n\
{chrom}\t5050\t5060\tp{c}\n{chrom}\t5050\t5060\tq{c}\n"
        ));
    }
    let a = indexed_pair(dir.path(), "a", &a_text);
    let b = indexed_pair(dir.path(), "b", &b_text);
    let cases: &[(&str, &[&str])] = &[
        (
            "closest",
            &["closest", "-a", "{A}", "-b", "{B}", "--distance"],
        ),
        (
            "closest -k 3",
            &["closest", "-a", "{A}", "-b", "{B}", "-k", "3"],
        ),
        (
            "pairs",
            &["intersect", "-a", "{A}", "-b", "{B}", "--mode", "pairs"],
        ),
        (
            "map",
            &["map", "-a", "{A}", "-b", "{B}", "-c", "4:collapse"],
        ),
    ];
    for (name, template) in cases {
        assert_indexed_equals_streaming(name, template, &a, &b);
    }
}
