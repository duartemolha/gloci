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
    // The indexed reader has no line numbers: it names the record by its
    // number within the chromosome, not by a BGZF offset.
    let stderr = String::from_utf8_lossy(&cmd.stderr);
    assert!(stderr.contains("record 1 of chr2"), "{stderr}");
    assert!(!stderr.contains("BGZF offset"), "{stderr}");

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
    // M4: it names the index, with no made-up line number.
    assert!(stderr.contains("b.bed.gz.tbi"), "{stderr}");
    assert!(!stderr.contains("line 1"), "{stderr}");
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

/// M5: a write error names the real destination, in indexed (chunked) and
/// streaming mode alike, not the placeholder "output".
#[cfg(target_os = "linux")]
#[test]
fn write_errors_name_the_destination() {
    let fx = Fixture::new();
    let (a, b) = (s(&fx.a_indexed), s(&fx.b_indexed));
    for cmd in ["intersect", "subtract", "closest", "map", "coverage"] {
        for t in ["1", "2"] {
            let mut args = vec![cmd, "-a", &a, "-b", &b, "-t", t, "-o", "/dev/full"];
            if cmd == "map" {
                args.extend(["-c", "4:collapse"]);
            }
            let o = gloci().args(&args).assert().code(1);
            let e = String::from_utf8(o.get_output().stderr.clone()).unwrap();
            assert!(e.contains("error writing /dev/full"), "{cmd} -t {t}: {e}");
        }
    }
}

/// Run `args` at `-t 1` and `-t 4` expecting failure; both must print the
/// same message, which must contain every `needles` entry.
fn assert_refuse_alike(args: &[&str], needles: &[&str]) {
    let run = |t: &str| {
        let mut a: Vec<&str> = args.to_vec();
        a.extend(["-t", t]);
        let o = gloci().args(&a).assert().failure().code(1);
        String::from_utf8_lossy(&o.get_output().stderr).into_owned()
    };
    let (one, four) = (run("1"), run("4"));
    assert_eq!(one, four, "-t 1 and -t 4 must refuse alike");
    for n in needles {
        assert!(one.contains(n), "missing {n:?} in: {one}");
    }
}

/// Indexed and streaming relations refuse an input that is not in table
/// order alike, whatever `--threads` is (a lex-sorted indexed BED; a VCF
/// whose header lists its contigs in another order than its records).
#[test]
fn indexed_and_streaming_refuse_alike() {
    let dir = tempdir();
    let bed = |c: &str| format!("{c}\t10\t20\n");
    let lex: String = ["chr1", "chr10", "chr2"].iter().map(|c| bed(c)).collect();
    let (a, _) = indexed_pair(dir.path(), "lexa", &lex);
    let (b, _) = indexed_pair(dir.path(), "lexb", &lex);
    assert_refuse_alike(
        &["intersect", "--a", &a, "--b", &b],
        &["lexb.bed.gz is not sorted", "chr10", "chr2", "sort-bed"],
    );

    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chr2,length=1000>\n##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tC\t.\t.\t.\nchr2\t10\t.\tA\tC\t.\t.\t.\n";
    let (va, _) = indexed_vcf_pair(dir.path(), "hdr", vcf);
    assert_refuse_alike(
        &["intersect", "--a", &va, "--b", &va],
        &["hdr.vcf.gz is not sorted"],
    );

    // The same VCF as an indexed BCF (records in file order chr1, chr2; the
    // dictionary lists chr2 first).
    // `gloci index` refuses such a file, so its CSI is built the way
    // bcftools would.
    let bytes = gloci_testkit::fixtures::vcf_to_bcf(vcf.as_bytes());
    let bcf = gloci_testkit::fixtures::write(dir.path(), "hdr.bcf", &bytes);
    gloci_testkit::fixtures::write(
        dir.path(),
        "hdr.bcf.csi",
        &gloci_testkit::fixtures::bcf_csi_any_order(&bytes),
    );
    assert_refuse_alike(
        &["intersect", "--a", &s(&bcf), "--b", &s(&bcf)],
        &["hdr.bcf is not sorted"],
    );
}

/// No false positive: a file in karyotype order is accepted under the
/// karyotype table, and indexed and streaming runs agree.
#[test]
fn karyotype_sorted_indexed_file_is_accepted() {
    let dir = tempdir();
    let bed = |c: &str| format!("{c}\t10\t20\n");
    let text: String = ["chr1", "chr2", "chr10", "chrX"]
        .iter()
        .map(|c| bed(c))
        .collect();
    let a = indexed_pair(dir.path(), "kara", &text);
    let b = indexed_pair(dir.path(), "karb", &text);
    assert_indexed_equals_streaming(
        "karyotype",
        &["intersect", "--a", "{A}", "--b", "{B}"],
        &a,
        &b,
    );
    let got = out(&["intersect", "--a", &a.0, "--b", &b.0, "-t", "4"]);
    assert_eq!(String::from_utf8_lossy(&got), text);
}

// ---- plain seekable files (D76) ----

/// Run `args` with `-vv` at `-t` and say whether it ran per chromosome.
fn ran_per_chromosome(args: &[&str], t: &str) -> bool {
    let mut a: Vec<&str> = args.to_vec();
    a.extend(["-t", t, "-vv"]);
    let o = gloci().args(&a).assert().success();
    String::from_utf8_lossy(&o.get_output().stderr).contains("running per chromosome")
}

/// A plain BED with a track line, comments and blank lines (also between
/// chromosomes and at the end), Windows line ends on some lines, and 8
/// chromosomes in karyotype order.
fn plain_text(seed: u64, with_noise: bool) -> String {
    let mut t = String::new();
    if with_noise {
        t.push_str("track name=x\n# a comment\n\n");
    }
    let mut state = seed;
    let mut next = |m: u64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) % m
    };
    for c in [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr10", "chrX", "chrY",
    ] {
        let mut starts: Vec<u64> = (0..40).map(|_| next(5000)).collect();
        starts.sort_unstable();
        for (i, st) in starts.iter().enumerate() {
            let len = 1 + next(120);
            let strand = ['+', '-', '.'][next(3) as usize];
            let eol = if with_noise && i % 7 == 3 {
                "\r\n"
            } else {
                "\n"
            };
            t.push_str(&format!(
                "{c}\t{st}\t{}\t{c}_{i}\t{}\t{strand}{eol}",
                st + len,
                next(900)
            ));
            if with_noise && i == 20 {
                t.push_str("\n# mid-file comment\n");
            }
        }
        if with_noise {
            t.push_str("\n#between chromosomes\n");
        }
    }
    t
}

const RELATIONS: &[&[&str]] = &[
    &["intersect"],
    &["intersect", "--mode", "pairs"],
    &["intersect", "--mode", "count"],
    &["subtract"],
    &["closest"],
    &["map", "-c", "4:collapse"],
    &["coverage"],
    &["coverage", "--hist"],
    &["intersect", "-s"],
    &["closest", "-k", "2"],
];

/// Plain BED A and B give, for every relation command, the bytes `-t 1`
/// streaming gives, at `-t 2` and `-t 8`, and run per chromosome.
#[test]
fn seekable_equals_streaming() {
    let dir = tempdir();
    for noise in [false, true] {
        let a = s(&plain_file(dir.path(), "a.bed", &plain_text(1, noise)));
        let b = s(&plain_file(dir.path(), "b.bed", &plain_text(2, noise)));
        for rel in RELATIONS {
            let mut args: Vec<&str> = vec![rel[0], "-a", &a, "-b", &b];
            args.extend(&rel[1..]);
            let mut one = args.clone();
            one.extend(["-t", "1"]);
            let want = out(&one);
            assert!(!want.is_empty(), "{rel:?}");
            for t in ["2", "8"] {
                let mut many = args.clone();
                many.extend(["-t", t]);
                assert_eq!(
                    String::from_utf8_lossy(&out(&many)),
                    String::from_utf8_lossy(&want),
                    "{rel:?} -t {t} noise={noise}"
                );
            }
            assert!(!ran_per_chromosome(&args, "1"), "{rel:?} streams at -t 1");
            assert!(ran_per_chromosome(&args, "2"), "{rel:?} -t 2");
        }
    }
}

/// Several B files, a B with chromosomes A lacks and the reverse, and
/// `--header` and `--names`-style labels behave as streaming does.
#[test]
fn seekable_several_b_and_missing_chromosomes() {
    let dir = tempdir();
    let a = s(&plain_file(dir.path(), "a.bed", &plain_text(3, true)));
    let only = |keep: &[&str], seed| -> String {
        plain_text(seed, false)
            .lines()
            .filter(|l| keep.iter().any(|k| l.split('\t').next() == Some(*k)))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let b1 = s(&plain_file(
        dir.path(),
        "b1.bed",
        &only(&["chr2", "chr10", "chrY"], 4),
    ));
    let b2 = s(&plain_file(
        dir.path(),
        "b2.bed",
        &only(&["chr1", "chr3", "chrX"], 5),
    ));
    for rel in [
        &["intersect", "--mode", "pairs"][..],
        &["coverage"],
        &["closest"],
        &["map", "-c", "4:count"],
    ] {
        let mut args: Vec<&str> = vec![rel[0], "-a", &a, "-b", &b1, "-b", &b2];
        args.extend(&rel[1..]);
        let mut one = args.clone();
        one.extend(["-t", "1"]);
        let want = out(&one);
        let mut four = args.clone();
        four.extend(["-t", "4"]);
        assert_eq!(out(&four), want, "{rel:?}");
        assert!(ran_per_chromosome(&args, "4"), "{rel:?}");
    }
}

/// An indexed input and a plain one, either way round.
#[test]
fn seekable_mixed_with_indexed() {
    let dir = tempdir();
    let (a_text, b_text) = (plain_text(6, false), plain_text(7, false));
    let (a_idx, a_plain) = indexed_pair(dir.path(), "ma", &a_text);
    let b_plain = s(&plain_file(dir.path(), "mb.bed", &b_text));
    let (b_idx, _) = indexed_pair(dir.path(), "mb", &b_text);
    let a_plain_uncompressed = s(&plain_file(dir.path(), "ma.bed", &a_text));
    let _ = a_plain;
    for rel in RELATIONS {
        for (a, b) in [(&a_idx, &b_plain), (&a_plain_uncompressed, &b_idx)] {
            let mut args: Vec<&str> = vec![rel[0], "-a", a, "-b", b];
            args.extend(&rel[1..]);
            let mut one = args.clone();
            one.extend(["-t", "1"]);
            let want = out(&one);
            let mut many = args.clone();
            many.extend(["-t", "4"]);
            assert_eq!(out(&many), want, "{rel:?} a={a}");
            assert!(ran_per_chromosome(&args, "4"), "{rel:?} a={a}");
        }
    }
}

/// Plain VCF A and B.
#[test]
fn seekable_vcf() {
    let dir = tempdir();
    let vcf = |seed: u64| -> String {
        let mut t = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n##contig=<ID=chr3,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        );
        let mut x = seed;
        for c in ["chr1", "chr2", "chr3"] {
            let mut ps: Vec<u64> = (0..30)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                    1 + (x >> 33) % 600
                })
                .collect();
            ps.sort_unstable();
            for p in ps {
                t.push_str(&format!("{c}\t{p}\t.\tAC\tA\t.\t.\t.\n"));
            }
        }
        t
    };
    let a = s(&plain_file(dir.path(), "a.vcf", &vcf(1)));
    let b = s(&plain_file(dir.path(), "b.vcf", &vcf(2)));
    for rel in [
        &["intersect"][..],
        &["intersect", "--mode", "pairs"],
        &["coverage"],
    ] {
        let mut args: Vec<&str> = vec![rel[0], "-a", &a, "-b", &b];
        args.extend(&rel[1..]);
        let mut one = args.clone();
        one.extend(["-t", "1"]);
        let want = out(&one);
        assert!(!want.is_empty(), "{rel:?}");
        let mut many = args.clone();
        many.extend(["-t", "4"]);
        assert_eq!(out(&many), want, "{rel:?}");
        assert!(ran_per_chromosome(&args, "4"), "{rel:?}");
    }
}

/// stdin, `/dev/stdin`, fd paths, gzip and GFF are never split: they stream.
#[test]
fn seekable_excludes_streams_and_other_formats() {
    let dir = tempdir();
    let text = plain_text(8, false);
    let a = s(&plain_file(dir.path(), "a.bed", &text));
    let b = s(&plain_file(dir.path(), "b.bed", &text));
    let want = out(&["intersect", "-a", &a, "-b", &b, "-t", "1"]);
    for path in ["-", "/dev/stdin"] {
        let o = gloci()
            .args(["intersect", "-a", path, "-b", &b, "-t", "4", "-vv"])
            .write_stdin(text.clone())
            .assert()
            .success();
        assert_eq!(o.get_output().stdout, want, "{path}");
        let err = String::from_utf8_lossy(&o.get_output().stderr);
        assert!(!err.contains("running per chromosome"), "{path}: {err}");
    }
    // gzip (not BGZF): streams
    let gz = s(&gloci_testkit::fixtures::write(
        dir.path(),
        "a.bed.gz",
        &gloci_testkit::fixtures::gzip(text.as_bytes()),
    ));
    assert!(!ran_per_chromosome(
        &["intersect", "-a", &gz, "-b", &b],
        "4"
    ));
    // GFF stays streaming (a ##FASTA section can follow its records)
    let gff = "##gff-version 3\nchr1\tsrc\tgene\t11\t20\t.\t+\t.\tID=g1\nchr2\tsrc\tgene\t11\t20\t.\t+\t.\tID=g2\n";
    let g = s(&plain_file(dir.path(), "a.gff", gff));
    assert!(!ran_per_chromosome(&["intersect", "-a", &g, "-b", &g], "4"));
    // `--unsorted` and `--sort` stream
    assert!(!ran_per_chromosome(
        &["intersect", "-a", &a, "-b", &b, "--unsorted"],
        "4"
    ));
}

/// Non-contiguous chromosomes, or an order the table rejects, are refused
/// as streaming refuses them, at every thread count.
#[test]
fn seekable_order_errors_are_streamings() {
    let dir = tempdir();
    let good = "chr1\t10\t20\nchr2\t10\t20\nchr3\t10\t20\n";
    let g = s(&plain_file(dir.path(), "good.bed", good));
    let cases: &[(&str, &str, &[&str])] = &[
        (
            "noncontig",
            "chr1\t10\t20\nchr2\t10\t20\nchr1\t30\t40\n",
            &["noncontig.bed is not sorted", "line 3"],
        ),
        (
            "backwards",
            "chr2\t10\t20\nchr1\t10\t20\n",
            &["backwards.bed is not sorted", "line 2"],
        ),
        (
            "lex",
            "chr1\t10\t20\nchr10\t10\t20\nchr2\t10\t20\n",
            &["lex.bed is not sorted", "chr10", "chr2"],
        ),
    ];
    for (name, text, needles) in cases {
        let bad = s(&plain_file(dir.path(), &format!("{name}.bed"), text));
        let o = gloci()
            .args(["intersect", "-a", &bad, "-b", &g, "-t", "4", "-vv"])
            .assert()
            .failure();
        let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
        assert!(!err.contains("running per chromosome"), "{name}: {err}");
        assert_refuse_alike(&["intersect", "-a", &bad, "-b", &g], needles);
        assert_refuse_alike(&["intersect", "-a", &g, "-b", &bad], needles);
        assert_refuse_alike(&["closest", "-a", &bad, "-b", &bad], needles);
    }
    // An order the table rejects only under -g: chr2 before chr1 in a genome
    // that lists chr1 first; and a chromosome the closed table lacks.
    let unknown = s(&plain_file(
        dir.path(),
        "unknown.bed",
        "chr1\t10\t20\nchrNope\t10\t20\n",
    ));
    assert_refuse_alike(
        &["intersect", "-a", &unknown, "-b", &g, "-g", "hg38"],
        &["chrNope"],
    );
    assert_refuse_alike(
        &["intersect", "-a", &g, "-b", &unknown, "-g", "hg38"],
        &["chrNope"],
    );
}

/// A within-chromosome error carries streaming's line number and text; when
/// several chromosomes fail, the one earliest in the file is reported.
#[test]
fn seekable_errors_inside_a_chromosome_match_streaming() {
    let dir = tempdir();
    let mut rows: Vec<String> = Vec::new();
    for c in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6"] {
        for i in 0..50 {
            rows.push(format!("{c}\t{}\t{}\n", i * 10, i * 10 + 5));
        }
    }
    let b = s(&plain_file(dir.path(), "b.bed", &rows.concat()));
    // A, with a blank and a comment line to move the line numbers.
    let with = |edits: &[(usize, &str)]| -> String {
        let mut r = rows.clone();
        for (i, t) in edits {
            r[*i] = t.to_string();
        }
        r.insert(10, "\n".to_string());
        r.insert(30, "# note\n".to_string());
        r.concat()
    };
    let cases: Vec<(&str, String, &[&str])> = vec![
        (
            "badnum",
            with(&[(120, "chr3\tNOPE\t5\n")]),
            &["badnum.bed", "start"],
        ),
        (
            "order",
            with(&[(210, "chr5\t3\t9\n")]),
            &["order.bed is not sorted"],
        ),
        (
            "two",
            with(&[(70, "chr2\t1\n"), (250, "chr6\tx\t5\n")]),
            &["two.bed", "expected at least 3"],
        ),
        ("cols", with(&[(160, "chr4 5 9\n")]), &["cols.bed"]),
    ];
    for (name, text, needles) in &cases {
        let a = s(&plain_file(dir.path(), &format!("{name}.bed"), text));
        for cmd in ["intersect", "coverage", "closest"] {
            assert_refuse_alike(&[cmd, "-a", &a, "-b", &b], needles);
            assert_refuse_alike(&[cmd, "-a", &b, "-b", &a], needles);
        }
        // The error came from a worker, not from a fallback to streaming.
        let o = gloci()
            .args(["intersect", "-a", &a, "-b", &b, "-t", "4", "-vv"])
            .assert()
            .failure();
        let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
        // (A line with no tab is a chromosome of its own, in the middle of
        // another: that is not contiguous, so streaming reports it.)
        assert_eq!(
            err.contains("running per chromosome"),
            *name != "cols",
            "{name}: {err}"
        );
    }
    // --skip-invalid: the same warning, with the same line numbers.
    let a = s(&plain_file(
        dir.path(),
        "skip.bed",
        &with(&[
            (70, "chr2\t1\n"),
            (130, "chr3\tx\t5\n"),
            (250, "chr6\tx\t5\n"),
        ]),
    ));
    let run = |t: &str| {
        let o = gloci()
            .args(["intersect", "-a", &a, "-b", &b, "--skip-invalid", "-t", t])
            .assert()
            .success();
        (
            o.get_output().stdout.clone(),
            String::from_utf8_lossy(&o.get_output().stderr).into_owned(),
        )
    };
    let (one, four) = (run("1"), run("4"));
    assert_eq!(one, four);
    assert!(one.1.contains("skipped 3 invalid records"), "{}", one.1);
}

/// `--chrom-map` renames chromosomes on plain inputs read per chromosome:
/// the output is streaming's, and the scan's contiguity check uses the
/// file's own names.
#[test]
fn seekable_chrom_map() {
    let dir = tempdir();
    let text = |names: &[&str]| -> String {
        names
            .iter()
            .flat_map(|c| (0..5).map(move |i| format!("{c}\t{}\t{}\tr{i}\n", i * 20, i * 20 + 30)))
            .collect()
    };
    let a = s(&plain_file(
        dir.path(),
        "a.bed",
        &text(&["1", "2", "3", "10", "X"]),
    ));
    let b = s(&plain_file(
        dir.path(),
        "b.bed",
        &text(&["chr1", "chr2", "chr10", "chrX"]),
    ));
    let map = s(&plain_file(
        dir.path(),
        "map.tsv",
        "1\tchr1\n2\tchr2\n3\tchr3\n10\tchr10\nX\tchrX\n",
    ));
    for rel in ["intersect", "coverage", "closest"] {
        let args = [rel, "-a", &a, "-b", &b, "--chrom-map", &map];
        let want = out(&[&args[..], &["-t", "1"]].concat());
        assert!(!want.is_empty());
        assert_eq!(out(&[&args[..], &["-t", "4"]].concat()), want, "{rel}");
        assert!(ran_per_chromosome(&args, "4"), "{rel}");
    }
    // Two file names renamed to one (`1` to `chr2`, where the file also has
    // `chr2`): not split, streamed as before.
    let clash_in = s(&plain_file(
        dir.path(),
        "clash.bed",
        "1\t0\t10\nchr2\t5\t15\nchr3\t0\t10\n",
    ));
    let clash = s(&plain_file(dir.path(), "clash.tsv", "1\tchr2\n"));
    let args = [
        "intersect",
        "-a",
        &clash_in,
        "-b",
        &clash_in,
        "--chrom-map",
        &clash,
    ];
    assert!(!ran_per_chromosome(&args, "4"));
    assert_eq!(
        out(&[&args[..], &["-t", "4"]].concat()),
        out(&[&args[..], &["-t", "1"]].concat())
    );
    // Names that are in file order only after renaming: streaming's verdict.
    let swap = s(&plain_file(dir.path(), "swap.tsv", "1\tchr9\n2\tchr1\n"));
    assert_refuse_alike(
        &["intersect", "-a", &a, "-b", &a, "--chrom-map", &swap],
        &["is not sorted"],
    );
}

/// `--region` on plain A, read per chromosome, selects what streaming does.
#[test]
fn seekable_with_regions() {
    let dir = tempdir();
    let a = s(&plain_file(dir.path(), "a.bed", &plain_text(9, true)));
    let b = s(&plain_file(dir.path(), "b.bed", &plain_text(10, false)));
    for rel in ["intersect", "closest", "coverage"] {
        let args = [
            rel,
            "-a",
            &a,
            "-b",
            &b,
            "--region",
            "chr3:100-900",
            "--region",
            "chrX:1-2000",
            "--region",
            "chr10:50-60",
        ];
        let want = out(&[&args[..], &["-t", "1"]].concat());
        assert!(!want.is_empty());
        assert_eq!(out(&[&args[..], &["-t", "4"]].concat()), want, "{rel}");
        assert!(ran_per_chromosome(&args, "4"), "{rel}");
    }
}

/// Records of unknown strand under `-s` are reported once per input, not
/// once per chunk, plain or indexed, as at `-t 1`.
#[test]
fn unknown_strand_warning_is_once_per_input() {
    let dir = tempdir();
    let text: String = ["chr1", "chr2", "chr3", "chr4"]
        .iter()
        .map(|c| format!("{c}\t1\t10\ta\t0\t.\n{c}\t20\t30\tb\t0\t+\n"))
        .collect();
    let a = s(&plain_file(dir.path(), "a.bed", &text));
    let (a_idx, _) = indexed_pair(dir.path(), "ia", &text);
    for a in [&a, &a_idx] {
        for t in ["1", "4"] {
            let o = gloci()
                .args(["intersect", "-a", a, "-b", a, "-s", "-t", t])
                .assert()
                .success();
            let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
            // A and B are the same file: one warning for each.
            assert_eq!(
                err.matches("unknown strand").count(),
                2,
                "{a} -t {t}: {err}"
            );
        }
    }
}
