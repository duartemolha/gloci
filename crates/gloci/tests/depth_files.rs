//! `gloci depth --out-prefix` (D173): the files match the mosdepth 0.3.10 outputs committed in
//! `tests/mosdepth/` (made by its `regen.sh`, in the comparison image), except where
//! `docs/differences.md` says mosdepth is wrong. Those cases pin gloci's own files in
//! `<case>/gloci/` (written under `GLOCI_BLESS=1`, then reviewed and committed).
mod common;
use common::*;
use gloci_testkit::differential::check_anchor;
use predicates::str::contains;
use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/samtools")
}

fn f(name: &str) -> String {
    inputs().join(name).display().to_string()
}

fn case_dir(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/mosdepth")
        .join(case)
}

fn bgzf_text(path: &Path) -> String {
    let raw = std::fs::read(path).unwrap();
    // A BGZF block: gzip magic, the extra-field flag and the `BC` subfield.
    assert!(
        raw.starts_with(&[0x1f, 0x8b, 8, 4]) && raw[12..14] == *b"BC",
        "{} is not BGZF",
        path.display()
    );
    let mut out = String::new();
    noodles_bgzf::io::Reader::new(&raw[..])
        .read_to_string(&mut out)
        .unwrap();
    out
}

/// The case table, shared with `regen.sh`: name, input, gloci's arguments (a `@` prefix is a
/// file of the samtools inputs), and the anchor of a documented difference.
struct Case {
    name: &'static str,
    input: &'static str,
    args: &'static [&'static str],
    anchor: Option<&'static str>,
}

const CASES: &[Case] = &[
    Case {
        name: "depth-window-100",
        input: "depth.bam",
        args: &["--by", "100"],
        anchor: None,
    },
    Case {
        name: "depth-targets-thresholds",
        input: "depth.bam",
        args: &["--by", "@targets.bed", "--thresholds", "1,2,10"],
        anchor: None,
    },
    Case {
        name: "pairs-window-thresholds",
        input: "mosdepth.bam",
        args: &["--by", "100", "--thresholds", "1,2,3"],
        anchor: None,
    },
    Case {
        name: "pairs-bed-mapq",
        input: "mosdepth.bam",
        args: &[
            "--min-mapq",
            "20",
            "--by",
            "@mosdepth.bed",
            "--thresholds",
            "2,1",
        ],
        anchor: None,
    },
    Case {
        name: "pairs-no-by",
        input: "mosdepth.bam",
        args: &[],
        anchor: None,
    },
    Case {
        name: "pairs-flags-zero",
        input: "mosdepth.bam",
        args: &["--exclude-flags", "0", "--by", "500"],
        anchor: None,
    },
    Case {
        name: "means-base-by-base",
        input: "mosdepth-means.bam",
        args: &["--by", "200"],
        anchor: None,
    },
    Case {
        name: "empty",
        input: "empty.bam",
        args: &["--by", "100"],
        anchor: None,
    },
    Case {
        name: "mates-contained",
        input: "mates.bam",
        args: &["--by", "100"],
        anchor: Some("mosdepth-contained-mate"),
    },
    Case {
        name: "stale-pair-entry",
        input: "mosdepth-stale.bam",
        args: &["--exclude-flags", "UNMAP,QCFAIL,DUP", "--by", "100"],
        anchor: Some("mosdepth-mate-table"),
    },
    Case {
        name: "bed-blank-line",
        input: "mosdepth.bam",
        args: &["--by", "@mosdepth-blank.bed"],
        anchor: Some("mosdepth-bed-blank-line"),
    },
];

fn resolve(args: &[&str]) -> Vec<String> {
    args.iter()
        .map(|a| match a.strip_prefix('@') {
            Some(file) => f(file),
            None => a.to_string(),
        })
        .collect()
}

/// The files of `dir` written for prefix `o`, by their mosdepth names without the prefix: text
/// as it is, each `.bed.gz` (checked to be BGZF) as its decompressed `.bed`; and the raw names.
fn written(dir: &Path) -> (BTreeMap<String, String>, Vec<String>) {
    let mut files = BTreeMap::new();
    let mut raw = Vec::new();
    for e in std::fs::read_dir(dir).unwrap() {
        let e = e.unwrap();
        let name = e.file_name().into_string().unwrap();
        let short = name.strip_prefix("o.").expect("every file starts with o.");
        raw.push(short.to_string());
        if short.ends_with(".csi") {
            continue;
        }
        let text = if let Some(bed) = short.strip_suffix(".bed.gz") {
            files.insert(format!("{bed}.bed"), bgzf_text(&e.path()));
            continue;
        } else {
            std::fs::read_to_string(e.path()).unwrap()
        };
        files.insert(short.to_string(), text);
    }
    raw.sort();
    (files, raw)
}

/// The files committed in `dir` (not its `gloci/` directory).
fn committed(dir: &Path) -> BTreeMap<String, String> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap())
        .filter(|e| e.file_type().unwrap().is_file())
        .map(|e| {
            (
                e.file_name().into_string().unwrap(),
                std::fs::read_to_string(e.path()).unwrap(),
            )
        })
        .collect()
}

/// The names `written` reports for a committed set: `regions.bed` is `o.regions.bed.gz` and its
/// index.
fn expected_raw(files: &BTreeMap<String, String>) -> Vec<String> {
    let mut raw = Vec::new();
    for n in files.keys() {
        if n.ends_with(".bed") {
            raw.push(format!("{n}.gz"));
            raw.push(format!("{n}.gz.csi"));
        } else {
            raw.push(n.clone());
        }
    }
    raw.sort();
    raw
}

fn run_case(input: &str, args: &[String]) -> (BTreeMap<String, String>, Vec<String>) {
    let tmp = tempdir();
    gloci()
        .args(["depth", input])
        .args(args)
        .arg("--out-prefix")
        .arg(tmp.path().join("o"))
        .assert()
        .success();
    written(tmp.path())
}

fn bless() -> bool {
    std::env::var_os("GLOCI_BLESS").is_some()
}

#[test]
fn files_match_mosdepth() {
    let doc = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/differences.md"),
    )
    .unwrap();
    for c in CASES {
        let (got, raw) = run_case(&f(c.input), &resolve(c.args));
        let want = committed(&case_dir(c.name));
        match c.anchor {
            None => {
                assert_eq!(raw, expected_raw(&want), "{}: files written", c.name);
                for (name, text) in &want {
                    assert_eq!(&got[name], text, "{}: {name}", c.name);
                }
            }
            Some(a) => {
                check_anchor(&doc, a).unwrap_or_else(|e| panic!("{}: {e}", c.name));
                assert_eq!(raw, expected_raw(&want), "{}: files written", c.name);
                assert_ne!(got, want, "{}: documented, so it must differ", c.name);
                let pinned = case_dir(c.name).join("gloci");
                if bless() {
                    std::fs::create_dir_all(&pinned).unwrap();
                    for (name, text) in &got {
                        std::fs::write(pinned.join(name), text).unwrap();
                    }
                }
                assert_eq!(got, committed(&pinned), "{}: gloci's pinned files", c.name);
            }
        }
    }
}

/// Other readers of the same alignments write the same files.
#[test]
fn cram_and_stdin_match_the_bam() {
    let want = committed(&case_dir("depth-window-100"));
    let (got, _) = run_case(&f("depth.cram"), &["--by".into(), "100".into()]);
    assert_eq!(got, want, "depth.cram");

    let tmp = tempdir();
    gloci()
        .args(["depth", "-", "--by", "100", "--out-prefix"])
        .arg(tmp.path().join("o"))
        .write_stdin(std::fs::read(f("depth.sam")).unwrap())
        .assert()
        .success();
    assert_eq!(written(tmp.path()).0, want, "depth.sam on stdin");
}

#[test]
fn files_have_a_csi_with_the_bed_preset() {
    let tmp = tempdir();
    let p = tmp.path().join("o");
    gloci()
        .args(["depth", &f("depth.bam"), "--by", "100", "--out-prefix"])
        .arg(&p)
        .assert()
        .success();
    assert!(tmp.path().join("o.regions.bed.gz.csi").exists());
    gloci()
        .arg("view")
        .arg(tmp.path().join("o.regions.bed.gz"))
        .args(["--region", "chr1:1-100"])
        .assert()
        .success()
        .stdout("chr1\t0\t100\t0.81\n");
}

#[test]
fn past_reference_end_is_cut_with_one_warning() {
    let tmp = tempdir();
    let o = gloci()
        .args([
            "depth",
            &f("depth-pastend.bam"),
            "--by",
            "100",
            "--out-prefix",
        ])
        .arg(tmp.path().join("o"))
        .assert()
        .success()
        .get_output()
        .clone();
    let err = String::from_utf8(o.stderr).unwrap();
    assert_eq!(err.matches("past the end").count(), 1, "{err}");
    let regions = bgzf_text(&tmp.path().join("o.regions.bed.gz"));
    assert!(regions.contains("chr2\t0\t100\t0.20\n"), "{regions}");
    assert!(regions.contains("chr2\t200\t300\t0.12\n"), "{regions}");
}

#[test]
fn usage_errors() {
    let bam = f("depth.bam");
    let cases: Vec<Vec<&str>> = vec![
        vec![&bam, &bam, "--by", "100", "--out-prefix"],
        vec![&bam, "--by", "0", "--out-prefix"],
        vec![&bam, "--by", "100", "--thresholds", "0", "--out-prefix"],
        vec![&bam, "--thresholds", "1", "--out-prefix"],
        vec![&bam, "--by", "100"],
        vec![&bam, "-o", "x.tsv", "--out-prefix"],
        vec![&bam, "--region", "chr1:1-10", "--out-prefix"],
        vec![&bam, "--min-baseq", "1", "--out-prefix"],
        vec![&bam, "--min-read-length", "1", "--out-prefix"],
        vec![&bam, "--require-flags", "1", "--out-prefix"],
        vec![&bam, "--summary", "--out-prefix"],
        vec![&bam, "--skip-mate-overlap", "--out-prefix"],
        vec![&bam, "--all-positions", "--out-prefix"],
        vec![&bam, "--all-references", "--out-prefix"],
        vec![&bam, "--count-deletions", "--out-prefix"],
        vec![&bam, "--column-names", "--out-prefix"],
        vec![&bam, "--regions-file", "x.bed", "--out-prefix"],
        vec![&bam, "--exclude-flags", "BOGUS", "--out-prefix"],
    ];
    for args in cases {
        let tmp = tempdir();
        let mut cmd = gloci();
        cmd.current_dir(tmp.path()).arg("depth").args(&args);
        // The last argument of most rows is the prefix flag still waiting for its value.
        if args.last() == Some(&"--out-prefix") {
            cmd.arg(tmp.path().join("o"));
        }
        cmd.assert().code(2);
        assert_eq!(
            std::fs::read_dir(tmp.path()).unwrap().count(),
            0,
            "{args:?} wrote files"
        );
    }
}

#[test]
fn bed_region_past_end_publishes_nothing() {
    let tmp = tempdir();
    let p = tmp.path().join("o");
    gloci()
        .args(["depth", &f("depth.bam"), "--by", "100", "--out-prefix"])
        .arg(&p)
        .assert()
        .success();
    let before = snapshot(tmp.path());
    let bad = plain_file(tmp.path(), "bad.bed", "chr1\t480\t501\n");
    gloci()
        .args(["depth", &f("depth.bam"), "--by"])
        .arg(&bad)
        .arg("--out-prefix")
        .arg(&p)
        .assert()
        .code(1)
        .stderr(contains("bad.bed"))
        .stderr(contains("line 1"))
        .stderr(contains("LN 500"));
    let mut after = snapshot(tmp.path());
    // Only bad.bed was added; every earlier file is byte-identical, and no temp file remains.
    assert!(after.remove("bad.bed").is_some());
    assert_eq!(after, before);
}

fn snapshot(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap())
        .map(|e| {
            (
                e.file_name().into_string().unwrap(),
                std::fs::read(e.path()).unwrap(),
            )
        })
        .collect()
}

#[test]
fn unwritable_prefix_fails_cleanly() {
    let tmp = tempdir();
    gloci()
        .args(["depth", &f("depth.bam"), "--by", "100", "--out-prefix"])
        .arg(tmp.path().join("missing-dir/o"))
        .assert()
        .code(1)
        .stderr(contains("o.mosdepth.global.dist.txt"));
    assert_eq!(std::fs::read_dir(tmp.path()).unwrap().count(), 0);
}

#[test]
fn threads_do_not_change_the_files() {
    let run = |t: &str| {
        let tmp = tempdir();
        gloci()
            .args(["depth", &f("mosdepth.bam"), "--by"])
            .arg(f("mosdepth.bed"))
            .args(["--thresholds", "1,2", "-t", t, "--out-prefix"])
            .arg(tmp.path().join("o"))
            .assert()
            .success();
        snapshot(tmp.path())
    };
    let one = run("1");
    assert!(one.contains_key("o.regions.bed.gz.csi"));
    assert_eq!(one, run("4"));
}

#[test]
fn by_windows_huge_header() {
    let tmp = tempdir();
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for i in 0..200_000 {
        sam.push_str(&format!("@SQ\tSN:r{i}\tLN:10\n"));
    }
    sam.push_str("q\t0\tr199999\t1\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n");
    let input = plain_file(tmp.path(), "big.sam", &sam);
    let start = std::time::Instant::now();
    gloci()
        .arg("depth")
        .arg(&input)
        .args(["--by", "5", "--out-prefix"])
        .arg(tmp.path().join("o"))
        .assert()
        .success();
    assert!(start.elapsed().as_secs() < 30, "{:?}", start.elapsed());
    let rows = bgzf_text(&tmp.path().join("o.regions.bed.gz"));
    assert_eq!(rows.lines().count(), 400_000);
}
