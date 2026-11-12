//! `cargo xtask paper`: tables, numbers and the Zenodo bundle (D98), checked
//! on a small fixture of results (three cases, five concordance rows, one of
//! them a failed run).
use std::collections::BTreeSet;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use xtask::diff::workspace_root;
use xtask::paper::bundle::bundle;
use xtask::paper::numbers::numbers;
use xtask::paper::tables::tables;

fn fixture() -> PathBuf {
    workspace_root().join("xtask/tests/fixtures/paper-results/0.0.0-fixture")
}

fn expected() -> PathBuf {
    workspace_root().join("xtask/tests/fixtures/paper-expected")
}

/// Every file under `dir`, as `/`-separated relative path.
fn files(dir: &Path) -> BTreeSet<String> {
    fn walk(base: &Path, d: &Path, out: &mut BTreeSet<String>) {
        for e in fs::read_dir(d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                walk(base, &p, out);
            } else {
                out.insert(
                    p.strip_prefix(base)
                        .unwrap()
                        .to_string_lossy()
                        .replace('\\', "/"),
                );
            }
        }
    }
    let mut out = BTreeSet::new();
    walk(dir, dir, &mut out);
    out
}

/// Tables go to `<tmp>/tables`, numbers beside them in `<tmp>`.
fn generate(tmp: &Path) {
    tables(&fixture(), &tmp.join("tables")).unwrap();
}

#[test]
fn tables_from_fixture_match_expected() {
    let tmp = tempfile::tempdir().unwrap();
    generate(tmp.path());
    if std::env::var_os("UPDATE_PAPER_EXPECTED").is_some() {
        let _ = fs::remove_dir_all(expected());
        fs::create_dir_all(expected().join("tables")).unwrap();
        for f in files(tmp.path()) {
            fs::copy(tmp.path().join(&f), expected().join(&f)).unwrap();
        }
    }
    assert_eq!(
        files(tmp.path()),
        files(&expected()),
        "set of generated files"
    );
    for f in files(&expected()) {
        let want = fs::read_to_string(expected().join(&f)).unwrap();
        let got = fs::read_to_string(tmp.path().join(&f)).unwrap();
        assert_eq!(got, want, "{f} differs (UPDATE_PAPER_EXPECTED=1 rewrites)");
    }
    let names: Vec<String> = files(&expected().join("tables")).into_iter().collect();
    assert_eq!(
        names,
        [
            "T1-datasets.csv",
            "T2-benchmarks.csv",
            "T3-speedups.csv",
            "T4-concordance.csv",
            "T5-differences.csv",
            "T6-command-coverage.csv",
            "T7-tests.csv",
            "T8-software.csv"
        ]
    );
}

#[test]
fn tables_are_deterministic() {
    let a = tempfile::tempdir().unwrap();
    let b = tempfile::tempdir().unwrap();
    generate(a.path());
    generate(b.path());
    assert_eq!(files(a.path()), files(b.path()));
    for f in files(a.path()) {
        assert_eq!(
            fs::read(a.path().join(&f)).unwrap(),
            fs::read(b.path().join(&f)).unwrap(),
            "{f}"
        );
    }
}

#[test]
fn numbers_tex_defines_every_key() {
    let tmp = tempfile::tempdir().unwrap();
    generate(tmp.path());
    let tex = fs::read_to_string(tmp.path().join("numbers.tex")).unwrap();
    let mut lines = tex.lines();
    assert_eq!(
        lines.next().unwrap(),
        r"\providecommand{\glnum}[1]{\csname glnum@#1\endcsname}"
    );
    let body: Vec<&str> = lines.collect();
    let map = numbers(&fixture()).unwrap();
    assert_eq!(body.len(), map.len());
    let mut sorted = body.clone();
    sorted.sort();
    assert_eq!(body, sorted, "one sorted line per key");
    for key in map.keys() {
        let prefix = format!(r"\expandafter\def\csname glnum@{key}\endcsname{{");
        assert!(
            body.iter()
                .any(|l| l.starts_with(&prefix) && l.ends_with('}')),
            "no definition for {key}"
        );
    }
    for key in [
        "bench.merge-exons.gloci.median_s",
        "bench.merge-exons.bedtools.max_rss_mb",
        "bench.merge-exons.speedup.bedtools",
        "concordance.cases",
        "concordance.compat.at_100",
        "concordance.default.mean_match_rate",
        "concordance.sort-exons.default.match_rate",
        "tests.golden_cases",
        "tests.golden.bedtools",
        "differences.anchors",
        "datasets.1kg-30x-chr20.records",
        "migration.bedtools.available",
        "migration.bedtools.total",
    ] {
        assert!(map.contains_key(key), "missing key {key}");
    }
    assert_eq!(map["migration.bedtools.available"], 1);
    assert_eq!(map["migration.bedtools.total"], 2);
    // The failed run is a case that ran: it counts, without a rate.
    assert_eq!(map["concordance.cases"], 3);
    assert_eq!(map["concordance.compat.at_100"], 2);
}

/// (path, mtime, uid, gid, mode, is_dir) of every entry, in archive order.
fn entries(bytes: &[u8]) -> Vec<(String, u64, u64, u64, u32, bool)> {
    let gz = flate2::read::GzDecoder::new(bytes);
    let mut ar = tar::Archive::new(gz);
    ar.entries()
        .unwrap()
        .map(|e| {
            let e = e.unwrap();
            let h = e.header();
            (
                e.path().unwrap().to_string_lossy().into_owned(),
                h.mtime().unwrap(),
                h.uid().unwrap(),
                h.gid().unwrap(),
                h.mode().unwrap(),
                h.entry_type().is_dir(),
            )
        })
        .collect()
}

#[test]
fn bundle_is_reproducible() {
    let tmp = tempfile::tempdir().unwrap();
    let (a, b) = (tmp.path().join("a.tar.gz"), tmp.path().join("b.tar.gz"));
    bundle(&fixture(), &a).unwrap();
    bundle(&fixture(), &b).unwrap();
    let bytes = fs::read(&a).unwrap();
    assert_eq!(bytes, fs::read(&b).unwrap(), "two runs are byte-identical");
    assert_eq!(&bytes[..2], &[0x1f, 0x8b]);
    assert_eq!(&bytes[4..8], &[0, 0, 0, 0], "gzip header mtime is 0");

    let es = entries(&bytes);
    let paths: Vec<&str> = es.iter().map(|e| e.0.as_str()).collect();
    let mut sorted = paths.clone();
    sorted.sort();
    assert_eq!(paths, sorted, "entries are sorted");
    // 2026-11-12T10:30:15Z
    let finished = 1_794_479_415;
    for (p, mtime, uid, gid, mode, dir) in &es {
        assert_eq!(*mtime, finished, "{p}: mtime is finished_at");
        assert_eq!((*uid, *gid), (0, 0), "{p}");
        assert_eq!(*mode, if *dir { 0o755 } else { 0o644 }, "{p}");
    }
    let top = paths[0].trim_end_matches('/').to_string();
    assert!(
        top.starts_with("gloci-paper-evidence-0.0.0-fixture"),
        "{top}"
    );
    let got_files: BTreeSet<String> = es
        .iter()
        .filter(|e| !e.5)
        .map(|e| e.0.strip_prefix(&format!("{top}/")).unwrap().to_string())
        .collect();

    let root = workspace_root();
    let mut want: BTreeSet<String> = files(&fixture())
        .into_iter()
        .map(|f| format!("results/0.0.0-fixture/{f}"))
        .collect();
    let tmp2 = tempfile::tempdir().unwrap();
    generate(tmp2.path());
    want.extend(files(tmp2.path()));
    for f in [
        "xtask/bench/datasets.toml",
        "paper/fetch-datasets.sh",
        "xtask/docker/Dockerfile",
        "xtask/docker/Dockerfile.bench",
        "CITATION.cff",
        "README.txt",
    ] {
        want.insert(f.to_string());
    }
    for e in fs::read_dir(root.join("xtask/paper")).unwrap() {
        let p = e.unwrap().path();
        if p.extension().is_some_and(|x| x == "toml") {
            want.insert(format!(
                "xtask/paper/{}",
                p.file_name().unwrap().to_string_lossy()
            ));
        }
    }
    assert_eq!(got_files, want);

    // Content of a repo file is the file itself; the tables are the ones
    // `tables` writes.
    let gz = flate2::read::GzDecoder::new(&bytes[..]);
    let mut ar = tar::Archive::new(gz);
    for e in ar.entries().unwrap() {
        let mut e = e.unwrap();
        let p = e.path().unwrap().to_string_lossy().into_owned();
        let rel = p.strip_prefix(&format!("{top}/")).unwrap_or("").to_string();
        if rel == "CITATION.cff" || rel == "tables/T4-concordance.csv" {
            let mut got = Vec::new();
            e.read_to_end(&mut got).unwrap();
            let src = if rel == "CITATION.cff" {
                root.join(&rel)
            } else {
                tmp2.path().join(&rel)
            };
            assert_eq!(got, fs::read(src).unwrap(), "{rel}");
        }
    }
}

/// With results committed under `paper/results/`, `paper/tables/` and
/// `paper/numbers.*` must be what `cargo xtask paper` writes from them.
#[test]
fn committed_tables_are_current() {
    let results = workspace_root().join("paper/results");
    let mut versions: Vec<PathBuf> = match fs::read_dir(&results) {
        Ok(rd) => rd
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.is_dir())
            .collect(),
        Err(_) => Vec::new(),
    };
    versions.sort();
    let Some(latest) = versions.pop() else {
        eprintln!(
            "skipped: paper/results is empty (run `cargo xtask paper run` on the reference machine)"
        );
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    tables(&latest, &tmp.path().join("tables")).unwrap();
    for f in files(tmp.path()) {
        assert_eq!(
            fs::read_to_string(tmp.path().join(&f)).unwrap(),
            fs::read_to_string(workspace_root().join("paper").join(&f)).unwrap_or_default(),
            "paper/{f} is stale: run `cargo xtask paper`"
        );
    }
}

/// `paper run` writes both files; here on this repository with no datasets
/// present (sizes are null) and no container (tool versions are read on the
/// host, or null).
#[test]
fn run_records_describe_this_repo_and_host() {
    use xtask::paper::{environment, snapshot};
    let root = workspace_root();
    let cat = xtask::paper::fetch::load_catalog().unwrap();
    let empty = tempfile::tempdir().unwrap();
    let snap = snapshot::snapshot(&root, &cat, empty.path(), "9.9.9", None).unwrap();
    let golden_dirs = fs::read_dir(root.join("crates/gloci/tests/golden"))
        .unwrap()
        .filter(|e| e.as_ref().unwrap().path().is_dir())
        .count() as u64;
    assert_eq!(snap["golden"]["total"], golden_dirs);
    let by_tool: u64 = snap["golden"]["by_tool"]
        .as_object()
        .unwrap()
        .values()
        .map(|v| v.as_u64().unwrap())
        .sum();
    assert_eq!(by_tool, golden_dirs);
    assert!(snap["golden"]["by_tool"]["bedtools"].as_u64().unwrap() > 0);
    assert!(snap["anchors"].as_array().unwrap().len() > 100);
    assert!(snap["migration"].as_array().unwrap().len() > 80);
    assert!(snap["tests"]["rust-crates"].as_u64().unwrap() > 100);
    assert!(snap["datasets"]["exons.bed"]["bytes"].is_null());
    assert!(snap["commit"].is_null());

    let cases = xtask::paper::bench::load_paper_cases().unwrap();
    let wheel = empty
        .path()
        .join("gloci-0.1.0-cp39-abi3-manylinux_2_28_x86_64.whl");
    fs::write(&wheel, b"not really a wheel").unwrap();
    let env = environment::environment(&environment::EnvArgs {
        version: "9.9.9",
        commit: environment::git_commit(&root),
        started_at: "2026-01-01T00:00:00Z",
        finished_at: "2026-01-01T00:00:01Z",
        image: None,
        gloci: None,
        wheel: Some(&wheel),
        cases: &cases,
        extra_tools: &["vcftools".to_string()],
    })
    .unwrap();
    if std::env::var_os("SHOW_PAPER_ENV").is_some() {
        eprintln!("{}", serde_json::to_string_pretty(&env).unwrap());
    }
    assert!(env["host"]["cores"].as_u64().is_some());
    assert!(env["docker"].is_null());
    assert!(env["tools"].as_object().unwrap().contains_key("bedtools"));
    assert!(env["tools"].as_object().unwrap().contains_key("gnu-sort"));
    assert_eq!(
        env["wheel"]["file_name"],
        "gloci-0.1.0-cp39-abi3-manylinux_2_28_x86_64.whl"
    );
    assert_eq!(env["wheel"]["sha256"].as_str().unwrap().len(), 64);
    // A missing wheel is refused.
    assert!(environment::wheel_info(&empty.path().join("nope.whl")).is_err());
}

/// A failed run's row (`,,,,,failed,`) has no counts: it loads, T4 shows it
/// with empty cells, and `numbers` has no rate for it.
#[test]
fn a_failed_concordance_row_has_no_counts() {
    let r = xtask::paper::data::load(&fixture()).unwrap();
    let f = r
        .concordance
        .iter()
        .find(|c| c.case == "merge-exons-failed")
        .expect("the fixture's failed row");
    assert_eq!(
        (f.n_gloci, f.n_reference, f.n_matching, f.match_rate),
        (None, None, None, None)
    );
    let tmp = tempfile::tempdir().unwrap();
    generate(tmp.path());
    let t4 = fs::read_to_string(tmp.path().join("tables/T4-concordance.csv")).unwrap();
    assert!(
        t4.lines()
            .any(|l| l == "merge-exons-failed,merge,bedtools,exons.bed,default,,,,,"),
        "{t4}"
    );
    let map = numbers(&fixture()).unwrap();
    assert!(!map.contains_key("concordance.merge-exons-failed.default.match_rate"));
    assert!(map.contains_key("concordance.sort-exons.default.match_rate"));
}
