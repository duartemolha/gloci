//! Paper datasets (D95): pins, derived inputs and the generated fetch script.
//! Nothing here touches the network: derived commands run on tiny local files.

use std::collections::HashSet;
use std::fs;

use xtask::bench::{Catalog, derived_command, ensure_derived, record_sha256};
use xtask::diff::workspace_root;
use xtask::paper::fetch::{SCRIPT_PATH, generate_script, load_catalog};

fn read(rel: &str) -> String {
    fs::read_to_string(workspace_root().join(rel)).unwrap()
}

#[test]
fn paper_datasets_have_provider_checksums() {
    let cat = load_catalog().unwrap();
    let paper: Vec<_> = cat.datasets.iter().filter(|d| d.paper).collect();
    for want in [
        "1kg-30x-chr20",
        "gencode-v46-gtf",
        "encode-ccre-v4",
        "encode-ctcf-k562-peaks",
        "encode-ctcf-k562-fc",
        "gnomad-v4.1-exomes-chr22",
        "encode-ctcf-k562-bam",
    ] {
        assert!(paper.iter().any(|d| d.name == want), "missing {want}");
    }
    for d in paper {
        assert!(d.url.starts_with("https://"), "{}: https only", d.name);
        let (alg, hex) = d.checksum().unwrap();
        assert_eq!(alg, "md5", "{}: use the provider's md5", d.name);
        assert!(
            hex.len() == 32 && hex.chars().all(|c| c.is_ascii_hexdigit()),
            "{}: bad md5",
            d.name
        );
    }
    // Derived digests are TO-PIN until recorded; a recorded one is sha256 hex.
    for d in &cat.derived {
        assert!(
            d.sha256.is_empty()
                || (d.sha256.len() == 64 && d.sha256.chars().all(|c| c.is_ascii_hexdigit())),
            "{}: sha256 must be empty (TO-PIN) or 64 hex digits",
            d.name
        );
    }
}

#[test]
fn derived_commands_reference_known_inputs() {
    let cat = load_catalog().unwrap();
    let mut names: HashSet<&str> = cat.datasets.iter().map(|d| d.name.as_str()).collect();
    let mut seen = HashSet::new();
    // A derived input may use only what comes before it.
    for d in &cat.derived {
        assert!(seen.insert(&d.name), "duplicate derived input {}", d.name);
        assert!(!d.from.is_empty(), "{}: no inputs", d.name);
        for f in &d.from {
            match f.strip_prefix("repo:") {
                Some(p) => assert!(
                    workspace_root().join(p).exists(),
                    "{}: repo file {p} is missing",
                    d.name
                ),
                None => assert!(
                    names.contains(f.as_str()),
                    "{}: unknown or later input {f}",
                    d.name
                ),
            }
        }
        assert!(d.command.contains("{out}"), "{}: no {{out}}", d.name);
        assert!(d.command.contains("{in}"), "{}: no {{in}}", d.name);
        let max = d.from.len();
        for i in (max + 1)..=9 {
            assert!(
                !d.command.contains(&format!("{{in{i}}}")),
                "{}: {{in{i}}} has no input",
                d.name
            );
        }
        names.insert(&d.name);
    }
}

#[test]
fn nightly_never_fetches_paper_data() {
    let cat = load_catalog().unwrap();
    let paper: HashSet<&str> = cat
        .datasets
        .iter()
        .filter(|d| d.paper)
        .map(|d| d.name.as_str())
        .collect();
    // 1kg-30x-chr20 was a nightly dataset before the paper used it.
    let cases = read("xtask/bench/cases.toml");
    for c in xtask::bench::parse_cases(&cases).unwrap() {
        if let Some(ds) = &c.dataset {
            assert!(
                !paper.contains(ds.as_str()) || ds == "1kg-30x-chr20",
                "nightly case {} uses paper dataset {ds}",
                c.name
            );
        }
        for d in &cat.derived {
            assert!(
                c.dataset.as_deref() != Some(d.name.as_str()),
                "nightly case {} uses derived input {}",
                c.name,
                d.name
            );
        }
    }
    let nightly = read(".github/workflows/nightly.yml");
    for d in &paper {
        if *d != "1kg-30x-chr20" {
            assert!(!nightly.contains(d), "nightly.yml mentions {d}");
        }
    }
}

#[test]
fn no_workflow_runs_paper_downloads() {
    let dir = workspace_root().join(".github/workflows");
    for e in fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        let text = fs::read_to_string(&p).unwrap();
        for bad in ["xtask paper", "fetch-datasets.sh"] {
            assert!(!text.contains(bad), "{} contains `{bad}`", p.display());
        }
    }
}

fn conda_pins(dockerfile: &str) -> Vec<(String, String)> {
    let want = [
        "bedtools", "bedops", "vcftools", "htslib", "bcftools", "samtools",
    ];
    dockerfile
        .split_whitespace()
        .filter_map(|w| w.split_once('='))
        .filter(|(n, _)| want.contains(n))
        .map(|(n, v)| (n.to_string(), v.to_string()))
        .collect()
}

#[test]
fn bench_and_diff_images_pin_same_reference_tools() {
    let bench = conda_pins(&read("xtask/docker/Dockerfile.bench"));
    let diff = conda_pins(&read("xtask/docker/Dockerfile"));
    assert_eq!(bench.len(), 6, "Dockerfile.bench pins: {bench:?}");
    assert_eq!(bench, diff, "the two images must pin the same versions");
    let b = read("xtask/docker/Dockerfile.bench");
    for tool in ["gawk=", "time=", "ucsc-bigwigtobedgraph="] {
        assert!(b.contains(tool), "Dockerfile.bench must pin {tool}");
    }
}

#[test]
fn fetch_script_is_current() {
    let script = generate_script(
        &load_catalog().unwrap(),
        &read("xtask/docker/Dockerfile.bench"),
    )
    .unwrap();
    assert_eq!(
        read(SCRIPT_PATH),
        script,
        "paper/fetch-datasets.sh is stale: run `cargo xtask paper fetch-script`"
    );
    assert!(script.starts_with("#!/bin/sh\n"));
}

fn tiny_catalog(sha: &str) -> Catalog {
    Catalog::parse(&format!(
        r#"
[[dataset]]
name = "x"
url = "https://example.invalid/x.txt"
md5 = "00000000000000000000000000000000"

[[derived]]
name = "up.txt"
from = ["repo:xtask/bench/datasets.toml"]
command = "head -n 1 {{in}} | tr a-z A-Z > {{out}} && echo idx > {{out}}.idx"
sha256 = "{sha}"
"#
    ))
    .unwrap()
}

#[test]
fn derived_input_is_built_checked_and_recorded() {
    let tmp = tempfile::tempdir().unwrap();
    let cache = tmp.path();

    // No digest yet: an error without --record that says how to fix it.
    let cat = tiny_catalog("");
    let err = ensure_derived(&cat, &cat.derived[0], cache, None, false).unwrap_err();
    assert!(err.to_string().contains("paper fetch --record"), "{err}");
    assert!(
        cache.join("up.txt.idx").exists(),
        "sidecar follows the output"
    );

    // A wrong digest is a mismatch.
    let cat = tiny_catalog(&"0".repeat(64));
    let err = ensure_derived(&cat, &cat.derived[0], cache, None, false).unwrap_err();
    assert!(err.to_string().contains("sha256 mismatch"), "{err}");

    // The right digest passes.
    let out = cache.join("up.txt");
    let hex = {
        use sha2::Digest;
        sha2::Sha256::digest(fs::read(&out).unwrap())
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>()
    };
    let cat = tiny_catalog(&hex);
    assert_eq!(
        ensure_derived(&cat, &cat.derived[0], cache, None, false).unwrap(),
        out
    );
}

#[test]
fn derived_command_and_record_helpers() {
    let c = derived_command(
        "cat {in} {in2} > {out}",
        &["a".to_string(), "b".to_string()],
        "o",
    );
    assert_eq!(c, "cat a b > o");
    let text = "[[derived]]\nname = \"a\"\nsha256 = \"\" # TO-PIN\n\n[[derived]]\nname = \"b\"\nsha256 = \"\"\n";
    let out = record_sha256(text, "b", "ff").unwrap();
    assert!(out.contains("name = \"b\"\nsha256 = \"ff\""));
    assert!(out.contains("sha256 = \"\" # TO-PIN"), "a is untouched");
    assert!(record_sha256(text, "zzz", "ff").is_err());
}

/// Run a generated script on a tiny local catalog: the dataset is pre-seeded
/// in the cache (no network), the commands are `tr` and friends.
fn run_script(
    only: Option<&str>,
    allow: bool,
    break_cmd: bool,
    tmp: &std::path::Path,
) -> std::process::Output {
    use md5::Digest;
    let md5 = md5::Md5::digest(b"hello\n")
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<String>();
    let second = if break_cmd {
        "zcat /nonexistent.gz | cat > {out}"
    } else {
        "tr a-z A-Z < {in} > {out}"
    };
    let cat = Catalog::parse(&format!(
        r#"
[[dataset]]
name = "x"
url = "https://example.invalid/x.txt"
md5 = "{md5}"
paper = true

[[derived]]
name = "up.txt"
from = ["x"]
command = "{second}"
sha256 = ""

[[derived]]
name = "twice.txt"
from = ["up.txt"]
command = "cat {{in}} {{in}} > {{out}}"
sha256 = ""
"#
    ))
    .unwrap();
    let script = generate_script(&cat, "").unwrap();
    fs::create_dir_all(tmp.join("paper")).unwrap();
    let path = tmp.join("paper/s.sh");
    fs::write(&path, script).unwrap();
    let cache = tmp.join("cache");
    fs::create_dir_all(&cache).unwrap();
    fs::write(cache.join("x.txt"), "hello\n").unwrap();
    let mut c = std::process::Command::new("sh");
    c.arg(&path).env("GLOCI_BENCH_CACHE", &cache);
    if allow {
        c.env("GLOCI_ALLOW_UNPINNED", "1");
    }
    if let Some(o) = only {
        c.arg(o);
    }
    c.output().unwrap()
}

#[test]
fn script_resolves_dependencies_and_rejects_unknown_names() {
    let tmp = tempfile::tempdir().unwrap();
    let o = run_script(Some("nope"), true, false, tmp.path());
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("unknown dataset"));

    let o = run_script(Some("twice.txt"), true, false, tmp.path());
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let cache = tmp.path().join("cache");
    assert_eq!(fs::read_to_string(cache.join("up.txt")).unwrap(), "HELLO\n");
    assert_eq!(
        fs::read_to_string(cache.join("twice.txt")).unwrap(),
        "HELLO\nHELLO\n"
    );
}

#[test]
fn script_fails_on_unpinned_digests_unless_allowed() {
    let tmp = tempfile::tempdir().unwrap();
    let o = run_script(None, false, false, tmp.path());
    assert_eq!(o.status.code(), Some(1));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("2 derived inputs UNVERIFIED"), "{err}");
    let o = run_script(None, true, false, tmp.path());
    assert!(o.status.success());
}

#[test]
fn script_runs_derive_commands_under_pipefail() {
    let tmp = tempfile::tempdir().unwrap();
    let o = run_script(None, true, true, tmp.path());
    assert!(
        !o.status.success(),
        "a failing pipe stage must fail the run"
    );
    assert!(
        !tmp.path().join("cache/up.txt").exists(),
        "no result is left behind"
    );
}

#[test]
fn derived_input_without_its_index_is_an_error() {
    let tmp = tempfile::tempdir().unwrap();
    let cat = Catalog::parse(
        r#"
[[derived]]
name = "a.bam"
from = ["repo:xtask/bench/datasets.toml"]
command = "cp {in} {out}"
sha256 = ""
"#,
    )
    .unwrap();
    let err = ensure_derived(&cat, &cat.derived[0], tmp.path(), None, true).unwrap_err();
    assert!(err.to_string().contains(".bai"), "{err}");
}
