//! Packaging invariants for the crates.io release (D79, D82).

use std::path::{Component, Path, PathBuf};
use std::process::Command;

const PUBLISHED: [&str; 6] = [
    "gloci-core",
    "gloci-io",
    "gloci-engine",
    "gloci-vcf",
    "gloci-ops",
    "gloci",
];

fn metadata() -> serde_json::Value {
    let root = xtask::diff::workspace_root();
    let out = Command::new(env!("CARGO"))
        .args(["metadata", "--format-version", "1", "--no-deps"])
        .current_dir(&root)
        .output()
        .unwrap();
    assert!(out.status.success(), "cargo metadata failed");
    serde_json::from_slice(&out.stdout).unwrap()
}

fn package<'a>(meta: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
    meta["packages"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == name)
        .unwrap_or_else(|| panic!("no package {name}"))
}

#[test]
fn published_crates_have_metadata() {
    let meta = metadata();
    for name in PUBLISHED {
        let p = package(&meta, name);
        assert!(
            p["description"].as_str().is_some_and(|d| !d.is_empty()),
            "{name}: description"
        );
        assert_eq!(p["license"], "MIT", "{name}: license");
        assert_eq!(
            p["repository"], "https://github.com/duartemolha/gloci",
            "{name}: repository"
        );
        assert!(p["publish"].is_null(), "{name} must be publishable");
    }
    for name in ["gloci-testkit", "xtask"] {
        let p = package(&meta, name);
        assert_eq!(
            p["publish"],
            serde_json::json!([]),
            "{name} must not publish"
        );
    }
}

#[test]
fn internal_dep_versions_match_workspace() {
    let meta = metadata();
    let version = package(&meta, "gloci")["version"].as_str().unwrap();
    let internal: Vec<&str> = meta["packages"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap())
        .collect();
    for p in meta["packages"].as_array().unwrap() {
        for d in p["dependencies"].as_array().unwrap() {
            let dn = d["name"].as_str().unwrap();
            if !internal.contains(&dn) {
                continue;
            }
            let req = d["req"].as_str().unwrap();
            assert_eq!(
                req,
                format!("^{version}"),
                "{} -> {dn}: path dependency needs version = \"{version}\"",
                p["name"]
            );
        }
    }
}

#[test]
fn licence_file_matches_metadata() {
    let root = xtask::diff::workspace_root();
    let text = std::fs::read_to_string(root.join("LICENSE")).unwrap();
    assert!(text.starts_with("MIT License"));
    assert!(text.contains("Copyright (c) 2026 Duarte Molha"));
    assert!(text.contains("Permission is hereby granted, free of charge"));
    let meta = metadata();
    assert_eq!(package(&meta, "gloci")["license"], "MIT");
}

/// Collapse `.` and `..` lexically, so a path that does not exist can still
/// be compared against a directory.
fn normalise(p: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for c in p.components() {
        match c {
            Component::ParentDir => {
                out.pop();
            }
            Component::CurDir => {}
            c => out.push(c),
        }
    }
    out
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            if p.file_name().is_some_and(|n| n != "target") {
                rust_files(&p, out);
            }
        } else if p.extension().is_some_and(|x| x == "rs") {
            out.push(p);
        }
    }
}

#[test]
fn no_include_escapes_its_crate() {
    let meta = metadata();
    let mut checked = 0;
    for name in PUBLISHED {
        let manifest = PathBuf::from(package(&meta, name)["manifest_path"].as_str().unwrap());
        let crate_dir = manifest.parent().unwrap();
        let mut files = Vec::new();
        rust_files(crate_dir, &mut files);
        for f in files {
            let text = std::fs::read_to_string(&f).unwrap();
            for mac in ["include_str!(", "include_bytes!("] {
                for (at, _) in text.match_indices(mac) {
                    let args = &text[at + mac.len()..];
                    // The first string literal is the path, or the start of a
                    // `concat!` of the path; either way its prefix decides
                    // which directory the include reaches into.
                    let q1 = args.find('"').expect("string literal in include");
                    let rest = &args[q1 + 1..];
                    let lit = &rest[..rest.find('"').unwrap()];
                    let target = normalise(&f.parent().unwrap().join(lit));
                    assert!(
                        target.starts_with(crate_dir),
                        "{}: {mac}{lit}) escapes crate {name}",
                        f.display()
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 0, "found no include macros; the scan is broken");
}

#[test]
fn every_workflow_action_is_pinned() {
    let root = xtask::diff::workspace_root();
    let dir = root.join(".github/workflows");
    let mut seen = 0;
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        if !p.extension().is_some_and(|x| x == "yml" || x == "yaml") {
            continue;
        }
        for (n, line) in std::fs::read_to_string(&p).unwrap().lines().enumerate() {
            let t = line.trim_start().trim_start_matches("- ").trim_start();
            let Some(spec) = t.strip_prefix("uses:") else {
                continue;
            };
            let spec = spec.split('#').next().unwrap().trim();
            if spec.starts_with("./") {
                continue;
            }
            seen += 1;
            let sha = spec.split_once('@').map(|(_, r)| r).unwrap_or("");
            assert!(
                sha.len() == 40 && sha.bytes().all(|b| b.is_ascii_hexdigit()),
                "{}:{}: action `{spec}` is not pinned to a full commit SHA",
                p.display(),
                n + 1
            );
        }
    }
    assert!(seen > 0, "no `uses:` lines found; the scan is broken");
}
