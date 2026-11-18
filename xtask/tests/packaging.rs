//! Packaging invariants for the crates.io release (D79, D82).

use std::path::{Component, Path, PathBuf};
use std::process::Command;

const PUBLISHED: [&str; 8] = [
    "gloci-core",
    "gloci-chains",
    "gloci-chains-t2t",
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
        if !PUBLISHED.contains(&p["name"].as_str().unwrap()) {
            continue;
        }
        for d in p["dependencies"].as_array().unwrap() {
            let dn = d["name"].as_str().unwrap();
            if !internal.contains(&dn) {
                continue;
            }
            if !PUBLISHED.contains(&dn) {
                continue; // unpublished crates are path-only (next test)
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
fn published_crates_do_not_version_unpublished_crates() {
    let meta = metadata();
    for name in PUBLISHED {
        for d in package(&meta, name)["dependencies"].as_array().unwrap() {
            let dn = d["name"].as_str().unwrap();
            if matches!(dn, "gloci-testkit" | "xtask") {
                assert_eq!(
                    d["req"], "*",
                    "{name} -> {dn} ({}): an unpublished crate must be a path-only dependency",
                    d["kind"]
                );
            }
        }
    }
}

#[test]
fn every_published_crate_ships_the_licence() {
    let root = xtask::diff::workspace_root();
    let want = std::fs::read(root.join("LICENSE")).unwrap();
    for name in PUBLISHED {
        let got = std::fs::read(root.join("crates").join(name).join("LICENSE"))
            .unwrap_or_else(|_| panic!("{name}: missing LICENSE"));
        assert!(got == want, "{name}: LICENSE differs from the root copy");
    }
}

/// Every release archive and the Python wheel carry the notices of the crates compiled in (D80).
#[test]
fn every_archive_ships_the_third_party_licences() {
    let root = xtask::diff::workspace_root();
    let text = std::fs::read_to_string(root.join("THIRD-PARTY-LICENSES")).unwrap();
    assert!(text.starts_with("THIRD-PARTY-LICENSES: "));
    assert!(text.contains("cargo xtask licenses"));
    let release = std::fs::read_to_string(root.join(".github/workflows/release.yml")).unwrap();
    assert!(
        release.contains("cp THIRD-PARTY-LICENSES \"dist/${STAGE}/\""),
        "the archive step must copy THIRD-PARTY-LICENSES next to LICENSE"
    );
    let ci = std::fs::read_to_string(root.join(".github/workflows/ci.yml")).unwrap();
    assert!(
        ci.contains("cargo xtask licenses --check"),
        "CI must check THIRD-PARTY-LICENSES"
    );
    let py = root.join("crates/gloci-py");
    let pyproject = std::fs::read_to_string(py.join("pyproject.toml")).unwrap();
    assert!(pyproject.contains("license-files = [\"LICENSE\", \"THIRD-PARTY-LICENSES\"]"));
    assert!(
        std::fs::read(py.join("THIRD-PARTY-LICENSES")).unwrap() == text.as_bytes(),
        "crates/gloci-py/THIRD-PARTY-LICENSES must be the root file (a symlink)"
    );
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

/// A lexical guard: it reads the first path literal of every include macro
/// and checks it stays inside the crate. It is backed by the real check, the
/// verified `cargo publish --dry-run` in CI, which builds each crate from its
/// tarball. Paths it cannot read as literals fail with the file and line.
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
            for (i, line) in text.lines().enumerate() {
                if line.trim_start().starts_with("//") {
                    continue;
                }
                for mac in ["include_str!(", "include_bytes!("] {
                    for (at, _) in line.match_indices(mac) {
                        let mut args = line[at + mac.len()..].trim_start();
                        if let Some(r) = args.strip_prefix("concat!(") {
                            args = r.trim_start();
                        }
                        let mut base = f.parent().unwrap().to_path_buf();
                        if let Some(r) = args.strip_prefix("env!(\"CARGO_MANIFEST_DIR\")") {
                            base = crate_dir.to_path_buf();
                            args = r.trim_start().trim_start_matches(',').trim_start();
                        }
                        let lit = args.strip_prefix('"').and_then(|r| r.split_once('"'));
                        let Some((lit, _)) = lit else {
                            panic!(
                                "{}:{}: include path is not a string literal; check by hand",
                                f.display(),
                                i + 1
                            );
                        };
                        let target = normalise(&base.join(lit.trim_start_matches('/')));
                        assert!(
                            target.starts_with(crate_dir),
                            "{}:{}: {mac}{lit}) escapes crate {name}",
                            f.display(),
                            i + 1
                        );
                        checked += 1;
                    }
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

/// The CI publish dry run and the release loop name every published crate, in dependency order,
/// so a crate added to `PUBLISHED` cannot be left out of either (D82).
#[test]
fn workflows_publish_every_published_crate_in_order() {
    let root = xtask::diff::workspace_root();
    let read = |f: &str| std::fs::read_to_string(root.join(".github/workflows").join(f)).unwrap();
    let ci = read("ci.yml");
    let dry_run = ci
        .lines()
        .find(|l| l.contains("cargo publish --dry-run"))
        .expect("ci.yml has a publish dry run");
    let flags: Vec<&str> = dry_run
        .split_whitespace()
        .skip_while(|w| *w != "-p")
        .collect::<Vec<_>>()
        .chunks(2)
        .map(|c| {
            assert_eq!(c[0], "-p", "ci.yml: {dry_run}");
            c[1]
        })
        .collect();
    assert_eq!(flags, PUBLISHED, "ci.yml's publish dry run");
    let release = read("release.yml");
    let looped: Vec<&str> = release
        .lines()
        .find_map(|l| l.trim().strip_prefix("for crate in "))
        .and_then(|l| l.strip_suffix("; do"))
        .expect("release.yml has a publish loop")
        .split_whitespace()
        .collect();
    assert_eq!(looped, PUBLISHED, "release.yml's publish loop");
}
