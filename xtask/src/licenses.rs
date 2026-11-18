//! THIRD-PARTY-LICENSES (D80): the notices of the crates compiled into the `gloci` binary and the
//! `gloci-py` extension, generated from `cargo metadata` and shipped in every release archive.
use anyhow::{Context, Result, bail};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Where the file is committed, relative to the repository.
pub const FILE: &str = "THIRD-PARTY-LICENSES";

/// The workspace packages whose normal dependencies ship: the binary and the Python extension.
const ROOTS: [&str; 2] = ["gloci", "gloci-py"];

/// File-name prefixes (compared in upper case) of the licence files of a package's source.
const LICENCE_PREFIXES: [&str; 5] = ["LICENSE", "LICENCE", "COPYING", "NOTICE", "UNLICENSE"];

/// Standard texts, for a package whose source has no licence file.
const STANDARD_TEXTS: [(&str, &str); 3] = [
    ("MIT", include_str!("../licenses/MIT.txt")),
    ("Apache-2.0", include_str!("../licenses/Apache-2.0.txt")),
    ("Zlib", include_str!("../licenses/Zlib.txt")),
];

/// A package compiled into a shipped artefact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Package {
    pub name: String,
    pub version: String,
    /// The SPDX expression of its manifest.
    pub license: String,
    pub repository: Option<String>,
    /// The directory holding its source.
    pub dir: PathBuf,
}

fn normal_dep(d: &Value) -> bool {
    d["dep_kinds"]
        .as_array()
        .is_some_and(|ks| ks.iter().any(|k| k["kind"].is_null()))
}

/// The external packages in the normal (not dev, not build) dependency closure of the
/// [`ROOTS`], from `cargo metadata --format-version 1` output without a platform filter, sorted
/// by name and version. Workspace members (the packages without a `source`) are left out: gloci's
/// own LICENSE covers them.
pub fn closure(meta: &Value) -> Result<Vec<Package>> {
    let packages: HashMap<&str, &Value> = meta["packages"]
        .as_array()
        .context("metadata has no packages")?
        .iter()
        .filter_map(|p| Some((p["id"].as_str()?, p)))
        .collect();
    let nodes: HashMap<&str, &Value> = meta["resolve"]["nodes"]
        .as_array()
        .context("metadata has no resolve graph (run it without --no-deps)")?
        .iter()
        .filter_map(|n| Some((n["id"].as_str()?, n)))
        .collect();
    let mut todo: Vec<&str> = meta["workspace_members"]
        .as_array()
        .context("metadata has no workspace_members")?
        .iter()
        .filter_map(Value::as_str)
        .filter(|id| {
            packages
                .get(id)
                .is_some_and(|p| ROOTS.contains(&p["name"].as_str().unwrap_or("")))
        })
        .collect();
    if todo.len() != ROOTS.len() {
        bail!("the workspace lacks one of {ROOTS:?}");
    }
    let mut seen: BTreeSet<&str> = todo.iter().copied().collect();
    while let Some(id) = todo.pop() {
        let node = nodes.get(id).with_context(|| format!("no node for {id}"))?;
        for d in node["deps"].as_array().into_iter().flatten() {
            let dep = d["pkg"].as_str().context("dependency without pkg")?;
            if normal_dep(d) && seen.insert(dep) {
                todo.push(dep);
            }
        }
    }
    let mut out = Vec::new();
    for id in seen {
        let p = packages
            .get(id)
            .with_context(|| format!("no package for {id}"))?;
        if p["source"].is_null() {
            continue;
        }
        let name = p["name"].as_str().unwrap_or_default().to_string();
        let manifest = p["manifest_path"].as_str().unwrap_or_default();
        out.push(Package {
            version: p["version"].as_str().unwrap_or_default().to_string(),
            license: p["license"]
                .as_str()
                .with_context(|| format!("{name} declares no licence expression"))?
                .to_string(),
            repository: p["repository"].as_str().map(str::to_string),
            dir: Path::new(manifest)
                .parent()
                .map(Path::to_path_buf)
                .unwrap_or_default(),
            name,
        });
    }
    out.sort_by(|a, b| (&a.name, &a.version).cmp(&(&b.name, &b.version)));
    Ok(out)
}

/// The texts of the licence files in `dir` (names starting with LICENSE, LICENCE, COPYING, NOTICE
/// or UNLICENSE, any case), in file-name order, joined by a blank line. `None` without any.
pub fn licence_files_text(dir: &Path) -> Result<Option<String>> {
    let mut names = Vec::new();
    for e in std::fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let e = e?;
        let name = e.file_name().to_string_lossy().into_owned();
        let upper = name.to_ascii_uppercase();
        if e.path().is_file() && LICENCE_PREFIXES.iter().any(|p| upper.starts_with(p)) {
            names.push(name);
        }
    }
    names.sort();
    let texts: Vec<String> = names
        .iter()
        .map(|n| {
            let bytes = std::fs::read(dir.join(n))?;
            Ok(clean(&String::from_utf8_lossy(&bytes)))
        })
        .collect::<Result<_>>()?;
    Ok((!texts.is_empty()).then(|| texts.join("\n\n")))
}

/// `text` with LF line endings, no trailing blanks on a line and no blank lines around it.
fn clean(text: &str) -> String {
    let lines: Vec<&str> = text.lines().map(str::trim_end).collect();
    lines.join("\n").trim_matches('\n').to_string()
}

/// The licence identifiers of an SPDX expression, in order, without duplicates.
pub fn spdx_ids(expr: &str) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for tok in expr
        .replace(['(', ')', '/'], " ")
        .split_whitespace()
        .filter(|t| !matches!(*t, "AND" | "OR" | "WITH"))
    {
        if !out.iter().any(|o| o == tok) {
            out.push(tok.to_string());
        }
    }
    out
}

/// The text of a package without a licence file: the standard text of each licence of its
/// expression that gloci carries (a copyleft alternative of an `OR` is not carried), then the
/// package's repository. An error when none is carried.
fn fallback_text(p: &Package) -> Result<String> {
    let mut parts = vec![format!(
        "{} {} has no licence file in its source. The standard text of its licence follows, \
         without the copyright line of its authors; its repository is {}.",
        p.name,
        p.version,
        p.repository.as_deref().unwrap_or("not declared")
    )];
    let mut found = false;
    for id in spdx_ids(&p.license) {
        if let Some((_, text)) = STANDARD_TEXTS.iter().find(|(n, _)| *n == id) {
            parts.push(clean(text));
            found = true;
        }
    }
    if !found {
        bail!(
            "{} {} ({}) has no licence file and xtask/licenses/ has no text for it",
            p.name,
            p.version,
            p.license
        );
    }
    Ok(parts.join("\n\n"))
}

/// The line MPL-2.0 requires of a distributor: where the source of the package is.
fn source_line(p: &Package) -> Option<String> {
    spdx_ids(&p.license).iter().any(|i| i == "MPL-2.0").then(|| {
        format!(
            "Source of {} {} (MPL-2.0, used unmodified): https://crates.io/crates/{}/{}, repository {}",
            p.name,
            p.version,
            p.name,
            p.version,
            p.repository.as_deref().unwrap_or("not declared")
        )
    })
}

const RULE: &str =
    "================================================================================";

/// The file's text. `texts` holds the licence-file text of each package by (name, version), none
/// for a package without one. Packages with an identical text share one entry, headed by a
/// `name version (licence)` line each.
pub fn render(
    packages: &[Package],
    texts: &BTreeMap<(String, String), Option<String>>,
) -> Result<String> {
    // text -> (packages, source lines)
    let mut groups: Vec<(String, Vec<&Package>)> = Vec::new();
    for p in packages {
        let text = match texts.get(&(p.name.clone(), p.version.clone())) {
            Some(Some(t)) => t.clone(),
            _ => fallback_text(p)?,
        };
        match groups.iter_mut().find(|(t, _)| *t == text) {
            Some((_, ps)) => ps.push(p),
            None => groups.push((text, vec![p])),
        }
    }
    let mut out = format!(
        "THIRD-PARTY-LICENSES: the licences of the {} crates compiled into gloci's binary and \
         Python extension.\n\
         Generated by `cargo xtask licenses` from Cargo.lock (CI runs `cargo xtask licenses --check`); \
         do not edit.\n\
         gloci itself is under the licence in LICENSE.\n",
        packages.len()
    );
    for (text, ps) in &groups {
        out.push_str(&format!("\n{RULE}\n"));
        for p in ps {
            out.push_str(&format!("{} {} ({})\n", p.name, p.version, p.license));
        }
        for line in ps.iter().filter_map(|p| source_line(p)) {
            out.push_str(&line);
            out.push('\n');
        }
        out.push_str(&format!("{RULE}\n\n{text}\n"));
    }
    Ok(out)
}

/// Generate the file's text for the workspace at `root`.
pub fn generate(root: &Path) -> Result<String> {
    let run = |args: &[&str]| -> Result<Vec<u8>> {
        let o = Command::new("cargo")
            .args(args)
            .current_dir(root)
            .output()
            .context("running cargo")?;
        if !o.status.success() {
            bail!(
                "cargo {} failed: {}",
                args.join(" "),
                String::from_utf8_lossy(&o.stderr)
            );
        }
        Ok(o.stdout)
    };
    run(&["fetch", "--locked"])?;
    let meta: Value =
        serde_json::from_slice(&run(&["metadata", "--format-version", "1", "--locked"])?)?;
    let packages = closure(&meta)?;
    let mut texts = BTreeMap::new();
    for p in &packages {
        texts.insert(
            (p.name.clone(), p.version.clone()),
            licence_files_text(&p.dir)?,
        );
    }
    render(&packages, &texts)
}

/// `cargo xtask licenses`: write the file. With `check`, compare it with a fresh generation and
/// return false (after saying so) when they differ.
pub fn run(root: &Path, check: bool) -> Result<bool> {
    let fresh = generate(root)?;
    let path = root.join(FILE);
    if check {
        let committed = std::fs::read_to_string(&path).unwrap_or_default();
        if committed != fresh {
            eprintln!("{FILE} is out of date: run `cargo xtask licenses` and commit the result");
            return Ok(false);
        }
        return Ok(true);
    }
    std::fs::write(&path, &fresh)?;
    eprintln!("wrote {}", path.display());
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn pkg(name: &str, license: &str, repo: Option<&str>) -> Package {
        Package {
            name: name.into(),
            version: "1.0.0".into(),
            license: license.into(),
            repository: repo.map(str::to_string),
            dir: PathBuf::new(),
        }
    }

    fn texts(items: &[(&str, Option<&str>)]) -> BTreeMap<(String, String), Option<String>> {
        items
            .iter()
            .map(|(n, t)| ((n.to_string(), "1.0.0".to_string()), t.map(str::to_string)))
            .collect()
    }

    #[test]
    fn identical_texts_share_one_entry() {
        let ps = [
            pkg("a", "MIT", None),
            pkg("b", "MIT", None),
            pkg("c", "Apache-2.0", None),
        ];
        let t = texts(&[
            ("a", Some("SAME")),
            ("b", Some("SAME")),
            ("c", Some("OTHER")),
        ]);
        let out = render(&ps, &t).unwrap();
        assert_eq!(out.matches("SAME").count(), 1);
        assert_eq!(out.matches("OTHER").count(), 1);
        assert!(out.contains("a 1.0.0 (MIT)\nb 1.0.0 (MIT)\n"), "{out}");
        assert!(out.contains("c 1.0.0 (Apache-2.0)\n"));
        assert!(out.contains("the licences of the 3 crates"));
    }

    #[test]
    fn a_package_without_a_licence_file_gets_the_standard_text_and_its_repository() {
        let ps = [pkg("n", "MIT OR Apache-2.0", Some("https://example.org/n"))];
        let out = render(&ps, &texts(&[("n", None)])).unwrap();
        assert!(out.contains("n 1.0.0 has no licence file"), "{out}");
        assert!(out.contains("https://example.org/n"));
        assert!(out.contains("Permission is hereby granted, free of charge"));
        assert!(out.contains("Apache License"));
    }

    #[test]
    fn a_copyleft_alternative_is_not_carried_and_an_unknown_licence_is_an_error() {
        let ps = [pkg("r", "MIT OR LGPL-2.1-or-later", None)];
        let out = render(&ps, &texts(&[("r", None)])).unwrap();
        assert!(out.contains("Permission is hereby granted"));
        assert!(!out.contains("GNU Lesser"));
        let e = render(&[pkg("x", "BSD-3-Clause", None)], &texts(&[("x", None)])).unwrap_err();
        assert!(e.to_string().contains("x 1.0.0 (BSD-3-Clause)"), "{e}");
    }

    #[test]
    fn an_mpl_package_names_where_its_source_is() {
        let ps = [pkg(
            "index_list",
            "MPL-2.0",
            Some("https://github.com/o/index_list"),
        )];
        let out = render(&ps, &texts(&[("index_list", Some("MPL TEXT"))])).unwrap();
        assert!(
            out.contains(
                "Source of index_list 1.0.0 (MPL-2.0, used unmodified): \
                 https://crates.io/crates/index_list/1.0.0, repository https://github.com/o/index_list"
            ),
            "{out}"
        );
        let none = render(&[pkg("m", "MIT", None)], &texts(&[("m", Some("T"))])).unwrap();
        assert!(!none.contains("Source of"));
    }

    #[test]
    fn spdx_expressions_split_into_identifiers() {
        assert_eq!(
            spdx_ids("(MIT OR Apache-2.0) AND Unicode-3.0"),
            ["MIT", "Apache-2.0", "Unicode-3.0"]
        );
        assert_eq!(spdx_ids("MIT/Apache-2.0"), ["MIT", "Apache-2.0"]);
        assert_eq!(
            spdx_ids("Apache-2.0 WITH LLVM-exception OR Apache-2.0"),
            ["Apache-2.0", "LLVM-exception"]
        );
    }

    fn entry(id: &str, name: &str, source: Option<&str>) -> Value {
        json!({"id": id, "name": name, "version": "1.0.0", "license": "MIT",
               "source": source, "manifest_path": format!("/src/{name}/Cargo.toml")})
    }

    fn dep(pkg: &str, kind: Option<&str>) -> Value {
        json!({"pkg": pkg, "dep_kinds": [{"kind": kind, "target": null}]})
    }

    #[test]
    fn the_closure_is_the_normal_dependencies_of_the_shipped_artefacts() {
        let reg = Some("registry+x");
        let meta = json!({
            "workspace_members": ["gloci", "gloci-py", "gloci-io", "xtask"],
            "packages": [
                entry("gloci", "gloci", None), entry("gloci-py", "gloci-py", None),
                entry("gloci-io", "gloci-io", None), entry("xtask", "xtask", None),
                entry("serde", "serde", reg), entry("itoa", "itoa", reg),
                entry("cc", "cc", reg), entry("tempfile", "tempfile", reg),
                entry("clap", "clap", reg), entry("pyo3", "pyo3", reg),
            ],
            "resolve": {"nodes": [
                {"id": "gloci", "deps": [dep("gloci-io", None), dep("clap", None),
                                          dep("tempfile", Some("dev")), dep("cc", Some("build"))]},
                {"id": "gloci-py", "deps": [dep("pyo3", None)]},
                {"id": "gloci-io", "deps": [dep("serde", None)]},
                {"id": "xtask", "deps": [dep("gloci", None), dep("itoa", None)]},
                {"id": "serde", "deps": []}, {"id": "itoa", "deps": []}, {"id": "cc", "deps": []},
                {"id": "tempfile", "deps": []}, {"id": "clap", "deps": []}, {"id": "pyo3", "deps": []},
            ]},
        });
        let names: Vec<String> = closure(&meta)
            .unwrap()
            .into_iter()
            .map(|p| p.name)
            .collect();
        // Not itoa (only xtask uses it), tempfile (dev), cc (build) or the workspace members.
        assert_eq!(names, ["clap", "pyo3", "serde"]);
    }

    #[test]
    fn licence_files_are_found_by_prefix_in_any_case() {
        let d = tempfile::tempdir().unwrap();
        for (n, t) in [
            ("LICENSE-MIT", "mit\r\ntext  \n"),
            ("license-apache", "apache"),
            ("NOTICE", "notice"),
            ("Cargo.toml", "x"),
            ("README.md", "y"),
        ] {
            std::fs::write(d.path().join(n), t).unwrap();
        }
        // Upper-case names sort first: LICENSE-MIT, NOTICE, license-apache.
        assert_eq!(
            licence_files_text(d.path()).unwrap().unwrap(),
            "mit\ntext\n\nnotice\n\napache"
        );
        let empty = tempfile::tempdir().unwrap();
        assert_eq!(licence_files_text(empty.path()).unwrap(), None);
    }
}
