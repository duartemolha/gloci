//! Release metadata checks (D81): one version, a CHANGELOG section for it,
//! a matching CITATION.cff, a LICENSE, and consistent internal dependencies.
use anyhow::{Result, anyhow, bail};
use std::path::Path;

/// One `## ` section of a Keep a Changelog file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Section {
    /// `None` for `[Unreleased]`.
    pub version: Option<String>,
    pub date: Option<String>,
    /// The lines under the heading, without leading or trailing blank lines.
    pub body: String,
}

fn parse_version(s: &str) -> Option<(u64, u64, u64)> {
    let mut it = s.split('.');
    let v = (
        it.next()?.parse().ok()?,
        it.next()?.parse().ok()?,
        it.next()?.parse().ok()?,
    );
    it.next().is_none().then_some(v)
}

fn is_link_def(l: &str) -> bool {
    l.strip_prefix('[')
        .and_then(|r| r.split_once("]:"))
        .is_some_and(|(name, rest)| {
            !name.is_empty() && !name.contains(']') && rest.starts_with(char::is_whitespace)
        })
}

fn is_iso_date(s: &str) -> bool {
    let b = s.as_bytes();
    if b.len() != 10 || b[4] != b'-' || b[7] != b'-' {
        return false;
    }
    let digits = |r: std::ops::Range<usize>| s[r].bytes().all(|c| c.is_ascii_digit());
    if !(digits(0..4) && digits(5..7) && digits(8..10)) {
        return false;
    }
    let (m, d): (u32, u32) = (s[5..7].parse().unwrap(), s[8..10].parse().unwrap());
    let y: u32 = s[0..4].parse().unwrap();
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    let max = match m {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if leap => 29,
        2 => 28,
        _ => return false,
    };
    (1..=max).contains(&d)
}

/// Parse the `## [..]` sections. `[Unreleased]` must come first, versions must
/// strictly decrease, and each date must be ISO (YYYY-MM-DD).
pub fn changelog_sections(text: &str) -> Result<Vec<Section>> {
    let mut out: Vec<(Option<String>, Option<String>, Vec<&str>)> = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if let Some(h) = line.strip_prefix("## ") {
            let h = h.trim();
            let rest = h.strip_prefix('[').and_then(|r| r.split_once(']')).ok_or_else(|| {
                anyhow!("CHANGELOG.md line {}: heading `{h}` must be `## [Unreleased]` or `## [X.Y.Z] - YYYY-MM-DD`", i + 1)
            })?;
            let (name, tail) = rest;
            if name == "Unreleased" {
                if !out.is_empty() {
                    bail!(
                        "CHANGELOG.md line {}: `[Unreleased]` must be the first section",
                        i + 1
                    );
                }
                if !tail.trim().is_empty() {
                    bail!("CHANGELOG.md line {}: `[Unreleased]` takes no date", i + 1);
                }
                out.push((None, None, Vec::new()));
                continue;
            }
            let ver = parse_version(name).ok_or_else(|| {
                anyhow!(
                    "CHANGELOG.md line {}: `{name}` is not a version X.Y.Z",
                    i + 1
                )
            })?;
            let date = tail
                .trim()
                .strip_prefix('-')
                .map(str::trim)
                .ok_or_else(|| {
                    anyhow!(
                        "CHANGELOG.md line {}: heading `{h}` needs ` - YYYY-MM-DD`",
                        i + 1
                    )
                })?;
            if !is_iso_date(date) {
                bail!(
                    "CHANGELOG.md line {}: date `{date}` is not ISO (YYYY-MM-DD)",
                    i + 1
                );
            }
            if let Some((Some(prev), _, _)) = out.last()
                && parse_version(prev).is_some_and(|p| p <= ver)
            {
                bail!(
                    "CHANGELOG.md line {}: version {name} must be lower than the section above ({prev}); newest first",
                    i + 1
                );
            }
            out.push((Some(name.to_string()), Some(date.to_string()), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.2.push(line);
        }
    }
    Ok(out
        .into_iter()
        .map(|(version, date, mut lines)| {
            // Keep a Changelog footer links (`[x]: https://...`) are not notes.
            while lines
                .last()
                .is_some_and(|l| l.trim().is_empty() || is_link_def(l))
            {
                lines.pop();
            }
            let start = lines.iter().position(|l| !l.trim().is_empty());
            let end = lines.iter().rposition(|l| !l.trim().is_empty());
            let body = match (start, end) {
                (Some(s), Some(e)) => lines[s..=e].join("\n"),
                _ => String::new(),
            };
            Section {
                version,
                date,
                body,
            }
        })
        .collect())
}

/// The body of `version`'s section, as the release notes.
pub fn release_notes(changelog: &str, version: &str) -> Result<String> {
    let sections = changelog_sections(changelog)?;
    sections
        .into_iter()
        .find(|s| s.version.as_deref() == Some(version))
        .map(|s| s.body)
        .ok_or_else(|| {
            anyhow!("CHANGELOG.md has no `## [{version}] - YYYY-MM-DD` section; add one before releasing")
        })
}

fn cff_field(cff: &str, key: &str) -> Option<String> {
    cff.lines().find_map(|l| {
        let v = l.strip_prefix(key)?.strip_prefix(':')?;
        // A YAML comment starts at a `#` after whitespace (outside quotes).
        let v = v.trim();
        let v = if v.starts_with(['"', '\'']) {
            v
        } else {
            v.split(" #").next().unwrap_or(v)
        };
        Some(v.trim().trim_matches(['"', '\'']).to_string())
    })
}

fn read_toml(path: &Path) -> Result<toml::Value> {
    let s = std::fs::read_to_string(path).map_err(|e| anyhow!("{}: {e}", path.display()))?;
    toml::from_str::<toml::Value>(&s).map_err(|e| anyhow!("{}: {e}", path.display()))
}

/// Internal dependencies (a `path` and a `version`) in one manifest whose
/// version differs from `version`.
fn dep_problems(manifest: &toml::Value, label: &str, version: &str, out: &mut Vec<String>) {
    let mut tables: Vec<(&str, &toml::Value)> = Vec::new();
    for k in ["dependencies", "dev-dependencies", "build-dependencies"] {
        if let Some(t) = manifest.get(k) {
            tables.push((k, t));
        }
    }
    if let Some(t) = manifest
        .get("workspace")
        .and_then(|w| w.get("dependencies"))
    {
        tables.push(("workspace.dependencies", t));
    }
    for (kind, t) in tables {
        let Some(t) = t.as_table() else { continue };
        for (name, dep) in t {
            let (Some(p), Some(v)) = (dep.get("path"), dep.get("version")) else {
                continue;
            };
            let _ = p;
            if v.as_str() != Some(version) {
                out.push(format!(
                    "{label}: [{kind}] `{name}` has version {v} but the workspace version is {version}; set it to \"{version}\""
                ));
            }
        }
    }
}

/// Check the repository's release metadata. Each problem says how to fix it.
pub fn release_check(root: &Path, tag: Option<&str>) -> Result<Vec<String>> {
    let mut problems = Vec::new();
    let ws = read_toml(&root.join("Cargo.toml"))?;
    let version = ws
        .get("workspace")
        .and_then(|w| w.get("package"))
        .and_then(|p| p.get("version"))
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow!("Cargo.toml has no [workspace.package] version"))?
        .to_string();

    if let Some(tag) = tag
        && tag != format!("v{version}")
    {
        problems.push(format!(
            "tag {tag} does not match the workspace version {version}; tag v{version} or bump the version"
        ));
    }

    // CHANGELOG
    let mut date: Option<String> = None;
    match std::fs::read_to_string(root.join("CHANGELOG.md")) {
        Err(e) => problems.push(format!("CHANGELOG.md: {e}; create it (Keep a Changelog format)")),
        Ok(text) => match changelog_sections(&text) {
            Err(e) => problems.push(format!("{e:#}")),
            Ok(sections) => match sections.iter().find(|s| s.version.as_deref() == Some(&version)) {
                None => problems.push(format!(
                    "CHANGELOG.md has no `## [{version}] - YYYY-MM-DD` section; add one for this release"
                )),
                Some(s) => {
                    date = s.date.clone();
                    if !s.body.lines().any(|l| {
                        let l = l.trim_start();
                        l.starts_with("- ") || l.starts_with("* ")
                    }) {
                        problems.push(format!(
                            "CHANGELOG.md section [{version}] has no bullet; list what the release contains"
                        ));
                    }
                }
            },
        },
    }

    // CITATION.cff
    match std::fs::read_to_string(root.join("CITATION.cff")) {
        Err(e) => problems.push(format!("CITATION.cff: {e}; create it")),
        Ok(cff) => {
            match cff_field(&cff, "version") {
                Some(v) if v == version => {}
                other => problems.push(format!(
                    "CITATION.cff version is {} but the workspace version is {version}; set `version: {version}`",
                    other.map_or("missing".into(), |v| format!("`{v}`"))
                )),
            }
            if let Some(d) = &date {
                match cff_field(&cff, "date-released") {
                    Some(v) if &v == d => {}
                    other => problems.push(format!(
                        "CITATION.cff date-released is {} but CHANGELOG.md dates [{version}] {d}; set `date-released: {d}`",
                        other.map_or("missing".into(), |v| format!("`{v}`"))
                    )),
                }
            }
        }
    }

    // LICENSE
    if !root.join("LICENSE").is_file() {
        problems.push("LICENSE is missing; add the MIT licence text as LICENSE".to_string());
    }

    // Internal dependency versions
    dep_problems(&ws, "Cargo.toml", &version, &mut problems);
    let mut manifests = vec![root.join("xtask/Cargo.toml")];
    if let Ok(rd) = std::fs::read_dir(root.join("crates")) {
        let mut dirs: Vec<_> = rd.filter_map(|e| e.ok()).map(|e| e.path()).collect();
        dirs.sort();
        manifests.extend(dirs.into_iter().map(|d| d.join("Cargo.toml")));
    }
    for m in manifests.into_iter().filter(|m| m.is_file()) {
        let label = m.strip_prefix(root).unwrap_or(&m).display().to_string();
        dep_problems(&read_toml(&m)?, &label, &version, &mut problems);
    }
    Ok(problems)
}
