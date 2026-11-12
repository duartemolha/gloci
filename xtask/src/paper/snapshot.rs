//! `repo-snapshot.json`: what the repository itself says at the commit of a
//! paper run: golden cases per tool, differences anchors, the migration
//! tables, test counts and the size of every dataset.
use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use serde_json::{Value, json};

use crate::bench::{Catalog, bench_cache_dir, input_path};

/// `(anchor, title)` of every `## Title {#anchor}` heading.
pub fn parse_anchors(md: &str) -> Vec<(String, String)> {
    md.lines()
        .filter_map(|l| {
            let h = l.strip_prefix("## ")?;
            let (title, rest) = h.rsplit_once(" {#")?;
            Some((
                rest.strip_suffix('}')?.to_string(),
                title.trim().to_string(),
            ))
        })
        .collect()
}

fn cells(line: &str) -> Vec<String> {
    let inner = line.trim().trim_start_matches('|').trim_end_matches('|');
    let mut out = vec![String::new()];
    let mut prev = ' ';
    for c in inner.chars() {
        if c == '|' && prev != '\\' {
            out.push(String::new());
        } else if let Some(cur) = out.last_mut() {
            cur.push(c);
        }
        prev = c;
    }
    out.iter().map(|c| c.trim().replace("\\|", "|")).collect()
}

fn tool_of_heading(h: &str) -> Option<&'static str> {
    let h = h.to_lowercase();
    [
        ("sorting, indexing", "htslib"),
        ("bedtools", "bedtools"),
        ("samtools", "samtools"),
        ("bedops", "bedops"),
        ("vcftools", "vcftools"),
        ("filter expressions", "bcftools"),
    ]
    .iter()
    .find(|(k, _)| h.starts_with(k))
    .map(|(_, t)| *t)
}

/// `(tool, old command, gloci command, status)` of every row of the tables
/// under "Command mapping" in docs/migration.md. Status is `unavailable` when
/// the gloci cell starts with "not available", `partial` when it says so
/// further on, else `available`.
pub fn parse_migration(md: &str) -> Vec<[String; 4]> {
    let mut rows = Vec::new();
    let mut in_mapping = false;
    let mut tool: Option<&str> = None;
    let mut table_line = 0usize;
    for l in md.lines() {
        if let Some(h) = l.strip_prefix("## ") {
            in_mapping = h.trim() == "Command mapping";
            tool = None;
            table_line = 0;
            continue;
        }
        if let Some(h) = l.strip_prefix("### ") {
            tool = tool_of_heading(h);
            table_line = 0;
            continue;
        }
        if !l.starts_with('|') {
            table_line = 0;
            continue;
        }
        table_line += 1;
        let (true, Some(t)) = (in_mapping, tool) else {
            continue;
        };
        // Line 1 is the header, line 2 the `---` rule.
        if table_line <= 2 {
            continue;
        }
        let c = cells(l);
        if c.len() < 2 {
            continue;
        }
        let g = c[1].to_lowercase();
        let status = if g.starts_with("not available") {
            "unavailable"
        } else if g.contains("not available") {
            "partial"
        } else {
            "available"
        };
        rows.push([
            t.to_string(),
            c[0].clone(),
            c[1].clone(),
            status.to_string(),
        ]);
    }
    rows
}

fn count_lines_matching(path: &Path, pat: &str) -> u64 {
    std::fs::read_to_string(path)
        .map(|t| t.lines().filter(|l| l.trim() == pat).count() as u64)
        .unwrap_or(0)
}

fn walk(dir: &Path, ext: &str, f: &mut dyn FnMut(&Path)) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    let mut v: Vec<_> = rd.filter_map(|e| e.ok().map(|e| e.path())).collect();
    v.sort();
    for p in v {
        if p.is_dir() {
            if p.file_name().is_some_and(|n| n == "target" || n == ".git") {
                continue;
            }
            walk(&p, ext, f);
        } else if p.extension().is_some_and(|e| e == ext) {
            f(&p);
        }
    }
}

/// Lines of a text file, skipping `#` comments and BED headers; BGZF/gzip
/// when the name says so. None for formats gloci's table does not count.
pub fn count_records(path: &Path) -> Option<u64> {
    let name = path.file_name()?.to_string_lossy().to_lowercase();
    let (inner, gz) = match name
        .strip_suffix(".gz")
        .or_else(|| name.strip_suffix(".bgz"))
    {
        Some(n) => (n.to_string(), true),
        None => (name.clone(), false),
    };
    if !["bed", "bedgraph", "vcf", "gtf"]
        .iter()
        .any(|e| inner.ends_with(&format!(".{e}")))
    {
        return None;
    }
    let f = std::fs::File::open(path).ok()?;
    let mut r: Box<dyn BufRead> = if gz {
        Box::new(std::io::BufReader::with_capacity(
            1 << 20,
            flate2::read::MultiGzDecoder::new(f),
        ))
    } else {
        Box::new(std::io::BufReader::with_capacity(1 << 20, f))
    };
    let (mut n, mut line) = (0u64, Vec::new());
    loop {
        line.clear();
        match r.read_until(b'\n', &mut line) {
            Ok(0) => return Some(n),
            Ok(_) => {
                if !(line.starts_with(b"#")
                    || line.starts_with(b"track")
                    || line.starts_with(b"browser")
                    || line == b"\n")
                {
                    n += 1;
                }
            }
            Err(_) => return None,
        }
    }
}

/// The snapshot of the repository at `root`, with dataset sizes read from
/// the bench cache under `data_dir`.
pub fn snapshot(
    root: &Path,
    cat: &Catalog,
    data_dir: &Path,
    version: &str,
    commit: Option<&str>,
) -> Result<Value> {
    // Golden cases: one directory each; the tool is the `-vs-<tool>` suffix.
    let golden_dir = root.join("crates/gloci/tests/golden");
    let mut by_tool: BTreeMap<String, u64> = BTreeMap::new();
    let mut golden_names = Vec::new();
    for e in std::fs::read_dir(&golden_dir)
        .with_context(|| format!("reading {}", golden_dir.display()))?
    {
        let e = e?;
        if !e.path().is_dir() {
            continue;
        }
        let name = e.file_name().to_string_lossy().into_owned();
        let tool = name.rsplit_once("-vs-").map_or("other", |(_, t)| t);
        *by_tool.entry(tool.to_string()).or_default() += 1;
        golden_names.push(name);
    }
    golden_names.sort();

    // Anchors, with the golden cases whose differential case names them.
    let cases_text = std::fs::read_to_string(root.join("xtask/diff/cases.toml"))?;
    let cases: toml::Value = toml::from_str(&cases_text)?;
    let case_list = cases
        .get("case")
        .and_then(|c| c.as_array())
        .cloned()
        .unwrap_or_default();
    let mut per_anchor: BTreeMap<String, u64> = BTreeMap::new();
    for c in &case_list {
        let (Some(name), Some(a)) = (
            c.get("name").and_then(|v| v.as_str()),
            c.get("expected_difference").and_then(|v| v.as_str()),
        ) else {
            continue;
        };
        if golden_names.binary_search(&name.to_string()).is_ok() {
            *per_anchor
                .entry(a.trim_start_matches('#').to_string())
                .or_default() += 1;
        }
    }
    let anchors: Vec<Value> = parse_anchors(&std::fs::read_to_string(root.join("docs/differences.md"))?)
        .into_iter()
        .map(|(a, title)| {
            json!({"anchor": a, "title": title, "golden_cases": per_anchor.get(&a).copied().unwrap_or(0)})
        })
        .collect();

    let migration: Vec<Value> = parse_migration(&std::fs::read_to_string(root.join("docs/migration.md"))?)
        .into_iter()
        .map(|[tool, old, new, status]| {
            json!({"tool": tool, "old_command": old, "gloci_command": new, "status": status})
        })
        .collect();

    // Tests: `#[test]` functions by place, plus the data-driven cases.
    let (mut crates, mut xt, mut py) = (0u64, 0u64, 0u64);
    walk(&root.join("crates"), "rs", &mut |p| {
        crates += count_lines_matching(p, "#[test]")
    });
    walk(&root.join("xtask"), "rs", &mut |p| {
        xt += count_lines_matching(p, "#[test]")
    });
    walk(&root.join("python"), "py", &mut |p| {
        py += std::fs::read_to_string(p)
            .map(|t| {
                t.lines()
                    .filter(|l| l.trim_start().starts_with("def test_"))
                    .count() as u64
            })
            .unwrap_or(0)
    });
    let tests = json!({
        "differential-cases": case_list.len(),
        "golden-cases": golden_names.len(),
        "python": py,
        "rust-crates": crates,
        "rust-xtask": xt,
    });

    // Datasets and derived inputs: size and record count when present.
    let cache = bench_cache_dir(data_dir);
    let mut datasets = serde_json::Map::new();
    let names = cat
        .datasets
        .iter()
        .filter(|d| d.paper)
        .map(|d| d.name.as_str())
        .chain(cat.derived.iter().map(|d| d.name.as_str()));
    for name in names {
        let path = input_path(cat, name, &cache)?;
        let bytes = std::fs::metadata(&path)
            .ok()
            .filter(|m| m.is_file())
            .map(|m| m.len());
        let records = bytes.and_then(|_| count_records(&path));
        datasets.insert(
            name.to_string(),
            json!({"bytes": bytes, "records": records}),
        );
    }

    Ok(json!({
        "gloci_version": version,
        "commit": commit,
        "golden": {"total": golden_names.len(), "by_tool": by_tool},
        "anchors": anchors,
        "migration": migration,
        "tests": tests,
        "datasets": datasets,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diff::workspace_root;

    #[test]
    fn migration_tables_parse() {
        let md = std::fs::read_to_string(workspace_root().join("docs/migration.md")).unwrap();
        let rows = parse_migration(&md);
        for t in [
            "htslib", "bedtools", "samtools", "bedops", "vcftools", "bcftools",
        ] {
            assert!(rows.iter().any(|r| r[0] == t), "no rows for {t}");
        }
        assert!(rows.len() > 80, "{} rows", rows.len());
        assert!(rows.iter().all(|r| !r[1].is_empty() && !r[2].is_empty()));
    }

    #[test]
    fn anchors_parse() {
        let md = std::fs::read_to_string(workspace_root().join("docs/differences.md")).unwrap();
        let a = parse_anchors(&md);
        assert!(
            a.iter()
                .any(|(k, t)| k == "sort-tie-order" && !t.is_empty())
        );
        assert!(a.len() > 100);
    }

    #[test]
    fn cells_split_on_unescaped_pipes() {
        assert_eq!(cells("| a \\| b | c |"), vec!["a | b", "c"]);
    }
}
