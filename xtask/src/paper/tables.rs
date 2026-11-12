//! `cargo xtask paper`: the CSV tables T1 to T8 and the LaTeX numbers file
//! (D98), computed from one results directory and the committed repository
//! files. Output is a pure function of those inputs.
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use anyhow::{Context, Result};
use serde_json::Value;

use super::data::{Results, csv_line};
use crate::bench::{Case, eval, medians};
use crate::diff::workspace_root;

pub const HEADERS: [(&str, &[&str]); 8] = [
    (
        "T1-datasets",
        &[
            "dataset",
            "description",
            "release",
            "url",
            "checksum_algorithm",
            "checksum",
            "bytes",
            "records",
        ],
    ),
    (
        "T2-benchmarks",
        &[
            "case",
            "label",
            "dataset",
            "tool",
            "role",
            "threads",
            "median_s",
            "mean_s",
            "stddev_s",
            "min_s",
            "max_s",
            "runs",
            "max_rss_mb",
        ],
    ),
    (
        "T3-speedups",
        &[
            "case",
            "label",
            "baseline",
            "speedup",
            "target",
            "target_met",
        ],
    ),
    (
        "T4-concordance",
        &[
            "case",
            "command",
            "reference_tool",
            "dataset",
            "mode",
            "n_gloci",
            "n_reference",
            "n_matching",
            "match_rate",
            "anchors",
        ],
    ),
    (
        "T5-differences",
        &["anchor", "title", "golden_cases", "concordance_cases"],
    ),
    (
        "T6-command-coverage",
        &["tool", "old_command", "gloci_command", "status"],
    ),
    ("T7-tests", &["category", "count"]),
    ("T8-software", &["tool", "version", "role"]),
];

/// Seconds `%.3f`, speedup `%.2f`, rate `%.6f`, MB `%.1f`.
pub fn secs(x: f64) -> String {
    format!("{x:.3}")
}
pub fn ratio(x: f64) -> String {
    format!("{x:.2}")
}
pub fn rate(x: f64) -> String {
    format!("{x:.6}")
}
pub fn mb(kb: u64) -> String {
    format!("{:.1}", kb as f64 / 1024.0)
}

fn opt<T>(x: Option<T>, f: impl Fn(T) -> String) -> String {
    x.map(f).unwrap_or_default()
}

fn text(v: &Value) -> String {
    match v {
        Value::Null => String::new(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// The paper's benchmark cases by name; empty when `bench.toml` is unreadable
/// is an error, since roles and targets come from it.
fn paper_cases() -> Result<BTreeMap<String, Case>> {
    Ok(crate::paper::bench::load_paper_cases()?
        .into_iter()
        .map(|c| (c.name.clone(), c))
        .collect())
}

/// A tool is a baseline when a target of its case names it.
fn role(case: Option<&Case>, tool: &str) -> &'static str {
    if tool == "gloci" {
        return "subject";
    }
    let named = case.is_some_and(|c| {
        c.targets
            .iter()
            .any(|t| t.split_whitespace().any(|w| w == tool))
    });
    if named { "baseline" } else { "yardstick" }
}

fn sorted(mut rows: Vec<Vec<String>>) -> Vec<Vec<String>> {
    rows.sort();
    rows
}

fn t1(r: &Results) -> Result<Vec<Vec<String>>> {
    let cat = super::fetch::load_catalog()?;
    let info = &r.snapshot["datasets"];
    let mut rows = Vec::new();
    for d in cat.datasets.iter().filter(|d| d.paper) {
        let (alg, sum) = d.checksum()?;
        rows.push(vec![
            d.name.clone(),
            d.description.clone().unwrap_or_default(),
            d.release.clone().unwrap_or_default(),
            d.url.clone(),
            alg.to_string(),
            sum.to_string(),
            text(&info[&d.name]["bytes"]),
            text(&info[&d.name]["records"]),
        ]);
    }
    Ok(sorted(rows))
}

fn dataset_of(case: Option<&Case>) -> String {
    case.and_then(|c| c.inputs.get("a").cloned())
        .or_else(|| case.and_then(|c| c.dataset.clone()))
        .unwrap_or_default()
}

fn t2(r: &Results, cases: &BTreeMap<String, Case>) -> Vec<Vec<String>> {
    let mut rows = Vec::new();
    for b in &r.bench {
        let c = cases.get(&b.case);
        let label = b
            .label
            .clone()
            .or_else(|| c.and_then(|c| c.label.clone()))
            .unwrap_or_else(|| b.case.clone());
        for t in &b.tools {
            rows.push(vec![
                b.case.clone(),
                label.clone(),
                dataset_of(c),
                t.tool.clone(),
                role(c, &t.tool).to_string(),
                opt(
                    b.threads
                        .or_else(|| c.and_then(|c| c.threads.map(u64::from))),
                    |x| x.to_string(),
                ),
                secs(t.median),
                opt(t.mean, secs),
                opt(t.stddev, secs),
                opt(t.min, secs),
                opt(t.max, secs),
                opt(t.runs, |x| x.to_string()),
                opt(t.rss_kb, mb),
            ]);
        }
    }
    sorted(rows)
}

fn t3(r: &Results, cases: &BTreeMap<String, Case>) -> Vec<Vec<String>> {
    let mut rows = Vec::new();
    for b in &r.bench {
        let Some(g) = b.tools.iter().find(|t| t.tool == "gloci") else {
            continue;
        };
        let c = cases.get(&b.case);
        let label = b
            .label
            .clone()
            .or_else(|| c.and_then(|c| c.label.clone()))
            .unwrap_or_else(|| b.case.clone());
        let pairs: Vec<(&str, f64)> = b
            .tools
            .iter()
            .map(|t| (t.tool.as_str(), t.median))
            .collect();
        let m = medians(&pairs);
        for t in b.tools.iter().filter(|t| t.tool != "gloci") {
            let named: Vec<&String> = c
                .map(|c| {
                    c.targets
                        .iter()
                        .filter(|x| x.split_whitespace().any(|w| w == t.tool))
                        .collect()
                })
                .unwrap_or_default();
            let met: Vec<bool> = named.iter().filter_map(|x| eval(x, &m).ok()).collect();
            rows.push(vec![
                b.case.clone(),
                label.clone(),
                t.tool.clone(),
                ratio(t.median / g.median),
                named
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join("; "),
                if met.is_empty() {
                    String::new()
                } else {
                    met.iter().all(|&x| x).to_string()
                },
            ]);
        }
    }
    sorted(rows)
}

fn t4(r: &Results) -> Vec<Vec<String>> {
    sorted(
        r.concordance
            .iter()
            .map(|c| {
                vec![
                    c.case.clone(),
                    c.command.clone(),
                    c.reference_tool.clone(),
                    c.dataset.clone(),
                    c.mode.clone(),
                    c.n_gloci.to_string(),
                    c.n_reference.to_string(),
                    c.n_matching.to_string(),
                    rate(c.match_rate),
                    c.anchors.join(";"),
                ]
            })
            .collect(),
    )
}

/// Distinct concordance cases that name `anchor`.
fn concordance_cases(r: &Results, anchor: &str) -> usize {
    r.concordance
        .iter()
        .filter(|c| c.anchors.iter().any(|a| a == anchor))
        .map(|c| c.case.as_str())
        .collect::<BTreeSet<_>>()
        .len()
}

fn t5(r: &Results) -> Vec<Vec<String>> {
    sorted(
        r.snapshot["anchors"]
            .as_array()
            .map(|a| a.as_slice())
            .unwrap_or_default()
            .iter()
            .map(|a| {
                let anchor = text(&a["anchor"]);
                vec![
                    anchor.clone(),
                    text(&a["title"]),
                    text(&a["golden_cases"]),
                    concordance_cases(r, &anchor).to_string(),
                ]
            })
            .collect(),
    )
}

fn t6(r: &Results) -> Vec<Vec<String>> {
    sorted(
        r.snapshot["migration"]
            .as_array()
            .map(|a| a.as_slice())
            .unwrap_or_default()
            .iter()
            .map(|m| {
                ["tool", "old_command", "gloci_command", "status"]
                    .iter()
                    .map(|k| text(&m[k]))
                    .collect()
            })
            .collect(),
    )
}

fn t7(r: &Results) -> Vec<Vec<String>> {
    sorted(
        r.snapshot["tests"]
            .as_object()
            .into_iter()
            .flatten()
            .map(|(k, v)| vec![k.clone(), text(v)])
            .collect(),
    )
}

/// A tool is a baseline when a target of any case names it.
fn software_role(cases: &BTreeMap<String, Case>, tool: &str) -> &'static str {
    match tool {
        "gloci" => "subject",
        "hyperfine" => "harness",
        "bgzip" | "tabix" => "helper",
        _ if cases.values().any(|c| role(Some(c), tool) == "baseline") => "baseline",
        _ => "yardstick",
    }
}

fn t8(r: &Results, cases: &BTreeMap<String, Case>) -> Vec<Vec<String>> {
    sorted(
        r.environment["tools"]
            .as_object()
            .into_iter()
            .flatten()
            .map(|(k, v)| vec![k.clone(), text(v), software_role(cases, k).to_string()])
            .collect(),
    )
}

/// Every output file, by path relative to the output root: the tables under
/// `tables/`, then `numbers.tex` and `numbers.json`.
pub fn render(results: &Path) -> Result<BTreeMap<String, String>> {
    let r = super::data::load(results)?;
    let cases = paper_cases().context("loading the paper cases (xtask/paper/bench.toml)")?;
    let bodies = [
        t1(&r)?,
        t2(&r, &cases),
        t3(&r, &cases),
        t4(&r),
        t5(&r),
        t6(&r),
        t7(&r),
        t8(&r, &cases),
    ];
    let mut out = BTreeMap::new();
    for ((name, header), rows) in HEADERS.iter().zip(bodies) {
        let mut s = csv_line(header);
        s.push('\n');
        for row in rows {
            s.push_str(&csv_line(&row));
            s.push('\n');
        }
        out.insert(format!("tables/{name}.csv"), s);
    }
    let map = super::numbers::numbers_of(&r)?;
    out.insert("numbers.tex".into(), super::numbers::numbers_tex(&map));
    out.insert(
        "numbers.json".into(),
        serde_json::to_string_pretty(&map)? + "\n",
    );
    Ok(out)
}

/// Write the tables into `out` (e.g. `paper/tables`) and `numbers.tex` and
/// `numbers.json` beside it.
pub fn tables(results: &Path, out: &Path) -> Result<()> {
    let parent = out.parent().filter(|p| !p.as_os_str().is_empty());
    let root = parent.unwrap_or(Path::new("."));
    for (rel, body) in render(results)? {
        let path = match rel.strip_prefix("tables/") {
            Some(name) => out.join(name),
            None => root.join(&rel),
        };
        if let Some(d) = path.parent() {
            std::fs::create_dir_all(d)?;
        }
        std::fs::write(&path, body).with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

/// The newest `paper/results/<version>` directory.
pub fn latest_results() -> Result<std::path::PathBuf> {
    let dir = workspace_root().join("paper/results");
    let mut v: Vec<_> = std::fs::read_dir(&dir)
        .with_context(|| {
            format!(
                "{} does not exist: run `cargo xtask paper run` on the reference machine first",
                dir.display()
            )
        })?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.is_dir())
        .collect();
    v.sort();
    v.pop().with_context(|| {
        format!(
            "{} is empty: run `cargo xtask paper run` first",
            dir.display()
        )
    })
}
