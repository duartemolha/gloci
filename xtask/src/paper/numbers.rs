//! The paper's numbers (D98): every figure the text quotes, keyed, and the
//! `numbers.tex` that defines one LaTeX macro value per key.
use std::collections::BTreeMap;
use std::path::Path;

use anyhow::Result;
use serde_json::{Number, Value};

use super::data::Results;
use super::tables::{mb, ratio, secs};

fn int(n: u64) -> Value {
    Value::from(n)
}

/// A float rounded exactly as its display is, so JSON and TeX agree.
fn float(shown: String) -> Value {
    Value::Number(Number::from_f64(shown.parse().unwrap_or(0.0)).unwrap_or_else(|| 0.into()))
}

pub fn numbers(results: &Path) -> Result<BTreeMap<String, Value>> {
    numbers_of(&super::data::load(results)?)
}

pub fn numbers_of(r: &Results) -> Result<BTreeMap<String, Value>> {
    let mut m = BTreeMap::new();
    for b in &r.bench {
        for t in &b.tools {
            let k = format!("bench.{}.{}", b.case, t.tool);
            m.insert(format!("{k}.median_s"), float(secs(t.median)));
            if let Some(kb) = t.rss_kb {
                m.insert(format!("{k}.max_rss_mb"), float(mb(kb)));
            }
        }
        if let Some(g) = b.tools.iter().find(|t| t.tool == "gloci") {
            for t in b.tools.iter().filter(|t| t.tool != "gloci") {
                let s = t.median / g.median;
                m.insert(
                    format!("bench.{}.speedup.{}", b.case, t.tool),
                    float(ratio(s)),
                );
            }
        }
    }
    if !r.concordance.is_empty() {
        let cases: std::collections::BTreeSet<&str> =
            r.concordance.iter().map(|c| c.case.as_str()).collect();
        m.insert("concordance.cases".into(), int(cases.len() as u64));
        let at_100 = r
            .concordance
            .iter()
            .filter(|c| c.mode == "compat" && c.match_rate.is_some_and(|r| r >= 1.0))
            .count();
        m.insert("concordance.compat.at_100".into(), int(at_100 as u64));
        let def: Vec<f64> = r
            .concordance
            .iter()
            .filter(|c| c.mode == "default")
            .filter_map(|c| c.match_rate)
            .collect();
        if !def.is_empty() {
            let mean = def.iter().sum::<f64>() / def.len() as f64;
            m.insert(
                "concordance.default.mean_match_rate".into(),
                float(format!("{mean:.6}")),
            );
        }
        // A failed run has no rate: no number for it.
        for c in &r.concordance {
            if let Some(rate) = c.match_rate {
                m.insert(
                    format!("concordance.{}.{}.match_rate", c.case, c.mode),
                    float(format!("{rate:.6}")),
                );
            }
        }
    }
    let golden = &r.snapshot["golden"];
    if let Some(n) = golden["total"].as_u64() {
        m.insert("tests.golden_cases".into(), int(n));
    }
    for (tool, n) in golden["by_tool"].as_object().into_iter().flatten() {
        if let Some(n) = n.as_u64() {
            m.insert(format!("tests.golden.{tool}"), int(n));
        }
    }
    if let Some(a) = r.snapshot["anchors"].as_array() {
        m.insert("differences.anchors".into(), int(a.len() as u64));
    }
    for (name, d) in r.snapshot["datasets"].as_object().into_iter().flatten() {
        if let Some(n) = d["records"].as_u64() {
            m.insert(format!("datasets.{name}.records"), int(n));
        }
    }
    let mut mig: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    for row in r.snapshot["migration"].as_array().into_iter().flatten() {
        let e = mig
            .entry(row["tool"].as_str().unwrap_or("").to_string())
            .or_default();
        e.1 += 1;
        if row["status"].as_str() != Some("unavailable") {
            e.0 += 1;
        }
    }
    for (tool, (avail, total)) in mig {
        m.insert(format!("migration.{tool}.available"), int(avail));
        m.insert(format!("migration.{tool}.total"), int(total));
    }
    Ok(m)
}

/// The text of one value, in the format its key calls for.
pub fn display(key: &str, v: &Value) -> String {
    match v.as_f64() {
        Some(x) if v.is_f64() => {
            let places = if key.ends_with(".median_s") {
                3
            } else if key.ends_with(".max_rss_mb") {
                1
            } else if key.contains(".speedup.") {
                2
            } else {
                6
            };
            format!("{x:.places$}")
        }
        _ => v.to_string(),
    }
}

pub fn numbers_tex(m: &BTreeMap<String, Value>) -> String {
    let mut lines: Vec<String> = m
        .iter()
        .map(|(k, v)| {
            format!(
                "\\expandafter\\def\\csname glnum@{k}\\endcsname{{{}}}\n",
                display(k, v)
            )
        })
        .collect();
    lines.sort();
    let mut s = String::from("\\providecommand{\\glnum}[1]{\\csname glnum@#1\\endcsname}\n");
    s.extend(lines);
    s
}
