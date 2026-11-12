//! What the paper tables read: one results directory
//! (`paper/results/<version>/`) with its benchmark JSONs, `concordance.csv`,
//! `environment.json` and `repo-snapshot.json`.
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail, ensure};
use serde_json::Value;

/// The columns `concordance.csv` has (written by `paper run`).
pub const CONCORDANCE_COLUMNS: [&str; 11] = [
    "case",
    "command",
    "reference_tool",
    "dataset",
    "mode",
    "n_gloci",
    "n_reference",
    "n_matching",
    "match_rate",
    "verdict",
    "anchors",
];

pub struct ToolRow {
    pub tool: String,
    pub median: f64,
    pub mean: Option<f64>,
    pub stddev: Option<f64>,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub runs: Option<u64>,
    /// Peak RSS in KiB; none when the measurement failed or is missing.
    pub rss_kb: Option<u64>,
}

pub struct BenchCase {
    pub case: String,
    pub label: Option<String>,
    pub threads: Option<u64>,
    pub tools: Vec<ToolRow>,
}

pub struct ConcRow {
    pub case: String,
    pub command: String,
    pub reference_tool: String,
    pub dataset: String,
    pub mode: String,
    pub n_gloci: u64,
    pub n_reference: u64,
    pub n_matching: u64,
    pub match_rate: f64,
    pub anchors: Vec<String>,
}

pub struct Results {
    pub version: String,
    pub bench: Vec<BenchCase>,
    pub concordance: Vec<ConcRow>,
    pub environment: Value,
    pub snapshot: Value,
}

fn read_json(path: &Path) -> Result<Value> {
    let text = std::fs::read_to_string(path).with_context(|| {
        format!(
            "reading {} (written by `cargo xtask paper run`)",
            path.display()
        )
    })?;
    serde_json::from_str(&text).with_context(|| format!("parsing {}", path.display()))
}

/// A hyperfine statistic: flat (1.x) or under `summary.time_wall_clock` (2.x).
fn stat(r: &Value, key: &str) -> Option<f64> {
    r[key]
        .as_f64()
        .or_else(|| r["summary"]["time_wall_clock"][key].as_f64())
}

pub fn parse_bench(case: &str, hf: &Value, rss: Option<&Value>) -> Result<BenchCase> {
    let rows = hf["results"]
        .as_array()
        .ok_or_else(|| anyhow!("{case}.json has no `results` array"))?;
    let rss_of = |tool: &str| -> Option<u64> {
        rss?["results"]
            .as_array()?
            .iter()
            .find(|r| r["command"].as_str() == Some(tool))?["max_rss_kb"]
            .as_u64()
    };
    let mut tools = Vec::new();
    for r in rows {
        let tool = r["name"]
            .as_str()
            .or_else(|| r["command"].as_str())
            .ok_or_else(|| anyhow!("{case}.json: a result has no command name"))?
            .to_string();
        let median =
            stat(r, "median").ok_or_else(|| anyhow!("{case}.json: `{tool}` has no median"))?;
        // hyperfine 2: `measurements`, `summary.<metric>.count`; 1.x: `times`.
        let runs = r["measurements"]
            .as_array()
            .map(|t| t.len() as u64)
            .or_else(|| r["summary"]["time_wall_clock"]["count"].as_u64())
            .or_else(|| r["times"].as_array().map(|t| t.len() as u64))
            .or_else(|| r["runs"].as_u64());
        tools.push(ToolRow {
            rss_kb: rss_of(&tool),
            mean: stat(r, "mean"),
            stddev: stat(r, "stddev"),
            min: stat(r, "min"),
            max: stat(r, "max"),
            runs,
            median,
            tool,
        });
    }
    Ok(BenchCase {
        case: case.to_string(),
        label: rss.and_then(|r| r["label"].as_str()).map(String::from),
        threads: rss.and_then(|r| r["threads"].as_u64()),
        tools,
    })
}

/// Split CSV text into records (RFC 4180 quoting).
pub fn csv_parse(text: &str) -> Vec<Vec<String>> {
    let mut rows = Vec::new();
    let mut row = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = text.chars().peekable();
    let mut any = false;
    while let Some(c) = chars.next() {
        any = true;
        if quoted {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    field.push('"');
                } else {
                    quoted = false;
                }
            } else {
                field.push(c);
            }
            continue;
        }
        match c {
            '"' => quoted = true,
            ',' => row.push(std::mem::take(&mut field)),
            '\r' => {}
            '\n' => {
                row.push(std::mem::take(&mut field));
                rows.push(std::mem::take(&mut row));
                any = false;
            }
            _ => field.push(c),
        }
    }
    if any {
        row.push(field);
        rows.push(row);
    }
    rows
}

/// One CSV line (no newline), quoting what needs it.
pub fn csv_line<S: AsRef<str>>(fields: &[S]) -> String {
    fields
        .iter()
        .map(|f| {
            let f = f.as_ref();
            if f.contains([',', '"', '\n', '\r']) {
                format!("\"{}\"", f.replace('"', "\"\""))
            } else {
                f.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn parse_concordance(text: &str) -> Result<Vec<ConcRow>> {
    let rows = csv_parse(text);
    let header = rows
        .first()
        .ok_or_else(|| anyhow!("concordance.csv is empty"))?;
    let idx = |name: &str| -> Result<usize> {
        header
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| anyhow!("concordance.csv has no column `{name}`"))
    };
    let cols: Vec<usize> = CONCORDANCE_COLUMNS
        .iter()
        .map(|c| idx(c))
        .collect::<Result<_>>()?;
    let mut out = Vec::new();
    for (n, r) in rows.iter().enumerate().skip(1) {
        ensure!(
            r.len() == header.len(),
            "concordance.csv line {}: {} fields, expected {}",
            n + 1,
            r.len(),
            header.len()
        );
        let get = |i: usize| r[cols[i]].as_str();
        let num = |i: usize| -> Result<u64> {
            get(i).parse().with_context(|| {
                format!(
                    "concordance.csv line {}: `{}` is not a count",
                    n + 1,
                    get(i)
                )
            })
        };
        out.push(ConcRow {
            case: get(0).into(),
            command: get(1).into(),
            reference_tool: get(2).into(),
            dataset: get(3).into(),
            mode: get(4).into(),
            n_gloci: num(5)?,
            n_reference: num(6)?,
            n_matching: num(7)?,
            match_rate: get(8).parse().with_context(|| {
                format!(
                    "concordance.csv line {}: bad match_rate `{}`",
                    n + 1,
                    get(8)
                )
            })?,
            anchors: get(10)
                .split([';', '|', ' '])
                .map(|a| a.trim().trim_start_matches('#'))
                .filter(|a| !a.is_empty())
                .map(String::from)
                .collect(),
        });
    }
    Ok(out)
}

pub fn load(dir: &Path) -> Result<Results> {
    ensure!(
        dir.is_dir(),
        "{} is not a results directory (`cargo xtask paper run` writes paper/results/<version>/)",
        dir.display()
    );
    let environment = read_json(&dir.join("environment.json"))?;
    let snapshot = read_json(&dir.join("repo-snapshot.json"))?;
    let version = environment["gloci_version"]
        .as_str()
        .map(String::from)
        .or_else(|| dir.file_name().map(|n| n.to_string_lossy().into_owned()))
        .ok_or_else(|| anyhow!("no gloci_version in environment.json"))?;
    let mut bench = Vec::new();
    let bench_dir = dir.join("bench");
    if bench_dir.is_dir() {
        let mut names: Vec<String> = std::fs::read_dir(&bench_dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter(|n| n.ends_with(".json") && !n.ends_with(".rss.json"))
            .collect();
        names.sort();
        for n in names {
            let case = n.trim_end_matches(".json");
            let hf = read_json(&bench_dir.join(&n))?;
            let rss_path = bench_dir.join(format!("{case}.rss.json"));
            let rss = if rss_path.exists() {
                Some(read_json(&rss_path)?)
            } else {
                None
            };
            bench.push(parse_bench(case, &hf, rss.as_ref())?);
        }
    }
    let conc_path = dir.join("concordance.csv");
    let concordance = if conc_path.exists() {
        parse_concordance(&std::fs::read_to_string(&conc_path)?)?
    } else {
        Vec::new()
    };
    if bench.is_empty() && concordance.is_empty() {
        bail!(
            "{} has neither benchmark results nor concordance.csv",
            dir.display()
        );
    }
    Ok(Results {
        version,
        bench,
        concordance,
        environment,
        snapshot,
    })
}

/// Seconds since the epoch of `YYYY-MM-DDTHH:MM:SSZ`.
pub fn parse_utc(s: &str) -> Result<u64> {
    let bad = || anyhow!("`{s}` is not a UTC time like 2026-11-12T10:30:15Z");
    let b = s.as_bytes();
    ensure!(b.len() == 20 && b[19] == b'Z', bad());
    let n = |a: usize, z: usize| -> Result<i64> { s[a..z].parse::<i64>().map_err(|_| bad()) };
    let (y, m, d) = (n(0, 4)?, n(5, 7)?, n(8, 10)?);
    let (hh, mm, ss) = (n(11, 13)?, n(14, 16)?, n(17, 19)?);
    ensure!(
        (1..=12).contains(&m) && (1..=31).contains(&d) && hh < 24 && mm < 60 && ss < 61,
        bad()
    );
    // Days from civil (Howard Hinnant).
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    let days = era * 146_097 + doe - 719_468;
    u64::try_from(days * 86_400 + hh * 3600 + mm * 60 + ss).map_err(|_| bad())
}

/// `secs` since the epoch as `YYYY-MM-DDTHH:MM:SSZ`.
pub fn format_utc(secs: u64) -> String {
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

pub fn now_utc() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    format_utc(secs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn utc_round_trip() {
        assert_eq!(parse_utc("2026-11-12T10:30:15Z").unwrap(), 1_794_479_415);
        assert_eq!(format_utc(1_794_479_415), "2026-11-12T10:30:15Z");
        assert_eq!(format_utc(0), "1970-01-01T00:00:00Z");
        assert_eq!(parse_utc("2024-02-29T23:59:59Z").unwrap(), 1_709_251_199);
        assert!(parse_utc("2026-11-12 10:30:15").is_err());
    }

    #[test]
    fn hyperfine2_file_has_runs() {
        let p = crate::diff::workspace_root()
            .join("xtask/tests/fixtures/paper-hyperfine2/bedcov-ctcf-peaks.json");
        let v: Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
        let b = parse_bench("bedcov-ctcf-peaks", &v, None).unwrap();
        assert_eq!(b.tools.len(), 2);
        for t in &b.tools {
            assert_eq!(t.runs, Some(3), "{}", t.tool);
            assert!(t.mean.is_some() && t.stddev.is_some() && t.min.is_some() && t.max.is_some());
        }
    }

    #[test]
    fn csv_quotes() {
        let l = csv_line(&["a,b", "c\"d", "e"]);
        assert_eq!(l, "\"a,b\",\"c\"\"d\",e");
        assert_eq!(
            csv_parse(&format!("{l}\nx,y,z\n")),
            vec![
                vec!["a,b".to_string(), "c\"d".into(), "e".into()],
                vec!["x".into(), "y".into(), "z".into()],
            ]
        );
    }
}
