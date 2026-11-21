//! `environment.json`: where and with what a paper run happened. A value that
//! cannot be read is `null`, never a guess.
use std::collections::BTreeSet;
use std::path::Path;
use std::process::Command;

use anyhow::{Result, ensure};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

use crate::bench::Case;

/// Stdout of a command, trimmed; None when it cannot run or fails.
fn capture(prog: &str, args: &[&str]) -> Option<String> {
    let out = Command::new(prog).args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
    (!s.is_empty()).then_some(s)
}

fn proc_field(file: &str, keys: &[&str]) -> Option<String> {
    let text = std::fs::read_to_string(file).ok()?;
    for key in keys {
        for l in text.lines() {
            if let Some((k, v)) = l.split_once(':')
                && k.trim() == *key
                && !matches!(v.trim(), "" | "-" | "unknown")
            {
                return Some(v.trim().to_string());
            }
        }
    }
    None
}

pub fn cpu_model() -> Option<String> {
    proc_field("/proc/cpuinfo", &["model name", "Model", "Hardware"])
        .or_else(|| {
            capture("lscpu", &[]).and_then(|t| {
                t.lines().find_map(|l| {
                    let (k, v) = l.split_once(':')?;
                    (k.trim() == "Model name" && !matches!(v.trim(), "" | "-" | "unknown"))
                        .then(|| v.trim().to_string())
                })
            })
        })
        .or_else(|| capture("sysctl", &["-n", "machdep.cpu.brand_string"]))
}

pub fn ram_bytes() -> Option<u64> {
    if let Some(kb) = proc_field("/proc/meminfo", &["MemTotal"]) {
        return kb
            .trim_end_matches("kB")
            .trim()
            .parse::<u64>()
            .ok()
            .map(|k| k * 1024);
    }
    capture("sysctl", &["-n", "hw.memsize"])?.parse().ok()
}

fn host() -> Value {
    json!({
        "cpu_model": cpu_model(),
        "cores": std::thread::available_parallelism().ok().map(|n| n.get()),
        "ram_bytes": ram_bytes(),
        "kernel": capture("uname", &["-r"]),
        "arch": capture("uname", &["-m"]),
    })
}

fn docker_info(image: &str) -> Value {
    let id = capture(
        "docker",
        &["image", "inspect", "--format", "{{.Id}}", image],
    );
    let digests = capture(
        "docker",
        &[
            "image",
            "inspect",
            "--format",
            "{{json .RepoDigests}}",
            image,
        ],
    )
    .and_then(|j| serde_json::from_str::<Vec<String>>(&j).ok())
    .filter(|v| !v.is_empty());
    json!({"image": image, "image_id": id, "repo_digests": digests})
}

/// The binary a case's tool name runs, and its version arguments.
fn version_command(tool: &str) -> (&str, &'static [&'static str]) {
    match tool {
        "gnu-sort" => ("sort", &["--version"]),
        "bedtk" => ("bedtk", &["version"]),
        _ => (tool, &["--version"]),
    }
}

/// A tool's version from its version command's output (either stream):
/// - a `version: X` line, shown after the first line (bedops and sort-bed
///   print their name, then a citation with digits, then `version:`);
/// - else the first line that has a digit;
/// - else the first non-empty line.
fn first_line(out: &std::process::Output) -> Option<String> {
    let both = [&out.stdout[..], &out.stderr[..]].concat();
    let text = String::from_utf8_lossy(&both);
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect();
    let first = *lines.first()?;
    let labelled = lines.iter().find_map(|l| {
        let (key, value) = l.split_once(':')?;
        key.trim()
            .eq_ignore_ascii_case("version")
            .then(|| format!("{first} {}", value.trim()))
    });
    labelled
        .or_else(|| {
            lines
                .iter()
                .find(|l| l.chars().any(|c| c.is_ascii_digit()))
                .map(|l| l.to_string())
        })
        .or_else(|| Some(first.to_string()))
}

/// `--user UID:GID` of the current user, so a version query does not run
/// the container as root; empty when `id` is not available.
fn docker_user() -> Vec<String> {
    match (capture("id", &["-u"]), capture("id", &["-g"])) {
        (Some(u), Some(g)) => vec!["--user".into(), format!("{u}:{g}")],
        _ => Vec::new(),
    }
}

fn tool_version(image: Option<&str>, tool: &str) -> Value {
    let (bin, args) = version_command(tool);
    let out = match image {
        Some(img) => Command::new("docker")
            .args(["run", "--rm", "--network", "none"])
            .args(docker_user().as_slice())
            .args([img, bin])
            .args(args)
            .output(),
        None => Command::new(bin).args(args).output(),
    };
    match out {
        Ok(o) if o.status.success() => first_line(&o).map_or(Value::Null, Value::from),
        _ => Value::Null,
    }
}

/// Names a case gives to another command of a tool already listed (`gloci`
/// is recorded apart): `vcftools-3runs` runs vcftools, `bcftools-stats`
/// bcftools.
fn is_alias(tool: &str) -> bool {
    tool == "gloci"
        || tool == "bcftools-stats"
        || tool.rsplit_once('-').is_some_and(|(_, s)| {
            s.strip_suffix("runs")
                .is_some_and(|n| n.parse::<u32>().is_ok())
        })
}

/// The tools the cases run (by the name the cases give them, so `gnu-sort`
/// stays `gnu-sort`), plus the ones the datasets need, without the aliases.
fn tool_names(cases: &[Case]) -> BTreeSet<String> {
    let mut s: BTreeSet<String> = ["hyperfine", "bgzip", "tabix"].map(String::from).into();
    for c in cases {
        for t in c.commands.keys().filter(|t| !is_alias(t)) {
            s.insert(t.clone());
        }
    }
    s
}

pub fn sha256_file(path: &Path) -> Result<String> {
    use std::io::Read;
    let mut h = Sha256::new();
    let mut f = std::fs::File::open(path)?;
    let mut buf = vec![0u8; 1 << 16];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(h.finalize().iter().map(|b| format!("{b:02x}")).collect())
}

/// `{file_name, sha256}` of the wheel; the file must exist.
pub fn wheel_info(wheel: &Path) -> Result<Value> {
    ensure!(
        wheel.is_file(),
        "--wheel {}: no such file (build one with `maturin build --release`)",
        wheel.display()
    );
    Ok(json!({
        "file_name": wheel.file_name().map(|n| n.to_string_lossy().into_owned()),
        "sha256": sha256_file(wheel)?,
    }))
}

pub struct EnvArgs<'a> {
    pub version: &'a str,
    pub commit: Option<String>,
    pub started_at: &'a str,
    pub finished_at: &'a str,
    /// Bench image; None for a host run.
    pub image: Option<&'a str>,
    pub gloci: Option<&'a Path>,
    pub wheel: Option<&'a Path>,
    pub cases: &'a [Case],
    /// Tools run outside the benchmark cases (the concordance run's
    /// reference tools), recorded too.
    pub extra_tools: &'a [String],
}

pub fn environment(a: &EnvArgs) -> Result<Value> {
    let mut tools = serde_json::Map::new();
    let gloci_version = a.gloci.and_then(|p| {
        let o = Command::new(p).arg("--version").output().ok()?;
        o.status.success().then(|| first_line(&o)).flatten()
    });
    tools.insert(
        "gloci".into(),
        gloci_version.map_or(Value::Null, Value::from),
    );
    let mut names = tool_names(a.cases);
    names.extend(a.extra_tools.iter().cloned());
    for t in names {
        tools.insert(t.clone(), tool_version(a.image, &t));
    }
    Ok(json!({
        "gloci_version": a.version,
        "gloci_commit": a.commit,
        "started_at": a.started_at,
        "finished_at": a.finished_at,
        "host": host(),
        "docker": a.image.map(docker_info),
        "wheel": a.wheel.map(wheel_info).transpose()?,
        "tools": tools,
    }))
}

pub fn git_commit(root: &Path) -> Option<String> {
    capture("git", &["-C", &root.to_string_lossy(), "rev-parse", "HEAD"])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::ExitStatusExt;

    fn out(s: &str) -> std::process::Output {
        std::process::Output {
            status: std::process::ExitStatus::from_raw(0),
            stdout: s.as_bytes().to_vec(),
            stderr: Vec::new(),
        }
    }

    #[test]
    fn version_line_has_a_digit() {
        assert_eq!(
            first_line(&out("bedops\n  version:  2.4.41 (typical)\n")).as_deref(),
            Some("bedops 2.4.41 (typical)")
        );
        // bedops and sort-bed print their citation (with digits) before the version.
        let bedops = "bedops\n  citation: http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract\n            https://doi.org/10.1093/bioinformatics/bts277\n  version:  2.4.41 (typical)\n  authors:  Shane Neph & Scott Kuehn\n";
        assert_eq!(
            first_line(&out(bedops)).as_deref(),
            Some("bedops 2.4.41 (typical)")
        );
        assert_eq!(
            first_line(&out("bedtools v2.31.1\n")).as_deref(),
            Some("bedtools v2.31.1")
        );
        assert_eq!(first_line(&out("tool\nother\n")).as_deref(), Some("tool"));
        assert_eq!(first_line(&out("\n")), None);
    }
}
