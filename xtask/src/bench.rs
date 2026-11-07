//! Benchmark runner: generates the input for each case, runs every tool under
//! hyperfine (inside the pinned reference-tool container by default), stores
//! the JSON results and evaluates the case's targets.
use anyhow::{Context, Result, anyhow, bail, ensure};
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;

pub const DEFAULT_IMAGE: &str = "gloci-diff:p0";

/// Median wall time per tool, in seconds.
pub struct Medians(HashMap<String, f64>);

pub fn medians(v: &[(&str, f64)]) -> Medians {
    Medians(v.iter().map(|(k, x)| (k.to_string(), *x)).collect())
}

/// Evaluate `NAME (<|<=) [FACTOR *] NAME`.
pub fn eval(expr: &str, m: &Medians) -> Result<bool> {
    let t: Vec<&str> = expr.split_whitespace().collect();
    let (lhs, op, factor, rhs) = match t.as_slice() {
        [l, op, r] => (*l, *op, 1.0, *r),
        [l, op, f, "*", r] => (
            *l,
            *op,
            f.parse::<f64>()
                .with_context(|| format!("bad factor in target `{expr}`"))?,
            *r,
        ),
        _ => bail!("bad target `{expr}`: expected NAME (<|<=) [FACTOR *] NAME"),
    };
    let get = |n: &str| {
        m.0.get(n)
            .copied()
            .ok_or_else(|| anyhow!("target `{expr}`: no median for `{n}`"))
    };
    let (a, b) = (get(lhs)?, factor * get(rhs)?);
    match op {
        "<" => Ok(a < b),
        "<=" => Ok(a <= b),
        _ => bail!("bad operator `{op}` in target `{expr}` (use < or <=)"),
    }
}

#[derive(Deserialize)]
struct CaseFile {
    case: Vec<Case>,
}

#[derive(Deserialize, Clone)]
pub struct Case {
    pub name: String,
    /// Arguments of `cargo xtask gen`, e.g. `bed --records 10000000 --seed 1`.
    pub setup: String,
    /// `bed` (default) or `bed.gz`: the extension of the input, which is
    /// BGZF-compressed with `gloci view` when it ends in `.gz`.
    #[serde(default = "default_ext")]
    pub input: String,
    /// Tool name -> shell command with `{gloci}`, `{in}` and `{out}`.
    pub commands: BTreeMap<String, String>,
    pub targets: Vec<String>,
}

fn default_ext() -> String {
    "bed".into()
}

pub fn parse_cases(text: &str) -> Result<Vec<Case>> {
    Ok(toml::from_str::<CaseFile>(text)?.case)
}

/// Replace the value after `--records` in a setup string.
pub fn override_records(setup: &str, records: u64) -> String {
    let mut t: Vec<String> = setup.split_whitespace().map(String::from).collect();
    if let Some(i) = t.iter().position(|x| x == "--records")
        && i + 1 < t.len()
    {
        t[i + 1] = records.to_string();
    }
    t.join(" ")
}

pub fn substitute(cmd: &str, gloci: &str, input: &str, output: &str) -> String {
    cmd.replace("{gloci}", gloci)
        .replace("{in}", input)
        .replace("{out}", output)
}

/// Read `results[].{command,median}` from a hyperfine `--export-json` file.
pub fn parse_hyperfine(json: &str) -> Result<Vec<(String, f64)>> {
    let v: serde_json::Value = serde_json::from_str(json)?;
    let rs = v["results"]
        .as_array()
        .ok_or_else(|| anyhow!("hyperfine JSON has no results"))?;
    rs.iter()
        .map(|r| {
            Ok((
                // `name` is the -n label; older hyperfine put it in `command`.
                r["name"]
                    .as_str()
                    .or_else(|| r["command"].as_str())
                    .ok_or_else(|| anyhow!("result without command"))?
                    .to_string(),
                // hyperfine >= 1.20 nests the statistics under `summary`.
                r["median"]
                    .as_f64()
                    .or_else(|| r["summary"]["time_wall_clock"]["median"].as_f64())
                    .ok_or_else(|| anyhow!("result without median"))?,
            ))
        })
        .collect()
}

pub struct Args {
    pub cases: PathBuf,
    pub out: PathBuf,
    pub check_targets: bool,
    /// gloci binary to benchmark; default: `cargo build --release` for the host.
    pub gloci: Option<PathBuf>,
    /// Run hyperfine on the host with host tools instead of in the container.
    pub no_docker: bool,
    pub image: String,
    /// Only run the case with this name.
    pub only: Option<String>,
    /// Override `--records` in every case's setup (reduced local runs).
    pub records: Option<u64>,
    pub data_dir: PathBuf,
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn run_ok(cmd: &mut Command) -> Result<()> {
    let st = cmd.status().with_context(|| format!("running {cmd:?}"))?;
    ensure!(st.success(), "{cmd:?} failed: {st}");
    Ok(())
}

fn id_of(flag: &str) -> Result<String> {
    let o = Command::new("id").arg(flag).output()?;
    Ok(String::from_utf8(o.stdout)?.trim().to_string())
}

/// Generate (or reuse) the input for a case; returns its path.
fn prepare_input(case: &Case, setup: &str, gloci: &Path, dir: &Path) -> Result<PathBuf> {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    setup.hash(&mut h);
    let stem = format!("{}-{:08x}", case.name, h.finish() as u32);
    let plain = dir.join(format!("{stem}.bed"));
    let want = dir.join(format!("{stem}.{}", case.input));
    if want.exists() {
        return Ok(want);
    }
    std::fs::create_dir_all(dir)?;
    if !plain.exists() {
        eprintln!("generating {}", plain.display());
        let tmp = plain.with_extension("bed.partial");
        let args: Vec<String> = setup.split_whitespace().map(String::from).collect();
        let mut f = std::fs::File::create(&tmp)?;
        crate::synth::run_args(&args, &mut f)?;
        std::fs::rename(&tmp, &plain)?;
    }
    if case.input.ends_with(".gz") {
        eprintln!("compressing to {}", want.display());
        // Write to a partial name first so an interrupted run never leaves a
        // truncated cached input.
        let tmp = dir.join(format!("{stem}.partial.{}", case.input));
        run_ok(
            Command::new(gloci)
                .arg("view")
                .arg(&plain)
                .arg("-o")
                .arg(&tmp),
        )?;
        std::fs::rename(&tmp, &want)?;
        std::fs::remove_file(&plain)?;
    }
    Ok(want)
}

pub fn run(a: &Args) -> Result<bool> {
    let root = workspace_root();
    let gloci = match &a.gloci {
        Some(p) => std::fs::canonicalize(p)
            .with_context(|| format!("gloci binary {} not found", p.display()))?,
        None => {
            run_ok(
                Command::new("cargo")
                    .args(["build", "--release", "-p", "gloci"])
                    .current_dir(&root),
            )?;
            root.join("target/release/gloci")
        }
    };
    let cases = parse_cases(&std::fs::read_to_string(&a.cases)?)?;
    let version = env!("CARGO_PKG_VERSION");
    let results = a.out.join(version);
    std::fs::create_dir_all(&results)?;
    let results = std::fs::canonicalize(&results)?;
    std::fs::create_dir_all(&a.data_dir)?;
    let data = std::fs::canonicalize(&a.data_dir)?;
    if !a.no_docker {
        let ok = Command::new("docker")
            .args(["image", "inspect", &a.image])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        ensure!(
            ok,
            "image {0} not found: docker build -t {0} xtask/docker",
            a.image
        );
    }
    let mut all_ok = true;
    let mut ran = 0;
    for case in cases
        .iter()
        .filter(|c| a.only.as_ref().is_none_or(|o| *o == c.name))
    {
        ran += 1;
        let setup = match a.records {
            Some(n) => override_records(&case.setup, n),
            None => case.setup.clone(),
        };
        let input = prepare_input(case, &setup, &gloci, &data)?;
        let out_name = format!("{}.out.{}", case.name, case.input);
        let json = results.join(format!("{}.json", case.name));
        let (gloci_cmd, in_arg, out_arg) = if a.no_docker {
            (
                gloci.to_string_lossy().into_owned(),
                input.to_string_lossy().into_owned(),
                data.join(&out_name).to_string_lossy().into_owned(),
            )
        } else {
            (
                "gloci".to_string(),
                format!("/data/{}", input.file_name().unwrap().to_string_lossy()),
                format!("/data/{out_name}"),
            )
        };
        let mut hf: Vec<String> = [
            "hyperfine",
            "--shell=default",
            "--warmup",
            "1",
            "--runs",
            "5",
            "--export-json",
        ]
        .map(String::from)
        .to_vec();
        hf.push(if a.no_docker {
            json.to_string_lossy().into_owned()
        } else {
            format!("/results/{}.json", case.name)
        });
        for (tool, cmd) in &case.commands {
            hf.push("-n".into());
            hf.push(tool.clone());
            hf.push(substitute(cmd, &gloci_cmd, &in_arg, &out_arg));
        }
        eprintln!("== {} ({} tools)", case.name, case.commands.len());
        let mut c;
        if a.no_docker {
            c = Command::new(&hf[0]);
            c.args(&hf[1..]);
        } else {
            c = Command::new("docker");
            c.args([
                "run",
                "--rm",
                "--user",
                &format!("{}:{}", id_of("-u")?, id_of("-g")?),
            ])
            .args(["-e", "HOME=/tmp"])
            .arg("-v")
            .arg(format!("{}:/data", data.display()))
            .arg("-v")
            .arg(format!("{}:/results", results.display()))
            .arg("-v")
            .arg(format!("{}:/usr/local/bin/gloci:ro", gloci.display()))
            .arg(&a.image)
            .args(&hf);
        }
        let st = c.status().context("running hyperfine")?;
        let _ = std::fs::remove_file(data.join(&out_name));
        ensure!(st.success(), "hyperfine failed for case {}", case.name);
        let med = parse_hyperfine(&std::fs::read_to_string(&json)?)?;
        let pairs: Vec<(&str, f64)> = med.iter().map(|(k, v)| (k.as_str(), *v)).collect();
        for (k, v) in &pairs {
            eprintln!("  {k:<12} median {v:.3} s");
        }
        let m = medians(&pairs);
        for t in &case.targets {
            let pass = eval(t, &m)?;
            eprintln!("  target `{t}`: {}", if pass { "ok" } else { "MISS" });
            all_ok &= pass;
        }
    }
    ensure!(ran > 0, "no case matched");
    Ok(all_ok || !a.check_targets)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn targets_evaluate() {
        let m = medians(&[("gloci", 1.0), ("sort-bed", 1.5), ("gnu-sort", 0.8)]);
        assert!(eval("gloci < sort-bed", &m).unwrap());
        assert!(eval("gloci <= 1.5 * gnu-sort", &m).unwrap());
        assert!(!eval("gloci < gnu-sort", &m).unwrap());
        assert!(eval("gloci < nosuch", &m).is_err());
    }

    #[test]
    fn bad_targets_are_errors() {
        let m = medians(&[("a", 1.0), ("b", 2.0)]);
        for bad in ["a > b", "a <", "a < x * b", "a < 2 b", "a < 2 * b * b", ""] {
            assert!(eval(bad, &m).is_err(), "{bad}");
        }
        assert!(eval("a <= 0.5 * b", &m).unwrap());
        assert!(!eval("a < 0.5 * b", &m).unwrap());
    }

    #[test]
    fn cases_file_parses_and_targets_name_commands() {
        let text =
            std::fs::read_to_string(workspace_root().join("xtask/bench/cases.toml")).unwrap();
        let cases = parse_cases(&text).unwrap();
        let names: Vec<_> = cases.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["sort-plain-10M", "sort-bgzf-10M"]);
        for c in &cases {
            for t in &c.targets {
                // every tool a target mentions has a command
                let m = medians(
                    &c.commands
                        .keys()
                        .map(|k| (k.as_str(), 1.0))
                        .collect::<Vec<_>>(),
                );
                eval(t, &m).unwrap_or_else(|e| panic!("{}: {e}", c.name));
            }
        }
    }

    #[test]
    fn hyperfine_json_and_substitution() {
        let j = r#"{"results":[{"command":"gloci","median":1.25},{"command":"sort-bed","median":2.0}]}"#;
        assert_eq!(
            parse_hyperfine(j).unwrap(),
            vec![("gloci".to_string(), 1.25), ("sort-bed".to_string(), 2.0)]
        );
        let new = r#"{"results":[{"command":"gloci","name":"gloci","summary":{"time_wall_clock":{"median":3.5}}}]}"#;
        assert_eq!(
            parse_hyperfine(new).unwrap(),
            vec![("gloci".to_string(), 3.5)]
        );
        assert_eq!(
            substitute("{gloci} sort {in} > {out}", "g", "i", "o"),
            "g sort i > o"
        );
        assert_eq!(
            override_records("bed --records 10000000 --seed 1", 5),
            "bed --records 5 --seed 1"
        );
    }
}
