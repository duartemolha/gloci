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

#[derive(Deserialize, Clone, Debug)]
pub struct Dataset {
    pub name: String,
    pub url: String,
    /// The checksum to verify: the provider's own when it publishes one
    /// (EBI publishes md5 manifests), else a sha256 computed on download.
    /// Exactly one of the two is given.
    #[serde(default)]
    pub sha256: Option<String>,
    #[serde(default)]
    pub md5: Option<String>,
    #[serde(default)]
    pub unpack: bool,
}

impl Dataset {
    /// The checksum algorithm and expected hex digest.
    pub fn checksum(&self) -> Result<(&'static str, &str)> {
        match (&self.sha256, &self.md5) {
            (Some(h), None) => Ok(("sha256", h)),
            (None, Some(h)) => Ok(("md5", h)),
            _ => bail!(
                "dataset `{}` needs exactly one of sha256 and md5",
                self.name
            ),
        }
    }
}

#[derive(Deserialize)]
struct DatasetsFile {
    dataset: Vec<Dataset>,
}

pub fn parse_datasets(text: &str) -> Result<Vec<Dataset>> {
    Ok(toml::from_str::<DatasetsFile>(text)?.dataset)
}

#[derive(Deserialize)]
struct CaseFile {
    case: Vec<Case>,
}

#[derive(Deserialize, Clone)]
pub struct Case {
    pub name: String,
    /// Arguments of `cargo xtask gen`, e.g. `bed --records 10000000 --seed 1`.
    #[serde(default)]
    pub setup: Option<String>,
    /// Dataset from `datasets.toml`, e.g. `1kg-30x-chr20`.
    #[serde(default)]
    pub dataset: Option<String>,
    /// `bed` (default), `bed.gz`, `vcf.gz`, etc.
    #[serde(default = "default_ext")]
    pub input: String,
    /// Tool name -> shell command with `{gloci}`, `{in}`, `{in_a}`, `{in_b}` and `{out}`.
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
    substitute_pair(cmd, gloci, input, input, output)
}

pub fn substitute_pair(cmd: &str, gloci: &str, in_a: &str, in_b: &str, output: &str) -> String {
    cmd.replace("{gloci}", gloci)
        .replace("{in_a}", in_a)
        .replace("{in_b}", in_b)
        .replace("{a}", in_a)
        .replace("{b}", in_b)
        .replace("{in}", in_a)
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

pub fn bench_cache_dir(data_dir: &Path) -> PathBuf {
    match std::env::var("GLOCI_BENCH_CACHE") {
        Ok(s) if !s.is_empty() => PathBuf::from(s),
        _ => data_dir.join("cache"),
    }
}

/// The hex digest of `path` with `algorithm` (`sha256` or `md5`).
fn digest_file(path: &Path, algorithm: &str) -> Result<String> {
    use sha2::Digest;
    fn hash<D: Digest + std::io::Write>(path: &Path, mut d: D) -> Result<String> {
        let mut f = std::fs::File::open(path)?;
        std::io::copy(&mut f, &mut d)?;
        Ok(d.finalize()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>())
    }
    match algorithm {
        "sha256" => hash(path, sha2::Sha256::new()),
        "md5" => hash(path, md5::Md5::new()),
        _ => bail!("unknown checksum algorithm {algorithm}"),
    }
}

pub fn ensure_dataset(ds: &Dataset, cache_dir: &Path) -> Result<PathBuf> {
    std::fs::create_dir_all(cache_dir)?;
    let filename = ds.url.split('/').next_back().unwrap_or(&ds.name);
    let target_path = cache_dir.join(filename);
    if !target_path.exists() {
        eprintln!("fetching dataset {} from {}", ds.name, ds.url);
        let tmp = cache_dir.join(format!("{filename}.tmp"));
        let st = Command::new("curl")
            .args(["-fsSL", "-o"])
            .arg(&tmp)
            .arg(&ds.url)
            .status()
            .with_context(|| format!("failed to run curl to download {}", ds.url))?;
        ensure!(st.success(), "curl failed to download {}", ds.url);
        std::fs::rename(&tmp, &target_path)?;
    }
    let (algorithm, want) = ds.checksum()?;
    let hash = digest_file(&target_path, algorithm)?;
    ensure!(
        hash.eq_ignore_ascii_case(want),
        "{algorithm} mismatch for {}: expected {want}, got {hash}; delete {} and fetch again",
        ds.name,
        target_path.display()
    );
    if ds.unpack {
        let unpacked_dir = cache_dir.join(&ds.name);
        if !unpacked_dir.exists() {
            eprintln!("unpacking {}", target_path.display());
            std::fs::create_dir_all(&unpacked_dir)?;
            let st = Command::new("tar")
                .args(["-xf"])
                .arg(&target_path)
                .args(["-C"])
                .arg(&unpacked_dir)
                .status()
                .with_context(|| format!("failed to unpack {}", target_path.display()))?;
            ensure!(
                st.success(),
                "tar failed to unpack {}",
                target_path.display()
            );
        }
        Ok(unpacked_dir)
    } else {
        Ok(target_path)
    }
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

struct PreparedInputs {
    in_a: PathBuf,
    in_b: Option<PathBuf>,
}

fn link_or_copy(src: &Path, dst: &Path) -> Result<()> {
    if dst.exists() {
        return Ok(());
    }
    if std::fs::hard_link(src, dst).is_ok() {
        return Ok(());
    }
    std::fs::copy(src, dst)?;
    Ok(())
}

/// Generate (or reuse) the input for a case; returns its paths.
fn prepare_input(
    case: &Case,
    setup: Option<&str>,
    gloci: &Path,
    dir: &Path,
) -> Result<PreparedInputs> {
    std::fs::create_dir_all(dir)?;

    if let Some(ref ds_name) = case.dataset {
        let text = std::fs::read_to_string(workspace_root().join("xtask/bench/datasets.toml"))
            .context("reading datasets.toml")?;
        let datasets = parse_datasets(&text)?;
        let ds = datasets
            .iter()
            .find(|d| d.name == *ds_name)
            .ok_or_else(|| anyhow!("unknown dataset `{ds_name}`"))?;
        let cache_dir = bench_cache_dir(dir);
        let path = ensure_dataset(ds, &cache_dir)?;
        if ds.unpack {
            let in_a = path.join("a.bed");
            let in_b = path.join("b.bed");
            ensure!(
                in_a.exists() && in_b.exists(),
                "dataset `{ds_name}` unpacked at {} is missing a.bed or b.bed",
                path.display()
            );
            let target_a = dir.join(format!("{}-a.bed", case.name));
            let target_b = dir.join(format!("{}-b.bed", case.name));
            link_or_copy(&in_a, &target_a)?;
            link_or_copy(&in_b, &target_b)?;
            return Ok(PreparedInputs {
                in_a: target_a,
                in_b: Some(target_b),
            });
        } else {
            let file_name = path
                .file_name()
                .ok_or_else(|| anyhow!("invalid dataset path {}", path.display()))?;
            let target = dir.join(file_name);
            link_or_copy(&path, &target)?;

            if target.to_string_lossy().ends_with(".gz") {
                let csi = PathBuf::from(format!("{}.csi", target.display()));
                let tbi = PathBuf::from(format!("{}.tbi", target.display()));
                if !csi.exists() && !tbi.exists() {
                    let cached_csi = PathBuf::from(format!("{}.csi", path.display()));
                    let cached_tbi = PathBuf::from(format!("{}.tbi", path.display()));
                    if cached_csi.exists() {
                        let _ = link_or_copy(&cached_csi, &csi);
                    } else if cached_tbi.exists() {
                        let _ = link_or_copy(&cached_tbi, &tbi);
                    } else {
                        run_ok(Command::new(gloci).arg("index").arg(&target))?;
                    }
                }
            }
            return Ok(PreparedInputs {
                in_a: target,
                in_b: None,
            });
        }
    }

    let setup_str =
        setup.ok_or_else(|| anyhow!("case {} has neither setup nor dataset", case.name))?;
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    setup_str.hash(&mut h);
    let stem = format!("{}-{:08x}", case.name, h.finish() as u32);

    let is_pair = setup_str.trim_start().starts_with("pair");
    let is_vcf = setup_str.trim_start().starts_with("vcf");
    let ext = if is_vcf { "vcf" } else { "bed" };

    if is_pair {
        let plain_a = dir.join(format!("{stem}_a.bed"));
        let plain_b = dir.join(format!("{stem}_b.bed"));
        let want_a = dir.join(format!("{stem}_a.{}", case.input));
        let want_b = dir.join(format!("{stem}_b.{}", case.input));
        if !want_a.exists() || !want_b.exists() {
            if !plain_a.exists() || !plain_b.exists() {
                eprintln!("generating {} and {}", plain_a.display(), plain_b.display());
                let tmp_a = plain_a.with_extension("bed.partial");
                let tmp_b = plain_b.with_extension("bed.partial");
                let mut args: Vec<String> =
                    setup_str.split_whitespace().map(String::from).collect();
                args.push("--out-a".into());
                args.push(tmp_a.to_string_lossy().into_owned());
                args.push("--out-b".into());
                args.push(tmp_b.to_string_lossy().into_owned());
                let mut dummy = Vec::new();
                crate::synth::run_args(&args, &mut dummy)?;
                std::fs::rename(&tmp_a, &plain_a)?;
                std::fs::rename(&tmp_b, &plain_b)?;
            }
            if case.input.ends_with(".gz") {
                eprintln!(
                    "compressing and indexing {} and {}",
                    want_a.display(),
                    want_b.display()
                );
                let tmp_a = dir.join(format!("{stem}_a.partial.{}", case.input));
                let tmp_b = dir.join(format!("{stem}_b.partial.{}", case.input));
                run_ok(
                    Command::new(gloci)
                        .arg("view")
                        .arg(&plain_a)
                        .arg("-o")
                        .arg(&tmp_a),
                )?;
                run_ok(
                    Command::new(gloci)
                        .arg("view")
                        .arg(&plain_b)
                        .arg("-o")
                        .arg(&tmp_b),
                )?;
                std::fs::rename(&tmp_a, &want_a)?;
                std::fs::rename(&tmp_b, &want_b)?;
                run_ok(Command::new(gloci).arg("index").arg(&want_a))?;
                run_ok(Command::new(gloci).arg("index").arg(&want_b))?;
                let _ = std::fs::remove_file(&plain_a);
                let _ = std::fs::remove_file(&plain_b);
            }
        }
        Ok(PreparedInputs {
            in_a: want_a,
            in_b: Some(want_b),
        })
    } else {
        let plain = dir.join(format!("{stem}.{ext}"));
        let want = dir.join(format!("{stem}.{}", case.input));
        if !want.exists() {
            if !plain.exists() {
                eprintln!("generating {}", plain.display());
                let tmp = plain.with_extension(format!("{ext}.partial"));
                let args: Vec<String> = setup_str.split_whitespace().map(String::from).collect();
                let mut f = std::fs::File::create(&tmp)?;
                crate::synth::run_args(&args, &mut f)?;
                std::fs::rename(&tmp, &plain)?;
            }
            if case.input.ends_with(".gz") {
                eprintln!("compressing and indexing {}", want.display());
                let tmp = dir.join(format!("{stem}.partial.{}", case.input));
                run_ok(
                    Command::new(gloci)
                        .arg("view")
                        .arg(&plain)
                        .arg("-o")
                        .arg(&tmp),
                )?;
                std::fs::rename(&tmp, &want)?;
                run_ok(Command::new(gloci).arg("index").arg(&want))?;
                let _ = std::fs::remove_file(&plain);
            }
        }
        Ok(PreparedInputs {
            in_a: want,
            in_b: None,
        })
    }
}

/// Remove a case's output and the side files commands write beside it
/// (`OUT.frq`, `OUT.genome`, ...).
fn remove_outputs(data: &Path, out_name: &str) {
    let _ = std::fs::remove_file(data.join(out_name));
    let side = format!("{out_name}.");
    if let Ok(rd) = std::fs::read_dir(data) {
        for e in rd.flatten() {
            if e.file_name().to_string_lossy().starts_with(&side) {
                let _ = std::fs::remove_file(e.path());
            }
        }
    }
}

fn tool_binary(tool: &str) -> &str {
    match tool {
        "gnu-sort" => "sort",
        "grit-genomics" => "grit",
        "vcfkit" => "vk",
        "vcftools-5runs" => "vcftools",
        "bcftools-stats" => "bcftools",
        _ => tool,
    }
}

fn check_tool_available(tool: &str, use_docker: bool, image: &str) -> bool {
    let bin = tool_binary(tool);
    if use_docker {
        Command::new("docker")
            .args(["run", "--rm", image, "which", bin])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    } else {
        Command::new("which")
            .arg(bin)
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }
}

fn run_native_bench(commands: &[(&str, String)], json_path: &Path) -> Result<Vec<(String, f64)>> {
    let mut results = Vec::new();
    for (tool, cmd) in commands {
        // 1 warmup run
        let _ = Command::new("sh").arg("-c").arg(cmd).status();
        // 3 timed runs
        let mut times = Vec::new();
        for _ in 0..3 {
            let t0 = std::time::Instant::now();
            let st = Command::new("sh")
                .arg("-c")
                .arg(cmd)
                .status()
                .with_context(|| format!("failed to execute command: {cmd}"))?;
            ensure!(st.success(), "command failed with status {st}: {cmd}");
            times.push(t0.elapsed().as_secs_f64());
        }
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = times[times.len() / 2];
        results.push((tool.to_string(), median));
    }
    let j = serde_json::json!({
        "results": results.iter().map(|(tool, med)| {
            serde_json::json!({
                "command": tool,
                "median": med,
            })
        }).collect::<Vec<_>>()
    });
    std::fs::write(json_path, serde_json::to_string_pretty(&j)?)?;
    Ok(results)
}

/// Cargo's target directory: `$CARGO_TARGET_DIR` (relative to the
/// workspace root when relative), else `target` in the workspace root.
pub fn target_dir() -> PathBuf {
    let root = workspace_root();
    match std::env::var_os("CARGO_TARGET_DIR") {
        Some(d) if !d.is_empty() => root.join(d),
        _ => root.join("target"),
    }
}

/// The gloci binary of `profile` (`release` or `debug`), built first with
/// cargo when missing.
fn built_gloci(profile: &str) -> Result<PathBuf> {
    let bin = target_dir().join(profile).join("gloci");
    if !bin.exists() {
        let mut cmd = Command::new("cargo");
        cmd.args(["build", "-p", "gloci"])
            .current_dir(workspace_root());
        if profile == "release" {
            cmd.arg("--release");
        }
        run_ok(&mut cmd)?;
    }
    ensure!(bin.exists(), "cargo built no gloci at {}", bin.display());
    Ok(bin)
}

pub fn resolve_gloci(explicit: Option<&Path>) -> Result<PathBuf> {
    match explicit {
        Some(p) => std::fs::canonicalize(p)
            .with_context(|| format!("gloci binary {} not found", p.display())),
        None => built_gloci("release"),
    }
}

pub fn run(a: &Args) -> Result<bool> {
    let gloci = resolve_gloci(a.gloci.as_deref())?;
    let cases = parse_cases(&std::fs::read_to_string(&a.cases)?)?;
    let version = env!("CARGO_PKG_VERSION");
    let results = a.out.join(version);
    std::fs::create_dir_all(&results)?;
    let results = std::fs::canonicalize(&results)?;
    std::fs::create_dir_all(&a.data_dir)?;
    let data = std::fs::canonicalize(&a.data_dir)?;

    let docker_available = if !a.no_docker {
        Command::new("docker")
            .args(["image", "inspect", &a.image])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    } else {
        false
    };
    let use_docker = !a.no_docker && docker_available;
    if !a.no_docker && !docker_available && a.check_targets {
        bail!(
            "image {0} not found: docker build -t {0} -f xtask/docker/Dockerfile.bench xtask/docker",
            a.image
        );
    } else if !a.no_docker && !docker_available {
        eprintln!(
            "Docker image {} not available; running with host tools",
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
        let setup = case.setup.as_ref().map(|s| match a.records {
            Some(n) => override_records(s, n),
            None => s.clone(),
        });
        let inputs = prepare_input(case, setup.as_deref(), &gloci, &data)?;
        let out_name = format!("{}.out.{}", case.name, case.input);
        let json = results.join(format!("{}.json", case.name));

        let (gloci_cmd, in_a_arg, in_b_arg, out_arg) = if use_docker {
            let in_a_file = inputs.in_a.file_name().unwrap().to_string_lossy();
            let in_b_file = inputs
                .in_b
                .as_ref()
                .map(|p| p.file_name().unwrap().to_string_lossy().into_owned());
            (
                "gloci".to_string(),
                format!("/data/{in_a_file}"),
                in_b_file.map(|b| format!("/data/{b}")).unwrap_or_default(),
                format!("/data/{out_name}"),
            )
        } else {
            (
                gloci.to_string_lossy().into_owned(),
                inputs.in_a.to_string_lossy().into_owned(),
                inputs
                    .in_b
                    .as_ref()
                    .map(|p| p.to_string_lossy().into_owned())
                    .unwrap_or_default(),
                data.join(&out_name).to_string_lossy().into_owned(),
            )
        };

        // Filter commands based on tool availability
        let mut active_commands: Vec<(String, String)> = Vec::new();
        for (tool, cmd) in &case.commands {
            if tool == "gloci" {
                active_commands.push((tool.clone(), cmd.clone()));
                continue;
            }
            let is_baseline = case
                .targets
                .iter()
                .any(|t| t.split_whitespace().any(|word| word == tool));
            if check_tool_available(tool, use_docker, &a.image) {
                active_commands.push((tool.clone(), cmd.clone()));
            } else if !is_baseline {
                eprintln!("  skipping yardstick {tool}: executable not found");
            } else if a.check_targets {
                bail!(
                    "binding baseline tool `{tool}` required by target in case `{}` is not available",
                    case.name
                );
            } else {
                eprintln!(
                    "  warning: baseline tool `{tool}` not found (skipped for local run without --check-targets)"
                );
            }
        }
        ensure!(
            !active_commands.is_empty(),
            "no commands available for case {}",
            case.name
        );

        let has_hyperfine = if use_docker {
            true
        } else {
            Command::new("which")
                .arg("hyperfine")
                .output()
                .map(|o| o.status.success())
                .unwrap_or(false)
        };

        if has_hyperfine {
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
            hf.push(if use_docker {
                format!("/results/{}.json", case.name)
            } else {
                json.to_string_lossy().into_owned()
            });
            for (tool, cmd) in &active_commands {
                hf.push("-n".into());
                hf.push(tool.clone());
                hf.push(substitute_pair(
                    cmd, &gloci_cmd, &in_a_arg, &in_b_arg, &out_arg,
                ));
            }
            eprintln!("== {} ({} tools)", case.name, active_commands.len());
            let mut c;
            if use_docker {
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
            } else {
                c = Command::new(&hf[0]);
                c.args(&hf[1..]);
            }
            let st = c.status().context("running hyperfine")?;
            remove_outputs(&data, &out_name);
            ensure!(st.success(), "hyperfine failed for case {}", case.name);
        } else {
            eprintln!(
                "== {} ({} tools, native timing runner)",
                case.name,
                active_commands.len()
            );
            let substituted: Vec<(&str, String)> = active_commands
                .iter()
                .map(|(tool, cmd)| {
                    (
                        tool.as_str(),
                        substitute_pair(cmd, &gloci_cmd, &in_a_arg, &in_b_arg, &out_arg),
                    )
                })
                .collect();
            run_native_bench(&substituted, &json)?;
            remove_outputs(&data, &out_name);
        }

        let med = parse_hyperfine(&std::fs::read_to_string(&json)?)?;
        let pairs: Vec<(&str, f64)> = med.iter().map(|(k, v)| (k.as_str(), *v)).collect();
        for (k, v) in &pairs {
            eprintln!("  {k:<12} median {v:.3} s");
        }
        let m = medians(&pairs);
        for t in &case.targets {
            match eval(t, &m) {
                Ok(pass) => {
                    eprintln!("  target `{t}`: {}", if pass { "ok" } else { "MISS" });
                    all_ok &= pass;
                }
                Err(e) => {
                    if a.check_targets {
                        return Err(e);
                    } else {
                        eprintln!("  target `{t}`: skipped (missing baseline in local run)");
                    }
                }
            }
        }
    }
    ensure!(ran > 0, "no case matched");
    Ok(all_ok || !a.check_targets)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digests() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("abc");
        std::fs::write(&p, b"abc").unwrap();
        assert_eq!(
            digest_file(&p, "md5").unwrap(),
            "900150983cd24fb0d6963f7d28e17f72"
        );
        assert_eq!(
            digest_file(&p, "sha256").unwrap(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    /// Every `gloci` command line in cases.toml parses with gloci's own clap
    /// parser, so a misspelled or removed flag fails `cargo test` instead of
    /// a nightly benchmark run.
    /// The words of a case's gloci command: from `{gloci}` to the first
    /// unquoted pipe, redirection or `&&`, with single quotes removed.
    fn gloci_words(cmd: &str) -> Vec<String> {
        let rest = &cmd[cmd.find("{gloci}").expect("a {gloci} word")..];
        let mut words = Vec::new();
        let (mut cur, mut have, mut quoted) = (String::new(), false, false);
        for ch in rest.chars() {
            match ch {
                '\'' => {
                    quoted = !quoted;
                    have = true;
                }
                '|' | '>' | '&' if !quoted => break,
                c if c.is_whitespace() && !quoted => {
                    if have {
                        words.push(std::mem::take(&mut cur));
                        have = false;
                    }
                }
                c => {
                    cur.push(c);
                    have = true;
                }
            }
        }
        if have {
            words.push(cur);
        }
        words
    }

    #[test]
    fn gloci_words_respect_quotes_and_prefixes() {
        assert_eq!(
            gloci_words(
                "printf x > {out}.g && {gloci} vcf filter --include 'Q>30 && D>1' {in} -o {out}"
            ),
            [
                "{gloci}",
                "vcf",
                "filter",
                "--include",
                "Q>30 && D>1",
                "{in}",
                "-o",
                "{out}"
            ]
        );
        assert_eq!(
            gloci_words("{gloci} sort {in} | wc"),
            ["{gloci}", "sort", "{in}"]
        );
    }

    #[test]
    fn every_gloci_bench_command_parses() {
        use clap::Parser;
        let text =
            std::fs::read_to_string(workspace_root().join("xtask/bench/cases.toml")).unwrap();
        let cases = parse_cases(&text).unwrap();
        let mut checked = 0;
        for c in &cases {
            let Some(cmd) = c.commands.get("gloci") else {
                continue;
            };
            let line = gloci_words(cmd).join(" ");
            let args: Vec<String> = gloci_words(cmd)
                .iter()
                .map(|w| {
                    w.replace("{gloci}", "gloci")
                        .replace("{in_a}", "a.bed")
                        .replace("{in_b}", "b.bed")
                        .replace("{in}", "in.bed")
                        .replace("{out}", "out.bed")
                })
                .collect();
            if let Err(e) = gloci::cli::Cli::try_parse_from(&args) {
                panic!("case `{}`: `{}` does not parse:\n{e}", c.name, line.trim());
            }
            checked += 1;
        }
        assert!(checked >= 25, "only {checked} gloci command lines found");
    }

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
        // Preserves P0 sort cases
        assert!(names.contains(&"sort-plain-10M"));
        assert!(names.contains(&"sort-bgzf-10M"));

        // Requires P1 cases
        let expected = [
            "intersect-plain-10M",
            "merge-plain-10M",
            "subtract-plain-10M",
            "closest-plain-10M",
            "map-plain-10M",
            "coverage-plain-10M",
            "intersect-bgzf-10M",
            "merge-bgzf-10M",
            "subtract-bgzf-10M",
            "closest-bgzf-10M",
            "map-bgzf-10M",
            "coverage-bgzf-10M",
            "vcf-filter-1kg-chr20",
            "vcf-filter-synthetic-2500",
            // P2 cases
            "partition-plain-10M",
            "genomecov-plain-10M",
            "cluster-plain-10M",
            "groupby-plain-10M",
            "shift-plain-10M",
            "vcf-stats-1kg-chr20",
            "vcf-stats-synthetic-2500",
            "vcf-filter-include-synthetic-2500",
            "convert-bcf-synthetic-100",
        ];
        for exp in expected {
            assert!(names.contains(&exp), "missing expected case `{exp}`");
        }

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
    fn datasets_file_parses() {
        let text =
            std::fs::read_to_string(workspace_root().join("xtask/bench/datasets.toml")).unwrap();
        let datasets = parse_datasets(&text).unwrap();
        assert!(datasets.iter().any(|d| d.name == "1kg-30x-chr20"));
        for d in &datasets {
            assert!(!d.name.is_empty());
            assert!(d.url.starts_with("https://") || d.url.starts_with("http://"));
            let (algorithm, hex) = d.checksum().unwrap();
            let len = if algorithm == "md5" { 32 } else { 64 };
            assert_eq!(hex.len(), len, "bad {algorithm} length for {}", d.name);
            assert!(
                hex.chars().all(|c| c.is_ascii_hexdigit()),
                "bad {algorithm} hex for {}",
                d.name
            );
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

    #[test]
    fn prepare_input_indexes_bgzf() {
        // A debug build (or $GLOCI_BIN): fast, and under any CARGO_TARGET_DIR.
        let gloci = match std::env::var_os("GLOCI_BIN") {
            Some(b) => PathBuf::from(b),
            None => built_gloci("debug").unwrap(),
        };
        let td = tempfile::tempdir().unwrap();
        let c = Case {
            name: "test-bgzf".into(),
            setup: Some("bed --records 10 --chroms 1 --seed 1 --sorted".into()),
            dataset: None,
            input: "bed.gz".into(),
            commands: BTreeMap::new(),
            targets: Vec::new(),
        };
        let inputs = prepare_input(&c, c.setup.as_deref(), &gloci, td.path()).unwrap();
        assert!(inputs.in_a.exists());
        let tbi = PathBuf::from(format!("{}.tbi", inputs.in_a.display()));
        let csi = PathBuf::from(format!("{}.csi", inputs.in_a.display()));
        assert!(
            tbi.exists() || csi.exists(),
            "index file must exist for prepared bgzf input at {}",
            inputs.in_a.display()
        );
    }

    /// A case's `vcf filter` command must keep records on its data, or it
    /// times nothing (the 1kg chr20 data has QUAL `.` throughout, so a
    /// `--min-qual` there dropped every record for every tool). Dry-runs
    /// gloci's command on the first 10,000 data lines of the dataset's
    /// committed head fixture (`xtask/bench/fixtures/DATASET.head.vcf`, which
    /// mimics its fields), or of a small generated VCF for setup cases.
    #[test]
    fn vcf_filter_cases_keep_records() {
        use std::io::Write;
        let root = workspace_root();
        let cases =
            parse_cases(&std::fs::read_to_string(root.join("xtask/bench/cases.toml")).unwrap())
                .unwrap();
        let gloci = match std::env::var_os("GLOCI_BIN") {
            Some(b) => PathBuf::from(b),
            None => built_gloci("debug").unwrap(),
        };
        let td = tempfile::tempdir().unwrap();
        let input = td.path().join("in.vcf");
        let mut checked = 0;
        for c in &cases {
            let Some(cmd) = c.commands.get("gloci") else {
                continue;
            };
            if !(cmd.contains("vcf filter")
                || cmd.contains("vcf stats")
                || (cmd.contains(" view ") && c.input.starts_with("vcf")))
            {
                continue;
            }
            let mut w = std::io::BufWriter::new(std::fs::File::create(&input).unwrap());
            if let Some(d) = &c.dataset {
                let fixture = root.join(format!("xtask/bench/fixtures/{d}.head.vcf"));
                let text = std::fs::read_to_string(&fixture).unwrap_or_else(|e| {
                    panic!(
                        "case `{}`: no head fixture {}: {e}",
                        c.name,
                        fixture.display()
                    )
                });
                let mut data = 0;
                for line in text.lines() {
                    if !line.starts_with('#') {
                        data += 1;
                        if data > 10_000 {
                            break;
                        }
                    }
                    writeln!(w, "{line}").unwrap();
                }
            } else {
                let setup = c.setup.as_deref().expect("a setup or a dataset");
                let t: Vec<&str> = setup.split_whitespace().collect();
                assert_eq!(
                    t.first(),
                    Some(&"vcf"),
                    "case `{}`: unsupported setup",
                    c.name
                );
                let get = |k: &str, d: u64| {
                    t.iter()
                        .position(|x| *x == k)
                        .and_then(|i| t.get(i + 1))
                        .map_or(d, |v| v.parse().unwrap())
                };
                crate::synth::vcf(
                    &crate::synth::VcfOpts {
                        records: get("--records", 1000).min(10_000),
                        samples: get("--samples", 10).min(100) as u32,
                        seed: get("--seed", 1),
                        chroms: get("--chroms", 1) as u32,
                    },
                    &mut w,
                )
                .unwrap();
            }
            w.flush().unwrap();
            drop(w);
            let out = td.path().join("out.vcf");
            let words = gloci_words(cmd);
            let line = words.join(" ");
            let args: Vec<String> = words
                .iter()
                .map(|w| {
                    w.replace("{gloci}", gloci.to_str().unwrap())
                        .replace("{in}", input.to_str().unwrap())
                        .replace("{out}", out.to_str().unwrap())
                })
                .collect();
            let mut command = Command::new(&args[0]);
            command.args(&args[1..]);
            // A command that writes to stdout ends in `> {out}`.
            if cmd[cmd.find("{gloci}").unwrap()..].contains("> {out}") {
                command.stdout(std::fs::File::create(&out).unwrap());
            }
            let r = command.output().unwrap();
            assert!(
                r.status.success(),
                "case `{}`: {}",
                c.name,
                String::from_utf8_lossy(&r.stderr)
            );
            if cmd.contains("vcf stats") {
                // Every requested statistic must have rows to time.
                for ext in ["frq", "lmiss", "imiss", "het", "TsTv.summary"] {
                    let f = format!("{}.{ext}", out.display());
                    let lines = std::fs::read_to_string(&f)
                        .unwrap_or_else(|e| panic!("case `{}`: no {f}: {e}", c.name))
                        .lines()
                        .count();
                    assert!(
                        lines >= 2,
                        "case `{}`: `{}` wrote no rows to {f}",
                        c.name,
                        line.trim()
                    );
                }
            } else if cmd.contains("vcf filter") {
                let kept = std::fs::read_to_string(&out)
                    .unwrap()
                    .lines()
                    .filter(|l| !l.starts_with('#'))
                    .count();
                assert!(
                    kept > 0,
                    "case `{}`: `{}` keeps no record of its input, so it would time nothing",
                    c.name,
                    line.trim()
                );
            } else {
                let len = std::fs::metadata(&out).map(|m| m.len()).unwrap_or(0);
                assert!(
                    len > 0,
                    "case `{}`: `{}` wrote nothing",
                    c.name,
                    line.trim()
                );
            }
            checked += 1;
        }
        assert!(checked >= 6, "only {checked} vcf cases checked");
    }
}
