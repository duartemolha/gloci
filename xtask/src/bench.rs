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
    /// A paper dataset (D95): fetched only by `cargo xtask paper fetch`
    /// or `paper/fetch-datasets.sh`, never by CI.
    #[serde(default)]
    pub paper: bool,
    /// What the dataset is and which release (paper table T1).
    #[serde(default)]
    pub description: Option<String>,
    #[serde(default)]
    pub release: Option<String>,
}

/// An input made from datasets or other derived inputs by a command run with
/// `LC_ALL=C` in the bench image (D95). The result is `<cache>/<name>`.
#[derive(Deserialize, Clone, Debug)]
pub struct Derived {
    pub name: String,
    /// Each entry names a dataset, a derived input, or `repo:<path>`.
    pub from: Vec<String>,
    /// A shell command. `{in}` is the first input, `{in2}`, `{in3}` the
    /// next ones, `{out}` the result.
    pub command: String,
    /// sha256 of the result; empty until recorded with `paper fetch --record`.
    #[serde(default)]
    pub sha256: String,
}

/// Every dataset and derived input of `datasets.toml`.
pub struct Catalog {
    pub datasets: Vec<Dataset>,
    pub derived: Vec<Derived>,
}

impl Catalog {
    pub fn parse(text: &str) -> Result<Catalog> {
        let f = toml::from_str::<DatasetsFile>(text)?;
        Ok(Catalog {
            datasets: f.dataset,
            derived: f.derived,
        })
    }
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
    #[serde(default)]
    dataset: Vec<Dataset>,
    #[serde(default)]
    derived: Vec<Derived>,
}

pub fn parse_datasets(text: &str) -> Result<Vec<Dataset>> {
    Ok(toml::from_str::<DatasetsFile>(text)?.dataset)
}

#[derive(Deserialize)]
struct CaseFile {
    case: Vec<Case>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct Case {
    pub name: String,
    /// Arguments of `cargo xtask gen`, e.g. `bed --records 10000000 --seed 1`.
    #[serde(default)]
    pub setup: Option<String>,
    /// Dataset from `datasets.toml`, e.g. `1kg-30x-chr20`.
    #[serde(default)]
    pub dataset: Option<String>,
    /// `bed` (default), `bed.gz`, `vcf.gz`, `bam`, `cram`, `fa`, etc.
    #[serde(default = "default_ext")]
    pub input: String,
    /// Tool name -> shell command with `{gloci}`, `{in}`, `{in_a}`, `{in_b}`
    /// and `{out}`; alignment and FASTA setups also give `{ref}` (the
    /// reference FASTA) and `{ins}` (every BAM of `--parts`), and `{in_b}` is
    /// their BED of target windows or intervals.
    pub commands: BTreeMap<String, String>,
    pub targets: Vec<String>,
    /// Timed runs per tool (after one warm-up), judged on the median; more for I/O-heavy
    /// cases (cat, maskfasta, sort output), whose wall clock is noisy on a slow disk.
    #[serde(default = "default_runs")]
    pub runs: u32,
    /// `windows`: `{in_b}` of an alignment setup is every 30th window of
    /// `gloci windows` over its reference, sized for about 10,000 targets
    /// (the `bedcov` cases), instead of the synthesised target intervals.
    /// `exome`: 200,000 targets of 100 bp spread evenly over the reference
    /// (fewer on a reference too short to keep them 200 bp apart).
    #[serde(default)]
    pub in_b: Option<String>,
    /// Real inputs (paper cases): `a`, `b`, `c` and `genome` name a dataset
    /// or derived input of `datasets.toml`, found in the bench cache and
    /// bound to `{in_a}`, `{in_b}`, `{in_c}` and `{genome}`.
    #[serde(default)]
    pub inputs: BTreeMap<String, String>,
    /// Row label of the paper table.
    #[serde(default)]
    pub label: Option<String>,
    /// Paper table the case belongs to.
    #[serde(default)]
    pub table: Option<String>,
    /// Threads of the gloci command (recorded with the results).
    #[serde(default)]
    pub threads: Option<u32>,
    /// Commands whose outputs must agree (`[[case.same]]`): the runner runs each once into
    /// its own file and compares them before it times anything ([`crate::verify`]).
    #[serde(default)]
    pub same: Vec<crate::verify::Same>,
    /// Why a command is compared with no other (`[case.not_compared]`: tool = "reason").
    #[serde(default)]
    pub not_compared: BTreeMap<String, String>,
}

/// When the runner checks that the commands of a case agree ([`crate::verify`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Verify {
    /// Before a case is timed; a case whose commands disagree is not timed, and the run
    /// fails.
    #[default]
    Before,
    /// Check every case and time nothing.
    Only,
    /// Time without checking.
    Skip,
}

/// Timed runs per tool unless a case sets `runs`.
pub const DEFAULT_RUNS: u32 = 5;

fn default_ext() -> String {
    "bed".into()
}

fn default_runs() -> u32 {
    DEFAULT_RUNS
}

pub fn parse_cases(text: &str) -> Result<Vec<Case>> {
    let cases = toml::from_str::<CaseFile>(text)?.case;
    for c in &cases {
        ensure!(c.runs >= 1, "{}: runs must be at least 1", c.name);
        crate::verify::check_schema(&c.name, &c.same, &c.not_compared, &c.commands)?;
    }
    Ok(cases)
}

/// Replace the value after `--records` in a setup string. A `bedpe` setup has
/// `--pairs` instead, and its B counts (`--b-pairs`, `--b-bed`) scale with it.
pub fn override_records(setup: &str, records: u64) -> String {
    let mut t: Vec<String> = setup.split_whitespace().map(String::from).collect();
    let value = |t: &[String], flag: &str| -> Option<(usize, u64)> {
        let i = t.iter().position(|x| x == flag)?;
        Some((i + 1, t.get(i + 1)?.parse().ok()?))
    };
    if let Some((i, _)) = value(&t, "--records") {
        t[i] = records.to_string();
    }
    if let Some((i, was)) = value(&t, "--pairs") {
        t[i] = records.to_string();
        for flag in ["--b-pairs", "--b-bed"] {
            if let Some((j, b)) = value(&t, flag) {
                t[j] = (b as u128 * u128::from(records) / u128::from(was.max(1))).to_string();
            }
        }
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

/// Replace `{in_c}`, `{genome}` and `{python}` (the placeholders of cases with
/// `inputs`). Run before `substitute_pair`.
pub fn substitute_extra(cmd: &str, in_c: &str, genome: &str, python: &str) -> String {
    cmd.replace("{in_c}", in_c)
        .replace("{genome}", genome)
        .replace("{python}", python)
}

/// Peak resident set size in KiB from the `-o` file of GNU `time -f %M`: the
/// last line that is a bare integer (a failing command adds a "Command exited
/// with non-zero status N" line before it).
pub fn parse_rss_kb(text: &str) -> Result<u64> {
    text.lines()
        .rev()
        .map(str::trim)
        .find(|l| !l.is_empty())
        .and_then(|l| l.parse::<u64>().ok())
        // A running process has pages: 0 is a broken `time` build (conda-forge's
        // reports 0 on aarch64), not a measurement.
        .filter(|&kb| kb > 0)
        .ok_or_else(|| anyhow!("no peak RSS (KiB) in GNU time output: {text:?}"))
}

/// The `<case>.rss.json` document: one row per command, `max_rss_kb` null
/// when the measurement failed.
pub fn rss_json(case: &Case, rows: &[(String, Option<u64>)]) -> serde_json::Value {
    serde_json::json!({
        "case": case.name,
        "label": case.label,
        "table": case.table,
        "threads": case.threads,
        "unit": "KiB",
        "results": rows.iter().map(|(tool, kb)| serde_json::json!({
            "command": tool,
            "max_rss_kb": kb,
        })).collect::<Vec<_>>(),
    })
}

/// Read the `max_rss_kb` of every command from a `<case>.rss.json`.
pub fn parse_rss_json(json: &str) -> Result<Vec<(String, Option<u64>)>> {
    let v: serde_json::Value = serde_json::from_str(json)?;
    let rows = v["results"]
        .as_array()
        .ok_or_else(|| anyhow!("rss json has no `results` array"))?;
    rows.iter()
        .map(|r| {
            let name = r["command"]
                .as_str()
                .ok_or_else(|| anyhow!("rss row without `command`"))?;
            Ok((name.to_string(), r["max_rss_kb"].as_u64()))
        })
        .collect()
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

/// Where an input of a derived command is: a dataset's or derived input's
/// cache path, or a file of the repository.
pub fn input_path(cat: &Catalog, from: &str, cache_dir: &Path) -> Result<PathBuf> {
    if let Some(p) = from.strip_prefix("repo:") {
        return Ok(workspace_root().join(p));
    }
    if let Some(d) = cat.datasets.iter().find(|d| d.name == from) {
        let filename = d.url.split('/').next_back().unwrap_or(&d.name);
        return Ok(cache_dir.join(filename));
    }
    if cat.derived.iter().any(|d| d.name == from) {
        return Ok(cache_dir.join(from));
    }
    bail!("`{from}` is neither a dataset, a derived input nor repo:<path>")
}

/// `s` as one single-quoted shell word.
pub fn sh_quote(s: &str) -> String {
    format!("'{}'", s.replace('\'', "'\\''"))
}

/// The index files a derived input must have beside it, by its name.
pub fn expected_sidecars(name: &str) -> &'static [&'static str] {
    if name.ends_with(".bam") {
        &["bai"]
    } else if name.ends_with(".cram") {
        &["crai"]
    } else if name.ends_with(".bcf") {
        &["csi"]
    } else if name.ends_with(".vcf.gz") || name.ends_with(".bed.gz") {
        &["tbi"]
    } else {
        &[]
    }
}

/// `command` with `{in}`, `{in2}`... and `{out}` replaced.
pub fn derived_command(command: &str, inputs: &[String], out: &str) -> String {
    let mut c = command.replace("{out}", out);
    for (i, p) in inputs.iter().enumerate().rev() {
        let key = if i == 0 {
            "{in}".to_string()
        } else {
            format!("{{in{}}}", i + 1)
        };
        c = c.replace(&key, p);
    }
    c
}

/// Set `sha256 = "<hex>"` in the `[[derived]]` block named `name`.
pub fn record_sha256(text: &str, name: &str, hex: &str) -> Result<String> {
    let mut out = Vec::new();
    let (mut in_block, mut done) = (false, false);
    for line in text.lines() {
        if line.starts_with("[[") {
            in_block = false;
        }
        if line.trim() == format!("name = \"{name}\"") {
            in_block = true;
        }
        if in_block && !done && line.trim_start().starts_with("sha256 =") {
            out.push(format!("sha256 = \"{hex}\""));
            done = true;
        } else {
            out.push(line.to_string());
        }
    }
    ensure!(
        done,
        "no `sha256 =` line in the derived block `{name}` of datasets.toml"
    );
    Ok(out.join("\n") + "\n")
}

/// Build the derived input `d` (and what it needs) in `cache_dir`.
///
/// The command runs with `LC_ALL=C` under `bash -o pipefail`, in `image`
/// (with the cache and the repository mounted) or on the host when `image`
/// is `None`. An empty recorded sha256 is an error unless `record` is set,
/// which writes the digest into datasets.toml.
pub fn ensure_derived(
    cat: &Catalog,
    d: &Derived,
    cache_dir: &Path,
    image: Option<&str>,
    record: bool,
) -> Result<PathBuf> {
    std::fs::create_dir_all(cache_dir)?;
    let out = cache_dir.join(&d.name);
    let mut inputs = Vec::new();
    for f in &d.from {
        if let Some(dep) = cat.derived.iter().find(|x| x.name == *f) {
            ensure_derived(cat, dep, cache_dir, image, record)?;
        } else if let Some(ds) = cat.datasets.iter().find(|x| x.name == *f) {
            ensure_dataset(ds, cache_dir)?;
        }
        inputs.push(input_path(cat, f, cache_dir)?);
    }
    let sidecar = |ext: &str| cache_dir.join(format!("{}.{ext}", d.name));
    let complete = out.exists()
        && expected_sidecars(&d.name)
            .iter()
            .all(|e| sidecar(e).exists());
    if !complete {
        eprintln!("deriving {}", d.name);
        // A stale partial result would make `bcftools index` or `tabix` refuse.
        let prefix = format!("{}.partial", d.name);
        for e in std::fs::read_dir(cache_dir)?.flatten() {
            if e.file_name().to_string_lossy().starts_with(&prefix) {
                std::fs::remove_file(e.path())?;
            }
        }
        let _ = std::fs::remove_file(&out);
        let tmp = cache_dir.join(&prefix);
        let ins: Vec<String> = inputs
            .iter()
            .map(|p| sh_quote(&p.display().to_string()))
            .collect();
        let status = match image {
            None => Command::new("bash")
                .args(["-o", "pipefail", "-c"])
                .arg(derived_command(
                    &d.command,
                    &ins,
                    &sh_quote(&tmp.display().to_string()),
                ))
                .env("LC_ALL", "C")
                .status(),
            Some(img) => {
                let root = workspace_root();
                let cache_abs = std::fs::canonicalize(cache_dir)?;
                let map = |s: &str| {
                    s.replace(&cache_abs.display().to_string(), "/cache")
                        .replace(&root.display().to_string(), "/repo")
                };
                let ins: Vec<String> = inputs
                    .iter()
                    .map(|p| sh_quote(&map(&p.display().to_string())))
                    .collect();
                let id = |flag: &str| -> Result<String> {
                    let o = Command::new("id")
                        .arg(flag)
                        .output()
                        .context("running id")?;
                    Ok(String::from_utf8_lossy(&o.stdout).trim().to_string())
                };
                let user = format!("{}:{}", id("-u")?, id("-g")?);
                Command::new("docker")
                    .args(["run", "--rm", "--user", &user, "-e", "LC_ALL=C", "-v"])
                    .arg(format!("{}:/cache", cache_abs.display()))
                    .arg("-v")
                    .arg(format!("{}:/repo:ro", root.display()))
                    .arg(img)
                    .args(["bash", "-o", "pipefail", "-c"])
                    .arg(derived_command(
                        &d.command,
                        &ins,
                        &sh_quote(&format!("/cache/{prefix}")),
                    ))
                    .status()
            }
        }
        .with_context(|| format!("failed to run the command of derived input {}", d.name))?;
        ensure!(
            status.success(),
            "the command of derived input {} failed ({status}); fix it in datasets.toml",
            d.name
        );
        // Index files first, the result last: an interrupted rename leaves no
        // main file, so the next run derives again.
        let pre = format!("{prefix}.");
        for e in std::fs::read_dir(cache_dir)?.flatten() {
            let n = e.file_name().to_string_lossy().into_owned();
            if let Some(ext) = n.strip_prefix(&pre) {
                std::fs::rename(e.path(), sidecar(ext))?;
            }
        }
        std::fs::rename(&tmp, &out)?;
        for e in expected_sidecars(&d.name) {
            ensure!(
                sidecar(e).exists(),
                "derived input {} did not write its .{e} index; fix its command in datasets.toml",
                d.name
            );
        }
    }
    let hex = digest_file(&out, "sha256")?;
    if d.sha256.is_empty() {
        ensure!(
            record,
            "derived input {} has no sha256 yet; run `cargo xtask paper fetch --record` once and commit datasets.toml",
            d.name
        );
        let path = workspace_root().join("xtask/bench/datasets.toml");
        let text = std::fs::read_to_string(&path)?;
        std::fs::write(&path, record_sha256(&text, &d.name, &hex)?)?;
        eprintln!("recorded sha256 of {}: {hex}", d.name);
    } else {
        ensure!(
            hex.eq_ignore_ascii_case(&d.sha256),
            "sha256 mismatch for derived input {}: expected {}, got {hex}; delete {} and derive again",
            d.name,
            d.sha256,
            out.display()
        );
    }
    Ok(out)
}

pub fn ensure_dataset(ds: &Dataset, cache_dir: &Path) -> Result<PathBuf> {
    std::fs::create_dir_all(cache_dir)?;
    let filename = ds.url.split('/').next_back().unwrap_or(&ds.name);
    let target_path = cache_dir.join(filename);
    if !target_path.exists() {
        eprintln!("fetching dataset {} from {}", ds.name, ds.url);
        let tmp = cache_dir.join(format!("{filename}.tmp"));
        let st = Command::new("curl")
            .args(["-fsSL", "--retry", "3", "--proto", "=https", "-o"])
            .arg(&tmp)
            .arg(&ds.url)
            .status()
            .with_context(|| format!("failed to run curl to download {}", ds.url))?;
        ensure!(st.success(), "curl failed to download {}", ds.url);
        std::fs::rename(&tmp, &target_path)?;
    }
    let (algorithm, want) = ds.checksum()?;
    let hash = digest_file(&target_path, algorithm)?;
    if !hash.eq_ignore_ascii_case(want) {
        let bad = cache_dir.join(format!("{filename}.bad"));
        std::fs::rename(&target_path, &bad)?;
        bail!(
            "{algorithm} mismatch for {}: expected {want}, got {hash}; moved the file to {} (delete it and fetch again)",
            ds.name,
            bad.display()
        );
    }
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
    /// Write results here instead of `<out>/<version>` (the paper run).
    pub results: Option<PathBuf>,
    /// A Python wheel of the version under test: installed with `pip install
    /// --no-deps` (in the container), and what the `gloci-py` commands import.
    pub wheel: Option<PathBuf>,
    /// Whether the outputs of each case are compared first.
    pub verify: Verify,
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
    /// The second input: B of a pair, or the BED of target windows or
    /// intervals of an alignment or FASTA setup.
    in_b: Option<PathBuf>,
    /// The reference FASTA of an alignment setup (`{ref}`).
    reference: Option<PathBuf>,
    /// Every file of a multi-part alignment setup (`{ins}`).
    ins: Vec<PathBuf>,
    /// The third input (`{in_c}`) of a case with `inputs`.
    in_c: Option<PathBuf>,
    /// The chromosome sizes file (`{genome}`) of a case with `inputs`.
    genome: Option<PathBuf>,
}

impl PreparedInputs {
    fn new(in_a: PathBuf, in_b: Option<PathBuf>) -> Self {
        PreparedInputs {
            in_a,
            in_b,
            reference: None,
            ins: Vec::new(),
            in_c: None,
            genome: None,
        }
    }
}

/// The value after `flag` in a setup's words, parsed.
fn setup_value<T: std::str::FromStr>(words: &[String], flag: &str) -> Option<T> {
    let i = words.iter().position(|w| w == flag)?;
    words.get(i + 1)?.parse().ok()
}

/// Replace (or append) `flag value` in a setup's words.
fn set_setup_value(words: &mut Vec<String>, flag: &str, value: String) {
    match words.iter().position(|w| w == flag) {
        Some(i) if i + 1 < words.len() => words[i + 1] = value,
        _ => {
            words.push(flag.to_string());
            words.push(value);
        }
    }
}

/// Run `tool args` (paths under `dir`) on the host when the tool is installed,
/// else in the benchmark image with `dir` mounted at `/data`.
fn run_tool(tool: &str, args: &[String], dir: &Path, image: Option<&str>) -> Result<()> {
    let host = Command::new("which")
        .arg(tool)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if host {
        return run_ok(Command::new(tool).args(args));
    }
    let image = image.ok_or_else(|| {
        anyhow!("{tool} is neither installed nor is the benchmark image available")
    })?;
    let mut c = Command::new("docker");
    c.args([
        "run",
        "--rm",
        "--user",
        &format!("{}:{}", id_of("-u")?, id_of("-g")?),
    ])
    .args(["-e", "HOME=/tmp", "-v"])
    .arg(format!("{}:/data", dir.display()))
    .arg(image)
    .arg(tool)
    .args(args.iter().map(|a| match Path::new(a).strip_prefix(dir) {
        Ok(rest) => format!("/data/{}", rest.display()),
        Err(_) => a.clone(),
    }));
    run_ok(&mut c)
}

/// Generate (or reuse) the inputs of a `sam ...` setup: the BAM (or CRAM) of
/// the case's `input`, its index, the reference FASTA and a BED of target
/// windows. `--parts N` makes N coordinate-sorted BAMs of one reference
/// (`ins`). The files are named by the setup alone, so cases with one setup
/// share them.
fn prepare_sam(
    case: &Case,
    setup: &str,
    gloci: &Path,
    dir: &Path,
    image: Option<&str>,
) -> Result<PreparedInputs> {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    setup.hash(&mut h);
    let stem = format!("sam-{:08x}", h.finish() as u32);
    let mut words: Vec<String> = setup.split_whitespace().map(String::from).collect();
    let parts: u64 = setup_value(&words, "--parts").unwrap_or(1).max(1);
    if let Some(i) = words.iter().position(|w| w == "--parts") {
        words.drain(i..(i + 2).min(words.len()));
    }
    let shuffled = words.iter().any(|w| w == "--shuffled");
    let by_name = words.iter().any(|w| w == "--by-name");
    ensure!(
        parts == 1 || !(shuffled || by_name),
        "--parts does not combine with --shuffled or --by-name"
    );
    let records: u64 = setup_value(&words, "--records").unwrap_or(1000);
    let depth: f64 = setup_value(&words, "--depth").unwrap_or(30.0);
    if parts > 1 {
        // Every part has the chromosome length of the whole.
        set_setup_value(&mut words, "--records", (records / parts).to_string());
        set_setup_value(&mut words, "--depth", (depth / parts as f64).to_string());
    }
    let reference = dir.join(format!("{stem}.ref.fa"));
    let targets = dir.join(format!("{stem}.targets.bed"));
    let bam_of = |k: u64| {
        if parts == 1 {
            dir.join(format!("{stem}.bam"))
        } else {
            dir.join(format!("{stem}.part{k}.bam"))
        }
    };
    let bams: Vec<PathBuf> = (0..parts).map(bam_of).collect();
    if !bams.iter().all(|b| b.exists()) || !reference.exists() || !targets.exists() {
        for (k, bam) in bams.iter().enumerate() {
            if bam.exists() && k > 0 {
                continue;
            }
            eprintln!("generating {}", bam.display());
            let mut args = words.clone();
            if parts > 1 {
                args.push("--part".into());
                args.push(k.to_string());
            }
            let (tmp_ref, tmp_targets) = (
                dir.join(format!("{stem}.ref.fa.partial")),
                dir.join(format!("{stem}.targets.bed.partial")),
            );
            if k == 0 {
                args.extend([
                    "--ref-out".into(),
                    tmp_ref.to_string_lossy().into_owned(),
                    "--targets-out".into(),
                    tmp_targets.to_string_lossy().into_owned(),
                ]);
            }
            let tmp_bam = dir.join(format!("{stem}.partial.bam"));
            // The SAM text streams into gloci view: it is never on disk.
            let mut child = Command::new(gloci)
                .args(["view", "-", "-o"])
                .arg(&tmp_bam)
                .stdin(std::process::Stdio::piped())
                .spawn()
                .context("running gloci view")?;
            let written = {
                let stdin = child.stdin.take().ok_or_else(|| anyhow!("no stdin"))?;
                let mut w = std::io::BufWriter::with_capacity(1 << 20, stdin);
                crate::synth::run_args(&args, &mut w)
            };
            let st = child.wait()?;
            written?;
            ensure!(st.success(), "gloci view of the generated SAM failed: {st}");
            std::fs::rename(&tmp_bam, bam)?;
            if k == 0 {
                std::fs::rename(&tmp_ref, &reference)?;
                std::fs::rename(&tmp_targets, &targets)?;
            }
            // Only coordinate-sorted BAMs have an index.
            if !shuffled && !by_name {
                run_ok(Command::new(gloci).arg("index").arg(bam))?;
            }
        }
    }
    if !PathBuf::from(format!("{}.fai", reference.display())).exists() {
        run_ok(Command::new(gloci).arg("index").arg(&reference))?;
    }
    let mut in_b = targets;
    let fai = PathBuf::from(format!("{}.fai", reference.display()));
    match case.in_b.as_deref() {
        Some("windows") => {
            in_b = dir.join(format!("{stem}.windows.bed"));
            if !in_b.exists() {
                write_window_targets(gloci, &fai, &in_b)?;
            }
        }
        Some("exome") => {
            in_b = dir.join(format!("{stem}.exome.bed"));
            if !in_b.exists() {
                write_exome_targets(&fai, &in_b)?;
            }
        }
        Some(other) => bail!("{}: unknown in_b `{other}`", case.name),
        None => {}
    }
    let mut in_a = bams[0].clone();
    if case.input == "cram" {
        let cram = dir.join(format!("{stem}.cram"));
        if !cram.exists() {
            eprintln!("writing {} with samtools", cram.display());
            let tmp = dir.join(format!("{stem}.partial.cram"));
            run_tool(
                "samtools",
                &[
                    "view".into(),
                    "-C".into(),
                    "-T".into(),
                    reference.to_string_lossy().into_owned(),
                    "-o".into(),
                    tmp.to_string_lossy().into_owned(),
                    in_a.to_string_lossy().into_owned(),
                ],
                dir,
                image,
            )?;
            std::fs::rename(&tmp, &cram)?;
        }
        if !PathBuf::from(format!("{}.crai", cram.display())).exists() {
            run_tool(
                "samtools",
                &["index".into(), cram.to_string_lossy().into_owned()],
                dir,
                image,
            )?;
        }
        in_a = cram;
    }
    Ok(PreparedInputs {
        in_a,
        in_b: Some(in_b),
        reference: Some(reference),
        ins: bams,
        in_c: None,
        genome: None,
    })
}

/// A BED of about 10,000 windows: every 30th window of `gloci windows` over
/// the `.fai`, in windows of `total length / 300,000` bases (at least 1 kb).
fn write_window_targets(gloci: &Path, fai: &Path, bed: &Path) -> Result<()> {
    let total: u64 = std::fs::read_to_string(fai)?
        .lines()
        .filter_map(|l| l.split('\t').nth(1)?.parse::<u64>().ok())
        .sum();
    let size = (total / 300_000).max(1000);
    let out = Command::new(gloci)
        .args(["windows", "-g"])
        .arg(fai)
        .args(["-w", &size.to_string()])
        .output()
        .context("running gloci windows")?;
    ensure!(
        out.status.success(),
        "gloci windows failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8(out.stdout)?;
    let mut kept = String::new();
    for l in text.lines().step_by(30) {
        kept.push_str(l);
        kept.push('\n');
    }
    ensure!(!kept.is_empty(), "no window targets for {}", fai.display());
    let tmp = bed.with_extension("bed.partial");
    std::fs::write(&tmp, kept)?;
    std::fs::rename(&tmp, bed)?;
    Ok(())
}

/// An exome-like BED: 200,000 targets of 100 bp, one every `total / 200,000` bases of the
/// reference in `.fai` order (at least 200 bases apart, so a short reference gets fewer),
/// each inside its chromosome.
fn write_exome_targets(fai: &Path, bed: &Path) -> Result<()> {
    let chroms: Vec<(String, u64)> = std::fs::read_to_string(fai)?
        .lines()
        .filter_map(|l| {
            let mut f = l.split('\t');
            Some((f.next()?.to_string(), f.next()?.parse::<u64>().ok()?))
        })
        .collect();
    let total: u64 = chroms.iter().map(|(_, n)| n).sum();
    let step = (total / 200_000).max(200);
    let mut kept = String::new();
    let mut n = 0u64;
    for (name, len) in &chroms {
        let mut pos = step / 2;
        while pos + 100 <= *len && n < 200_000 {
            kept.push_str(&format!("{name}\t{pos}\t{}\n", pos + 100));
            n += 1;
            pos += step;
        }
    }
    ensure!(n > 0, "no exome targets for {}", fai.display());
    let tmp = bed.with_extension("bed.partial");
    std::fs::write(&tmp, kept)?;
    std::fs::rename(&tmp, bed)?;
    Ok(())
}

/// Generate (or reuse) the FASTA, its `.fai` and the BED of intervals of a
/// `fasta ...` setup.
fn prepare_fasta(setup: &str, gloci: &Path, dir: &Path) -> Result<PreparedInputs> {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    setup.hash(&mut h);
    let stem = format!("fa-{:08x}", h.finish() as u32);
    let (fa, bed) = (
        dir.join(format!("{stem}.fa")),
        dir.join(format!("{stem}.bed")),
    );
    if !fa.exists() || !bed.exists() {
        eprintln!("generating {}", fa.display());
        let (tmp_fa, tmp_bed) = (
            dir.join(format!("{stem}.fa.partial")),
            dir.join(format!("{stem}.bed.partial")),
        );
        let mut args: Vec<String> = setup.split_whitespace().map(String::from).collect();
        args.push("--bed-out".into());
        args.push(tmp_bed.to_string_lossy().into_owned());
        let mut f = std::fs::File::create(&tmp_fa)?;
        crate::synth::run_args(&args, &mut f)?;
        std::fs::rename(&tmp_fa, &fa)?;
        std::fs::rename(&tmp_bed, &bed)?;
    }
    if !PathBuf::from(format!("{}.fai", fa.display())).exists() {
        run_ok(Command::new(gloci).arg("index").arg(&fa))?;
    }
    Ok(PreparedInputs {
        reference: Some(fa.clone()),
        ..PreparedInputs::new(fa, Some(bed))
    })
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
/// Inputs of a case with `inputs`: the datasets and derived inputs of
/// `datasets.toml` where `paper fetch` left them (the bench cache, which for
/// docker runs must lie under the data directory).
fn prepare_real_inputs(case: &Case, dir: &Path) -> Result<PreparedInputs> {
    let text = std::fs::read_to_string(workspace_root().join("xtask/bench/datasets.toml"))
        .context("reading datasets.toml")?;
    let cat = Catalog::parse(&text)?;
    let cache = bench_cache_dir(dir);
    let mut found: BTreeMap<&str, PathBuf> = BTreeMap::new();
    for (slot, name) in &case.inputs {
        ensure!(
            matches!(slot.as_str(), "a" | "b" | "c" | "genome"),
            "case {}: unknown input slot `{slot}` (use a, b, c or genome)",
            case.name
        );
        let p = input_path(&cat, name, &cache)?;
        ensure!(
            p.exists(),
            "case {}: input `{name}` is missing at {}: run `cargo xtask paper fetch --data-dir {}`",
            case.name,
            p.display(),
            dir.display()
        );
        found.insert(slot.as_str(), p);
    }
    let in_a = found
        .remove("a")
        .ok_or_else(|| anyhow!("case {}: `inputs` has no `a`", case.name))?;
    let mut prepared = PreparedInputs::new(in_a, found.remove("b"));
    prepared.in_c = found.remove("c");
    prepared.genome = found.remove("genome");
    Ok(prepared)
}

fn prepare_input(
    case: &Case,
    setup: Option<&str>,
    gloci: &Path,
    dir: &Path,
    image: Option<&str>,
) -> Result<PreparedInputs> {
    std::fs::create_dir_all(dir)?;

    if !case.inputs.is_empty() {
        return prepare_real_inputs(case, dir);
    }

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
            return Ok(PreparedInputs::new(target_a, Some(target_b)));
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
            return Ok(PreparedInputs::new(target, None));
        }
    }

    let setup_str =
        setup.ok_or_else(|| anyhow!("case {} has neither setup nor dataset", case.name))?;
    match setup_str.split_whitespace().next() {
        Some("sam") => return prepare_sam(case, setup_str, gloci, dir, image),
        Some("fasta") => return prepare_fasta(setup_str, gloci, dir),
        _ => {}
    }
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    setup_str.hash(&mut h);
    let stem = format!("{}-{:08x}", case.name, h.finish() as u32);

    let kind = setup_str.split_whitespace().next().unwrap_or("");
    // Two-file setups: `pair` (two BEDs), `vcf-pair` (two VCFs) and `bedpe`
    // with --b-pairs (two BEDPEs) or --b-bed (a BEDPE and a BED).
    let (ext_a, ext_b) = match kind {
        "pair" => (Some("bed"), "bed"),
        "vcf-pair" => (Some("vcf"), "vcf"),
        "bedpe" if setup_str.contains("--b-pairs") => (Some("bedpe"), "bedpe"),
        "bedpe" if setup_str.contains("--b-bed") => (Some("bedpe"), "bed"),
        _ => (None, ""),
    };
    let is_vcf = kind == "vcf";
    let ext = if is_vcf { "vcf" } else { "bed" };

    if let Some(ext_a) = ext_a {
        let gz = if case.input.ends_with(".gz") {
            ".gz"
        } else {
            ""
        };
        let plain_a = dir.join(format!("{stem}_a.{ext_a}"));
        let plain_b = dir.join(format!("{stem}_b.{ext_b}"));
        let want_a = dir.join(format!("{stem}_a.{ext_a}{gz}"));
        let want_b = dir.join(format!("{stem}_b.{ext_b}{gz}"));
        if !want_a.exists() || !want_b.exists() {
            if !plain_a.exists() || !plain_b.exists() {
                eprintln!("generating {} and {}", plain_a.display(), plain_b.display());
                let tmp_a = plain_a.with_extension(format!("{ext_a}.partial"));
                let tmp_b = plain_b.with_extension(format!("{ext_b}.partial"));
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
                let tmp_a = dir.join(format!("{stem}_a.partial.{ext_a}{gz}"));
                let tmp_b = dir.join(format!("{stem}_b.partial.{ext_b}{gz}"));
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
        Ok(PreparedInputs::new(want_a, Some(want_b)))
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
                // An index needs sorted records: the input of a sort case is not sorted,
                // and `gloci index` refuses it.
                let sorted = setup_str.split_whitespace().any(|w| w == "--sorted");
                eprintln!(
                    "compressing{} {}",
                    if sorted { " and indexing" } else { "" },
                    want.display()
                );
                let tmp = dir.join(format!("{stem}.partial.{}", case.input));
                run_ok(
                    Command::new(gloci)
                        .arg("view")
                        .arg(&plain)
                        .arg("-o")
                        .arg(&tmp),
                )?;
                std::fs::rename(&tmp, &want)?;
                if sorted {
                    run_ok(Command::new(gloci).arg("index").arg(&want))?;
                }
                let _ = std::fs::remove_file(&plain);
            }
        }
        Ok(PreparedInputs::new(want, None))
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

/// `docker run` of the bench image as the current user, with the data,
/// results and gloci binary mounted; the caller appends the command.
/// The Python benchmark scripts, as the container sees them.
const PY_DIR_IN_CONTAINER: &str = "/bench-py";

/// `pip install --no-deps` of `wheel` into `target` (D96: the wheel brings
/// its own dependencies, which are none beyond Python).
pub fn pip_install_args(wheel: &str, target: &str) -> Vec<String> {
    [
        "python3",
        "-m",
        "pip",
        "install",
        "--no-deps",
        "--no-compile",
        "--target",
        target,
        wheel,
    ]
    .map(String::from)
    .to_vec()
}

/// What `{python}` stands for: the interpreter, with the directory the wheel
/// was installed into on its path.
pub fn python_with_wheel(python: &str, pylib: Option<&str>) -> String {
    match pylib {
        Some(dir) => format!("PYTHONPATH={dir} {python}"),
        None => python.to_string(),
    }
}

/// Install `wheel` into `<results>/.pylib`, in the container when there is
/// one. Returns that directory as the commands see it.
fn install_wheel(
    wheel: &Path,
    results: &Path,
    docker: Option<(&str, &Path, &Path)>,
    python: &str,
) -> Result<String> {
    let wheel = std::fs::canonicalize(wheel)
        .with_context(|| format!("--wheel {}: no such file", wheel.display()))?;
    let name = wheel
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .context("--wheel: no file name")?;
    let host_dir = results.join(".pylib");
    let _ = std::fs::remove_dir_all(&host_dir);
    match docker {
        Some((image, data, results_mount)) => {
            let in_c = "/results/.pylib";
            let mut c = docker_run(
                image,
                data,
                results_mount,
                &[
                    "-v".to_string(),
                    format!("{}:/wheel/{name}:ro", wheel.display()),
                ],
            )?;
            c.args(pip_install_args(&format!("/wheel/{name}"), in_c));
            run_ok(&mut c).context(
                "installing --wheel in the bench image (is it a wheel for the image's architecture and Python 3.12?)",
            )?;
            Ok(in_c.to_string())
        }
        None => {
            let mut args = pip_install_args(&wheel.to_string_lossy(), &host_dir.to_string_lossy());
            args[0] = python.to_string();
            run_ok(Command::new(&args[0]).args(&args[1..]))
                .context("installing --wheel with the host's pip")?;
            Ok(host_dir.to_string_lossy().into_owned())
        }
    }
}

/// The Python module a benchmark tool name stands for, when it is one.
fn python_module(tool: &str) -> Option<&'static str> {
    match tool {
        "pybedtools" => Some("pybedtools"),
        "polars-bio" => Some("polars_bio"),
        "pyranges1" => Some("pyranges1"),
        "bioframe" => Some("bioframe"),
        _ => None,
    }
}

fn docker_run(image: &str, data: &Path, results: &Path, extra: &[String]) -> Result<Command> {
    let mut c = Command::new("docker");
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
    .args(extra)
    .arg(image);
    Ok(c)
}

/// Peak RSS (KiB) of one run of `cmd` under GNU `time -f %M`. The command's
/// own failure does not matter here (hyperfine already ran it); a missing
/// `time` or an unreadable report is an `Err`. For a pipeline, `time` reports
/// the largest single process, not the sum.
fn measure_rss(
    cmd: &str,
    docker: Option<(&str, &Path, &Path, &[String])>,
    results: &Path,
    tag: &str,
) -> Result<u64> {
    let file = results.join(format!("{tag}.rss.txt"));
    let _ = std::fs::remove_file(&file);
    let mut c = match docker {
        Some((image, data, results, extra)) => {
            let mut c = docker_run(image, data, results, extra)?;
            c.args([
                "/usr/bin/time",
                "-f",
                "%M",
                "-o",
                &format!("/results/{tag}.rss.txt"),
            ]);
            c
        }
        None => {
            let mut c = Command::new("/usr/bin/time");
            c.args(["-f", "%M", "-o"]).arg(&file);
            c
        }
    };
    c.args(["sh", "-c", cmd])
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null());
    c.status()
        .context("running /usr/bin/time (GNU time not installed?)")?;
    let text = std::fs::read_to_string(&file).context("GNU time wrote no report")?;
    let _ = std::fs::remove_file(&file);
    parse_rss_kb(&text)
}

fn tool_binary(tool: &str) -> &str {
    match tool {
        "gnu-sort" => "sort",
        "grit-genomics" => "grit",
        "vcfkit" => "vk",
        "vcftools-5runs" | "vcftools-3runs" => "vcftools",
        "samtools-t8" | "samtools-all" => "samtools",
        "bcftools-stats" => "bcftools",
        _ => tool,
    }
}

fn check_tool_available(tool: &str, use_docker: bool, image: &str) -> bool {
    // A Python library is available when `python3` can import it.
    if let Some(module) = python_module(tool) {
        let py = std::env::var("GLOCI_PYTHON").unwrap_or_else(|_| "python3".into());
        let code = format!("import {module}");
        let mut c = if use_docker {
            let mut c = Command::new("docker");
            c.args(["run", "--rm", image, "python3"]);
            c
        } else {
            Command::new(py)
        };
        return c
            .args(["-c", &code])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
    }
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

fn run_native_bench(
    commands: &[(&str, String)],
    runs: u32,
    json_path: &Path,
) -> Result<Vec<(String, f64)>> {
    let mut results = Vec::new();
    for (tool, cmd) in commands {
        // 1 warmup run
        let _ = Command::new("sh").arg("-c").arg(cmd).status();
        // `runs` timed runs, judged on the median
        let mut times = Vec::new();
        for _ in 0..runs {
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
    let cases = parse_cases(&std::fs::read_to_string(&a.cases)?)?;
    run_cases(a, cases)
}

/// Run `cases` (those of `a.only` when set) and write the results to
/// `a.results`, or `<a.out>/<version>`.
pub fn run_cases(a: &Args, cases: Vec<Case>) -> Result<bool> {
    let gloci = resolve_gloci(a.gloci.as_deref())?;
    let version = env!("CARGO_PKG_VERSION");
    let results = a.results.clone().unwrap_or_else(|| a.out.join(version));
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

    // What the container is given besides the data and results: the binary
    // under test and the Python benchmark scripts.
    let py_dir = workspace_root().join("xtask/bench/py");
    let mounts: Vec<String> = vec![
        "-v".into(),
        format!("{}:/usr/local/bin/gloci:ro", gloci.display()),
        "-v".into(),
        format!("{}:{PY_DIR_IN_CONTAINER}:ro", py_dir.display()),
    ];
    let python = std::env::var("GLOCI_PYTHON").unwrap_or_else(|_| "python3".into());
    let wants_py = cases
        .iter()
        .filter(|c| a.only.as_ref().is_none_or(|o| *o == c.name))
        .any(|c| c.commands.contains_key("gloci-py"));
    let pylib = match (&a.wheel, wants_py) {
        (Some(w), true) => Some(install_wheel(
            w,
            &results,
            use_docker.then_some((a.image.as_str(), data.as_path(), results.as_path())),
            &python,
        )?),
        _ => None,
    };
    let python_cmd = python_with_wheel(&python, pylib.as_deref());
    let py_arg = if use_docker {
        PY_DIR_IN_CONTAINER.to_string()
    } else {
        py_dir.to_string_lossy().into_owned()
    };

    let mut all_ok = true;
    let mut ran = 0;
    // Cases whose commands do not give the same result.
    let mut disagree: Vec<String> = Vec::new();
    for case in cases
        .iter()
        .filter(|c| a.only.as_ref().is_none_or(|o| *o == c.name))
    {
        ran += 1;
        if a.verify == Verify::Only && case.same.is_empty() {
            // Nothing to compare, so nothing to generate or run either (and the Python
            // cases, which have no group, need a wheel that a comparison run has no use for).
            eprintln!("== {}: no `same` group, nothing compared", case.name);
            continue;
        }
        let setup = case.setup.as_ref().map(|s| match a.records {
            Some(n) => override_records(s, n),
            None => s.clone(),
        });
        let inputs = prepare_input(
            case,
            setup.as_deref(),
            &gloci,
            &data,
            use_docker.then_some(a.image.as_str()),
        )?;
        let out_name = format!("{}.out.{}", case.name, case.input);
        let json = results.join(format!("{}.json", case.name));

        // A path as the commands see it: under /data in the container.
        let shown = |p: &Path| -> Result<String> {
            if !use_docker {
                return Ok(p.to_string_lossy().into_owned());
            }
            match p.strip_prefix(&data) {
                Ok(rel) => Ok(format!("/data/{}", rel.display())),
                Err(_) => {
                    ensure!(
                        case.inputs.is_empty(),
                        "case {}: input {} is outside the data directory {}, which is all the container sees (unset GLOCI_BENCH_CACHE)",
                        case.name,
                        p.display(),
                        data.display()
                    );
                    Ok(format!(
                        "/data/{}",
                        p.file_name().unwrap().to_string_lossy()
                    ))
                }
            }
        };
        let opt = |p: &Option<PathBuf>| -> Result<String> {
            p.as_deref()
                .map(&shown)
                .transpose()
                .map(Option::unwrap_or_default)
        };
        let gloci_cmd = if use_docker {
            "gloci".to_string()
        } else {
            gloci.to_string_lossy().into_owned()
        };
        let in_a_arg = shown(&inputs.in_a)?;
        let in_b_arg = opt(&inputs.in_b)?;
        let in_c_arg = opt(&inputs.in_c)?;
        let genome_arg = opt(&inputs.genome)?;
        let out_arg = if use_docker {
            format!("/data/{out_name}")
        } else {
            data.join(&out_name).to_string_lossy().into_owned()
        };
        let ref_arg = opt(&inputs.reference)?;
        let ins_arg = inputs
            .ins
            .iter()
            .map(|p| shown(p))
            .collect::<Result<Vec<_>>>()?
            .join(" ");
        let fill_out = |cmd: &str, out: &str| {
            substitute_pair(
                &substitute_extra(
                    &cmd.replace("{ref}", &ref_arg)
                        .replace("{ins}", &ins_arg)
                        .replace("{py}", &py_arg),
                    &in_c_arg,
                    &genome_arg,
                    &python_cmd,
                ),
                &gloci_cmd,
                &in_a_arg,
                &in_b_arg,
                out,
            )
        };
        let fill = |cmd: &str| fill_out(cmd, &out_arg);

        // Filter commands based on tool availability
        let mut active_commands: Vec<(String, String)> = Vec::new();
        for (tool, cmd) in &case.commands {
            // `gloci-py` runs the wheel: without `--wheel` there is nothing to time.
            if tool == "gloci-py" {
                if pylib.is_some() {
                    active_commands.push((tool.clone(), cmd.clone()));
                } else {
                    eprintln!("  skipping gloci-py: no --wheel");
                }
                continue;
            }
            // `gloci` and its variants (`gloci-t8`) run the binary under test.
            if tool == "gloci" || tool.starts_with("gloci-") {
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

        // Do the commands agree? Each command of a `same` group is run once into its own
        // file and the files are compared; a case that fails is not timed.
        if a.verify != Verify::Skip {
            if case.same.is_empty() {
                eprintln!("== {}: no `same` group, nothing compared", case.name);
            } else {
                eprintln!("== {} (comparing the outputs)", case.name);
                let docker = use_docker.then_some((
                    a.image.as_str(),
                    data.as_path(),
                    results.as_path(),
                    mounts.as_slice(),
                ));
                let reports = verify_case(case, &active_commands, &data, &fill_out, docker)?;
                let mut text = String::new();
                for r in &reports {
                    eprintln!("  {}", r.line());
                    text.push_str(&r.line());
                    text.push('\n');
                }
                for (tool, why) in &case.not_compared {
                    text.push_str(&format!("not compared: {tool}: {why}\n"));
                }
                std::fs::write(results.join(format!("{}.verify.txt", case.name)), text)?;
                if reports.iter().any(|r| !r.ok()) {
                    eprintln!("  {}: not timed, its commands disagree", case.name);
                    disagree.push(case.name.clone());
                    continue;
                }
            }
        }
        if a.verify == Verify::Only {
            continue;
        }

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
            let runs = case.runs.to_string();
            let mut hf: Vec<String> = [
                "hyperfine",
                "--shell=default",
                "--warmup",
                "1",
                "--runs",
                &runs,
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
                hf.push(fill(cmd));
            }
            eprintln!("== {} ({} tools)", case.name, active_commands.len());
            let mut c;
            if use_docker {
                c = docker_run(&a.image, &data, &results, &mounts)?;
                c.args(&hf);
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
                .map(|(tool, cmd)| (tool.as_str(), fill(cmd)))
                .collect();
            run_native_bench(&substituted, case.runs, &json)?;
            remove_outputs(&data, &out_name);
        }

        // Peak RSS: one more run of every command under GNU time (D96).
        let mut rss_rows: Vec<(String, Option<u64>)> = Vec::new();
        for (tool, cmd) in &active_commands {
            remove_outputs(&data, &out_name);
            let tag = format!("{}.{tool}", case.name);
            let docker = use_docker.then_some((
                a.image.as_str(),
                data.as_path(),
                results.as_path(),
                mounts.as_slice(),
            ));
            let kb = match measure_rss(&fill(cmd), docker, &results, &tag) {
                Ok(kb) => {
                    eprintln!("  {tool:<12} peak RSS {:.1} MiB", kb as f64 / 1024.0);
                    Some(kb)
                }
                Err(e) => {
                    eprintln!("  warning: no peak RSS for {tool}: {e:#}");
                    None
                }
            };
            rss_rows.push((tool.clone(), kb));
        }
        remove_outputs(&data, &out_name);
        std::fs::write(
            results.join(format!("{}.rss.json", case.name)),
            serde_json::to_string_pretty(&rss_json(case, &rss_rows))?,
        )?;

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
    ensure!(
        disagree.is_empty(),
        "the commands of {} do not give the same result (see above, and <case>.verify.txt in {}); \
fix the commands or the case's `same` groups, or pass --no-verify to time them anyway",
        disagree.join(", "),
        results.display()
    );
    Ok(all_ok || !a.check_targets)
}

/// `cmd` through `sh -c`, in the bench container when `docker` is given.
fn shell(cmd: &str, docker: Option<(&str, &Path, &Path, &[String])>) -> Result<Command> {
    let mut c = match docker {
        Some((image, data, results, extra)) => docker_run(image, data, results, extra)?,
        None => Command::new("sh"),
    };
    if docker.is_some() {
        c.arg("sh");
    }
    c.args(["-c", cmd]);
    Ok(c)
}

/// Run `cmd` ([`shell`]) and return whether it succeeded and what it wrote to standard
/// error.
fn run_shell(cmd: &str, docker: Option<(&str, &Path, &Path, &[String])>) -> Result<(bool, String)> {
    let o = shell(cmd, docker)?
        .stdout(std::process::Stdio::null())
        .output()
        .with_context(|| format!("running {cmd}"))?;
    Ok((
        o.status.success(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    ))
}

/// Run every command of the `same` groups of `case` once, each into its own file
/// (`<case>.verify.<tool>.<ext>` in the data directory), and compare the files group by
/// group. `active` are the commands whose tools are installed; `fill_out` fills a command's
/// placeholders with a given `{out}`.
fn verify_case(
    case: &Case,
    active: &[(String, String)],
    data: &Path,
    fill_out: &dyn Fn(&str, &str) -> String,
    docker: Option<(&str, &Path, &Path, &[String])>,
) -> Result<Vec<crate::verify::GroupReport>> {
    let out_name = |tool: &str| format!("{}.verify.{tool}.{}", case.name, case.input);
    let shown = |name: &str| match docker {
        Some(_) => format!("/data/{name}"),
        None => data.join(name).to_string_lossy().into_owned(),
    };
    let mut tools: Vec<&str> = case
        .same
        .iter()
        .flat_map(|g| g.tools.iter().map(String::as_str))
        .collect();
    tools.sort_unstable();
    tools.dedup();
    let mut outputs: BTreeMap<String, PathBuf> = BTreeMap::new();
    for tool in tools {
        let Some((_, cmd)) = active.iter().find(|(t, _)| t == tool) else {
            continue;
        };
        let name = out_name(tool);
        remove_outputs(data, &name);
        let cmd = fill_out(cmd, &shown(&name));
        let (ok, stderr) = run_shell(&cmd, docker)?;
        ensure!(
            ok,
            "case {}: `{tool}` failed when run to compare its output: {cmd}\n{}",
            case.name,
            stderr.trim_end()
        );
        outputs.insert(tool.to_string(), data.join(&name));
    }
    let mut reports = Vec::new();
    for g in &case.same {
        // A binary output is compared as the records a reference tool prints for it, read
        // from that tool as it prints them.
        let output = |tool: &str| -> Result<Option<crate::verify::Output>> {
            let Some(p) = outputs.get(tool) else {
                return Ok(None);
            };
            Ok(Some(match crate::verify::decoder(&g.compare) {
                Some(decoder) if p.exists() => {
                    let name = p.file_name().unwrap().to_string_lossy();
                    crate::verify::Output::Records(shell(
                        &format!("{decoder} {}", shown(&name)),
                        docker,
                    )?)
                }
                _ => crate::verify::Output::File(p.clone()),
            }))
        };
        reports.push(crate::verify::compare_group(g, &output)?);
    }
    for tool in outputs.keys() {
        remove_outputs(data, &out_name(tool));
    }
    Ok(reports)
}

/// The words of the `{gloci}` command in `cmd`, up to the first unquoted
/// pipe, redirection or `&`.
pub fn gloci_words(cmd: &str) -> Vec<String> {
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

    /// The words of a case's gloci command: from `{gloci}` to the first
    /// unquoted pipe, redirection or `&&`, with single quotes removed.
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
    fn runs_default_and_per_case() {
        let root = crate::diff::workspace_root();
        let text = std::fs::read_to_string(root.join("xtask/bench/cases.toml")).unwrap();
        let cases = parse_cases(&text).unwrap();
        let runs = |n: &str| cases.iter().find(|c| c.name == n).unwrap().runs;
        assert_eq!(runs("sort-plain-10M"), DEFAULT_RUNS);
        for n in [
            "cat-bam-8x",
            "maskfasta-3G",
            "sort-bam-50M",
            "sort-name-bam-50M",
            "merge-bam-8x",
        ] {
            assert_eq!(runs(n), 10, "{n}");
        }
        let zero = "[[case]]\nname = \"z\"\nruns = 0\ntargets = []\n[case.commands]\n";
        assert!(parse_cases(zero).is_err());
    }

    /// Every `gloci` command line in cases.toml parses with gloci's own clap
    /// parser, so a misspelled or removed flag fails `cargo test` instead of
    /// a nightly benchmark run.
    #[test]
    fn every_gloci_bench_command_parses() {
        use clap::Parser;
        let text =
            std::fs::read_to_string(workspace_root().join("xtask/bench/cases.toml")).unwrap();
        let cases = parse_cases(&text).unwrap();
        let mut checked = 0;
        for c in &cases {
            // `gloci` and its variants (`gloci-t8`), not the Python library (`gloci-py`).
            for (_, cmd) in c
                .commands
                .iter()
                .filter(|(k, _)| k.starts_with("gloci") && *k != "gloci-py")
            {
                let line = gloci_words(cmd).join(" ");
                let args: Vec<String> = gloci_words(cmd)
                    .iter()
                    .flat_map(|w| {
                        let w = w
                            .replace("{gloci}", "gloci")
                            .replace("{in_a}", "a.bed")
                            .replace("{in_b}", "b.bed")
                            .replace("{ref}", "ref.fa")
                            .replace("{in}", "in.bed")
                            .replace("{out}", "out.bed");
                        if w == "{ins}" {
                            vec!["x.bam".to_string(), "y.bam".to_string()]
                        } else {
                            vec![w]
                        }
                    })
                    .collect();
                if let Err(e) = gloci::cli::Cli::try_parse_from(&args) {
                    panic!("case `{}`: `{}` does not parse:\n{e}", c.name, line.trim());
                }
                checked += 1;
            }
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
            // P3a cases
            "convert-bam-50M",
            "genomecov-bam-50M",
            "coverage-bam-targets",
            "intersect-bam-any",
            "coverage-cram-targets",
            "getfasta-1M",
            "nuc-1M",
            "maskfasta-3G",
            "index-bam-50M",
            "sort-bam-50M",
            "sort-name-bam-50M",
            "merge-bam-8x",
            "cat-bam-8x",
            // P3b cases
            "jaccard-plain-10M",
            "fisher-plain-10M",
            "reldist-plain-10M",
            "sample-plain-10M",
            "shuffle-plain-10M",
            "random-10M",
            "pair-intersect-1M-vs-10M",
            "pair-pair-1M",
            "convert-bedpe-50M",
            "vcf-compare-synthetic-2500",
            // P3c cases
            "flagstat-bam-50M",
            "idxstats-bam-200k-sq",
            "view-count-bam-50M",
            "depth-bam-50M",
            "depth-bam-50M-t8",
            "bedcov-bam-10k-targets",
            "bedcov-bam-200k-exome",
            "faidx-10k-regions",
            "dict-3G",
            "reheader-bam-50M",
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

    /// The input of a sort case is not sorted: it is compressed and not indexed (`gloci index`
    /// refuses it, which stopped the whole run at `sort-bgzf-10M`).
    #[test]
    fn prepare_input_does_not_index_an_unsorted_bgzf_input() {
        let gloci = match std::env::var_os("GLOCI_BIN") {
            Some(b) => PathBuf::from(b),
            None => built_gloci("debug").unwrap(),
        };
        let td = tempfile::tempdir().unwrap();
        let c = Case {
            name: "test-sort-bgzf".into(),
            setup: Some("bed --records 200 --chroms 3 --seed 1 --profile uniform".into()),
            dataset: None,
            input: "bed.gz".into(),
            commands: BTreeMap::new(),
            targets: Vec::new(),
            runs: DEFAULT_RUNS,
            in_b: None,
            inputs: BTreeMap::new(),
            label: None,
            table: None,
            threads: None,
            same: Vec::new(),
            not_compared: BTreeMap::new(),
        };
        let inputs = prepare_input(&c, c.setup.as_deref(), &gloci, td.path(), None).unwrap();
        assert!(inputs.in_a.exists());
        for ext in ["tbi", "csi"] {
            assert!(!PathBuf::from(format!("{}.{ext}", inputs.in_a.display())).exists());
        }
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
            runs: DEFAULT_RUNS,
            in_b: None,
            inputs: BTreeMap::new(),
            label: None,
            table: None,
            threads: None,
            same: Vec::new(),
            not_compared: BTreeMap::new(),
        };
        let inputs = prepare_input(&c, c.setup.as_deref(), &gloci, td.path(), None).unwrap();
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

    /// Every P3a and P3c case's gloci command runs on a tiny generated version of its
    /// own setup (a SAM, or a FASTA with its BED), succeeds and writes
    /// output. A CRAM case reads the BAM here (samtools makes the CRAM only
    /// in the benchmark image); `{ref}` is still the generated reference.
    #[test]
    fn gloci_bench_commands_run_on_small_data() {
        use std::process::Stdio;
        let gloci = match std::env::var_os("GLOCI_BIN") {
            Some(b) => PathBuf::from(b),
            None => built_gloci("debug").unwrap(),
        };
        let text =
            std::fs::read_to_string(workspace_root().join("xtask/bench/cases.toml")).unwrap();
        let td = tempfile::tempdir().unwrap();
        let mut checked = 0;
        for c in parse_cases(&text).unwrap() {
            let Some(setup) = c.setup.as_deref() else {
                continue;
            };
            let kind = setup.split_whitespace().next();
            if kind != Some("sam") && kind != Some("fasta") {
                continue;
            }
            let mut words: Vec<String> = setup.split_whitespace().map(String::from).collect();
            set_setup_value(&mut words, "--records", "4000".into());
            if kind == Some("fasta") {
                set_setup_value(&mut words, "--chrom-len", "60000".into());
            }
            if setup_value::<u64>(&words, "--chroms").unwrap_or(0) > 1000 {
                // A 200k-contig reference is 10 Gbases: keep the shape, not the size.
                set_setup_value(&mut words, "--chroms", "300".into());
            }
            if c.name.starts_with("faidx-") {
                set_setup_value(&mut words, "--records", "300".into());
            }
            let small = Case {
                input: if c.input == "cram" {
                    "bam".into()
                } else {
                    c.input.clone()
                },
                ..c.clone()
            };
            let inputs =
                prepare_input(&small, Some(&words.join(" ")), &gloci, td.path(), None).unwrap();
            let out = td.path().join(format!("{}.out.{}", c.name, c.input));
            let cmd = c.commands.get("gloci").expect("a gloci command");
            let path = |p: &Path| p.to_string_lossy().into_owned();
            let args: Vec<String> = gloci_words(cmd)
                .iter()
                .flat_map(|w| {
                    if w == "{ins}" {
                        return inputs.ins.iter().map(|p| path(p)).collect();
                    }
                    vec![
                        w.replace("{gloci}", &path(&gloci))
                            .replace(
                                "{in_b}",
                                &inputs.in_b.as_deref().map(path).unwrap_or_default(),
                            )
                            .replace(
                                "{ref}",
                                &inputs.reference.as_deref().map(path).unwrap_or_default(),
                            )
                            .replace("{in}", &path(&inputs.in_a))
                            .replace("{out}", &path(&out)),
                    ]
                })
                .collect();
            let mut command = Command::new(&args[0]);
            command.args(&args[1..]).stderr(Stdio::piped());
            if cmd[cmd.find("{gloci}").unwrap()..].contains("> {out}") {
                command.stdout(std::fs::File::create(&out).unwrap());
            }
            let r = command.output().unwrap();
            assert!(
                r.status.success(),
                "case `{}`: `{}`: {}",
                c.name,
                args.join(" "),
                String::from_utf8_lossy(&r.stderr)
            );
            // The index case rewrites the input's index; the others write out.
            let wrote = if c.name.starts_with("index-") {
                std::fs::metadata(format!("{}.bai", inputs.in_a.display()))
            } else {
                std::fs::metadata(&out)
            };
            assert!(
                wrote.map(|m| m.len()).unwrap_or(0) > 0,
                "case `{}`: `{}` wrote nothing",
                c.name,
                args.join(" ")
            );
            checked += 1;
        }
        assert_eq!(checked, 24, "P3a and P3c cases and convert-bedpe run");
    }
    /// The P3b cases' gloci commands (shell lines, some writing a genome file
    /// first) run on a small version of their own setup and write output.
    /// `convert-bedpe-50M` is a `sam` setup, run by the P3a test above.
    #[test]
    fn p3b_gloci_bench_commands_run_on_small_data() {
        let gloci = match std::env::var_os("GLOCI_BIN") {
            Some(b) => PathBuf::from(b),
            None => built_gloci("debug").unwrap(),
        };
        let text =
            std::fs::read_to_string(workspace_root().join("xtask/bench/cases.toml")).unwrap();
        let td = tempfile::tempdir().unwrap();
        let mut checked = 0;
        for c in parse_cases(&text).unwrap() {
            let Some(setup) = c.setup.as_deref() else {
                continue;
            };
            let kind = setup.split_whitespace().next().unwrap();
            let p3b = [
                "jaccard-",
                "fisher-",
                "reldist-",
                "sample-",
                "shuffle-",
                "random-",
                "pair-",
                "vcf-compare-",
            ];
            if !p3b.iter().any(|p| c.name.starts_with(p)) {
                continue;
            }
            let small = override_records(setup, 3000);
            // random reads no input; keep its one record.
            let small = if c.name.starts_with("random-") {
                setup.to_string()
            } else {
                small
            };
            let small = if kind == "vcf-pair" {
                // 2500 samples of 3000 records is a big file for a test.
                let mut w: Vec<String> = small.split_whitespace().map(String::from).collect();
                set_setup_value(&mut w, "--samples", "40".into());
                w.join(" ")
            } else {
                small
            };
            let case = Case {
                input: c.input.clone(),
                ..c.clone()
            };
            let inputs = prepare_input(&case, Some(&small), &gloci, td.path(), None).unwrap();
            let out = td.path().join(format!("{}.out.{}", c.name, c.input));
            let cmd = c.commands.get("gloci").expect("a gloci command");
            let path = |p: &Path| p.to_string_lossy().into_owned();
            // -n of sample and random, scaled to the small data.
            let cmd = cmd
                .replace("-n 1000000 ", "-n 1000 ")
                .replace("-n 10000000 ", "-n 1000 ");
            let line = substitute_pair(
                &cmd,
                &path(&gloci),
                &path(&inputs.in_a),
                &inputs.in_b.as_deref().map(path).unwrap_or_default(),
                &path(&out),
            );
            let r = Command::new("sh").arg("-c").arg(&line).output().unwrap();
            assert!(
                r.status.success(),
                "case `{}`: `{line}`: {}",
                c.name,
                String::from_utf8_lossy(&r.stderr)
            );
            let file = if c.name.starts_with("vcf-compare-") {
                format!("{}.diff.sites", out.display())
            } else {
                path(&out)
            };
            assert!(
                std::fs::metadata(&file).map(|m| m.len()).unwrap_or(0) > 0,
                "case `{}`: `{line}` wrote nothing to {file}",
                c.name
            );
            if c.name.starts_with("vcf-compare-") {
                for ext in ["diff.indv", "diff.discordance_matrix"] {
                    assert!(
                        std::fs::metadata(format!("{}.{ext}", out.display())).is_ok(),
                        "case `{}`: no .{ext}",
                        c.name
                    );
                }
            }
            checked += 1;
        }
        assert_eq!(checked, 9, "P3b cases run");
    }
}
