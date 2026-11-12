//! `cargo xtask paper run --bench-only`: the paper's real-data benchmarks
//! (D96), run by the ordinary bench runner with peak RSS in every run.
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail, ensure};
use serde::Deserialize;

use crate::bench::{Case, parse_cases};
use crate::diff::workspace_root;

/// Where the paper's cases are, relative to the repository.
pub const BENCH_PATH: &str = "xtask/paper/bench.toml";

#[derive(Deserialize)]
struct PaperBenchFile {
    /// Names of cases of `xtask/bench/cases.toml` that the paper reports too.
    #[serde(default)]
    from_nightly: Vec<String>,
    #[serde(default)]
    case: Vec<toml::Value>,
}

/// The paper's cases: the `from_nightly` ones, cloned from `nightly` by name,
/// then those of the file itself.
pub fn parse_paper_cases(text: &str, nightly: &[Case]) -> Result<Vec<Case>> {
    let f: PaperBenchFile = toml::from_str(text).context("parsing paper bench.toml")?;
    let mut out = Vec::new();
    for name in &f.from_nightly {
        let c = nightly
            .iter()
            .find(|c| &c.name == name)
            .with_context(|| format!("from_nightly: no case `{name}` in cases.toml"))?;
        out.push(c.clone());
    }
    let own = toml::to_string(&toml::Value::Table(
        [("case".to_string(), toml::Value::Array(f.case))]
            .into_iter()
            .collect(),
    ))?;
    out.extend(parse_cases(&own)?);
    let mut seen = std::collections::HashSet::new();
    for c in &out {
        ensure!(seen.insert(c.name.clone()), "duplicate case `{}`", c.name);
    }
    Ok(out)
}

pub fn load_paper_cases() -> Result<Vec<Case>> {
    let root = workspace_root();
    let nightly = parse_cases(&std::fs::read_to_string(
        root.join("xtask/bench/cases.toml"),
    )?)?;
    parse_paper_cases(&std::fs::read_to_string(root.join(BENCH_PATH))?, &nightly)
}

pub struct PaperBenchArgs {
    pub only: Option<String>,
    pub gloci: Option<PathBuf>,
    pub no_docker: bool,
    pub image: String,
    pub data_dir: PathBuf,
    pub out: PathBuf,
    /// Timed runs per command instead of each case's own (smoke runs).
    pub runs: Option<u32>,
    /// A Python wheel of the same version: its name and sha256 are recorded
    /// in environment.json (the Python benchmarks come later).
    pub wheel: Option<PathBuf>,
}

pub fn run(a: &PaperBenchArgs) -> Result<bool> {
    let mut cases = load_paper_cases()?;
    if let Some(w) = &a.wheel {
        crate::paper::environment::wheel_info(w)?;
    }
    if let Some(n) = a.runs {
        ensure!(n >= 1, "--runs must be at least 1");
        for c in &mut cases {
            c.runs = n;
        }
    }
    if let Some(o) = &a.only
        && !cases.iter().any(|c| &c.name == o)
    {
        bail!("--only {o}: not a paper benchmark (see {BENCH_PATH})");
    }
    let version_dir: PathBuf = Path::new(&a.out).join(env!("CARGO_PKG_VERSION"));
    let results = version_dir.join("bench");
    std::fs::create_dir_all(&results)?;
    let started_at = crate::paper::data::now_utc();
    let ok = crate::bench::run_cases(
        &crate::bench::Args {
            cases: workspace_root().join(BENCH_PATH),
            out: a.out.clone(),
            check_targets: false,
            gloci: a.gloci.clone(),
            no_docker: a.no_docker,
            image: a.image.clone(),
            only: a.only.clone(),
            records: None,
            data_dir: a.data_dir.clone(),
            results: Some(results),
        },
        cases.clone(),
    )?;
    write_run_records(a, &version_dir, &cases, &started_at)?;
    Ok(ok)
}

/// `environment.json` and `repo-snapshot.json` next to the results.
pub fn write_run_records(
    a: &PaperBenchArgs,
    dir: &Path,
    cases: &[Case],
    started_at: &str,
) -> Result<()> {
    use crate::paper::{environment, snapshot};
    let root = workspace_root();
    let commit = environment::git_commit(&root);
    let version = env!("CARGO_PKG_VERSION");
    let gloci = a.gloci.clone().or_else(|| {
        let p = crate::bench::target_dir().join("release/gloci");
        p.is_file().then_some(p)
    });
    let env = environment::environment(&environment::EnvArgs {
        version,
        commit: commit.clone(),
        started_at,
        finished_at: &crate::paper::data::now_utc(),
        image: (!a.no_docker).then_some(a.image.as_str()),
        gloci: gloci.as_deref(),
        wheel: a.wheel.as_deref(),
        cases,
        extra_tools: &[],
    })?;
    std::fs::write(
        dir.join("environment.json"),
        serde_json::to_string_pretty(&env)? + "\n",
    )?;
    let cat = crate::paper::fetch::load_catalog()?;
    let snap = snapshot::snapshot(&root, &cat, &a.data_dir, version, commit.as_deref())?;
    std::fs::write(
        dir.join("repo-snapshot.json"),
        serde_json::to_string_pretty(&snap)? + "\n",
    )?;
    Ok(())
}
