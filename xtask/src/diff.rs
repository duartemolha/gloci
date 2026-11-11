//! Differential harness: runs gloci and the pinned reference tools (bedtools,
//! bedops, GNU sort, htslib, bcftools) on the same inputs inside the
//! `gloci-diff:p0` container and compares the outputs.
//!
//! Every case in `xtask/diff/cases.toml` either must match, or names a
//! documented difference (an anchor in `docs/differences.md`) and must differ.
//! `--update-golden` stores the reference output of each single-invocation case
//! under `crates/gloci/tests/golden/<case>/`, which plain `cargo test` checks
//! without Docker (`crates/gloci/tests/golden.rs`).
//!
//! Commands run with the case directory as working directory: the input is
//! `input.<ext>` and every `aux` file keeps its basename. `gloci` is on PATH.
use anyhow::{Context, Result, anyhow, bail, ensure};
pub use gloci_testkit::differential::{Verdict, compare};
use serde::Deserialize;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub use crate::bench::DEFAULT_IMAGE;

#[derive(Deserialize, Default, Clone, Debug)]
pub struct Input {
    /// Arguments of `cargo xtask gen`: `["bed", "--records", "1500"]`.
    #[serde(rename = "gen")]
    pub generate: Option<Vec<String>>,
    /// A committed file, relative to the workspace root.
    pub file: Option<String>,
}

#[derive(Deserialize, Default, Clone, Debug)]
pub struct Case {
    pub name: String,
    #[serde(default)]
    pub input: Input,
    /// Extra files (workspace-relative) copied next to the input by basename.
    #[serde(default)]
    pub aux: Vec<String>,
    /// Shell command; stdout is compared.
    #[serde(default)]
    pub gloci: String,
    #[serde(default)]
    pub reference: String,
    #[serde(default = "exact")]
    pub compare: String,
    /// Anchor in docs/differences.md of the documented difference.
    pub expected_difference: Option<String>,
}

fn exact() -> String {
    "exact".into()
}

#[derive(Deserialize)]
struct CaseFile {
    case: Vec<Case>,
}

pub fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

pub fn load_cases(path: &Path) -> Result<Vec<Case>> {
    let text = std::fs::read_to_string(path).with_context(|| path.display().to_string())?;
    let cases = toml::from_str::<CaseFile>(&text)?.case;
    for c in &cases {
        ensure!(
            matches!(
                c.compare.as_str(),
                "exact"
                    | "lines-unordered"
                    | "vcf-records"
                    | "exit-status"
                    | "numeric"
                    | "sam"
                    | "vcf"
                    | "bam"
                    | "cram"
            ),
            "{}: unknown compare mode `{}`",
            c.name,
            c.compare
        );
        ensure!(
            c.input.generate.is_some() != c.input.file.is_some(),
            "{}: input needs exactly one of gen or file",
            c.name
        );
        ensure!(
            !c.gloci.is_empty() && !c.reference.is_empty(),
            "{}: gloci and reference are required",
            c.name
        );
    }
    Ok(cases)
}

/// A documented difference must name an anchor that exists.
pub fn validate(c: &Case, doc: &str) -> Result<(), String> {
    match &c.expected_difference {
        Some(a) => gloci_testkit::differential::check_anchor(doc, a),
        None => Ok(()),
    }
}

pub fn verdict(c: &Case, equal: bool) -> Verdict {
    gloci_testkit::differential::verdict(c.expected_difference.as_deref(), equal)
}

/// The gloci arguments when the case is a single plain `gloci ARGS...`
/// invocation (no shell syntax), which the golden test can replay.
pub fn golden_args(c: &Case) -> Option<String> {
    let rest = c.gloci.trim().strip_prefix("gloci ")?;
    if !rest.contains('\'') {
        if rest.contains(['\n', ';', '&', '|', '<', '>', '$', '"', '`', '\\']) {
            return None;
        }
        return Some(rest.split_whitespace().collect::<Vec<_>>().join(" "));
    }
    // Single-quoted arguments (filter expressions) are plain text to the shell;
    // outside the quotes no shell syntax is allowed.
    let mut quoted = false;
    for ch in rest.chars() {
        match (ch, quoted) {
            ('\'', _) => quoted = !quoted,
            ('\n', _) => return None,
            (';' | '&' | '|' | '<' | '>' | '$' | '"' | '`' | '\\', false) => return None,
            _ => {}
        }
    }
    gloci_testkit::differential::split_args(rest)?;
    Some(rest.to_string())
}

fn input_ext(c: &Case) -> String {
    match (&c.input.generate, &c.input.file) {
        (Some(g), _) => g[0].clone(),
        // everything after the first dot: `bed.gz` for `x.bed.gz`
        (_, Some(f)) => Path::new(f)
            .file_name()
            .and_then(|n| n.to_str())
            .and_then(|n| n.split_once('.'))
            .map(|(_, e)| e.to_string())
            .unwrap_or_default(),
        _ => unreachable!(),
    }
}

/// `(file name, bytes)` of the input and aux files of a case.
pub fn case_files(c: &Case, root: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let mut v = Vec::new();
    let bytes = match (&c.input.generate, &c.input.file) {
        (Some(g), _) => {
            let mut b = Vec::new();
            crate::synth::run_args(g, &mut b)?;
            b
        }
        (_, Some(f)) => std::fs::read(root.join(f)).with_context(|| format!("{}: {f}", c.name))?,
        _ => unreachable!(),
    };
    v.push((format!("input.{}", input_ext(c)), bytes));
    for a in &c.aux {
        let name = Path::new(a)
            .file_name()
            .ok_or_else(|| anyhow!("{}: bad aux {a}", c.name))?
            .to_string_lossy()
            .into_owned();
        v.push((
            name,
            std::fs::read(root.join(a)).with_context(|| format!("{}: {a}", c.name))?,
        ));
    }
    Ok(v)
}

pub struct Args {
    pub image: String,
    pub cases: PathBuf,
    pub only: Option<String>,
    pub update_golden: bool,
    /// gloci binary to mount (default: `cargo build` for the host).
    pub gloci: Option<PathBuf>,
}

/// Run `script` with `sh -c` inside the container: the case directory is
/// mounted at /work (the working directory) and gloci at /usr/local/bin/gloci.
fn docker_sh(image: &str, gloci: &Path, dir: &Path, script: &str) -> Result<Output> {
    let id = |f: &str| -> Result<String> {
        let o = Command::new("id").arg(f).output()?;
        Ok(String::from_utf8(o.stdout)?.trim().to_string())
    };
    Command::new("docker")
        .args(["run", "--rm", "--network", "none", "--user"])
        .arg(format!("{}:{}", id("-u")?, id("-g")?))
        .args(["-e", "HOME=/tmp", "-w", "/work", "-v"])
        .arg(format!("{}:/work", dir.display()))
        .arg("-v")
        .arg(format!("{}:/usr/local/bin/gloci:ro", gloci.display()))
        .arg(image)
        .args(["sh", "-c", script])
        .output()
        .context("running docker")
}

fn stderr_of(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).trim().to_string()
}

/// The bytes to store for a golden output: `old`, the file already there, when it compares
/// equal to `new` under the case's compare mode, so that a regeneration does not churn files
/// whose only change is in lines the mode ignores (bcftools' date-stamped
/// `##bcftools_*Command` headers); else `new`.
pub fn kept_golden<'a>(compare_mode: &str, old: Option<&'a [u8]>, new: &'a [u8]) -> &'a [u8] {
    match old {
        Some(o) if compare(compare_mode, o, new) => o,
        _ => new,
    }
}

fn write_golden(
    c: &Case,
    root: &Path,
    args: &str,
    expected: &[u8],
    expected_gloci: Option<&[u8]>,
) -> Result<()> {
    let dir = root.join("crates/gloci/tests/golden").join(&c.name);
    let old_expected = std::fs::read(dir.join("expected")).ok();
    let old_gloci = std::fs::read(dir.join("expected_gloci")).ok();
    let expected = kept_golden(&c.compare, old_expected.as_deref(), expected);
    let expected_gloci = expected_gloci.map(|g| kept_golden(&c.compare, old_gloci.as_deref(), g));
    if dir.exists() {
        std::fs::remove_dir_all(&dir)?;
    }
    std::fs::create_dir_all(&dir)?;
    for (n, b) in case_files(c, root)? {
        std::fs::write(dir.join(n), b)?;
    }
    std::fs::write(dir.join("args"), format!("{args}\n"))?;
    std::fs::write(dir.join("compare"), format!("{}\n", c.compare))?;
    if let Some(a) = &c.expected_difference {
        std::fs::write(dir.join("expected_difference"), format!("{a}\n"))?;
    }
    std::fs::write(dir.join("expected"), expected)?;
    // gloci's own current output, so a regression to a different wrong output
    // cannot pass as "still differs".
    if let Some(g) = expected_gloci {
        std::fs::write(dir.join("expected_gloci"), g)?;
    }
    Ok(())
}

/// Returns true when every case passed.
pub fn run(a: &Args) -> Result<bool> {
    let root = workspace_root();
    let gloci = match &a.gloci {
        Some(p) => std::fs::canonicalize(p)
            .with_context(|| format!("gloci binary {} not found", p.display()))?,
        None => {
            let st = Command::new("cargo")
                .args(["build", "-p", "gloci"])
                .current_dir(&root)
                .status()?;
            ensure!(st.success(), "cargo build failed");
            root.join("target/debug/gloci")
        }
    };
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
    let doc = std::fs::read_to_string(root.join("docs/differences.md"))?;
    let cases = load_cases(&a.cases)?;
    let work = root.join("target/diff-work");
    let mut all_ok = true;
    let mut ran = 0;
    for c in cases
        .iter()
        .filter(|c| a.only.as_ref().is_none_or(|o| *o == c.name))
    {
        if let Err(e) = validate(c, &doc) {
            eprintln!("FAIL  {}: {e}", c.name);
            all_ok = false;
            continue;
        }
        let golden = golden_args(c);
        if a.update_golden && golden.is_none() {
            eprintln!(
                "skip  {} (not a single gloci invocation; live only)",
                c.name
            );
            continue;
        }
        ran += 1;
        let dir = work.join(&c.name);
        if dir.exists() {
            std::fs::remove_dir_all(&dir)?;
        }
        std::fs::create_dir_all(&dir)?;
        for (n, b) in case_files(c, &root)? {
            std::fs::write(dir.join(n), b)?;
        }
        let outcome = c.compare == "exit-status";
        let observe = |o: &Output| {
            gloci_testkit::differential::observe(&c.compare, o.status.code(), &o.stdout)
        };
        let reference = docker_sh(&a.image, &gloci, &dir, &c.reference)?;
        if !outcome && !reference.status.success() {
            bail!(
                "{}: reference command failed: {}",
                c.name,
                stderr_of(&reference)
            );
        }
        let reference_out = observe(&reference);
        if a.update_golden {
            let mine_out = if c.expected_difference.is_some() {
                let mine = docker_sh(&a.image, &gloci, &dir, &c.gloci)?;
                if !outcome && !mine.status.success() {
                    bail!("{}: gloci failed: {}", c.name, stderr_of(&mine));
                }
                Some(observe(&mine))
            } else {
                None
            };
            write_golden(
                c,
                &root,
                golden.as_deref().unwrap(),
                &reference_out,
                mine_out.as_deref(),
            )?;
            eprintln!("wrote {}", c.name);
            continue;
        }
        let mine = docker_sh(&a.image, &gloci, &dir, &c.gloci)?;
        if !outcome && !mine.status.success() {
            eprintln!("FAIL  {}: gloci failed: {}", c.name, stderr_of(&mine));
            all_ok = false;
            continue;
        }
        let mine_out = observe(&mine);
        let equal = compare(&c.compare, &mine_out, &reference_out);
        let v = verdict(c, equal);
        let word = match v {
            Verdict::Pass if c.expected_difference.is_some() => "ok    (differs as documented)",
            Verdict::Pass => "ok    (matches)",
            Verdict::Mismatch => "FAIL  mismatch",
            Verdict::StaleDifference => "FAIL  stale expected_difference: outputs are now equal",
        };
        eprintln!("{word}  {}", c.name);
        if v != Verdict::Pass {
            all_ok = false;
            let _ = std::fs::write(dir.join("gloci.out"), &mine_out);
            let _ = std::fs::write(dir.join("reference.out"), &reference_out);
            eprintln!("      outputs in {}", dir.display());
        }
        if v == Verdict::Pass {
            let _ = std::fs::remove_dir_all(&dir);
        }
    }
    ensure!(ran > 0, "no case matched");
    Ok(all_ok)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(expected: Option<&str>) -> Case {
        Case {
            name: "c".into(),
            expected_difference: expected.map(String::from),
            ..Case::default()
        }
    }

    #[test]
    fn compare_modes() {
        assert!(compare("exact", b"a\nb\n", b"a\nb\n"));
        assert!(
            compare("lines-unordered", b"a\nb\n", b"b\na\n")
                && !compare("exact", b"a\nb\n", b"b\na\n")
        );
        assert!(compare(
            "vcf-records",
            b"##x=1\n#CHROM\nr1\n",
            b"##y=2\n#CHROM\nr1\n"
        ));
        assert!(!compare("vcf-records", b"#h\nr1\nr2\n", b"#h\nr2\nr1\n"));
        assert!(compare(
            "vcf",
            b"##a=<ID=x>\n##bcftools_viewVersion=1\n##bcftools_viewCommand=v\n#CHROM\nr1\n",
            b"##a=<ID=x>\n##glociCommand=g\n#CHROM\nr1\n"
        ));
        assert!(!compare(
            "vcf",
            b"##a=<ID=x>\n#CHROM\nr1\n",
            b"##a=<ID=x,IDX=0>\n#CHROM\nr1\n"
        ));
        assert!(compare("numeric", b"a\t2.5\n", b"a\t2.5000001\n"));
        assert!(!compare("numeric", b"a\t2.5\n", b"a\t2.6\n"));
    }

    #[test]
    fn regenerating_keeps_a_golden_equal_under_its_compare_mode() {
        let old = b"##bcftools_isecCommand=isec; Date=Mon Nov 10\n#CHROM\nr1\n".as_slice();
        let new = b"##bcftools_isecCommand=isec; Date=Tue Nov 11\n#CHROM\nr1\n".as_slice();
        assert_eq!(kept_golden("vcf-records", Some(old), new), old);
        // a real change, or a mode that sees the header, rewrites it
        let changed = b"##x\n#CHROM\nr2\n".as_slice();
        assert_eq!(kept_golden("vcf-records", Some(old), changed), changed);
        assert_eq!(kept_golden("exact", Some(old), new), new);
        assert_eq!(kept_golden("exact", None, new), new);
    }

    #[test]
    fn expected_difference_rules() {
        let doc = "## View region footprints {#view-region-footprints}\n";
        assert!(validate(&case(Some("view-region-footprints")), doc).is_ok());
        assert!(validate(&case(Some("#view-region-footprints")), doc).is_ok());
        assert!(validate(&case(Some("nope")), doc).is_err());
        assert!(validate(&case(None), doc).is_ok());
        assert_eq!(
            verdict(&case(Some("view-region-footprints")), true),
            Verdict::StaleDifference
        );
        assert_eq!(
            verdict(&case(Some("view-region-footprints")), false),
            Verdict::Pass
        );
        assert_eq!(verdict(&case(None), false), Verdict::Mismatch);
        assert_eq!(verdict(&case(None), true), Verdict::Pass);
    }

    #[test]
    fn manifest_parses_and_all_anchors_exist() {
        let root = workspace_root();
        let cases = load_cases(&root.join("xtask/diff/cases.toml")).unwrap();
        let doc = std::fs::read_to_string(root.join("docs/differences.md")).unwrap();
        let names: Vec<_> = cases.iter().map(|c| c.name.as_str()).collect();
        for want in [
            "sort-bed-vs-sort-bed",
            "sort-bed-vs-bedtools-sort",
            "sort-bed-vs-gnu-sort",
            "sort-bed-karyotype-default-vs-bedtools-sort",
            "sort-bed-vs-bedtools-sort-genome",
            "sort-vcf-vs-bcftools-sort",
            "view-region-ref-inclusive-vs-tabix",
            "view-region-default-vs-tabix",
            "index-interop-tabix-queries",
        ] {
            assert!(names.contains(&want), "missing case {want}");
        }
        let mut seen = std::collections::HashSet::new();
        for c in &cases {
            assert!(seen.insert(&c.name), "duplicate case {}", c.name);
            validate(c, &doc).unwrap_or_else(|e| panic!("{}: {e}", c.name));
        }
    }

    #[test]
    fn every_documented_difference_is_pinned_by_a_case() {
        let root = workspace_root();
        let cases = load_cases(&root.join("xtask/diff/cases.toml")).unwrap();
        let doc = std::fs::read_to_string(root.join("docs/differences.md")).unwrap();
        for a in gloci_testkit::differential::anchors(&doc) {
            assert!(
                cases.iter().any(|c| c
                    .expected_difference
                    .as_deref()
                    .is_some_and(|d| d.trim_start_matches('#') == a)),
                "docs/differences.md section {{#{a}}} is not pinned by any case"
            );
        }
    }

    #[test]
    fn golden_dirs_match_the_manifest() {
        let root = workspace_root();
        let cases = load_cases(&root.join("xtask/diff/cases.toml")).unwrap();
        let golden = root.join("crates/gloci/tests/golden");
        let mut eligible = 0;
        for c in &cases {
            let Some(args) = golden_args(c) else { continue };
            eligible += 1;
            let dir = golden.join(&c.name);
            let read = |f: &str| {
                std::fs::read(dir.join(f)).unwrap_or_else(|e| {
                    panic!(
                        "{}: {f}: {e} (run cargo xtask diff --update-golden)",
                        c.name
                    )
                })
            };
            assert_eq!(read("args"), format!("{args}\n").into_bytes(), "{}", c.name);
            assert_eq!(read("compare"), format!("{}\n", c.compare).into_bytes());
            for (name, bytes) in case_files(c, &root).unwrap() {
                assert_eq!(read(&name), bytes, "{}: {name} is stale", c.name);
            }
            // an empty reference output is meaningful only for a documented
            // difference (the reference tool misses what gloci finds)
            assert!(
                c.expected_difference.is_some() || !read("expected").is_empty(),
                "{}: empty expected output makes the case vacuous",
                c.name
            );
            assert_eq!(
                dir.join("expected_gloci").exists(),
                c.expected_difference.is_some(),
                "{}: expected_gloci must exist exactly for documented differences",
                c.name
            );
        }
        assert!(eligible >= 8, "only {eligible} golden cases");
        // no orphan directories
        for e in std::fs::read_dir(&golden).unwrap() {
            let n = e.unwrap().file_name().to_string_lossy().into_owned();
            assert!(
                cases
                    .iter()
                    .any(|c| c.name == n && golden_args(c).is_some()),
                "orphan golden dir {n}"
            );
        }
    }
}
