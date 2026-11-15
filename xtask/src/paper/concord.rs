//! `cargo xtask paper run --concordance-only`: gloci's output on the paper
//! datasets against the reference tools', record by record (D97).
//!
//! Each case of `xtask/paper/concordance.toml` runs a reference command once
//! and gloci once or twice in the bench image, with the dataset cache mounted
//! read-only: the default run (with `--chrom-order lex` added to every gloci
//! subcommand that takes it, as the derived inputs are sorted with
//! `LC_ALL=C`) and, when the case has `compat` flags, the compat run. Each run
//! is compared with the reference output by
//! [`gloci_testkit::differential::match_stats`] and gets a [`Verdict`].
//! `<out>/<version>/concordance.csv` gets one row per run, and
//! `<out>/<version>/concordance/<case>.<mode>.txt` the commands, the counts
//! and the first 10 unmatched lines of each side.
//!
//! Reference outputs are kept in `target/concordance-work` and reused while
//! their key is the same: the image Id (or the host tool's `--version`), the
//! pinned digest of every input and the command ([`reference_stamp`]), so a
//! re-run after a gloci fix does not wait for vcftools again. Delete the
//! directory to run them afresh.
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{Context, Result, bail, ensure};
use gloci_testkit::differential::{
    CONCORDANCE_MODES, MatchStats, Residual, line_key, match_stats_with_residual,
};
use serde::Deserialize;

use crate::bench::{Catalog, bench_cache_dir, input_path, resolve_gloci, target_dir};
use crate::diff::workspace_root;

/// Where the cases are, relative to the repository.
pub const CASES_PATH: &str = "xtask/paper/concordance.toml";

/// The columns of `concordance.csv`.
pub const CSV_HEADER: &str = "case,command,reference_tool,dataset,mode,n_gloci,n_reference,n_matching,match_rate,verdict,anchors";

/// Commands of docs/migration.md without a concordance case, and why.
pub const NO_CONCORDANCE_CASE: &[(&str, &str)] = &[
    (
        "index",
        "writes a binary index, not records; covered by diff cases only (tabix, samtools and bcftools interop)",
    ),
    (
        "vcf norm",
        "needs a reference FASTA, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "vcf merge",
        "needs a second call set with other samples, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "vcf compare",
        "needs two call sets of the same samples; the 1kg halves are disjoint regions; covered by diff cases only",
    ),
    (
        "getfasta",
        "needs a reference FASTA, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "maskfasta",
        "needs a reference FASTA, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "nuc",
        "needs a reference FASTA, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "faidx",
        "needs a reference FASTA, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "dict",
        "needs a reference FASTA, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "pair-intersect",
        "needs BEDPE pairs, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "pair-pair",
        "needs BEDPE pairs, which no paper dataset provides; covered by diff cases only",
    ),
    (
        "shuffle",
        "random output from another generator than bedtools' (differences.md#rng-stream): no record-level reference",
    ),
    (
        "random",
        "random output from another generator than bedtools' (differences.md#rng-stream): no record-level reference",
    ),
    (
        "sample",
        "random output from another generator than bedtools' (differences.md#rng-stream): no record-level reference",
    ),
    (
        "flags",
        "converts FLAG numbers given on the command line; there is no dataset to run it on; covered by diff cases only",
    ),
    (
        "check",
        "reports pass or fail per file (exit status); covered by diff cases only",
    ),
    (
        "samples",
        "lists header read groups; the ENCODE BAM has none; covered by diff cases only",
    ),
    (
        "reheader",
        "rewrites a header without touching records; covered by diff cases only",
    ),
    (
        "igv",
        "writes a script for IGV from the feature list, with no computation to compare at scale; covered by diff cases only",
    ),
    (
        "links",
        "writes a page of browser links from the feature list, with no computation to compare at scale; covered by diff cases only",
    ),
    (
        "cat",
        "needs several BAM files of one reference; the paper has one; covered by diff cases only",
    ),
];

#[derive(Deserialize, Clone, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct ConcordCase {
    pub name: String,
    /// The gloci command as docs/migration.md maps it, such as `vcf filter`.
    pub command: String,
    pub reference_tool: String,
    /// `a`, `b`, `c` and `genome`: datasets or derived inputs of
    /// `datasets.toml`, bound to `{in_a}`, `{in_b}`, `{in_c}` and `{genome}`.
    pub inputs: BTreeMap<String, String>,
    /// Shell commands that write their result to `{out}`; `{gloci}` is the
    /// binary under test.
    pub gloci: String,
    pub reference: String,
    #[serde(default = "exact")]
    pub compare: String,
    /// Flags the compat run adds after the first gloci subcommand.
    #[serde(default)]
    pub compat: Option<String>,
    /// The documented difference the default run may show.
    #[serde(default)]
    pub expected_difference: Option<String>,
    /// The documented difference the compat run may show.
    #[serde(default)]
    pub compat_difference: Option<String>,
    /// `numeric` only: the reference prints this many fixed decimals (as
    /// bedtools coverage prints 7 and bedmap 6), which widens the tolerance
    /// of its numbers to half a unit of the last one.
    #[serde(default)]
    pub ref_decimals: Option<u32>,
    /// The lowest rate the default run's documented difference explains:
    /// the rate of the run that was checked line by line, truncated to six
    /// decimals. Below it the run is a mismatch.
    #[serde(default)]
    pub min_rate: Option<f64>,
    /// The same for the compat run.
    #[serde(default)]
    pub compat_min_rate: Option<f64>,
    /// The default run's documented difference is the order alone: the run
    /// passes only with every line matched (rate 1) in another order.
    #[serde(default)]
    pub order_only: bool,
    /// The same for the compat run.
    #[serde(default)]
    pub compat_order_only: bool,
    /// Why an empty reference output is right for this case. Without it, a
    /// reference that writes nothing fails the run (a tool that wrote
    /// nothing must not count as a match).
    #[serde(default)]
    pub allow_empty: Option<String>,
}

/// What a run of a case may show: its documented difference, with the
/// lowest rate that difference explains or that it is about order alone.
#[derive(Debug, Clone, Copy, Default)]
pub struct Expect<'a> {
    pub anchor: Option<&'a str>,
    pub min_rate: Option<f64>,
    pub order_only: bool,
    /// An empty reference output is allowed.
    pub allow_empty: bool,
}

impl ConcordCase {
    /// What the default (`compat` false) or compat run may show.
    pub fn expect(&self, compat: bool) -> Expect<'_> {
        if compat {
            Expect {
                anchor: self.compat_difference.as_deref(),
                min_rate: self.compat_min_rate,
                order_only: self.compat_order_only,
                allow_empty: self.allow_empty.is_some(),
            }
        } else {
            Expect {
                anchor: self.expected_difference.as_deref(),
                min_rate: self.min_rate,
                order_only: self.order_only,
                allow_empty: self.allow_empty.is_some(),
            }
        }
    }
}

fn exact() -> String {
    "exact".into()
}

#[derive(Deserialize)]
struct CaseFile {
    #[serde(default)]
    case: Vec<ConcordCase>,
}

const SLOTS: [(&str, &str); 4] = [
    ("a", "{in_a}"),
    ("b", "{in_b}"),
    ("c", "{in_c}"),
    ("genome", "{genome}"),
];

/// The clap command of a subcommand path such as `vcf filter`.
fn subcommand(path: &str) -> Option<clap::Command> {
    let mut cmd = gloci::command();
    for w in path.split_whitespace() {
        cmd = cmd.find_subcommand(w)?.clone();
    }
    (cmd.get_name() != "gloci").then_some(cmd)
}

fn has_subcommands(cmd: &clap::Command) -> bool {
    cmd.get_subcommands().any(|s| s.get_name() != "help")
}

pub fn parse_cases(text: &str) -> Result<Vec<ConcordCase>> {
    let cases = toml::from_str::<CaseFile>(text)
        .context("parsing concordance.toml")?
        .case;
    let mut seen = BTreeSet::new();
    for c in &cases {
        let n = &c.name;
        ensure!(!n.is_empty(), "a case has no name");
        ensure!(seen.insert(n.clone()), "duplicate case `{n}`");
        ensure!(
            subcommand(&c.command).is_some_and(|s| !has_subcommands(&s)),
            "{n}: command `{}` is not a gloci subcommand",
            c.command
        );
        ensure!(
            CONCORDANCE_MODES.contains(&c.compare.as_str()),
            "{n}: compare `{}` is not one of {CONCORDANCE_MODES:?}",
            c.compare
        );
        ensure!(!c.inputs.is_empty(), "{n}: no inputs");
        for slot in c.inputs.keys() {
            ensure!(
                SLOTS.iter().any(|(s, _)| s == slot),
                "{n}: unknown input slot `{slot}` (use a, b, c or genome)"
            );
        }
        ensure!(c.gloci.contains("{gloci}"), "{n}: gloci has no {{gloci}}");
        for (what, cmd) in [("gloci", &c.gloci), ("reference", &c.reference)] {
            ensure!(
                cmd.contains("{out}"),
                "{n}: {what} must write its result to {{out}}"
            );
            for (slot, ph) in SLOTS {
                ensure!(
                    !cmd.contains(ph) || c.inputs.contains_key(slot),
                    "{n}: {what} uses {ph} but inputs has no `{slot}`"
                );
            }
        }
        ensure!(
            c.compat_difference.is_none() || c.compat.is_some(),
            "{n}: compat_difference without compat flags"
        );
        ensure!(
            c.compat.as_deref().is_none_or(|f| !f.trim().is_empty()),
            "{n}: empty compat flags"
        );
        ensure!(
            c.ref_decimals.is_none() || c.compare == "numeric",
            "{n}: ref_decimals applies to compare = \"numeric\" only"
        );
        for (what, anchor, min_rate, order_only) in [
            ("", &c.expected_difference, c.min_rate, c.order_only),
            (
                "compat_",
                &c.compat_difference,
                c.compat_min_rate,
                c.compat_order_only,
            ),
        ] {
            ensure!(
                anchor.is_some() || (min_rate.is_none() && !order_only),
                "{n}: {what}min_rate and {what}order_only need {what}difference"
            );
            ensure!(
                !(order_only && min_rate.is_some()),
                "{n}: {what}order_only runs match every line, so no {what}min_rate"
            );
            ensure!(
                min_rate.is_none_or(|r| (0.0..1.0).contains(&r)),
                "{n}: {what}min_rate must be at least 0 and below 1"
            );
        }
        ensure!(
            c.allow_empty
                .as_deref()
                .is_none_or(|r| !r.trim().is_empty()),
            "{n}: allow_empty needs a reason"
        );
        for (what, cmd) in [("gloci", &c.gloci), ("reference", &c.reference)] {
            ensure!(
                !cmd.contains("<(") && !cmd.contains(">("),
                "{n}: {what} uses a process substitution, whose failure the run cannot see: \
                 write the intermediate to a file ({{out}}.x) in the same command instead"
            );
        }
    }
    Ok(cases)
}

pub fn load_cases() -> Result<Vec<ConcordCase>> {
    let p = workspace_root().join(CASES_PATH);
    parse_cases(&std::fs::read_to_string(&p).with_context(|| p.display().to_string())?)
}

/// The words of the subcommand path at the start of `rest` (`sort`, or
/// `vcf filter`), and the bytes they take, leading spaces included.
fn subcommand_path(rest: &str) -> (usize, String) {
    let mut at = 0;
    let mut words: Vec<&str> = Vec::new();
    loop {
        let skip = rest[at..].len() - rest[at..].trim_start().len();
        let start = at + skip;
        let end = rest[start..]
            .find(char::is_whitespace)
            .map_or(rest.len(), |e| start + e);
        let w = &rest[start..end];
        if w.is_empty() {
            break;
        }
        let mut path = words.clone();
        path.push(w);
        let Some(cmd) = subcommand(&path.join(" ")) else {
            break;
        };
        words.push(w);
        at = end;
        if !has_subcommands(&cmd) {
            break;
        }
    }
    (at, words.join(" "))
}

fn takes_chrom_order(path: &str) -> bool {
    subcommand(path).is_some_and(|c| {
        c.get_arguments()
            .any(|a| a.get_long() == Some("chrom-order"))
    })
}

/// The gloci command template of a run: `--chrom-order lex` after every
/// gloci subcommand that takes it and does not set it, and, in the compat
/// run, the case's compat flags after the first subcommand.
pub fn gloci_command(c: &ConcordCase, compat: bool) -> String {
    const G: &str = "{gloci}";
    let mut out = String::new();
    let mut rest = c.gloci.as_str();
    let mut first = true;
    while let Some(i) = rest.find(G) {
        out.push_str(&rest[..i + G.len()]);
        rest = &rest[i + G.len()..];
        let (len, path) = subcommand_path(rest);
        out.push_str(&rest[..len]);
        rest = &rest[len..];
        let segment = &rest[..rest.find(G).unwrap_or(rest.len())];
        if takes_chrom_order(&path) && !segment.contains("--chrom-order") {
            out.push_str(" --chrom-order lex");
        }
        if compat
            && first
            && let Some(f) = &c.compat
        {
            out.push(' ');
            out.push_str(f.trim());
        }
        first = false;
    }
    out.push_str(rest);
    out
}

/// What the placeholders of a command stand for.
pub struct Slots {
    pub gloci: String,
    pub in_a: String,
    pub in_b: String,
    pub in_c: String,
    pub genome: String,
    pub out: String,
}

pub fn fill(cmd: &str, s: &Slots) -> String {
    cmd.replace("{gloci}", &s.gloci)
        .replace("{in_a}", &s.in_a)
        .replace("{in_b}", &s.in_b)
        .replace("{in_c}", &s.in_c)
        .replace("{genome}", &s.genome)
        .replace("{out}", &s.out)
}

/// The argument lists of every `gloci_word` invocation in a shell command:
/// words split as the shell would (single and double quotes kept together),
/// each list ending at an unquoted `|`, `&`, `;`, `<` or `>` (a file
/// descriptor number before a redirection is dropped).
pub fn gloci_invocations(cmd: &str, gloci_word: &str) -> Vec<Vec<String>> {
    let mut segments: Vec<Vec<String>> = vec![Vec::new()];
    let mut word = String::new();
    let (mut have, mut quote) = (false, None::<char>);
    for ch in cmd.chars() {
        match (ch, quote) {
            ('\'' | '"', None) => {
                quote = Some(ch);
                have = true;
            }
            (c, Some(q)) if c == q => quote = None,
            (c, Some(_)) => word.push(c),
            (c, None) if c.is_whitespace() => {
                if have {
                    segments.last_mut().unwrap().push(std::mem::take(&mut word));
                    have = false;
                }
            }
            ('|' | '&' | ';' | '<' | '>', None) => {
                let fd = matches!(ch, '<' | '>')
                    && !word.is_empty()
                    && word.chars().all(|d| d.is_ascii_digit());
                if have && !fd {
                    segments.last_mut().unwrap().push(std::mem::take(&mut word));
                }
                word.clear();
                have = false;
                segments.push(Vec::new());
            }
            (c, None) => {
                word.push(c);
                have = true;
            }
        }
    }
    if have {
        segments.last_mut().unwrap().push(word);
    }
    segments
        .into_iter()
        .filter(|s| s.first().is_some_and(|w| w == gloci_word))
        .collect()
}

/// The gloci commands the mapping tables of docs/migration.md give for an
/// available tool option (rows that say "not yet available" or "not planned"
/// are skipped), as subcommand paths such as `vcf filter`.
pub fn migration_commands(md: &str) -> Vec<String> {
    let tables = md.split_once("## Command mapping\n").map_or("", |(_, t)| t);
    let tables = tables
        .split_once("\n## Other behaviour changes")
        .map_or(tables, |(t, _)| t);
    let mut out = BTreeSet::new();
    for line in tables.lines().filter(|l| l.starts_with('|')) {
        let cells = split_cells(line);
        let Some(right) = cells.get(2).map(|c| c.trim()) else {
            continue;
        };
        if right.starts_with("not yet available") || right.starts_with("not planned") {
            continue;
        }
        for (i, code) in right.split('`').enumerate() {
            if i % 2 == 0 {
                continue;
            }
            let mut rest = code;
            while let Some(at) = rest.find("gloci ") {
                let before = rest[..at].chars().next_back();
                rest = &rest[at + "gloci ".len()..];
                if before.is_some_and(|c| c.is_alphanumeric() || c == '-') {
                    continue;
                }
                let (len, path) = subcommand_path(rest);
                if len > 0 && subcommand(&path).is_some_and(|c| !has_subcommands(&c)) {
                    out.insert(path);
                }
            }
        }
    }
    out.into_iter().collect()
}

/// The cells of a Markdown table row, split at unescaped `|`.
fn split_cells(line: &str) -> Vec<String> {
    let mut cells = vec![String::new()];
    let mut escaped = false;
    for ch in line.chars() {
        match ch {
            '|' if !escaped => cells.push(String::new()),
            c => {
                escaped = c == '\\';
                cells.last_mut().unwrap().push(c);
            }
        }
        if ch == '|' {
            escaped = false;
        }
    }
    cells
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// The same records in the same order.
    Match,
    /// A match where the case names a documented difference: on this data
    /// the difference never fires.
    AnchorUnused,
    /// The outputs differ and the case names the documented difference
    /// (whose sample of unmatched lines must show that difference).
    Documented,
    /// The same records in another order, with no documented difference.
    Reordered,
    /// Different records, with no documented difference.
    Mismatch,
    /// A command failed or wrote no output.
    Failed,
}

impl Verdict {
    pub fn as_str(&self) -> &'static str {
        match self {
            Verdict::Match => "match",
            Verdict::AnchorUnused => "anchor-unused",
            Verdict::Documented => "documented",
            Verdict::Reordered => "reordered",
            Verdict::Mismatch => "mismatch",
            Verdict::Failed => "failed",
        }
    }

    pub fn passes(&self) -> bool {
        matches!(
            self,
            Verdict::Match | Verdict::AnchorUnused | Verdict::Documented
        )
    }
}

/// The verdict of a run: `stats` is `None` when a command failed; `expect`
/// is the documented difference the run may show.
///
/// - A failed command, or an empty reference output where the case does not
///   allow one: `Failed`.
/// - The same lines in the same order: `Match`, or `AnchorUnused` when a
///   difference is documented.
/// - Different with nothing documented: `Reordered` (every line matched) or
///   `Mismatch`.
/// - An order-only difference documents a run whose lines all match in
///   another order; anything else is `Mismatch`.
/// - Any other documented difference documents a run whose rate is at least
///   its `min_rate`; below it, or without one, `Mismatch`. A run whose lines
///   all match in another order is `Reordered`: its difference is not about
///   order.
pub fn verdict(stats: Option<&MatchStats>, expect: &Expect) -> Verdict {
    let Some(s) = stats else {
        return Verdict::Failed;
    };
    if s.n_reference == 0 && !expect.allow_empty {
        return Verdict::Failed;
    }
    let all = s.rate() == 1.0;
    match (all && s.same_order, expect.anchor.is_some()) {
        (true, false) => Verdict::Match,
        (true, true) => Verdict::AnchorUnused,
        (false, false) if all => Verdict::Reordered,
        (false, false) => Verdict::Mismatch,
        (false, true) if expect.order_only => {
            if all {
                Verdict::Documented
            } else {
                Verdict::Mismatch
            }
        }
        (false, true) if all => Verdict::Reordered,
        (false, true) => match expect.min_rate {
            Some(m) if s.rate() >= m => Verdict::Documented,
            _ => Verdict::Mismatch,
        },
    }
}

/// The match rate with six decimals, truncated from the counts (exact
/// integer arithmetic), so a run that is not a full match never reads
/// `1.000000`.
pub fn rate_text(s: &MatchStats) -> String {
    let n = s.n_gloci.max(s.n_reference);
    if n == 0 {
        return "1.000000".to_string();
    }
    let micro = u128::from(s.n_matching) * 1_000_000 / u128::from(n);
    format!("{}.{:06}", micro / 1_000_000, micro % 1_000_000)
}

fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// The CSV text with the rows of the cases in `ran` replaced by `new_rows`:
/// rows of other cases of `order` are kept, rows of cases no longer in
/// `order` dropped, and rows sorted by `order`, then default before compat.
pub fn merge_csv(
    old: Option<&str>,
    new_rows: &[String],
    ran: &[String],
    order: &[String],
) -> String {
    let case_of = |row: &str| row.split(',').next().unwrap_or("").to_string();
    let mode_rank = |row: &str| match row.split(',').nth(4) {
        Some("default") => 0,
        _ => 1,
    };
    let mut rows: Vec<String> = old
        .unwrap_or("")
        .lines()
        .skip(1)
        .filter(|r| !r.is_empty())
        .filter(|r| {
            let c = case_of(r);
            !ran.contains(&c) && order.contains(&c)
        })
        .map(String::from)
        .collect();
    rows.extend(new_rows.iter().cloned());
    rows.sort_by_key(|r| {
        let c = case_of(r);
        (order.iter().position(|o| *o == c), mode_rank(r))
    });
    let mut out = format!("{CSV_HEADER}\n");
    for r in rows {
        out.push_str(&r);
        out.push('\n');
    }
    out
}

pub struct ConcordArgs {
    pub only: Option<String>,
    pub gloci: Option<PathBuf>,
    pub no_docker: bool,
    pub image: String,
    pub data_dir: PathBuf,
    /// Results go to `<out>/<version>`.
    pub out: PathBuf,
}

/// How commands run: in the bench image, or on the host.
struct Runner {
    docker: Option<DockerRun>,
    work: PathBuf,
}

struct DockerRun {
    image: String,
    cache: PathBuf,
    gloci: PathBuf,
    user: String,
}

/// The result of one command: whether it succeeded, its stderr, its time.
struct Ran {
    ok: bool,
    stderr: String,
    secs: f64,
}

impl Runner {
    /// Run `cmd` under `bash -o pipefail` with `LC_ALL=C`; its stderr goes to
    /// `<work>/<tag>.err`.
    fn run(&self, cmd: &str, tag: &str) -> Result<Ran> {
        let err_path = self.work.join(format!("{tag}.err"));
        let err = std::fs::File::create(&err_path)?;
        let mut c = match &self.docker {
            Some(d) => {
                let mut c = Command::new("docker");
                c.args(["run", "--rm", "--network", "none", "--name"])
                    .arg(format!("gloci-concord-{}-{tag}", std::process::id()))
                    .args(["--user", &d.user, "-e", "HOME=/tmp", "-e", "LC_ALL=C"])
                    .arg("-v")
                    .arg(format!("{}:/cache:ro", d.cache.display()))
                    .arg("-v")
                    .arg(format!("{}:/work", self.work.display()))
                    .arg("-v")
                    .arg(format!("{}:/usr/local/bin/gloci:ro", d.gloci.display()))
                    .args(["-w", "/work", &d.image, "bash", "-o", "pipefail", "-c", cmd]);
                c
            }
            None => host_shell(cmd, &self.work),
        };
        let t0 = std::time::Instant::now();
        let st = c
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(err)
            .status()
            .with_context(|| format!("running {cmd}"))?;
        let secs = t0.elapsed().as_secs_f64();
        let text = std::fs::read_to_string(&err_path).unwrap_or_default();
        let tail: Vec<&str> = text.lines().rev().take(15).collect();
        let stderr = tail.into_iter().rev().collect::<Vec<_>>().join("\n");
        Ok(Ran {
            ok: st.success(),
            stderr,
            secs,
        })
    }

    /// A work file as the commands see it.
    fn shown(&self, name: &str) -> String {
        match self.docker {
            Some(_) => format!("/work/{name}"),
            None => self.work.join(name).to_string_lossy().into_owned(),
        }
    }
}

/// `cmd` as the host runs it: `bash -o pipefail` with `LC_ALL=C` in `dir`, as
/// the image does.
fn host_shell(cmd: &str, dir: &Path) -> Command {
    let mut c = Command::new("bash");
    c.args(["-o", "pipefail", "-c", cmd])
        .env("LC_ALL", "C")
        .current_dir(dir);
    c
}

/// Whether `cmd` succeeds when run as the commands of a case run on the host.
pub fn host_command_succeeds(cmd: &str, dir: &Path) -> Result<bool> {
    Ok(host_shell(cmd, dir)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .with_context(|| format!("running {cmd}"))?
        .success())
}

/// Remove `name` and the files commands write beside it (`name.p.frq`, ...).
fn remove_outputs(dir: &Path, name: &str) {
    let _ = std::fs::remove_file(dir.join(name));
    let side = format!("{name}.");
    if let Ok(rd) = std::fs::read_dir(dir) {
        for e in rd.flatten() {
            if e.file_name().to_string_lossy().starts_with(&side) {
                let _ = std::fs::remove_file(e.path());
            }
        }
    }
}

/// A reader of an output, decompressing gzip (BGZF included) by its magic.
fn open_output(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut f = std::fs::File::open(path).with_context(|| path.display().to_string())?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    let f = std::fs::File::open(path)?;
    Ok(if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::with_capacity(
            1 << 20,
            flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f)),
        ))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    })
}

const SHOWN_BYTES: usize = 400;

fn shown_line(line: &[u8]) -> String {
    let s = String::from_utf8_lossy(line);
    if s.len() > SHOWN_BYTES {
        let mut cut = SHOWN_BYTES;
        while !s.is_char_boundary(cut) {
            cut -= 1;
        }
        format!("{}... ({} bytes)", &s[..cut], line.len())
    } else {
        s.into_owned()
    }
}

/// Up to `n` lines of `path` that `residual` left unmatched on this side
/// (`gloci_side`), each surplus line once.
fn unmatched_lines(
    mode: &str,
    path: &Path,
    residual: &mut Residual,
    gloci_side: bool,
    n: usize,
) -> Result<Vec<String>> {
    let mut out = Vec::new();
    if n == 0 {
        return Ok(out);
    }
    let mut r = open_output(path)?;
    let mut index = 0u64;
    gloci_testkit::differential::for_each_line_key(mode, &mut *r, |_, line| {
        if residual.take(gloci_side, index, line) {
            out.push(if mode == "numeric" {
                format!("line {}: {}", index + 1, shown_line(line))
            } else {
                shown_line(line)
            });
        }
        index += 1;
        out.len() < n
    })?;
    Ok(out)
}

/// The first line (1-based, among the lines the mode keeps) where the two
/// outputs differ, and the line of each side there.
fn first_order_difference(mode: &str, a: &Path, b: &Path) -> Result<Option<(u64, String, String)>> {
    let (mut ra, mut rb) = (open_output(a)?, open_output(b)?);
    let next = |r: &mut dyn BufRead| -> Result<Option<(u64, Vec<u8>)>> {
        let mut buf = Vec::new();
        loop {
            buf.clear();
            if r.read_until(b'\n', &mut buf)? == 0 {
                return Ok(None);
            }
            let line = buf.strip_suffix(b"\n").unwrap_or(&buf).to_vec();
            if let Some(k) = line_key(mode, &line) {
                return Ok(Some((k, line)));
            }
        }
    };
    let mut i = 0u64;
    loop {
        i += 1;
        match (next(&mut *ra)?, next(&mut *rb)?) {
            (None, None) => return Ok(None),
            (Some((ka, _)), Some((kb, _))) if ka == kb => {}
            (x, y) => {
                let s = |v: Option<(u64, Vec<u8>)>| {
                    v.map_or("(end)".to_string(), |(_, l)| shown_line(&l))
                };
                return Ok(Some((i, s(x), s(y))));
            }
        }
    }
}

struct RunReport<'a> {
    case: &'a ConcordCase,
    mode: &'a str,
    gloci_cmd: &'a str,
    reference_cmd: &'a str,
    expect: Expect<'a>,
    stats: Option<MatchStats>,
    verdict: Verdict,
    times: String,
    body: String,
}

fn write_report(dir: &Path, r: &RunReport) -> Result<()> {
    let mut t = String::new();
    t.push_str(&format!("case: {}\nmode: {}\n", r.case.name, r.mode));
    t.push_str(&format!(
        "command: gloci {} (against {})\ncompare: {}\n",
        r.case.command, r.case.reference_tool, r.case.compare
    ));
    t.push_str(&format!(
        "gloci: {}\nreference: {}\n",
        r.gloci_cmd, r.reference_cmd
    ));
    if let Some(a) = r.expect.anchor {
        t.push_str(&format!("documented difference: {a}"));
        if r.expect.order_only {
            t.push_str(" (order only: every line must match)");
        } else if let Some(m) = r.expect.min_rate {
            t.push_str(&format!(" (min_rate {m:.6})"));
        }
        t.push('\n');
    }
    t.push_str(&r.times);
    if let Some(s) = &r.stats {
        t.push_str(&format!(
            "n_gloci: {}\nn_reference: {}\nn_matching: {}\nmatch_rate: {}\nsame order: {}\n",
            s.n_gloci,
            s.n_reference,
            s.n_matching,
            rate_text(s),
            if s.same_order { "yes" } else { "no" }
        ));
    }
    t.push_str(&format!("verdict: {}\n", r.verdict.as_str()));
    t.push_str(&r.body);
    let path = dir.join(format!("{}.{}.txt", r.case.name, r.mode));
    std::fs::write(&path, t).with_context(|| path.display().to_string())
}

/// One row of `concordance.csv`: `match_rate` is a decimal from 0 to 1, six
/// places truncated ([`rate_text`]), empty with the counts when a command
/// failed; `anchors` holds the
/// run's documented difference as a plain name (no `#`; names would be
/// joined with `;`).
pub fn csv_row(
    c: &ConcordCase,
    mode: &str,
    stats: Option<&MatchStats>,
    v: Verdict,
    anchor: Option<&str>,
) -> String {
    let dataset = c.inputs.values().cloned().collect::<Vec<_>>().join("+");
    let nums = match stats {
        Some(s) => format!(
            "{},{},{},{}",
            s.n_gloci,
            s.n_reference,
            s.n_matching,
            rate_text(s)
        ),
        None => ",,,".to_string(),
    };
    [
        csv_field(&c.name),
        csv_field(&c.command),
        csv_field(&c.reference_tool),
        csv_field(&dataset),
        mode.to_string(),
        nums,
        v.as_str().to_string(),
        csv_field(anchor.map_or("", |a| a.trim_start_matches('#'))),
    ]
    .join(",")
}

/// The digest `datasets.toml` pins for an input: a derived input's sha256,
/// or a dataset's own checksum, as `algorithm:hex`; `None` when unpinned.
pub fn pinned_digest(cat: &Catalog, name: &str) -> Option<String> {
    if let Some(d) = cat.derived.iter().find(|d| d.name == name) {
        return (!d.sha256.is_empty()).then(|| format!("sha256:{}", d.sha256));
    }
    let d = cat.datasets.iter().find(|d| d.name == name)?;
    let (alg, hex) = d.checksum().ok()?;
    Some(format!("{alg}:{hex}"))
}

/// The key a cached reference output is kept under: the tools (the image
/// and its Id, or the host tool's version), every input with its digest, and
/// the command. The output is reused only while all of them are the same.
pub fn reference_stamp(tools: &str, inputs: &[(String, String, String)], cmd: &str) -> String {
    let mut s = format!("tools: {tools}\n");
    for (slot, name, digest) in inputs {
        s.push_str(&format!("input {slot}: {name} {digest}\n"));
    }
    s.push_str(&format!("command: {cmd}\n"));
    s
}

/// The Id of a docker image (`docker image inspect --format {{.Id}}`).
fn image_id(image: &str) -> Result<String> {
    let o = Command::new("docker")
        .args(["image", "inspect", "--format", "{{.Id}}", image])
        .output()
        .context("running docker image inspect")?;
    ensure!(o.status.success(), "docker image inspect {image} failed");
    Ok(String::from_utf8(o.stdout)?.trim().to_string())
}

/// The first line of `<tool> --version` on the host (the reference tool of a
/// `--no-docker` run; `reference_tool` names a binary: sort-bed, bedops,
/// bedtools, tabix, vcftools, bcftools, samtools), or a note that it has none.
fn host_tool_version(bin: &str) -> String {
    Command::new(bin)
        .arg("--version")
        .output()
        .ok()
        .map(|o| {
            // The first line with a digit (bedops prints its name first).
            let text = String::from_utf8_lossy(&[o.stdout, o.stderr].concat()).into_owned();
            let lines: Vec<&str> = text
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty())
                .collect();
            lines
                .iter()
                .find(|l| l.chars().any(|c| c.is_ascii_digit()))
                .or(lines.first())
                .map_or(String::new(), |l| l.to_string())
        })
        .filter(|v| !v.is_empty())
        .map_or_else(|| format!("{bin} (no --version)"), |v| format!("{bin} {v}"))
}

fn id_of(flag: &str) -> Result<String> {
    let o = Command::new("id").arg(flag).output()?;
    Ok(String::from_utf8(o.stdout)?.trim().to_string())
}

/// Run the concordance cases (one with `--only`); true when every run passes.
pub fn run(a: &ConcordArgs) -> Result<bool> {
    let started_at = crate::paper::data::now_utc();
    let cases = load_cases()?;
    if let Some(o) = &a.only {
        ensure!(
            cases.iter().any(|c| &c.name == o),
            "--only {o}: not a concordance case (see {CASES_PATH})"
        );
    }
    let gloci = resolve_gloci(a.gloci.as_deref())?;
    let cat: Catalog = crate::paper::fetch::load_catalog()?;
    let cache_dir = bench_cache_dir(&a.data_dir);
    ensure!(
        cache_dir.is_dir(),
        "no dataset cache at {}: run `cargo xtask paper fetch --data-dir {}`",
        cache_dir.display(),
        a.data_dir.display()
    );
    let cache = std::fs::canonicalize(&cache_dir)?;
    let results = a.out.join(env!("CARGO_PKG_VERSION"));
    let detail = results.join("concordance");
    std::fs::create_dir_all(&detail)?;
    let results = std::fs::canonicalize(&results)?;
    let detail = std::fs::canonicalize(&detail)?;
    let work = target_dir().join("concordance-work");
    std::fs::create_dir_all(&work)?;
    let work = std::fs::canonicalize(&work)?;
    let docker = if a.no_docker {
        None
    } else {
        let ok = Command::new("docker")
            .args(["image", "inspect", &a.image])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if !ok {
            bail!(
                "image {0} not found: docker build -t {0} -f xtask/docker/Dockerfile.bench xtask/docker (or pass --no-docker to use host tools)",
                a.image
            );
        }
        Some(DockerRun {
            image: a.image.clone(),
            cache: cache.clone(),
            gloci: gloci.clone(),
            user: format!("{}:{}", id_of("-u")?, id_of("-g")?),
        })
    };
    let runner = Runner { docker, work };
    let image_key = match &runner.docker {
        Some(d) => Some(format!("image {} {}", d.image, image_id(&d.image)?)),
        None => None,
    };

    let mut rows = Vec::new();
    let mut ran = Vec::new();
    let mut all_ok = true;
    for c in cases
        .iter()
        .filter(|c| a.only.as_ref().is_none_or(|o| *o == c.name))
    {
        ran.push(c.name.clone());
        eprintln!(
            "== {} (gloci {} against {})",
            c.name, c.command, c.reference_tool
        );
        let mut slots = Slots {
            gloci: match runner.docker {
                Some(_) => "gloci".into(),
                None => gloci.to_string_lossy().into_owned(),
            },
            in_a: String::new(),
            in_b: String::new(),
            in_c: String::new(),
            genome: String::new(),
            out: String::new(),
        };
        for (slot, name) in &c.inputs {
            let p = input_path(&cat, name, &cache)?;
            ensure!(
                p.exists(),
                "{}: input `{name}` is missing at {}: run `cargo xtask paper fetch --data-dir {}`",
                c.name,
                p.display(),
                a.data_dir.display()
            );
            let shown = match runner.docker {
                Some(_) => {
                    let rel = p.strip_prefix(&cache).with_context(|| {
                        format!(
                            "{}: {} is outside the cache {}",
                            c.name,
                            p.display(),
                            cache.display()
                        )
                    })?;
                    format!("/cache/{}", rel.display())
                }
                None => p.to_string_lossy().into_owned(),
            };
            match slot.as_str() {
                "a" => slots.in_a = shown,
                "b" => slots.in_b = shown,
                "c" => slots.in_c = shown,
                _ => slots.genome = shown,
            }
        }

        // The reference output, reused while its command and image are the same.
        let ref_name = format!("{}.ref.out", c.name);
        slots.out = runner.shown(&ref_name);
        let ref_cmd = fill(&c.reference, &slots);
        let tools = match &image_key {
            Some(k) => k.clone(),
            None => format!("host {}", host_tool_version(&c.reference_tool)),
        };
        let digests: Vec<(String, String, String)> = c
            .inputs
            .iter()
            .map(|(slot, name)| {
                let digest = match pinned_digest(&cat, name) {
                    Some(d) => d,
                    None => {
                        let p = input_path(&cat, name, &cache).unwrap_or_default();
                        let m = std::fs::metadata(&p).ok();
                        format!(
                            "unpinned size {} modified {:?}",
                            m.as_ref().map_or(0, |m| m.len()),
                            m.and_then(|m| m.modified().ok())
                        )
                    }
                };
                (slot.clone(), name.clone(), digest)
            })
            .collect();
        let stamp = reference_stamp(&tools, &digests, &ref_cmd);
        let stamp_path = runner.work.join(format!("{}.ref.cmd", c.name));
        let ref_path = runner.work.join(&ref_name);
        let reused =
            ref_path.exists() && std::fs::read_to_string(&stamp_path).is_ok_and(|s| s == stamp);
        let (ref_ok, ref_err, ref_time) = if reused {
            eprintln!("  reference output reused ({})", ref_path.display());
            (true, String::new(), "reference: reused\n".to_string())
        } else {
            let _ = std::fs::remove_file(&stamp_path);
            remove_outputs(&runner.work, &ref_name);
            let r = runner.run(&ref_cmd, &format!("{}.ref", c.name))?;
            let ok = r.ok && ref_path.exists();
            eprintln!(
                "  reference {:.1} s{}",
                r.secs,
                if ok { "" } else { " FAILED" }
            );
            if ok {
                std::fs::write(&stamp_path, &stamp)?;
            }
            (ok, r.stderr, format!("reference: {:.1} s\n", r.secs))
        };

        let modes: &[(&str, bool)] = if c.compat.is_some() {
            &[("default", false), ("compat", true)]
        } else {
            &[("default", false)]
        };
        for &(mode, compat) in modes {
            let expect = c.expect(compat);
            let anchor = expect.anchor;
            let out_name = format!("{}.{mode}.out", c.name);
            slots.out = runner.shown(&out_name);
            let gloci_cmd = fill(&gloci_command(c, compat), &slots);
            let mut times = ref_time.clone();
            let mut body = String::new();
            let mut stats = None;
            if !ref_ok {
                body.push_str(&format!("## reference failed\n{ref_err}\n"));
            } else {
                remove_outputs(&runner.work, &out_name);
                let r = runner.run(&gloci_cmd, &format!("{}.{mode}", c.name))?;
                let out_path = runner.work.join(&out_name);
                times.push_str(&format!("gloci: {:.1} s\n", r.secs));
                if !r.ok || !out_path.exists() {
                    body.push_str(&format!("## gloci failed\n{}\n", r.stderr));
                } else {
                    let (s, mut residual) = match_stats_with_residual(
                        &c.compare,
                        c.ref_decimals,
                        &mut *open_output(&out_path)?,
                        &mut *open_output(&ref_path)?,
                    )
                    .with_context(|| format!("{}: comparing outputs", c.name))?;
                    if s.rate() < 1.0 {
                        let g = unmatched_lines(&c.compare, &out_path, &mut residual, true, 10)?;
                        let rf = unmatched_lines(&c.compare, &ref_path, &mut residual, false, 10)?;
                        body.push_str(&format!(
                            "## only in gloci ({}, first {})\n",
                            s.n_gloci - s.n_matching,
                            g.len()
                        ));
                        for l in g {
                            body.push_str(&l);
                            body.push('\n');
                        }
                        body.push_str(&format!(
                            "## only in the reference ({}, first {})\n",
                            s.n_reference - s.n_matching,
                            rf.len()
                        ));
                        for l in rf {
                            body.push_str(&l);
                            body.push('\n');
                        }
                    } else if !s.same_order
                        && let Some((i, g, rf)) =
                            first_order_difference(&c.compare, &out_path, &ref_path)?
                    {
                        body.push_str(&format!(
                            "## first difference in order, line {i}\ngloci:     {g}\nreference: {rf}\n"
                        ));
                    }
                    stats = Some(s);
                }
                remove_outputs(&runner.work, &out_name);
            }
            let v = verdict(stats.as_ref(), &expect);
            all_ok &= v.passes();
            let summary = match &stats {
                Some(s) => format!(
                    "{}/{} gloci, {} reference, rate {}{}",
                    s.n_matching,
                    s.n_gloci,
                    s.n_reference,
                    rate_text(s),
                    if s.same_order { "" } else { ", order differs" }
                ),
                None => "no comparison".to_string(),
            };
            eprintln!("  {mode:<8} {summary}: {}", v.as_str());
            if v == Verdict::AnchorUnused {
                eprintln!("  anchor unused: {} ({})", anchor.unwrap_or(""), c.name);
            }
            write_report(
                &detail,
                &RunReport {
                    case: c,
                    mode,
                    gloci_cmd: &gloci_cmd,
                    reference_cmd: &ref_cmd,
                    expect,
                    stats,
                    verdict: v,
                    times,
                    body,
                },
            )?;
            rows.push(csv_row(c, mode, stats.as_ref(), v, anchor));
        }
    }

    let csv_path = results.join("concordance.csv");
    let old = std::fs::read_to_string(&csv_path).ok();
    let order: Vec<String> = cases.iter().map(|c| c.name.clone()).collect();
    let text = merge_csv(old.as_deref(), &rows, &ran, &order);
    let tmp = results.join("concordance.csv.tmp");
    std::fs::File::create(&tmp)?.write_all(text.as_bytes())?;
    std::fs::rename(&tmp, &csv_path)?;
    eprintln!("wrote {}", csv_path.display());
    write_run_records(a, &results, &gloci, &cases, &started_at)?;
    Ok(all_ok)
}

/// `environment.json` and `repo-snapshot.json` beside `concordance.csv`, which
/// the paper tables read, when the directory has none yet: a benchmark run
/// into the same directory writes its own, which also list the tools and are
/// kept.
pub fn write_run_records(
    a: &ConcordArgs,
    dir: &Path,
    gloci: &Path,
    cases: &[ConcordCase],
    started_at: &str,
) -> Result<()> {
    use crate::paper::{environment, snapshot};
    let (env_path, snap_path) = (dir.join("environment.json"), dir.join("repo-snapshot.json"));
    if env_path.exists() && snap_path.exists() {
        return Ok(());
    }
    let root = workspace_root();
    let commit = environment::git_commit(&root);
    let version = env!("CARGO_PKG_VERSION");
    let tools: Vec<String> = cases
        .iter()
        .map(|c| c.reference_tool.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if !env_path.exists() {
        let env = environment::environment(&environment::EnvArgs {
            version,
            commit: commit.clone(),
            started_at,
            finished_at: &crate::paper::data::now_utc(),
            image: (!a.no_docker).then_some(a.image.as_str()),
            gloci: Some(gloci),
            wheel: None,
            cases: &[],
            extra_tools: &tools,
        })?;
        std::fs::write(&env_path, serde_json::to_string_pretty(&env)? + "\n")?;
    }
    if !snap_path.exists() {
        let cat = crate::paper::fetch::load_catalog()?;
        let snap = snapshot::snapshot(&root, &cat, &a.data_dir, version, commit.as_deref())?;
        std::fs::write(&snap_path, serde_json::to_string_pretty(&snap)? + "\n")?;
    }
    Ok(())
}
