//! Do the commands of a benchmark case give the same result?
//!
//! A benchmark that times commands without looking at what they write can time the wrong
//! thing: on 2026-11-13 the paper cases had a `bedmap` that printed `NAN` for every row in a
//! second, a `bedops --partition` fed input in an order it does not read, and tools compared
//! on work that was not the same. So a case says which of its commands must agree
//! (`[[case.same]]`) and how, the runner runs each of those commands once into its own file
//! before it times anything, and this module compares the files. A command that is compared
//! with no other says why (`[case.not_compared]`).
//!
//! ```toml
//! [[case.same]]
//! tools = ["gloci", "bedtools", "bedops"]   # the first is compared with each of the others
//! compare = "exact"                         # the default
//! ```
//!
//! | `compare` | the outputs agree when |
//! |---|---|
//! | `exact` | they are the same bytes |
//! | `unordered` | they hold the same lines, in any order |
//! | `distinct` | they hold the same distinct lines, in any order and however often (`bedtools pairtopair` writes a pair that hits in both orientations twice) |
//! | `numeric` | line by line, the text fields are the same and the numbers agree: within a relative 1e-5, or, with `decimals = N`, within half a unit of the N-th decimal of the other tool's number |
//! | `vcf-records` | the lines that do not start with `#` are the same bytes |
//! | `vcf-sites` | those lines have the same CHROM, POS, REF and ALT, whatever the case of the letters |
//! | `bam-records` | `samtools view` prints the same records for both |
//! | `bcf-records` | `bcftools view -H` prints the same records for both |
//!
//! `columns = N` compares the first N tab-separated columns only (`exact`, `unordered` and
//! `distinct`).
//! `missing = "NAN"` (`numeric`) reads a field that is exactly that, in the outputs after the
//! first, as `.`: bedmap prints `NAN` for an interval nothing maps to.
//! Output that is gzip or BGZF is read through its compression. The records of a BAM or BCF
//! are read from the decoder's standard output and never written down: those of a
//! 50-million-read BAM are 12 GB of text.

use anyhow::{Context, Result, bail, ensure};
use gloci_testkit::differential::match_stats_with_residual;
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};

/// Commands of a case whose outputs must agree: the first with each of the others.
#[derive(Deserialize, Clone, Debug, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Same {
    pub tools: Vec<String>,
    #[serde(default = "exact")]
    pub compare: String,
    /// Compare the first N tab-separated columns only.
    #[serde(default)]
    pub columns: Option<usize>,
    /// `numeric`: the tools after the first print this many decimals.
    #[serde(default)]
    pub decimals: Option<u32>,
    /// `numeric`: what the tools after the first print for a value that is not there
    /// (bedmap's `NAN`), read as the first tool's `.`.
    #[serde(default)]
    pub missing: Option<String>,
}

fn exact() -> String {
    "exact".into()
}

pub const MODES: &[&str] = &[
    "exact",
    "unordered",
    "distinct",
    "numeric",
    "vcf-records",
    "vcf-sites",
    "bam-records",
    "bcf-records",
];

/// The command that prints the records of a binary output as text, for the modes that
/// compare decoded records; `None` for the modes that read the output as it is.
pub fn decoder(mode: &str) -> Option<&'static str> {
    match mode {
        "bam-records" => Some("samtools view"),
        "bcf-records" => Some("bcftools view -H"),
        _ => None,
    }
}

/// Check the `same` groups and `not_compared` reasons of the case `name` against its
/// commands.
pub fn check_schema(
    name: &str,
    same: &[Same],
    not_compared: &BTreeMap<String, String>,
    commands: &BTreeMap<String, String>,
) -> Result<()> {
    for g in same {
        ensure!(
            g.tools.len() >= 2,
            "{name}: a `same` group needs at least two tools, found {:?}",
            g.tools
        );
        ensure!(
            MODES.contains(&g.compare.as_str()),
            "{name}: unknown compare mode `{}` (use one of {})",
            g.compare,
            MODES.join(", ")
        );
        for t in &g.tools {
            ensure!(
                commands.contains_key(t),
                "{name}: `same` names `{t}`, which is not a command of the case"
            );
        }
        let mut seen = g.tools.clone();
        seen.sort();
        seen.dedup();
        ensure!(
            seen.len() == g.tools.len(),
            "{name}: a `same` group names a tool twice: {:?}",
            g.tools
        );
        ensure!(
            g.columns.is_none() || matches!(g.compare.as_str(), "exact" | "unordered" | "distinct"),
            "{name}: `columns` applies to compare = \"exact\", \"unordered\" and \"distinct\" only"
        );
        ensure!(g.columns != Some(0), "{name}: `columns` must be at least 1");
        ensure!(
            g.decimals.is_none() || g.compare == "numeric",
            "{name}: `decimals` applies to compare = \"numeric\" only"
        );
        ensure!(
            g.missing.is_none() || g.compare == "numeric",
            "{name}: `missing` applies to compare = \"numeric\" only"
        );
    }
    for (t, why) in not_compared {
        ensure!(
            commands.contains_key(t),
            "{name}: `not_compared` names `{t}`, which is not a command of the case"
        );
        ensure!(
            !why.trim().is_empty(),
            "{name}: `not_compared.{t}` needs a reason"
        );
        ensure!(
            !same.iter().any(|g| g.tools.contains(t)),
            "{name}: `{t}` is in a `same` group and in `not_compared`"
        );
    }
    Ok(())
}

/// The commands of a case that must be accounted for and are not: every `gloci` command
/// (but the Python library's, which the Python parity tests cover) and every tool a target
/// names, when it is in no `same` group and has no `not_compared` reason.
pub fn unaccounted(
    same: &[Same],
    not_compared: &BTreeMap<String, String>,
    commands: &BTreeMap<String, String>,
    targets: &[String],
) -> Vec<String> {
    let named_by_target = |tool: &str| {
        targets
            .iter()
            .any(|t| t.split_whitespace().any(|word| word == tool))
    };
    commands
        .keys()
        .filter(|t| t.as_str() != "gloci-py")
        .filter(|t| *t == "gloci" || t.starts_with("gloci-") || named_by_target(t))
        .filter(|t| !same.iter().any(|g| g.tools.contains(t)) && !not_compared.contains_key(*t))
        .cloned()
        .collect()
}

fn open(path: &Path) -> Result<Box<dyn BufRead>> {
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

/// What a mode compares of one line (without its newline), written to `out`; false for a
/// line the mode skips.
fn view(mode: &str, columns: Option<usize>, line: &[u8], out: &mut Vec<u8>) -> bool {
    out.clear();
    match mode {
        "vcf-records" | "vcf-sites" if line.starts_with(b"#") => return false,
        "vcf-sites" => {
            for (i, f) in line.split(|&b| b == b'\t').enumerate().take(5) {
                if matches!(i, 0 | 1 | 3 | 4) {
                    if i > 0 {
                        out.push(b'\t');
                    }
                    out.extend(f.iter().map(u8::to_ascii_uppercase));
                }
            }
        }
        _ => match columns {
            Some(n) => {
                let end = line
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| **b == b'\t')
                    .nth(n - 1)
                    .map_or(line.len(), |(i, _)| i);
                out.extend_from_slice(&line[..end]);
            }
            None => out.extend_from_slice(line),
        },
    }
    true
}

/// The next line of `r` that the mode keeps, as the mode compares it, in `out`; false at
/// the end.
fn next_view(
    mode: &str,
    columns: Option<usize>,
    r: &mut dyn BufRead,
    raw: &mut Vec<u8>,
    out: &mut Vec<u8>,
) -> Result<bool> {
    loop {
        raw.clear();
        if r.read_until(b'\n', raw)? == 0 {
            return Ok(false);
        }
        let line = raw.strip_suffix(b"\n").unwrap_or(raw);
        if view(mode, columns, line, out) {
            return Ok(true);
        }
    }
}

fn shown(line: &[u8]) -> String {
    let s = String::from_utf8_lossy(line);
    let mut cut = s.len().min(160);
    while !s.is_char_boundary(cut) {
        cut -= 1;
    }
    if cut < s.len() {
        format!("{}... ({} bytes)", &s[..cut], s.len())
    } else {
        s.into_owned()
    }
}

/// A tool's output with every field that is exactly `missing` replaced by `.`.
struct MissingAsDot<'a> {
    inner: &'a mut dyn BufRead,
    missing: Vec<u8>,
    raw: Vec<u8>,
    line: Vec<u8>,
    at: usize,
}

impl<'a> MissingAsDot<'a> {
    fn new(inner: &'a mut dyn BufRead, missing: &[u8]) -> MissingAsDot<'a> {
        MissingAsDot {
            inner,
            missing: missing.to_vec(),
            raw: Vec::new(),
            line: Vec::new(),
            at: 0,
        }
    }
}

impl Read for MissingAsDot<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let have = self.fill_buf()?;
        let n = have.len().min(buf.len());
        buf[..n].copy_from_slice(&have[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for MissingAsDot<'_> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.at == self.line.len() {
            self.raw.clear();
            self.line.clear();
            self.at = 0;
            if self.inner.read_until(b'\n', &mut self.raw)? > 0 {
                let (body, eol) = match self.raw.strip_suffix(b"\n") {
                    Some(b) => (b, &b"\n"[..]),
                    None => (&self.raw[..], &b""[..]),
                };
                for (i, f) in body.split(|&b| b == b'\t').enumerate() {
                    if i > 0 {
                        self.line.push(b'\t');
                    }
                    if f == self.missing {
                        self.line.push(b'.');
                    } else {
                        self.line.extend_from_slice(f);
                    }
                }
                self.line.extend_from_slice(eol);
            }
        }
        Ok(&self.line[self.at..])
    }

    fn consume(&mut self, n: usize) {
        self.at = (self.at + n).min(self.line.len());
    }
}

/// Whether the outputs `a` (the group's first tool) and `b` agree under the group `g`:
/// `Ok(None)` when they do, else what differs first.
pub fn differ(g: &Same, a: &Path, b: &Path) -> Result<Option<String>> {
    differ_streams(g, &mut *open(a)?, &mut *open(b)?)
}

/// [`differ`] for two outputs that are already open.
fn differ_streams(g: &Same, ra: &mut dyn BufRead, rb: &mut dyn BufRead) -> Result<Option<String>> {
    match g.compare.as_str() {
        "numeric" => {
            let mut dotted;
            let rb: &mut dyn BufRead = match &g.missing {
                Some(m) => {
                    dotted = MissingAsDot::new(rb, m.as_bytes());
                    &mut dotted
                }
                None => rb,
            };
            let (s, _) = match_stats_with_residual("numeric", g.decimals, &mut *ra, &mut *rb)?;
            Ok(
                (s.n_matching != s.n_gloci || s.n_gloci != s.n_reference).then(|| {
                    format!(
                        "{} and {} lines, of which {} agree in place",
                        s.n_gloci, s.n_reference, s.n_matching
                    )
                }),
            )
        }
        "unordered" => {
            let mut counts: HashMap<Vec<u8>, i64> = HashMap::new();
            let (mut raw, mut v) = (Vec::new(), Vec::new());
            let (mut na, mut nb) = (0u64, 0u64);
            while next_view("unordered", g.columns, &mut *ra, &mut raw, &mut v)? {
                *counts.entry(v.clone()).or_default() += 1;
                na += 1;
            }
            while next_view("unordered", g.columns, &mut *rb, &mut raw, &mut v)? {
                *counts.entry(v.clone()).or_default() -= 1;
                nb += 1;
            }
            let odd: Vec<(&Vec<u8>, &i64)> = counts.iter().filter(|(_, c)| **c != 0).collect();
            Ok(odd.first().map(|(line, c)| {
                format!(
                    "{na} and {nb} lines; {} differ, such as this one, {} more often in the {}: {}",
                    odd.len(),
                    c.abs(),
                    if **c > 0 { "first" } else { "second" },
                    shown(line)
                )
            }))
        }
        "distinct" => {
            // 1: in the first output, 2: in the second, 3: in both.
            let mut seen: HashMap<Vec<u8>, u8> = HashMap::new();
            let (mut raw, mut v) = (Vec::new(), Vec::new());
            let (mut na, mut nb) = (0u64, 0u64);
            while next_view("distinct", g.columns, &mut *ra, &mut raw, &mut v)? {
                *seen.entry(v.clone()).or_default() |= 1;
                na += 1;
            }
            while next_view("distinct", g.columns, &mut *rb, &mut raw, &mut v)? {
                *seen.entry(v.clone()).or_default() |= 2;
                nb += 1;
            }
            let odd: Vec<(&Vec<u8>, &u8)> = seen.iter().filter(|(_, s)| **s != 3).collect();
            Ok(odd.first().map(|(line, s)| {
                format!(
                    "{na} and {nb} lines; {} distinct lines are in one output only, such as \
                     this one, in the {} only: {}",
                    odd.len(),
                    if **s == 1 { "first" } else { "second" },
                    shown(line)
                )
            }))
        }
        mode => {
            // exact, vcf-records, vcf-sites, and the decoded records of bam-records and
            // bcf-records: line by line, in order.
            let mode = if decoder(mode).is_some() {
                "exact"
            } else {
                mode
            };
            let (mut raw_a, mut raw_b, mut va, mut vb) =
                (Vec::new(), Vec::new(), Vec::new(), Vec::new());
            let mut n = 0u64;
            loop {
                let more_a = next_view(mode, g.columns, &mut *ra, &mut raw_a, &mut va)?;
                let more_b = next_view(mode, g.columns, &mut *rb, &mut raw_b, &mut vb)?;
                n += 1;
                match (more_a, more_b) {
                    (false, false) => return Ok(None),
                    (true, true) if va == vb => {}
                    (true, true) => {
                        return Ok(Some(format!(
                            "line {n} differs:\n      {}\n      {}",
                            shown(&va),
                            shown(&vb)
                        )));
                    }
                    (true, false) => {
                        return Ok(Some(format!(
                            "the second output ends after {} lines; the first goes on: {}",
                            n - 1,
                            shown(&va)
                        )));
                    }
                    (false, true) => {
                        return Ok(Some(format!(
                            "the first output ends after {} lines; the second goes on: {}",
                            n - 1,
                            shown(&vb)
                        )));
                    }
                }
            }
        }
    }
}

/// How the comparison of one group went.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupReport {
    pub tools: Vec<String>,
    pub compare: String,
    /// Tools of the group that were not run (not installed).
    pub missing: Vec<String>,
    /// `(tool, what differs)` for each tool that disagrees with the first.
    pub different: Vec<(String, String)>,
}

impl GroupReport {
    pub fn ok(&self) -> bool {
        self.different.is_empty()
    }

    /// One line for the log: `same (exact): gloci = bedtools = bedops`.
    pub fn line(&self) -> String {
        let ran: Vec<&str> = self
            .tools
            .iter()
            .filter(|t| !self.missing.contains(t))
            .map(String::as_str)
            .collect();
        let mut s = if ran.len() < 2 {
            format!("not checked ({}): {}", self.compare, self.tools.join(", "))
        } else if self.ok() {
            format!("same ({}): {}", self.compare, ran.join(" = "))
        } else {
            format!("DIFFERENT ({}): {}", self.compare, ran.join(", "))
        };
        if !self.missing.is_empty() {
            s.push_str(&format!(" [not installed: {}]", self.missing.join(", ")));
        }
        for (t, why) in &self.different {
            s.push_str(&format!("\n    {} against {t}: {why}", ran[0]));
        }
        s
    }
}

/// Where the output of a tool is read from.
#[derive(Debug)]
pub enum Output {
    /// A file, read through gzip or BGZF when it is compressed.
    File(PathBuf),
    /// The records of a binary output: the standard output of the command that prints them
    /// ([`decoder`]). Nothing is written to disk.
    Records(Command),
}

/// An output being read, and the command that prints it when there is one.
struct Stream {
    reader: Box<dyn BufRead>,
    child: Option<(Child, String)>,
}

impl Stream {
    fn open(o: Output) -> Result<Stream> {
        match o {
            Output::File(p) => Ok(Stream {
                reader: open(&p)?,
                child: None,
            }),
            Output::Records(mut c) => {
                let what = format!("{c:?}");
                let mut child = c
                    .stdin(Stdio::null())
                    .stdout(Stdio::piped())
                    .spawn()
                    .with_context(|| format!("running {what}"))?;
                let out = child.stdout.take().expect("a piped stdout");
                Ok(Stream {
                    reader: Box::new(BufReader::with_capacity(1 << 20, out)),
                    child: Some((child, what)),
                })
            }
        }
    }

    /// Wait for the command behind the stream, and fail when it failed: an output that
    /// could not be read must not pass for one that agrees. What is left of its output is
    /// read to the end first, so that the command ends by itself (killing `docker run`
    /// would leave its container running).
    fn finish(mut self) -> Result<()> {
        if let Some((mut child, what)) = self.child.take() {
            std::io::copy(&mut self.reader, &mut std::io::sink())?;
            let status = child.wait()?;
            ensure!(
                status.success(),
                "cannot read an output to compare it: {what} ({status})"
            );
        }
        Ok(())
    }
}

/// Compare the outputs of the group `g`: `output(tool)` is where a tool's output is read
/// from, or `None` for a tool that was not run.
pub fn compare_group(
    g: &Same,
    output: &dyn Fn(&str) -> Result<Option<Output>>,
) -> Result<GroupReport> {
    let mut report = GroupReport {
        tools: g.tools.clone(),
        compare: g.compare.clone(),
        missing: Vec::new(),
        different: Vec::new(),
    };
    let mut ran: Vec<&str> = Vec::new();
    for t in &g.tools {
        match output(t)? {
            Some(_) => ran.push(t),
            None => report.missing.push(t.clone()),
        }
    }
    let Some((first, rest)) = ran.split_first() else {
        return Ok(report);
    };
    let wrote = |o: &Output| !matches!(o, Output::File(p) if !p.exists());
    for t in rest {
        // The first output is read again for each of the others.
        let (a, b) = (
            output(first)?.expect("a tool that ran"),
            output(t)?.expect("a tool that ran"),
        );
        if !wrote(&a) {
            bail!("{first} wrote no output ({a:?})");
        }
        if !wrote(&b) {
            report
                .different
                .push((t.to_string(), "it wrote no output".into()));
            continue;
        }
        let (mut sa, mut sb) = (Stream::open(a)?, Stream::open(b)?);
        let why = differ_streams(g, &mut *sa.reader, &mut *sb.reader)?;
        sa.finish()?;
        sb.finish()?;
        if let Some(why) = why {
            report.different.push((t.to_string(), why));
        }
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(compare: &str, columns: Option<usize>, decimals: Option<u32>) -> Same {
        Same {
            tools: vec!["a".into(), "b".into()],
            compare: compare.into(),
            columns,
            decimals,
            missing: None,
        }
    }

    fn cmp(g: &Same, a: &[u8], b: &[u8]) -> Option<String> {
        let d = tempfile::tempdir().unwrap();
        let (pa, pb) = (d.path().join("a"), d.path().join("b"));
        std::fs::write(&pa, a).unwrap();
        std::fs::write(&pb, b).unwrap();
        differ(g, &pa, &pb).unwrap()
    }

    #[test]
    fn exact_is_bytes_in_order() {
        let g = group("exact", None, None);
        assert_eq!(cmp(&g, b"a\tb\nc\n", b"a\tb\nc\n"), None);
        let why = cmp(&g, b"a\nb\n", b"b\na\n").unwrap();
        assert!(why.starts_with("line 1 differs"), "{why}");
        let why = cmp(&g, b"a\nb\n", b"a\n").unwrap();
        assert!(
            why.contains("the second output ends after 1 lines"),
            "{why}"
        );
        let why = cmp(&g, b"a\n", b"a\nb\n").unwrap();
        assert!(why.contains("the first output ends after 1 lines"), "{why}");
        // The bedmap fault of 2026-11-13: the same number of rows, other values.
        let why = cmp(&g, b"x\t0.403605\ny\t0\n", b"x\t0.000000\ny\tNAN\n").unwrap();
        assert!(why.starts_with("line 1 differs"), "{why}");
    }

    #[test]
    fn unordered_is_the_same_lines() {
        let g = group("unordered", None, None);
        assert_eq!(cmp(&g, b"a\nb\nb\n", b"b\na\nb\n"), None);
        let why = cmp(&g, b"a\nb\nb\n", b"b\na\na\n").unwrap();
        assert!(why.contains("3 and 3 lines; 2 differ"), "{why}");
    }

    #[test]
    fn distinct_ignores_order_and_repeats() {
        // bedtools pairtopair writes a pair that hits in both orientations twice.
        let g = group("distinct", None, None);
        assert_eq!(cmp(&g, b"a\nb\n", b"b\na\nb\n"), None);
        let why = cmp(&g, b"a\nb\n", b"b\nb\nc\n").unwrap();
        assert!(
            why.contains("2 and 3 lines; 2 distinct lines are in one output only"),
            "{why}"
        );
    }

    #[test]
    fn columns_cut_each_line() {
        // multiinter writes columns gloci does not.
        let g = group("exact", Some(3), None);
        assert_eq!(cmp(&g, b"c\t1\t2\n", b"c\t1\t2\t1\t1,2\n"), None);
        assert!(cmp(&g, b"c\t1\t2\n", b"c\t1\t3\tx\n").is_some());
        // A line with fewer columns is compared whole.
        assert_eq!(cmp(&g, b"c\t1\n", b"c\t1\n"), None);
    }

    #[test]
    fn numeric_allows_for_printed_decimals() {
        let g = group("numeric", None, None);
        assert_eq!(cmp(&g, b"a\t0.5\n", b"a\t0.5000001\n"), None);
        assert!(cmp(&g, b"a\t0.5\n", b"a\t0.6\n").is_some());
        assert!(cmp(&g, b"a\t0.5\n", b"b\t0.5\n").is_some());
        // bedmap prints six decimals: 0.0000123 is 0.000012 there.
        assert!(cmp(&g, b"a\t0.0000123\n", b"a\t0.000012\n").is_some());
        let g6 = group("numeric", None, Some(6));
        assert_eq!(cmp(&g6, b"a\t0.0000123\n", b"a\t0.000012\n"), None);
        assert!(cmp(&g6, b"a\t0.0000123\n", b"a\t0.000014\n").is_some());
        let why = cmp(&g, b"a\t1\nb\t2\n", b"a\t1\n").unwrap();
        assert_eq!(why, "2 and 1 lines, of which 1 agree in place");
        // bedmap prints NAN where nothing maps; gloci and bedtools print `.`.
        assert!(cmp(&g6, b"a\t.\nb\t2\n", b"a\tNAN\nb\t2.000000\n").is_some());
        let nan = Same {
            missing: Some("NAN".into()),
            ..group("numeric", None, Some(6))
        };
        assert_eq!(cmp(&nan, b"a\t.\nb\t2\n", b"a\tNAN\nb\t2.000000\n"), None);
        // Only a whole field counts, and only in the outputs after the first.
        assert!(cmp(&nan, b"a\t.\n", b"a\tNANO\n").is_some());
        assert!(cmp(&nan, b"a\tNAN\n", b"a\t.\n").is_some());
    }

    #[test]
    fn vcf_modes_skip_the_header() {
        let g = group("vcf-records", None, None);
        assert_eq!(
            cmp(
                &g,
                b"##a\n#CHROM\nchr1\t5\t.\tA\tC\n",
                b"##b\n##c\n#CHROM\nchr1\t5\t.\tA\tC\n"
            ),
            None
        );
        assert!(
            cmp(
                &g,
                b"#h\nchr1\t5\t.\tA\tC\tx\n",
                b"#h\nchr1\t5\t.\tA\tC\ty\n"
            )
            .is_some()
        );
        // Sites: CHROM, POS, REF and ALT, in any letter case (vcftools upper-cases ALT).
        let g = group("vcf-sites", None, None);
        assert_eq!(
            cmp(
                &g,
                b"#h\nchr1\t5\tid\tA\t<DEL:ME:Alu>\t9\tPASS\tX=1\n",
                b"chr1\t5\t.\tA\t<DEL:ME:ALU>\t.\t.\tY=2\n"
            ),
            None
        );
        assert!(cmp(&g, b"chr1\t5\t.\tA\tC\n", b"chr1\t6\t.\tA\tC\n").is_some());
    }

    #[test]
    fn compressed_outputs_are_read_through() {
        use std::io::Write;
        let gz = |b: &[u8]| {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            e.write_all(b).unwrap();
            e.finish().unwrap()
        };
        let g = group("exact", None, None);
        assert_eq!(cmp(&g, &gz(b"a\nb\n"), b"a\nb\n"), None);
        assert!(cmp(&g, &gz(b"a\nb\n"), &gz(b"a\nc\n")).is_some());
    }

    fn commands(names: &[&str]) -> BTreeMap<String, String> {
        names.iter().map(|n| (n.to_string(), "x".into())).collect()
    }

    #[test]
    fn schema_is_checked() {
        let cmds = commands(&["gloci", "bedtools", "bedops"]);
        let none = BTreeMap::new();
        let ok = |g: Same| check_schema("c", &[g], &none, &cmds);
        assert!(
            ok(Same {
                tools: vec!["gloci".into(), "bedtools".into()],
                ..group("exact", None, None)
            })
            .is_ok()
        );
        let e = |g: Same| ok(g).unwrap_err().to_string();
        assert!(
            e(Same {
                tools: vec!["gloci".into()],
                ..group("exact", None, None)
            })
            .contains("at least two")
        );
        assert!(
            e(Same {
                tools: vec!["gloci".into(), "gia".into()],
                ..group("exact", None, None)
            })
            .contains("not a command")
        );
        assert!(
            e(Same {
                tools: vec!["gloci".into(), "gloci".into()],
                ..group("exact", None, None)
            })
            .contains("twice")
        );
        let both = |c: &str, col, dec| Same {
            tools: vec!["gloci".into(), "bedops".into()],
            ..group(c, col, dec)
        };
        assert!(e(both("fuzzy", None, None)).contains("unknown compare mode"));
        assert!(e(both("numeric", Some(3), None)).contains("`columns` applies"));
        assert!(e(both("exact", None, Some(6))).contains("`decimals` applies"));
        assert!(e(both("exact", Some(0), None)).contains("at least 1"));
        let why: BTreeMap<String, String> = [("bedops".to_string(), " ".to_string())].into();
        assert!(
            check_schema("c", &[], &why, &cmds)
                .unwrap_err()
                .to_string()
                .contains("needs a reason")
        );
        let why: BTreeMap<String, String> = [("bedops".to_string(), "merged".to_string())].into();
        assert!(
            check_schema("c", &[both("exact", None, None)], &why, &cmds)
                .unwrap_err()
                .to_string()
                .contains("in a `same` group and in `not_compared`")
        );
    }

    #[test]
    fn unaccounted_names_gloci_commands_and_baselines() {
        let cmds = commands(&["gloci", "gloci-t8", "gloci-py", "bedtools", "bedops", "gia"]);
        let targets = vec![
            "gloci <= 0.5 * bedtools".to_string(),
            "gloci-py < pybedtools".into(),
        ];
        let none = BTreeMap::new();
        // gia is a yardstick (no target names it) and gloci-py is the Python tests' business.
        assert_eq!(
            unaccounted(&[], &none, &cmds, &targets),
            ["bedtools", "gloci", "gloci-t8"]
        );
        let g = Same {
            tools: vec!["gloci".into(), "bedtools".into()],
            ..group("exact", None, None)
        };
        let why: BTreeMap<String, String> =
            [("gloci-t8".to_string(), "the same command".to_string())].into();
        assert!(unaccounted(&[g], &why, &cmds, &targets).is_empty());
    }

    #[test]
    fn a_group_reports_what_ran_and_what_differs() {
        let d = tempfile::tempdir().unwrap();
        for (n, text) in [("gloci", "a\n"), ("bedtools", "a\n"), ("bedops", "b\n")] {
            std::fs::write(d.path().join(n), text).unwrap();
        }
        let g = Same {
            tools: vec![
                "gloci".into(),
                "bedtools".into(),
                "bedops".into(),
                "gia".into(),
            ],
            ..group("exact", None, None)
        };
        let out = |t: &str| Ok((t != "gia").then(|| Output::File(d.path().join(t))));
        let r = compare_group(&g, &out).unwrap();
        assert!(!r.ok());
        assert_eq!(r.missing, ["gia"]);
        assert_eq!(r.different.len(), 1);
        assert_eq!(r.different[0].0, "bedops");
        let line = r.line();
        assert!(
            line.starts_with("DIFFERENT (exact): gloci, bedtools, bedops [not installed: gia]"),
            "{line}"
        );
        assert!(
            line.contains("gloci against bedops: line 1 differs"),
            "{line}"
        );
        // Only one tool ran: nothing was compared.
        let alone = |t: &str| Ok((t == "gloci").then(|| Output::File(d.path().join(t))));
        let r = compare_group(&g, &alone).unwrap();
        assert!(
            r.ok() && r.line().starts_with("not checked (exact)"),
            "{}",
            r.line()
        );
    }

    fn sh(script: &str) -> Output {
        let mut c = Command::new("sh");
        c.args(["-c", script]);
        Output::Records(c)
    }

    #[test]
    fn decoded_records_are_read_from_the_decoder() {
        let g = Same {
            tools: vec!["gloci".into(), "samtools".into()],
            ..group("bam-records", None, None)
        };
        let print = |a: &'static str, b: &'static str| {
            move |t: &str| Ok(Some(sh(if t == "gloci" { a } else { b })))
        };
        let r = compare_group(&g, &print("printf 'r1\\nr2\\n'", "printf 'r1\\nr2\\n'")).unwrap();
        assert!(r.ok(), "{}", r.line());
        let r = compare_group(&g, &print("printf 'r1\\nr2\\n'", "printf 'r1\\nr3\\n'")).unwrap();
        assert!(r.line().contains("line 2 differs"), "{}", r.line());
        // More than a pipe holds after the difference: the decoders still end.
        let long = "printf 'x\\n'; yes r | head -n 200000";
        let r = compare_group(&g, &print("printf 'y\\n'; yes r | head -n 200000", long)).unwrap();
        assert!(r.line().contains("line 1 differs"), "{}", r.line());
    }

    #[test]
    fn a_decoder_that_fails_is_an_error_not_agreement() {
        // Two truncated files that print the same records before the decoder gives up.
        let g = Same {
            tools: vec!["gloci".into(), "samtools".into()],
            ..group("bam-records", None, None)
        };
        let out = |_: &str| Ok(Some(sh("printf 'r1\\n'; exit 3")));
        let e = compare_group(&g, &out).unwrap_err().to_string();
        assert!(e.contains("cannot read an output to compare it"), "{e}");
    }
}
