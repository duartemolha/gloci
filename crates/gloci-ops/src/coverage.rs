//! `coverage`: how much of each A record the B records cover
//! (bedtools coverage and annotate).
//!
//! A reporting layer over the relation scaffold (`relate`). Every accepted
//! hit (D36) is clipped to A and adds 1 to the depth over the bases it
//! shares with A. With `--split`, A and B are their BED12 blocks, depth is
//! counted over the bases shared by an A block and a B block, and A's
//! length is the sum of its blocks. A point A (D56) has no bases: it
//! reports the number of B overlapping it, 0 bases, length 0 and fraction
//! `.`; `--mean` prints `.` and `--hist`/`--per-base` print no lines.
//!
//! One TSV row per A is written for the summary and `--mean` reports, and
//! one row per depth (`--hist`) or per base (`--per-base`). Fractions and
//! means print in Rust's shortest round-trip form (D46). With several `-b`
//! files only the summary report is allowed, and it gives one covered
//! fraction per file (D57).

use crate::intersect::Text;
use crate::relate::{
    Acceptor, RelateInputs, RelateOptions, RelationFiles, RunOptions, run_files, run_relation_in,
};
use gloci_core::ChromTable;
use gloci_core::Interval;
use gloci_engine::{ARec, Hit, Neighbours, SweepSummary, Visitor};
use gloci_io::{IoError, RecordWriter, UsageError};
use std::collections::BTreeMap;
use std::io::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageReport {
    /// Hits, covered bases, A length, covered fraction.
    Summary,
    /// One row per depth present in A, with `all` rows after the last A
    /// (`write_hist_totals`).
    Hist,
    /// One row per base of A: its 1-based offset and depth.
    PerBase,
    /// Mean depth over A.
    Mean,
}

#[derive(Debug, Clone)]
pub struct CoverageOptions {
    pub report: CoverageReport,
    /// `window` must be `None`: coverage has no window flags.
    pub relate: RelateOptions,
}

/// Bases per depth, and the total length, over the non-point A records
/// reported with `--hist`. Chunks are summed with `merge`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CoverageTotals {
    depth_bases: BTreeMap<u64, u64>,
    length: u64,
}

impl CoverageTotals {
    pub fn merge(&mut self, other: &CoverageTotals) {
        for (d, b) in &other.depth_bases {
            *self.depth_bases.entry(*d).or_default() += b;
        }
        self.length += other.length;
    }

    /// Bases of A at `depth`, over all A.
    pub fn bases_at(&self, depth: u64) -> u64 {
        self.depth_bases.get(&depth).copied().unwrap_or(0)
    }

    /// The summed length of all A.
    pub fn length(&self) -> u64 {
        self.length
    }
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Option checks that need only the number of `-b` files.
pub fn validate(o: &CoverageOptions, n_b: usize) -> Result<(), IoError> {
    if n_b == 0 {
        return Err(usage("coverage needs at least one -b file"));
    }
    if n_b > 1 && o.report != CoverageReport::Summary {
        return Err(usage(
            "with several -b files coverage writes one covered fraction per file; \
--hist, --per-base and --mean take a single -b (run coverage once per file)",
        ));
    }
    if o.relate.window.is_some() {
        return Err(usage("coverage has no window options"));
    }
    Ok(())
}

/// The names of the appended columns, for a VCF A's column-name line.
pub fn column_names(o: &CoverageOptions, n_b: usize) -> Vec<String> {
    let strs = |v: &[&str]| v.iter().map(|s| (*s).to_string()).collect();
    match o.report {
        CoverageReport::Summary if n_b > 1 => (0..n_b)
            .map(|i| match o.relate.labels.get(i) {
                Some(l) => format!("covered_fraction_{l}"),
                None => format!("covered_fraction_b{}", i + 1),
            })
            .collect(),
        CoverageReport::Summary => strs(&["hits", "covered_bases", "a_length", "covered_fraction"]),
        CoverageReport::Hist => strs(&["depth", "bases", "a_length", "fraction"]),
        CoverageReport::PerBase => strs(&["position", "depth"]),
        CoverageReport::Mean => strs(&["mean_depth"]),
    }
}

/// Write the `all` rows of `--hist`: depth, bases, total length, fraction.
pub fn write_hist_totals(t: &CoverageTotals, w: &mut dyn RecordWriter) -> Result<(), IoError> {
    if t.length == 0 {
        return Ok(());
    }
    let mut row = Vec::new();
    for (d, b) in &t.depth_bases {
        row.clear();
        write!(row, "all\t{d}\t{b}\t{}\t{}", t.length, ratio(*b, t.length)).expect("vec write");
        write_row(w, &row)?;
    }
    Ok(())
}

fn write_row(w: &mut dyn RecordWriter, row: &[u8]) -> Result<(), IoError> {
    w.write_record(row).map_err(|source| IoError::Write {
        path: "output".to_string(),
        source,
    })
}

fn ratio(num: u64, den: u64) -> String {
    (num as f64 / den as f64).to_string()
}

/// Segments `(start, end, depth)` with depth above 0, in order, from the
/// pieces `(start, end)` (each adds 1 over its bases).
fn profile(pieces: &[(u64, u64)]) -> Vec<(u64, u64, u64)> {
    let mut edges: Vec<(u64, i64)> = Vec::with_capacity(pieces.len() * 2);
    for &(s, e) in pieces {
        edges.push((s, 1));
        edges.push((e, -1));
    }
    edges.sort_unstable();
    let mut out: Vec<(u64, u64, u64)> = Vec::new();
    let (mut depth, mut at) = (0i64, 0u64);
    let mut i = 0;
    while i < edges.len() {
        let pos = edges[i].0;
        if depth > 0 && pos > at {
            out.push((at, pos, depth as u64));
        }
        while i < edges.len() && edges[i].0 == pos {
            depth += edges[i].1;
            i += 1;
        }
        at = pos;
    }
    out
}

fn union_len(pieces: &[(u64, u64)]) -> u64 {
    profile(pieces).iter().map(|(s, e, _)| e - s).sum()
}

struct Reporter<'o, 'w> {
    report: CoverageReport,
    split: bool,
    n_b: usize,
    acceptor: Acceptor<'o>,
    w: &'w mut dyn RecordWriter,
    a: Text,
    row: Vec<u8>,
    totals: CoverageTotals,
}

impl Reporter<'_, '_> {
    fn emit(&mut self) -> Result<(), IoError> {
        write_row(self.w, &self.row)
    }
}

impl Visitor for Reporter<'_, '_> {
    fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], _: &Neighbours<'_>) -> Result<(), IoError> {
        let a_blocks: Vec<Interval> = if self.split {
            self.acceptor.a_blocks(a)?.to_vec()
        } else {
            vec![a.meta.span]
        };
        let length: u64 = a_blocks.iter().map(Interval::len).sum();
        let mut count = 0u64;
        // Pieces of A covered by an accepted hit, per file.
        let mut pieces: Vec<Vec<(u64, u64)>> = vec![Vec::new(); self.n_b];
        for h in hits {
            if !self.acceptor.accepted(a, h)? {
                continue;
            }
            count += 1;
            if length == 0 {
                continue;
            }
            let b_blocks = if self.split {
                self.acceptor.hit_blocks(h)?
            } else {
                vec![h.span]
            };
            let dst = &mut pieces[usize::from(h.file)];
            for ab in &a_blocks {
                for bb in &b_blocks {
                    let (s, e) = (ab.start.max(bb.start), ab.end.min(bb.end));
                    if s < e {
                        dst.push((s, e));
                    }
                }
            }
        }
        let line = self.a.of(a.raw, a.meta.loc)?;
        let mut row = std::mem::take(&mut self.row);
        row.clear();
        row.extend_from_slice(line);
        self.row = row;
        let head = self.row.len();
        match self.report {
            CoverageReport::Summary if self.n_b > 1 => {
                for p in &pieces {
                    self.row.push(b'\t');
                    if length == 0 {
                        self.row.push(b'.');
                    } else {
                        let r = ratio(union_len(p), length);
                        self.row.extend_from_slice(r.as_bytes());
                    }
                }
                self.emit()
            }
            CoverageReport::Summary => {
                let covered = union_len(&pieces[0]);
                let frac = if length == 0 {
                    ".".to_string()
                } else {
                    ratio(covered, length)
                };
                write!(self.row, "\t{count}\t{covered}\t{length}\t{frac}").expect("vec write");
                self.emit()
            }
            CoverageReport::Mean => {
                if length == 0 {
                    self.row.extend_from_slice(b"\t.");
                } else {
                    let total: u64 = pieces[0].iter().map(|(s, e)| e - s).sum();
                    write!(self.row, "\t{}", ratio(total, length)).expect("vec write");
                }
                self.emit()
            }
            CoverageReport::Hist => {
                if length == 0 {
                    return Ok(());
                }
                let mut by_depth: BTreeMap<u64, u64> = BTreeMap::new();
                let mut covered = 0;
                for (s, e, d) in profile(&pieces[0]) {
                    *by_depth.entry(d).or_default() += e - s;
                    covered += e - s;
                }
                if covered < length {
                    by_depth.insert(0, length - covered);
                }
                for (d, b) in by_depth {
                    self.row.truncate(head);
                    write!(self.row, "\t{d}\t{b}\t{length}\t{}", ratio(b, length))
                        .expect("vec write");
                    self.emit()?;
                    *self.totals.depth_bases.entry(d).or_default() += b;
                }
                self.totals.length += length;
                Ok(())
            }
            CoverageReport::PerBase => {
                let segs = profile(&pieces[0]);
                let (mut offset, mut k) = (0u64, 0usize);
                for blk in &a_blocks {
                    for pos in blk.start..blk.end {
                        while k < segs.len() && segs[k].1 <= pos {
                            k += 1;
                        }
                        let depth = match segs.get(k) {
                            Some(&(s, _, d)) if s <= pos => d,
                            _ => 0,
                        };
                        offset += 1;
                        self.row.truncate(head);
                        write!(self.row, "\t{offset}\t{depth}").expect("vec write");
                        self.emit()?;
                    }
                }
                Ok(())
            }
        }
    }
}

/// Run `coverage` over `inputs`, writing TSV rows to `w` (`open_report`
/// with `ReportOutput::Tsv`, columns from `column_names`), and with
/// `--hist` the `all` rows last. Nothing is finished here; the caller calls
/// `finish`.
pub fn coverage(
    inputs: RelateInputs,
    o: &CoverageOptions,
    w: &mut dyn RecordWriter,
) -> Result<(SweepSummary, CoverageTotals), IoError> {
    validate(o, inputs.b.len())?;
    let (summary, totals) = coverage_in(inputs, None, o, w)?;
    if o.report == CoverageReport::Hist {
        write_hist_totals(&totals, w)?;
    }
    Ok((summary, totals))
}

/// Run `coverage` over the files a command opened: per chromosome when
/// every input is indexed (D41), else streamed. Writes as [`coverage`],
/// the `--hist` totals summed over every chromosome.
pub fn coverage_files(
    files: RelationFiles,
    o: &CoverageOptions,
    w: &mut dyn RecordWriter,
    run: &RunOptions,
) -> Result<(SweepSummary, CoverageTotals), IoError> {
    validate(o, files.b_count())?;
    let (summary, parts) = run_files(files, &o.relate, run, w, &|inputs, chroms, w| {
        coverage_in(inputs, chroms, o, w)
    })?;
    let mut totals = CoverageTotals::default();
    for t in &parts {
        totals.merge(t);
    }
    if o.report == CoverageReport::Hist {
        write_hist_totals(&totals, w)?;
    }
    Ok((summary, totals))
}

fn coverage_in(
    inputs: RelateInputs,
    chroms: Option<ChromTable>,
    o: &CoverageOptions,
    w: &mut dyn RecordWriter,
) -> Result<(SweepSummary, CoverageTotals), IoError> {
    let n_b = inputs.b.len();
    let a = Text::new(inputs.a.as_ref())?;
    let acceptor = Acceptor::new(&o.relate, &inputs);
    let mut reporter = Reporter {
        report: o.report,
        split: o.relate.split,
        n_b,
        acceptor,
        w,
        a,
        row: Vec::new(),
        totals: CoverageTotals::default(),
    };
    let summary = run_relation_in(inputs, chroms, &o.relate, None, &mut reporter)?;
    Ok((summary, reporter.totals))
}
