//! The relation scaffold shared by `intersect`, `merge`, `subtract`,
//! `closest`, `map` and `coverage`: it runs the engine over streams
//! (`run_relation`), decides which candidate hits count (`criteria`), and
//! writes results (`out`, `rewrite`).
//!
//! Operations take readers (spec §8.1). [`RelationFiles`] carries the
//! readers a command opened together with their paths, so that the one
//! entry point per command can run per chromosome on indexed inputs (D41)
//! and plain seekable files (D76) (`indexed`), or stream, with identical
//! output.

mod criteria;
mod files;
mod indexed;
mod out;
mod peek;
mod rewrite;

pub use criteria::Acceptor;
pub use files::{InputPaths, RelationFiles, RunOptions, runs_per_chromosome};
pub(crate) use files::{b_width, relation_chroms, run_files};

pub use out::{ReportOutput, open_report};
pub use rewrite::{b_placeholder, rewrite_coords};

use crate::skipped::{Skipped, log_skipped};
use crate::sort::{SortSettings, sort_readers};
use gloci_core::{ChromOrder, ChromTable, OverlapCriteria, StrandRule, Window};
use gloci_engine::{
    Inputs, NearestConfig, SweepConfig, SweepSummary, Visitor, sweep_sorted, sweep_unsorted,
};
use gloci_io::{
    Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordReader, UsageError,
};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

impl RelateOptions {
    /// Do the commands compare blocks rather than spans?
    pub fn splits(&self) -> bool {
        self.split || self.split_deletions
    }
}

#[derive(Debug, Clone)]
pub struct RelateOptions {
    pub criteria: OverlapCriteria,
    /// Padding of B (`-w`, `--left`, `--right`, `--strand-window`). Cannot
    /// be combined with thresholds (D37); the command layer enforces that.
    pub window: Option<Window>,
    /// `--split`: compare BED12 blocks and alignment CIGAR segments instead of spans.
    pub split: bool,
    /// `--split-deletions`: alignments are also cut at `D` (D105), with or without `--split`;
    /// BED12 blocks need `--split`.
    pub split_deletions: bool,
    pub unsorted: bool,
    /// `--sort`: sort A and then each B on its own (so the B files stay
    /// distinct) with the external sorter before the sorted sweep. The
    /// readers in [`RelateInputs`] must have been opened with
    /// [`crate::sort_input_options`]. Cannot be combined with `unsorted`.
    pub sort: Option<SortSettings>,
    /// One label per B input (`--names`), for commands that print them.
    pub labels: Vec<String>,
    pub genome: Option<String>,
    pub chrom_order: ChromOrder,
    /// Write header lines (`--header`).
    pub header: bool,
}

pub struct RelateInputs {
    pub a: Box<dyn RecordReader>,
    pub b: Vec<Box<dyn RecordReader>>,
}

thread_local! {
    /// Set while this thread runs one chunk of a per-chromosome run, whose
    /// readers the driver already watches for unknown strands (once per
    /// input, not once per chunk).
    static IN_CHUNK: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Marks the current thread as running a chunk of a per-chromosome run until
/// dropped.
struct ChunkScope;

impl ChunkScope {
    fn enter() -> ChunkScope {
        IN_CHUNK.with(|c| c.set(true));
        ChunkScope
    }
}

impl Drop for ChunkScope {
    fn drop(&mut self) {
        IN_CHUNK.with(|c| c.set(false));
    }
}

/// Passes records through and notes whether any had an unknown strand.
struct StrandWatch {
    inner: Box<dyn RecordReader>,
    unknown: Arc<AtomicBool>,
}

impl RecordReader for StrandWatch {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn compression(&self) -> Compression {
        self.inner.compression()
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        let before = batch.recs.len();
        let more = self.inner.read_batch(chroms, batch, max)?;
        if batch.recs[before..]
            .iter()
            .any(|m| m.strand == gloci_core::Strand::Unknown)
        {
            self.unknown.store(true, Ordering::Relaxed);
        }
        Ok(more)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn consume_to_end(&self) -> bool {
        self.inner.consume_to_end()
    }
}

/// Passes records through and records distinct chromosome names seen in first-seen order.
struct ChromWatch {
    inner: Box<dyn RecordReader>,
    last_chrom: Option<gloci_core::ChromId>,
    /// Ids already in `seen` (one table, so an id is one name): an A that
    /// comes back to a chromosome, or a file with many contigs, costs O(1).
    in_seen: gloci_io::naming::SeenChroms,
    seen: Arc<Mutex<Vec<Vec<u8>>>>,
}

impl RecordReader for ChromWatch {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn compression(&self) -> Compression {
        self.inner.compression()
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        let before = batch.recs.len();
        let more = self.inner.read_batch(chroms, batch, max)?;
        for m in &batch.recs[before..] {
            if self.last_chrom != Some(m.chrom) {
                self.last_chrom = Some(m.chrom);
                if self.in_seen.first_time(m.chrom) {
                    self.seen
                        .lock()
                        .unwrap()
                        .push(chroms.name(m.chrom).to_vec());
                }
            }
        }
        Ok(more)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn consume_to_end(&self) -> bool {
        self.inner.consume_to_end()
    }
}

/// Errors read from a sorted temporary file count data records, not lines.
fn sorted_loc(e: IoError) -> IoError {
    e.with_sorted_input_loc()
}

fn has_threshold(c: &OverlapCriteria) -> bool {
    c.min_frac_a.is_some()
        || c.min_frac_b.is_some()
        || c.reciprocal
        || c.either
        || c.min_bases.is_some()
}

/// Run the sweep over `inputs` and call `v` for every A record with its
/// candidate hits (D36). A window together with any threshold is a usage
/// error (D37). The chromosome table is built per D51. With `-s`
/// or `-S`, warns once per input that has records of unknown strand (D54);
/// the sorted sweep checks only the records it reads.
pub fn run_relation(
    inputs: RelateInputs,
    o: &RelateOptions,
    nearest: Option<NearestConfig>,
    v: &mut dyn Visitor,
) -> Result<SweepSummary, IoError> {
    let summary = run_relation_in(inputs, None, o, nearest, v)?;
    report_skipped(&summary);
    Ok(summary)
}

/// Warn about the records each input skipped under `--skip-invalid`
/// (`summary.skipped`): the count and the first 10 examples.
pub(crate) fn report_skipped(summary: &SweepSummary) {
    log_skipped(&summary.skipped);
}

/// [`run_relation`] with the chromosome table given when something was
/// already read with it (a record read ahead); `None` builds it (D51).
pub(crate) fn run_relation_in(
    inputs: RelateInputs,
    chroms: Option<ChromTable>,
    o: &RelateOptions,
    nearest: Option<NearestConfig>,
    v: &mut dyn Visitor,
) -> Result<SweepSummary, IoError> {
    if o.window.is_some() && has_threshold(&o.criteria) {
        return Err(UsageError(
            "-w/--left/--right cannot be combined with -f, -F, -r, -e or --min-bases (D37)".into(),
        )
        .into());
    }
    if o.sort.is_some() && o.unsorted {
        return Err(UsageError("--sort and --unsorted cannot be combined".into()).into());
    }
    let mut chroms = match chroms {
        Some(c) => c,
        None => relation_chroms(&inputs, o)?,
    };
    let RelateInputs { a, b } = inputs;
    // Messages that are not about a record of the sorted stream name the
    // user's files, not the temporary sorted ones.
    let names: Vec<String> = std::iter::once(a.name())
        .chain(b.iter().map(|r| r.name()))
        .map(str::to_string)
        .collect();
    let name_a = names[0].clone();
    let names_b: Vec<String> = names[1..].to_vec();
    let (a, b) = match &o.sort {
        None => (a, b),
        Some(settings) => {
            let a = sort_readers(vec![a], &mut chroms, settings).map_err(sorted_loc)?;
            let b = b
                .into_iter()
                .map(|r| sort_readers(vec![r], &mut chroms, settings))
                .collect::<Result<Vec<_>, _>>()
                .map_err(sorted_loc)?;
            (a, b)
        }
    };
    let in_chunk = IN_CHUNK.with(|c| c.get());
    let watching = o.criteria.strand != StrandRule::Any && !in_chunk;
    let mut watched: Vec<(String, Arc<AtomicBool>)> = Vec::new();
    let seen_a = Arc::new(Mutex::new(Vec::new()));
    let seen_bs: Vec<Arc<Mutex<Vec<Vec<u8>>>>> = (0..b.len())
        .map(|_| Arc::new(Mutex::new(Vec::new())))
        .collect();
    let mut orig_names = names.into_iter();
    // A `--sort`ed stream has nothing to report: `sort_readers` reported
    // the inputs it read.
    let mut skipped = Skipped::default();
    let mut wrap =
        |r: Box<dyn RecordReader>, seen: Arc<Mutex<Vec<Vec<u8>>>>| -> Box<dyn RecordReader> {
            let name = orig_names.next().unwrap_or_default();
            let r = skipped.watch_as(&name, r);
            let r: Box<dyn RecordReader> = Box::new(ChromWatch {
                inner: r,
                last_chrom: None,
                in_seen: Default::default(),
                seen,
            });
            if !watching {
                return r;
            }
            let flag = Arc::new(AtomicBool::new(false));
            watched.push((name, flag.clone()));
            Box::new(StrandWatch {
                inner: r,
                unknown: flag,
            })
        };
    let inputs = Inputs {
        a: wrap(a, seen_a.clone()),
        b: b.into_iter()
            .zip(&seen_bs)
            .map(|(r, seen)| wrap(r, seen.clone()))
            .collect(),
    };
    let cfg = SweepConfig {
        window: o.window,
        inclusive_points: o.criteria.inclusive_points,
        nearest,
    };
    let mut summary = if o.unsorted {
        sweep_unsorted(inputs, &mut chroms, &cfg, v)?
    } else {
        sweep_sorted(inputs, &mut chroms, &cfg, v).map_err(sorted_loc)?
    };
    // The sweep has dropped the readers.
    summary.skipped = skipped.take();
    skipped.log_alignment_counts();
    for (name, flag) in &watched {
        if flag.load(Ordering::Relaxed) {
            log::warn!("{name} has records with unknown strand; they never match under -s or -S");
        }
    }
    let a_chroms = seen_a.lock().unwrap().clone();
    // (A chunk of a per-chromosome run sees only its own chromosomes: the
    // run checks the whole input once.)
    if !a_chroms.is_empty() && !in_chunk {
        if names_b.len() == 1 {
            if !summary.shared_chroms {
                let b_chroms = seen_bs[0].lock().unwrap().clone();
                if !b_chroms.is_empty() {
                    log::warn!(
                        "{}",
                        gloci_io::naming::format_disjoint_chroms_warning(
                            &name_a,
                            &a_chroms,
                            &names_b[0],
                            &b_chroms
                        )
                    );
                }
            }
        } else {
            for (b_name, b_seen_arc) in names_b.iter().zip(&seen_bs) {
                let b_chroms = b_seen_arc.lock().unwrap().clone();
                if !b_chroms.is_empty() {
                    let b_set: std::collections::HashSet<&Vec<u8>> = b_chroms.iter().collect();
                    let shares = a_chroms.iter().any(|a| b_set.contains(a));
                    if !shares {
                        log::warn!(
                            "{}",
                            gloci_io::naming::format_disjoint_chroms_warning(
                                &name_a, &a_chroms, b_name, &b_chroms
                            )
                        );
                    }
                }
            }
        }
    }
    Ok(summary)
}
