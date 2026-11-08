//! The relation scaffold shared by `intersect`, `merge`, `subtract`,
//! `closest`, `map` and `coverage`: it runs the engine over streams
//! (`run_relation`), decides which candidate hits count (`criteria`), and
//! writes results (`out`, `rewrite`).
//!
//! Operations take readers, never paths (spec §8.1). A path-aware driver
//! for per-chromosome parallelism comes later and is a separate path.

mod criteria;
mod out;
mod rewrite;

pub use criteria::Acceptor;
pub use out::{ReportOutput, open_report};
pub use rewrite::{b_placeholder, rewrite_coords};

use crate::sort::{SortSettings, sort_readers};
use gloci_core::{ChromOrder, ChromTable, OverlapCriteria, StrandRule, Window};
use gloci_engine::{
    Inputs, NearestConfig, SweepConfig, SweepSummary, Visitor, sweep_sorted, sweep_unsorted,
};
use gloci_io::{
    Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordReader, UsageError,
    chrom_table_for,
};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

#[derive(Debug, Clone)]
pub struct RelateOptions {
    pub criteria: OverlapCriteria,
    /// Padding of B (`-w`, `--left`, `--right`, `--strand-window`). Cannot
    /// be combined with thresholds (D37); the command layer enforces that.
    pub window: Option<Window>,
    /// `--split`: compare BED12 blocks instead of spans.
    pub split: bool,
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
    if o.window.is_some() && has_threshold(&o.criteria) {
        return Err(UsageError(
            "-w/--left/--right cannot be combined with -f, -F, -r, -e or --min-bases (D37)".into(),
        )
        .into());
    }
    if o.sort.is_some() && o.unsorted {
        return Err(UsageError("--sort and --unsorted cannot be combined".into()).into());
    }
    let RelateInputs { a, b } = inputs;
    let mut chroms = {
        let all: Vec<&dyn RecordReader> = std::iter::once(a.as_ref())
            .chain(b.iter().map(|r| r.as_ref()))
            .collect();
        chrom_table_for(o.genome.as_deref(), &all, o.chrom_order)?
    };
    // Messages that are not about a record of the sorted stream name the
    // user's files, not the temporary sorted ones.
    let names: Vec<String> = std::iter::once(a.name())
        .chain(b.iter().map(|r| r.name()))
        .map(str::to_string)
        .collect();
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
    let watching = o.criteria.strand != StrandRule::Any;
    let mut watched: Vec<(String, Arc<AtomicBool>)> = Vec::new();
    let mut names = names.into_iter();
    let mut wrap = |r: Box<dyn RecordReader>| -> Box<dyn RecordReader> {
        let name = names.next().unwrap_or_default();
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
        a: wrap(a),
        b: b.into_iter().map(&mut wrap).collect(),
    };
    let cfg = SweepConfig {
        window: o.window,
        inclusive_points: o.criteria.inclusive_points,
        nearest,
    };
    let summary = if o.unsorted {
        sweep_unsorted(inputs, &mut chroms, &cfg, v)?
    } else {
        sweep_sorted(inputs, &mut chroms, &cfg, v).map_err(sorted_loc)?
    };
    for (name, flag) in &watched {
        if flag.load(Ordering::Relaxed) {
            log::warn!("{name} has records with unknown strand; they never match under -s or -S");
        }
    }
    Ok(summary)
}
