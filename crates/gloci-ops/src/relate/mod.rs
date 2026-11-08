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

pub use criteria::accepted;
pub use out::{ReportOutput, open_report};
pub use rewrite::{b_placeholder, rewrite_coords};

use gloci_core::{ChromOrder, ChromTable, OverlapCriteria, StrandRule, Window};
use gloci_engine::{
    Inputs, NearestConfig, SweepConfig, SweepSummary, Visitor, sweep_sorted, sweep_unsorted,
};
use gloci_io::{
    Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordReader, chrom_table_for,
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

/// Run the sweep over `inputs` and call `v` for every A record with its
/// candidate hits (D36). The chromosome table is built per D51. With `-s`
/// or `-S`, warns once per input that has records of unknown strand (D54);
/// the sorted sweep checks only the records it reads.
pub fn run_relation(
    inputs: RelateInputs,
    o: &RelateOptions,
    nearest: Option<NearestConfig>,
    v: &mut dyn Visitor,
) -> Result<SweepSummary, IoError> {
    let RelateInputs { a, b } = inputs;
    let mut chroms = {
        let all: Vec<&dyn RecordReader> = std::iter::once(a.as_ref())
            .chain(b.iter().map(|r| r.as_ref()))
            .collect();
        chrom_table_for(o.genome.as_deref(), &all, o.chrom_order)?
    };
    let watching = o.criteria.strand != StrandRule::Any;
    let mut watched: Vec<(String, Arc<AtomicBool>)> = Vec::new();
    let mut wrap = |r: Box<dyn RecordReader>| -> Box<dyn RecordReader> {
        if !watching {
            return r;
        }
        let flag = Arc::new(AtomicBool::new(false));
        watched.push((r.name().to_string(), flag.clone()));
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
        sweep_sorted(inputs, &mut chroms, &cfg, v)?
    };
    for (name, flag) in &watched {
        if flag.load(Ordering::Relaxed) {
            log::warn!("{name} has records with unknown strand; they never match under -s or -S");
        }
    }
    Ok(summary)
}
