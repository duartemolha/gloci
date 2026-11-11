//! Commands over BEDPE pairs (D146): `pair-intersect` (D154) and
//! `pair-pair` (D155). Each end of a pair is queried against an in-memory
//! [`IntervalIndex`](gloci_engine::IntervalIndex).
//!
//! A streams in batches, each batch is worked on by one of `threads`
//! workers, and the finished batches are written in A's order, so the output
//! does not depend on the thread count. At most `2 * threads` batches are in
//! flight. That driver is [`run_a`], shared by both commands.

mod intersect;
mod pair_pair;

pub use intersect::{PairIntersectOptions, PairType, pair_intersect, validate};
pub use pair_pair::{PairPairOptions, PairPairType, pair_pair};

use gloci_core::ChromTable;
use gloci_io::{BedpeBatch, BedpeReader, IoError};
use std::collections::BTreeMap;
use std::io::Write;
use std::num::NonZeroUsize;
use std::sync::mpsc::{channel, sync_channel};
use std::sync::{Arc, Mutex};

/// What a run saw.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PairSummary {
    /// Pairs read from A.
    pub a_pairs: u64,
    /// Pairs a span type left out because an end is absent (one or both) or
    /// the ends are on different chromosomes.
    pub span_skipped: u64,
}

/// What one batch saw.
#[derive(Default)]
pub(crate) struct Counts {
    pub pairs: u64,
    pub span_skipped: u64,
    /// An end of unknown strand met a strand rule.
    pub a_unknown: bool,
}

/// The per-batch work of a pair command: rows for a batch of A's pairs.
pub(crate) trait PairWorker: Sync {
    fn batch(&self, batch: &BedpeBatch, out: &mut Vec<u8>) -> Result<Counts, IoError>;
}

/// The result of one batch, by sequence number.
struct Done {
    out: Vec<u8>,
    counts: Result<Counts, IoError>,
}

fn werr(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".to_string(),
        source,
    }
}

/// Read all of A in batches, run `worker` over them (in parallel when
/// `threads > 1`) and write the rows in A's order. Returns the summary and
/// whether an end of unknown strand met a strand rule. `label` names the
/// command in internal errors.
pub(crate) fn run_a(
    a: &mut BedpeReader,
    chroms: &mut ChromTable,
    worker: &dyn PairWorker,
    threads: NonZeroUsize,
    out: &mut dyn Write,
    label: &str,
) -> Result<(PairSummary, bool), IoError> {
    let mut summary = PairSummary::default();
    let mut a_unknown = false;
    let mut take = |c: &Counts| {
        summary.a_pairs += c.pairs;
        summary.span_skipped += c.span_skipped;
        a_unknown |= c.a_unknown;
    };
    let mut max = 1usize;
    let mut next = |a: &mut BedpeReader,
                    chroms: &mut ChromTable,
                    batch: &mut BedpeBatch|
     -> Result<bool, IoError> {
        batch.clear();
        let more = a.read_batch(chroms, batch, max)?;
        if !batch.is_empty() {
            let per_rec = batch.buf.len() / batch.len() + 1;
            max = ((4usize << 20) / per_rec).clamp(1, 16384);
        }
        Ok(more)
    };
    if threads.get() == 1 {
        let (mut batch, mut buf) = (BedpeBatch::default(), Vec::new());
        loop {
            let more = next(a, chroms, &mut batch)?;
            buf.clear();
            take(&worker.batch(&batch, &mut buf)?);
            out.write_all(&buf).map_err(werr)?;
            if !more {
                break;
            }
        }
    } else {
        parallel(a, chroms, worker, threads, out, &mut next, &mut take, label)?;
    }
    if let Some(s) = a.diagnostics().summary(a.name()) {
        log::warn!("{s}");
    }
    Ok((summary, a_unknown))
}

type NextBatch<'f> =
    dyn FnMut(&mut BedpeReader, &mut ChromTable, &mut BedpeBatch) -> Result<bool, IoError> + 'f;

/// Batches of A across workers, written in order; at most `2 * threads`
/// batches in flight.
#[allow(clippy::too_many_arguments)]
fn parallel(
    a: &mut BedpeReader,
    chroms: &mut ChromTable,
    worker: &dyn PairWorker,
    threads: NonZeroUsize,
    out: &mut dyn Write,
    next: &mut NextBatch<'_>,
    take: &mut dyn FnMut(&Counts),
    label: &str,
) -> Result<(), IoError> {
    let n = threads.get();
    let limit = 2 * n;
    let died = || IoError::Internal(format!("a {label} worker died"));
    std::thread::scope(|scope| {
        let (job_tx, job_rx) = sync_channel::<(u64, BedpeBatch)>(limit);
        let job_rx = Arc::new(Mutex::new(job_rx));
        let (done_tx, done_rx) = channel::<(u64, Done)>();
        for _ in 0..n {
            let (job_rx, done_tx) = (job_rx.clone(), done_tx.clone());
            scope.spawn(move || {
                loop {
                    let job = job_rx.lock().unwrap_or_else(|e| e.into_inner()).recv();
                    let Ok((seq, batch)) = job else { break };
                    let mut buf = Vec::new();
                    let counts = worker.batch(&batch, &mut buf);
                    if done_tx.send((seq, Done { out: buf, counts })).is_err() {
                        break;
                    }
                }
            });
        }
        drop(done_tx);
        let mut pending: BTreeMap<u64, Done> = BTreeMap::new();
        let (mut sent, mut written) = (0u64, 0u64);
        let mut flush = |pending: &mut BTreeMap<u64, Done>, written: &mut u64| {
            while let Some(d) = pending.remove(written) {
                take(&d.counts?);
                out.write_all(&d.out).map_err(werr)?;
                *written += 1;
            }
            Ok::<(), IoError>(())
        };
        let mut result = Ok(());
        loop {
            let mut batch = BedpeBatch::default();
            let more = match next(a, chroms, &mut batch) {
                Ok(m) => m,
                Err(e) => {
                    result = Err(e);
                    break;
                }
            };
            while sent - written >= limit as u64 && result.is_ok() {
                match done_rx.recv() {
                    Ok((seq, d)) => {
                        pending.insert(seq, d);
                        result = flush(&mut pending, &mut written);
                    }
                    Err(_) => result = Err(died()),
                }
            }
            if result.is_err() {
                break;
            }
            if !batch.is_empty() {
                if job_tx.send((sent, batch)).is_err() {
                    result = Err(died());
                    break;
                }
                sent += 1;
            }
            if !more {
                break;
            }
        }
        drop(job_tx);
        while result.is_ok() && written < sent {
            match done_rx.recv() {
                Ok((seq, d)) => {
                    pending.insert(seq, d);
                    result = flush(&mut pending, &mut written);
                }
                Err(_) => result = Err(died()),
            }
        }
        result
    })
}
