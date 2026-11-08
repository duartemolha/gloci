//! The sorted sweep (spec §8.2). A and every B advance together in
//! `(chromosome, key)` order; each B input keeps an active set of the
//! records that can still overlap the current or a later A.
//!
//! - **Loading.** A B is loaded while `key <= a.span.end + pad`, where
//!   `pad = max(left, right)`. Since `span.start >= key`, that loads every
//!   B whose padded span can reach A.
//! - **Eviction.** A B leaves when `span.end + pad < a.key`. A's key, not
//!   its span start: a later A's span can start before this one's (a
//!   trimmed VCF deletion), but never before this key. A B that already
//!   meets this when read is never loaded.
//! - **Matching.** Hits are the active B where
//!   `overlaps(a.span, window.apply(b.span, a.strand), inclusive_points)`,
//!   by file, then B input order (D38).
//! - **Chromosomes.** When A moves to a new chromosome, the active sets
//!   are cleared and each B reads past records on earlier chromosomes;
//!   a record on a later one waits. Reading stops when A is exhausted.

use crate::active::ActiveSet;
use crate::hit::{ARec, Hit, Neighbours, Visitor};
use gloci_core::{ChromId, ChromTable, StrandRule, Window, overlaps, thousands};
use gloci_io::{
    IoError, OrderGuard, RecordBatch, RecordMeta, RecordReader, UsageError, next_batch_size,
};
use std::cmp::Ordering;

/// What `closest` asks of the sweep. Task 6 implements it; until then the
/// neighbours stay empty.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NearestConfig {
    pub k: usize,
    pub max_distance: Option<u64>,
    pub strand: StrandRule,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SweepConfig {
    /// Padding applied to every B span before testing (`-w`, `--left`,
    /// `--right`, `--strand-window`).
    pub window: Option<Window>,
    /// `--inclusive-points`: a point on an interval's edge overlaps it.
    pub inclusive_points: bool,
    /// Also report the nearest B on each side (`closest`).
    pub nearest: Option<NearestConfig>,
}

/// The streams of a relation: A, and the B inputs in command-line order.
pub struct Inputs {
    pub a: Box<dyn RecordReader>,
    pub b: Vec<Box<dyn RecordReader>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SweepSummary {
    pub a_records: u64,
    /// Records read from each B. Reading stops when A is exhausted, so
    /// this can be less than the file holds.
    pub b_records: Vec<u64>,
    /// At least one chromosome had both A and B records.
    pub shared_chroms: bool,
    /// The largest active set (all B inputs) any A was matched against,
    /// for tests: it is bounded by the overlap depth, not the file size.
    #[doc(hidden)]
    pub peak_active: usize,
}

/// Run the sweep: call `v` for every A record, in A's input order, with
/// its candidate hits (D36, D38). Every input is order-checked against
/// `chroms` (D51), and the first violation is an error naming that input.
pub fn sweep_sorted(
    inputs: Inputs,
    chroms: &mut ChromTable,
    cfg: &SweepConfig,
    v: &mut dyn Visitor,
) -> Result<SweepSummary, IoError> {
    let Inputs { a: mut reader, b } = inputs;
    let max_files = usize::from(u16::MAX) + 1;
    if b.len() > max_files {
        return Err(UsageError(format!(
            "too many B inputs: {} (at most {})",
            thousands(b.len() as u64),
            thousands(max_files as u64)
        ))
        .into());
    }
    let pad = cfg.window.map_or(0, |w| w.left.max(w.right));
    let mut streams: Vec<BStream> = b.into_iter().map(BStream::new).collect();
    let mut active: Vec<ActiveSet> = streams.iter().map(|_| ActiveSet::default()).collect();
    let mut guard = OrderGuard::new(reader.name(), reader.format());
    let mut summary = SweepSummary {
        a_records: 0,
        b_records: Vec::new(),
        shared_chroms: false,
        peak_active: 0,
    };
    // Ruling R-P1: `cfg.nearest` is accepted, and the neighbours stay
    // empty until Task 6.
    let near = Neighbours::default();
    let mut spare: Vec<Hit<'static>> = Vec::new();
    let mut chrom = None;
    let (mut batch, mut max) = (RecordBatch::default(), 1);
    loop {
        batch.clear();
        let more = reader.read_batch(chroms, &mut batch, max)?;
        guard.check(chroms, &batch)?;
        max = next_batch_size(&batch, max);
        summary.a_records += batch.len() as u64;
        for (i, a) in batch.recs.iter().enumerate() {
            if chrom != Some(a.chrom) {
                chrom = Some(a.chrom);
                for (s, set) in streams.iter_mut().zip(&mut active) {
                    set.clear();
                    summary.shared_chroms |= s.enter(chroms, a.chrom)?;
                }
            }
            let lim = a.span.end.saturating_add(pad);
            for (s, set) in streams.iter_mut().zip(&mut active) {
                set.compact();
                s.load(chroms, a, lim, pad, set)?;
            }
            let mut hits = recycle(std::mem::take(&mut spare));
            let mut live = 0;
            for (f, set) in active.iter_mut().enumerate() {
                let f = f as u16;
                let ip = cfg.inclusive_points;
                live += match cfg.window {
                    None => set.scan(f, lim, a.key, 0, |b| overlaps(a.span, b, ip), &mut hits),
                    Some(w) => set.scan(
                        f,
                        lim,
                        a.key,
                        pad,
                        |b| overlaps(a.span, w.apply(b, a.strand), ip),
                        &mut hits,
                    ),
                };
            }
            summary.peak_active = summary.peak_active.max(live);
            let rec = ARec {
                meta: *a,
                raw: batch.raw(i),
            };
            v.visit(&rec, &hits, &near)?;
            spare = recycle(hits);
        }
        if !more {
            break;
        }
    }
    summary.b_records = streams.iter().map(|s| s.read).collect();
    Ok(summary)
}

/// Reuse the allocation of a hit vector for hits that borrow something
/// else: the vector is emptied, and collecting an empty `into_iter` of a
/// same-layout type reuses its buffer.
fn recycle<'x, 'y>(mut v: Vec<Hit<'x>>) -> Vec<Hit<'y>> {
    v.clear();
    v.into_iter()
        .map(|_| unreachable!("the vector is empty"))
        .collect()
}

/// Where a B stream's next record is, relative to A's chromosome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Place {
    /// On A's chromosome.
    Same,
    /// On a later chromosome, or the stream is exhausted.
    Past,
}

/// One B input: its reader, the batch being consumed and its order guard.
struct BStream {
    reader: Box<dyn RecordReader>,
    guard: OrderGuard,
    batch: RecordBatch,
    /// The next record of `batch` to consume.
    pos: usize,
    max: usize,
    more: bool,
    place: Place,
    read: u64,
}

impl BStream {
    fn new(reader: Box<dyn RecordReader>) -> BStream {
        BStream {
            guard: OrderGuard::new(reader.name(), reader.format()),
            reader,
            batch: RecordBatch::default(),
            pos: 0,
            max: 1,
            more: true,
            place: Place::Past,
            read: 0,
        }
    }

    /// The index in `batch` of the next record, or `None` at the end.
    #[inline]
    fn peek(&mut self, chroms: &mut ChromTable) -> Result<Option<usize>, IoError> {
        if self.pos < self.batch.recs.len() {
            return Ok(Some(self.pos));
        }
        self.refill(chroms)
    }

    /// Read (and order-check) batches until one has a record.
    #[cold]
    fn refill(&mut self, chroms: &mut ChromTable) -> Result<Option<usize>, IoError> {
        while self.pos == self.batch.recs.len() {
            if !self.more {
                return Ok(None);
            }
            self.batch.clear();
            self.pos = 0;
            self.more = self.reader.read_batch(chroms, &mut self.batch, self.max)?;
            self.guard.check(chroms, &self.batch)?;
            self.max = next_batch_size(&self.batch, self.max);
            self.read += self.batch.len() as u64;
        }
        Ok(Some(self.pos))
    }

    /// A moved to chromosome `c`: read past the records on chromosomes
    /// ranked before it. Returns whether the next record is on `c`.
    fn enter(&mut self, chroms: &mut ChromTable, c: ChromId) -> Result<bool, IoError> {
        // The chromosome last compared with `c`, and whether it is before.
        let mut last: Option<(ChromId, bool)> = None;
        while let Some(i) = self.peek(chroms)? {
            let b = self.batch.recs[i].chrom;
            if b == c {
                self.place = Place::Same;
                return Ok(true);
            }
            let before = match last {
                Some((id, before)) if id == b => before,
                _ => {
                    let before = chroms.cmp(b, c) == Ordering::Less;
                    last = Some((b, before));
                    before
                }
            };
            if !before {
                break;
            }
            self.pos += 1;
        }
        self.place = Place::Past;
        Ok(false)
    }

    /// Load into `set` the records on A's chromosome with `key <= lim`,
    /// except those already ending before A's key (with the pad). A record
    /// on another chromosome is a later one: the order guard has checked it.
    fn load(
        &mut self,
        chroms: &mut ChromTable,
        a: &RecordMeta,
        lim: u64,
        pad: u64,
        set: &mut ActiveSet,
    ) -> Result<(), IoError> {
        if self.place != Place::Same {
            return Ok(());
        }
        while let Some(i) = self.peek(chroms)? {
            let m = &self.batch.recs[i];
            if m.chrom != a.chrom {
                self.place = Place::Past;
                return Ok(());
            }
            if m.key > lim {
                return Ok(());
            }
            if m.span.end.saturating_add(pad) >= a.key {
                set.push(m, self.batch.raw(i));
            }
            self.pos += 1;
        }
        self.place = Place::Past;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recycle_keeps_the_allocation() {
        let v: Vec<Hit<'static>> = Vec::with_capacity(64);
        let w: Vec<Hit<'_>> = recycle(v);
        assert!(w.is_empty());
        assert!(w.capacity() >= 64);
    }
}
