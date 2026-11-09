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
//! - **Neighbours** (`closest`). Evicted and skipped B go to a left pool
//!   per B input, and B are read ahead of the load limit for the right
//!   side; see `nearest`.

use crate::active::ActiveSet;
use crate::hit::{ARec, Hit, Neighbours, Visitor, recycle};
use crate::nearest::{LeftPool, LookAhead, NearScan, Target, Want, keep_nearest};
use gloci_core::{ChromId, ChromTable, StrandRule, Window, overlaps, thousands};
use gloci_io::{
    IoError, OrderGuard, RecordBatch, RecordMeta, RecordReader, UsageError, next_batch_size,
};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// What `closest` asks of the sweep: the `k` nearest B on each side of
/// each A, per B file ([`Neighbours`]). `k = 0` asks for none.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NearestConfig {
    pub k: usize,
    /// Only B at most this far from A (spec §4.4 distance).
    pub max_distance: Option<u64>,
    /// Only B on the same (or opposite) strand as A (D55).
    pub strand: StrandRule,
    /// With `strand` `Any`: keep the `k` nearest (plus ties) of `-` B and,
    /// separately, of the rest (`+` and unknown, which D6 counts as `+`),
    /// and report both groups together. `closest --signed b` needs this:
    /// whether a B is upstream there depends on its own strand. Under a
    /// strand rule only one B strand counts, so it changes nothing.
    pub split_by_b_strand: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SweepConfig {
    /// Padding applied to every B span before testing (`-w`, `--left`,
    /// `--right`, `--strand-window`).
    pub window: Option<Window>,
    /// `--inclusive-points`: a point on an interval's edge overlaps it.
    pub inclusive_points: bool,
    /// Also report the nearest B on each side of each A (`closest`).
    pub nearest: Option<NearestConfig>,
}

/// The streams of a relation: A, and the B inputs in command-line order.
pub struct Inputs {
    pub a: Box<dyn RecordReader>,
    pub b: Vec<Box<dyn RecordReader>>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SweepSummary {
    pub a_records: u64,
    /// Records read from each B. The sorted sweep stops reading when A is
    /// exhausted, so this can be less than the file holds; unsorted mode
    /// reads every record.
    pub b_records: Vec<u64>,
    /// At least one chromosome had both A and B records.
    pub shared_chroms: bool,
    /// Records skipped under `--skip-invalid`: each input's name and
    /// diagnostics, A first. The engine leaves this empty; the relation
    /// layer fills it once the readers are done.
    pub skipped: Vec<(String, gloci_io::Diagnostics)>,
    #[doc(hidden)]
    pub stats: SweepStats,
}

impl SweepSummary {
    pub fn merge(&mut self, other: &SweepSummary) {
        self.a_records += other.a_records;
        if self.b_records.len() < other.b_records.len() {
            self.b_records.resize(other.b_records.len(), 0);
        }
        for (b, ob) in self.b_records.iter_mut().zip(&other.b_records) {
            *b += ob;
        }
        self.shared_chroms |= other.shared_chroms;
        if self.skipped.is_empty() {
            self.skipped = other.skipped.clone();
        } else {
            for ((_, d), (_, od)) in self.skipped.iter_mut().zip(&other.skipped) {
                d.absorb(od.clone());
            }
        }
        self.stats.peak_active = self.stats.peak_active.max(other.stats.peak_active);
        self.stats.peak_pool = self.stats.peak_pool.max(other.stats.peak_pool);
        self.stats.peak_pool_bytes = self.stats.peak_pool_bytes.max(other.stats.peak_pool_bytes);
        self.stats.lapper_indexes += other.stats.lapper_indexes;
        self.stats.seek_queries += other.stats.seek_queries;
    }
}

/// What the sweep did, for tests. The sorted sweep fills the peak sizes
/// (memory is bounded by the overlap depth and `k`, spec §8.2, not by the
/// file size). `sweep_unsorted` leaves those 0 and fills the rest.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SweepStats {
    /// The largest active set (all B inputs) any A was matched against.
    pub peak_active: usize,
    /// The most B the left pools (all B inputs) held after any A.
    pub peak_pool: usize,
    /// The most raw bytes each B input's left pool held at any moment,
    /// summed over the inputs: a bound on what they held at once.
    pub peak_pool_bytes: usize,
    /// Unsorted mode: the (B input, chromosome) indexes built (D61).
    pub lapper_indexes: usize,
    /// Unsorted mode: the queries answered through lapper's `seek` (made
    /// while A was still in order).
    pub seek_queries: u64,
}

/// Run the sweep: call `v` for every A record, in A's input order, with
/// its candidate hits (D36, D38) and, when `cfg.nearest` is set, its
/// neighbours. Every input is order-checked against `chroms` (D51), and
/// the first violation is an error naming that input.
pub fn sweep_sorted(
    inputs: Inputs,
    chroms: &mut ChromTable,
    cfg: &SweepConfig,
    v: &mut dyn Visitor,
) -> Result<SweepSummary, IoError> {
    let Inputs { a: mut reader, b } = inputs;
    check_b_count(b.len())?;
    let pad = cfg.window.map_or(0, |w| w.left.max(w.right));
    let mut streams: Vec<BStream> = b.into_iter().map(BStream::new).collect();
    let mut active: Vec<ActiveSet> = streams.iter().map(|_| ActiveSet::default()).collect();
    let nearest = cfg.nearest.filter(|n| n.k > 0);
    let mut pools: Vec<LeftPool> = match &nearest {
        Some(n) => streams.iter().map(|_| LeftPool::new(n)).collect(),
        None => Vec::new(),
    };
    let mut right_bounds = vec![0; streams.len()];
    let mut starts: [BinaryHeap<u64>; 2] = Default::default();
    let mut guard = OrderGuard::new(reader.name(), reader.format());
    let mut summary = SweepSummary {
        a_records: 0,
        b_records: Vec::new(),
        shared_chroms: false,
        skipped: Vec::new(),
        stats: SweepStats::default(),
    };
    let mut spare: Vec<Hit<'static>> = Vec::new();
    let mut spare_near: Neighbours<'static> = Neighbours::default();
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
                pools.iter_mut().for_each(LeftPool::clear);
            }
            let lim = a.span.end.saturating_add(pad);
            let target = nearest.map(|n| Target::new(&n, a.span, a.strand, cfg.inclusive_points));
            for (f, (s, set)) in streams.iter_mut().zip(&mut active).enumerate() {
                set.compact();
                let mut pool = pools.get_mut(f);
                if let Some(p) = pool.as_deref_mut() {
                    p.expire(a.key);
                }
                s.load(chroms, a, lim, pad, set, pool)?;
                right_bounds[f] = match target {
                    // Nothing past the load limit is looked at.
                    Some(t) if t.want != Want::Nothing => {
                        s.look_ahead(chroms, a, lim, LookAhead::new(t, &mut starts), set)?
                    }
                    _ => lim,
                };
            }
            let mut hits = recycle(std::mem::take(&mut spare));
            let mut near = std::mem::take(&mut spare_near).recycle();
            let mut live = 0;
            match target {
                None => {
                    for (f, set) in active.iter_mut().enumerate() {
                        live += scan(set, f as u16, a, lim, pad, cfg, &mut hits, None);
                    }
                }
                Some(t) => {
                    let mut pooled = 0;
                    for (f, (set, pool)) in active.iter_mut().zip(pools.iter_mut()).enumerate() {
                        let (left, right) = (near.left.len(), near.right.len());
                        let mut ns = NearScan {
                            t,
                            a_key: a.key,
                            right_bound: right_bounds[f],
                            pool: &mut *pool,
                            left: &mut near.left,
                            right: &mut near.right,
                        };
                        let f = f as u16;
                        live += scan(set, f, a, lim, pad, cfg, &mut hits, Some(&mut ns));
                        // The pool is done changing for this A; its members
                        // are borrowed until the visit ends.
                        let pool: &LeftPool = pool;
                        pool.extend(f, &t, &mut near.left);
                        keep_nearest(&mut near.left, left, &t);
                        keep_nearest(&mut near.right, right, &t);
                        near.end_file();
                        pooled += pool.len();
                    }
                    summary.stats.peak_pool = summary.stats.peak_pool.max(pooled);
                }
            }
            summary.stats.peak_active = summary.stats.peak_active.max(live);
            let rec = ARec {
                meta: *a,
                raw: batch.raw(i),
            };
            v.visit(&rec, &hits, &near)?;
            spare = recycle(hits);
            spare_near = near.recycle();
        }
        if !more {
            break;
        }
    }
    summary.b_records = streams.iter().map(|s| s.read).collect();
    summary.stats.peak_pool_bytes = pools.iter().map(LeftPool::peak_bytes).sum();
    Ok(summary)
}

/// A hit's `file` is a `u16`: at most 65,536 B inputs.
pub(crate) fn check_b_count(n: usize) -> Result<(), IoError> {
    let max_files = usize::from(u16::MAX) + 1;
    if n > max_files {
        return Err(UsageError(format!(
            "too many B inputs: {} (at most {})",
            thousands(n as u64),
            thousands(max_files as u64)
        ))
        .into());
    }
    Ok(())
}

/// Match A against one B input's active set: `ActiveSet::scan` with the
/// overlap test of `cfg` (its window applied relative to A's strand).
#[allow(clippy::too_many_arguments)]
fn scan<'h>(
    set: &'h mut ActiveSet,
    f: u16,
    a: &RecordMeta,
    lim: u64,
    pad: u64,
    cfg: &SweepConfig,
    hits: &mut Vec<Hit<'h>>,
    near: Option<&mut NearScan<'_, 'h>>,
) -> usize {
    let ip = cfg.inclusive_points;
    match cfg.window {
        None => set.scan(f, lim, a.key, 0, |b| overlaps(a.span, b, ip), hits, near),
        Some(w) => set.scan(
            f,
            lim,
            a.key,
            pad,
            |b| overlaps(a.span, w.apply(b, a.strand), ip),
            hits,
            near,
        ),
    }
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
    /// except those already ending before A's key (with the pad), which go
    /// to `pool` when there is one. A record on another chromosome is a
    /// later one: the order guard has checked it.
    fn load(
        &mut self,
        chroms: &mut ChromTable,
        a: &RecordMeta,
        lim: u64,
        pad: u64,
        set: &mut ActiveSet,
        mut pool: Option<&mut LeftPool>,
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
            } else if let Some(p) = pool.as_deref_mut() {
                p.offer(m.span, m.strand, m.loc, self.batch.raw(i), a.key);
            }
            self.pos += 1;
        }
        self.place = Place::Past;
        Ok(())
    }

    /// After `load`: read records past the load limit into `set` until no
    /// later one can be among the nearest right of A, counting those
    /// already there first (see `nearest`). Returns the largest key that
    /// still can be.
    fn look_ahead(
        &mut self,
        chroms: &mut ChromTable,
        a: &RecordMeta,
        lim: u64,
        mut ahead: LookAhead<'_>,
        set: &mut ActiveSet,
    ) -> Result<u64, IoError> {
        for (key, start, strand) in set.past(lim) {
            if key > ahead.bound() {
                return Ok(ahead.bound());
            }
            ahead.add(start, strand);
        }
        if self.place != Place::Same {
            return Ok(ahead.bound());
        }
        while let Some(i) = self.peek(chroms)? {
            let m = &self.batch.recs[i];
            if m.chrom != a.chrom {
                self.place = Place::Past;
                return Ok(ahead.bound());
            }
            if m.key > ahead.bound() {
                return Ok(ahead.bound());
            }
            // `load` took every record up to the limit, so this one starts
            // after A's end and cannot be evicted for this A.
            debug_assert!(m.key > lim);
            set.push(m, self.batch.raw(i));
            ahead.add(m.span.start, m.strand);
            self.pos += 1;
        }
        self.place = Place::Past;
        Ok(ahead.bound())
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
