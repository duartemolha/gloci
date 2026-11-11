//! Twin of `crate::merge`: [`MergedReader`] re-implements its clustering (no aggregation, pulled
//! instead of pushed). Change one, change both: `tests/jaccard_props.rs`
//! (`merged_reader_matches_merge`) fails when they diverge.
//!
//! Stream adapters for the interval statistics: [`MergedReader`] turns a
//! sorted input into its merged clusters (as `gloci merge` writes them),
//! and [`Tally`] counts what flows through.

use crate::merge::{joins, slot_of};
use gloci_core::{ChromId, ChromTable, Interval, Strand};
use gloci_io::{
    Compression, Diagnostics, Format, Header, IoError, OrderGuard, RecordBatch, RecordMeta,
    RecordReader, next_batch_size,
};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

/// A footprint waiting for its turn in start order.
struct Pending {
    start: u64,
    end: u64,
    slot: usize,
    seq: u64,
}

impl Pending {
    fn key(&self) -> (u64, u64) {
        (self.start, self.seq)
    }
}

impl PartialEq for Pending {
    fn eq(&self, o: &Self) -> bool {
        self.key() == o.key()
    }
}
impl Eq for Pending {}
impl PartialOrd for Pending {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Pending {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.key().cmp(&o.key())
    }
}

struct Open {
    start: u64,
    end: u64,
}

/// A finished cluster.
#[derive(Clone, Copy)]
struct Cluster {
    chrom: ChromId,
    start: u64,
    end: u64,
    strand: Strand,
}

/// A finished cluster, ordered by (start, end, strand slot, creation).
struct Closed {
    order: (u64, u64, usize, u64),
    cluster: Cluster,
}

impl PartialEq for Closed {
    fn eq(&self, o: &Self) -> bool {
        self.order == o.order
    }
}
impl Eq for Closed {}
impl PartialOrd for Closed {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Closed {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.order.cmp(&o.order)
    }
}

fn strand_of_slot(per_strand: bool, slot: usize) -> Strand {
    match (per_strand, slot) {
        (true, 0) => Strand::Plus,
        (true, 1) => Strand::Minus,
        _ => Strand::Unknown,
    }
}

/// A sorted input read as the clusters `gloci merge` would write (records
/// that overlap or are bookended join; with `by_strand`, `+`, `-` and unknown
/// strands cluster apart). Each record is `chrom\tstart\tend` with the
/// cluster's strand, in start order, so it is itself a sorted input. State is
/// the open clusters and the closed ones not yet in order, never the input.
///
/// The input is order-checked here (a violation is an order error with the
/// `--sort` hint), because the merged stream hides it from the sweep.
pub struct MergedReader {
    inner: Box<dyn RecordReader>,
    by_strand: bool,
    guard: OrderGuard,
    batch: RecordBatch,
    pos: usize,
    more: bool,
    max: usize,
    heap: BinaryHeap<Reverse<Pending>>,
    open: Vec<Option<Open>>,
    closed: BinaryHeap<Reverse<Closed>>,
    ready: VecDeque<Cluster>,
    chrom: Option<ChromId>,
    seq: u64,
    created: u64,
    emitted: u64,
    peak_held: usize,
    done: bool,
    header: Header,
    diagnostics: Diagnostics,
}

impl MergedReader {
    pub fn new(inner: Box<dyn RecordReader>, by_strand: bool) -> Self {
        let slots = if by_strand { 3 } else { 1 };
        MergedReader {
            guard: OrderGuard::new(inner.name(), inner.format()),
            inner,
            by_strand,
            batch: RecordBatch::default(),
            pos: 0,
            more: true,
            max: 1,
            heap: BinaryHeap::new(),
            open: (0..slots).map(|_| None).collect(),
            closed: BinaryHeap::new(),
            ready: VecDeque::new(),
            chrom: None,
            seq: 0,
            created: 0,
            emitted: 0,
            peak_held: 0,
            done: false,
            header: Header::None,
            diagnostics: Diagnostics::default(),
        }
    }

    /// The most footprints and clusters held at once (for the
    /// memory-flatness test).
    #[doc(hidden)]
    pub fn peak_held(&self) -> usize {
        self.peak_held
    }

    /// Process the next input record, or finish at the end of the input.
    fn step(&mut self, chroms: &mut ChromTable) -> Result<(), IoError> {
        if self.pos >= self.batch.len() {
            self.batch.clear();
            self.pos = 0;
            while self.more && self.batch.is_empty() {
                self.more = self.inner.read_batch(chroms, &mut self.batch, self.max)?;
                self.guard.check(chroms, &self.batch)?;
                self.max = next_batch_size(&self.batch, self.max);
            }
            if self.batch.is_empty() {
                self.finish_chrom();
                self.done = true;
                return Ok(());
            }
        }
        let m = self.batch.recs[self.pos];
        self.pos += 1;
        self.record(&m);
        let held = self.heap.len() + self.closed.len() + self.ready.len();
        self.peak_held = self.peak_held.max(held);
        Ok(())
    }

    fn record(&mut self, m: &RecordMeta) {
        if self.chrom != Some(m.chrom) {
            self.finish_chrom();
            self.chrom = Some(m.chrom);
        }
        // A span never starts before its key, so everything starting at or
        // before this key is final in start order.
        while self.heap.peek().is_some_and(|p| p.0.start <= m.key) {
            let Reverse(p) = self.heap.pop().expect("peeked");
            self.add(p);
        }
        self.heap.push(Reverse(Pending {
            start: m.span.start,
            end: m.span.end,
            slot: slot_of(self.by_strand, m.strand),
            seq: self.seq,
        }));
        self.seq += 1;
        // An open cluster that ends before this key can never grow: close it
        // now, or one stale cluster of another strand would hold every
        // closed cluster back until the chromosome ends.
        for slot in 0..self.open.len() {
            if self.open[slot].as_ref().is_some_and(|c| c.end < m.key) {
                self.close(slot);
            }
        }
        let horizon = self
            .open
            .iter()
            .flatten()
            .map(|c| c.start)
            .min()
            .map_or(m.key, |s| s.min(m.key));
        self.flush(Some(horizon));
    }

    fn add(&mut self, p: Pending) {
        let slot = p.slot;
        let joins = self.open[slot]
            .as_ref()
            .is_some_and(|c| joins(p.start, c.end, 0));
        if !joins {
            self.close(slot);
            self.open[slot] = Some(Open {
                start: p.start,
                end: p.end,
            });
        }
        let c = self.open[slot].as_mut().expect("opened");
        c.end = c.end.max(p.end);
    }

    fn close(&mut self, slot: usize) {
        let (Some(c), Some(chrom)) = (self.open[slot].take(), self.chrom) else {
            return;
        };
        self.created += 1;
        self.closed.push(Reverse(Closed {
            order: (c.start, c.end, slot, self.created),
            cluster: Cluster {
                chrom,
                start: c.start,
                end: c.end,
                strand: strand_of_slot(self.by_strand, slot),
            },
        }));
    }

    /// Move closed clusters that no later cluster can precede (those that
    /// start before `horizon`; all of them when `None`) to the ready queue.
    fn flush(&mut self, horizon: Option<u64>) {
        while let Some(Reverse(top)) = self.closed.peek() {
            if horizon.is_some_and(|h| top.order.0 >= h) {
                break;
            }
            let Reverse(c) = self.closed.pop().expect("peeked");
            self.ready.push_back(c.cluster);
        }
    }

    fn finish_chrom(&mut self) {
        while let Some(Reverse(p)) = self.heap.pop() {
            self.add(p);
        }
        for slot in 0..self.open.len() {
            self.close(slot);
        }
        self.flush(None);
    }
}

impl RecordReader for MergedReader {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        Format::Bed
    }
    fn compression(&self) -> Compression {
        Compression::None
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        while self.ready.is_empty() && !self.done {
            self.step(chroms)?;
        }
        for _ in 0..max.max(1) {
            let Some(c) = self.ready.pop_front() else {
                break;
            };
            let off = batch.buf.len();
            batch.buf.extend_from_slice(chroms.name(c.chrom));
            batch.buf.push(b'\t');
            batch.buf.extend_from_slice(c.start.to_string().as_bytes());
            batch.buf.push(b'\t');
            batch.buf.extend_from_slice(c.end.to_string().as_bytes());
            self.emitted += 1;
            batch.recs.push(RecordMeta {
                chrom: c.chrom,
                key: c.start,
                span: Interval {
                    start: c.start,
                    end: c.end,
                },
                strand: c.strand,
                loc: self.emitted,
                off,
                len: u32::try_from(batch.buf.len() - off).unwrap_or(u32::MAX),
            });
        }
        Ok(!(self.done && self.ready.is_empty()))
    }
    /// The skipped records of the input are reported by whoever watches the
    /// input itself.
    fn diagnostics(&self) -> &Diagnostics {
        &self.diagnostics
    }
    fn consume_to_end(&self) -> bool {
        true
    }
}

/// What a [`Tally`] saw.
#[derive(Debug, Default)]
pub struct TallyCounts {
    pub records: AtomicU64,
    pub bases: AtomicU64,
}

/// Passes records through and counts them and their bases. It asks the
/// sweep to read the input to its end, so the totals are complete.
pub struct Tally {
    inner: Box<dyn RecordReader>,
    counts: Arc<TallyCounts>,
}

impl Tally {
    pub fn new(inner: Box<dyn RecordReader>) -> (Tally, Arc<TallyCounts>) {
        let counts = Arc::new(TallyCounts::default());
        (
            Tally {
                inner,
                counts: counts.clone(),
            },
            counts,
        )
    }
}

impl RecordReader for Tally {
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
        let new = &batch.recs[before..];
        let bases: u64 = new.iter().map(|m| m.span.end - m.span.start).sum();
        self.counts
            .records
            .fetch_add(new.len() as u64, Ordering::Relaxed);
        self.counts.bases.fetch_add(bases, Ordering::Relaxed);
        Ok(more)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn consume_to_end(&self) -> bool {
        true
    }
    fn whole_records(&self) -> bool {
        self.inner.whole_records()
    }
}

/// Drops the records of unknown strand (`-s`/`-S` can never match them) and
/// notes that it did.
pub struct DropUnknown {
    inner: Box<dyn RecordReader>,
    scratch: RecordBatch,
    dropped: Arc<std::sync::atomic::AtomicBool>,
}

impl DropUnknown {
    pub fn new(inner: Box<dyn RecordReader>) -> (DropUnknown, Arc<std::sync::atomic::AtomicBool>) {
        let dropped = Arc::new(std::sync::atomic::AtomicBool::new(false));
        (
            DropUnknown {
                inner,
                scratch: RecordBatch::default(),
                dropped: dropped.clone(),
            },
            dropped,
        )
    }
}

impl RecordReader for DropUnknown {
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
        self.scratch.clear();
        let more = self.inner.read_batch(chroms, &mut self.scratch, max)?;
        for (i, m) in self.scratch.recs.iter().enumerate() {
            if m.strand == Strand::Unknown {
                self.dropped.store(true, Ordering::Relaxed);
                continue;
            }
            let off = batch.buf.len();
            batch.buf.extend_from_slice(self.scratch.raw(i));
            batch.recs.push(RecordMeta { off, ..*m });
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
