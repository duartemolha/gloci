//! Unsorted mode (`--unsorted`, D39, D40, spec §8.2). Every B input is
//! loaded into memory and indexed per chromosome; then A is streamed in any
//! order, and each A queries the indexes. Neither side is order-checked.
//! The visitor sees exactly what the sorted sweep shows it for each A: the
//! same hits in the same order (D38) and the same [`Neighbours`].
//!
//! - **Index.** One [`Lapper<u64, u32>`] per B input and chromosome,
//!   holding positions in that chromosome's B table. Lapper's intervals
//!   are half-open and an empty one never matches, so B `[start, end)` is
//!   stored as `[start, end + 1)` (saturating): the closed range
//!   `[start, end]`, which a point never empties. A B starting at
//!   `u64::MAX` has no such range and is kept in a side list.
//! - **Candidates.** An A queries the closed range
//!   `[a.start - pad - 1, a.end + pad + 1]` (saturating), which takes in
//!   every B whose padded span can overlap it. Each candidate is re-tested
//!   with the sorted sweep's rule, and the hits are sorted into B input
//!   order (D38).
//! - **Sorted A.** While A arrives grouped by chromosome with query starts
//!   non-decreasing within each, lapper is queried through `seek` with one
//!   cursor per index, which resumes where the previous query began. From
//!   the first A out of that order on, plain `find`.
//! - **Neighbours** (`closest`). Each chromosome also keeps its B sorted by
//!   end and by start, in the strand buckets of the left pools. The left
//!   side walks back from the last B ending at or before A's start, the
//!   right side forward from the first B starting at or after A's end, each
//!   until it holds the `k` nearest and their ties. Distances never
//!   decrease along either walk (gloci-core `distance` is the gap plus 1).
//! - **Memory** is proportional to the B inputs (D39).

use crate::hit::{ARec, Hit, Neighbours, Visitor};
use crate::nearest::{Buckets, Target, keep_nearest};
use crate::sweep::{Inputs, NearestConfig, SweepConfig, SweepStats, SweepSummary, check_b_count};
use gloci_core::{ChromId, ChromTable, Interval, Strand, overlaps, thousands};
use gloci_io::{IoError, RecordBatch, RecordReader, UsageError, next_batch_size};
use rust_lapper::{Interval as LapInterval, Lapper};

type Index = Lapper<u64, u32>;

/// Run unsorted mode: load every B, then call `v` for every A record, in
/// A's input order, with its candidate hits (D36, D38) and, when
/// `cfg.nearest` is set, its neighbours, exactly as [`sweep_sorted`] would.
/// A's records and B's are interned in `chroms`; neither is order-checked.
///
/// [`sweep_sorted`]: crate::sweep_sorted
pub fn sweep_unsorted(
    inputs: Inputs,
    chroms: &mut ChromTable,
    cfg: &SweepConfig,
    v: &mut dyn Visitor,
) -> Result<SweepSummary, IoError> {
    let Inputs { a: mut reader, b } = inputs;
    check_b_count(b.len())?;
    let pad = cfg.window.map_or(0, |w| w.left.max(w.right));
    let nearest = cfg.nearest.filter(|n| n.k > 0);
    let mut stats = SweepStats::default();
    let mut b_records = Vec::with_capacity(b.len());
    let mut loaded = Vec::with_capacity(b.len());
    for r in b {
        let (input, n) = BInput::load(r, chroms, nearest.as_ref(), &mut stats)?;
        loaded.push(input);
        b_records.push(n);
    }
    // B is fixed from here on; hits borrow it.
    let loaded = &loaded[..];
    let mut summary = SweepSummary {
        a_records: 0,
        b_records,
        shared_chroms: false,
        skipped: Vec::new(),
        stats,
    };
    let mut order = AOrder::default();
    // Per B input: the seek cursor of the current chromosome.
    let mut cursors = vec![0usize; loaded.len()];
    let mut hits: Vec<Hit<'_>> = Vec::new();
    let mut near: Neighbours<'_> = Neighbours::default();
    let mut cand: Vec<u32> = Vec::new();
    let (mut batch, mut max) = (RecordBatch::default(), 1);
    loop {
        batch.clear();
        let more = reader.read_batch(chroms, &mut batch, max)?;
        max = next_batch_size(&batch, max);
        summary.a_records += batch.len() as u64;
        for (i, a) in batch.recs.iter().enumerate() {
            let lo = a.span.start.saturating_sub(pad).saturating_sub(1);
            let hi = a.span.end.saturating_add(pad).saturating_add(1);
            // The query start is what `seek` needs in order.
            let seek = match order.next(a.chrom, lo) {
                Step::Same => true,
                Step::Entered => {
                    cursors.fill(0);
                    true
                }
                Step::OutOfOrder | Step::Unsorted => false,
            };
            let target = nearest.map(|n| Target::new(&n, a.span, a.strand, cfg.inclusive_points));
            let ip = cfg.inclusive_points;
            let is_hit = |b: Interval| match cfg.window {
                None => overlaps(a.span, b, ip),
                Some(w) => overlaps(a.span, w.apply(b, a.strand), ip),
            };
            hits.clear();
            near = std::mem::take(&mut near).recycle();
            for (f, (input, cursor)) in loaded.iter().zip(&mut cursors).enumerate() {
                let f = f as u16;
                let Some(c) = input.chrom(a.chrom) else {
                    if target.is_some() {
                        near.end_file();
                    }
                    continue;
                };
                summary.shared_chroms = true;
                c.candidates(lo, hi, seek.then_some(&mut *cursor), &mut cand);
                summary.stats.seek_queries += u64::from(seek);
                cand.retain(|&j| is_hit(c.recs[j as usize].span));
                cand.sort_unstable();
                hits.extend(
                    cand.iter()
                        .map(|&j| c.recs[j as usize].hit(f, &input.arena)),
                );
                if let Some(t) = &target {
                    c.neighbours(&input.arena, f, t, &mut near);
                }
            }
            let rec = ARec {
                meta: *a,
                raw: batch.raw(i),
            };
            v.visit(&rec, &hits, &near)?;
        }
        if !more {
            break;
        }
    }
    Ok(summary)
}

/// Whether A has arrived grouped by chromosome with non-decreasing keys
/// (the query starts) within each chromosome, so far.
struct AOrder {
    sorted: bool,
    last: Option<(ChromId, u64)>,
    /// By chromosome id: has A been on it?
    seen: Vec<bool>,
}

impl Default for AOrder {
    fn default() -> AOrder {
        AOrder {
            sorted: true,
            last: None,
            seen: Vec::new(),
        }
    }
}

/// What the next A record does to the order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Step {
    /// Still in order, on the same chromosome.
    Same,
    /// Still in order, on a chromosome A has not been on before.
    Entered,
    /// The first record out of order.
    OutOfOrder,
    /// A was already out of order.
    Unsorted,
}

impl AOrder {
    fn next(&mut self, chrom: ChromId, key: u64) -> Step {
        if !self.sorted {
            return Step::Unsorted;
        }
        let step = match self.last {
            Some((c, k)) if c == chrom => {
                if key >= k {
                    Step::Same
                } else {
                    Step::OutOfOrder
                }
            }
            _ => {
                let i = chrom.0 as usize;
                if self.seen.len() <= i {
                    self.seen.resize(i + 1, false);
                }
                if std::mem::replace(&mut self.seen[i], true) {
                    Step::OutOfOrder
                } else {
                    Step::Entered
                }
            }
        };
        self.sorted = step != Step::OutOfOrder;
        self.last = Some((chrom, key));
        step
    }
}

/// A loaded B record; its raw bytes are `arena[off..off + len]` in its
/// input's arena.
#[derive(Debug, Clone, Copy)]
struct BRec {
    span: Interval,
    loc: u64,
    off: usize,
    len: u32,
    strand: Strand,
}

impl BRec {
    fn hit<'b>(&self, file: u16, arena: &'b [u8]) -> Hit<'b> {
        Hit {
            file,
            span: self.span,
            strand: self.strand,
            raw: &arena[self.off..self.off + self.len as usize],
            loc: self.loc,
        }
    }
}

/// One loaded B input: every record's raw bytes, and its records by
/// chromosome.
struct BInput {
    arena: Vec<u8>,
    /// By chromosome id; `None` where the input has no records.
    chroms: Vec<Option<ChromB>>,
}

impl BInput {
    /// Read every record of `reader`. Returns the input and its record
    /// count.
    fn load(
        mut reader: Box<dyn RecordReader>,
        chroms: &mut ChromTable,
        nearest: Option<&NearestConfig>,
        stats: &mut SweepStats,
    ) -> Result<(BInput, u64), IoError> {
        let mut arena = Vec::new();
        let mut by_chrom: Vec<Vec<BRec>> = Vec::new();
        let (mut batch, mut max, mut read) = (RecordBatch::default(), 1, 0);
        loop {
            batch.clear();
            let more = reader.read_batch(chroms, &mut batch, max)?;
            max = next_batch_size(&batch, max);
            read += batch.len() as u64;
            for (i, m) in batch.recs.iter().enumerate() {
                let c = m.chrom.0 as usize;
                if by_chrom.len() <= c {
                    by_chrom.resize_with(c + 1, Vec::new);
                }
                let off = arena.len();
                arena.extend_from_slice(batch.raw(i));
                by_chrom[c].push(BRec {
                    span: m.span,
                    loc: m.loc,
                    off,
                    len: m.len,
                    strand: m.strand,
                });
            }
            if !more {
                break;
            }
        }
        // Memory is the cost of this mode (D39): drop the growth slack.
        arena.shrink_to_fit();
        let mut per_chrom = Vec::with_capacity(by_chrom.len());
        for (c, mut recs) in by_chrom.into_iter().enumerate() {
            if recs.is_empty() {
                per_chrom.push(None);
                continue;
            }
            if recs.len() >= u32::MAX as usize {
                return Err(UsageError(format!(
                    "{}: {} records on {} (--unsorted indexes at most {} per chromosome); \
                     sort the inputs and run without --unsorted",
                    reader.name(),
                    thousands(recs.len() as u64),
                    String::from_utf8_lossy(chroms.name(ChromId(c as u32))),
                    thousands(u64::from(u32::MAX) - 1)
                ))
                .into());
            }
            recs.shrink_to_fit();
            per_chrom.push(Some(ChromB::new(recs, nearest, stats)));
        }
        Ok((
            BInput {
                arena,
                chroms: per_chrom,
            },
            read,
        ))
    }

    fn chrom(&self, c: ChromId) -> Option<&ChromB> {
        self.chroms.get(c.0 as usize)?.as_ref()
    }
}

/// One B input's records on one chromosome, in input order, and their
/// index, whose ids are positions in `recs`.
struct ChromB {
    recs: Vec<BRec>,
    index: Index,
    /// Positions of the B starting at `u64::MAX`, which `index` cannot hold.
    at_max: Vec<u32>,
    /// With a nearest request: per bucket ([`Buckets`]), the positions of
    /// its B sorted by end, and sorted by start.
    by_end: [Vec<u32>; 2],
    by_start: [Vec<u32>; 2],
}

impl ChromB {
    fn new(recs: Vec<BRec>, nearest: Option<&NearestConfig>, stats: &mut SweepStats) -> ChromB {
        stats.lapper_indexes += 1;
        let (mut items, mut at_max) = (Vec::with_capacity(recs.len()), Vec::new());
        for (i, r) in recs.iter().enumerate() {
            if r.span.start == u64::MAX {
                at_max.push(i as u32);
            } else {
                items.push(LapInterval {
                    start: r.span.start,
                    stop: r.span.end.saturating_add(1),
                    val: i as u32,
                });
            }
        }
        let index = Lapper::new(items);
        let (mut by_end, mut by_start): ([Vec<u32>; 2], [Vec<u32>; 2]) = Default::default();
        if let Some(n) = nearest {
            let buckets = Buckets::new(n);
            for (i, r) in recs.iter().enumerate() {
                if let Some(g) = buckets.of(r.strand) {
                    by_end[g].push(i as u32);
                }
            }
            by_start.clone_from(&by_end);
            for g in 0..2 {
                by_end[g].sort_unstable_by_key(|&i| (recs[i as usize].span.end, i));
                by_start[g].sort_unstable_by_key(|&i| (recs[i as usize].span.start, i));
            }
        }
        ChromB {
            recs,
            index,
            at_max,
            by_end,
            by_start,
        }
    }

    /// Replace `out` with the positions of the B whose closed range
    /// `[start, end]` meets `[lo, hi]`. With a `cursor`, which must be
    /// this index's, queries through `seek`; the caller guarantees the
    /// query starts do not decrease between calls on one cursor.
    fn candidates(&self, lo: u64, hi: u64, cursor: Option<&mut usize>, out: &mut Vec<u32>) {
        out.clear();
        // Lapper takes `[lo, hi + 1)`; `hi + 1` saturates only for B
        // starting at `u64::MAX`, which are the side list.
        let stop = hi.saturating_add(1);
        match cursor {
            Some(c) => out.extend(self.index.seek(lo, stop, c).map(|iv| iv.val)),
            None => out.extend(self.index.find(lo, stop).map(|iv| iv.val)),
        }
        if hi == u64::MAX {
            out.extend_from_slice(&self.at_max);
        }
    }

    /// Push A's neighbours among this chromosome's B (as hits of `file`)
    /// and close the file's runs: the `k` nearest plus ties on each side,
    /// in each bucket A wants, sorted as [`Neighbours`] states.
    fn neighbours<'b>(&self, arena: &'b [u8], file: u16, t: &Target, near: &mut Neighbours<'b>) {
        let (left, right) = (near.left.len(), near.right.len());
        for g in t.want.buckets() {
            // Left: B ending at or before A's start, latest end first.
            let order = &self.by_end[g];
            let p = order.partition_point(|&i| self.recs[i as usize].span.end <= t.span.start);
            self.take_nearest(order[..p].iter().rev(), arena, file, t, &mut near.left);
            // Right: B starting at or after A's end, earliest start first.
            let order = &self.by_start[g];
            let p = order.partition_point(|&i| self.recs[i as usize].span.start < t.span.end);
            self.take_nearest(order[p..].iter(), arena, file, t, &mut near.right);
        }
        keep_nearest(&mut near.left, left, t);
        keep_nearest(&mut near.right, right, t);
        near.end_file();
    }

    /// Push the B of `order`, whose distances from A never decrease, that
    /// do not overlap A and are within `max_distance`, until the `k`
    /// nearest and every B tied with the k-th are in.
    fn take_nearest<'b, 'o>(
        &self,
        order: impl Iterator<Item = &'o u32>,
        arena: &'b [u8],
        file: u16,
        t: &Target,
        out: &mut Vec<Hit<'b>>,
    ) {
        let (mut kept, mut kth) = (0, None);
        for &i in order {
            let b = &self.recs[i as usize];
            // Only a B touching A's edge can overlap it here.
            if t.overlaps(b.span) {
                continue;
            }
            let d = t.distance(b.span);
            if !t.within(b.span) || kth.is_some_and(|k| d > k) {
                break;
            }
            out.push(b.hit(file, arena));
            kept += 1;
            if kept == t.k {
                kth = Some(d);
            }
        }
    }
}
