//! Unsorted mode (`--unsorted`, D39, D40, spec §8.2). Every B input is
//! loaded into memory and indexed per chromosome; then A is streamed in any
//! order, and each A queries the indexes. Neither side is order-checked.
//! The visitor sees exactly what the sorted sweep shows it for each A: the
//! same hits in the same order (D38) and the same [`Neighbours`].
//!
//! - **Index.** One per B input and chromosome: a coitree when every B end
//!   there, padded by the window, is below `i32::MAX - 2`; an [`ITree64`]
//!   otherwise. Both hold B `[start, end)` as the closed range
//!   `[start, end]`, so a point is never an empty interval.
//! - **Candidates.** An A queries `[a.start - pad - 1, a.end + pad + 1]`
//!   (saturating), which takes in every B whose padded span can overlap
//!   it. Each candidate is re-tested with the sorted sweep's rule, and the
//!   hits are sorted into B input order (D38).
//! - **Sorted A.** While A arrives grouped by chromosome with keys
//!   non-decreasing within each, coitrees are queried through a sorted
//!   querent (one per tree), which reuses the previous query's results. A
//!   querent's state is per tree, so only the order within a chromosome
//!   matters to it. From the first A out of that order on, plain queries.
//! - **Neighbours** (`closest`). Each chromosome also keeps its B sorted by
//!   end and by start, in the strand buckets of the left pools. The left
//!   side walks back from the last B ending at or before A's start, the
//!   right side forward from the first B starting at or after A's end, each
//!   until it holds the `k` nearest and their ties. Distances never
//!   decrease along either walk (gloci-core `distance` is the gap plus 1).
//! - **Memory** is proportional to the B inputs (D39).

use crate::hit::{ARec, Hit, Neighbours, Visitor};
use crate::itree64::ITree64;
use crate::nearest::{Buckets, Target, keep_nearest};
use crate::sweep::{Inputs, NearestConfig, SweepConfig, SweepStats, SweepSummary, check_b_count};
use coitrees::{COITree, COITreeSortedQuerent, GenericInterval, IntervalTree, SortedQuerent};
use gloci_core::{ChromId, ChromTable, Interval, Strand, overlaps, thousands};
use gloci_io::{IoError, RecordBatch, RecordReader, UsageError, next_batch_size};

/// Coitrees index a chromosome only while every padded B end there is
/// below this, and queries are clamped to it. The SIMD (NEON, AVX2)
/// coitrees store `first - 1` and `last + 1`, and their sorted querents
/// return a chunk's empty padding lanes (id 0) for a query ending at
/// `i32::MAX`.
const COI_LIMIT: u64 = (i32::MAX - 2) as u64;

type Coi = COITree<u32, u32>;
type Querent<'t> = COITreeSortedQuerent<'t, u32, u32>;

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
        let (input, n) = BInput::load(r, chroms, pad, nearest.as_ref(), &mut stats)?;
        loaded.push(input);
        b_records.push(n);
    }
    // B is fixed from here on; hits and querents borrow it.
    let loaded = &loaded[..];
    let mut summary = SweepSummary {
        a_records: 0,
        b_records,
        shared_chroms: false,
        skipped: Vec::new(),
        stats,
    };
    let mut order = AOrder::default();
    let mut querents: Vec<Option<Querent<'_>>> = loaded.iter().map(|_| None).collect();
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
            match order.next(a.chrom, a.key) {
                Step::Same => {}
                Step::Entered => {
                    for (q, input) in querents.iter_mut().zip(loaded) {
                        *q = input.chrom(a.chrom).and_then(ChromB::querent);
                    }
                }
                Step::OutOfOrder => querents.iter_mut().for_each(|q| *q = None),
                Step::Unsorted => {}
            }
            let lo = a.span.start.saturating_sub(pad).saturating_sub(1);
            let hi = a.span.end.saturating_add(pad).saturating_add(1);
            let target = nearest.map(|n| Target::new(&n, a.span, a.strand, cfg.inclusive_points));
            let ip = cfg.inclusive_points;
            let is_hit = |b: Interval| match cfg.window {
                None => overlaps(a.span, b, ip),
                Some(w) => overlaps(a.span, w.apply(b, a.strand), ip),
            };
            hits.clear();
            near = std::mem::take(&mut near).recycle();
            for (f, (input, q)) in loaded.iter().zip(&mut querents).enumerate() {
                let f = f as u16;
                let Some(c) = input.chrom(a.chrom) else {
                    if target.is_some() {
                        near.end_file();
                    }
                    continue;
                };
                summary.shared_chroms = true;
                let by_querent = c.candidates(lo, hi, q.as_mut(), &mut cand);
                summary.stats.querent_queries += u64::from(by_querent);
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
/// within each chromosome, so far.
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
        pad: u64,
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
            per_chrom.push(Some(ChromB::new(recs, pad, nearest, stats)));
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
    /// With a nearest request: per bucket ([`Buckets`]), the positions of
    /// its B sorted by end, and sorted by start.
    by_end: [Vec<u32>; 2],
    by_start: [Vec<u32>; 2],
}

enum Index {
    Coi(Coi),
    Wide(ITree64),
}

impl ChromB {
    fn new(
        recs: Vec<BRec>,
        pad: u64,
        nearest: Option<&NearestConfig>,
        stats: &mut SweepStats,
    ) -> ChromB {
        let index = if recs
            .iter()
            .all(|r| r.span.end.saturating_add(pad) < COI_LIMIT)
        {
            stats.coitrees_indexes += 1;
            let nodes: Vec<coitrees::Interval<u32>> = recs
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    coitrees::Interval::new(r.span.start as i32, r.span.end as i32, i as u32)
                })
                .collect();
            Index::Coi(Coi::new(&nodes))
        } else {
            stats.itree64_indexes += 1;
            let items = recs
                .iter()
                .enumerate()
                .map(|(i, r)| (r.span, i as u32))
                .collect();
            Index::Wide(ITree64::new(items))
        };
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
            by_end,
            by_start,
        }
    }

    /// A sorted querent over this chromosome's coitree, if it has one.
    fn querent(&self) -> Option<Querent<'_>> {
        match &self.index {
            Index::Coi(t) => Some(Querent::new(t)),
            Index::Wide(_) => None,
        }
    }

    /// Replace `out` with the positions of the B whose closed range
    /// `[start, end]` meets `[lo, hi]`. `querent`, when given, must be
    /// this chromosome's. Returns whether it answered.
    fn candidates(
        &self,
        lo: u64,
        hi: u64,
        querent: Option<&mut Querent<'_>>,
        out: &mut Vec<u32>,
    ) -> bool {
        out.clear();
        match &self.index {
            Index::Wide(t) => t.query(lo, hi, out),
            Index::Coi(t) => {
                // Every B coordinate here is below the limit, so clamping
                // the query to it changes no answer.
                let (lo, hi) = (lo.min(COI_LIMIT) as i32, hi.min(COI_LIMIT) as i32);
                match querent {
                    Some(q) => {
                        querent_query(q, lo, hi, out);
                        return true;
                    }
                    None => tree_query(t, lo, hi, out),
                }
            }
        }
        false
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

/// Push the ids of the intervals of `tree` meeting `[lo, hi]`. Generic so
/// it works with every coitrees variant (`COITree` is the AVX2, NEON or
/// basic tree, by target).
fn tree_query<'t, T: IntervalTree<'t, Metadata = u32>>(
    tree: &'t T,
    lo: i32,
    hi: i32,
    out: &mut Vec<u32>,
) {
    tree.query(lo, hi, |n| out.push(*GenericInterval::<u32>::metadata(n)));
}

/// [`tree_query`] through a sorted querent.
fn querent_query<'t, Q: SortedQuerent<'t, Metadata = u32>>(
    q: &mut Q,
    lo: i32,
    hi: i32,
    out: &mut Vec<u32>,
) {
    q.query(lo, hi, |n| out.push(*GenericInterval::<u32>::metadata(n)));
}

#[cfg(test)]
mod tests {
    use super::*;
    use coitrees::{BasicCOITree, BasicSortedQuerent};

    /// The basic coitree (x86-64 without AVX2) answers like this target's
    /// `COITree` through the same generic helpers, so both build and agree.
    #[test]
    fn every_coitree_variant_gives_the_same_candidates() {
        let nodes: Vec<coitrees::Interval<u32>> = (0..500u32)
            .map(|i| {
                let first = (i * 37 % 1000) as i32;
                coitrees::Interval::new(first, first + (i % 50) as i32, i)
            })
            .collect();
        let basic = BasicCOITree::<u32, u32>::new(&nodes);
        let here = Coi::new(&nodes);
        let mut bq = BasicSortedQuerent::new(&basic);
        let mut hq = Querent::new(&here);
        for (lo, hi) in [
            (0, 0),
            (5, 9),
            (3, 400),
            (390, 410),
            (100, 120),
            (999, 1100),
        ] {
            let mut want: Vec<u32> = nodes
                .iter()
                .filter(|n| n.first <= hi && n.last >= lo)
                .map(|n| n.metadata)
                .collect();
            want.sort_unstable();
            let run = |f: &mut dyn FnMut(&mut Vec<u32>)| {
                let mut v = Vec::new();
                f(&mut v);
                v.sort_unstable();
                v
            };
            assert_eq!(run(&mut |v| tree_query(&basic, lo, hi, v)), want);
            assert_eq!(run(&mut |v| tree_query(&here, lo, hi, v)), want);
            assert_eq!(run(&mut |v| querent_query(&mut bq, lo, hi, v)), want);
            assert_eq!(run(&mut |v| querent_query(&mut hq, lo, hi, v)), want);
        }
    }
}
