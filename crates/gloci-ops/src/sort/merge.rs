//! k-way merges: of spilled runs (with cascading passes, D24) and of
//! already sorted inputs (`sort -m`).

use super::Tie;
use super::spill::{RecHead, Run, RunReader, SpillDir};
use gloci_core::ChromTable;
use gloci_io::{IoError, OrderGuard, RecordBatch, RecordReader};
use std::cmp::Ordering;
use std::collections::binary_heap::PeekMut;
use std::collections::{BinaryHeap, VecDeque};
use std::io;
use std::path::Path;

/// The most runs merged at once (D24).
pub(crate) const MAX_FAN_IN: usize = 128;

/// Estimated buffer bytes per open input in `merge_sorted`.
pub const MERGE_BUFFER_PER_INPUT: usize = 384 << 10;

/// How many inputs to merge at once under memory budget `max_mem`, capped at 128.
pub fn merge_fan_in(max_mem: u64) -> usize {
    (max_mem as usize / MERGE_BUFFER_PER_INPUT).clamp(2, MAX_FAN_IN)
}

/// How many runs to merge at once: up to [`MAX_FAN_IN`], fewer when the
/// longest record times the fan-in would exceed the budget (every open run
/// holds its current record), but at least 2. Depends only on the budget and
/// the input, and the output does not depend on it.
pub(crate) fn fan_in(budget: usize, max_rec: usize) -> usize {
    (budget / max_rec.max(1)).clamp(2, MAX_FAN_IN)
}

/// The current record of one source in a merge. `BinaryHeap` is a max-heap,
/// so the order is reversed: the smallest record is on top.
struct Head<const FULL: bool> {
    rank: u32,
    rec: RecHead,
    src: usize,
    bytes: Vec<u8>,
}

impl<const FULL: bool> Head<FULL> {
    fn cmp_key(&self, other: &Self) -> Ordering {
        let ord = (self.rank, self.rec.key, self.rec.end).cmp(&(
            other.rank,
            other.rec.key,
            other.rec.end,
        ));
        let ord = if FULL {
            ord.then_with(|| self.bytes.cmp(&other.bytes))
        } else {
            ord
        };
        ord.then(self.rec.seq.cmp(&other.rec.seq))
    }
}

impl<const FULL: bool> Ord for Head<FULL> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_key(other).reverse()
    }
}

impl<const FULL: bool> PartialOrd for Head<FULL> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const FULL: bool> PartialEq for Head<FULL> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<const FULL: bool> Eq for Head<FULL> {}

fn rank_of(ranks: &[u32], chrom: u32) -> Result<u32, IoError> {
    ranks
        .get(chrom as usize)
        .copied()
        .ok_or_else(|| IoError::Read {
            path: "sort spill file".into(),
            source: io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unknown chromosome id {chrom}"),
            ),
        })
}

/// Receives merged records: another run (intermediate pass) or the output.
pub(crate) type Sink<'a> = dyn FnMut(&RecHead, &[u8]) -> Result<(), IoError> + 'a;

/// Merge `runs` in `(rank, key, end, tiebreak)` order into `sink`.
pub(crate) fn merge_runs<const FULL: bool>(
    runs: &[Run],
    ranks: &[u32],
    sink: &mut Sink,
) -> Result<(), IoError> {
    let mut readers = runs
        .iter()
        .map(RunReader::open)
        .collect::<Result<Vec<_>, _>>()?;
    let mut heap = BinaryHeap::with_capacity(readers.len());
    for (src, r) in readers.iter_mut().enumerate() {
        let mut rec = RecHead::default();
        let mut bytes = Vec::new();
        if r.next(&mut rec, &mut bytes)? {
            let rank = rank_of(ranks, rec.chrom)?;
            heap.push(Head::<FULL> {
                rank,
                rec,
                src,
                bytes,
            });
        }
    }
    loop {
        let Some(mut top) = heap.peek_mut() else {
            break;
        };
        sink(&top.rec, &top.bytes)?;
        let h = &mut *top;
        if readers[h.src].next(&mut h.rec, &mut h.bytes)? {
            // Replaced in place; dropping the PeekMut sifts it down.
            h.rank = rank_of(ranks, h.rec.chrom)?;
        } else {
            PeekMut::pop(top);
        }
    }
    Ok(())
}

/// Bytes of records to read from one input at a time in `merge_sorted`.
const BATCH_BYTES: usize = 256 << 10;
const MAX_BATCH: usize = 4096;

/// One input of `merge_sorted`: its current batch and position.
struct Input<'a> {
    reader: &'a mut dyn RecordReader,
    guard: OrderGuard,
    batch: RecordBatch,
    pos: usize,
    more: bool,
    /// Records to ask for per batch; adapted to the record length.
    max: usize,
}

impl Input<'_> {
    /// Load the next non-empty batch, checking its order; false at the end.
    fn refill(&mut self, chroms: &mut ChromTable) -> Result<bool, IoError> {
        self.pos = 0;
        loop {
            self.batch.clear();
            if !self.more {
                return Ok(false);
            }
            self.more = self.reader.read_batch(chroms, &mut self.batch, self.max)?;
            self.guard.check(chroms, &self.batch)?;
            if !self.batch.is_empty() {
                let per_rec = self.batch.buf.len() / self.batch.len() + 1;
                self.max = (BATCH_BYTES / per_rec).clamp(1, MAX_BATCH);
                return Ok(true);
            }
        }
    }

    fn head<const FULL: bool>(&self, src: usize, base_seq: u64, ranks: &[u32]) -> Head<FULL> {
        let m = &self.batch.recs[self.pos];
        Head {
            rank: ranks[m.chrom.0 as usize],
            rec: RecHead {
                chrom: m.chrom.0,
                key: m.key,
                end: m.span.end,
                // Ties between inputs go to the earlier input.
                seq: base_seq + src as u64,
            },
            src,
            bytes: if FULL {
                self.batch.raw(self.pos).to_vec()
            } else {
                Vec::new()
            },
        }
    }

    /// Refresh `h` from the current record (`h.src` is unchanged).
    fn update<const FULL: bool>(&self, h: &mut Head<FULL>, ranks: &[u32]) {
        let m = &self.batch.recs[self.pos];
        h.rank = ranks[m.chrom.0 as usize];
        h.rec.chrom = m.chrom.0;
        h.rec.key = m.key;
        h.rec.end = m.span.end;
        if FULL {
            h.bytes.clear();
            h.bytes.extend_from_slice(self.batch.raw(self.pos));
        }
    }
}

/// Cascade-merge `runs` down to `out`, merging up to `fan` runs at once.
pub(crate) fn cascade_merge_runs<const FULL: bool>(
    runs: Vec<Run>,
    spill: &mut SpillDir,
    ranks: &[u32],
    fan: usize,
    sink: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    let mut count = 0u64;
    let mut runs: VecDeque<Run> = runs.into();
    while runs.len() > fan {
        let k = fan.min(runs.len() - fan + 1);
        let group: Vec<Run> = runs.drain(..k).collect();
        let mut w = spill.create_run()?;
        merge_runs::<FULL>(&group, ranks, &mut |h, raw| w.write(h, raw))?;
        runs.push_back(w.finish()?);
        group.into_iter().for_each(Run::remove);
    }
    merge_runs::<FULL>(runs.make_contiguous(), ranks, &mut |_, raw| {
        sink(raw).map_err(super::out_err)?;
        count += 1;
        Ok(())
    })?;
    Ok(count)
}

/// Merge inputs in batches of up to `fan_in` into LZ4 spill runs, closing
/// each batch's readers before opening the next batch, then cascade-merge the runs.
pub fn merge_sorted_batches<F>(
    num_inputs: usize,
    chroms: &mut ChromTable,
    tie: Tie,
    max_mem: u64,
    tmp_dir: &Path,
    mut open_batch: F,
    out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> Result<u64, IoError>
where
    F: FnMut(usize, usize) -> Result<Vec<Box<dyn RecordReader>>, IoError>,
{
    let fan = merge_fan_in(max_mem);
    if num_inputs <= fan {
        let mut readers = open_batch(0, num_inputs)?;
        let res = merge_sorted(&mut readers, chroms, tie, tmp_dir, out);
        for r in &readers {
            if let Some(s) = r.diagnostics().summary(r.name()) {
                log::warn!("{s}");
            }
        }
        return res;
    }
    let mut spill = SpillDir::create(tmp_dir)?;
    let mut runs = Vec::new();
    let mut start = 0;
    while start < num_inputs {
        let count = fan.min(num_inputs - start);
        let mut readers = open_batch(start, count)?;
        let mut w = spill.create_run()?;
        let mut sink = |h: &RecHead, raw: &[u8]| w.write(h, raw).map_err(io::Error::other);
        match tie {
            Tie::InputOrder => {
                merge_inputs::<false>(&mut readers, chroms, start as u64, &mut sink)?;
            }
            Tie::FullLine => {
                merge_inputs::<true>(&mut readers, chroms, start as u64, &mut sink)?;
            }
        }
        runs.push(w.finish()?);
        for r in &readers {
            if let Some(s) = r.diagnostics().summary(r.name()) {
                log::warn!("{s}");
            }
        }
        drop(readers);
        start += count;
    }
    let ranks = chroms.ranks();
    let total = match tie {
        Tie::InputOrder => cascade_merge_runs::<false>(runs, &mut spill, &ranks, fan, out)?,
        Tie::FullLine => cascade_merge_runs::<true>(runs, &mut spill, &ranks, fan, out)?,
    };
    spill.close();
    Ok(total)
}

/// k-way merge of already sorted inputs (`sort -m`), in `(chromosome order,
/// key, end, tiebreak)` order; ties go to the earlier input, or with
/// [`Tie::FullLine`] to the smaller line first. Every input is checked by its
/// own [`OrderGuard`] as it is read, batch by batch. More than 128 inputs
/// are merged in rounds through spill files in `tmp_dir` (`--tmp-dir`).
/// Returns the number of records written.
pub fn merge_sorted(
    readers: &mut [Box<dyn RecordReader>],
    chroms: &mut ChromTable,
    tie: Tie,
    tmp_dir: &Path,
    out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    if readers.len() <= MAX_FAN_IN {
        return match tie {
            Tie::InputOrder => merge_inputs::<false>(readers, chroms, 0, &mut |_, raw| out(raw)),
            Tie::FullLine => merge_inputs::<true>(readers, chroms, 0, &mut |_, raw| out(raw)),
        };
    }
    let fan = MAX_FAN_IN;
    let mut spill = SpillDir::create(tmp_dir)?;
    let mut runs = Vec::new();
    let mut start = 0;
    while start < readers.len() {
        let count = fan.min(readers.len() - start);
        let mut w = spill.create_run()?;
        let mut sink = |h: &RecHead, raw: &[u8]| w.write(h, raw).map_err(io::Error::other);
        let chunk = &mut readers[start..start + count];
        match tie {
            Tie::InputOrder => {
                merge_inputs::<false>(chunk, chroms, start as u64, &mut sink)?;
            }
            Tie::FullLine => {
                merge_inputs::<true>(chunk, chroms, start as u64, &mut sink)?;
            }
        }
        runs.push(w.finish()?);
        start += count;
    }
    let ranks = chroms.ranks();
    let total = match tie {
        Tie::InputOrder => cascade_merge_runs::<false>(runs, &mut spill, &ranks, fan, out)?,
        Tie::FullLine => cascade_merge_runs::<true>(runs, &mut spill, &ranks, fan, out)?,
    };
    spill.close();
    Ok(total)
}

fn merge_inputs<const FULL: bool>(
    readers: &mut [Box<dyn RecordReader>],
    chroms: &mut ChromTable,
    base_seq: u64,
    sink: &mut dyn FnMut(&RecHead, &[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    let mut inputs: Vec<Input> = readers
        .iter_mut()
        .map(|r| Input {
            guard: OrderGuard::new(r.name(), r.format()),
            reader: r.as_mut(),
            batch: RecordBatch::default(),
            pos: 0,
            more: true,
            // The first batch learns the record length.
            max: 1,
        })
        .collect();
    let mut live = Vec::with_capacity(inputs.len());
    for (src, input) in inputs.iter_mut().enumerate() {
        if input.refill(chroms)? {
            live.push(src);
        }
    }
    // Ranks are computed after reading, and again whenever reading interns a
    // new chromosome (which can shift the ranks of others).
    let mut ranks = chroms.ranks();
    let mut heap: BinaryHeap<Head<FULL>> = live
        .into_iter()
        .map(|src| inputs[src].head(src, base_seq, &ranks))
        .collect();
    let mut written = 0u64;
    loop {
        let Some(mut top) = heap.peek_mut() else {
            break;
        };
        let input = &mut inputs[top.src];
        sink(&top.rec, input.batch.raw(input.pos)).map_err(super::out_err)?;
        written += 1;
        input.pos += 1;
        if input.pos < input.batch.len() {
            input.update(&mut *top, &ranks);
            continue;
        }
        let src = top.src;
        PeekMut::pop(top);
        let input = &mut inputs[src];
        if input.refill(chroms)? {
            if ranks.len() != chroms.len() {
                ranks = chroms.ranks();
                let mut heads = std::mem::take(&mut heap).into_vec();
                for h in &mut heads {
                    h.rank = ranks[h.rec.chrom as usize];
                }
                heap = BinaryHeap::from(heads);
            }
            heap.push(input.head(src, base_seq, &ranks));
        }
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fan_in_bounds() {
        assert_eq!(fan_in(2 << 30, 100), MAX_FAN_IN);
        assert_eq!(fan_in(64 << 10, 40), MAX_FAN_IN);
        assert_eq!(fan_in(16 << 20, 3 << 20), 5);
        assert_eq!(fan_in(1 << 20, 3 << 20), 2);
        assert_eq!(fan_in(0, 0), 2);
    }
}
