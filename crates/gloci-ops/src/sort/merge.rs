//! k-way merges: of spilled runs (with cascading passes, D24) and of
//! already sorted inputs (`sort -m`).

use super::Tie;
use super::spill::{RecHead, Run, RunReader};
use gloci_core::ChromTable;
use gloci_io::{IoError, OrderGuard, RecordBatch, RecordReader};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::binary_heap::PeekMut;
use std::io;

/// The most runs merged at once (D24).
pub(crate) const MAX_FAN_IN: usize = 128;

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

    fn head<const FULL: bool>(&self, src: usize, ranks: &[u32]) -> Head<FULL> {
        let m = &self.batch.recs[self.pos];
        Head {
            rank: ranks[m.chrom.0 as usize],
            rec: RecHead {
                chrom: m.chrom.0,
                key: m.key,
                end: m.span.end,
                // Ties between inputs go to the earlier input.
                seq: src as u64,
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

/// k-way merge of already sorted inputs (`sort -m`), in `(chromosome order,
/// key, end, tiebreak)` order; ties go to the earlier input, or with
/// [`Tie::FullLine`] to the smaller line first. Every input is checked by its
/// own [`OrderGuard`] as it is read, batch by batch. Returns the number of
/// records written.
pub fn merge_sorted(
    readers: &mut [Box<dyn RecordReader>],
    chroms: &mut ChromTable,
    tie: Tie,
    out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    match tie {
        Tie::InputOrder => merge_inputs::<false>(readers, chroms, out),
        Tie::FullLine => merge_inputs::<true>(readers, chroms, out),
    }
}

fn merge_inputs<const FULL: bool>(
    readers: &mut [Box<dyn RecordReader>],
    chroms: &mut ChromTable,
    out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
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
        .map(|src| inputs[src].head(src, &ranks))
        .collect();
    let mut written = 0u64;
    loop {
        let Some(mut top) = heap.peek_mut() else {
            break;
        };
        let input = &mut inputs[top.src];
        out(input.batch.raw(input.pos)).map_err(super::out_err)?;
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
            heap.push(input.head(src, &ranks));
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
