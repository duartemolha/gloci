//! k-way merge of sorted record streams, shared by `merge` and `complement`.

use gloci_core::ChromTable;
use gloci_io::{IoError, OrderGuard, RecordBatch, RecordMeta, RecordReader};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Bytes of records to read from one input at a time.
const BATCH_BYTES: usize = 256 << 10;
const MAX_BATCH: usize = 4096;

/// Callback of [`kmerge`]: input index, record position data, raw bytes.
pub type Visit<'a> = dyn FnMut(usize, &RecordMeta, &[u8]) -> Result<(), IoError> + 'a;
/// Callback of [`kmerge_with`]: [`Visit`] plus the chromosome table.
pub type VisitWith<'a> =
    dyn FnMut(usize, &RecordMeta, &[u8], &ChromTable) -> Result<(), IoError> + 'a;

/// One input: its current batch and position.
struct Input {
    reader: Box<dyn RecordReader>,
    guard: OrderGuard,
    batch: RecordBatch,
    pos: usize,
    more: bool,
    /// Records to ask for per batch; adapted to the record length.
    max: usize,
}

impl Input {
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

    fn head(&self, src: usize, ranks: &[u32]) -> (u32, u64, usize) {
        let m = &self.batch.recs[self.pos];
        (ranks[m.chrom.0 as usize], m.key, src)
    }
}

/// Merge sorted inputs into one stream, calling `f(input index, meta, raw)`
/// for each record in `(chromosome order, key, input index, input order)`
/// order. Every input is order-checked against `chroms` as it is read; a
/// violation is an `OrderError` with the sort hint.
pub fn kmerge(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    f: &mut Visit,
) -> Result<(), IoError> {
    kmerge_with(inputs, chroms, &mut |i, m, raw, _| f(i, m, raw))
}

/// [`kmerge`] whose callback also gets the chromosome table as read so far,
/// to turn `meta.chrom` into a name.
pub fn kmerge_with(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    f: &mut VisitWith,
) -> Result<(), IoError> {
    let mut inputs: Vec<Input> = inputs
        .into_iter()
        .map(|r| Input {
            guard: OrderGuard::new(r.name(), r.format()),
            reader: r,
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
    // Ranks are recomputed whenever reading interns a new chromosome, which
    // can shift the ranks of others.
    let mut ranks = chroms.ranks();
    let mut heap: BinaryHeap<Reverse<(u32, u64, usize)>> = live
        .into_iter()
        .map(|src| Reverse(inputs[src].head(src, &ranks)))
        .collect();
    while let Some(Reverse((_, _, src))) = heap.pop() {
        let input = &mut inputs[src];
        f(
            src,
            &input.batch.recs[input.pos],
            input.batch.raw(input.pos),
            chroms,
        )?;
        input.pos += 1;
        let more = input.pos < input.batch.len() || input.refill(chroms)?;
        if !more {
            continue;
        }
        if ranks.len() != chroms.len() {
            ranks = chroms.ranks();
            // Stale ranks of the other heads: rebuild from each input's
            // current record.
            let others: Vec<usize> = heap.drain().map(|Reverse((_, _, s))| s).collect();
            for s in others {
                heap.push(Reverse(inputs[s].head(s, &ranks)));
            }
        }
        heap.push(Reverse(inputs[src].head(src, &ranks)));
    }
    Ok(())
}
