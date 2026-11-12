//! k-way merges: of spilled runs (with cascading passes, D24) and of
//! already sorted inputs (`sort -m`).

use super::Tie;
use super::order::{
    M_COORD, M_FULL, M_INTERVAL, M_NAME, M_NAME_BYTES, coordinate_key, keeps_bytes, name_cmp_raw,
};
use super::spill::{RecHead, Run, RunReader, SpillDir};
use crate::aln_sort::{BamRecords, cmp as aln_cmp, order_name};
use gloci_core::ChromTable;
use gloci_io::{AlnSortOrder, AlnSource, IoError, Loc, OrderGuard, RecordBatch, RecordReader};
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
struct Head<const M: u8> {
    rank: u32,
    rec: RecHead,
    src: usize,
    bytes: Vec<u8>,
}

impl<const M: u8> Head<M> {
    fn cmp_key(&self, other: &Self) -> Ordering {
        let pos = |a: &Self, b: &Self| {
            (a.rank, a.rec.key, a.rec.end).cmp(&(b.rank, b.rec.key, b.rec.end))
        };
        let ord = match M {
            M_INTERVAL | M_COORD => pos(self, other),
            M_FULL => pos(self, other).then_with(|| self.bytes.cmp(&other.bytes)),
            _ => name_cmp_raw(&self.bytes, &other.bytes, M == M_NAME),
        };
        ord.then(self.rec.seq.cmp(&other.rec.seq))
    }
}

impl<const M: u8> Ord for Head<M> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_key(other).reverse()
    }
}

impl<const M: u8> PartialOrd for Head<M> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const M: u8> PartialEq for Head<M> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<const M: u8> Eq for Head<M> {}

fn rank_of<const M: u8>(ranks: &[u32], chrom: u32) -> Result<u32, IoError> {
    if M >= M_COORD {
        // Alignment orders compare reference ids as they are.
        return Ok(chrom);
    }
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
pub(crate) fn merge_runs<const M: u8>(
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
            let rank = rank_of::<M>(ranks, rec.chrom)?;
            heap.push(Head::<M> {
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
        gloci_io::cancel::check()?;
        sink(&top.rec, &top.bytes)?;
        let h = &mut *top;
        if readers[h.src].next(&mut h.rec, &mut h.bytes)? {
            // Replaced in place; dropping the PeekMut sifts it down.
            h.rank = rank_of::<M>(ranks, h.rec.chrom)?;
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

    fn head<const M: u8>(&self, src: usize, base_seq: u64, ranks: &[u32]) -> Head<M> {
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
            bytes: if keeps_bytes(M) {
                self.batch.raw(self.pos).to_vec()
            } else {
                Vec::new()
            },
        }
    }

    /// Refresh `h` from the current record (`h.src` is unchanged).
    fn update<const M: u8>(&self, h: &mut Head<M>, ranks: &[u32]) {
        let m = &self.batch.recs[self.pos];
        h.rank = ranks[m.chrom.0 as usize];
        h.rec.chrom = m.chrom.0;
        h.rec.key = m.key;
        h.rec.end = m.span.end;
        if keeps_bytes(M) {
            h.bytes.clear();
            h.bytes.extend_from_slice(self.batch.raw(self.pos));
        }
    }
}

/// Cascade-merge `runs` down to `out`, merging up to `fan` runs at once.
pub(crate) fn cascade_merge_runs<const M: u8>(
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
        merge_runs::<M>(&group, ranks, &mut |h, raw| w.write(h, raw))?;
        runs.push_back(w.finish()?);
        group.into_iter().for_each(Run::remove);
    }
    merge_runs::<M>(runs.make_contiguous(), ranks, &mut |_, raw| {
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
                merge_inputs::<M_INTERVAL>(&mut readers, chroms, start as u64, &mut sink)?;
            }
            Tie::FullLine => {
                merge_inputs::<M_FULL>(&mut readers, chroms, start as u64, &mut sink)?;
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
        Tie::InputOrder => cascade_merge_runs::<M_INTERVAL>(runs, &mut spill, &ranks, fan, out)?,
        Tie::FullLine => cascade_merge_runs::<M_FULL>(runs, &mut spill, &ranks, fan, out)?,
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
            Tie::InputOrder => {
                merge_inputs::<M_INTERVAL>(readers, chroms, 0, &mut |_, raw| out(raw))
            }
            Tie::FullLine => merge_inputs::<M_FULL>(readers, chroms, 0, &mut |_, raw| out(raw)),
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
                merge_inputs::<M_INTERVAL>(chunk, chroms, start as u64, &mut sink)?;
            }
            Tie::FullLine => {
                merge_inputs::<M_FULL>(chunk, chroms, start as u64, &mut sink)?;
            }
        }
        runs.push(w.finish()?);
        start += count;
    }
    let ranks = chroms.ranks();
    let total = match tie {
        Tie::InputOrder => cascade_merge_runs::<M_INTERVAL>(runs, &mut spill, &ranks, fan, out)?,
        Tie::FullLine => cascade_merge_runs::<M_FULL>(runs, &mut spill, &ranks, fan, out)?,
    };
    spill.close();
    Ok(total)
}

fn merge_inputs<const M: u8>(
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
    let mut heap: BinaryHeap<Head<M>> = live
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

/// The sort key of an alignment record in mode `M`: `(tid, key)` for coordinate order, and
/// nothing for the name orders, which compare the bytes.
fn aln_key<const M: u8>(raw: &[u8]) -> (u32, u64) {
    if M == M_COORD {
        coordinate_key(raw, u32::MAX)
    } else {
        (0, 0)
    }
}

fn aln_head<const M: u8>(raw: &[u8], src: usize, seq: u64) -> Head<M> {
    let (rank, key) = aln_key::<M>(raw);
    Head {
        rank,
        rec: RecHead {
            chrom: rank,
            key,
            end: 0,
            seq,
        },
        src,
        bytes: raw.to_vec(),
    }
}

/// k-way merge of alignment inputs, each already in `order`, holding one record per input. Ties
/// go to the earlier input (`base_seq + index`), then to record order, as samtools' `heap_lt`.
/// Each input is checked as it is read; a record that comes before its predecessor is an error
/// naming the file and the record.
fn merge_alignment_inputs<const M: u8>(
    inputs: &mut [BamRecords],
    order: AlnSortOrder,
    base_seq: u64,
    sink: &mut dyn FnMut(&RecHead, &[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    let mut heap: BinaryHeap<Head<M>> = BinaryHeap::with_capacity(inputs.len());
    // The number of the record each input is at.
    let mut at = vec![0u64; inputs.len()];
    for (src, input) in inputs.iter_mut().enumerate() {
        if let Some(n) = input.next()? {
            at[src] = n;
            heap.push(aln_head::<M>(input.bam(), src, base_seq + src as u64));
        }
    }
    let mut written = 0u64;
    loop {
        let Some(mut top) = heap.peek_mut() else {
            break;
        };
        sink(&top.rec, &top.bytes).map_err(super::out_err)?;
        written += 1;
        let src = top.src;
        let input = &mut inputs[src];
        let Some(n) = input.next()? else {
            PeekMut::pop(top);
            continue;
        };
        if aln_cmp(order, &top.bytes, input.bam()) == Ordering::Greater {
            return Err(unsorted_input(
                input.src.name(),
                order,
                (n, input.bam()),
                (at[src], &top.bytes),
            ));
        }
        at[src] = n;
        let (rank, key) = aln_key::<M>(input.bam());
        let h = &mut *top;
        h.rank = rank;
        h.rec.chrom = rank;
        h.rec.key = key;
        h.bytes.clear();
        h.bytes.extend_from_slice(input.bam());
        // Dropping the PeekMut sifts the head down.
    }
    Ok(written)
}

fn unsorted_input(
    name: &str,
    order: AlnSortOrder,
    (n, cur): (u64, &[u8]),
    (prev_n, prev): (u64, &[u8]),
) -> IoError {
    let label = |raw: &[u8]| String::from_utf8_lossy(super::qname_of(raw)).into_owned();
    IoError::Parse {
        path: name.to_string(),
        loc: Loc::Record(n),
        msg: format!(
            "the file is not in {} order: this record ({}) comes after record {prev_n} ({}); \
sort it first with: gloci sort {name} -o sorted.bam{}",
            order_name(order),
            label(cur),
            label(prev),
            match order {
                AlnSortOrder::Coordinate => "",
                AlnSortOrder::Name => " -n",
                AlnSortOrder::NameBytes => " -N",
            }
        ),
    }
}

/// Opens alignment input `i` (counted from 0).
pub(crate) type OpenAln<'a> = dyn FnMut(usize) -> Result<Box<dyn AlnSource>, IoError> + 'a;

/// Merge the `n` alignment inputs that `open` opens, each in `order`, into `out` as BAM record
/// bytes. At most `fan` are open at a time: more are merged in rounds, each round's inputs
/// opened when it starts and closed when it ends, into spill runs in `tmp_dir` that are then
/// cascade-merged. Returns the records written.
pub(crate) fn merge_alignment_sources(
    n: usize,
    open: &mut OpenAln,
    order: AlnSortOrder,
    fan: usize,
    tmp_dir: &Path,
    out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    match order {
        AlnSortOrder::Coordinate => {
            merge_alignments_in::<M_COORD>(n, open, order, fan, tmp_dir, out)
        }
        AlnSortOrder::Name => merge_alignments_in::<M_NAME>(n, open, order, fan, tmp_dir, out),
        AlnSortOrder::NameBytes => {
            merge_alignments_in::<M_NAME_BYTES>(n, open, order, fan, tmp_dir, out)
        }
    }
}

fn open_round(open: &mut OpenAln, start: usize, count: usize) -> Result<Vec<BamRecords>, IoError> {
    (start..start + count)
        .map(|i| open(i).map(BamRecords::new))
        .collect()
}

fn merge_alignments_in<const M: u8>(
    n: usize,
    open: &mut OpenAln,
    order: AlnSortOrder,
    fan: usize,
    tmp_dir: &Path,
    out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
) -> Result<u64, IoError> {
    let fan = fan.clamp(2, MAX_FAN_IN);
    if n <= fan {
        let mut inputs = open_round(open, 0, n)?;
        let total = merge_alignment_inputs::<M>(&mut inputs, order, 0, &mut |_, raw| out(raw))?;
        inputs.into_iter().for_each(BamRecords::finish);
        return Ok(total);
    }
    let mut spill = SpillDir::create(tmp_dir)?;
    let mut runs = Vec::new();
    let mut start = 0;
    while start < n {
        let count = fan.min(n - start);
        let mut inputs = open_round(open, start, count)?;
        let mut w = spill.create_run()?;
        merge_alignment_inputs::<M>(&mut inputs, order, start as u64, &mut |h, raw| {
            w.write(h, raw).map_err(io::Error::other)
        })?;
        runs.push(w.finish()?);
        inputs.into_iter().for_each(BamRecords::finish);
        start += count;
    }
    let total = cascade_merge_runs::<M>(runs, &mut spill, &[], fan, out)?;
    spill.close();
    Ok(total)
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
