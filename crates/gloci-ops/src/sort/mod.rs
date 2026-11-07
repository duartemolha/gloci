//! External merge sort with bounded memory and a deterministic order
//! (spec §7.4, D13, D24).
//!
//! Records are ordered by `(chromosome rank, key, end, tiebreak)`. The
//! tiebreak is the global input sequence number ([`Tie::InputOrder`], a
//! stable sort) or the raw bytes and then the sequence number
//! ([`Tie::FullLine`], `--full-key`). Sequence numbers are unique, so the
//! order is total and the output is the same for every thread count.
//!
//! Records are copied into an in-memory chunk bounded by `max_mem`. If the
//! whole input fits in one chunk it is sorted and written without touching
//! disk. Otherwise each full chunk is sorted and spilled as an LZ4 run, and
//! the runs are merged up to 128 at a time, in cascading passes when there
//! are more (D24).
//!
//! Memory: the chunk (record bytes plus a 40-byte entry per record) stays
//! within `max_mem`, except that a single record larger than `max_mem` is
//! sorted in a chunk of its own. The chunk is freed before merging. Each run
//! open in a merge holds its current record, so the fan-in drops below 128
//! when 128 of the longest records would not fit in `max_mem` (very long VCF
//! lines); it also holds about 136 KiB of LZ4 and read buffers, a fixed
//! overhead of at most about 17 MiB. Chunk boundaries and the fan-in depend
//! only on `max_mem` and the input, never on the thread count, and the output
//! depends on neither.

mod chunk;
mod merge;
mod spill;

pub use merge::merge_sorted;

use chunk::Chunk;
use gloci_core::ChromTable;
use gloci_io::{IoError, RecordBatch};
use spill::{RecHead, Run, SpillDir};
use std::collections::VecDeque;
use std::io;
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// How records with equal `(chromosome, key, end)` are ordered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tie {
    /// Input order (a stable sort).
    InputOrder,
    /// By the whole raw record, then input order (`--full-key`).
    FullLine,
}

/// Which records to write, comparing whole raw records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dedup {
    All,
    /// One copy of each distinct record (`uniq`); implies [`Tie::FullLine`].
    Unique,
    /// One copy of each record that occurs more than once (`uniq -d`);
    /// implies [`Tie::FullLine`].
    Duplicates,
}

#[derive(Debug, Clone)]
pub struct SortConfig {
    /// Memory budget in bytes for records held in memory. The library
    /// accepts any value; a record larger than the budget is sorted in a
    /// chunk of its own.
    pub max_mem: u64,
    /// Where the spill directory is created.
    pub tmp_dir: PathBuf,
    /// More than 1 sorts chunks in parallel on the global rayon pool. Never
    /// changes the output.
    pub threads: NonZeroUsize,
    pub tie: Tie,
    pub dedup: Dedup,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SortStats {
    /// Records written to the output.
    pub records: u64,
    /// Chunks spilled to disk (0 when the input fit in memory).
    pub spills: u32,
}

/// Name used for errors returned by the output callback.
const OUTPUT: &str = "output";

/// An error from the output callback, its `io::ErrorKind` kept (so a broken
/// pipe stays recognisable).
pub(crate) fn out_err(source: io::Error) -> IoError {
    IoError::Write {
        path: OUTPUT.to_string(),
        source,
    }
}

/// Sorts records pushed batch by batch; see the module docs.
pub struct Sorter {
    tie: Tie,
    dedup: Dedup,
    parallel: bool,
    budget: usize,
    chunk: Chunk,
    seq: u64,
    max_rec: usize,
    spills: u32,
    runs: Vec<Run>,
    spill: SpillDir,
}

impl Sorter {
    /// Creates the spill directory under `cfg.tmp_dir` up front, so an
    /// unusable temp directory fails before any input is read.
    pub fn new(cfg: SortConfig) -> Result<Sorter, IoError> {
        let budget = usize::try_from(cfg.max_mem).unwrap_or(usize::MAX);
        let tie = match cfg.dedup {
            Dedup::All => cfg.tie,
            Dedup::Unique | Dedup::Duplicates => Tie::FullLine,
        };
        Ok(Sorter {
            tie,
            dedup: cfg.dedup,
            parallel: cfg.threads.get() > 1,
            budget,
            chunk: Chunk::new(budget),
            seq: 0,
            max_rec: 0,
            spills: 0,
            runs: Vec::new(),
            spill: SpillDir::create(&cfg.tmp_dir)?,
        })
    }

    /// Add a batch; `chroms` is the table its `ChromId`s come from (shared
    /// by all inputs). May sort and spill a full chunk.
    pub fn push(&mut self, chroms: &ChromTable, batch: &RecordBatch) -> Result<(), IoError> {
        for m in &batch.recs {
            let raw = &batch.buf[m.off..m.off + m.len as usize];
            self.max_rec = self.max_rec.max(raw.len());
            let (chrom, end) = (m.chrom.0, m.span.end);
            if !self.chunk.try_push(chrom, m.key, end, self.seq, raw) {
                self.spill_chunk(chroms)?;
                let pushed = self.chunk.try_push(chrom, m.key, end, self.seq, raw);
                // Never drop a record silently.
                assert!(pushed, "an empty chunk accepts any record");
            }
            self.seq += 1;
        }
        Ok(())
    }

    /// Sort the chunk and write it as a run.
    fn spill_chunk(&mut self, chroms: &ChromTable) -> Result<(), IoError> {
        self.chunk.sort(&chroms.ranks(), self.tie, self.parallel);
        let mut w = self.spill.create_run()?;
        for e in &self.chunk.entries {
            let h = RecHead {
                chrom: e.chrom,
                key: e.key,
                end: e.end,
                seq: e.seq,
            };
            w.write(&h, self.chunk.raw(e))?;
        }
        self.runs.push(w.finish()?);
        self.spills += 1;
        self.chunk.clear();
        Ok(())
    }

    /// Write every record in order to `out` and remove the spill directory.
    pub fn finish(
        mut self,
        chroms: &ChromTable,
        out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
    ) -> Result<SortStats, IoError> {
        let mut emit = Emitter::new(self.dedup, out);
        if self.runs.is_empty() {
            self.chunk.sort(&chroms.ranks(), self.tie, self.parallel);
            for e in &self.chunk.entries {
                emit.emit(self.chunk.raw(e))?;
            }
        } else {
            if !self.chunk.is_empty() {
                self.spill_chunk(chroms)?;
            }
            self.chunk.release();
            // Every chromosome of the runs is in `chroms` by now.
            let ranks = chroms.ranks();
            match self.tie {
                Tie::InputOrder => self.merge_all::<false>(&ranks, &mut emit)?,
                Tie::FullLine => self.merge_all::<true>(&ranks, &mut emit)?,
            }
        }
        let stats = SortStats {
            records: emit.written,
            spills: self.spills,
        };
        self.spill.close();
        Ok(stats)
    }

    /// Merge the runs into `emit`, first merging groups of runs into new
    /// runs while there are more than the fan-in.
    fn merge_all<const FULL: bool>(
        &mut self,
        ranks: &[u32],
        emit: &mut Emitter,
    ) -> Result<(), IoError> {
        let fan = merge::fan_in(self.budget, self.max_rec);
        let mut runs: VecDeque<Run> = std::mem::take(&mut self.runs).into();
        while runs.len() > fan {
            // Merge just enough runs that the rest can follow `fan` at a time.
            let k = fan.min(runs.len() - fan + 1);
            let group: Vec<Run> = runs.drain(..k).collect();
            let mut w = self.spill.create_run()?;
            merge::merge_runs::<FULL>(&group, ranks, &mut |h, raw| w.write(h, raw))?;
            runs.push_back(w.finish()?);
            group.into_iter().for_each(Run::remove);
        }
        merge::merge_runs::<FULL>(runs.make_contiguous(), ranks, &mut |_, raw| emit.emit(raw))
    }
}

/// Writes sorted records to the output, applying `--unique` /
/// `--duplicates` to whole raw records (identical records are adjacent in
/// the full-line order).
struct Emitter<'a> {
    dedup: Dedup,
    out: &'a mut dyn FnMut(&[u8]) -> io::Result<()>,
    prev: Vec<u8>,
    /// Copies of `prev` seen so far (0 before the first record).
    copies: u64,
    written: u64,
}

impl<'a> Emitter<'a> {
    fn new(dedup: Dedup, out: &'a mut dyn FnMut(&[u8]) -> io::Result<()>) -> Self {
        Emitter {
            dedup,
            out,
            prev: Vec::new(),
            copies: 0,
            written: 0,
        }
    }

    fn write(&mut self, rec: &[u8]) -> Result<(), IoError> {
        (self.out)(rec).map_err(out_err)?;
        self.written += 1;
        Ok(())
    }

    fn emit(&mut self, rec: &[u8]) -> Result<(), IoError> {
        if self.dedup == Dedup::All {
            return self.write(rec);
        }
        if self.copies > 0 && rec == self.prev.as_slice() {
            self.copies += 1;
            if self.dedup == Dedup::Duplicates && self.copies == 2 {
                return self.write(rec);
            }
            return Ok(());
        }
        self.prev.clear();
        self.prev.extend_from_slice(rec);
        self.copies = 1;
        if self.dedup == Dedup::Unique {
            return self.write(rec);
        }
        Ok(())
    }
}
