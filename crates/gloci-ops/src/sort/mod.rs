//! External merge sort with bounded memory and a deterministic order
//! (spec §7.4, D13, D24).
//!
//! Records are ordered by `(chromosome rank, key, end, tiebreak)`. The
//! tiebreak is input order ([`Tie::InputOrder`], a stable sort) or the raw
//! bytes and then input order ([`Tie::FullLine`], `--full-key`). Input order
//! is a record's offset in its chunk, and in a spilled run the chunk's byte
//! position plus that offset; it is total, so the output is the same for
//! every thread count.
//!
//! Records are copied into an in-memory chunk bounded by `max_mem`. If the
//! whole input fits in one chunk it is sorted and written without touching
//! disk. Otherwise each full chunk is sorted and spilled as an LZ4 run, and
//! the runs are merged up to 128 at a time, in cascading passes when there
//! are more (D24).
//!
//! Memory: the chunk (record bytes plus a 32-byte entry per record) stays
//! within `max_mem`, except that a single record larger than `max_mem` is
//! sorted in a chunk of its own. The chunk is freed before merging. Each run
//! open in a merge holds its current record, so the fan-in drops below 128
//! when 128 of the longest records would not fit in `max_mem` (very long VCF
//! lines); it also holds about 136 KiB of LZ4 and read buffers, a fixed
//! overhead of at most about 17 MiB. Chunk boundaries and the fan-in depend
//! only on `max_mem` and the input, never on the thread count, and the output
//! depends on neither.

mod chunk;
mod inputs;
mod merge;
mod order;
mod spill;

pub(crate) use inputs::coordinate_header;
pub use inputs::{
    HeaderRuleCmd, MIN_MAX_MEM, SortSettings, check_definitions, check_definitions_for,
    check_definitions_pair, check_samples, check_samples_for, check_samples_pair,
    sort_input_options, sort_readers, sort_together,
};
pub(crate) use merge::{MAX_FAN_IN, merge_alignment_sources};
pub use merge::{MERGE_BUFFER_PER_INPUT, merge_fan_in, merge_sorted, merge_sorted_batches};
pub use order::{Order, coordinate_key, name_cmp_raw, name_flag_rank, qname_of, strnum_cmp};

use chunk::Chunk;
use gloci_core::ChromTable;
use gloci_io::{IoError, RecordBatch, RecordMeta};
use order::{M_COORD, M_FULL, M_INTERVAL, M_NAME, M_NAME_BYTES};
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
    /// More than 1 sorts chunks in parallel on the current rayon pool (the run's own pool,
    /// which the command runs on). Never changes the output.
    pub threads: NonZeroUsize,
    /// How records are ordered; `Order::Interval(Tie::InputOrder)` for intervals.
    pub order: Order,
    /// Only for intervals; alignment orders take [`Dedup::All`].
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
    order: Order,
    dedup: Dedup,
    parallel: bool,
    budget: usize,
    chunk: Chunk,
    /// The input-order tiebreak of the current chunk's first record in the runs it spills:
    /// the record bytes of every earlier chunk plus one per chunk. A spilled record's tiebreak
    /// is this plus its arena offset, so tiebreaks rise in input order across runs and the
    /// runs of different chunks never share one (a run's own records are never compared).
    spill_base: u64,
    max_rec: usize,
    spills: u32,
    runs: Vec<Run>,
    spill: SpillDir,
    /// The run this sorter belongs to (D91).
    cancel: gloci_io::cancel::Watch,
}

impl Sorter {
    /// Creates the spill directory under `cfg.tmp_dir` up front, so an
    /// unusable temp directory fails before any input is read.
    pub fn new(cfg: SortConfig) -> Result<Sorter, IoError> {
        let budget = usize::try_from(cfg.max_mem).unwrap_or(usize::MAX);
        let order = match (cfg.dedup, cfg.order) {
            (Dedup::All, o) => o,
            (_, Order::Interval(_)) => Order::Interval(Tie::FullLine),
            (_, Order::Aln(_)) => {
                return Err(IoError::Internal(
                    "--unique and --duplicates do not apply to alignment sorting".into(),
                ));
            }
        };
        Ok(Sorter {
            order,
            dedup: cfg.dedup,
            parallel: cfg.threads.get() > 1,
            budget,
            chunk: Chunk::new(budget),
            spill_base: 0,
            max_rec: 0,
            spills: 0,
            runs: Vec::new(),
            spill: SpillDir::create(&cfg.tmp_dir)?,
            cancel: gloci_io::cancel::Watch::current(),
        })
    }

    /// Add a batch; `chroms` is the table its `ChromId`s come from (shared
    /// by all inputs). May sort and spill a full chunk.
    pub fn push(&mut self, chroms: &ChromTable, batch: &RecordBatch) -> Result<(), IoError> {
        self.push_keyed(chroms, batch, |m| (m.key, m.span.end))
    }

    /// [`Sorter::push`] with each record's `(key, end)` given by `keys` instead of its own
    /// `key` and `span.end`. `--sort` of alignments keys them by `(pos << 1) | reverse` with
    /// end 0, so they come out in `gloci sort`'s coordinate order (D121).
    pub fn push_keyed(
        &mut self,
        chroms: &ChromTable,
        batch: &RecordBatch,
        keys: impl Fn(&RecordMeta) -> (u64, u64),
    ) -> Result<(), IoError> {
        for m in &batch.recs {
            let raw = &batch.buf[m.off..m.off + m.len as usize];
            self.max_rec = self.max_rec.max(raw.len());
            let chrom = m.chrom.0;
            let (key, end) = keys(m);
            if !self.chunk.try_push(chrom, key, end, raw) {
                self.spill_chunk(&chroms.ranks())?;
                // Never drop a record silently.
                if !self.chunk.try_push(chrom, key, end, raw) {
                    return Err(IoError::Internal(format!(
                        "an empty sort chunk rejected a record of {} bytes",
                        raw.len()
                    )));
                }
            }
        }
        Ok(())
    }

    /// Add one BAM record (no `block_size`) of an alignment sort ([`Order::Aln`]); the input
    /// sequence number is the last tie. A record too short for its fixed fields is an error.
    pub fn push_raw(&mut self, raw: &[u8]) -> Result<(), IoError> {
        if !matches!(self.order, Order::Aln(_)) {
            return Err(IoError::Internal(
                "push_raw is for alignment sorting".into(),
            ));
        }
        if raw.len() < 32 {
            return Err(IoError::Internal(format!(
                "an alignment record of {} bytes is shorter than its 32 fixed bytes",
                raw.len()
            )));
        }
        self.max_rec = self.max_rec.max(raw.len());
        // Unplaced reads sort after every reference: u32::MAX is above every refID, as
        // `sort --check` and `sort -m` key them (see `coordinate_key`).
        let (chrom, key) = if self.order.mode() == M_COORD {
            coordinate_key(raw, u32::MAX)
        } else {
            (0, 0)
        };
        if !self.chunk.try_push(chrom, key, 0, raw) {
            self.spill_chunk(&[])?;
            if !self.chunk.try_push(chrom, key, 0, raw) {
                return Err(IoError::Internal(format!(
                    "an empty sort chunk rejected a record of {} bytes",
                    raw.len()
                )));
            }
        }
        Ok(())
    }

    /// Sort the chunk and write it as a run.
    fn spill_chunk(&mut self, ranks: &[u32]) -> Result<(), IoError> {
        self.cancel.check()?;
        self.chunk.sort(ranks, self.order, self.parallel);
        let mut w = self.spill.create_run()?;
        for e in &self.chunk.entries {
            self.cancel.check()?;
            let h = RecHead {
                chrom: e.chrom,
                key: e.key,
                end: e.end,
                seq: self.spill_base + e.off as u64,
            };
            w.write(&h, self.chunk.raw(e))?;
        }
        self.runs.push(w.finish()?);
        self.spills += 1;
        self.spill_base += self.chunk.arena.len() as u64 + 1;
        self.chunk.clear();
        Ok(())
    }

    /// Write every record in order to `out` and remove the spill directory.
    pub fn finish(
        self,
        chroms: &ChromTable,
        out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
    ) -> Result<SortStats, IoError> {
        // Every chromosome of the runs is in `chroms` by now.
        self.finish_ranked(&chroms.ranks(), out)
    }

    /// As [`Sorter::finish`] for an alignment sort, which has no chromosome table.
    pub fn finish_raw(
        self,
        out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
    ) -> Result<SortStats, IoError> {
        self.finish_ranked(&[], out)
    }

    fn finish_ranked(
        mut self,
        ranks: &[u32],
        out: &mut dyn FnMut(&[u8]) -> io::Result<()>,
    ) -> Result<SortStats, IoError> {
        let mut emit = Emitter::new(self.dedup, out);
        if self.runs.is_empty() {
            self.cancel.check()?;
            self.chunk.sort(ranks, self.order, self.parallel);
            for e in &self.chunk.entries {
                self.cancel.check()?;
                emit.emit(self.chunk.raw(e))?;
            }
        } else {
            if !self.chunk.is_empty() {
                self.spill_chunk(ranks)?;
            }
            self.chunk.release();
            match self.order.mode() {
                M_INTERVAL => self.merge_all::<M_INTERVAL>(ranks, &mut emit)?,
                M_FULL => self.merge_all::<M_FULL>(ranks, &mut emit)?,
                M_COORD => self.merge_all::<M_COORD>(ranks, &mut emit)?,
                M_NAME => self.merge_all::<M_NAME>(ranks, &mut emit)?,
                M_NAME_BYTES => self.merge_all::<M_NAME_BYTES>(ranks, &mut emit)?,
                m => {
                    return Err(IoError::Internal(format!("unknown sort mode {m}")));
                }
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
    fn merge_all<const M: u8>(&mut self, ranks: &[u32], emit: &mut Emitter) -> Result<(), IoError> {
        let fan = merge::fan_in(self.budget, self.max_rec);
        let mut runs: VecDeque<Run> = std::mem::take(&mut self.runs).into();
        while runs.len() > fan {
            // Merge just enough runs that the rest can follow `fan` at a time.
            let k = fan.min(runs.len() - fan + 1);
            let group: Vec<Run> = runs.drain(..k).collect();
            let mut w = self.spill.create_run()?;
            merge::merge_runs::<M>(&group, ranks, &mut |h, raw| w.write(h, raw))?;
            runs.push_back(w.finish()?);
            group.into_iter().for_each(Run::remove);
        }
        merge::merge_runs::<M>(runs.make_contiguous(), ranks, &mut |_, raw| emit.emit(raw))
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
