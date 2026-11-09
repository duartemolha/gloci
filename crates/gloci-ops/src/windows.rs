//! `windows`: cut each chromosome (`-g`) or each interval (`-i`) into
//! windows (replaces `bedtools makewindows` and `bedops --chop`).
//!
//! Semantics follow bedtools: `Size` windows start every `step` bases and the
//! last one is truncated at the interval end (a step larger than the size
//! leaves gaps, a smaller one overlaps); `Count` gives `count` windows of
//! `floor(len / count)` bases, the last one taking the remainder. An interval
//! shorter than `count` gets no windows (bedtools warns and skips it too).
//! The id column is the 1-based window number, counted from the end with
//! `reverse`. Points have no bases and give no windows.

use gloci_core::{ChromTable, Genome};
use gloci_io::{IoError, RecordBatch, RecordReader, RecordWriter, UsageError, next_batch_size};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Split {
    /// Windows of `size` bases starting every `step` bases.
    Size { size: u64, step: u64 },
    /// This many windows per interval.
    Count(u64),
}

#[derive(Debug, Clone, Copy)]
pub struct WindowsOptions {
    pub split: Split,
    /// Number the windows from the end of the interval.
    pub reverse: bool,
    /// Append the 1-based window number.
    pub id_column: bool,
}

impl WindowsOptions {
    fn check(&self) -> Result<(), UsageError> {
        match self.split {
            Split::Size { size: 0, .. } => {
                Err(UsageError("-w: the window size must be at least 1".into()))
            }
            Split::Size { step: 0, .. } => Err(UsageError("--step must be at least 1".into())),
            Split::Count(0) => Err(UsageError("-n: the window count must be at least 1".into())),
            _ => Ok(()),
        }
    }
}

struct Cutter<'a> {
    o: &'a WindowsOptions,
    w: &'a mut dyn RecordWriter,
    row: Vec<u8>,
    /// Intervals shorter than the window count: how many, and the first.
    short: u64,
    first_short: Option<String>,
}

impl Cutter<'_> {
    /// One summary warning for the intervals `-n` had to skip.
    fn finish(&self) {
        if let (Some(first), Split::Count(c)) = (&self.first_short, self.o.split) {
            log::warn!(
                "{} interval(s) shorter than the {c} windows requested were skipped (first: {first})",
                self.short
            );
        }
    }

    fn emit(&mut self, name: &[u8], s: u64, e: u64, id: u64) -> Result<(), IoError> {
        self.row.clear();
        self.row.extend_from_slice(name);
        for v in [s, e] {
            self.row.push(b'\t');
            self.row.extend_from_slice(v.to_string().as_bytes());
        }
        if self.o.id_column {
            self.row.push(b'\t');
            self.row.extend_from_slice(id.to_string().as_bytes());
        }
        self.w
            .write_record(&self.row)
            .map_err(|source| IoError::Write {
                path: "output".into(),
                source,
            })
    }

    /// Cut `[start, end)` of the chromosome `name`.
    fn cut(&mut self, name: &[u8], start: u64, end: u64) -> Result<(), IoError> {
        let len = end - start;
        if len == 0 {
            return Ok(());
        }
        let (n, size, step) = match self.o.split {
            Split::Size { size, step } => (len.div_ceil(step), size, step),
            Split::Count(c) => {
                if len < c {
                    self.short += 1;
                    if self.first_short.is_none() {
                        self.first_short =
                            Some(format!("{}:{start}-{end}", String::from_utf8_lossy(name)));
                    }
                    return Ok(());
                }
                (c, len / c, len / c)
            }
        };
        for i in 0..n {
            let s = start + i * step;
            let e = if matches!(self.o.split, Split::Count(_)) && i + 1 == n {
                end
            } else {
                s.saturating_add(size).min(end)
            };
            let id = if self.o.reverse { n - i } else { i + 1 };
            self.emit(name, s, e, id)?;
        }
        Ok(())
    }
}

/// Windows over every chromosome of `genome`, in genome order.
pub fn windows_genome(
    genome: &Genome,
    o: &WindowsOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    o.check()?;
    let mut c = Cutter {
        o,
        w,
        row: Vec::new(),
        short: 0,
        first_short: None,
    };
    for i in 0..genome.len() {
        let Some(len) = genome.length(i) else {
            return Err(UsageError(format!(
                "-g: the genome has no length for {}; windows needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
                String::from_utf8_lossy(genome.name(i))
            ))
            .into());
        };
        if len == 0 {
            log::warn!(
                "{} has length 0 in the genome; no windows written for it",
                String::from_utf8_lossy(genome.name(i))
            );
        }
        c.cut(genome.name(i), 0, len)?;
    }
    c.finish();
    Ok(())
}

/// Windows over the footprint of every record of `input`, in input order
/// (the input needs no sorting).
pub fn windows_intervals(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &WindowsOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    o.check()?;
    let mut c = Cutter {
        o,
        w,
        row: Vec::new(),
        short: 0,
        first_short: None,
    };
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for m in &batch.recs {
            c.cut(chroms.name(m.chrom), m.span.start, m.span.end)?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            c.finish();
            crate::skipped::log_reader(input.as_ref());
            return Ok(());
        }
    }
}
