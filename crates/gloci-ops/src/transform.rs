//! `slop` and `flank`: grow each feature, or write the strips beside it
//! (replace `bedtools slop`/`flank` and bedops `--range`, D44, D49).
//!
//! Records stream in input order and need no sorted input, so no order
//! check applies to them. Output keeps the input format; only the
//! coordinate columns change. When a transform's output goes out of order
//! (`-s` or `--pct` in slop, any flank) one warning names the first
//! offending output record (input that is itself unsorted never triggers it).
//!
//! With a genome, a record that straddles the chromosome end is clipped, but
//! one that lies wholly beyond it is an error (bedtools writes nonsense for
//! it): the genome almost certainly is not the input's build.

use crate::relate::rewrite_coords;
use gloci_core::{ChromId, ChromTable, Fraction, Genome, Interval, Strand};
use gloci_io::{
    Format, IoError, Loc, RecordBatch, RecordReader, RecordWriter, UsageError, next_batch_size,
};

/// How far to pad: a number of bases, or a fraction of the feature length
/// (exact decimal, floored).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Amount {
    Bases(u64),
    Fraction(Fraction),
}

impl Amount {
    fn of(&self, len: u64) -> u64 {
        match self {
            Amount::Bases(n) => *n,
            Amount::Fraction(f) => f.floor_mul(len),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pad {
    pub left: Amount,
    pub right: Amount,
    /// Swap `left` and `right` for `-` strand records; an unknown strand
    /// counts as `+` (D49).
    pub strand_aware: bool,
}

/// `slop`/`flank` take BED and GFF/GTF only.
pub fn check_format(command: &str, format: Format) -> Result<(), UsageError> {
    match format {
        Format::Vcf | Format::Bcf => Err(UsageError(format!(
            "{command} works on BED and GFF/GTF input; a VCF/BCF record has no \
footprint to pad (convert it first, e.g. with gloci view -O bed)"
        ))),
        _ => Ok(()),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Slop,
    Flank,
}

struct Out<'a> {
    w: &'a mut dyn RecordWriter,
    format: Format,
    buf: Vec<u8>,
    /// Chromosome, input start and output start of the previous record.
    last: Option<(ChromId, u64, u64)>,
    lines: u64,
    warned: bool,
}

impl Out<'_> {
    fn emit(
        &mut self,
        chrom: ChromId,
        in_start: u64,
        raw: &[u8],
        iv: Interval,
    ) -> Result<(), IoError> {
        self.lines += 1;
        // Warn only when the transform itself breaks the order: same
        // chromosome, input in order, output start going down.
        if !self.warned
            && self
                .last
                .is_some_and(|(c, i, o)| c == chrom && i <= in_start && o > iv.start)
        {
            self.warned = true;
            log::warn!(
                "output record {} is out of order; the output of this command is not sorted \
(pipe the output through gloci sort)",
                self.lines
            );
        }
        self.last = Some((chrom, in_start, iv.start));
        self.buf.clear();
        rewrite_coords(self.format, raw, iv, &mut self.buf);
        self.w
            .write_record(&self.buf)
            .map_err(|source| IoError::Write {
                path: "output".into(),
                source,
            })
    }
}

fn run(
    mode: Mode,
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    pad: &Pad,
    genome: Option<&Genome>,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    let format = input.format();
    let name = input.name().to_string();
    check_format(if mode == Mode::Slop { "slop" } else { "flank" }, format)?;
    if let Some(g) = genome
        && let Some(i) = (0..g.len()).find(|&i| g.length(i).is_none())
    {
        return Err(UsageError(format!(
            "-g: the genome has no length for {}; clipping needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
            String::from_utf8_lossy(g.name(i))
        ))
        .into());
    }
    let mut out = Out {
        w,
        format,
        buf: Vec::new(),
        last: None,
        lines: 0,
        warned: false,
    };
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let raw = batch.raw(i);
            let len = m.span.len();
            let (mut l, mut r) = (pad.left.of(len), pad.right.of(len));
            if pad.strand_aware && m.strand == Strand::Minus {
                std::mem::swap(&mut l, &mut r);
            }
            // A closed table: the chromosome id is the genome index.
            let limit = genome.map_or(u64::MAX, |g| {
                g.length(m.chrom.0 as usize).expect("lengths checked")
            });
            let (s, e) = (m.span.start, m.span.end);
            if s > limit || (s == limit && e > s) {
                return Err(IoError::Parse {
                    path: name,
                    loc: Loc::Line(m.loc),
                    msg: format!(
                        "{}:{s}-{e} lies beyond the end of {} (length {limit}); \
check that -g matches the genome build of the input",
                        String::from_utf8_lossy(chroms.name(m.chrom)),
                        String::from_utf8_lossy(chroms.name(m.chrom)),
                    ),
                });
            }
            let iv = |a: u64, b: u64| Interval {
                start: a.min(limit),
                end: b.min(limit),
            };
            match mode {
                Mode::Slop => out.emit(
                    m.chrom,
                    s,
                    raw,
                    iv(s.saturating_sub(l), e.saturating_add(r)),
                )?,
                Mode::Flank => {
                    let left = iv(s.saturating_sub(l), s);
                    if !left.is_point() {
                        out.emit(m.chrom, s, raw, left)?;
                    }
                    let right = iv(e, e.saturating_add(r));
                    if !right.is_point() {
                        out.emit(m.chrom, s, raw, right)?;
                    }
                }
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            crate::skipped::log_reader(input.as_ref());
            return Ok(());
        }
    }
}

/// Write each record with `pad` added to its sides, clipped at 0 and, with a
/// genome, at the chromosome length (`[start - L, end + R)`).
pub fn slop(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    pad: &Pad,
    genome: Option<&Genome>,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    run(Mode::Slop, input, chroms, pad, genome, w)
}

/// Write the strips beside each record, left `[start - L, start)` then right
/// `[end, end + R)`, clipped like [`slop`]; a strip clipped to nothing is
/// not written.
pub fn flank(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    pad: &Pad,
    genome: Option<&Genome>,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    run(Mode::Flank, input, chroms, pad, genome, w)
}
