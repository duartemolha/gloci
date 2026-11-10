//! `genomecov`: the depth of coverage of a genome by the features of one
//! sorted BED-family input (replaces `bedtools genomecov -i`; BAM is P3).
//!
//! The input is swept once. Each feature adds one to the depth over its
//! bases (over each BED12 block with `split`, or over its single 5' or 3'
//! base). The sweep keeps only the depth-change events of the features
//! seen so far that are not yet behind the current start, so memory follows
//! the nesting depth, not the file size. Depth runs are streamed out as they
//! close: bedGraph rows and per-base rows are written immediately, and the
//! histogram keeps one counter per depth value (bounded by the greatest
//! depth, or by `max_depth`).
//!
//! The runs cover whole chromosomes from base 0 to the chromosome length
//! (`-g` is required), in genome order, chromosomes without features
//! included. Points cover no bases (D68) and are ignored.

use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable, Genome, Interval, Strand};
use gloci_engine::kmerge_with;
use gloci_io::bed::bed12_blocks;
use gloci_io::{Format, IoError, Loc, RecordMeta, RecordReader, RecordWriter, UsageError};
use std::collections::BTreeMap;
use std::io::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GenomecovMode {
    /// Per chromosome and for the whole genome: bases at each depth.
    #[default]
    Hist,
    /// Runs of equal depth above zero.
    BedGraph,
    /// Runs of equal depth including zero.
    BedGraphAll,
    /// One row per base (1-based position).
    PerBase,
}

/// Which part of each feature is counted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Piece {
    #[default]
    Whole,
    /// The first base in the feature's direction (strand-aware).
    FivePrime,
    /// The last base in the feature's direction (strand-aware).
    ThreePrime,
}

#[derive(Debug, Clone, Default)]
pub struct GenomecovOptions {
    pub mode: GenomecovMode,
    /// Multiply every reported depth (not for `Hist`).
    pub scale: Option<f64>,
    /// Count only features on this strand (`Plus` or `Minus`).
    pub strand: Option<Strand>,
    /// Count each BED12 block instead of the whole feature.
    pub split: bool,
    pub piece: Piece,
    /// `Hist` only: depths at or above this go into its bin.
    pub max_depth: Option<u64>,
}

pub fn validate(o: &GenomecovOptions) -> Result<(), UsageError> {
    let u = |m: &str| Err(UsageError(m.to_string()));
    if o.scale.is_some() && o.mode == GenomecovMode::Hist {
        return u(
            "--scale needs --bedgraph, --bga or --per-base: the histogram counts bases, \
which a scale factor would make meaningless; pick one of those modes",
        );
    }
    if let Some(s) = o.scale
        && !s.is_finite()
    {
        return u("--scale must be a finite number, for example 0.5 or 1e6");
    }
    if o.max_depth.is_some() && o.mode != GenomecovMode::Hist {
        return u(
            "--max-depth only applies to the histogram (the default mode); \
drop it, or drop --bedgraph/--bga/--per-base",
        );
    }
    if o.split && o.piece != Piece::Whole {
        return u("--split cannot be combined with --5p or --3p; use one of them");
    }
    if !matches!(o.strand, None | Some(Strand::Plus | Strand::Minus)) {
        return u("--strand must be + or -");
    }
    Ok(())
}

/// A number as C++ `ostream` prints a `double` by default (`%g`, six
/// significant digits), which is how bedtools writes scaled depths and
/// fractions: `0.6`, `1`, `0.0704225`, `1e-07`, `1e+08`.
pub fn format_g6(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    if !x.is_finite() {
        return if x.is_nan() {
            "nan".into()
        } else if x < 0.0 {
            "-inf".into()
        } else {
            "inf".into()
        };
    }
    let sci = format!("{x:.5e}");
    let (mant, exp) = sci.split_once('e').expect("exponent form");
    let exp: i32 = exp.parse().expect("exponent");
    let trim = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if !(-4..6).contains(&exp) {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mant), exp.abs())
    } else {
        trim(&format!("{x:.*}", (5 - exp) as usize))
    }
}

struct Cov<'a> {
    o: &'a GenomecovOptions,
    genome: &'a Genome,
    w: &'a mut dyn RecordWriter,
    chrom: Option<ChromId>,
    len: u64,
    /// Depth changes not yet swept: position -> net change.
    events: BTreeMap<u64, i64>,
    depth: i64,
    /// Start of the stretch at the current depth.
    at: u64,
    /// The run being built: runs of equal depth are joined.
    run: Option<(u64, u64, u64)>,
    /// Bases per depth in this chromosome / the whole genome (`Hist`).
    bins: BTreeMap<u64, u64>,
    total_bins: BTreeMap<u64, u64>,
    total_len: u64,
    next_genome: usize,
    row: Vec<u8>,
    peak_events: usize,
    pieces: Vec<Interval>,
}

impl Cov<'_> {
    fn name(&self) -> &[u8] {
        self.genome
            .name(self.chrom.expect("open chromosome").0 as usize)
    }

    fn emit(&mut self) -> Result<(), IoError> {
        let r = self
            .w
            .write_record(&self.row)
            .map_err(|source| IoError::Write {
                path: "output".into(),
                source,
            });
        self.row.clear();
        r
    }

    fn depth_text(&mut self, d: u64) {
        match self.o.scale {
            None => write!(self.row, "{d}").expect("vec write"),
            Some(s) => write!(self.row, "{}", format_g6(d as f64 * s)).expect("vec write"),
        }
    }

    fn write_run(&mut self, s: u64, e: u64, d: u64) -> Result<(), IoError> {
        match self.o.mode {
            GenomecovMode::Hist => {
                let d = self.o.max_depth.map_or(d, |m| d.min(m));
                *self.bins.entry(d).or_default() += e - s;
            }
            GenomecovMode::BedGraph | GenomecovMode::BedGraphAll => {
                if d == 0 && self.o.mode == GenomecovMode::BedGraph {
                    return Ok(());
                }
                self.row.clear();
                let name = self.name().to_vec();
                self.row.extend_from_slice(&name);
                write!(self.row, "\t{s}\t{e}\t").expect("vec write");
                self.depth_text(d);
                self.emit()?;
            }
            GenomecovMode::PerBase => {
                let name = self.name().to_vec();
                for pos in s + 1..=e {
                    self.row.clear();
                    self.row.extend_from_slice(&name);
                    write!(self.row, "\t{pos}\t").expect("vec write");
                    self.depth_text(d);
                    self.emit()?;
                }
            }
        }
        Ok(())
    }

    fn flush_run(&mut self) -> Result<(), IoError> {
        if let Some((s, e, d)) = self.run.take() {
            self.write_run(s, e, d)?;
        }
        Ok(())
    }

    fn add_run(&mut self, s: u64, e: u64, d: u64) -> Result<(), IoError> {
        if e <= s {
            return Ok(());
        }
        if let Some(r) = self.run.as_mut()
            && r.1 == s
            && r.2 == d
        {
            r.1 = e;
            return Ok(());
        }
        self.flush_run()?;
        self.run = Some((s, e, d));
        Ok(())
    }

    /// Sweep every event before `limit`; later features cannot change them.
    fn sweep_before(&mut self, limit: u64) -> Result<(), IoError> {
        while let Some(entry) = self.events.first_entry() {
            let pos = *entry.key();
            if pos >= limit {
                break;
            }
            let delta = entry.remove();
            self.add_run(self.at, pos, self.depth as u64)?;
            self.depth += delta;
            self.at = pos;
        }
        Ok(())
    }

    fn begin(&mut self, i: usize) {
        self.chrom = Some(ChromId(i as u32));
        self.len = self.genome.length(i).expect("lengths checked");
        self.depth = 0;
        self.at = 0;
        self.bins.clear();
    }

    fn finish_chrom(&mut self) -> Result<(), IoError> {
        if self.chrom.is_none() {
            return Ok(());
        }
        self.sweep_before(u64::MAX)?;
        self.add_run(self.at, self.len, self.depth as u64)?;
        self.flush_run()?;
        if self.o.mode == GenomecovMode::Hist && self.len > 0 {
            let name = self.name().to_vec();
            let bins = std::mem::take(&mut self.bins);
            for (&d, &b) in &bins {
                self.row.clear();
                self.row.extend_from_slice(&name);
                write!(
                    self.row,
                    "\t{d}\t{b}\t{}\t{}",
                    self.len,
                    format_g6(b as f64 / self.len as f64)
                )
                .expect("vec write");
                self.emit()?;
                *self.total_bins.entry(d).or_default() += b;
            }
            self.total_len += self.len;
        }
        self.chrom = None;
        Ok(())
    }

    /// Chromosomes of the genome before `upto` that have no features.
    fn untouched_before(&mut self, upto: usize) -> Result<(), IoError> {
        while self.next_genome < upto {
            let i = self.next_genome;
            self.next_genome += 1;
            self.begin(i);
            self.finish_chrom()?;
        }
        Ok(())
    }

    fn pieces_of(&mut self, m: &RecordMeta, raw: &[u8], path: &str) -> Result<(), IoError> {
        self.pieces.clear();
        let Interval { start, end } = m.span;
        if self.o.split {
            match bed12_blocks(raw) {
                Ok(Some(b)) => self
                    .pieces
                    .extend(b.into_iter().filter(|b| b.start < b.end)),
                Ok(None) => self.pieces.push(m.span),
                Err(msg) => {
                    return Err(IoError::Parse {
                        path: path.to_string(),
                        loc: Loc::Line(m.loc),
                        msg: format!("invalid BED12 blocks: {msg}"),
                    });
                }
            }
            return Ok(());
        }
        // An unknown strand counts as plus, as in slop.
        let minus = m.strand == Strand::Minus;
        match self.o.piece {
            Piece::Whole => self.pieces.push(m.span),
            Piece::FivePrime | Piece::ThreePrime => {
                let first = (self.o.piece == Piece::FivePrime) != minus;
                let p = if first { start } else { end - 1 };
                self.pieces.push(Interval {
                    start: p,
                    end: p + 1,
                });
            }
        }
        Ok(())
    }

    fn record(&mut self, m: &RecordMeta, raw: &[u8], path: &str, bed: bool) -> Result<(), IoError> {
        if m.span.start == m.span.end {
            return Ok(());
        }
        if let Some(want) = self.o.strand {
            if bed && m.strand == Strand::Unknown && raw.split(|&b| b == b'\t').nth(5).is_none() {
                return Err(IoError::Parse {
                    path: path.to_string(),
                    loc: Loc::Line(m.loc),
                    msg: "--strand needs the strand in column 6 but this record has fewer \
than 6 columns; add the strand column or drop --strand"
                        .into(),
                });
            }
            if m.strand != want {
                return Ok(());
            }
        }
        if self.chrom != Some(m.chrom) {
            self.finish_chrom()?;
            let i = m.chrom.0 as usize;
            // A closed genome table: the id is the genome index.
            self.untouched_before(i)?;
            self.next_genome = i + 1;
            self.begin(i);
        }
        if m.span.end > self.len {
            let name = String::from_utf8_lossy(self.name()).into_owned();
            return Err(IoError::Parse {
                path: path.to_string(),
                loc: Loc::Line(m.loc),
                msg: format!(
                    "{name}:{}-{} extends past the end of {name} (length {}); \
check that -g matches the genome build of the input",
                    m.span.start, m.span.end, self.len
                ),
            });
        }
        self.sweep_before(m.key)?;
        self.pieces_of(m, raw, path)?;
        for p in std::mem::take(&mut self.pieces) {
            *self.events.entry(p.start).or_default() += 1;
            *self.events.entry(p.end).or_default() -= 1;
            self.pieces.push(p);
        }
        self.peak_events = self.peak_events.max(self.events.len());
        Ok(())
    }

    fn finish(&mut self) -> Result<(), IoError> {
        self.finish_chrom()?;
        self.untouched_before(self.genome.len())?;
        if self.o.mode == GenomecovMode::Hist && self.total_len > 0 {
            let bins = std::mem::take(&mut self.total_bins);
            for (&d, &b) in &bins {
                self.row.clear();
                write!(
                    self.row,
                    "genome\t{d}\t{b}\t{}\t{}",
                    self.total_len,
                    format_g6(b as f64 / self.total_len as f64)
                )
                .expect("vec write");
                self.emit()?;
            }
        }
        Ok(())
    }
}

/// Write the coverage of `genome` by `input` (sorted in `chroms` order, a
/// closed table built from `genome`, so chromosome ids are genome indexes).
/// Rows per `o.mode`; see the module documentation.
pub fn genomecov(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    genome: &Genome,
    o: &GenomecovOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    genomecov_counting(input, chroms, genome, o, w)
        .map(|_| ())
        .map_err(|e| e.with_sorted_input_loc())
}

/// [`genomecov`], returning the most pending depth-change events ever held
/// (for the memory-flatness test).
#[doc(hidden)]
pub fn genomecov_counting(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    genome: &Genome,
    o: &GenomecovOptions,
    w: &mut dyn RecordWriter,
) -> Result<usize, IoError> {
    validate(o)?;
    if let Some(i) = (0..genome.len()).find(|&i| genome.length(i).is_none()) {
        return Err(UsageError(format!(
            "-g: the genome has no length for {}; genomecov needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
            String::from_utf8_lossy(genome.name(i))
        ))
        .into());
    }
    let path = input.name().to_string();
    let bed = matches!(input.format(), Format::Bed | Format::BedGraph | Format::Wig);
    let mut skipped = Skipped::default();
    let input = skipped.watch(input);
    let mut s = Cov {
        o,
        genome,
        w,
        chrom: None,
        len: 0,
        events: BTreeMap::new(),
        depth: 0,
        at: 0,
        run: None,
        bins: BTreeMap::new(),
        total_bins: BTreeMap::new(),
        total_len: 0,
        next_genome: 0,
        row: Vec::new(),
        peak_events: 0,
        pieces: Vec::new(),
    };
    kmerge_with(vec![input], chroms, &mut |_, meta, raw, _| {
        s.record(meta, raw, &path, bed)
    })?;
    s.finish()?;
    log_skipped(&skipped.take());
    Ok(s.peak_events)
}

#[cfg(test)]
mod tests {
    use super::format_g6;

    #[test]
    fn g6_matches_cpp_default_stream() {
        for (x, want) in [
            (0.0, "0"),
            (1.0, "1"),
            (0.6, "0.6"),
            (0.676056338, "0.676056"),
            (0.0704225352, "0.0704225"),
            (1.5, "1.5"),
            (7.5, "7.5"),
            (1e-7, "1e-07"),
            (3e-7, "3e-07"),
            (1e8, "1e+08"),
            (123456.0, "123456"),
            (1234567.0, "1.23457e+06"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (-3.0, "-3"),
            (0.3333333333, "0.333333"),
            (999999.5, "1e+06"),
        ] {
            assert_eq!(format_g6(x), want, "{x}");
        }
    }
}
