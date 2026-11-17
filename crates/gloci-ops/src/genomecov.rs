//! `genomecov`: the depth of coverage of a genome by the features of one
//! sorted BED-family or alignment input (replaces `bedtools genomecov -i` and
//! `-ibam`).
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
//! (`-g`, or an alignment header's lengths), in genome order, chromosomes
//! without features included. Points cover no bases (D68) and are ignored.
//!
//! An alignment counts its span (D101), or with `split`/`split_deletions` its
//! blocks (D105, the shared `AlnBlocks`). An indexed BAM at more than one
//! thread runs per chromosome (D117): chunks of whole genome chromosomes are
//! swept in parallel and written in genome order, and the `genome` rows of
//! the histogram are summed from the chunks, so the output is the same at
//! every thread count.

use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable, Genome, Interval, Strand};
use gloci_engine::{all_indexed, kmerge_with, run_chunks};
use gloci_io::bed::bed12_blocks;
use gloci_io::{
    AlnBlocks, Diagnostics, Format, IndexedFile, IoError, Loc, ReaderOptions, RecordMeta,
    RecordReader, RecordWriter, UsageError,
};
use std::collections::BTreeMap;
use std::io::Write;
use std::num::NonZeroUsize;
use std::ops::Range;
use std::path::PathBuf;

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
    /// Count each BED12 block, or each alignment block between `N` operations, instead of the
    /// whole feature.
    pub split: bool,
    /// Alignments only: also cut the blocks at deletions (D105).
    pub split_deletions: bool,
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
    if (o.split || o.split_deletions) && o.piece != Piece::Whole {
        return u(
            "--split and --split-deletions cannot be combined with --5p or --3p; use one of them",
        );
    }
    if !matches!(o.strand, None | Some(Strand::Plus | Strand::Minus)) {
        return u("--strand must be + or -");
    }
    Ok(())
}

/// The `--hist` fraction as bedtools computes it: both counts are rounded to
/// `float`, divided as `float`, and the result is printed with six
/// significant digits. (A division in double differs in the last digit for
/// some values: 162111521 / 207915645 is 0.779698 as float, 0.779699 as
/// double.)
fn format_fraction(bases: u64, len: u64) -> String {
    format_g6(f64::from(bases as f32 / len as f32))
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

impl Cov<'_> {
    /// Where record number `n` of the input is: a line of text, or a record of a binary file.
    fn loc(&self, n: u64) -> Loc {
        if self.format.is_binary() {
            Loc::Record(n)
        } else {
            Loc::Line(n)
        }
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
    /// The input's format, for the place of an error.
    format: Format,
    /// The chromosome of the last record read, and how many records of it were read: an
    /// alignment error names its place so (the same at every thread count).
    seen_on: Option<ChromId>,
    seen: u64,
    /// Where an alignment's blocks come from, when it has any to count.
    blocks: Option<AlnBlocks>,
}

impl<'a> Cov<'a> {
    fn name(&self) -> &'a [u8] {
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
                let name = self.name();
                self.row.extend_from_slice(name);
                write!(self.row, "\t{s}\t{e}\t").expect("vec write");
                self.depth_text(d);
                self.emit()?;
            }
            GenomecovMode::PerBase => {
                let name = self.name();
                for pos in s + 1..=e {
                    self.row.clear();
                    self.row.extend_from_slice(name);
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
            let name = self.name();
            let bins = std::mem::take(&mut self.bins);
            for (&d, &b) in &bins {
                self.row.clear();
                self.row.extend_from_slice(name);
                write!(
                    self.row,
                    "\t{d}\t{b}\t{}\t{}",
                    self.len,
                    format_fraction(b, self.len)
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
        if let Some(b) = &self.blocks
            && (self.o.split || self.o.split_deletions)
        {
            self.pieces.extend(b.get(raw, m.loc, m.span)?);
            return Ok(());
        }
        if self.o.split {
            match bed12_blocks(raw) {
                Ok(Some(b)) => self
                    .pieces
                    .extend(b.into_iter().filter(|b| b.start < b.end)),
                Ok(None) => self.pieces.push(m.span),
                Err(msg) => {
                    return Err(IoError::Parse {
                        path: path.to_string(),
                        loc: self.loc(m.loc),
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
        if self.seen_on != Some(m.chrom) {
            self.seen_on = Some(m.chrom);
            self.seen = 0;
        }
        self.seen += 1;
        if m.span.start == m.span.end {
            return Ok(());
        }
        if let Some(want) = self.o.strand {
            if bed && m.strand == Strand::Unknown && raw.split(|&b| b == b'\t').nth(5).is_none() {
                return Err(IoError::Parse {
                    path: path.to_string(),
                    loc: self.loc(m.loc),
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
            let msg = format!(
                "{name}:{}-{} extends past the end of {name} (length {}); \
check that -g matches the genome build of the input",
                m.span.start, m.span.end, self.len
            );
            if self.format.is_alignment() {
                // The reader's own record number differs between a streamed and a
                // per-chromosome run, so the place is counted here, among the reads that
                // reach the sweep (D101's place, minus the reads left out), the same at
                // every thread count.
                return Err(IoError::Detect {
                    path: path.to_string(),
                    msg: format!(
                        "{path} record {} of {name} (counting only the reads genomecov counts): {msg}",
                        self.seen
                    ),
                });
            }
            return Err(IoError::Parse {
                path: path.to_string(),
                loc: self.loc(m.loc),
                msg,
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

    /// Close the open chromosome and cover the genome chromosomes up to `end` that had no
    /// features.
    fn finish_to(&mut self, end: usize) -> Result<(), IoError> {
        self.finish_chrom()?;
        self.untouched_before(end)
    }
}

/// The `genome` rows of the histogram, from the bases per depth over all chromosomes.
fn write_genome_rows(
    w: &mut dyn RecordWriter,
    bins: &BTreeMap<u64, u64>,
    total_len: u64,
) -> Result<(), IoError> {
    if total_len == 0 {
        return Ok(());
    }
    let mut row = Vec::new();
    for (&d, &b) in bins {
        row.clear();
        write!(
            row,
            "genome\t{d}\t{b}\t{total_len}\t{}",
            format_fraction(b, total_len)
        )
        .expect("vec write");
        w.write_record(&row).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })?;
    }
    Ok(())
}

/// What sweeping a run of genome chromosomes leaves for the caller to combine.
struct Swept {
    total_bins: BTreeMap<u64, u64>,
    total_len: u64,
    peak_events: usize,
    skipped: Vec<(String, Diagnostics)>,
}

/// Sweep `input` (sorted in `chroms` order, a closed table built from `genome`, so chromosome
/// ids are genome indexes) over the genome chromosomes `range`, writing their rows to `w`.
/// The `genome` histogram rows are not written: the totals come back.
fn sweep(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    genome: &Genome,
    o: &GenomecovOptions,
    range: Range<usize>,
    w: &mut dyn RecordWriter,
) -> Result<Swept, Box<(IoError, Diagnostics)>> {
    let path = input.name().to_string();
    let format = input.format();
    let bed = matches!(
        format,
        Format::Bed | Format::BedGraph | Format::Wig | Format::BigWig | Format::BigBed
    );
    let blocks = if o.split || o.split_deletions {
        AlnBlocks::for_reader(input.as_ref(), o.split_deletions)
    } else {
        None
    };
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
        next_genome: range.start,
        row: Vec::new(),
        peak_events: 0,
        pieces: Vec::new(),
        format,
        seen_on: None,
        seen: 0,
        blocks,
    };
    let done = kmerge_with(vec![input], chroms, &mut |_, meta, raw, _| {
        s.record(meta, raw, &path, bed)
    })
    .and_then(|_| s.finish_to(range.end));
    let (total_bins, total_len, peak_events) = (s.total_bins, s.total_len, s.peak_events);
    // The reader is gone: what it noted (even on the way to an error) is kept.
    let diag = skipped.take().into_iter().next().map(|(_, d)| d);
    match done {
        Ok(()) => Ok(Swept {
            total_bins,
            total_len,
            peak_events,
            skipped: skipped_named(&path, diag),
        }),
        Err(e) => Err(Box::new((e, diag.unwrap_or_default()))),
    }
}

fn skipped_named(path: &str, d: Option<Diagnostics>) -> Vec<(String, Diagnostics)> {
    d.map(|d| (path.to_string(), d)).into_iter().collect()
}

fn check_genome(
    genome: &Genome,
    o: &GenomecovOptions,
    input: &dyn RecordReader,
) -> Result<(), IoError> {
    validate(o)?;
    if o.split_deletions && !input.format().is_alignment() {
        return Err(UsageError(
            "--split-deletions cuts alignments at deletions, and the input is not BAM, CRAM or \
SAM; drop it, or use --split for BED12 blocks"
                .into(),
        )
        .into());
    }
    if let Some(i) = (0..genome.len()).find(|&i| genome.length(i).is_none()) {
        return Err(UsageError(format!(
            "-g: the genome has no length for {}; genomecov needs chromosome lengths \
(use a built-in genome, a .fai or .dict file, or a BAM file)",
            String::from_utf8_lossy(genome.name(i))
        ))
        .into());
    }
    Ok(())
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
    check_genome(genome, o, input.as_ref())?;
    let alignment = input.format().is_alignment();
    let swept = sweep(input, chroms, genome, o, 0..genome.len(), w).map_err(|f| f.0)?;
    if o.mode == GenomecovMode::Hist {
        write_genome_rows(w, &swept.total_bins, swept.total_len)?;
    }
    if alignment && let Some((name, d)) = swept.skipped.first() {
        d.log_alignment_counts(name);
    }
    log_skipped(&swept.skipped);
    Ok(swept.peak_events)
}

/// How an indexed BAM is read per chromosome: its path and reader options, and the run's
/// threads and spill directory.
pub struct ChunkRun {
    pub path: PathBuf,
    pub ro: ReaderOptions,
    pub threads: NonZeroUsize,
    pub tmp_dir: PathBuf,
}

/// Chunks per thread, as the relations.
const CHUNKS_PER_THREAD: usize = 4;

/// One chunk of a per-chromosome run: the genome chromosomes `range`, of which `names` hold
/// reads.
#[derive(Clone)]
struct Job {
    range: Range<usize>,
    names: Vec<Vec<u8>>,
}

/// The BAM or CRAM `reader` reads, ready to be read a chromosome at a time, or `None` when it must
/// stream: not an indexed BAM or CRAM, or one streaming reports errors for (a chromosome the genome
/// lacks, or an order the genome rejects).
fn open_chunked(
    reader: &dyn RecordReader,
    chroms: &ChromTable,
    run: &ChunkRun,
) -> Result<Option<IndexedFile>, IoError> {
    let Some(path) = run.path.to_str() else {
        return Ok(None);
    };
    if !matches!(reader.format(), Format::Bam | Format::Cram) || !all_indexed(&[path]) {
        return Ok(None);
    }
    let mut file = IndexedFile::open(&run.path, reader.format(), reader.header().clone(), &run.ro)?;
    // An alignment reader's one-time warnings are printed once, by this run.
    file.defer_warnings();
    let mut probe = chroms.clone();
    if file
        .non_empty_reference_names()
        .iter()
        .any(|n| probe.intern(n).is_err())
    {
        log::debug!("the input has a chromosome the genome lacks; streaming");
        return Ok(None);
    }
    let known: Vec<Vec<u8>> = file
        .file_order_reference_names()
        .into_iter()
        .filter(|n| chroms.cmp_names(n, n).is_some())
        .collect();
    if known
        .windows(2)
        .any(|w| chroms.cmp_names(&w[0], &w[1]) == Some(std::cmp::Ordering::Greater))
    {
        log::debug!("the input is out of genome order; streaming");
        return Ok(None);
    }
    file.info().warn_if_stale(file.name());
    Ok(Some(file))
}

/// [`genomecov`] for an input that may be an indexed BAM: with `chunked` and a BAM that has an
/// index, chunks of whole genome chromosomes are swept in parallel and written in genome order
/// (D117), with the same output as streaming; anything else streams. `input` is the opened
/// reader (its header and format are used).
pub fn genomecov_files(
    input: Box<dyn RecordReader>,
    chunked: Option<&ChunkRun>,
    chroms: &mut ChromTable,
    genome: &Genome,
    o: &GenomecovOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    let file = match chunked {
        Some(run) if run.threads.get() > 1 => {
            check_genome(genome, o, input.as_ref())?;
            open_chunked(input.as_ref(), chroms, run)?.map(|f| (f, run))
        }
        _ => None,
    };
    let Some((file, run)) = file else {
        return genomecov(input, chroms, genome, o, w);
    };
    let names: Vec<Vec<u8>> = (0..genome.len()).map(|i| genome.name(i).to_vec()).collect();
    let in_file: Vec<bool> = names.iter().map(|n| file.weight(n) > 0).collect();
    let weights: Vec<u64> = names.iter().map(|n| file.weight(n) + 1).collect();
    let mut start = 0;
    let jobs: Vec<Job> = crate::relate::group_chromosomes(
        names.clone(),
        &weights,
        CHUNKS_PER_THREAD * run.threads.get(),
    )
    .into_iter()
    .map(|g| {
        let range = start..start + g.len();
        start = range.end;
        let names = range
            .clone()
            .filter(|&i| in_file[i])
            .map(|i| names[i].clone())
            .collect();
        Job { range, names }
    })
    .collect();
    let table = chroms.clone();
    // What the chunks that failed had noted: a streamed run has printed the one-time warnings
    // met before its error, so the run prints these.
    let failed: std::sync::Mutex<Vec<Diagnostics>> = std::sync::Mutex::new(Vec::new());
    let step = |job: &Job, sink: &mut gloci_engine::ChunkSink| -> Result<Swept, IoError> {
        let reader = file.chromosomes_reader(&job.names)?;
        let mut chroms = table.clone();
        sweep(reader, &mut chroms, genome, o, job.range.clone(), sink).map_err(|f| {
            let (e, d) = *f;
            failed.lock().unwrap_or_else(|p| p.into_inner()).push(d);
            e
        })
    };
    log::debug!("genomecov running per chromosome: {} chunk(s)", jobs.len());
    let swept = match run_chunks(jobs, run.threads, &run.tmp_dir, &step, w) {
        Ok(s) => s,
        Err(e) => {
            for d in failed.into_inner().unwrap_or_else(|p| p.into_inner()) {
                d.log_deferred_warnings();
            }
            return Err(e.with_sorted_input_loc());
        }
    };

    let mut total_bins: BTreeMap<u64, u64> = BTreeMap::new();
    let mut total_len = 0;
    let mut skipped: Vec<(String, Diagnostics)> = Vec::new();
    for s in swept {
        for (d, b) in s.total_bins {
            *total_bins.entry(d).or_default() += b;
        }
        total_len += s.total_len;
        if skipped.is_empty() {
            skipped = s.skipped;
        } else {
            for ((_, d), (_, od)) in skipped.iter_mut().zip(s.skipped) {
                d.absorb(od);
            }
        }
    }
    if o.mode == GenomecovMode::Hist {
        write_genome_rows(w, &total_bins, total_len)?;
    }
    // What streaming reports once for the whole input, from the chunks together: the reads
    // with no place at the end of the file are in no chromosome's chunk, the index counts them.
    if let Some((name, d)) = skipped.first_mut() {
        d.unmapped += file.unplaced_reads();
        d.log_deferred_warnings();
        d.log_alignment_counts(name);
    }
    log_skipped(&skipped);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{format_fraction, format_g6};

    #[test]
    fn fractions_divide_in_single_precision_like_bedtools() {
        // Printed by bedtools genomecov 2.31.1 (golden genomecov-hist-float-fraction-vs-bedtools).
        assert_eq!(format_fraction(162_111_521, 207_915_645), "0.779698");
        assert_eq!(format_fraction(143_954_154, 166_830_933), "0.862875");
        assert_eq!(format_fraction(234_958_625, 250_000_000), "0.939834");
        assert_eq!(format_fraction(3, 5), "0.6");
        assert_eq!(format_fraction(7, 7), "1");
        assert_eq!(format_fraction(0, 9), "0");
    }

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
