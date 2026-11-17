//! `shuffle`: every record moved to a uniformly random place (D152, D153),
//! the equivalent of `bedtools shuffle`.
//!
//! A record of length `L` goes to any `(chromosome, start)` that keeps the
//! whole record inside the allowed set (`--include` regions, else the genome,
//! less `--exclude`), each equally likely ([`Allowed`]). Data record *j* draws
//! from stream *j* of the seed, retries included (D144), so the output is the
//! same at any thread count. Records keep their other columns and their
//! strand; a record with no place to go is an error unless
//! `--skip-unplaceable` says to drop it.

pub mod allowed;
mod bedpe;

pub use allowed::Allowed;
pub use bedpe::{refuse_flags as refuse_bedpe_flags, shuffle_bedpe};

use crate::rng::Stream;
use gloci_core::{ChromId, ChromTable, Genome, Interval};
use gloci_io::{
    Format, IoError, Loc, ReaderOptions, RecordBatch, RecordReader, RecordWriter, UsageError,
    next_batch_size, open_reader,
};
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::num::NonZeroUsize;

#[derive(Debug, Clone)]
pub struct ShuffleOptions {
    pub seed: u64,
    /// A file of regions the records must stay inside (default: the genome).
    pub include: Option<String>,
    /// A file of regions the records must not overlap.
    pub exclude: Option<String>,
    /// Keep each record on its chromosome.
    pub same_chrom: bool,
    /// Pick the chromosome uniformly among those with room, then the place.
    pub chrom_first: bool,
    /// No two output records overlap.
    pub no_overlapping: bool,
    /// Draws per record under `--no-overlapping` (and for BEDPE `--exclude`).
    pub max_tries: u32,
    /// Drop records with no place instead of failing.
    pub skip_unplaceable: bool,
    /// `--chrom-map`: renames the chromosomes of `include` and `exclude` as it
    /// does those of the input.
    pub chrom_map: Option<std::sync::Arc<gloci_io::ChromMap>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ShuffleSummary {
    /// Records read.
    pub records: u64,
    /// Records that had no place (skipped).
    pub unplaceable: u64,
}

/// Unplaceable records listed in the warning.
const LISTED: usize = 10;

/// Every chromosome of the genome needs a length, and gets an id.
pub(crate) fn prepare_genome(genome: &Genome, chroms: &mut ChromTable) -> Result<(), IoError> {
    for g in 0..genome.len() {
        if genome.length(g).is_none() {
            return Err(UsageError(format!(
                "-g: the genome has no length for {}; shuffle needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
                String::from_utf8_lossy(genome.name(g))
            ))
            .into());
        }
        chroms.intern(genome.name(g)).map_err(|e| {
            IoError::Usage(UsageError(format!(
                "-g: cannot use chromosome {} here: {e}",
                String::from_utf8_lossy(genome.name(g))
            )))
        })?;
    }
    Ok(())
}

/// Read the regions of `path` (BED, GFF or VCF; any order).
pub(crate) fn load_regions(
    path: &str,
    chroms: &mut ChromTable,
    threads: NonZeroUsize,
    chrom_map: Option<std::sync::Arc<gloci_io::ChromMap>>,
) -> Result<Vec<(ChromId, Interval)>, IoError> {
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map,
        threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: gloci_io::InvalidPolicy::Fail,
    };
    let mut reader = open_reader(path, &ro)?;
    let mut out = Vec::new();
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = reader.read_batch(chroms, &mut batch, max)?;
        out.extend(batch.recs.iter().map(|m| (m.chrom, m.span)));
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    Ok(out)
}

/// Warn about regions on chromosomes the genome does not have.
pub(crate) fn warn_unknown(
    regions: &[(ChromId, Interval)],
    known: impl Fn(ChromId) -> bool,
    flag: &str,
    path: &str,
) {
    let n = regions.iter().filter(|(c, _)| !known(*c)).count();
    if n > 0 {
        log::warn!(
            "{n} of the {} regions in {path} ({flag}) are on chromosomes that are not in -g \
and were ignored",
            regions.len()
        );
    }
}

pub(crate) fn warn_unplaceable(summary: &ShuffleSummary, listed: &[String]) {
    let n = summary.unplaceable;
    if n == 0 {
        return;
    }
    let (noun, verb) = if n == 1 {
        ("record", "was")
    } else {
        ("records", "were")
    };
    let more = if n as usize > listed.len() {
        format!(", and {} more", n as usize - listed.len())
    } else {
        String::new()
    };
    log::warn!(
        "{n} {noun} could not be placed and {verb} skipped ({}{more})",
        listed.join("; ")
    );
}

/// The refusal for an input that is not intervals.
fn refuse_format(name: &str, f: Format) -> Result<(), UsageError> {
    match f {
        Format::Bed
        | Format::BedGraph
        | Format::Gff3
        | Format::Gtf
        | Format::Wig
        | Format::BigWig
        | Format::BigBed => Ok(()),
        Format::Vcf | Format::Bcf => Err(UsageError(format!(
            "{name} is {}; shuffle moves intervals, and a variant cannot be moved without \
changing its REF; convert it to intervals first: gloci convert {name} -o x.bed",
            f.name()
        ))),
        f if f.is_alignment() => Err(UsageError(format!(
            "{name} is {} alignments; shuffle reads intervals; convert them first: \
gloci convert {name} -o x.bed",
            f.name()
        ))),
        f => Err(UsageError(format!(
            "{name} is {}; shuffle reads BED, bedGraph and GFF/GTF (and BEDPE with --format bedpe)",
            f.name()
        ))),
    }
}

/// `raw` moved to `place`: chromosome name and coordinates replaced.
fn moved(
    format: Format,
    raw: &[u8],
    name: &[u8],
    place: (ChromId, u64),
    len: u64,
    out: &mut Vec<u8>,
) {
    let iv = Interval {
        start: place.1,
        end: place.1 + len,
    };
    let mut tmp = Vec::with_capacity(raw.len() + 8);
    crate::rewrite_coords(format, raw, iv, &mut tmp);
    out.clear();
    out.extend_from_slice(name);
    if let Some(i) = tmp.iter().position(|&b| b == b'\t') {
        out.extend_from_slice(&tmp[i..]);
    }
}

/// Places already taken, for `--no-overlapping`.
struct Taken {
    by_chrom: HashMap<ChromId, BTreeSet<(u64, u64)>>,
}

impl Taken {
    /// Whether `[start, end)` overlaps a taken interval. Touching intervals
    /// and points on an edge do not overlap.
    fn overlaps(&self, c: ChromId, start: u64, end: u64) -> bool {
        // Taken intervals are disjoint, so those starting before `end` end in
        // the same order and the last one reaches furthest.
        self.by_chrom
            .get(&c)
            .and_then(|s| s.range(..(end, 0)).next_back())
            .is_some_and(|&(_, e)| e > start)
    }

    fn insert(&mut self, c: ChromId, start: u64, end: u64) {
        self.by_chrom.entry(c).or_default().insert((start, end));
    }
}

/// One draw for a record of `len` bases.
fn draw(
    allowed: &Allowed,
    o: &ShuffleOptions,
    len: u64,
    chrom: ChromId,
    rng: &mut Stream,
) -> Option<(ChromId, u64)> {
    if o.chrom_first {
        allowed.place_chrom_first(len, rng)
    } else {
        allowed.place(len, o.same_chrom.then_some(chrom), rng)
    }
}

/// Why a record of `len` bases has no place, and what to change.
fn no_room(
    o: &ShuffleOptions,
    chroms: &ChromTable,
    chrom: ChromId,
    len: u64,
    overlapped: bool,
) -> String {
    if overlapped {
        return format!(
            "no place for this record (length {len}) in {} draws that avoids the records \
already placed; raise --max-tries, shrink --exclude or widen --include, or pass \
--skip-unplaceable to drop such records",
            o.max_tries
        );
    }
    let where_ = if o.same_chrom {
        format!(" on {}", String::from_utf8_lossy(chroms.name(chrom)))
    } else {
        String::new()
    };
    format!(
        "no place for this record (length {len}): no allowed region{where_} is {len} bases or \
longer; shrink --exclude, widen --include{}, or pass --skip-unplaceable to drop such records",
        if o.same_chrom {
            ", drop --same-chrom"
        } else {
            ""
        }
    )
}

/// Write every record of `reader` moved to a random place, in input order, and
/// return how many were read and how many had no place.
///
/// Record *j* (from 0, counting skipped ones) uses stream *j* of `o.seed`.
/// The genome must carry chromosome lengths. `threads` only speeds the
/// drawing up; the output never depends on it.
pub fn shuffle(
    mut reader: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    genome: &Genome,
    o: &ShuffleOptions,
    threads: NonZeroUsize,
    w: &mut dyn RecordWriter,
) -> Result<ShuffleSummary, IoError> {
    refuse_format(reader.name(), reader.format())?;
    if o.same_chrom && o.chrom_first {
        return Err(UsageError(
            "--same-chrom and --chrom-first contradict each other: the first keeps each record \
on its chromosome, the second picks a chromosome at random"
                .into(),
        )
        .into());
    }
    prepare_genome(genome, chroms)?;
    let include = match &o.include {
        Some(p) => Some((load_regions(p, chroms, threads, o.chrom_map.clone())?, p)),
        None => None,
    };
    let exclude = match &o.exclude {
        Some(p) => Some((load_regions(p, chroms, threads, o.chrom_map.clone())?, p)),
        None => None,
    };
    let allowed = Allowed::build(
        genome,
        chroms,
        include.as_ref().map(|(r, _)| r.as_slice()),
        exclude.as_ref().map_or(&[], |(r, _)| r.as_slice()),
    );
    for (regions, flag) in [(&include, "--include"), (&exclude, "--exclude")] {
        if let Some((r, p)) = regions {
            warn_unknown(r, |c| allowed.knows(c), flag, p);
        }
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .build()
        .map_err(|e| IoError::Internal(format!("cannot start the thread pool: {e}")))?;
    let write_err = |source| IoError::Write {
        path: "output".into(),
        source,
    };
    let format = reader.format();
    let name = reader.name().to_string();
    let mut summary = ShuffleSummary::default();
    let mut listed: Vec<String> = Vec::new();
    let mut taken = Taken {
        by_chrom: HashMap::new(),
    };
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    let mut max = 1;
    let mut j0 = 0u64;
    loop {
        batch.clear();
        let more = reader.read_batch(chroms, &mut batch, max)?;
        // Each record's place, or why it has none (true: it kept overlapping).
        let places: Vec<Result<(ChromId, u64), bool>> = if o.no_overlapping {
            batch
                .recs
                .iter()
                .enumerate()
                .map(|(i, m)| {
                    let len = m.span.len();
                    let mut rng = Stream::new(o.seed, j0 + i as u64);
                    let mut tries = 0u32;
                    loop {
                        let Some((c, s)) = draw(&allowed, o, len, m.chrom, &mut rng) else {
                            return Err(false);
                        };
                        if !taken.overlaps(c, s, s + len) {
                            taken.insert(c, s, s + len);
                            return Ok((c, s));
                        }
                        tries += 1;
                        if tries >= o.max_tries {
                            return Err(true);
                        }
                    }
                })
                .collect()
        } else {
            let one = |(i, m): (usize, &gloci_io::RecordMeta)| {
                let mut rng = Stream::new(o.seed, j0 + i as u64);
                draw(&allowed, o, m.span.len(), m.chrom, &mut rng).ok_or(false)
            };
            if threads.get() > 1 && batch.len() > 256 {
                pool.install(|| batch.recs.par_iter().enumerate().map(one).collect())
            } else {
                batch.recs.iter().enumerate().map(one).collect()
            }
        };
        for (i, (m, place)) in batch.recs.iter().zip(places).enumerate() {
            summary.records += 1;
            match place {
                Ok(place) => {
                    moved(
                        format,
                        batch.raw(i),
                        chroms.name(place.0),
                        place,
                        m.span.len(),
                        &mut out,
                    );
                    w.write_record(&out).map_err(write_err)?;
                }
                Err(overlapped) => {
                    let len = m.span.len();
                    if !o.skip_unplaceable {
                        return Err(IoError::Parse {
                            path: name,
                            loc: Loc::Line(m.loc),
                            msg: no_room(o, chroms, m.chrom, len, overlapped),
                        });
                    }
                    summary.unplaceable += 1;
                    if listed.len() < LISTED {
                        listed.push(format!("{name} line {} length {len}", m.loc));
                    }
                }
            }
        }
        j0 += batch.len() as u64;
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    crate::skipped::log_reader(reader.as_ref());
    warn_unplaceable(&summary, &listed);
    Ok(summary)
}
