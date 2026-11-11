//! `shuffle --format bedpe` (D153).
//!
//! A pair whose ends are on one chromosome is placed as its outer span, so the
//! ends keep their lengths and the gap between them; a pair on two chromosomes
//! places each end independently; an absent end (`. -1 -1`) stays absent.
//! With `--exclude` every end must avoid the excluded regions (the gap between
//! the ends of a one-chromosome pair may cross them, as in bedtools); a place
//! that does not is redrawn up to `--max-tries` times. `--include`,
//! `--same-chrom`, `--chrom-first` and `--no-overlapping` are usage errors
//! (bedtools ignores them).

use super::allowed::Regions;
use super::{
    Allowed, ShuffleOptions, ShuffleSummary, load_exclude, prepare_genome, warn_unplaceable,
};
use crate::rng::Stream;
use gloci_core::{ChromId, ChromTable, Genome};
use gloci_io::{BedpeBatch, BedpeReader, IoError, Loc, UsageError};
use std::io::Write;

/// Where one end went: chromosome, start, end.
type Moved = Option<(ChromId, u64, u64)>;

fn refuse_flags(o: &ShuffleOptions) -> Result<(), UsageError> {
    let flag = if o.include.is_some() {
        "--include"
    } else if o.same_chrom {
        "--same-chrom"
    } else if o.chrom_first {
        "--chrom-first"
    } else if o.no_overlapping {
        "--no-overlapping"
    } else {
        return Ok(());
    };
    Err(UsageError(format!(
        "{flag} cannot be used with --format bedpe (bedtools ignores it silently, so gloci \
refuses it); BEDPE shuffle takes only -g, --exclude, --seed, --max-tries and --skip-unplaceable"
    )))
}

/// Place `len` bases anywhere, redrawing while the result overlaps an excluded region.
fn place_end(
    allowed: &Allowed,
    excl: &Regions,
    o: &ShuffleOptions,
    len: u64,
    rng: &mut Stream,
) -> Result<(ChromId, u64), bool> {
    let tries = if excl.is_empty() { 1 } else { o.max_tries };
    for _ in 0..tries {
        let (c, s) = allowed.place(len, None, rng).ok_or(false)?;
        if !excl.overlaps(c, s, s + len) {
            return Ok((c, s));
        }
    }
    Err(true)
}

/// The new ends of one pair, or why it has no place (true: it kept landing on
/// an excluded region).
fn place_pair(
    allowed: &Allowed,
    excl: &Regions,
    o: &ShuffleOptions,
    ends: &[Option<gloci_io::BedpeEnd>; 2],
    rng: &mut Stream,
) -> Result<[Moved; 2], bool> {
    match ends {
        [Some(a), Some(b)] if a.chrom == b.chrom => {
            let lo = a.span.start.min(b.span.start);
            let hi = a.span.end.max(b.span.end);
            let tries = if excl.is_empty() { 1 } else { o.max_tries };
            for _ in 0..tries {
                let (c, s) = allowed.place(hi - lo, None, rng).ok_or(false)?;
                let shift = |e: &gloci_io::BedpeEnd| {
                    let start = s + (e.span.start - lo);
                    (c, start, start + e.span.len())
                };
                let (na, nb) = (shift(a), shift(b));
                if !excl.overlaps(c, na.1, na.2) && !excl.overlaps(c, nb.1, nb.2) {
                    return Ok([Some(na), Some(nb)]);
                }
            }
            Err(true)
        }
        _ => {
            let mut out = [None, None];
            for (slot, e) in out.iter_mut().zip(ends) {
                if let Some(e) = e {
                    let len = e.span.len();
                    let (c, s) = place_end(allowed, excl, o, len, rng)?;
                    *slot = Some((c, s, s + len));
                }
            }
            Ok(out)
        }
    }
}

/// `raw` with the chromosome and coordinates of the moved ends replaced.
fn write_moved(raw: &[u8], chroms: &ChromTable, moved: &[Moved; 2], out: &mut Vec<u8>) {
    for (i, f) in raw.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        match (i, moved) {
            (0..=2, [Some(m), _]) | (3..=5, [_, Some(m)]) => match i % 3 {
                0 => out.extend_from_slice(chroms.name(m.0)),
                1 => out.extend_from_slice(m.1.to_string().as_bytes()),
                _ => out.extend_from_slice(m.2.to_string().as_bytes()),
            },
            _ => out.extend_from_slice(f),
        }
    }
    out.push(b'\n');
}

/// Write every pair of `reader` moved to a random place, in input order (the
/// header lines are the caller's to write).
///
/// Pair *j* (from 0, counting skipped ones) uses stream *j* of `o.seed`.
pub fn shuffle_bedpe(
    mut reader: BedpeReader,
    chroms: &mut ChromTable,
    genome: &Genome,
    o: &ShuffleOptions,
    w: &mut dyn Write,
) -> Result<ShuffleSummary, IoError> {
    refuse_flags(o)?;
    prepare_genome(genome, chroms)?;
    let threads = std::num::NonZeroUsize::MIN;
    let excl = load_exclude(o.exclude.as_deref(), chroms, threads)?;
    let allowed = Allowed::build(genome, chroms, None, &[]);
    let name = reader.name().to_string();
    let mut summary = ShuffleSummary::default();
    let mut listed: Vec<String> = Vec::new();
    let mut batch = BedpeBatch::default();
    let mut out = Vec::new();
    let mut max = 1;
    let mut j = 0u64;
    loop {
        batch.clear();
        let more = reader.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            summary.records += 1;
            let mut rng = Stream::new(o.seed, j);
            j += 1;
            out.clear();
            match place_pair(&allowed, &excl, o, &m.ends, &mut rng) {
                Ok(moved) => write_moved(batch.raw(i), chroms, &moved, &mut out),
                Err(overlapped) => {
                    let lens: Vec<String> = m
                        .ends
                        .iter()
                        .flatten()
                        .map(|e| e.span.len().to_string())
                        .collect();
                    if !o.skip_unplaceable {
                        let why = if overlapped {
                            format!(
                                "no place for this pair in {} draws that keep its ends off the \
--exclude regions; raise --max-tries, shrink --exclude",
                                o.max_tries
                            )
                        } else {
                            "no chromosome is long enough for this pair; use a genome with a \
longer chromosome"
                                .to_string()
                        };
                        return Err(IoError::Parse {
                            path: name,
                            loc: Loc::Line(m.loc),
                            msg: format!(
                                "{why}, or pass --skip-unplaceable to drop such pairs (end \
lengths {})",
                                lens.join(", ")
                            ),
                        });
                    }
                    summary.unplaceable += 1;
                    if listed.len() < super::LISTED {
                        listed.push(format!("{name} line {} lengths {}", m.loc, lens.join(", ")));
                    }
                    continue;
                }
            }
            w.write_all(&out).map_err(|source| IoError::Write {
                path: "output".into(),
                source,
            })?;
        }
        let per = batch.buf.len() / batch.len().max(1) + 1;
        max = (4_194_304 / per).clamp(1, 16384);
        if !more {
            break;
        }
    }
    if let Some(s) = reader.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    warn_unplaceable(&summary, &listed);
    Ok(summary)
}
