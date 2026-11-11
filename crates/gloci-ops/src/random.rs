//! `random`: intervals placed uniformly at random over a genome (D150),
//! the equivalent of `bedtools random`.
//!
//! Every valid start of every chromosome is equally likely: chromosome `c`
//! has `len_c - l + 1` valid starts (none when it is shorter than `l`), and
//! interval *i* draws one position over all of them with prefix sums and a
//! binary search, then a coin for the strand. Interval *i* uses stream *i*
//! of the seed (D144), so the output is the same at any thread count.

use crate::rng::Stream;
use gloci_core::Genome;
use gloci_io::{IoError, RecordWriter, UsageError};
use rayon::prelude::*;
use std::num::NonZeroUsize;

#[derive(Debug, Clone, Copy)]
pub struct RandomOptions {
    /// Interval length in bases (at least 1).
    pub length: u64,
    /// How many intervals.
    pub count: u64,
    pub seed: u64,
}

/// Intervals per parallel work item, and per written batch.
const CHUNK: u64 = 4096;
const CHUNKS_PER_BATCH: u64 = 64;

/// Valid starts of every chromosome, as running totals.
struct Placement {
    /// `prefix[c]` is the number of valid starts before chromosome `c`
    /// (`prefix.len() == chromosomes + 1`).
    prefix: Vec<u64>,
}

impl Placement {
    fn new(genome: &Genome, length: u64) -> Result<Placement, UsageError> {
        let mut prefix = Vec::with_capacity(genome.len() + 1);
        prefix.push(0u64);
        let mut total = 0u64;
        let mut longest: Option<(usize, u64)> = None;
        for i in 0..genome.len() {
            let Some(len) = genome.length(i) else {
                return Err(UsageError(format!(
                    "-g: the genome has no length for {}; random needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
                    String::from_utf8_lossy(genome.name(i))
                )));
            };
            if longest.is_none_or(|(_, l)| len > l) {
                longest = Some((i, len));
            }
            if len >= length {
                total = total.checked_add(len - length + 1).ok_or_else(|| {
                    UsageError("-g: the chromosome lengths add up to more than 2^64 bases".into())
                })?;
            }
            prefix.push(total);
        }
        if total == 0 {
            return Err(match longest {
                Some((i, len)) => UsageError(format!(
                    "-l {length} is longer than every chromosome; the longest is {} ({len} bp), \
so use -l {len} or less",
                    String::from_utf8_lossy(genome.name(i))
                )),
                None => UsageError("-g: the genome has no chromosomes".into()),
            });
        }
        Ok(Placement { prefix })
    }

    fn total(&self) -> u64 {
        self.prefix.last().copied().unwrap_or(0)
    }

    /// The chromosome and start of valid position `r` (`r < total`).
    fn locate(&self, r: u64) -> (usize, u64) {
        // The last chromosome whose prefix is <= r (empty ones share a prefix
        // with their successor, so the partition point skips them).
        let c = self.prefix.partition_point(|&p| p <= r) - 1;
        (c, r - self.prefix[c])
    }
}

fn push_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(v.to_string().as_bytes());
}

/// Write `o.count` random intervals over `genome` as BED6 rows: chromosome,
/// start, end, id (1-based), length and strand, in id order (not sorted).
pub fn random_intervals(
    genome: &Genome,
    o: &RandomOptions,
    threads: NonZeroUsize,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    if o.length == 0 {
        return Err(UsageError("-l: the interval length must be at least 1".into()).into());
    }
    let placement = Placement::new(genome, o.length)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .build()
        .map_err(|e| IoError::Internal(format!("cannot start the thread pool: {e}")))?;
    let write_err = |source| IoError::Write {
        path: "output".into(),
        source,
    };
    let per_batch = CHUNK * CHUNKS_PER_BATCH;
    let mut done = 0u64;
    while done < o.count {
        let batch_end = o.count.min(done + per_batch);
        let chunks: Vec<(u64, u64)> = (done..batch_end)
            .step_by(CHUNK as usize)
            .map(|s| (s, batch_end.min(s + CHUNK)))
            .collect();
        let rows: Vec<Vec<Vec<u8>>> = pool.install(|| {
            chunks
                .par_iter()
                .map(|&(from, to)| {
                    let mut out = Vec::with_capacity((to - from) as usize);
                    for i in from..to {
                        let mut s = Stream::new(o.seed, i);
                        let (c, start) = placement.locate(s.below(placement.total()));
                        let plus = s.coin();
                        let mut row = Vec::with_capacity(48);
                        row.extend_from_slice(genome.name(c));
                        row.push(b'\t');
                        push_u64(&mut row, start);
                        row.push(b'\t');
                        push_u64(&mut row, start + o.length);
                        row.push(b'\t');
                        push_u64(&mut row, i + 1);
                        row.push(b'\t');
                        push_u64(&mut row, o.length);
                        row.extend_from_slice(if plus { b"\t+" } else { b"\t-" });
                        out.push(row);
                    }
                    out
                })
                .collect()
        });
        for row in rows.iter().flatten() {
            w.write_record(row).map_err(write_err)?;
        }
        done = batch_end;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome(sizes: &[(&str, u64)]) -> Genome {
        Genome::new(
            sizes
                .iter()
                .map(|(n, l)| (n.as_bytes().to_vec(), Some(*l)))
                .collect(),
        )
        .expect("genome")
    }

    #[test]
    fn locate_skips_short_chromosomes() {
        let g = genome(&[("a", 5), ("b", 100), ("c", 12)]);
        let p = Placement::new(&g, 10).expect("placement");
        // a: none, b: 91 starts, c: 3 starts.
        assert_eq!(p.total(), 94);
        assert_eq!(p.locate(0), (1, 0));
        assert_eq!(p.locate(90), (1, 90));
        assert_eq!(p.locate(91), (2, 0));
        assert_eq!(p.locate(93), (2, 2));
    }

    #[test]
    fn too_long_names_the_longest() {
        let g = genome(&[("chr1", 100), ("chr3", 700)]);
        let e = Placement::new(&g, 800).err().expect("error").0;
        assert!(e.contains("chr3 (700 bp)"), "{e}");
    }
}
