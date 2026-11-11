//! `gloci nuc`: base composition (and optionally the sequence and a pattern
//! count) under each interval (replaces `bedtools nuc`, D112).
//!
//! The intervals are read in input order and each sequence is streamed from
//! the indexed FASTA in chunks of at most [`gloci_io::FETCH_CHUNK`], so memory
//! does not depend on the interval length. bedtools prints the counts before
//! the `-seq` column, so with `--seq` one pass counts and a second pass
//! writes the sequence; a record is never held whole.
//!
//! Where bedtools is matched: counts are case-insensitive with `U` counted as
//! `T` and everything else as other; the percentages divide by the interval
//! length (so `N` and other bases lower them) in f32 and print with six
//! decimals; under `-s` a `-` record is reverse-complemented first, so the
//! counts, the sequence and the pattern count are all of the reverse
//! complement (a `U` there has become `A`); pattern matches overlap and are
//! case-sensitive unless `ignore_case`. The header line is bedtools' own,
//! written for the first record even when that record is skipped.

use crate::getfasta::{outside, plural, write_err};
use crate::seq::{SeqSummary, revcomp_in_place};
use gloci_core::{ChromTable, Strand};
use gloci_io::{
    BcfText, IndexedFasta, IoError, Loc, RecordBatch, RecordReader, column_count, next_batch_size,
};
use std::io::Write;

#[derive(Debug, Clone)]
pub struct NucOptions {
    /// Reverse-complement `-` records before counting.
    pub strand: bool,
    /// Add the sequence column.
    pub seq: bool,
    /// Add a column counting this pattern (overlapping matches).
    pub pattern: Option<Vec<u8>>,
    /// Match the pattern without regard to case.
    pub ignore_case: bool,
    /// `--skip-invalid`: skip (and count) an interval outside the FASTA
    /// instead of failing.
    pub skip_outside: bool,
}

/// Base counts of a sequence, case-insensitive; `U` counts as `T`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BaseCounts {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub other: u64,
}

impl BaseCounts {
    pub fn add(&mut self, chunk: &[u8]) {
        let mut k = [0u64; 6];
        for &b in chunk {
            k[CLASS[b as usize] as usize] += 1;
        }
        self.a += k[0];
        self.c += k[1];
        self.g += k[2];
        self.t += k[3];
        self.n += k[4];
        self.other += k[5];
    }

    pub fn total(&self) -> u64 {
        self.a + self.c + self.g + self.t + self.n + self.other
    }
}

const CLASS: [u8; 256] = {
    let mut t = [5u8; 256];
    let list: [(u8, u8); 6] = [
        (b'A', 0),
        (b'C', 1),
        (b'G', 2),
        (b'T', 3),
        (b'U', 3),
        (b'N', 4),
    ];
    let mut i = 0;
    while i < list.len() {
        t[list[i].0 as usize] = list[i].1;
        t[(list[i].0 + 32) as usize] = list[i].1;
        i += 1;
    }
    t
};

/// `count / len` in f32, printed with six decimals as bedtools' `%f` of a
/// float does (the f32 value is what is rounded, not the f64 quotient).
pub fn pct_f32(count: u64, len: u64) -> String {
    format!("{:.6}", count as f32 / len as f32)
}

/// Counts overlapping occurrences of a pattern across a stream of chunks,
/// keeping only the last `len - 1` bytes between chunks.
struct PatternCounter {
    pat: Vec<u8>,
    ignore_case: bool,
    carry: Vec<u8>,
    buf: Vec<u8>,
    count: u64,
}

impl PatternCounter {
    fn new(pat: &[u8], ignore_case: bool) -> Self {
        PatternCounter {
            pat: if ignore_case {
                pat.to_ascii_lowercase()
            } else {
                pat.to_vec()
            },
            ignore_case,
            carry: Vec::new(),
            buf: Vec::new(),
            count: 0,
        }
    }

    fn reset(&mut self) {
        self.carry.clear();
        self.count = 0;
    }

    fn feed(&mut self, chunk: &[u8]) {
        let plen = self.pat.len();
        if plen == 0 {
            return;
        }
        self.buf.clear();
        self.buf.extend_from_slice(&self.carry);
        self.buf.extend_from_slice(chunk);
        if self.buf.len() >= plen {
            let first = self.pat[0];
            let ic = self.ignore_case;
            for w in self.buf.windows(plen) {
                let head = if ic { w[0].to_ascii_lowercase() } else { w[0] };
                if head != first {
                    continue;
                }
                let hit = if ic {
                    w.iter()
                        .zip(&self.pat)
                        .all(|(a, b)| a.to_ascii_lowercase() == *b)
                } else {
                    w == self.pat.as_slice()
                };
                self.count += u64::from(hit);
            }
        }
        // A match inside the carry alone is impossible (it is shorter than
        // the pattern), so nothing is counted twice.
        let keep = (plen - 1).min(self.buf.len());
        self.carry.clear();
        self.carry
            .extend_from_slice(&self.buf[self.buf.len() - keep..]);
    }
}

fn header_line(k: usize, o: &NucOptions) -> String {
    let mut h = String::new();
    for i in 1..=k {
        h.push_str(&format!("{i}_usercol\t"));
    }
    let names = [
        "pct_at", "pct_gc", "num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len",
    ];
    let mut n = k;
    for name in names {
        n += 1;
        h.push_str(&format!("{n}_{name}\t"));
    }
    if o.seq {
        n += 1;
        h.push_str(&format!("{n}_seq\t"));
    }
    if o.pattern.is_some() {
        n += 1;
        h.push_str(&format!("{n}_user_patt_count\t"));
    }
    h.pop();
    format!("#{h}\n")
}

/// Write the base composition of every interval of `input`.
pub fn nuc(
    mut input: Box<dyn RecordReader>,
    fasta: &mut IndexedFasta,
    chroms: &mut ChromTable,
    o: &NucOptions,
    out: &mut dyn Write,
) -> Result<SeqSummary, IoError> {
    let name = input.name().to_string();
    let mut bcf = BcfText::for_reader(input.as_ref())?;
    let binary = bcf.is_some();
    let mut sum = SeqSummary::default();
    let mut first_outside: Option<String> = None;
    let mut n_outside = 0u64;
    let mut batch = RecordBatch::default();
    let mut scratch: Vec<u8> = Vec::new();
    let mut row: Vec<u8> = Vec::new();
    let mut counter = o
        .pattern
        .as_ref()
        .map(|p| PatternCounter::new(p, o.ignore_case));
    let mut header_done = false;
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let at = |msg: String| IoError::Parse {
                path: name.clone(),
                loc: if binary {
                    Loc::Record(m.loc)
                } else {
                    Loc::Line(m.loc)
                },
                msg,
            };
            let line = match &mut bcf {
                Some(b) => b.render(batch.raw(i), m.loc)?,
                None => batch.raw(i),
            };
            if !header_done {
                header_done = true;
                out.write_all(header_line(column_count(line), o).as_bytes())
                    .map_err(write_err)?;
            }
            if m.span.is_point() {
                sum.points += 1;
                continue;
            }
            let chrom = chroms.name(m.chrom);
            if let Some(why) = outside(fasta, chrom, m.span) {
                first_outside.get_or_insert_with(|| why.clone());
                if !o.skip_outside {
                    return Err(at(format!(
                        "{why}; fix the interval or the FASTA, or pass --skip-invalid to skip \
intervals outside the FASTA"
                    )));
                }
                n_outside += 1;
                sum.skipped += 1;
                continue;
            }
            let minus = o.strand && m.strand == Strand::Minus;
            // Pass 1: counts (and the pattern) of the sequence as written.
            let mut counts = BaseCounts::default();
            if let Some(c) = counter.as_mut() {
                c.reset();
            }
            let mut count_chunk = |chunk: &[u8]| -> Result<(), IoError> {
                let seen = if minus {
                    scratch.clear();
                    scratch.extend_from_slice(chunk);
                    revcomp_in_place(&mut scratch, false);
                    scratch.as_slice()
                } else {
                    chunk
                };
                counts.add(seen);
                if let Some(c) = counter.as_mut() {
                    c.feed(seen);
                }
                Ok(())
            };
            if minus {
                fasta.fetch_rev(chrom, m.span, &mut count_chunk)?;
            } else {
                fasta.fetch(chrom, m.span, &mut count_chunk)?;
            }
            let len = counts.total();
            row.clear();
            row.extend_from_slice(line);
            row.extend_from_slice(
                format!(
                    "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    pct_f32(counts.a + counts.t, len),
                    pct_f32(counts.c + counts.g, len),
                    counts.a,
                    counts.c,
                    counts.g,
                    counts.t,
                    counts.n,
                    counts.other,
                    len
                )
                .as_bytes(),
            );
            out.write_all(&row).map_err(write_err)?;
            if o.seq {
                // Pass 2: the sequence, in the same order as pass 1.
                out.write_all(b"\t").map_err(write_err)?;
                let mut emit = |chunk: &[u8]| -> Result<(), IoError> {
                    if minus {
                        scratch.clear();
                        scratch.extend_from_slice(chunk);
                        revcomp_in_place(&mut scratch, false);
                        out.write_all(&scratch)
                    } else {
                        out.write_all(chunk)
                    }
                    .map_err(write_err)
                };
                if minus {
                    fasta.fetch_rev(chrom, m.span, &mut emit)?;
                } else {
                    fasta.fetch(chrom, m.span, &mut emit)?;
                }
            }
            if let Some(c) = &counter {
                out.write_all(format!("\t{}", c.count).as_bytes())
                    .map_err(write_err)?;
            }
            out.write_all(b"\n").map_err(write_err)?;
            sum.written += 1;
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    if let Some(s) = input.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    if sum.points > 0 {
        log::warn!(
            "skipped {} zero-length interval{} (they have no sequence to count)",
            sum.points,
            plural(sum.points)
        );
    }
    if let Some(first) = first_outside {
        log::warn!(
            "skipped {n_outside} interval{} outside the FASTA (first: {first})",
            plural(n_outside)
        );
    }
    Ok(sum)
}
