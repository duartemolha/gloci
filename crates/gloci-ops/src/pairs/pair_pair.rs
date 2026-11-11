//! `pair-pair` (D155, bedtools pairtopair): the pairs of A, by how their
//! ends overlap the ends of the pairs of B.
//!
//! B's end 1 and end 2 each get an [`IntervalIndex`] (a B record's id is its
//! input position in both). Orientation 1 tests A1-B1 and A2-B2, orientation
//! 2 tests A1-B2 and A2-B1. A B is a hit of `both` when one orientation hits
//! on both ends, of `either` when any of the four tests hits. A hit B is
//! written once, as its line reads, in B's input order. The batching and
//! thread handling are `pairs::run_a`'s.

use super::{Counts, PairSummary, PairWorker, run_a};
use gloci_core::{ChromId, ChromTable, Interval, OverlapCriteria, Strand, StrandRule, Window};
use gloci_engine::{IndexEntry, IntervalIndex};
use gloci_io::{BedpeBatch, BedpeMeta, BedpeReader, IoError, Loc, bedpe_name};
use std::io::Write;
use std::num::NonZeroUsize;

/// What `pair-pair` reports (bedtools `-type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairPairType {
    /// A alone, when no end of A overlaps an end of B.
    Neither,
    /// A with each B that has an end overlapping an end of A.
    Either,
    /// A with each B whose two ends overlap the two ends of A.
    Both,
    /// A alone, when no B has both ends overlapping.
    NotBoth,
}

impl PairPairType {
    /// Whether the report is A's pairs alone, which are BEDPE lines.
    pub fn writes_a_alone(self) -> bool {
        matches!(self, PairPairType::Neither | PairPairType::NotBoth)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PairPairOptions {
    pub kind: PairPairType,
    pub criteria: OverlapCriteria,
    /// B's ends are padded by this before they are tested.
    pub window: Option<Window>,
    /// Drop a B named like A (column 7).
    pub different_names: bool,
}

type End = Option<(Interval, Strand)>;

/// B's records after loading.
struct BStore {
    /// Per end: the index over that end of every record.
    index: [IntervalIndex; 2],
    /// Per id: both ends.
    meta: Vec<[End; 2]>,
    /// Per id: offset and length of its line in `text`.
    at: Vec<(usize, u32)>,
    text: Vec<u8>,
    unknown_strand: bool,
}

impl BStore {
    fn line(&self, id: u32) -> &[u8] {
        let (off, len) = self.at[id as usize];
        &self.text[off..off + len as usize]
    }
}

fn load_b(
    mut b: BedpeReader,
    chroms: &mut ChromTable,
    different_names: bool,
) -> Result<BStore, IoError> {
    let mut store = BStore {
        index: [IntervalIndex::default(), IntervalIndex::default()],
        meta: Vec::new(),
        at: Vec::new(),
        text: Vec::new(),
        unknown_strand: false,
    };
    let mut entries: [Vec<IndexEntry>; 2] = [Vec::new(), Vec::new()];
    let (mut batch, mut max) = (BedpeBatch::default(), 1usize);
    loop {
        batch.clear();
        let more = b.read_batch(chroms, &mut batch, max)?;
        if !batch.is_empty() {
            max = ((4usize << 20) / (batch.buf.len() / batch.len() + 1)).clamp(1, 16384);
        }
        for (i, m) in batch.recs.iter().enumerate() {
            let id = u32::try_from(store.meta.len()).map_err(|_| {
                IoError::from(gloci_io::UsageError(
                    "B has more than 4 billion records; split it".into(),
                ))
            })?;
            let line = batch.raw(i);
            if different_names && bedpe_name(line).is_none() {
                return Err(no_name(b.name(), m.loc, "B"));
            }
            store.at.push((store.text.len(), line.len() as u32));
            store.text.extend_from_slice(line);
            let mut ends = [None, None];
            for (k, e) in m.ends.iter().enumerate() {
                if let Some(e) = e {
                    ends[k] = Some((e.span, e.strand));
                    store.unknown_strand |= e.strand == Strand::Unknown;
                    entries[k].push(IndexEntry {
                        chrom: e.chrom,
                        span: e.span,
                        id,
                    });
                }
            }
            store.meta.push(ends);
        }
        if !more {
            break;
        }
    }
    if let Some(s) = b.diagnostics().summary(b.name()) {
        log::warn!("{s}");
    }
    let [e1, e2] = entries;
    store.index = [IntervalIndex::build(e1), IntervalIndex::build(e2)];
    Ok(store)
}

fn no_name(path: &str, line: u64, which: &str) -> IoError {
    IoError::Parse {
        path: path.to_string(),
        loc: Loc::Line(line),
        msg: format!(
            "{which} pair has no name (column 7), which --require-different-names compares; \
give every pair a name or drop the option"
        ),
    }
}

#[derive(Default)]
struct Scratch {
    cand: Vec<u32>,
    /// `h[a_end][b_end]`: the B records whose `b_end` meets A's `a_end`.
    h: [[Vec<u32>; 2]; 2],
    both: Vec<u32>,
    any: Vec<u32>,
}

struct Worker<'a> {
    a_name: &'a str,
    b: &'a BStore,
    o: &'a PairPairOptions,
}

/// The ids in both of two ascending lists, appended to `out`.
fn intersect_into(x: &[u32], y: &[u32], out: &mut Vec<u32>) {
    let (mut i, mut j) = (0, 0);
    while i < x.len() && j < y.len() {
        match x[i].cmp(&y[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(x[i]);
                i += 1;
                j += 1;
            }
        }
    }
}

impl Worker<'_> {
    /// The B records whose end `b_end` meets `span` of `chrom`, ascending.
    fn hits(
        &self,
        chrom: ChromId,
        span: Interval,
        strand: Strand,
        b_end: usize,
        s: &mut Scratch,
    ) -> Vec<u32> {
        let c = &self.o.criteria;
        self.b.index[b_end].query(
            chrom,
            span,
            self.o.window,
            strand,
            c.inclusive_points,
            &mut s.cand,
        );
        s.cand
            .iter()
            .copied()
            .filter(|&id| {
                let (bs, bstrand) =
                    self.b.meta[id as usize][b_end].expect("indexed end is present");
                let bs = self.o.window.map_or(bs, |w| w.apply(bs, strand));
                c.accepts(span, strand, bs, bstrand)
            })
            .collect()
    }

    fn pair(
        &self,
        raw: &[u8],
        m: &BedpeMeta,
        s: &mut Scratch,
        out: &mut Vec<u8>,
        n: &mut Counts,
    ) -> Result<(), IoError> {
        let strand_rule = self.o.criteria.strand != StrandRule::Any;
        let a_name = if self.o.different_names {
            Some(bedpe_name(raw).ok_or_else(|| no_name(self.a_name, m.loc, "A"))?)
        } else {
            None
        };
        for (k, end) in m.ends.iter().enumerate() {
            for j in 0..2 {
                s.h[k][j] = match end {
                    Some(e) => {
                        n.a_unknown |= strand_rule && e.strand == Strand::Unknown;
                        let mut ids = self.hits(e.chrom, e.span, e.strand, j, s);
                        if let Some(name) = a_name {
                            ids.retain(|&id| bedpe_name(self.b.line(id)) != Some(name));
                        }
                        ids
                    }
                    None => Vec::new(),
                };
            }
        }
        s.both.clear();
        intersect_into(&s.h[0][0], &s.h[1][1], &mut s.both);
        intersect_into(&s.h[0][1], &s.h[1][0], &mut s.both);
        s.both.sort_unstable();
        s.both.dedup();
        let ids: &[u32] = match self.o.kind {
            PairPairType::Both => &s.both,
            _ => {
                s.any.clear();
                for h in s.h.iter().flatten() {
                    s.any.extend_from_slice(h);
                }
                s.any.sort_unstable();
                s.any.dedup();
                &s.any
            }
        };
        match self.o.kind {
            PairPairType::Both | PairPairType::Either => {
                for &id in ids {
                    out.extend_from_slice(raw);
                    out.push(b'\t');
                    out.extend_from_slice(self.b.line(id));
                    out.push(b'\n');
                }
            }
            PairPairType::Neither if ids.is_empty() => {
                out.extend_from_slice(raw);
                out.push(b'\n');
            }
            PairPairType::NotBoth if s.both.is_empty() => {
                out.extend_from_slice(raw);
                out.push(b'\n');
            }
            _ => {}
        }
        Ok(())
    }
}

impl PairWorker for Worker<'_> {
    fn batch(&self, batch: &BedpeBatch, out: &mut Vec<u8>) -> Result<Counts, IoError> {
        let mut n = Counts::default();
        let mut s = Scratch::default();
        for (i, m) in batch.recs.iter().enumerate() {
            n.pairs += 1;
            self.pair(batch.raw(i), m, &mut s, out, &mut n)?;
        }
        Ok(n)
    }
}

/// Run `pair-pair`: A's pairs against B's, written to `out` (TSV rows for
/// A+B, BEDPE lines for A alone). `chroms` interns both inputs' chromosomes.
/// With `threads > 1` the output is byte-identical to `threads == 1`.
pub fn pair_pair(
    mut a: BedpeReader,
    b: BedpeReader,
    chroms: &mut ChromTable,
    o: &PairPairOptions,
    threads: NonZeroUsize,
    out: &mut dyn Write,
) -> Result<PairSummary, IoError> {
    let store = load_b(b, chroms, o.different_names)?;
    let a_name = a.name().to_string();
    let worker = Worker {
        a_name: &a_name,
        b: &store,
        o,
    };
    let (summary, a_unknown) = run_a(&mut a, chroms, &worker, threads, out, "pair-pair")?;
    if o.criteria.strand != StrandRule::Any {
        if a_unknown {
            log::warn!("{a_name} has pair ends of unknown strand; they never match under -s or -S");
        }
        if store.unknown_strand {
            log::warn!("B has pair ends of unknown strand; they never match under -s or -S");
        }
    }
    Ok(summary)
}
