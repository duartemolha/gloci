//! `stats reldist`: the distribution of the relative distance of each A
//! midpoint to the nearest pair of B midpoints (replaces `bedtools reldist`).

use crate::intersect::Text;
use crate::relate::{RelateInputs, RelateOptions, relation_chroms};
use crate::skipped::{Skipped, log_skipped};
use crate::sort::{SortSettings, sort_readers};
use gloci_core::{ChromId, ChromOrder, ChromTable, Interval};
use gloci_io::{
    IoError, OrderGuard, RecordBatch, RecordMeta, RecordReader, RecordWriter, UsageError,
};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};

/// Records to request per read.
const READ_MAX: usize = 8192;

#[derive(Debug, Clone)]
pub struct ReldistOptions {
    /// Print each qualifying A record with its relative distance instead of the summary.
    pub detail: bool,
    /// Accept inputs in any order: all B midpoints are held in memory.
    pub unsorted: bool,
    /// `--sort`: sort A and B before reading them.
    pub sort: Option<SortSettings>,
    pub genome: Option<String>,
    pub chrom_order: ChromOrder,
}

/// The summary: how many A records fall in each bin of `⌊rel·100⌋`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReldistSummary {
    pub bins: BTreeMap<u64, u64>,
    pub total: u64,
}

/// The midpoint of an A interval: `⌊(start + end) / 2⌋`.
fn a_mid(s: Interval) -> i64 {
    (s.start + (s.end - s.start) / 2) as i64
}

/// The midpoint of a B interval: `⌊(start + end) / 2⌋` without overflow. bedtools adds the
/// two ends in a 32-bit `int`, which wraps past 2^31 (`#reldist-midpoint-overflow`).
fn b_mid(s: Interval) -> i64 {
    (s.start + (s.end - s.start) / 2) as i64
}

/// The relative distance of `m` to the B midpoints `mids` (sorted): `None`
/// when bedtools skips it (a chromosome with fewer than two midpoints, `m`
/// before the first or after the last midpoint).
fn rel_dist(mids: &[i64], m: i64) -> Option<f64> {
    if mids.len() < 2 {
        return None;
    }
    let i = mids.partition_point(|&x| x < m);
    let right = *mids.get(i)?;
    if right == m {
        return Some(0.0);
    }
    let left = *mids.get(i.checked_sub(1)?)?;
    let (dl, dr) = (m - left, right - m);
    Some(dl.min(dr) as f64 / (right - left) as f64)
}

struct Stream {
    r: Box<dyn RecordReader>,
    guard: Option<OrderGuard>,
    batch: RecordBatch,
    pos: usize,
    eof: bool,
}

impl Stream {
    fn new(r: Box<dyn RecordReader>, ordered: bool) -> Stream {
        let guard = ordered.then(|| OrderGuard::new(r.name(), r.format()));
        Stream {
            r,
            guard,
            batch: RecordBatch::default(),
            pos: 0,
            eof: false,
        }
    }

    /// Make the record at `pos` available; false at the end of the input.
    fn fill(&mut self, chroms: &mut ChromTable) -> Result<bool, IoError> {
        while self.pos >= self.batch.len() {
            if self.eof {
                return Ok(false);
            }
            self.batch.clear();
            self.pos = 0;
            if !self.r.read_batch(chroms, &mut self.batch, READ_MAX)? {
                self.eof = true;
            }
            if let Some(g) = &mut self.guard {
                g.check(chroms, &self.batch)?;
            }
        }
        Ok(true)
    }

    fn meta(&self) -> RecordMeta {
        self.batch.recs[self.pos]
    }
}

/// B midpoints of chromosome `c`, read from the sorted stream `b` (the
/// records of earlier chromosomes are dropped).
fn load_chrom(
    b: &mut Stream,
    c: ChromId,
    chroms: &mut ChromTable,
    mids: &mut Vec<i64>,
) -> Result<(), IoError> {
    mids.clear();
    while b.fill(chroms)? {
        let m = b.meta();
        match chroms.cmp(m.chrom, c) {
            Ordering::Less => {}
            Ordering::Equal => mids.push(b_mid(m.span)),
            Ordering::Greater => break,
        }
        b.pos += 1;
    }
    mids.sort_unstable();
    Ok(())
}

/// Bin the relative distance of every A record, writing `--detail` rows to `detail`.
pub fn reldist(
    inputs: RelateInputs,
    o: &ReldistOptions,
    mut detail: Option<&mut dyn RecordWriter>,
) -> Result<ReldistSummary, IoError> {
    let RelateInputs { a, b } = inputs;
    if b.len() != 1 {
        return Err(UsageError(format!(
            "reldist compares two files, A and one B, but {} B files were given; \
run it once per B file",
            b.len()
        ))
        .into());
    }
    if o.sort.is_some() && o.unsorted {
        return Err(UsageError("--sort and --unsorted cannot be combined".into()).into());
    }
    let relate = RelateOptions {
        criteria: gloci_core::OverlapCriteria::default(),
        window: None,
        split: false,
        split_deletions: false,
        unsorted: o.unsorted,
        sort: None,
        labels: Vec::new(),
        genome: o.genome.clone(),
        chrom_order: o.chrom_order,
        header: false,
    };
    let inputs = RelateInputs { a, b };
    let mut chroms = relation_chroms(&inputs, &relate)?;
    let RelateInputs { a, b } = inputs;
    let mut skipped = Skipped::default();
    let mut prepare = |r: Box<dyn RecordReader>| -> Result<Box<dyn RecordReader>, IoError> {
        match &o.sort {
            Some(s) => {
                sort_readers(vec![r], &mut chroms, s).map_err(IoError::with_sorted_input_loc)
            }
            None => Ok(skipped.watch(r)),
        }
    };
    let a = prepare(a)?;
    let b = prepare(b.into_iter().next().expect("one B"))?;
    let mut text = Text::new(a.as_ref())?;
    // A sorted stream is in the table's order; otherwise the inputs are checked.
    let ordered = !o.unsorted && o.sort.is_none();
    let mut a = Stream::new(a, ordered);
    let mut b = Stream::new(b, ordered);
    let mut summary = ReldistSummary::default();
    let mut row = Vec::new();
    let mut mids: Vec<i64> = Vec::new();
    let mut by_chrom: HashMap<ChromId, Vec<i64>> = HashMap::new();
    if o.unsorted {
        while b.fill(&mut chroms)? {
            let m = b.meta();
            by_chrom.entry(m.chrom).or_default().push(b_mid(m.span));
            b.pos += 1;
        }
        by_chrom.values_mut().for_each(|v| v.sort_unstable());
    }
    let mut cur: Option<ChromId> = None;
    while a.fill(&mut chroms)? {
        let m = a.meta();
        if !o.unsorted && cur != Some(m.chrom) {
            cur = Some(m.chrom);
            load_chrom(&mut b, m.chrom, &mut chroms, &mut mids)?;
        }
        let list: &[i64] = if o.unsorted {
            by_chrom.get(&m.chrom).map_or(&[], Vec::as_slice)
        } else {
            &mids
        };
        let rel = rel_dist(list, a_mid(m.span));
        if let Some(rel) = rel {
            *summary
                .bins
                .entry((rel * 100.0).floor() as u64)
                .or_insert(0) += 1;
            summary.total += 1;
            if let Some(w) = detail.as_deref_mut() {
                row.clear();
                row.extend_from_slice(text.of(a.batch.raw(a.pos), m.loc)?);
                row.push(b'\t');
                gloci_core::cfmt::write_fixed(&mut row, rel, 3);
                w.write_record(&row).map_err(|source| IoError::Write {
                    path: "output".to_string(),
                    source,
                })?;
            }
        }
        a.pos += 1;
    }
    // Read B to its end, so that a bad record after A's last chromosome is still found.
    if !o.unsorted {
        while b.fill(&mut chroms)? {
            b.pos = b.batch.len();
        }
    }
    drop((a, b));
    log_skipped(&skipped.take());
    skipped.log_alignment_counts();
    Ok(summary)
}

/// The summary table, as bedtools prints it.
pub fn write_reldist(s: &ReldistSummary, out: &mut Vec<u8>) {
    out.extend_from_slice(b"reldist\tcount\ttotal\tfraction\n");
    for (&bin, &count) in &s.bins {
        gloci_core::cfmt::write_fixed(out, bin as f64 / 100.0, 2);
        out.push(b'\t');
        gloci_core::cfmt::write_u64(out, count);
        out.push(b'\t');
        gloci_core::cfmt::write_u64(out, s.total);
        out.push(b'\t');
        gloci_core::cfmt::write_fixed(out, count as f64 / s.total as f64, 3);
        out.push(b'\n');
    }
}
