//! `merge`: cluster the footprints of sorted inputs and write one BED row
//! per cluster (replaces `bedtools merge` and `bedops -m`).

use crate::agg::{Accumulators, AggSpec, Op, ResolvedAgg, resolve};
use gloci_core::{ChromId, ChromTable, Strand};
use gloci_engine::kmerge_with;
use gloci_io::{Format, IoError, OutFormat, RecordMeta, RecordReader, RecordWriter, UsageError};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(Debug, Clone)]
pub struct MergeOptions {
    /// Records join the open cluster when `start <= cluster_end + max_gap`
    /// (signed): `0` joins overlapping and bookended features, `-N` needs at
    /// least N shared bases.
    pub max_gap: i64,
    /// Cluster each strand separately (`+`, `-` and unknown).
    pub per_strand: bool,
    /// Add a column with the number of records in each cluster.
    pub count: bool,
    pub aggs: Vec<AggSpec>,
    /// Value written for an aggregation with no values.
    pub empty: String,
    /// Separator of list-valued aggregations (`collapse`, `distinct`).
    pub delim: String,
}

/// Reject option and input combinations that cannot work, before any input
/// is read.
pub fn validate(formats: &[Format], o: &MergeOptions) -> Result<(), UsageError> {
    if let Some(a) = o.aggs.iter().find(|a| a.op == Op::Bases) {
        return Err(UsageError(format!(
            "-c {}: bases measures overlap with another input, and merge has none; \
use count for the number of records",
            a.label
        )));
    }
    if !o.aggs.is_empty() {
        if formats.contains(&Format::Bcf) {
            return Err(UsageError(
                "-c cannot read BCF records; convert first: gloci view FILE -O vcf".into(),
            ));
        }
        let family = |f: Format| OutFormat::family(f);
        if formats.iter().any(|&f| family(f) != family(formats[0])) {
            return Err(UsageError(
                "-c needs inputs of one format, because columns are read by position".into(),
            ));
        }
    }
    Ok(())
}

/// A footprint waiting for its turn in start order.
struct Pending {
    start: u64,
    end: u64,
    slot: usize,
    seq: u64,
    input: usize,
    loc: u64,
    /// The record's text, kept only when there are aggregations.
    raw: Vec<u8>,
}

impl Pending {
    fn key(&self) -> (u64, u64) {
        (self.start, self.seq)
    }
}

impl PartialEq for Pending {
    fn eq(&self, o: &Self) -> bool {
        self.key() == o.key()
    }
}
impl Eq for Pending {}
impl PartialOrd for Pending {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Pending {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.key().cmp(&o.key())
    }
}

struct Open {
    start: u64,
    end: u64,
    count: u64,
}

/// A finished cluster row, ordered by (start, end, strand slot, creation).
struct Closed {
    order: (u64, u64, usize, u64),
    row: Vec<u8>,
}

impl PartialEq for Closed {
    fn eq(&self, o: &Self) -> bool {
        self.order == o.order
    }
}
impl Eq for Closed {}
impl PartialOrd for Closed {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Closed {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.order.cmp(&o.order)
    }
}

fn slot_of(per_strand: bool, s: Strand) -> usize {
    if !per_strand {
        return 0;
    }
    match s {
        Strand::Plus => 0,
        Strand::Minus => 1,
        Strand::Unknown => 2,
    }
}

struct Merger<'a> {
    o: &'a MergeOptions,
    names: Vec<String>,
    w: &'a mut dyn RecordWriter,
    heap: BinaryHeap<Reverse<Pending>>,
    open: Vec<Option<Open>>,
    accs: Vec<Accumulators>,
    closed: BinaryHeap<Reverse<Closed>>,
    chrom: Option<ChromId>,
    chrom_name: Vec<u8>,
    seq: u64,
    created: u64,
    keep_raw: bool,
    /// The most closed clusters held back at once.
    peak_closed: usize,
}

impl Merger<'_> {
    fn record(
        &mut self,
        input: usize,
        m: &RecordMeta,
        raw: &[u8],
        chroms: &ChromTable,
    ) -> Result<(), IoError> {
        if self.chrom != Some(m.chrom) {
            self.finish_chrom()?;
            self.chrom = Some(m.chrom);
            self.chrom_name = chroms.name(m.chrom).to_vec();
        }
        // A span never starts before its key, so everything starting at or
        // before this key is final in start order.
        while self.heap.peek().is_some_and(|p| p.0.start <= m.key) {
            let Reverse(p) = self.heap.pop().expect("peeked");
            self.add(p)?;
        }
        let slot = slot_of(self.o.per_strand, m.strand);
        self.heap.push(Reverse(Pending {
            start: m.span.start,
            end: m.span.end,
            slot,
            seq: self.seq,
            input,
            loc: m.loc,
            raw: if self.keep_raw {
                raw.to_vec()
            } else {
                Vec::new()
            },
        }));
        self.seq += 1;
        // No later span starts before this key, so an open cluster that ends
        // more than the gap before it can never grow: close it now, or one
        // stale cluster of another strand would hold every closed cluster
        // back until the chromosome ends.
        for slot in 0..self.open.len() {
            if self.open[slot]
                .as_ref()
                .is_some_and(|c| i128::from(c.end) + i128::from(self.o.max_gap) < i128::from(m.key))
            {
                self.close(slot);
            }
        }
        let horizon = self
            .open
            .iter()
            .flatten()
            .map(|c| c.start)
            .min()
            .map_or(m.key, |s| s.min(m.key));
        self.flush(Some(horizon))
    }

    fn add(&mut self, p: Pending) -> Result<(), IoError> {
        let slot = p.slot;
        let joins = self.open[slot]
            .as_ref()
            .is_some_and(|c| i128::from(p.start) <= i128::from(c.end) + i128::from(self.o.max_gap));
        if !joins {
            self.close(slot);
            self.open[slot] = Some(Open {
                start: p.start,
                end: p.end,
                count: 0,
            });
        }
        let c = self.open[slot].as_mut().expect("opened");
        c.end = c.end.max(p.end);
        c.count += 1;
        if self.keep_raw {
            self.accs[slot].add(&p.raw, 0, &self.names[p.input], p.loc)?;
        }
        Ok(())
    }

    fn close(&mut self, slot: usize) {
        let Some(c) = self.open[slot].take() else {
            return;
        };
        let mut row = self.chrom_name.clone();
        row.push(b'\t');
        row.extend_from_slice(c.start.to_string().as_bytes());
        row.push(b'\t');
        row.extend_from_slice(c.end.to_string().as_bytes());
        if self.o.count {
            row.push(b'\t');
            row.extend_from_slice(c.count.to_string().as_bytes());
        }
        self.accs[slot].finish_into(&mut row);
        self.created += 1;
        self.closed.push(Reverse(Closed {
            order: (c.start, c.end, slot, self.created),
            row,
        }));
        self.peak_closed = self.peak_closed.max(self.closed.len());
    }

    /// Write closed clusters that no later cluster can precede: those that
    /// start before `horizon` (all of them when `None`).
    fn flush(&mut self, horizon: Option<u64>) -> Result<(), IoError> {
        while let Some(Reverse(top)) = self.closed.peek() {
            if horizon.is_some_and(|h| top.order.0 >= h) {
                break;
            }
            let Reverse(c) = self.closed.pop().expect("peeked");
            self.w
                .write_record(&c.row)
                .map_err(|source| IoError::Write {
                    path: "output".into(),
                    source,
                })?;
        }
        Ok(())
    }

    fn finish_chrom(&mut self) -> Result<(), IoError> {
        while let Some(Reverse(p)) = self.heap.pop() {
            self.add(p)?;
        }
        for slot in 0..self.open.len() {
            self.close(slot);
        }
        self.flush(None)
    }
}

/// Merge `inputs` (each sorted in `chroms` order) and write one row per
/// cluster: `chrom start end`, the record count with `o.count`, then one
/// column per aggregation.
pub fn merge(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    o: &MergeOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    merge_counting(inputs, chroms, o, w).map(|_| ())
}

/// [`merge`], returning the most finished clusters ever held back at once
/// (for the memory-flatness test).
#[doc(hidden)]
pub fn merge_counting(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    o: &MergeOptions,
    w: &mut dyn RecordWriter,
) -> Result<usize, IoError> {
    let formats: Vec<Format> = inputs.iter().map(|r| r.format()).collect();
    validate(&formats, o)?;
    let names: Vec<String> = inputs.iter().map(|r| r.name().to_string()).collect();
    let resolved: Vec<ResolvedAgg> = match inputs.first() {
        Some(first) if !o.aggs.is_empty() => resolve(&o.aggs, first.format(), first.header())?,
        _ => Vec::new(),
    };
    let slots = if o.per_strand { 3 } else { 1 };
    let accs = (0..slots)
        .map(|_| Accumulators::new(resolved.clone(), o.empty.clone(), o.delim.clone()))
        .collect();
    let mut unknown_strand = vec![false; inputs.len()];
    let mut m = Merger {
        o,
        names: names.clone(),
        w,
        heap: BinaryHeap::new(),
        open: (0..slots).map(|_| None).collect(),
        accs,
        closed: BinaryHeap::new(),
        chrom: None,
        chrom_name: Vec::new(),
        seq: 0,
        created: 0,
        keep_raw: !o.aggs.is_empty(),
        peak_closed: 0,
    };
    kmerge_with(inputs, chroms, &mut |input, meta, raw, chroms| {
        if meta.strand == Strand::Unknown {
            unknown_strand[input] = true;
        }
        m.record(input, meta, raw, chroms)
    })?;
    m.finish_chrom()?;
    if o.per_strand {
        for (name, _) in names.iter().zip(&unknown_strand).filter(|(_, u)| **u) {
            log::warn!(
                "{name}: some records have no strand; with -s they are merged together, \
apart from + and - records"
            );
        }
    }
    Ok(m.peak_closed)
}
