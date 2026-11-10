//! `cluster`: write each record with the id of the cluster it belongs to
//! (replaces `bedtools cluster`).
//!
//! Clusters are those of `merge` (same gap and strand rules, shared with it),
//! but every input record is kept: the output is the input line, a tab and the
//! cluster id. Ids start at 1 and increase through the file, across
//! chromosomes; with strands, a new cluster takes the next id when its first
//! record is read.

use crate::merge::{joins, slot_of};
use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable, Strand};
use gloci_engine::kmerge_with;
use gloci_io::{Format, IoError, RecordMeta, RecordReader, RecordWriter, UsageError};
use std::cell::Cell;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};
use std::rc::Rc;

#[derive(Debug, Clone)]
pub struct ClusterOptions {
    /// As `merge --max-gap`: a record joins the open cluster when
    /// `start <= cluster_end + max_gap` (signed).
    pub max_gap: i64,
    /// Cluster each strand separately (`+`, `-` and unknown).
    pub per_strand: bool,
}

/// Reject input formats whose records are not text lines.
pub fn validate(format: Format) -> Result<(), UsageError> {
    if format == Format::Bcf {
        return Err(UsageError(
            "cluster appends a column to each input line, which BCF records are not; convert \
             first: gloci view FILE -O vcf | gloci cluster -"
                .into(),
        ));
    }
    Ok(())
}

/// A cluster; its id is taken when its first record (in input order) is
/// written (0 until then).
type Group = Rc<Cell<u64>>;

/// The open cluster of one strand slot.
struct Open {
    end: u64,
    group: Group,
}

/// A record waiting to be written in input order.
struct Held {
    raw: Vec<u8>,
    start: u64,
    end: u64,
    slot: usize,
    group: Option<Group>,
}

/// Groups records exactly as `merge` does (footprints taken in start order,
/// through a window that holds a record until nothing can start before it),
/// and writes them in input order with ids in order of first appearance.
struct Clusterer<'a> {
    o: &'a ClusterOptions,
    w: &'a mut dyn RecordWriter,
    /// (start, seq) of the records not yet grouped.
    heap: BinaryHeap<Reverse<(u64, u64)>>,
    /// The records from sequence number `base` on, until written.
    queue: VecDeque<Held>,
    base: u64,
    open: [Option<Open>; 3],
    chrom: Option<ChromId>,
    next_id: u64,
    row: Vec<u8>,
}

impl Clusterer<'_> {
    fn record(&mut self, m: &RecordMeta, raw: &[u8]) -> Result<(), IoError> {
        if self.chrom != Some(m.chrom) {
            self.finish_chrom()?;
            self.chrom = Some(m.chrom);
        }
        // A span never starts before its key, so everything starting at or
        // before this key is final in start order.
        while self.heap.peek().is_some_and(|p| p.0.0 <= m.key) {
            let Reverse((_, seq)) = self.heap.pop().expect("peeked");
            self.add(seq);
        }
        let seq = self.base + self.queue.len() as u64;
        self.queue.push_back(Held {
            raw: raw.to_vec(),
            start: m.span.start,
            end: m.span.end,
            slot: slot_of(self.o.per_strand, m.strand),
            group: None,
        });
        self.heap.push(Reverse((m.span.start, seq)));
        self.flush()
    }

    /// Put record `seq` (start order) into the open cluster of its strand, or
    /// open a new one.
    fn add(&mut self, seq: u64) {
        let idx = (seq - self.base) as usize;
        let Held {
            start, end, slot, ..
        } = self.queue[idx];
        let group = match &mut self.open[slot] {
            Some(c) if joins(start, c.end, self.o.max_gap) => {
                c.end = c.end.max(end);
                c.group.clone()
            }
            open => {
                let group: Group = Rc::new(Cell::new(0));
                *open = Some(Open {
                    end,
                    group: group.clone(),
                });
                group
            }
        };
        self.queue[idx].group = Some(group);
    }

    /// Write the records at the front whose cluster is known.
    fn flush(&mut self) -> Result<(), IoError> {
        while self.queue.front().is_some_and(|h| h.group.is_some()) {
            let h = self.queue.pop_front().expect("peeked");
            self.base += 1;
            let group = h.group.expect("checked");
            if group.get() == 0 {
                self.next_id += 1;
                group.set(self.next_id);
            }
            self.row.clear();
            self.row.extend_from_slice(&h.raw);
            self.row.push(b'\t');
            self.row
                .extend_from_slice(group.get().to_string().as_bytes());
            self.w
                .write_record(&self.row)
                .map_err(|source| IoError::Write {
                    path: "output".into(),
                    source,
                })?;
        }
        Ok(())
    }

    fn finish_chrom(&mut self) -> Result<(), IoError> {
        while let Some(Reverse((_, seq))) = self.heap.pop() {
            self.add(seq);
        }
        self.open = [None, None, None];
        self.flush()
    }
}

/// Write every record of the sorted `input` followed by a tab and its
/// cluster id.
pub fn cluster(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &ClusterOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    validate(input.format())?;
    let mut skipped = Skipped::default();
    let input = skipped.watch(input);
    let mut unknown_strand = false;
    let mut c = Clusterer {
        o,
        w,
        heap: BinaryHeap::new(),
        queue: VecDeque::new(),
        base: 0,
        open: [None, None, None],
        chrom: None,
        next_id: 0,
        row: Vec::new(),
    };
    kmerge_with(vec![input], chroms, &mut |_, m, raw, _| {
        unknown_strand |= m.strand == Strand::Unknown;
        c.record(m, raw)
    })
    .map_err(|e| e.with_sorted_input_loc())?;
    c.finish_chrom().map_err(|e| e.with_sorted_input_loc())?;
    log_skipped(&skipped.take());
    if o.per_strand && unknown_strand {
        log::warn!(
            "some records have no strand; with -s they are clustered together, apart from + and - \
records"
        );
    }
    Ok(())
}
