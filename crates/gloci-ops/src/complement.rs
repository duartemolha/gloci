//! `complement`: the bases no feature covers (replaces `bedtools complement`
//! and `bedops -c`).
//!
//! With a genome, every genome chromosome is reported in genome order, with
//! its uncovered gaps over `[0, length)`; chromosomes without features give
//! one full-length row. Without one, only the gaps between features of the
//! same chromosome. Points cover no bases and are ignored (D48).

use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable, Genome};
use gloci_engine::kmerge_with;
use gloci_io::{IoError, Loc, RecordMeta, RecordReader, RecordWriter, UsageError};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// A footprint waiting for its turn in start order (a span never starts
/// before its record's key, but the key is what the inputs are sorted by).
#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct Pending {
    start: u64,
    end: u64,
    input: usize,
    loc: u64,
}

struct Comp<'a> {
    genome: Option<&'a Genome>,
    names: Vec<String>,
    binary: Vec<bool>,
    w: &'a mut dyn RecordWriter,
    heap: BinaryHeap<Reverse<Pending>>,
    chrom: Option<ChromId>,
    /// End of the covered prefix of the current chromosome; `None` until its
    /// first feature when there is no genome (no leading gap).
    cursor: Option<u64>,
    /// With a genome: the next genome chromosome not yet written.
    next_genome: usize,
}

impl Comp<'_> {
    fn row(&mut self, name: &[u8], start: u64, end: u64) -> Result<(), IoError> {
        let mut row = name.to_vec();
        row.push(b'\t');
        row.extend_from_slice(start.to_string().as_bytes());
        row.push(b'\t');
        row.extend_from_slice(end.to_string().as_bytes());
        self.w.write_record(&row).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })
    }

    /// Write full-length rows for genome chromosomes before index `upto`.
    fn untouched_before(&mut self, upto: usize) -> Result<(), IoError> {
        let Some(g) = self.genome else {
            return Ok(());
        };
        while self.next_genome < upto {
            let i = self.next_genome;
            self.next_genome += 1;
            let len = g.length(i).expect("lengths checked");
            if len > 0 {
                self.row(g.name(i), 0, len)?;
            }
        }
        Ok(())
    }

    fn record(&mut self, input: usize, m: &RecordMeta, chroms: &ChromTable) -> Result<(), IoError> {
        if m.span.start == m.span.end {
            return Ok(());
        }
        if self.chrom != Some(m.chrom) {
            self.finish_chrom(chroms)?;
            self.chrom = Some(m.chrom);
            if self.genome.is_some() {
                // A closed genome table: the id is the genome index.
                self.untouched_before(m.chrom.0 as usize)?;
                self.cursor = Some(0);
            } else {
                self.cursor = None;
            }
        }
        while self.heap.peek().is_some_and(|p| p.0.start <= m.key) {
            let Reverse(p) = self.heap.pop().expect("peeked");
            self.add(p, chroms)?;
        }
        self.heap.push(Reverse(Pending {
            start: m.span.start,
            end: m.span.end,
            input,
            loc: m.loc,
        }));
        Ok(())
    }

    fn add(&mut self, p: Pending, chroms: &ChromTable) -> Result<(), IoError> {
        let chrom = self.chrom.expect("a chromosome is open");
        if let Some(g) = self.genome {
            let len = g.length(chrom.0 as usize).expect("lengths checked");
            if p.end > len {
                return Err(IoError::Parse {
                    path: self.names[p.input].clone(),
                    loc: if self.binary[p.input] {
                        Loc::Record(p.loc)
                    } else {
                        Loc::Line(p.loc)
                    },
                    msg: format!(
                        "{}:{}-{} extends past the end of {} (length {len}); \
check that -g matches the genome build of the input",
                        String::from_utf8_lossy(chroms.name(chrom)),
                        p.start,
                        p.end,
                        String::from_utf8_lossy(chroms.name(chrom)),
                    ),
                });
            }
        }
        if let Some(c) = self.cursor
            && p.start > c
        {
            self.row(chroms.name(chrom), c, p.start)?;
        }
        self.cursor = Some(self.cursor.map_or(p.end, |c| c.max(p.end)));
        Ok(())
    }

    fn finish_chrom(&mut self, chroms: &ChromTable) -> Result<(), IoError> {
        let Some(chrom) = self.chrom else {
            return Ok(());
        };
        while let Some(Reverse(p)) = self.heap.pop() {
            self.add(p, chroms)?;
        }
        if let Some(g) = self.genome {
            let i = chrom.0 as usize;
            let len = g.length(i).expect("lengths checked");
            let c = self.cursor.unwrap_or(0);
            if c < len {
                self.row(chroms.name(chrom), c, len)?;
            }
            self.next_genome = i + 1;
        }
        self.chrom = None;
        Ok(())
    }
}

/// Write the complement of the union of `inputs` (each sorted in `chroms`
/// order) as BED3 rows. `genome` must be the genome `chroms` was built from
/// (a closed table, so chromosome ids are genome indexes).
pub fn complement(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    genome: Option<&Genome>,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    complement_impl(inputs, chroms, genome, w).map_err(|e| e.with_sorted_input_loc())
}

fn complement_impl(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    genome: Option<&Genome>,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    if let Some(g) = genome
        && let Some(i) = (0..g.len()).find(|&i| g.length(i).is_none())
    {
        return Err(UsageError(format!(
            "-g: the genome has no length for {}; complement needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
            String::from_utf8_lossy(g.name(i))
        ))
        .into());
    }
    let names = inputs
        .iter()
        .map(|r| r.name().to_string())
        .collect::<Vec<_>>();
    let mut skipped = Skipped::default();
    let inputs: Vec<Box<dyn RecordReader>> = inputs.into_iter().map(|r| skipped.watch(r)).collect();
    let mut chroms_seen: Vec<Vec<Vec<u8>>> = vec![Vec::new(); inputs.len()];
    let mut last_chrom: Vec<Option<ChromId>> = vec![None; inputs.len()];
    let mut marks = vec![gloci_io::naming::SeenChroms::default(); inputs.len()];
    let mut c = Comp {
        genome,
        names: names.clone(),
        binary: inputs.iter().map(|r| r.format().is_binary()).collect(),
        w,
        heap: BinaryHeap::new(),
        chrom: None,
        cursor: None,
        next_genome: 0,
    };
    kmerge_with(inputs, chroms, &mut |input, meta, _raw, chroms| {
        if last_chrom[input] != Some(meta.chrom) {
            last_chrom[input] = Some(meta.chrom);
            if marks[input].first_time(meta.chrom) {
                chroms_seen[input].push(chroms.name(meta.chrom).to_vec());
            }
        }
        c.record(input, meta, chroms)
    })?;
    c.finish_chrom(chroms)?;
    if let Some(g) = genome {
        c.untouched_before(g.len())?;
    }
    log_skipped(&skipped.take());
    if let Some(w) = gloci_io::naming::nary_naming_warning(&names, &chroms_seen) {
        log::warn!("{w}");
    }
    Ok(())
}
