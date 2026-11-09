//! `partition`: split the footprints of sorted inputs at every boundary into
//! disjoint segments and report which inputs cover each (replaces `bedops -p`,
//! `-i` and `-s`, `bedtools multiinter` and `unionbedg`).
//!
//! The inputs are k-way merged ([`kmerge_with`], which also checks their
//! order) and swept once. The sweep holds only the ends of the intervals that
//! are open at the current position, so memory follows the nesting depth, not
//! the file size. Points cover no bases: they add no segment and no boundary.

use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable, Genome};
use gloci_engine::kmerge_with;
use gloci_io::{Format, IoError, Loc, RecordMeta, RecordReader, RecordWriter, UsageError};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Default)]
pub struct PartitionOptions {
    /// One label per input, written in the membership column.
    pub names: Vec<String>,
    /// Keep segments covered by at least this many inputs (default 1, or 0
    /// when `empty` is set, so that the uncovered segments are kept).
    pub min_files: Option<usize>,
    /// Keep segments covered by at most this many inputs.
    pub max_files: Option<usize>,
    /// Keep segments covered by every input.
    pub all: bool,
    /// Join neighbouring kept segments with the same output columns.
    pub merge_adjacent: bool,
    /// Write the count and labels after the coordinates (ignored with
    /// `bedgraph`).
    pub membership: bool,
    /// Inputs are bedGraphs: write one value column per input instead of the
    /// membership.
    pub bedgraph: bool,
    /// With `bedgraph`: also write the segments no input covers, with this
    /// value in every column. With a genome, whole chromosomes are covered.
    pub empty: Option<String>,
}

/// Value of an input in a bedGraph segment it does not cover, when
/// `empty` is not given.
const ABSENT: &str = "0";

/// Reject option and input combinations that cannot work, before any input
/// is read.
pub fn validate(formats: &[Format], o: &PartitionOptions) -> Result<(), UsageError> {
    if o.names.len() != formats.len() {
        return Err(UsageError(format!(
            "--names gives {} labels for {} input files; give one label per file",
            o.names.len(),
            formats.len()
        )));
    }
    for (i, n) in o.names.iter().enumerate() {
        if n.is_empty() || n.contains([',', '\t', '\n']) {
            return Err(UsageError(format!(
                "label '{n}' cannot be used: labels are joined with commas, so they must be \
non-empty and contain no comma or tab; rename it with --names"
            )));
        }
        if o.names[..i].contains(n) {
            return Err(UsageError(format!(
                "two inputs are labelled '{n}', which would make the membership ambiguous; \
give distinct labels with --names"
            )));
        }
    }
    if o.bedgraph
        && let Some(f) = formats
            .iter()
            .find(|f| !matches!(f, Format::Bed | Format::BedGraph | Format::Wig))
    {
        return Err(UsageError(format!(
            "--bedgraph needs bedGraph inputs (chrom, start, end, value) but one is {}; \
drop --bedgraph for a membership table",
            f.name()
        )));
    }
    if o.empty.is_some() && !o.bedgraph {
        return Err(UsageError(
            "--empty needs --bedgraph; the membership table has no value columns to fill".into(),
        ));
    }
    if o.min_files == Some(0) && o.empty.is_none() {
        return Err(UsageError(
            "--min-files 0 would keep segments no file covers, which are written only with \
--empty; add --empty VALUE (with --bedgraph), or use --min-files 1 or drop it"
                .into(),
        ));
    }
    if let (Some(lo), Some(hi)) = (o.min_files, o.max_files)
        && lo > hi
    {
        return Err(UsageError(format!(
            "--min-files {lo} is above --max-files {hi}, so no segment can pass; \
swap or change them"
        )));
    }
    Ok(())
}

struct Pending {
    chrom: ChromId,
    name: Vec<u8>,
    start: u64,
    end: u64,
    payload: Vec<u8>,
}

struct Sweep<'a> {
    o: &'a PartitionOptions,
    names: Vec<String>,
    genome: Option<&'a Genome>,
    w: &'a mut dyn RecordWriter,
    lo: usize,
    hi: usize,
    chrom: Option<ChromId>,
    chrom_name: Vec<u8>,
    /// Ends of the open intervals.
    ends: BinaryHeap<Reverse<(u64, usize)>>,
    /// Open intervals per input.
    cover: Vec<u32>,
    /// Inputs with at least one open interval.
    active: usize,
    /// bedGraph value of each input's open interval.
    values: Vec<Vec<u8>>,
    /// Start of the segment being swept; `None` until the first feature of
    /// a chromosome when nothing says where the chromosome starts.
    cursor: Option<u64>,
    pending: Option<Pending>,
    scratch: Vec<u8>,
    /// Output row buffer, reused.
    row: Vec<u8>,
    /// With a genome: the next genome chromosome not yet visited.
    next_genome: usize,
    peak_open: usize,
}

impl Sweep<'_> {
    fn fill_genome(&self) -> bool {
        self.o.empty.is_some() && self.genome.is_some()
    }

    fn write_row(
        &mut self,
        name: &[u8],
        start: u64,
        end: u64,
        payload: &[u8],
    ) -> Result<(), IoError> {
        let mut row = std::mem::take(&mut self.row);
        row.clear();
        row.extend_from_slice(name);
        row.push(b'\t');
        row.extend_from_slice(start.to_string().as_bytes());
        row.push(b'\t');
        row.extend_from_slice(end.to_string().as_bytes());
        if self.o.bedgraph || self.o.membership {
            row.push(b'\t');
            row.extend_from_slice(payload);
        }
        let r = self.w.write_record(&row).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        });
        self.row = row;
        r
    }

    fn flush_pending(&mut self) -> Result<(), IoError> {
        if let Some(p) = self.pending.take() {
            self.write_row(&p.name, p.start, p.end, &p.payload)?;
        }
        Ok(())
    }

    /// Offer the segment `[start, end)` whose columns are in `self.scratch`.
    fn put(&mut self, start: u64, end: u64) -> Result<(), IoError> {
        let chrom = self.chrom.expect("a chromosome is open");
        if self.o.merge_adjacent
            && let Some(p) = self.pending.as_mut()
            && p.chrom == chrom
            && p.end == start
            && p.payload == self.scratch
        {
            p.end = end;
            return Ok(());
        }
        self.flush_pending()?;
        if self.o.merge_adjacent {
            self.pending = Some(Pending {
                chrom,
                name: self.chrom_name.clone(),
                start,
                end,
                payload: self.scratch.clone(),
            });
            Ok(())
        } else {
            let name = std::mem::take(&mut self.chrom_name);
            let payload = std::mem::take(&mut self.scratch);
            let r = self.write_row(&name, start, end, &payload);
            self.scratch = payload;
            self.chrom_name = name;
            r
        }
    }

    /// Sweep the segment from the cursor to `end` with the current coverage.
    fn segment(&mut self, end: u64) -> Result<(), IoError> {
        let Some(start) = self.cursor else {
            return Ok(());
        };
        if end <= start {
            return Ok(());
        }
        self.cursor = Some(end);
        let n = self.active;
        if n == 0 && self.o.empty.is_none() {
            return Ok(());
        }
        if n < self.lo || n > self.hi {
            return Ok(());
        }
        self.scratch.clear();
        if self.o.bedgraph {
            let absent = self.o.empty.as_deref().unwrap_or(ABSENT).as_bytes();
            for i in 0..self.cover.len() {
                if i > 0 {
                    self.scratch.push(b'\t');
                }
                if self.cover[i] > 0 {
                    self.scratch.extend_from_slice(&self.values[i]);
                } else {
                    self.scratch.extend_from_slice(absent);
                }
            }
        } else if self.o.membership {
            self.scratch.extend_from_slice(n.to_string().as_bytes());
            self.scratch.push(b'\t');
            let mut first = true;
            for (i, c) in self.cover.iter().enumerate() {
                if *c > 0 {
                    if !first {
                        self.scratch.push(b',');
                    }
                    first = false;
                    self.scratch.extend_from_slice(self.names[i].as_bytes());
                }
            }
        }
        self.put(start, end)
    }

    fn close_end(&mut self, end: u64, input: usize) -> Result<(), IoError> {
        self.segment(end)?;
        self.cover[input] -= 1;
        if self.cover[input] == 0 {
            self.active -= 1;
        }
        Ok(())
    }

    /// Finish every segment that ends at or before `pos`, then the stretch
    /// up to `pos`.
    fn advance_to(&mut self, pos: u64) -> Result<(), IoError> {
        while let Some(&Reverse((e, i))) = self.ends.peek() {
            if e > pos {
                break;
            }
            self.ends.pop();
            self.close_end(e, i)?;
        }
        self.segment(pos)
    }

    fn gap_row(&mut self, chrom: ChromId, name: &[u8], len: u64) -> Result<(), IoError> {
        let saved = (self.chrom, std::mem::take(&mut self.chrom_name));
        self.chrom = Some(chrom);
        self.chrom_name = name.to_vec();
        self.cursor = Some(0);
        let r = self.segment(len);
        self.cursor = None;
        self.chrom = saved.0;
        self.chrom_name = saved.1;
        r
    }

    /// With a genome: whole-chromosome rows for the chromosomes before `upto`.
    fn untouched_before(&mut self, upto: usize) -> Result<(), IoError> {
        let Some(g) = self.genome else {
            return Ok(());
        };
        if !self.fill_genome() {
            return Ok(());
        }
        while self.next_genome < upto {
            let i = self.next_genome;
            self.next_genome += 1;
            let len = g.length(i).expect("lengths checked");
            // The previous chromosome's pending row must go out first.
            self.flush_pending()?;
            self.gap_row(ChromId(i as u32), g.name(i), len)?;
        }
        self.flush_pending()
    }

    fn record(
        &mut self,
        input: usize,
        m: &RecordMeta,
        raw: &[u8],
        chroms: &ChromTable,
        inputs: &[(String, bool)],
    ) -> Result<(), IoError> {
        if m.span.start == m.span.end {
            return Ok(());
        }
        if self.chrom != Some(m.chrom) {
            self.finish_chrom()?;
            self.chrom = Some(m.chrom);
            self.chrom_name = chroms.name(m.chrom).to_vec();
            if self.genome.is_some() {
                // A closed genome table: the id is the genome index.
                self.untouched_before(m.chrom.0 as usize)?;
                self.next_genome = m.chrom.0 as usize + 1;
            }
            self.cursor = self.fill_genome().then_some(0);
        }
        let (start, end) = (m.span.start, m.span.end);
        if self.fill_genome() {
            let g = self.genome.expect("genome");
            let len = g.length(m.chrom.0 as usize).expect("lengths checked");
            if end > len {
                return Err(IoError::Parse {
                    path: inputs[input].0.clone(),
                    loc: Loc::Line(m.loc),
                    msg: format!(
                        "{}:{start}-{end} extends past the end of {} (length {len}); \
check that -g matches the genome build of the input",
                        String::from_utf8_lossy(&self.chrom_name),
                        String::from_utf8_lossy(&self.chrom_name),
                    ),
                });
            }
        }
        self.advance_to(start)?;
        if self.cursor.is_none() {
            self.cursor = Some(start);
        }
        if self.o.bedgraph {
            let at = |msg: String| IoError::Parse {
                path: inputs[input].0.clone(),
                loc: Loc::Line(m.loc),
                msg,
            };
            if self.cover[input] > 0 {
                return Err(at(format!(
                    "{}:{start}-{end} overlaps an earlier interval of the same file; \
--bedgraph needs non-overlapping intervals in each file; resolve the overlap, \
for example by summing with `gloci partition --bedgraph` on files split so that \
none overlaps, or drop --bedgraph",
                    String::from_utf8_lossy(&self.chrom_name)
                )));
            }
            let Some(v) = raw.split(|&b| b == b'\t').nth(3) else {
                return Err(at(
                    "no value in column 4; --bedgraph needs bedGraph input (chrom, start, end, \
value)"
                        .into(),
                ));
            };
            self.values[input].clear();
            self.values[input].extend_from_slice(v);
        }
        self.cover[input] += 1;
        if self.cover[input] == 1 {
            self.active += 1;
        }
        self.ends.push(Reverse((end, input)));
        self.peak_open = self.peak_open.max(self.ends.len());
        Ok(())
    }

    fn finish_chrom(&mut self) -> Result<(), IoError> {
        let Some(chrom) = self.chrom else {
            return Ok(());
        };
        while let Some(Reverse((e, i))) = self.ends.pop() {
            self.close_end(e, i)?;
        }
        if self.fill_genome() {
            let g = self.genome.expect("genome");
            let len = g.length(chrom.0 as usize).expect("lengths checked");
            self.segment(len)?;
        }
        self.flush_pending()?;
        self.chrom = None;
        self.cursor = None;
        Ok(())
    }
}

/// Partition `inputs` (each sorted in `chroms` order) and write one row per
/// segment: `chrom start end`, then the count and labels of the covering
/// inputs (with `o.membership`) or one value per input (with `o.bedgraph`).
/// `genome`, when given, must be the genome `chroms` was built from (a closed
/// table, so chromosome ids are genome indexes); with `o.empty` it also
/// supplies the chromosome lengths.
pub fn partition(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    genome: Option<&Genome>,
    o: &PartitionOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    partition_counting(inputs, chroms, genome, o, w)
        .map(|_| ())
        .map_err(|e| match e.with_sorted_input_loc() {
            IoError::Parse { path, loc, msg } if path.starts_with(gloci_io::SORTED_INPUTS) => {
                IoError::Parse {
                    path,
                    loc,
                    msg: format!(
                        "{msg} (the record number counts records after --sort, not source lines)"
                    ),
                }
            }
            e => e,
        })
}

/// [`partition`], returning the most intervals ever open at once (for the
/// memory-flatness test).
#[doc(hidden)]
pub fn partition_counting(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    genome: Option<&Genome>,
    o: &PartitionOptions,
    w: &mut dyn RecordWriter,
) -> Result<usize, IoError> {
    let formats: Vec<Format> = inputs.iter().map(|r| r.format()).collect();
    validate(&formats, o)?;
    if o.empty.is_some()
        && let Some(g) = genome
        && let Some(i) = (0..g.len()).find(|&i| g.length(i).is_none())
    {
        return Err(UsageError(format!(
            "-g: the genome has no length for {}; --empty with -g needs chromosome lengths \
(use a built-in genome, a .genome or .fai file, or a .dict file)",
            String::from_utf8_lossy(g.name(i))
        ))
        .into());
    }
    let n = inputs.len();
    let info: Vec<(String, bool)> = inputs
        .iter()
        .map(|r| (r.name().to_string(), r.format() == Format::Bcf))
        .collect();
    let mut skipped = Skipped::default();
    let inputs: Vec<Box<dyn RecordReader>> = inputs.into_iter().map(|r| skipped.watch(r)).collect();
    let mut chroms_seen: Vec<Vec<Vec<u8>>> = vec![Vec::new(); n];
    let mut last_chrom: Vec<Option<ChromId>> = vec![None; n];
    let (lo, hi) = if o.all {
        (n, n)
    } else {
        (
            o.min_files.unwrap_or(usize::from(o.empty.is_none())),
            o.max_files.unwrap_or(usize::MAX),
        )
    };
    let mut s = Sweep {
        o,
        names: o.names.clone(),
        genome,
        w,
        lo,
        hi,
        chrom: None,
        chrom_name: Vec::new(),
        ends: BinaryHeap::new(),
        cover: vec![0; n],
        active: 0,
        values: vec![Vec::new(); n],
        cursor: None,
        pending: None,
        scratch: Vec::new(),
        row: Vec::new(),
        next_genome: 0,
        peak_open: 0,
    };
    kmerge_with(inputs, chroms, &mut |input, meta, raw, chroms| {
        if last_chrom[input] != Some(meta.chrom) {
            last_chrom[input] = Some(meta.chrom);
            let c_name = chroms.name(meta.chrom);
            if !chroms_seen[input].iter().any(|x| x.as_slice() == c_name) {
                chroms_seen[input].push(c_name.to_vec());
            }
        }
        s.record(input, meta, raw, chroms, &info)
    })?;
    s.finish_chrom()?;
    if let Some(g) = genome {
        s.untouched_before(g.len())?;
    }
    log_skipped(&skipped.take());
    let names: Vec<String> = info.iter().map(|i| i.0.clone()).collect();
    if let Some(w) = gloci_io::naming::nary_naming_warning(&names, &chroms_seen) {
        log::warn!("{w}");
    }
    Ok(s.peak_open)
}
