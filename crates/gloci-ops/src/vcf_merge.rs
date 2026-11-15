//! `gloci vcf merge`: a k-way streaming merge of VCF/BCF files sorted in one
//! chromosome order. Each input is read a batch at a time and checked for
//! order; at each position the records of every input at it are paired into
//! rows (`gloci_vcf::merge::rows`) and each row is merged into one record.
//! Nothing is indexed, and only the records of one position are held.

use crate::vcf_split::vcf_rewrite_output;
use gloci_core::{ChromId, ChromTable};
use gloci_io::{
    BcfText, Header, IoError, OrderGuard, OrderOptions, OutFormat, RecordBatch, RecordReader,
    RecordWriter, TextToBcf, VcfHeader, next_batch_size,
};
use gloci_vcf::merge::{Merger, PairIn, stage_position};
use std::cmp::Ordering;

pub const COMMAND: &str = "vcf merge";

/// Counters of one run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MergeSummary {
    pub records_in: u64,
    pub records_out: u64,
    /// Records written that combine records of more than one input.
    pub records_merged: u64,
    /// Records skipped by `--apply-filters`.
    pub records_filtered: u64,
    /// The most records held at once.
    pub peak_buffered: u64,
}

/// The output format (`requested`, else the first input's) for the merged header text.
pub fn vcf_merge_output(
    first: &dyn RecordReader,
    requested: Option<OutFormat>,
    header_text: &[u8],
) -> Result<(OutFormat, Header), IoError> {
    let (to, _) = vcf_rewrite_output(COMMAND, first, requested)?;
    let header = VcfHeader::parse(header_text.to_vec()).map_err(|e| IoError::Detect {
        path: first.name().to_string(),
        msg: e,
    })?;
    Ok((to, Header::Vcf(header)))
}

/// `--apply-filters LIST`: keep a record whose FILTER has any listed name.
#[derive(Clone, Debug)]
pub struct FilterList {
    names: Vec<Vec<u8>>,
}

impl FilterList {
    pub fn parse(list: &str) -> FilterList {
        FilterList {
            names: list
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.as_bytes().to_vec())
                .collect(),
        }
    }

    fn keeps(&self, line: &[u8]) -> bool {
        let filter = line.split(|&b| b == b'\t').nth(6).unwrap_or(b".");
        if filter == b"." {
            return self.names.iter().any(|n| n == b".");
        }
        filter
            .split(|&b| b == b';')
            .any(|f| self.names.iter().any(|n| n == f))
    }
}

struct Stream {
    reader: Box<dyn RecordReader>,
    name: String,
    batch: RecordBatch,
    i: usize,
    more: bool,
    max: usize,
    guard: OrderGuard,
    bcf: Option<BcfText>,
    filter: Option<FilterList>,
    skipped: u64,
}

impl Stream {
    fn new(reader: Box<dyn RecordReader>, filter: Option<FilterList>) -> Result<Stream, IoError> {
        let name = reader.name().to_string();
        // vcf merge has --chrom-order but no -g.
        let guard = OrderGuard::new(&name, reader.format()).with_options(OrderOptions::CHROM_ORDER);
        let bcf = BcfText::for_reader(reader.as_ref())?;
        Ok(Stream {
            reader,
            name,
            batch: RecordBatch::default(),
            i: 0,
            more: true,
            max: 1,
            guard,
            bcf,
            filter,
            skipped: 0,
        })
    }

    /// The `(chrom, key)` of the next record, reading ahead when needed and
    /// skipping the records `--apply-filters` drops.
    fn peek(&mut self, chroms: &mut ChromTable) -> Result<Option<(ChromId, u64)>, IoError> {
        loop {
            while self.i >= self.batch.len() {
                if !self.more {
                    return Ok(None);
                }
                self.batch.clear();
                self.i = 0;
                self.more = self.reader.read_batch(chroms, &mut self.batch, self.max)?;
                self.guard.check(chroms, &self.batch)?;
                self.max = next_batch_size(&self.batch, self.max);
            }
            if let Some(f) = &self.filter {
                let loc = self.batch.recs[self.i].loc;
                let raw = self.batch.raw(self.i);
                let keep = match self.bcf.as_mut() {
                    Some(t) => f.keeps(t.render(raw, loc)?),
                    None => f.keeps(raw),
                };
                if !keep {
                    self.i += 1;
                    self.skipped += 1;
                    continue;
                }
            }
            let m = &self.batch.recs[self.i];
            return Ok(Some((m.chrom, m.key)));
        }
    }

    /// Take the next record as a VCF line.
    fn take(&mut self) -> Result<Vec<u8>, IoError> {
        let loc = self.batch.recs[self.i].loc;
        let raw = self.batch.raw(self.i);
        let line = match self.bcf.as_mut() {
            Some(t) => t.render(raw, loc)?.to_vec(),
            None => raw.to_vec(),
        };
        self.i += 1;
        Ok(line)
    }
}

fn data_err(msg: String) -> IoError {
    IoError::Detect {
        path: COMMAND.to_string(),
        msg,
    }
}

/// Merge `inputs` into `w` (opened with [`vcf_merge_output`]'s format and header).
pub fn vcf_merge(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    merger: &mut Merger,
    filter: Option<FilterList>,
    out_header: &VcfHeader,
    w: &mut dyn RecordWriter,
) -> Result<MergeSummary, IoError> {
    let mut streams = inputs
        .into_iter()
        .map(|r| Stream::new(r, filter.clone()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut text_to_bcf = if w.is_bcf() {
        Some(TextToBcf::new(COMMAND, out_header)?)
    } else {
        None
    };
    let write_err = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    let n = streams.len();
    let mut sum = MergeSummary::default();
    let (mut line, mut enc) = (Vec::new(), Vec::new());
    let mut heads: Vec<Option<(ChromId, u64)>> = vec![None; n];
    loop {
        let mut min: Option<(ChromId, u64)> = None;
        for (i, s) in streams.iter_mut().enumerate() {
            heads[i] = s.peek(chroms)?;
            if let Some(h) = heads[i]
                && min.is_none_or(|m| chroms.cmp(h.0, m.0).then(h.1.cmp(&m.1)) == Ordering::Less)
            {
                min = Some(h);
            }
        }
        let Some(min) = min else { break };
        // Every record at the position, per input.
        let mut lines: Vec<Vec<Vec<u8>>> = vec![Vec::new(); n];
        let mut held = 0u64;
        for (i, s) in streams.iter_mut().enumerate() {
            while s.peek(chroms)? == Some(min) {
                lines[i].push(s.take()?);
                sum.records_in += 1;
                held += 1;
            }
            let _ = &s.name;
        }
        sum.peak_buffered = sum.peak_buffered.max(held);
        let mut pairs: Vec<Vec<PairIn>> = Vec::with_capacity(n);
        for (i, ls) in lines.iter().enumerate() {
            let mut v = Vec::with_capacity(ls.len());
            for l in ls {
                v.push(
                    merger
                        .pair_in(l)
                        .map_err(|e| data_err(format!("{}: {e}", streams[i].name)))?,
                );
            }
            pairs.push(v);
        }
        let staged = stage_position(&pairs, merger.options().mode).map_err(|e| {
            data_err(format!(
                "{e} at {}:{}",
                String::from_utf8_lossy(chroms.name(min.0)),
                min.1 + 1
            ))
        })?;
        for row in staged {
            let members: Vec<Option<&[u8]>> = row
                .members
                .iter()
                .enumerate()
                .map(|(i, r)| r.map(|k| lines[i][k].as_slice()))
                .collect();
            line.clear();
            merger
                .merge_row(&members, &row.alleles, &mut line)
                .map_err(data_err)?;
            if members.iter().flatten().count() > 1 {
                sum.records_merged += 1;
            }
            sum.records_out += 1;
            match text_to_bcf.as_mut() {
                Some(t) => {
                    enc.clear();
                    t.encode(&line, 0, &mut enc)?;
                    w.write_record(&enc).map_err(write_err)?;
                }
                None => w.write_record(&line).map_err(write_err)?,
            }
        }
    }
    sum.records_filtered = streams.iter().map(|s| s.skipped).sum();
    for s in &streams {
        if let Some(m) = s.reader.diagnostics().summary(&s.name) {
            log::warn!("{m}");
        }
    }
    Ok(sum)
}
