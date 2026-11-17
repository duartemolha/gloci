//! `gloci liftover`: move records to another assembly through a chain file, by CrossMap
//! 0.7.4's rules (D174).

mod bed;
mod names;

pub use names::{ChromStyle, source_name, target_name};

use gloci_core::{ChromTable, Interval};
use gloci_io::{
    ChainIndex, Format, IoError, RecordBatch, RecordReader, RecordWriter, Segment, SourceRef,
    UsageError, next_batch_size,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Default)]
pub struct LiftoverOptions {
    pub chrom_style: ChromStyle,
    /// Skip records the op finds invalid and count them (D174, R6).
    pub skip_invalid: bool,
}

/// Why a record is left unmapped (D174, R5 to R10).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Unmapped {
    OffChain,
    BlockSplit,
    BlocksApart,
    BadBlocks,
    SeveralPieces,
    Partial,
    Strand,
    LowRatio,
    CrossChroms,
    MixedStrands,
}

impl Unmapped {
    pub fn describe(self) -> &'static str {
        match self {
            Unmapped::OffChain => "on no chain",
            Unmapped::BlockSplit => "a block in several pieces",
            Unmapped::BlocksApart => "blocks on different chromosomes or strands",
            Unmapped::BadBlocks => "blocks that come out unsorted or overlapping",
            Unmapped::SeveralPieces => "a feature in several pieces",
            Unmapped::Partial => "a feature only partly on a chain",
            Unmapped::Strand => "a strand that cannot be mapped",
            Unmapped::LowRatio => "too little of it on the chain",
            Unmapped::CrossChroms => "pieces on different chromosomes",
            Unmapped::MixedStrands => "pieces on different strands",
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LiftoverCounts {
    pub records: u64,
    pub written: u64,
    /// Records written as 2 or more records.
    pub split: u64,
    /// Records whose pieces' source union is shorter than the record (never a zero-length one).
    pub shortened: u64,
    /// Records skipped by the op under `skip_invalid`.
    pub invalid: u64,
    pub unmapped: BTreeMap<Unmapped, u64>,
}

impl LiftoverCounts {
    pub fn unmapped_total(&self) -> u64 {
        self.unmapped.values().sum()
    }
}

/// The bases of `iv` that some piece covers (pieces are sorted by source start).
fn covered(pieces: &[Segment]) -> u64 {
    let (mut total, mut reach) = (0, 0);
    for p in pieces {
        let start = p.source.start.max(reach);
        if p.source.end > start {
            total += p.source.end - start;
        }
        reach = reach.max(p.source.end);
    }
    total
}

fn refuse(format: Format) -> UsageError {
    UsageError(match format {
        Format::Vcf | Format::Bcf => "VCF liftover is not supported yet (CrossMap vcf)".to_string(),
        Format::Wig => {
            "WIG liftover is not supported; convert it first with gloci convert x.wig -o x.bedgraph"
                .to_string()
        }
        Format::Bedpe => "BEDPE liftover is not supported".to_string(),
        Format::Gff3 | Format::Gtf => "GFF/GTF liftover is not supported yet".to_string(),
        f => format!("liftover works on BED and bedGraph input, not {f:?}"),
    })
}

/// Lift every record of `input` through `chain`, writing mapped records to `out` and, as read,
/// unmapped ones to `unmapped` (D174).
pub fn liftover(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    chain: &ChainIndex,
    o: &LiftoverOptions,
    out: &mut dyn RecordWriter,
    mut unmapped: Option<&mut dyn RecordWriter>,
) -> Result<LiftoverCounts, IoError> {
    let format = input.format();
    if !matches!(format, Format::Bed | Format::BedGraph) {
        return Err(refuse(format).into());
    }
    let write_err = |source| IoError::Write {
        path: "output".into(),
        source,
    };
    let mut counts = LiftoverCounts::default();
    let mut sources: Vec<Option<Option<SourceRef>>> = Vec::new();
    let mut pieces: Vec<Segment> = Vec::new();
    let mut buf = Vec::new();
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let raw = batch.raw(i);
            counts.records += 1;
            let id = m.chrom.0 as usize;
            if sources.len() <= id {
                sources.resize(id + 1, None);
            }
            let name = String::from_utf8_lossy(chroms.name(m.chrom));
            let src = sources[id].get_or_insert_with(|| source_name(chain, &name));
            pieces.clear();
            if let Some(src) = src {
                chain.segments(
                    src,
                    Interval {
                        start: m.span.start,
                        end: m.span.end,
                    },
                    &mut pieces,
                );
            }
            if pieces.is_empty() {
                *counts.unmapped.entry(Unmapped::OffChain).or_default() += 1;
                if let Some(u) = unmapped.as_deref_mut() {
                    u.write_record(raw).map_err(write_err)?;
                }
                continue;
            }
            for seg in &pieces {
                let target = target_name(o.chrom_style, &name, &chain.chain(seg.chain).target);
                bed::rewrite_piece(raw, &target, seg, &mut buf);
                out.write_record(&buf).map_err(write_err)?;
            }
            counts.written += pieces.len() as u64;
            if pieces.len() > 1 {
                counts.split += 1;
            }
            if m.span.len() > 0 && covered(&pieces) < m.span.len() {
                counts.shortened += 1;
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            crate::skipped::log_reader(input.as_ref());
            return Ok(counts);
        }
    }
}
