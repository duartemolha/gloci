//! `gloci liftover`: move records to another assembly through a chain file, by CrossMap
//! 0.7.4's rules (D174).

mod bed;
mod gff;
mod names;

pub use names::{ChromStyle, source_name, target_name};

use gloci_core::{ChromTable, Interval};
use gloci_io::{
    ChainIndex, Format, IoError, Loc, RecordBatch, RecordReader, RecordWriter, Segment, SourceRef,
    UsageError, next_batch_size,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Default)]
pub struct LiftoverOptions {
    pub chrom_style: ChromStyle,
    /// Skip records the op finds invalid and count them (D174, R6).
    pub skip_invalid: bool,
    /// Lift BED12 records by columns 1 to 3 only, like any other BED record (CrossMap
    /// `--naive-bed-parsing`; D174, R8).
    pub ignore_blocks: bool,
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
    /// A predicate that reads after a count: "7 on no chain".
    pub fn describe(self) -> &'static str {
        match self {
            Unmapped::OffChain => "on no chain",
            Unmapped::BlockSplit => "with a block not in exactly one piece",
            Unmapped::BlocksApart => "with blocks on different chromosomes or strands",
            Unmapped::BadBlocks => "with blocks or a thick range that do not fit after lifting",
            Unmapped::SeveralPieces => "in several pieces",
            Unmapped::Partial => "only partly on a chain (also written, shortened)",
            Unmapped::Strand => "with a strand other than +, - or .",
            Unmapped::LowRatio => "with a mapped fraction below --min-ratio",
            Unmapped::CrossChroms => "with pieces on different chromosomes",
            Unmapped::MixedStrands => "with pieces on different strands",
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
    /// The first skipped invalid record: its line and why.
    pub first_invalid: Option<(u64, String)>,
    pub unmapped: BTreeMap<Unmapped, u64>,
}

impl LiftoverCounts {
    pub fn unmapped_total(&self) -> u64 {
        self.unmapped.values().sum()
    }
}

/// The bases of the record that some piece covers (pieces are sorted by source start).
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

/// An input chromosome as the record loop needs it.
struct Input {
    src: Option<SourceRef>,
    chr: bool,
}

/// The styled target name per (chain, input starts with `chr`), built once.
struct Styler {
    style: ChromStyle,
    names: Vec<Option<String>>,
}

impl Styler {
    fn new(style: ChromStyle, chain: &ChainIndex) -> Styler {
        Styler {
            style,
            names: vec![None; chain.chains() * 2],
        }
    }

    fn name(&mut self, chain: &ChainIndex, id: u32, chr: bool) -> &str {
        let style = self.style;
        self.names[id as usize * 2 + usize::from(chr)].get_or_insert_with(|| {
            target_name(style, if chr { "chr" } else { "" }, &chain.chain(id).target)
        })
    }

    fn same_name(&mut self, chain: &ChainIndex, a: u32, b: u32, chr: bool) -> bool {
        self.name(chain, a, chr);
        self.name(chain, b, chr);
        let at = |id: u32| &self.names[id as usize * 2 + usize::from(chr)];
        at(a) == at(b)
    }
}

/// Count a record that did not lift and, when asked, write it as read.
fn leave_unmapped(
    counts: &mut LiftoverCounts,
    unmapped: &mut Option<&mut dyn RecordWriter>,
    why: Unmapped,
    raw: &[u8],
) -> Result<(), IoError> {
    *counts.unmapped.entry(why).or_default() += 1;
    if let Some(u) = unmapped.as_deref_mut() {
        u.write_record(raw)
            .map_err(|e| write_err("unmapped output", e))?;
    }
    Ok(())
}

fn write_err(path: &str, source: std::io::Error) -> IoError {
    IoError::Write {
        path: path.into(),
        source,
    }
}

fn refuse(format: Format) -> UsageError {
    UsageError(match format {
        Format::Vcf | Format::Bcf => "VCF liftover is not supported yet (CrossMap vcf)".to_string(),
        Format::Wig => {
            "WIG liftover is not supported; convert it first with gloci convert x.wig -o x.bedgraph"
                .to_string()
        }
        Format::Bedpe => "BEDPE liftover is not supported".to_string(),
        f => format!(
            "liftover works on BED, bedGraph, GFF and GTF input, not {}",
            f.name()
        ),
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
    if !matches!(
        format,
        Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf
    ) {
        return Err(refuse(format).into());
    }
    let is_gff = matches!(format, Format::Gff3 | Format::Gtf);
    if is_gff && o.ignore_blocks {
        return Err(UsageError(
            "--ignore-blocks applies to BED12 records, not GFF or GTF features".into(),
        )
        .into());
    }
    let mut counts = LiftoverCounts::default();
    // Per input chromosome: its chain source, its name and whether it starts with `chr`.
    let mut sources: Vec<Option<Input>> = Vec::new();
    let mut styler = Styler::new(o.chrom_style, chain);
    let mut scratch = bed::BlockScratch::default();
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
                sources.resize_with(id + 1, || None);
            }
            let from = sources[id].get_or_insert_with(|| {
                let name = String::from_utf8_lossy(chroms.name(m.chrom)).into_owned();
                Input {
                    src: source_name(chain, &name),
                    chr: name.starts_with("chr"),
                }
            });
            if is_gff {
                match gff::lift(
                    raw,
                    m.span,
                    chain,
                    from.src.as_ref(),
                    &mut styler,
                    from.chr,
                    &mut pieces,
                    &mut buf,
                ) {
                    gff::Lift::Written { partial } => {
                        out.write_record(&buf).map_err(|e| write_err("output", e))?;
                        counts.written += 1;
                        if partial {
                            counts.shortened += 1;
                            leave_unmapped(&mut counts, &mut unmapped, Unmapped::Partial, raw)?;
                        }
                    }
                    gff::Lift::Unmapped(why) => {
                        leave_unmapped(&mut counts, &mut unmapped, why, raw)?;
                    }
                }
                continue;
            }
            if format == Format::Bed
                && !o.ignore_blocks
                && raw.split(|&b| b == b'\t').nth(11).is_some()
            {
                let lifted = bed::lift_blocks(
                    raw,
                    m.span,
                    chain,
                    from.src.as_ref(),
                    &mut styler,
                    from.chr,
                    &mut scratch,
                    &mut buf,
                );
                match lifted {
                    Ok(bed::Lift::Written { shortened }) => {
                        out.write_record(&buf).map_err(|e| write_err("output", e))?;
                        counts.written += 1;
                        counts.shortened += u64::from(shortened);
                    }
                    Ok(bed::Lift::Unmapped(why)) => {
                        leave_unmapped(&mut counts, &mut unmapped, why, raw)?;
                    }
                    Err(msg) if o.skip_invalid => {
                        counts.invalid += 1;
                        counts.first_invalid.get_or_insert((m.loc, msg));
                    }
                    Err(msg) => {
                        return Err(IoError::Parse {
                            path: input.name().to_string(),
                            loc: Loc::Line(m.loc),
                            msg: format!("{msg}; use --ignore-blocks to lift columns 1 to 3 only"),
                        });
                    }
                }
                continue;
            }
            pieces.clear();
            if let Some(src) = &from.src {
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
                leave_unmapped(&mut counts, &mut unmapped, Unmapped::OffChain, raw)?;
                continue;
            }
            for seg in &pieces {
                let target = styler.name(chain, seg.chain, from.chr);
                bed::rewrite_piece(raw, target, seg, &mut buf);
                out.write_record(&buf).map_err(|e| write_err("output", e))?;
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
