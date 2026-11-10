//! Which candidate hits count (D36, D37).

use super::{RelateInputs, RelateOptions};
use gloci_core::{
    Interval, OverlapCriteria, Strand, StrandRule, blocks_overlap, blocks_overlap_len, overlap_len,
};
use gloci_engine::{ARec, Hit};
use gloci_io::{AlnBlocks, Format, IoError, Loc, RecordReader, bed12_blocks};

/// Where an input's blocks come from.
enum BlockSource {
    /// BED12 columns, used under `--split`.
    Bed12 { name: String },
    /// CIGAR segments (D105), used under `--split` or `--split-deletions`.
    Alignment(AlnBlocks),
    /// One block, the span: VCF, BCF, GFF, GTF.
    Span,
}

impl BlockSource {
    fn of(r: &dyn RecordReader, split_deletions: bool) -> BlockSource {
        match r.format() {
            Format::Bed | Format::BedGraph | Format::Wig => BlockSource::Bed12 {
                name: r.name().to_string(),
            },
            f if f.is_alignment() => AlnBlocks::for_reader(r, split_deletions)
                .map_or(BlockSource::Span, BlockSource::Alignment),
            _ => BlockSource::Span,
        }
    }

    /// Blocks of a record. Only BED-family records with 12 or more columns have blocks, and
    /// an alignment's are its CIGAR segments; every other record (VCF, BCF, GFF, BED below
    /// BED12) is one block, its span. `split` says whether BED12 blocks count: alignments
    /// always split, because the caller only asks when `--split` or `--split-deletions` is on.
    /// A BED12 line whose blocks do not parse is a data error naming the file, line and
    /// reason.
    fn blocks(
        &self,
        split: bool,
        loc: u64,
        raw: &[u8],
        span: Interval,
    ) -> Result<Vec<Interval>, IoError> {
        match self {
            BlockSource::Span => Ok(vec![span]),
            BlockSource::Alignment(b) => b.get(raw, loc, span),
            BlockSource::Bed12 { .. } if !split => Ok(vec![span]),
            BlockSource::Bed12 { name } => match bed12_blocks(raw) {
                Ok(Some(b)) => Ok(b),
                Ok(None) => Ok(vec![span]),
                Err(msg) => Err(IoError::Parse {
                    path: name.clone(),
                    loc: Loc::Line(loc),
                    msg: format!("invalid BED12 blocks: {msg}"),
                }),
            },
        }
    }
}

fn strand_ok(rule: StrandRule, a: Strand, b: Strand) -> bool {
    match rule {
        StrandRule::Any => true,
        rule => {
            a != Strand::Unknown && b != Strand::Unknown && ((a == b) == (rule == StrandRule::Same))
        }
    }
}

/// Decides which candidate hits count. Knows where each input's blocks come
/// from, so `--split` parses blocks only from BED12 records and alignments,
/// and keeps A's blocks across the hits of one A record.
pub struct Acceptor<'o> {
    o: &'o RelateOptions,
    /// A first, then each B.
    inputs: Vec<BlockSource>,
    a_blocks: Option<(u64, Vec<Interval>)>,
}

impl<'o> Acceptor<'o> {
    pub fn new(o: &'o RelateOptions, inputs: &RelateInputs) -> Acceptor<'o> {
        let info = |r: &dyn RecordReader| BlockSource::of(r, o.split_deletions);
        let mut v = vec![info(inputs.a.as_ref())];
        v.extend(inputs.b.iter().map(|r| info(r.as_ref())));
        Acceptor {
            o,
            inputs: v,
            a_blocks: None,
        }
    }

    /// Does `hit` pass the criteria for `a`?
    ///
    /// - **Window.** The engine already found the hit by its padded span,
    ///   and thresholds cannot be combined with a window (D37), so only the
    ///   strand rule is left to apply (with `--split`, the padded blocks
    ///   must meet).
    /// - **Otherwise** the thresholds are measured against B's own span.
    ///   With `--split`, A and B are their blocks: they must share a base
    ///   (or point), and thresholds use the total block lengths and the
    ///   total overlap, `blocks_overlap_len`.
    pub fn accepted(&mut self, a: &ARec, hit: &Hit) -> Result<bool, IoError> {
        let o = self.o;
        let c = &o.criteria;
        let (a_span, b_span) = (a.meta.span, hit.span);
        if !o.splits() {
            return Ok(match o.window {
                Some(_) => strand_ok(c.strand, a.meta.strand, hit.strand),
                None => c.accepts(a_span, a.meta.strand, b_span, hit.strand),
            });
        }
        self.load_a_blocks(a)?;
        let ab = &self.a_blocks.as_ref().expect("set above").1;
        let mut bb = self.hit_blocks(hit)?;
        if let Some(w) = o.window {
            for x in &mut bb {
                *x = w.apply(*x, a.meta.strand);
            }
            return Ok(blocks_overlap(ab, &bb, c.inclusive_points)
                && strand_ok(c.strand, a.meta.strand, hit.strand));
        }
        if !blocks_overlap(ab, &bb, c.inclusive_points) {
            return Ok(false);
        }
        let total = |v: &[Interval]| v.iter().map(Interval::len).sum::<u64>();
        let (alen, blen, ov) = (total(ab), total(&bb), blocks_overlap_len(ab, &bb));
        // `accepts` measures two intervals, so give it two with these sizes
        // and this overlap. Overlapping blocks that share no base are
        // points (D34), and those thresholds count as met.
        let (sa, sb, inclusive) = if alen == 0 || blen == 0 || ov == 0 {
            (Interval::point(0), Interval::point(0), true)
        } else {
            (
                Interval {
                    start: 0,
                    end: alen,
                },
                Interval {
                    start: alen - ov,
                    end: alen - ov + blen,
                },
                c.inclusive_points,
            )
        };
        let c = OverlapCriteria {
            inclusive_points: inclusive,
            ..*c
        };
        Ok(c.accepts(sa, a.meta.strand, sb, hit.strand))
    }

    /// Bases shared by `a` and `hit` (unpadded): the span overlap, or with
    /// `--split` the overlap of the blocks. Zero for points.
    pub fn overlap_len(&mut self, a: &ARec, hit: &Hit) -> Result<u64, IoError> {
        if !self.o.splits() {
            return Ok(overlap_len(a.meta.span, hit.span));
        }
        self.load_a_blocks(a)?;
        let ab = &self.a_blocks.as_ref().expect("set above").1;
        let bb = self.hit_blocks(hit)?;
        Ok(blocks_overlap_len(ab, &bb))
    }

    /// A's blocks under `--split` (one block, its span, for a record
    /// without BED12 blocks), kept for the hits of this A.
    pub fn a_blocks(&mut self, a: &ARec) -> Result<&[Interval], IoError> {
        self.load_a_blocks(a)?;
        Ok(&self.a_blocks.as_ref().expect("set above").1)
    }

    /// The blocks of `hit` (B's own, unpadded): its BED12 blocks under
    /// `--split` when it has them, its CIGAR segments when it is an alignment, else its
    /// span as one block.
    pub fn hit_blocks(&self, hit: &Hit) -> Result<Vec<Interval>, IoError> {
        let src = self
            .inputs
            .get(usize::from(hit.file) + 1)
            .ok_or_else(|| IoError::Internal(format!("hit from unknown B input {}", hit.file)))?;
        src.blocks(self.o.split, hit.loc, hit.raw, hit.span)
    }

    fn load_a_blocks(&mut self, a: &ARec) -> Result<(), IoError> {
        if self.a_blocks.as_ref().is_none_or(|(l, _)| *l != a.meta.loc) {
            let b = self.inputs[0].blocks(self.o.split, a.meta.loc, a.raw, a.meta.span)?;
            self.a_blocks = Some((a.meta.loc, b));
        }
        Ok(())
    }
}
