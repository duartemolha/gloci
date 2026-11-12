//! What a region query through a BAM index makes of one record (D107, D116): the span a
//! region is tested against, and the two ways records are kept.

use super::place::{IntervalRules, Verdict};
use super::{AlignOptions, Aln, AlnHeader, AlnRef};
use crate::diag::Diagnostics;
use crate::record::Parsed;
use gloci_core::{ChromId, ChromTable, Interval};
use std::sync::Arc;

const FLAG_UNMAPPED: u16 = 0x4;

/// The span of a record for a region query and for the index (D116): where the CIGAR puts
/// it, not cut at the end of its reference (that cut is for interval commands), or one base
/// at its position for an unmapped read that has one, as samtools and htslib place it.
/// `None` for a record with no place.
pub fn alignment_query_span(a: &AlnRef<'_>) -> Option<(usize, Interval)> {
    let rid = a.ref_id()?;
    let pos = a.pos()?;
    let iv = if a.flag() & FLAG_UNMAPPED != 0 {
        Interval {
            start: pos,
            end: pos + 1,
        }
    } else {
        a.span(u64::MAX)?.0
    };
    (iv.start < iv.end).then_some((rid, iv))
}

/// Parses the records an index-backed BAM reader hands over.
pub(crate) struct BamRegionParse {
    hdr: Arc<AlnHeader>,
    opts: AlignOptions,
    /// Keep every record, placed unmapped reads included, with the span a region is tested
    /// against (`view`); otherwise the records interval commands read (D101).
    all: bool,
    rules: IntervalRules,
    scratch: Vec<u8>,
}

impl BamRegionParse {
    pub(crate) fn new(
        hdr: Arc<AlnHeader>,
        opts: AlignOptions,
        all: bool,
        defer_warnings: bool,
    ) -> BamRegionParse {
        BamRegionParse {
            rules: IntervalRules::new(hdr.refs.len()).deferring(defer_warnings),
            hdr,
            opts,
            all,
            scratch: Vec::new(),
        }
    }

    /// Check the record `buf[off..]` (no `block_size`). `Ok(None)`: it is left out (and
    /// counted in `diag`); the caller drops it. A kept record of an interval reader is
    /// rewritten as the slim layout unless the options ask for whole records.
    pub(crate) fn parse(
        &mut self,
        file: &str,
        buf: &mut Vec<u8>,
        off: usize,
        chroms: &mut ChromTable,
        diag: &mut Diagnostics,
    ) -> Result<Option<Parsed>, String> {
        let a = Aln::parse_framed(&buf[off..], &self.hdr)?;
        if self.all {
            if self.opts.drops(a.flag(), a.mapq()) {
                diag.filtered += 1;
                return Ok(None);
            }
            let Some((rid, iv)) = alignment_query_span(&AlnRef::Bam(a)) else {
                return Ok(None);
            };
            return Ok(Some((ChromId(rid as u32), iv.start, iv, a.strand())));
        }
        let p = match self
            .rules
            .judge(file, &AlnRef::Bam(a), &self.hdr, &self.opts, chroms, diag)
        {
            Verdict::Keep(p) => p,
            Verdict::Skip => return Ok(None),
            Verdict::Bad(msg) => return Err(msg),
        };
        if !self.opts.full_records {
            self.scratch.clear();
            a.write_slim(self.opts.score_tag, &mut self.scratch);
            buf.truncate(off);
            buf.extend_from_slice(&self.scratch);
        }
        Ok(Some((p.chrom, p.pos, p.span, p.strand)))
    }
}

/// The index span of a BAM record, 0-based half-open, and whether it is mapped; `None` for
/// a record that has no place (an unplaced read, counted in the index's `n_no_coor`).
pub(crate) fn index_span(a: &Aln<'_>) -> Option<(usize, Interval, bool)> {
    let (rid, iv) = alignment_query_span(&AlnRef::Bam(*a))?;
    Some((rid, iv, a.flag() & FLAG_UNMAPPED == 0))
}
