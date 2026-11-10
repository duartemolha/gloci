//! What an interval command makes of one alignment (D101, D102): skipped (unmapped, filtered,
//! off the end of its reference), or kept as a chromosome, a span and a strand. Shared by the
//! BAM and SAM readers, so both skip, count and warn alike.

use super::{AlignOptions, AlnHeader, AlnRef};
use crate::diag::Diagnostics;
use crate::record::intern_chrom;
use gloci_core::{ChromId, ChromTable, Interval, Strand};

const FLAG_UNMAPPED: u16 = 0x4;

/// A kept alignment as an interval.
pub(crate) struct Placed {
    pub chrom: ChromId,
    pub pos: u64,
    pub span: Interval,
    pub strand: Strand,
}

pub(crate) enum Verdict {
    Keep(Placed),
    /// Left out, and counted in the diagnostics.
    Skip,
    /// The chromosome cannot be used (a closed `-g` table without it).
    Bad(String),
}

pub(crate) struct IntervalRules {
    /// Reference id to chromosome id, filled as references are first met.
    ids: Vec<Option<ChromId>>,
    clip_warned: bool,
    off_ref_warned: bool,
}

impl IntervalRules {
    pub fn new(refs: usize) -> IntervalRules {
        IntervalRules {
            ids: vec![None; refs],
            clip_warned: false,
            off_ref_warned: false,
        }
    }

    /// Judge `a`, a record of the input `file`, against `hdr`.
    pub fn judge(
        &mut self,
        file: &str,
        a: &AlnRef<'_>,
        hdr: &AlnHeader,
        opts: &AlignOptions,
        chroms: &mut ChromTable,
        diag: &mut Diagnostics,
    ) -> Verdict {
        let flag = a.flag();
        if flag & FLAG_UNMAPPED != 0 {
            diag.unmapped += 1;
            return Verdict::Skip;
        }
        if opts.drops(flag, a.mapq()) {
            diag.filtered += 1;
            return Verdict::Skip;
        }
        // `parse` guarantees a mapped read has a reference and a position.
        let (Some(rid), Some(pos)) = (a.ref_id(), a.pos()) else {
            return Verdict::Skip;
        };
        let (name, len) = &hdr.refs[rid];
        let Some((span, clipped)) = a.span(*len) else {
            return Verdict::Skip;
        };
        if span.start >= span.end {
            // The read starts at or past the end of its reference: it covers no base.
            diag.off_reference += 1;
            if !self.off_ref_warned {
                self.off_ref_warned = true;
                log::warn!(
                    "{file}: read {} starts at {} on {}, which is {len} bases long; reads that \
start past the end of their reference cover no base and are skipped",
                    String::from_utf8_lossy(a.name()),
                    pos + 1,
                    String::from_utf8_lossy(name),
                );
            }
            return Verdict::Skip;
        }
        if clipped && !self.clip_warned {
            self.clip_warned = true;
            log::warn!(
                "{file}: read {} on {} extends past the end of the reference ({len} bases); \
reads like it are cut at the reference length",
                String::from_utf8_lossy(a.name()),
                String::from_utf8_lossy(name),
            );
        }
        let chrom = match self.ids[rid] {
            Some(c) => c,
            None => match intern_chrom(chroms, name) {
                Ok(c) => {
                    self.ids[rid] = Some(c);
                    c
                }
                Err(msg) => return Verdict::Bad(msg),
            },
        };
        Verdict::Keep(Placed {
            chrom,
            pos,
            span,
            strand: a.strand(),
        })
    }
}
