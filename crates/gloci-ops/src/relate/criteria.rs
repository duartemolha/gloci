//! Which candidate hits count (D36, D37).

use super::RelateOptions;
use gloci_core::{
    Interval, OverlapCriteria, Strand, StrandRule, blocks_overlap, blocks_overlap_len,
};
use gloci_engine::{ARec, Hit};
use gloci_io::bed12_blocks;

/// Blocks of a record: its BED12 blocks, else its span as the one block.
/// A BED12 line whose blocks do not parse also counts as one block.
fn blocks(raw: &[u8], span: Interval) -> Vec<Interval> {
    match bed12_blocks(raw) {
        Ok(Some(b)) => b,
        _ => vec![span],
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

/// Does `hit` pass the criteria for `a`?
///
/// - **Window.** The engine already found the hit by its padded span, and
///   thresholds cannot be combined with a window (D37), so only the strand
///   rule is left to apply (with `--split`, the padded blocks must meet).
/// - **Otherwise** the thresholds are measured against B's own span. With
///   `--split`, A and B are their blocks: they must share a base (or
///   point), and thresholds use the total block lengths and the total
///   overlap, `blocks_overlap_len`.
pub fn accepted(a: &ARec, hit: &Hit, o: &RelateOptions) -> bool {
    let c = &o.criteria;
    let (a_span, b_span) = (a.meta.span, hit.span);
    if !o.split {
        return match o.window {
            Some(_) => strand_ok(c.strand, a.meta.strand, hit.strand),
            None => c.accepts(a_span, a.meta.strand, b_span, hit.strand),
        };
    }
    let ab = blocks(a.raw, a_span);
    let mut bb = blocks(hit.raw, b_span);
    if let Some(w) = o.window {
        for x in &mut bb {
            *x = w.apply(*x, a.meta.strand);
        }
        return blocks_overlap(&ab, &bb, c.inclusive_points)
            && strand_ok(c.strand, a.meta.strand, hit.strand);
    }
    if !blocks_overlap(&ab, &bb, c.inclusive_points) {
        return false;
    }
    let total = |v: &[Interval]| v.iter().map(Interval::len).sum::<u64>();
    let (alen, blen, ov) = (total(&ab), total(&bb), blocks_overlap_len(&ab, &bb));
    // `accepts` measures two intervals, so give it two with these sizes and
    // this overlap. Overlapping blocks that share no base are points
    // (D34), and those thresholds count as met.
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
    c.accepts(sa, a.meta.strand, sb, hit.strand)
}
