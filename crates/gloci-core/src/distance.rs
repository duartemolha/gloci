//! Distance between intervals (spec 4.4).

use crate::interval::{Interval, Strand};
use crate::overlap::overlaps;

/// 0 when overlapping, otherwise the gap in bases plus 1 (bookended = 1).
pub fn distance(a: Interval, b: Interval, inclusive_points: bool) -> u64 {
    if overlaps(a, b, inclusive_points) {
        return 0;
    }
    let gap = b
        .start
        .saturating_sub(a.end)
        .max(a.start.saturating_sub(b.end));
    gap.saturating_add(1)
}

/// Reference frame for signed distances (`--signed ref|a|b`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SignFrame {
    Ref,
    A,
    B,
}

/// Signed distance; unknown strand counts as `+` (D6).
pub fn signed_distance(
    a: Interval,
    a_strand: Strand,
    b: Interval,
    b_strand: Strand,
    frame: SignFrame,
    inclusive_points: bool,
) -> i64 {
    let d = distance(a, b, inclusive_points);
    if d == 0 {
        return 0;
    }
    let d = i64::try_from(d).unwrap_or(i64::MAX);
    let b_lower = b.end <= a.start;
    let neg = match frame {
        SignFrame::Ref => b_lower,
        SignFrame::A => b_lower ^ (a_strand == Strand::Minus),
        SignFrame::B => !b_lower ^ (b_strand == Strand::Minus),
    };
    if neg { -d } else { d }
}
