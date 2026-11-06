//! Overlap rules (spec 4.3): point semantics, exact fractions, strand, windows.

use crate::error::FractionError;
use crate::interval::{Interval, Strand};

/// Do `a` and `b` overlap? A point is the interval `[p, p)`.
///
/// Two non-empty intervals must share a base. A point overlaps a non-empty
/// interval when strictly inside it (`s < p < e`), or `s <= p <= e` with
/// `inclusive_points`. Two points overlap when equal.
pub fn overlaps(a: Interval, b: Interval, inclusive_points: bool) -> bool {
    match (a.is_point(), b.is_point()) {
        (true, true) => a.start == b.start,
        (true, false) => point_in(a.start, b, inclusive_points),
        (false, true) => point_in(b.start, a, inclusive_points),
        (false, false) => a.start < b.end && b.start < a.end,
    }
}

fn point_in(p: u64, iv: Interval, inclusive: bool) -> bool {
    if inclusive {
        iv.start <= p && p <= iv.end
    } else {
        iv.start < p && p < iv.end
    }
}

/// Shared bases of `a` and `b`; 0 if either is a point.
pub fn overlap_len(a: Interval, b: Interval) -> u64 {
    if a.is_point() || b.is_point() {
        return 0;
    }
    a.end.min(b.end).saturating_sub(a.start.max(b.start))
}

/// An exact decimal fraction in (0, 1], stored as `num / 10^digits`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fraction {
    num: u64,
    digits: u32,
}

impl Fraction {
    pub fn parse(s: &str) -> Result<Fraction, FractionError> {
        let syntax = || FractionError::Syntax(s.to_string());
        let (int, frac) = match s.split_once('.') {
            Some((i, f)) => (i, f),
            None => (s, ""),
        };
        if int.is_empty() && frac.is_empty() {
            return Err(syntax());
        }
        if !int.bytes().chain(frac.bytes()).all(|b| b.is_ascii_digit()) {
            return Err(syntax());
        }
        if frac.len() > 18 {
            return Err(FractionError::TooPrecise(s.to_string()));
        }
        let digits = frac.len() as u32;
        let scale = 10u128.pow(digits);
        let out_of_range = || FractionError::OutOfRange(s.to_string());
        let mut num: u128 = 0;
        for b in int.bytes().chain(frac.bytes()) {
            num = num
                .checked_mul(10)
                .and_then(|n| n.checked_add(u128::from(b - b'0')))
                .filter(|&n| n <= scale)
                .ok_or_else(out_of_range)?;
        }
        if num == 0 || num > scale {
            return Err(out_of_range());
        }
        Ok(Fraction {
            num: num as u64,
            digits,
        })
    }

    /// Is `part / whole >= self`? Exact integer arithmetic.
    pub fn met(&self, part: u64, whole: u64) -> bool {
        u128::from(part) * 10u128.pow(self.digits) >= u128::from(self.num) * u128::from(whole)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StrandRule {
    #[default]
    Any,
    Same,
    Opposite,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OverlapCriteria {
    pub min_frac_a: Option<Fraction>,
    pub min_frac_b: Option<Fraction>,
    pub reciprocal: bool,
    pub either: bool,
    pub min_bases: Option<u64>,
    pub strand: StrandRule,
    pub inclusive_points: bool,
}

impl OverlapCriteria {
    pub fn accepts(&self, a: Interval, a_strand: Strand, b: Interval, b_strand: Strand) -> bool {
        if !overlaps(a, b, self.inclusive_points) {
            return false;
        }
        let ov = overlap_len(a, b);

        let min_b = if self.reciprocal && self.min_frac_b.is_none() {
            self.min_frac_a
        } else {
            self.min_frac_b
        };
        // A threshold on a point side counts as met. With `either`, only the
        // thresholds actually given take part, so one threshold reduces to itself.
        let fa = self.min_frac_a.map(|f| a.is_point() || f.met(ov, a.len()));
        let fb = min_b.map(|f| b.is_point() || f.met(ov, b.len()));
        let frac_ok = match (fa, fb) {
            (None, None) => true,
            (Some(x), None) | (None, Some(x)) => x,
            (Some(x), Some(y)) => {
                if self.either {
                    x || y
                } else {
                    x && y
                }
            }
        };

        let bases_ok = self
            .min_bases
            .is_none_or(|n| a.is_point() || b.is_point() || ov >= n);

        let strand_ok = match self.strand {
            StrandRule::Any => true,
            rule => {
                a_strand != Strand::Unknown
                    && b_strand != Strand::Unknown
                    && ((a_strand == b_strand) == (rule == StrandRule::Same))
            }
        };

        frac_ok && bases_ok && strand_ok
    }
}

/// Padding applied to B before testing (`-w`, `--left`, `--right`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Window {
    pub left: u64,
    pub right: u64,
    pub strand_aware: bool,
}

impl Window {
    pub fn apply(&self, b: Interval, a_strand: Strand) -> Interval {
        let (l, r) = if self.strand_aware && a_strand == Strand::Minus {
            (self.right, self.left)
        } else {
            (self.left, self.right)
        };
        Interval {
            start: b.start.saturating_sub(l),
            end: b.end.saturating_add(r),
        }
    }
}
