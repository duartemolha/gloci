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

/// An exact decimal fraction in (0, 1], stored as `num / scale` (`scale` is a
/// power of ten). Equality compares by value: `1` equals `1.0`, `50%` equals `0.5`.
#[derive(Clone, Copy, Debug)]
pub struct Fraction {
    num: u64,
    scale: u128,
}

impl PartialEq for Fraction {
    fn eq(&self, other: &Self) -> bool {
        u128::from(self.num) * other.scale == u128::from(other.num) * self.scale
    }
}

impl Eq for Fraction {}

impl Fraction {
    pub fn parse(s: &str) -> Result<Fraction, FractionError> {
        Self::parse_inner(s, true)
    }

    /// Like [`Fraction::parse`] for a scale factor: any non-negative decimal
    /// (`0`, `0.5`, `2`, `150%`) that fits `u64` once scaled. Used by `slop`
    /// and `flank --pct`, where `-b 2` means twice the feature length.
    pub fn parse_ratio(s: &str) -> Result<Fraction, FractionError> {
        Self::parse_inner(s, false)
    }

    /// `floor(self * len)`, exact, saturating at `u64::MAX`.
    pub fn floor_mul(&self, len: u64) -> u64 {
        let v = u128::from(self.num) * u128::from(len) / self.scale;
        u64::try_from(v).unwrap_or(u64::MAX)
    }

    fn parse_inner(s: &str, unit: bool) -> Result<Fraction, FractionError> {
        let syntax = || FractionError::Syntax(s.to_string());
        // A trailing `%` divides by 100: "50%" is 0.5.
        let (body, extra_digits) = match s.strip_suffix('%') {
            Some(b) => (b, 2),
            None => (s, 0),
        };
        let (int, frac) = match body.split_once('.') {
            Some((i, f)) => (i, f),
            None => (body, ""),
        };
        if int.is_empty() && frac.is_empty() {
            return Err(syntax());
        }
        if !int.bytes().chain(frac.bytes()).all(|b| b.is_ascii_digit()) {
            return Err(syntax());
        }
        if frac.len() + extra_digits > 18 {
            return Err(FractionError::TooPrecise(s.to_string()));
        }
        let digits = (frac.len() + extra_digits) as u32;
        let scale = 10u128.pow(digits);
        let out_of_range = || {
            if unit {
                FractionError::OutOfRange(s.to_string())
            } else {
                FractionError::RatioRange(s.to_string())
            }
        };
        let mut num: u128 = 0;
        for b in int.bytes().chain(frac.bytes()) {
            num = num
                .checked_mul(10)
                .and_then(|n| n.checked_add(u128::from(b - b'0')))
                .filter(|&n| {
                    if unit {
                        n <= scale
                    } else {
                        n <= u128::from(u64::MAX)
                    }
                })
                .ok_or_else(out_of_range)?;
        }
        if unit && (num == 0 || num > scale) {
            return Err(out_of_range());
        }
        Ok(Fraction {
            num: num as u64,
            scale,
        })
    }

    /// Is `part / whole >= self`? Exact integer arithmetic.
    pub fn met(&self, part: u64, whole: u64) -> bool {
        u128::from(part) * self.scale >= u128::from(self.num) * u128::from(whole)
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

        // D34 (spec 4.3, amended 2026-11-08): every threshold counts as met when
        // either side is a point. `-r` with one fraction applies it to both sides.
        let point = a.is_point() || b.is_point();
        let (min_a, min_b) = if self.reciprocal {
            let one = self.min_frac_a.or(self.min_frac_b);
            (self.min_frac_a.or(one), self.min_frac_b.or(one))
        } else {
            (self.min_frac_a, self.min_frac_b)
        };
        // With `either`, only the thresholds actually given take part, so one
        // threshold reduces to itself.
        let fa = min_a.map(|f| point || f.met(ov, a.len()));
        let fb = min_b.map(|f| point || f.met(ov, b.len()));
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

        let bases_ok = self.min_bases.is_none_or(|n| point || ov >= n);

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

/// Do any block of `a` and any block of `b` overlap (spec 4.3 point rules)?
pub fn blocks_overlap(a: &[Interval], b: &[Interval], inclusive_points: bool) -> bool {
    a.iter()
        .any(|&x| b.iter().any(|&y| overlaps(x, y, inclusive_points)))
}

/// Total shared bases between two block lists, each sorted and disjoint.
pub fn blocks_overlap_len(a: &[Interval], b: &[Interval]) -> u64 {
    let (mut i, mut j, mut total) = (0, 0, 0u64);
    while i < a.len() && j < b.len() {
        total += overlap_len(a[i], b[j]);
        if a[i].end <= b[j].end {
            i += 1;
        } else {
            j += 1;
        }
    }
    total
}
