//! Half-open 0-based intervals, strands and coordinate parsing.

use crate::error::CoordError;

/// A 0-based, half-open interval `[start, end)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn new(start: u64, end: u64) -> Result<Interval, CoordError> {
        if start > end {
            return Err(CoordError::StartAfterEnd { start, end });
        }
        Ok(Interval { start, end })
    }

    /// The zero-length interval `[p, p)`.
    pub fn point(p: u64) -> Interval {
        Interval { start: p, end: p }
    }

    pub fn is_point(&self) -> bool {
        self.start == self.end
    }

    #[allow(clippy::len_without_is_empty)] // emptiness is spelled `is_point`
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    /// `+`, `-`, `.`; GFF `?` maps to Unknown.
    pub fn parse(field: &[u8]) -> Option<Strand> {
        match field {
            b"+" => Some(Strand::Plus),
            b"-" => Some(Strand::Minus),
            b"." | b"?" => Some(Strand::Unknown),
            _ => None,
        }
    }
}

/// Parse an unsigned decimal coordinate (leading zeros allowed, no sign).
pub fn parse_coord(field: &[u8]) -> Result<u64, CoordError> {
    let shown = || String::from_utf8_lossy(field).into_owned();
    let digits = field.strip_prefix(b"-").unwrap_or(field);
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(CoordError::NotInteger(shown()));
    }
    if digits.len() != field.len() {
        return Err(CoordError::Negative(shown()));
    }
    let mut v: u64 = 0;
    for &d in digits {
        v = v
            .checked_mul(10)
            .and_then(|v| v.checked_add(u64::from(d - b'0')))
            .ok_or_else(|| CoordError::TooLarge(shown()))?;
    }
    Ok(v)
}

/// Convert 1-based closed `[s, e]` to `[s-1, e)`.
pub fn one_based_closed(s: u64, e: u64) -> Result<Interval, CoordError> {
    if s < 1 {
        return Err(CoordError::StartBelowOne);
    }
    Interval::new(s - 1, e)
}
