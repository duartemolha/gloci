//! Region strings (decision D8).

use crate::error::{CoordError, RegionError};
use crate::interval::{Interval, one_based_closed, parse_coord};

/// A 0-based region; `end: None` means to the chromosome end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub chrom: Vec<u8>,
    pub start: u64,
    pub end: Option<u64>,
}

impl Region {
    /// Accepts `chr`, `chr:s`, `chr:s-`, `chr:s-e` and `{name}:s-e` (1-based closed,
    /// commas allowed). A whole-string known chromosome name wins; if its prefix before
    /// the last `:` is also a known name the string is ambiguous.
    pub fn parse(s: &str, is_known: impl Fn(&[u8]) -> bool) -> Result<Region, RegionError> {
        if s.is_empty() {
            return Err(RegionError::Empty);
        }
        if let Some(rest) = s.strip_prefix('{') {
            let close = rest
                .find('}')
                .ok_or_else(|| RegionError::Unterminated(s.into()))?;
            let name = &rest[..close];
            if name.is_empty() {
                return Err(RegionError::MissingChrom(s.into()));
            }
            return match rest[close + 1..].strip_prefix(':') {
                Some(range) => parse_range(s, name, range),
                None if rest[close + 1..].is_empty() => Ok(whole(name)),
                None => Err(RegionError::Unterminated(s.into())),
            };
        }
        let split = s.rfind(':');
        if is_known(s.as_bytes()) {
            if let Some(i) = split {
                if is_known(&s.as_bytes()[..i]) && parse_range(s, &s[..i], &s[i + 1..]).is_ok() {
                    return Err(RegionError::Ambiguous(s.into()));
                }
            }
            return Ok(whole(s));
        }
        match split {
            None => Ok(whole(s)),
            Some(0) => Err(RegionError::MissingChrom(s.into())),
            Some(i) => parse_range(s, &s[..i], &s[i + 1..]),
        }
    }

    /// The region as an interval; an open end becomes `chrom_len`, or `u64::MAX` if unknown.
    pub fn interval(&self, chrom_len: Option<u64>) -> Interval {
        let limit = chrom_len.unwrap_or(u64::MAX);
        let end = self.end.unwrap_or(limit).min(limit).max(self.start);
        Interval {
            start: self.start,
            end,
        }
    }
}

fn whole(name: &str) -> Region {
    Region {
        chrom: name.as_bytes().to_vec(),
        start: 0,
        end: None,
    }
}

fn parse_range(full: &str, chrom: &str, range: &str) -> Result<Region, RegionError> {
    let wrap = |source: CoordError| RegionError::Coord {
        region: full.into(),
        source,
    };
    let num = |t: &str| parse_coord(t.replace(',', "").as_bytes()).map_err(wrap);
    let (start, end) = match range.split_once('-') {
        None => (num(range)?, None),
        Some((a, "")) => (num(a)?, None),
        Some((a, b)) => (num(a)?, Some(num(b)?)),
    };
    let (start, end) = match end {
        Some(e) => {
            let iv = one_based_closed(start, e).map_err(wrap)?;
            (iv.start, Some(iv.end))
        }
        None => (one_based_closed(start, start).map_err(wrap)?.start, None),
    };
    Ok(Region {
        chrom: chrom.as_bytes().to_vec(),
        start,
        end,
    })
}
