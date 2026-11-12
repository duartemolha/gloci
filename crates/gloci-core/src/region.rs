//! Region strings (decision D8).

use crate::error::{CoordError, RegionError};
use crate::interval::{Interval, one_based_closed, parse_coord};
use crate::overlap::overlaps;

/// A 0-based region; `end: None` means to the chromosome end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub chrom: Vec<u8>,
    pub start: u64,
    pub end: Option<u64>,
}

impl Region {
    /// Accepts `chr`, `chr:s`, `chr:s-`, `chr:-e` (= `chr:1-e`, as htslib), `chr:s-e` and `{name}:s-e` (1-based closed,
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
            if let Some(i) = split
                && is_known(&s.as_bytes()[..i])
                && parse_range(s, &s[..i], &s[i + 1..]).is_ok()
            {
                return Err(RegionError::Ambiguous(s.into()));
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

/// The regions requested on one chromosome, sorted and without duplicates
/// (D10). A span matches when it overlaps one of them under §4.3. They are
/// not merged into their union: a point between two bookended regions is in
/// neither.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RegionSet {
    regions: Vec<Interval>,
    /// `reach[i]`: the largest end among `regions[..=i]`.
    reach: Vec<u64>,
}

impl RegionSet {
    pub fn new(mut regions: Vec<Interval>) -> RegionSet {
        regions.sort_unstable();
        regions.dedup();
        let mut reach = Vec::with_capacity(regions.len());
        let mut m = 0;
        for r in &regions {
            m = r.end.max(m);
            reach.push(m);
        }
        RegionSet { regions, reach }
    }

    /// The regions, sorted by start.
    pub fn regions(&self) -> &[Interval] {
        &self.regions
    }

    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    /// The largest region end (0 when there are none). No span that starts
    /// after it can match.
    pub fn end(&self) -> u64 {
        self.reach.last().copied().unwrap_or(0)
    }

    /// Does `span` overlap one of the regions under §4.3?
    pub fn matches(&self, span: Interval, inclusive_points: bool) -> bool {
        let n = self.regions.partition_point(|r| r.start <= span.end);
        for i in (0..n).rev() {
            if self.reach[i] < span.start {
                break;
            }
            if overlaps(span, self.regions[i], inclusive_points) {
                return true;
            }
        }
        false
    }
}

fn whole(name: &str) -> Region {
    Region {
        chrom: name.as_bytes().to_vec(),
        start: 0,
        end: None,
    }
}

fn validate_commas(s: &str) -> bool {
    if !s.contains(',') {
        return true;
    }
    let mut parts = s.split(',');
    let Some(first) = parts.next() else {
        return false;
    };
    if first.is_empty() || first.len() > 3 || !first.chars().all(|c| c.is_ascii_digit()) {
        return false;
    }
    let mut count = 0;
    for part in parts {
        count += 1;
        if part.len() != 3 || !part.chars().all(|c| c.is_ascii_digit()) {
            return false;
        }
    }
    count > 0
}

fn parse_range(full: &str, chrom: &str, range: &str) -> Result<Region, RegionError> {
    let wrap = |source: CoordError| RegionError::Coord {
        region: full.into(),
        source,
    };
    let num = |t: &str| {
        if !validate_commas(t) {
            return Err(wrap(CoordError::MalformedCommas(t.into())));
        }
        parse_coord(t.replace(',', "").as_bytes()).map_err(wrap)
    };
    let (start, end) = match range.split_once('-') {
        None => (num(range)?, None),
        // htslib: `chr:-E` is `chr:1-E`.
        // (`chr:-0` stays an error: it would be an empty region.)
        Some(("", b)) if !b.is_empty() => {
            let e = num(b)?;
            if e == 0 {
                return Err(wrap(
                    one_based_closed(0, 0).expect_err("start 0 is invalid"),
                ));
            }
            (1, Some(e))
        }
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
