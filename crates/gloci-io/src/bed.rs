//! BED and bedGraph reader (bedGraph is BED with a value column).

use crate::chrom_map::ChromMap;
use crate::diag::InvalidPolicy;
use crate::error::IoError;
use crate::input::Input;
use crate::record::{Class, Parsed, TextKind, TextReader, coord, intern_chrom};
use gloci_core::{ChromTable, Interval, Strand};
use std::sync::Arc;

#[derive(Default, Clone)]
pub(crate) struct Bed;
pub(crate) type BedReader = TextReader<Bed>;

impl TextReader<Bed> {
    pub(crate) fn open_bed(
        input: Input,
        policy: InvalidPolicy,
        map: Option<Arc<ChromMap>>,
    ) -> Result<Self, IoError> {
        let mut r = Self::open(input, policy, Bed)?;
        r.set_chrom_map(map);
        Ok(r)
    }
}

fn is_word(line: &[u8], w: &[u8]) -> bool {
    line.starts_with(w) && line.get(w.len()).is_none_or(|b| *b == b' ' || *b == b'\t')
}

pub(crate) fn is_header_line(line: &[u8]) -> bool {
    line.first() == Some(&b'#') || is_word(line, b"track") || is_word(line, b"browser")
}

/// Is `line` (without its line terminator) a BED record, rather than a
/// blank, `#`, `track` or `browser` line (D14)?
pub fn is_bed_record_line(line: &[u8]) -> bool {
    !line.is_empty() && !is_header_line(line)
}

impl TextKind for Bed {
    fn classify(line: &[u8], _seen: bool) -> Class {
        if line.is_empty() {
            Class::Skip
        } else if is_header_line(line) {
            Class::Header
        } else {
            Class::Record
        }
    }

    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let mut f = line.split(|&b| b == b'\t');
        let (Some(c), Some(s), Some(e)) = (f.next(), f.next(), f.next()) else {
            let n = line.split(|&b| b == b'\t').count();
            let hint = if memchr::memchr(b' ', line).is_some() {
                " (the line contains spaces; BED must be tab-separated)"
            } else {
                ""
            };
            return Err(format!(
                "expected at least 3 tab-separated columns, found {n}{hint}"
            ));
        };
        let start = coord("start", s)?;
        let end = coord("end", e)?;
        let span = Interval::new(start, end).map_err(|e| e.to_string())?;
        let chrom = intern_chrom(chroms, c)?;
        // Column 6 counts only as exactly `+`, `-` or `.` (D54).
        let strand = match f.nth(2) {
            Some(b"+") => Strand::Plus,
            Some(b"-") => Strand::Minus,
            _ => Strand::Unknown,
        };
        Ok((chrom, start, span, strand))
    }
}

/// Blocks of a BED12 line as absolute intervals; `None` below 12 columns.
pub fn bed12_blocks(line: &[u8]) -> Result<Option<Vec<Interval>>, String> {
    let f: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
    if f.len() < 12 {
        return Ok(None);
    }
    let start = coord("start", f[1])?;
    let end = coord("end", f[2])?;
    let count = coord("blockCount", f[9])? as usize;
    let list = |name: &str, v: &[u8]| -> Result<Vec<u64>, String> {
        let v = v.strip_suffix(b",").unwrap_or(v);
        if v.is_empty() {
            return Ok(Vec::new());
        }
        v.split(|&b| b == b',').map(|x| coord(name, x)).collect()
    };
    let sizes = list("blockSizes", f[10])?;
    let starts = list("blockStarts", f[11])?;
    if sizes.len() != count || starts.len() != count {
        return Err(format!(
            "blockCount is {count} but there are {} blockSizes and {} blockStarts",
            sizes.len(),
            starts.len()
        ));
    }
    if count == 0 {
        return Err("blockCount is 0".into());
    }
    if starts[0] != 0 {
        return Err("the first block must start at offset 0".into());
    }
    let mut blocks: Vec<Interval> = Vec::with_capacity(count);
    for (&rel, &size) in starts.iter().zip(&sizes) {
        let s = start.checked_add(rel).ok_or("block start overflows")?;
        let e = s.checked_add(size).ok_or("block end overflows")?;
        if blocks.last().is_some_and(|p| s < p.end) {
            return Err("blocks are unsorted or overlap".into());
        }
        blocks.push(Interval { start: s, end: e });
    }
    if blocks.last().map(|b| b.end) != Some(end) {
        return Err(format!(
            "the last block ends at {} but chromEnd is {end}",
            blocks.last().map_or(0, |b| b.end)
        ));
    }
    Ok(Some(blocks))
}
