//! The BED interval path (D174, R5): one output record per piece. bedGraph takes it too.

use super::{Styler, Unmapped};
use gloci_core::{Interval, parse_coord};
use gloci_io::{ChainIndex, Segment, SourceRef, bed12_blocks};

/// Write `raw` as the record of one piece into `buf` (cleared first): columns 1-3 become `name`
/// and the piece's target interval; column 6, when exactly `+` or `-`, is flipped by a reverse
/// piece; every other column is copied byte for byte. Columns split on tabs only.
pub(crate) fn rewrite_piece(raw: &[u8], name: &str, seg: &Segment, buf: &mut Vec<u8>) {
    buf.clear();
    let mut itoa = itoa::Buffer::new();
    for (i, f) in raw.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            buf.push(b'\t');
        }
        match i {
            0 => buf.extend_from_slice(name.as_bytes()),
            1 => buf.extend_from_slice(itoa.format(seg.target.start).as_bytes()),
            2 => buf.extend_from_slice(itoa.format(seg.target.end).as_bytes()),
            5 if f == b"+" || f == b"-" => {
                let minus = f == b"-";
                buf.push(if seg.reverse == minus { b'+' } else { b'-' });
            }
            _ => buf.extend_from_slice(f),
        }
    }
}

/// The one piece that stands for a record cut into several (D174, R10): CrossMap `region`'s
/// merge. The pieces must cover `min_ratio` of the record, share a target chromosome and a mapped
/// strand; the result spans the first start to the last end.
pub(crate) fn merge_pieces(
    pieces: &[Segment],
    len: u64,
    min_ratio: f64,
    chain: &ChainIndex,
    styler: &mut Styler,
    chr: bool,
) -> Result<Segment, Unmapped> {
    let mapped: u64 = pieces.iter().map(|p| p.source.len()).sum();
    let ratio = if len == 0 {
        1.0
    } else {
        mapped as f64 / len as f64
    };
    if ratio < min_ratio {
        return Err(Unmapped::LowRatio);
    }
    let first = pieces[0];
    if pieces
        .iter()
        .any(|p| p.chain != first.chain && !styler.same_name(chain, first.chain, p.chain, chr))
    {
        return Err(Unmapped::CrossChroms);
    }
    if pieces.iter().any(|p| p.reverse != first.reverse) {
        return Err(Unmapped::MixedStrands);
    }
    let mut merged = first;
    merged.target.start = pieces.iter().map(|p| p.target.start).min().unwrap_or(0);
    merged.target.end = pieces.iter().map(|p| p.target.end).max().unwrap_or(0);
    Ok(merged)
}

/// What became of a BED12 record on the block path.
pub(crate) enum Lift {
    /// Written into the buffer; `shortened` when a block lost bases no chain block covers.
    Written {
        shortened: bool,
    },
    Unmapped(Unmapped),
}

/// Scratch space the block path reuses from record to record.
#[derive(Default)]
pub(crate) struct BlockScratch {
    pieces: Vec<Segment>,
    one: Vec<Segment>,
}

/// The first block `i` for which `hit(src)` holds, else the first for which `edge(src)` does.
fn first_block(
    pieces: &[Segment],
    hit: impl Fn(&Interval) -> bool,
    edge: impl Fn(&Interval) -> bool,
) -> Option<usize> {
    pieces
        .iter()
        .position(|p| hit(&p.source))
        .or_else(|| pieces.iter().position(|p| edge(&p.source)))
}

fn number(raw: &[u8], col: usize, what: &str) -> Result<u64, String> {
    let f = raw.split(|&b| b == b'\t').nth(col).unwrap_or_default();
    parse_coord(f).map_err(|e| format!("{what}: {e}"))
}

fn join_numbers(buf: &mut Vec<u8>, values: impl Iterator<Item = u64>) {
    let mut itoa = itoa::Buffer::new();
    for (i, v) in values.enumerate() {
        if i > 0 {
            buf.push(b',');
        }
        buf.extend_from_slice(itoa.format(v).as_bytes());
    }
}

/// Lift one record of 12 or more columns block by block (D174, R8) into `buf`. An `Err` is an
/// invalid line: its blocks or thick range cannot be read.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lift_blocks(
    raw: &[u8],
    span: Interval,
    chain: &ChainIndex,
    src: Option<&SourceRef>,
    styler: &mut Styler,
    chr: bool,
    scratch: &mut BlockScratch,
    buf: &mut Vec<u8>,
) -> Result<Lift, String> {
    let blocks = bed12_blocks(raw)?.ok_or("fewer than 12 columns")?;
    let thick_start = number(raw, 6, "thickStart")?;
    let thick_end = number(raw, 7, "thickEnd")?;
    for (what, v) in [("thickStart", thick_start), ("thickEnd", thick_end)] {
        if v < span.start || v > span.end {
            return Err(format!(
                "{what} {v} lies outside chromStart {} and chromEnd {}",
                span.start, span.end
            ));
        }
    }
    let minus = raw.split(|&b| b == b'\t').nth(5) == Some(b"-");

    let BlockScratch { pieces, one } = scratch;
    pieces.clear();
    let mut shortened = false;
    // A chromosome the chain lacks, or no block on any chain, is an off-chain record; a
    // block that is not in exactly one piece, with another on a chain, is a split one.
    let Some(src) = src else {
        return Ok(Lift::Unmapped(Unmapped::OffChain));
    };
    let (mut any, mut split) = (false, false);
    for b in &blocks {
        chain.segments(src, *b, one);
        any |= !one.is_empty();
        match one.as_slice() {
            [p] => {
                shortened |= p.source.len() < b.len();
                pieces.push(*p);
            }
            _ => split = true,
        }
    }
    if !any {
        return Ok(Lift::Unmapped(Unmapped::OffChain));
    }
    if split {
        return Ok(Lift::Unmapped(Unmapped::BlockSplit));
    }
    let first = pieces[0];
    for p in pieces.iter() {
        let same_name =
            p.chain == first.chain || styler.same_name(chain, first.chain, p.chain, chr);
        if !same_name || (p.reverse == minus) != (first.reverse == minus) {
            return Ok(Lift::Unmapped(Unmapped::BlocksApart));
        }
    }

    let cs = first.target.start;
    let ce = pieces[pieces.len() - 1].target.end;
    let found = |i: Option<usize>| i.map(|i| &pieces[i]);
    let ts_block = found(first_block(
        pieces,
        |s| s.start <= thick_start && thick_start < s.end,
        |s| thick_start == s.end,
    ));
    let te_block = found(first_block(
        pieces,
        |s| s.start < thick_end && thick_end <= s.end,
        |s| thick_end == s.start,
    ));
    let fwd = |p: &Segment, x: u64| i128::from(p.target.start) + i128::from(x - p.source.start);
    let mirror = |p: &Segment, x: u64| i128::from(p.target.end) - i128::from(x - p.source.start);
    let mapped = match (ts_block, te_block) {
        (Some(s), Some(e)) if !first.reverse => Some((fwd(s, thick_start), fwd(e, thick_end))),
        (Some(s), Some(e)) => Some((mirror(e, thick_end), mirror(s, thick_start))),
        _ => None,
    };
    let (ts, te) = mapped.unwrap_or((
        i128::from(cs) + i128::from(thick_start) - i128::from(span.start),
        i128::from(ce) - (i128::from(span.end) - i128::from(thick_end)),
    ));
    let ascending = pieces
        .windows(2)
        .all(|w| w[0].target.end <= w[1].target.start);
    if cs > ce || ts > te || i128::from(cs) > ts || te > i128::from(ce) || !ascending {
        return Ok(Lift::Unmapped(Unmapped::BadBlocks));
    }

    buf.clear();
    let mut itoa = itoa::Buffer::new();
    let name = styler.name(chain, first.chain, chr);
    for (i, f) in raw.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            buf.push(b'\t');
        }
        match i {
            0 => buf.extend_from_slice(name.as_bytes()),
            1 => buf.extend_from_slice(itoa.format(cs).as_bytes()),
            2 => buf.extend_from_slice(itoa.format(ce).as_bytes()),
            5 if f == b"+" || f == b"-" => {
                buf.push(if first.reverse == minus { b'+' } else { b'-' });
            }
            6 => buf.extend_from_slice(itoa.format(ts).as_bytes()),
            9 => buf.extend_from_slice(itoa.format(pieces.len()).as_bytes()),
            7 => buf.extend_from_slice(itoa.format(te).as_bytes()),
            10 => join_numbers(buf, pieces.iter().map(|p| p.target.len())),
            11 => join_numbers(buf, pieces.iter().map(|p| p.target.start - cs)),
            _ => buf.extend_from_slice(f),
        }
    }
    Ok(Lift::Written { shortened })
}
