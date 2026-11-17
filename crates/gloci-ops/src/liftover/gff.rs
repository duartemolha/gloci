//! GFF3 and GTF features (D174, R9): a feature lifts only when it lies in exactly one piece.

use super::{Styler, Unmapped};
use gloci_core::Interval;
use gloci_io::{ChainIndex, Segment, SourceRef};

/// What became of a feature.
pub(crate) enum Lift {
    /// Written into the buffer; `partial` when the piece is shorter than the feature, which
    /// CrossMap also writes to the unmapped output.
    Written {
        partial: bool,
    },
    Unmapped(Unmapped),
}

/// Lift one feature into `buf`: columns 1, 4, 5 and 7 are replaced, the rest is copied.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lift(
    raw: &[u8],
    span: Interval,
    chain: &ChainIndex,
    src: Option<&SourceRef>,
    styler: &mut Styler,
    chr: bool,
    pieces: &mut Vec<Segment>,
    buf: &mut Vec<u8>,
) -> Lift {
    let strand = raw.split(|&b| b == b'\t').nth(6).unwrap_or_default();
    if !matches!(strand, b"+" | b"-" | b".") {
        return Lift::Unmapped(Unmapped::Strand);
    }
    pieces.clear();
    if let Some(src) = src {
        chain.segments(src, span, pieces);
    }
    let seg = match pieces.as_slice() {
        [] => return Lift::Unmapped(Unmapped::OffChain),
        [one] => one,
        _ => return Lift::Unmapped(Unmapped::SeveralPieces),
    };
    let name = styler.name(chain, seg.chain, chr);
    buf.clear();
    let mut itoa = itoa::Buffer::new();
    for (i, f) in raw.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            buf.push(b'\t');
        }
        match i {
            0 => buf.extend_from_slice(name.as_bytes()),
            3 => buf.extend_from_slice(itoa.format(seg.target.start + 1).as_bytes()),
            4 => buf.extend_from_slice(itoa.format(seg.target.end).as_bytes()),
            6 if strand != b"." => {
                buf.push(if seg.reverse == (strand == b"-") {
                    b'+'
                } else {
                    b'-'
                });
            }
            _ => buf.extend_from_slice(f),
        }
    }
    Lift::Written {
        partial: seg.source.len() < span.len(),
    }
}
