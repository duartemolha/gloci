//! The BED interval path (D174, R5): one output record per piece. bedGraph takes it too.

use gloci_io::Segment;

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
