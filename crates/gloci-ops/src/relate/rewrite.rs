//! Output helpers for commands that change or invent coordinates.

use gloci_core::Interval;
use gloci_io::Format;

/// Copy `raw` with its coordinates replaced by `iv` (D44): BED columns 2-3,
/// GFF/GTF columns 4-5 as 1-based closed. Other columns, including BED12
/// thick and block columns, are left alone. VCF/BCF records are copied
/// unchanged; callers refuse to clip them.
pub fn rewrite_coords(format: Format, raw: &[u8], iv: Interval, out: &mut Vec<u8>) {
    let (first, start, end) = match format {
        Format::Bed | Format::BedGraph => (1, iv.start, iv.end),
        Format::Gff3 | Format::Gtf => (3, iv.start + 1, iv.end),
        Format::Vcf | Format::Bcf => {
            out.extend_from_slice(raw);
            return;
        }
    };
    for (i, f) in raw.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        if i == first {
            out.extend_from_slice(start.to_string().as_bytes());
        } else if i == first + 1 {
            out.extend_from_slice(end.to_string().as_bytes());
        } else {
            out.extend_from_slice(f);
        }
    }
}

/// The columns written for a missing B (D45): for BED-family B, `.`, `-1`,
/// `-1`, then `.` for each further column (`b_columns` is B's column
/// count); for other formats, `.` for every column.
pub fn b_placeholder(b_format: Format, b_columns: usize, out: &mut Vec<u8>) {
    let mut fields: Vec<&str> = Vec::new();
    match b_format {
        Format::Bed | Format::BedGraph => {
            fields.extend([".", "-1", "-1"]);
            fields.resize(b_columns.max(3), ".");
        }
        _ => fields.resize(b_columns, "."),
    }
    out.extend_from_slice(fields.join("\t").as_bytes());
}
