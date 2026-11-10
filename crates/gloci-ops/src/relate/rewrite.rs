//! Output helpers for commands that change or invent coordinates.

use gloci_core::Interval;
use gloci_io::Format;

/// Copy `raw` with its coordinates replaced by `iv` (D44): BED columns 2-3,
/// GFF/GTF columns 4-5 as 1-based closed. Other columns, including BED12
/// thick and block columns, are left alone. VCF/BCF records are copied
/// unchanged; callers refuse to clip them.
pub fn rewrite_coords(format: Format, raw: &[u8], iv: Interval, out: &mut Vec<u8>) {
    let (first, start, end) = match format {
        Format::Bed | Format::BedGraph | Format::Wig => (1, iv.start, iv.end),
        Format::Gff3 | Format::Gtf => (3, iv.start + 1, iv.end),
        Format::Vcf
        | Format::Bcf
        | Format::Bam
        | Format::Cram
        | Format::Sam
        | Format::Fasta
        | Format::Psl
        | Format::Rmsk => {
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

/// The columns written for a missing B (D45, as bedtools writes them,
/// verified against 2.31.1): `b_columns` is B's column count.
///
/// - **BED-family:** `.`, `-1`, `-1`, then `.` for each further column,
///   except that the score of a 5- or 6-column B is `-1`.
/// - **GFF/GTF:** `.` for seqid, source and type, `-1` for start and end,
///   then `.`.
/// - **VCF/BCF:** `.` for the chromosome, `-1` for POS, then `.`.
///
/// This corrects D45's first wording, which had `.` for every column of a
/// VCF/GFF B and for the score of a BED5/BED6 B.
pub fn b_placeholder(b_format: Format, b_columns: usize, out: &mut Vec<u8>) {
    let mut fields: Vec<&str> = match b_format {
        // An alignment is printed as its BED6 row (D104).
        Format::Bed
        | Format::BedGraph
        | Format::Wig
        | Format::Bam
        | Format::Cram
        | Format::Sam => {
            let mut f = vec![".", "-1", "-1"];
            f.resize(b_columns.max(3), ".");
            if matches!(b_columns, 5 | 6) {
                f[4] = "-1";
            }
            f
        }
        Format::Gff3 | Format::Gtf => {
            let mut f = vec![".", ".", ".", "-1", "-1"];
            f.resize(b_columns.max(5), ".");
            f
        }
        Format::Vcf | Format::Bcf | Format::Fasta | Format::Psl | Format::Rmsk => {
            let mut f = vec![".", "-1"];
            f.resize(b_columns.max(2), ".");
            f
        }
    };
    if fields.is_empty() {
        fields.push(".");
    }
    out.extend_from_slice(fields.join("\t").as_bytes());
}
