//! The `-H` line of a resolved query, as bcftools 1.21 writes it (D172).

use super::{Field, Piece, Query};

fn name(f: &Field) -> &[u8] {
    match f {
        Field::Chrom => b"CHROM",
        Field::Pos => b"POS",
        Field::Pos0 => b"POS0",
        Field::End => b"END",
        Field::End0 => b"END0",
        Field::Id => b"ID",
        Field::Ref => b"REF",
        Field::Alt(_) => b"ALT",
        Field::FirstAlt => b"FIRST_ALT",
        Field::Qual => b"QUAL",
        Field::Filter => b"FILTER",
        Field::Type => b"TYPE",
        Field::Info => b"INFO",
        Field::Line => b"LINE",
        Field::Sample => b"SAMPLE",
        Field::Gt | Field::Tgt => b"GT",
        Field::InfoTag(t) | Field::FormatTag(t) => &t.name,
    }
}

fn holds_newline(pieces: &[Piece]) -> bool {
    pieces
        .iter()
        .any(|p| matches!(p, Piece::Text(t) if t.contains(&b'\n')))
}

impl Query {
    /// Append the -H line for printed samples `names`; `filter_reads_samples` drops the `SAMPLE:` prefixes (D172).
    pub fn write_header(&self, names: &[&[u8]], filter_reads_samples: bool, out: &mut Vec<u8>) {
        let newline_block = |p: &Piece| matches!(p, Piece::Block(b) if holds_newline(b));
        let any_newline = self.pieces.iter().any(newline_block);
        let prefix = !filter_reads_samples && !any_newline;
        out.push(b'#');
        let mut n = 0;
        let mut once = false;
        for piece in &self.pieces {
            match piece {
                Piece::Text(t) => out.extend_from_slice(t),
                Piece::Field(f) => {
                    n += 1;
                    out.extend_from_slice(format!("[{n}]").as_bytes());
                    out.extend_from_slice(name(f));
                }
                Piece::Block(inner) => {
                    once |= holds_newline(inner);
                    let shown = if once {
                        &names[..names.len().min(1)]
                    } else {
                        names
                    };
                    for sample in shown {
                        for p in inner {
                            match p {
                                Piece::Text(t) => {
                                    out.extend(t.iter().filter(|&&b| b != b'\n'));
                                }
                                Piece::Field(f) => {
                                    n += 1;
                                    out.extend_from_slice(format!("[{n}]").as_bytes());
                                    if prefix {
                                        out.extend_from_slice(sample);
                                        out.push(b':');
                                    }
                                    out.extend_from_slice(name(f));
                                }
                                Piece::Block(_) => {}
                            }
                        }
                    }
                }
            }
        }
        if any_newline {
            out.push(b'\n');
        }
    }
}
