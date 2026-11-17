//! Writing one record of a resolved query (D172).

use super::value::{end_of, write_float, write_tgt, write_typed, write_types};
use super::{Field, Piece, Query, SampleCols, Tag, TagType};
use crate::{format_cell, format_index, split_columns};

/// One sample of a `[ ]` block.
struct Sample<'a> {
    name: &'a [u8],
    cell: &'a [u8],
    /// Index of the piece's FORMAT key in the line's FORMAT column.
    key: Option<usize>,
}

impl Query {
    /// Append the output of one VCF text line (no newline; at least 8 columns).
    /// Err names the tag and, for a FORMAT value, the sample: `INFO/I value "abc" is not an Integer`.
    pub fn write_record(
        &self,
        line: &[u8],
        samples: &SampleCols<'_>,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        let mut cols = Vec::new();
        split_columns(line, &mut cols)?;
        let start = out.len();
        let result = self.write_pieces(&self.pieces, line, &cols, samples, out);
        if result.is_err() {
            out.truncate(start);
        }
        result
    }

    fn write_pieces(
        &self,
        pieces: &[Piece],
        line: &[u8],
        cols: &[&[u8]],
        samples: &SampleCols<'_>,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        for piece in pieces {
            match piece {
                Piece::Text(t) => out.extend_from_slice(t),
                Piece::Field(f) => write_field(f, line, cols, None, out)?,
                Piece::Block(inner) => write_block(inner, line, cols, samples, out)?,
            }
        }
        Ok(())
    }
}

fn write_block(
    inner: &[Piece],
    line: &[u8],
    cols: &[&[u8]],
    samples: &SampleCols<'_>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let format = cols.get(8).copied().unwrap_or(b".");
    let keys: Vec<Option<usize>> = inner
        .iter()
        .map(|p| match p {
            Piece::Field(Field::Gt | Field::Tgt) => format_index(format, b"GT"),
            Piece::Field(Field::FormatTag(t)) => format_index(format, &t.name),
            _ => None,
        })
        .collect();
    for &pos in samples.order {
        let header_index = samples.kept[pos];
        if let Some(pass) = samples.pass
            && !pass.get(header_index).copied().unwrap_or(true)
        {
            continue;
        }
        let cell = cols.get(9 + pos).copied().unwrap_or(b".");
        for (piece, &key) in inner.iter().zip(&keys) {
            match piece {
                Piece::Text(t) => out.extend_from_slice(t),
                Piece::Field(f) => {
                    let sample = Sample {
                        name: &samples.names[header_index],
                        cell,
                        key,
                    };
                    write_field(f, line, cols, Some(&sample), out)?;
                }
                // A block never nests (the syntax refuses it).
                Piece::Block(_) => {}
            }
        }
    }
    Ok(())
}

fn number(col: &[u8], what: &str) -> Result<i64, String> {
    std::str::from_utf8(col)
        .ok()
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| {
            format!(
                "{what} value \"{}\" is not an Integer",
                String::from_utf8_lossy(col)
            )
        })
}

fn push_i64(out: &mut Vec<u8>, v: i64) {
    out.extend_from_slice(itoa::Buffer::new().format(v).as_bytes());
}

fn write_field(
    f: &Field,
    line: &[u8],
    cols: &[&[u8]],
    sample: Option<&Sample<'_>>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    match f {
        Field::Chrom => out.extend_from_slice(cols[0]),
        Field::Pos => out.extend_from_slice(cols[1]),
        Field::Pos0 => push_i64(out, number(cols[1], "POS")? - 1),
        Field::End => push_i64(out, end_of(number(cols[1], "POS")?, cols[3], cols[7])),
        Field::End0 => push_i64(out, end_of(number(cols[1], "POS")?, cols[3], cols[7]) - 1),
        Field::Id => out.extend_from_slice(cols[2]),
        Field::Ref => out.extend_from_slice(cols[3]),
        Field::Alt(None) => out.extend_from_slice(cols[4]),
        Field::Alt(Some(n)) => match cols[4].split(|&b| b == b',').nth(*n) {
            Some(a) if !a.is_empty() => out.extend_from_slice(a),
            _ => out.push(b'.'),
        },
        Field::FirstAlt => {
            out.extend_from_slice(cols[4].split(|&b| b == b',').next().unwrap_or(b"."))
        }
        Field::Qual => {
            if cols[5] == b"." {
                out.push(b'.');
            } else {
                write_float(cols[5], out).ok_or_else(|| {
                    format!(
                        "QUAL value \"{}\" is not a Float",
                        String::from_utf8_lossy(cols[5])
                    )
                })?;
            }
        }
        Field::Filter => out.extend_from_slice(cols[6]),
        Field::Type => write_types(cols[3], cols[4], out),
        Field::Info => out.extend_from_slice(cols[7]),
        Field::InfoTag(tag) => write_info_tag(tag, cols[7], out)?,
        Field::Line => out.extend_from_slice(line),
        Field::Sample => out.extend_from_slice(sample.map_or(b".", |s| s.name)),
        Field::Gt => out.extend_from_slice(sample.map_or(b".", |s| format_cell(s.cell, s.key))),
        Field::Tgt => {
            if let Some(s) = sample {
                let gt = format_cell(s.cell, s.key);
                let gt = if gt == b"." { b"".as_slice() } else { gt };
                write_tgt(gt, cols[3], cols[4], out).map_err(|e| {
                    format!(
                        "FORMAT/GT of sample {}: {e}",
                        String::from_utf8_lossy(s.name)
                    )
                })?;
            }
        }
        Field::FormatTag(tag) => {
            if let Some(s) = sample {
                let value = format_cell(s.cell, s.key);
                write_typed(value, tag.ty, tag.sub, out).map_err(|bad| {
                    format!(
                        "FORMAT/{} value \"{bad}\" of sample {} is not {}",
                        String::from_utf8_lossy(&tag.name),
                        String::from_utf8_lossy(s.name),
                        article(tag.ty)
                    )
                })?;
            }
        }
    }
    Ok(())
}

fn article(ty: TagType) -> &'static str {
    match ty {
        TagType::Integer => "an Integer",
        _ => "a Float",
    }
}

fn write_info_tag(tag: &Tag, info: &[u8], out: &mut Vec<u8>) -> Result<(), String> {
    let mut found: Option<Option<&[u8]>> = None;
    if info != b"." {
        for item in info.split(|&b| b == b';') {
            let (name, value) = match item.iter().position(|&b| b == b'=') {
                Some(eq) => (&item[..eq], Some(&item[eq + 1..])),
                None => (item, None),
            };
            if name == tag.name.as_slice() {
                found = Some(value);
                break;
            }
        }
    }
    if tag.ty == TagType::Flag {
        out.push(if found.is_some() { b'1' } else { b'.' });
        return Ok(());
    }
    match found {
        Some(Some(v)) => write_typed(v, tag.ty, tag.sub, out).map_err(|bad| {
            format!(
                "INFO/{} value \"{bad}\" is not {}",
                String::from_utf8_lossy(&tag.name),
                article(tag.ty)
            )
        }),
        _ => {
            out.push(b'.');
            Ok(())
        }
    }
}
