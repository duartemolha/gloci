//! GFF3/GTF to the 10-column BED of bedops `gff2bed`: chrom, start, end,
//! name, score, strand, source, type, phase, attributes.

use super::{ConvertOptions, Sink, push_num};
use gloci_core::ChromTable;
use gloci_io::{Format, IoError, RecordBatch, RecordReader, field, next_batch_size};

/// The value of attribute `key` in column 9, or `None` when it is absent or
/// empty. GFF3 attributes are `key=value;...`; GTF ones are
/// `key "value"; ...` (a quoted value may hold `;`).
pub(super) fn attr_value<'a>(attrs: &'a [u8], key: &[u8], gtf: bool) -> Option<&'a [u8]> {
    let mut rest = attrs;
    while !rest.is_empty() {
        // One attribute, up to the `;` outside quotes.
        let mut quoted = false;
        let mut skip = false;
        let mut end = rest.len();
        for (i, &b) in rest.iter().enumerate() {
            if skip {
                skip = false;
                continue;
            }
            match b {
                // A backslash escapes the next byte in a GTF value; the value
                // is written as it is, backslash included.
                b'\\' if gtf => skip = true,
                b'"' if gtf => quoted = !quoted,
                b';' if !quoted => {
                    end = i;
                    break;
                }
                _ => {}
            }
        }
        let item = rest[..end].trim_ascii();
        rest = rest.get(end + 1..).unwrap_or(b"");
        let split = if gtf {
            item.iter().position(|b| b.is_ascii_whitespace())
        } else {
            item.iter().position(|&b| b == b'=')
        };
        let Some(i) = split else { continue };
        if &item[..i] != key {
            continue;
        }
        let v = item[i + 1..].trim_ascii();
        let v = if gtf {
            v.strip_prefix(b"\"")
                .and_then(|v| v.strip_suffix(b"\""))
                .unwrap_or(v)
        } else {
            v
        };
        return (!v.is_empty()).then_some(v);
    }
    None
}

pub(super) fn convert(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &ConvertOptions,
    sink: &mut Sink,
) -> Result<(), IoError> {
    let gtf = input.format() == Format::Gtf;
    let key = o
        .name_attr
        .as_deref()
        .unwrap_or(if gtf { "gene_id" } else { "ID" })
        .as_bytes();
    let mut batch = RecordBatch::default();
    let mut row = Vec::new();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let raw = batch.raw(i);
            let col = |n| field(raw, n).unwrap_or(b".");
            if !o.features.is_empty() && !o.features.iter().any(|f| f.as_bytes() == col(2)) {
                continue;
            }
            row.clear();
            row.extend_from_slice(col(0));
            row.push(b'\t');
            push_num(&mut row, m.span.start);
            row.push(b'\t');
            push_num(&mut row, m.span.end);
            row.push(b'\t');
            row.extend_from_slice(attr_value(col(8), key, gtf).unwrap_or(b"."));
            // score, strand, source, type, phase, attributes
            for n in [5, 6, 1, 2, 7, 8] {
                row.push(b'\t');
                row.extend_from_slice(col(n));
            }
            sink.row(m.chrom, m.key, m.span.start, &row)?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::attr_value;

    #[test]
    fn gff3_attributes() {
        let a = b"ID=g1;Name=Foo bar;Empty=;Parent=a,b";
        assert_eq!(attr_value(a, b"ID", false), Some(&b"g1"[..]));
        assert_eq!(attr_value(a, b"Name", false), Some(&b"Foo bar"[..]));
        assert_eq!(attr_value(a, b"Empty", false), None);
        assert_eq!(attr_value(a, b"Parent", false), Some(&b"a,b"[..]));
        assert_eq!(attr_value(a, b"Missing", false), None);
        // A prefix of a key is not the key.
        assert_eq!(attr_value(b"IDx=1; ID=2", b"ID", false), Some(&b"2"[..]));
    }

    #[test]
    fn gtf_attributes() {
        let a = b"gene_id \"g1\"; gene_name \"A;B\"; tag basic; empty \"\";";
        assert_eq!(attr_value(a, b"gene_id", true), Some(&b"g1"[..]));
        assert_eq!(attr_value(a, b"gene_name", true), Some(&b"A;B"[..]));
        assert_eq!(attr_value(a, b"tag", true), Some(&b"basic"[..]));
        assert_eq!(attr_value(a, b"empty", true), None);
        assert_eq!(attr_value(a, b"gene", true), None);
        // An escaped quote does not end the value.
        let q = br#"gene_name "A\"B;C"; gene_id "g""#;
        assert_eq!(attr_value(q, b"gene_name", true), Some(&br#"A\"B;C"#[..]));
        assert_eq!(attr_value(q, b"gene_id", true), Some(&b"g"[..]));
    }
}
