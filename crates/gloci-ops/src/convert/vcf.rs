//! VCF/BCF to BED: chrom, start, end, ID, REF, ALT, QUAL, FILTER.
//!
//! Without `--split-alleles` the footprint is the record's, under the
//! reader's `--vcf-span`. With it every row is one ALT's allele footprint
//! whatever the span (spec 4.2: one interval per ALT).

use super::{ConvertOptions, Sink, push_num};
use gloci_core::{ChromTable, Variant, VariantKind, alt_footprints, alt_kinds};
use gloci_io::{BcfText, IoError, Loc, RecordBatch, RecordReader, field, next_batch_size};

pub(super) fn convert(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &ConvertOptions,
    sink: &mut Sink,
) -> Result<(), IoError> {
    let name = input.name().to_string();
    let mut bcf = BcfText::for_reader(input.as_ref())?;
    let binary = bcf.is_some();
    let mut batch = RecordBatch::default();
    let mut row = Vec::new();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let line = match &mut bcf {
                Some(b) => b.render(batch.raw(i), m.loc)?,
                None => batch.raw(i),
            };
            let bad = |msg: String| IoError::Parse {
                path: name.clone(),
                loc: if binary {
                    Loc::Record(m.loc)
                } else {
                    Loc::Line(m.loc)
                },
                msg,
            };
            let col = |n| field(line, n).unwrap_or(b".");
            let (reference, alts) = (col(3), col(4));
            // The reader checked the record: POS is a number, END parses.
            let pos = std::str::from_utf8(col(1))
                .ok()
                .and_then(|p| p.parse::<u64>().ok())
                .ok_or_else(|| bad("POS is not an integer".into()))?;
            let (mut end, mut svlen) = (None, None);
            if alts.contains(&b'<') {
                for kv in col(7).split(|&b| b == b';') {
                    if end.is_none()
                        && let Some(v) = kv.strip_prefix(b"END=")
                    {
                        end = std::str::from_utf8(v)
                            .ok()
                            .and_then(|v| v.parse::<u64>().ok());
                    }
                    if svlen.is_none()
                        && let Some(v) = kv.strip_prefix(b"SVLEN=")
                    {
                        svlen = Some(v);
                    }
                }
            }
            let v = Variant {
                pos,
                reference,
                alts,
                end,
                svlen,
            };
            let wanted = |k: VariantKind| o.types.as_ref().is_none_or(|t| t.contains(&k));
            let mut write = |sink: &mut Sink, start: u64, end: u64, alt: &[u8]| {
                row.clear();
                row.extend_from_slice(col(0));
                row.push(b'\t');
                push_num(&mut row, start);
                row.push(b'\t');
                push_num(&mut row, end);
                for v in [col(2), reference, alt, col(5), col(6)] {
                    row.push(b'\t');
                    row.extend_from_slice(v);
                }
                sink.row(m.chrom, m.key, start, &row)
            };
            if !o.split_alleles {
                if o.types.is_some() && !alt_kinds(&v).into_iter().any(wanted) {
                    continue;
                }
                write(sink, m.span.start, m.span.end, alts)?;
                continue;
            }
            let kinds = alt_kinds(&v);
            let prints = alt_footprints(&v).map_err(|e| bad(e.to_string()))?;
            let alt_list: Vec<&[u8]> = alts.split(|&b| b == b',').collect();
            if prints.iter().all(Option::is_none) {
                // `*`, `<NON_REF>` or `.` alone: nothing to split, keep the
                // record whole.
                if kinds.iter().copied().any(wanted) {
                    write(sink, m.span.start, m.span.end, alts)?;
                }
                continue;
            }
            for ((fp, kind), alt) in prints.iter().zip(&kinds).zip(&alt_list) {
                let Some(fp) = fp else { continue };
                if !wanted(*kind) {
                    continue;
                }
                write(sink, fp.start, fp.end, alt)?;
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}
