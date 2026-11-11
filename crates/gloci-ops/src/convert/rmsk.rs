//! RepeatMasker `.out` to the BED of bedops `rmsk2bed`: sequence, begin-1,
//! end, repeat name, score, strand (`C` is `-`), divergence, deletions,
//! insertions, `(left)` of the query, class/family, the three repeat
//! positions exactly as written (parentheses kept), then the ID and any
//! `*` flag.

use super::{Sink, push_num, strand_byte};
use gloci_core::ChromTable;
use gloci_io::{IoError, RecordBatch, RecordReader, next_batch_size, rmsk_fields};

/// RepeatMasker columns (0-based) in the order `rmsk2bed` writes them after
/// the strand; the fields past column 13 (ID, `*`) follow when present.
const TAIL: [usize; 8] = [1, 2, 3, 7, 10, 11, 12, 13];

pub(super) fn convert(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    sink: &mut Sink,
) -> Result<(), IoError> {
    let mut batch = RecordBatch::default();
    let mut row = Vec::new();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let f: Vec<&[u8]> = rmsk_fields(batch.raw(i)).collect();
            let col = |n: usize| f.get(n).copied().unwrap_or(b".");
            row.clear();
            row.extend_from_slice(chroms.name(m.chrom));
            row.push(b'\t');
            push_num(&mut row, m.span.start);
            row.push(b'\t');
            push_num(&mut row, m.span.end);
            row.push(b'\t');
            row.extend_from_slice(col(9));
            row.push(b'\t');
            row.extend_from_slice(col(0));
            row.push(b'\t');
            row.push(strand_byte(m.strand));
            for n in TAIL {
                row.push(b'\t');
                row.extend_from_slice(col(n));
            }
            // The ID and the `*` flag, when the line has them.
            for extra in f.iter().skip(14) {
                row.push(b'\t');
                row.extend_from_slice(extra);
            }
            sink.row(m.chrom, 0, m.span.start, &row)?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}
