//! PSL to the 21-column BED of bedops `psl2bed`: tName, tStart, tEnd, qName,
//! matches, strand, then qSize, misMatches, repMatches, nCount, qNumInsert,
//! qBaseInsert, tNumInsert, tBaseInsert, qStart, qEnd, tSize, blockCount,
//! blockSizes, qStarts, tStarts.

use super::{Sink, push_num};
use gloci_core::ChromTable;
use gloci_io::{IoError, RecordBatch, RecordReader, field, next_batch_size};

/// PSL columns (0-based) in the order `psl2bed` writes them after the strand.
const TAIL: [usize; 15] = [10, 1, 2, 3, 4, 5, 6, 7, 11, 12, 14, 17, 18, 19, 20];

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
            let raw = batch.raw(i);
            let col = |n| field(raw, n).unwrap_or(b".");
            row.clear();
            // The chromosome as the reader named it (after --chrom-map).
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
            // As written, so a translated search keeps its two characters.
            row.extend_from_slice(col(8));
            for n in TAIL {
                row.push(b'\t');
                row.extend_from_slice(col(n));
            }
            // BLAT writes hits in query order: any descent is worth a warning.
            sink.row(m.chrom, 0, m.span.start, &row)?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}
