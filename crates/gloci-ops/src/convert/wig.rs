//! WIG to bedGraph: one row per WIG data line (D63), values as written.

use super::Sink;
use gloci_core::ChromTable;
use gloci_io::{IoError, RecordBatch, RecordReader, next_batch_size};

pub(super) fn convert(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    sink: &mut Sink,
) -> Result<(), IoError> {
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            // The reader's records are already `chrom start end value`.
            sink.row(m.chrom, m.key, m.span.start, batch.raw(i))?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}
