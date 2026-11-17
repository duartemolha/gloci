//! WIG, bigWig and bigBed: one row per record of the reader (D63, D171), already bedGraph or
//! BED text, values as written.

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
            // The reader's records are already `chrom start end value` (bigBed: the BED line).
            sink.row(m.chrom, m.key, m.span.start, batch.raw(i))?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            crate::skipped::log_reader(input.as_ref());
            return Ok(());
        }
    }
}
