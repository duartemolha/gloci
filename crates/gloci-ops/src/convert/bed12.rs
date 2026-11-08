//! BED12 to BED6: one row per block, keeping name, score and strand.

use super::{Sink, push_num};
use gloci_core::ChromTable;
use gloci_io::{
    IoError, Loc, RecordBatch, RecordReader, bed12_blocks, column_count, field, next_batch_size,
};

pub(super) fn convert(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    sink: &mut Sink,
) -> Result<(), IoError> {
    let name = input.name().to_string();
    let mut batch = RecordBatch::default();
    let mut row = Vec::new();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let raw = batch.raw(i);
            let fail = |msg: String| IoError::Parse {
                path: name.clone(),
                loc: Loc::Line(m.loc),
                msg,
            };
            let blocks = bed12_blocks(raw).map_err(|e| fail(format!("invalid BED12: {e}")))?;
            let Some(blocks) = blocks else {
                return Err(fail(format!(
                    "--blocks needs BED12 records with 12 columns; this one has {}",
                    column_count(raw)
                )));
            };
            let col = |n| field(raw, n).unwrap_or(b".");
            for b in blocks {
                row.clear();
                row.extend_from_slice(col(0));
                row.push(b'\t');
                push_num(&mut row, b.start);
                row.push(b'\t');
                push_num(&mut row, b.end);
                for n in 3..6 {
                    row.push(b'\t');
                    row.extend_from_slice(col(n));
                }
                sink.row(m.chrom, m.key, b.start, &row)?;
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}
