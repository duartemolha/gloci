//! Reading ahead to an input's first record without reading anything twice.

use gloci_core::ChromTable;
use gloci_io::{
    Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordMeta, RecordReader,
    column_count,
};

/// A reader whose first record was read ahead and is handed back by the
/// first `read_batch`, so the input (a pipe, say) is read only once. The
/// look-ahead used the run's chromosome table, so `read_batch` must be given
/// that same table.
struct PeekReader {
    inner: Box<dyn RecordReader>,
    /// The batch read ahead and whether more input followed it.
    ahead: Option<(RecordBatch, bool)>,
}

impl RecordReader for PeekReader {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn compression(&self) -> Compression {
        self.inner.compression()
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        match self.ahead.take() {
            Some((ahead, more)) => {
                let base = batch.buf.len();
                batch.buf.extend_from_slice(&ahead.buf);
                batch.recs.extend(ahead.recs.iter().map(|m| RecordMeta {
                    off: m.off + base,
                    ..*m
                }));
                Ok(more)
            }
            None => self.inner.read_batch(chroms, batch, max),
        }
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
}

/// Read `r` up to its first record with `chroms` (the table the run will
/// use). Returns a reader that yields every record of `r`, that first one
/// included, and the column count of that record (`None` when `r` has no
/// record). Invalid records skipped on the way stay counted in `r`'s
/// diagnostics; an error is returned as reading would return it.
pub(crate) fn peek_first(
    mut r: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
) -> Result<(Box<dyn RecordReader>, Option<usize>), IoError> {
    let mut batch = RecordBatch::default();
    loop {
        let more = r.read_batch(chroms, &mut batch, 1)?;
        if !batch.is_empty() || !more {
            let columns = (!batch.is_empty()).then(|| column_count(batch.raw(0)));
            let reader = PeekReader {
                inner: r,
                ahead: Some((batch, more)),
            };
            return Ok((Box::new(reader), columns));
        }
    }
}
