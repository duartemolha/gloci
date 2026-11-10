//! BAM, CRAM and SAM to BED (D106, bedtools `bamtobed`): input order is kept, unmapped reads
//! are skipped by the reader.

use super::{ConvertOptions, Sink};
use gloci_io::{AlnBed, AlnBedLayout, AlnBedOptions, IoError, RecordBatch, RecordReader, field};

impl ConvertOptions {
    /// What the alignment flags ask `AlnBed` for.
    fn aln_bed(&self) -> AlnBedOptions {
        let blocks = self.blocks || self.split_deletions;
        AlnBedOptions {
            layout: if self.bed12 {
                AlnBedLayout::Bed12
            } else if blocks {
                AlnBedLayout::Blocks
            } else {
                AlnBedLayout::Bed6
            },
            score_tag: self.score_tag,
            cigar: self.cigar,
            split_deletions: self.split_deletions,
            color: self.color.clone().unwrap_or_else(|| "255,0,0".into()),
        }
    }
}

pub(super) fn convert(
    mut input: Box<dyn RecordReader>,
    chroms: &mut gloci_core::ChromTable,
    o: &ConvertOptions,
    sink: &mut Sink,
) -> Result<(), IoError> {
    let name = input.name().to_string();
    let Some(mut bed) = AlnBed::for_reader(&*input, o.aln_bed()) else {
        return Err(IoError::Internal(format!(
            "{name} is not an alignment reader"
        )));
    };
    if o.blocks || o.split_deletions {
        // A read's blocks can start after the next read does: the rows of a split sorted BAM
        // are out of order by nature, which the order warning would blame on the input.
        sink.warned = true;
    }
    let mut batch = RecordBatch::default();
    let mut rows = Vec::new();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            rows.clear();
            bed.write(batch.raw(i), m.loc, &mut rows)?;
            for row in rows.split(|&b| b == b'\n').filter(|r| !r.is_empty()) {
                let start = field(row, 1)
                    .and_then(|f| std::str::from_utf8(f).ok())
                    .and_then(|f| f.parse().ok())
                    .unwrap_or(m.span.start);
                sink.row(m.chrom, m.key, start, row)?;
            }
        }
        max = gloci_io::next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    let d = input.diagnostics();
    d.log_alignment_counts(&name);
    if let Some(s) = d.summary(&name) {
        log::warn!("{s}");
    }
    Ok(())
}
