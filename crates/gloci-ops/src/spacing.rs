//! `spacing`: write each record with its distance to the previous record
//! (replaces `bedtools spacing`).
//!
//! The previous record is the one before it in the file on the same
//! chromosome (not the one that reaches furthest). The distance is the bases
//! between them, `0` when bookended, `-1` when the footprints overlap, and `.`
//! for the first record of a chromosome.

use crate::cluster::validate;
use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable};
use gloci_engine::kmerge_with;
use gloci_io::{IoError, RecordReader, RecordWriter};

/// The text of the spacing column for a record at `start` after one ending at
/// `prev_end` (`None` for the first record of a chromosome).
pub fn spacing_text(prev_end: Option<u64>, start: u64) -> String {
    match prev_end {
        None => ".".into(),
        Some(e) if start < e => "-1".into(),
        Some(e) => (start - e).to_string(),
    }
}

/// Write every record of the sorted `input` followed by a tab and its spacing.
pub fn spacing(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    validate(input.format())?;
    let mut skipped = Skipped::default();
    let input = skipped.watch(input);
    let mut chrom: Option<ChromId> = None;
    let mut prev_end: Option<u64> = None;
    let mut row = Vec::new();
    kmerge_with(vec![input], chroms, &mut |_, m, raw, _| {
        if chrom != Some(m.chrom) {
            chrom = Some(m.chrom);
            prev_end = None;
        }
        row.clear();
        row.extend_from_slice(raw);
        row.push(b'\t');
        row.extend_from_slice(spacing_text(prev_end, m.span.start).as_bytes());
        prev_end = Some(m.span.end);
        w.write_record(&row).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })
    })
    .map_err(|e| e.with_sorted_input_loc())?;
    log_skipped(&skipped.take());
    Ok(())
}
