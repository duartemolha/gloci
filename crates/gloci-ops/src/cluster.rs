//! `cluster`: write each record with the id of the cluster it belongs to
//! (replaces `bedtools cluster`).
//!
//! Clusters are those of `merge` (same gap and strand rules, shared with it),
//! but every input record is kept: the output is the input line, a tab and the
//! cluster id. Ids start at 1 and increase through the file, across
//! chromosomes; with strands, a new cluster takes the next id when its first
//! record is read.

use crate::merge::{joins, slot_of};
use crate::skipped::{Skipped, log_skipped};
use gloci_core::{ChromId, ChromTable, Strand};
use gloci_engine::kmerge_with;
use gloci_io::{Format, IoError, RecordReader, RecordWriter, UsageError};

#[derive(Debug, Clone)]
pub struct ClusterOptions {
    /// As `merge --max-gap`: a record joins the open cluster when
    /// `start <= cluster_end + max_gap` (signed).
    pub max_gap: i64,
    /// Cluster each strand separately (`+`, `-` and unknown).
    pub per_strand: bool,
}

/// Reject input formats whose records are not text lines.
pub fn validate(format: Format) -> Result<(), UsageError> {
    if format == Format::Bcf {
        return Err(UsageError(
            "cluster appends a column to each input line, which BCF records are not; convert \
             first: gloci view FILE -O vcf | gloci cluster -"
                .into(),
        ));
    }
    Ok(())
}

/// The open cluster of one strand slot.
#[derive(Clone, Copy)]
struct Open {
    end: u64,
    id: u64,
}

/// Write every record of the sorted `input` followed by a tab and its
/// cluster id.
pub fn cluster(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &ClusterOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    validate(input.format())?;
    let mut skipped = Skipped::default();
    let input = skipped.watch(input);
    let mut chrom: Option<ChromId> = None;
    let mut open: [Option<Open>; 3] = [None; 3];
    let mut next_id = 0u64;
    let mut unknown_strand = false;
    let mut row = Vec::new();
    kmerge_with(vec![input], chroms, &mut |_, m, raw, _| {
        if chrom != Some(m.chrom) {
            chrom = Some(m.chrom);
            open = [None; 3];
        }
        unknown_strand |= m.strand == Strand::Unknown;
        let slot = &mut open[slot_of(o.per_strand, m.strand)];
        match slot {
            Some(c) if joins(m.span.start, c.end, o.max_gap) => c.end = c.end.max(m.span.end),
            _ => {
                next_id += 1;
                *slot = Some(Open {
                    end: m.span.end,
                    id: next_id,
                });
            }
        }
        let id = slot.expect("set above").id;
        row.clear();
        row.extend_from_slice(raw);
        row.push(b'\t');
        row.extend_from_slice(id.to_string().as_bytes());
        w.write_record(&row).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })
    })
    .map_err(|e| e.with_sorted_input_loc())?;
    log_skipped(&skipped.take());
    if o.per_strand && unknown_strand {
        log::warn!(
            "some records have no strand; with -s they are clustered together, apart from + and - \
records"
        );
    }
    Ok(())
}
