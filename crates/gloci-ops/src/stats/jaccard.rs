//! `stats jaccard`: the Jaccard similarity of two interval sets (replaces
//! `bedtools jaccard`).

use super::merged::{DropUnknown, MergedReader, Tally};
use crate::relate::{RelateInputs, RelateOptions, relation_chroms, run_relation_in};
use crate::skipped::{Skipped, log_skipped};
use crate::sort::{SortSettings, sort_readers};
use gloci_core::{ChromOrder, OverlapCriteria, StrandRule, overlap_len};
use gloci_engine::{ARec, Hit, Neighbours, Visitor};
use gloci_io::{IoError, UsageError};
use std::sync::atomic::Ordering;

#[derive(Debug, Clone)]
pub struct JaccardOptions {
    pub criteria: OverlapCriteria,
    /// `--sort`: sort A and B before merging them.
    pub sort: Option<SortSettings>,
    pub genome: Option<String>,
    pub chrom_order: ChromOrder,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JaccardResult {
    /// Bases shared by the merged intervals of A and B.
    pub intersection: u64,
    /// Bases covered by either.
    pub union: u64,
    /// Pairs of merged intervals that overlap and pass the criteria.
    pub n_intersections: u64,
}

impl JaccardResult {
    /// The statistic, divided in single precision as bedtools does.
    pub fn jaccard(&self) -> f32 {
        self.intersection as f32 / self.union as f32
    }
}

/// Sums the overlaps of merged A and merged B pairs.
struct Sum<'a> {
    criteria: &'a OverlapCriteria,
    intersection: u64,
    pairs: u64,
}

impl Visitor for Sum<'_> {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        _near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        for h in hits {
            if self
                .criteria
                .accepts(a.meta.span, a.meta.strand, h.span, h.strand)
            {
                self.intersection += overlap_len(a.meta.span, h.span);
                self.pairs += 1;
            }
        }
        Ok(())
    }
}

/// Merge A and B (overlapping and bookended intervals join; per strand with
/// `-s`/`-S`) and compare the merged sets: the intersection is the sum of
/// overlaps of the merged pairs the criteria accept, the union what is left
/// of the two totals. Both inputs are read to the end.
pub fn jaccard(inputs: RelateInputs, o: &JaccardOptions) -> Result<JaccardResult, IoError> {
    let RelateInputs { a, b } = inputs;
    if b.len() != 1 {
        return Err(UsageError(format!(
            "jaccard compares two files, A and one B, but {} B files were given; \
run it once per B file",
            b.len()
        ))
        .into());
    }
    let relate = RelateOptions {
        criteria: o.criteria,
        window: None,
        split: false,
        split_deletions: false,
        unsorted: false,
        sort: None,
        labels: Vec::new(),
        genome: o.genome.clone(),
        chrom_order: o.chrom_order,
        header: false,
    };
    let inputs = RelateInputs { a, b };
    let mut chroms = relation_chroms(&inputs, &relate)?;
    let RelateInputs { a, b } = inputs;
    let by_strand = o.criteria.strand != StrandRule::Any;
    let mut skipped = Skipped::default();
    let prepare = |r: Box<dyn gloci_io::RecordReader>,
                   chroms: &mut gloci_core::ChromTable,
                   skipped: &mut Skipped|
     -> Result<Box<dyn gloci_io::RecordReader>, IoError> {
        match &o.sort {
            Some(s) => sort_readers(vec![r], chroms, s).map_err(IoError::with_sorted_input_loc),
            None => Ok(skipped.watch(r)),
        }
    };
    let a = prepare(a, &mut chroms, &mut skipped)?;
    let b = prepare(
        b.into_iter().next().expect("one B"),
        &mut chroms,
        &mut skipped,
    )?;
    let (name_a, name_b) = (a.name().to_string(), b.name().to_string());
    // Under -s/-S a cluster of unknown strand can never match, so it takes no
    // part in the tally (bedtools leaves such records out too).
    let wrap = |r: Box<dyn gloci_io::RecordReader>| -> (Box<dyn gloci_io::RecordReader>, _) {
        let merged = Box::new(MergedReader::new(r, by_strand));
        if by_strand {
            let (d, flag) = DropUnknown::new(merged);
            (Box::new(d), Some(flag))
        } else {
            (merged, None)
        }
    };
    let (a, flag_a) = wrap(a);
    let (b, flag_b) = wrap(b);
    let (a, a_counts) = Tally::new(a);
    let (b, b_counts) = Tally::new(b);
    let mut sum = Sum {
        criteria: &o.criteria,
        intersection: 0,
        pairs: 0,
    };
    run_relation_in(
        RelateInputs {
            a: Box::new(a),
            b: vec![Box::new(b)],
        },
        Some(chroms),
        &relate,
        None,
        &mut sum,
    )?;
    for (name, flag) in [(name_a, flag_a), (name_b, flag_b)] {
        if flag.is_some_and(|f| f.load(Ordering::Relaxed)) {
            log::warn!(
                "{name} has records with unknown strand; under -s and -S they are left out \
of the statistic (they can never match)"
            );
        }
    }
    log_skipped(&skipped.take());
    skipped.log_alignment_counts();
    let total = a_counts
        .bases
        .load(Ordering::Relaxed)
        .saturating_add(b_counts.bases.load(Ordering::Relaxed));
    Ok(JaccardResult {
        intersection: sum.intersection,
        union: total.saturating_sub(sum.intersection),
        n_intersections: sum.pairs,
    })
}

/// The header and the one row, as bedtools prints them.
pub fn write_jaccard(r: &JaccardResult, out: &mut Vec<u8>) {
    out.extend_from_slice(b"intersection\tunion\tjaccard\tn_intersections\n");
    out.extend_from_slice(format!("{}\t{}\t", r.intersection, r.union).as_bytes());
    gloci_core::cfmt::write_g6(out, f64::from(r.jaccard()));
    out.extend_from_slice(format!("\t{}\n", r.n_intersections).as_bytes());
}
