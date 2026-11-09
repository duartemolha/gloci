//! `gloci shift`: move each feature along its chromosome (replaces
//! `bedtools shift`, D70). The engine, clipping, GFF rewriting and the order
//! warning are those of `slop`/`flank` (see [`crate::transform`]).
//!
//! The shifted interval is clipped like `slop` output: each coordinate is
//! held at 0 and, with a genome, at the chromosome length, so a feature
//! pushed over an edge shrinks. One pushed wholly over an edge would be
//! empty; it is dropped with a warning (bedtools writes a bogus interval).

use crate::transform::{Mode, run};
use gloci_core::{ChromTable, Fraction, Genome, Interval, Strand};
use gloci_io::{IoError, RecordReader, RecordWriter};

/// How far to move: signed bases, or a signed fraction of the feature length
/// (exact decimal, rounded toward the start of the chromosome like
/// bedtools: floor of the signed product).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShiftBy {
    Bases(i64),
    Fraction { negative: bool, fraction: Fraction },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Shift {
    pub by: ShiftBy,
    /// Negate the shift for `-` strand records; an unknown strand counts as `+`.
    pub strand_aware: bool,
}

impl Shift {
    /// The signed distance for a feature of `len` bases.
    fn delta(&self, len: u64, strand: Strand) -> i128 {
        let d = match self.by {
            ShiftBy::Bases(n) => i128::from(n),
            ShiftBy::Fraction {
                negative: false,
                fraction,
            } => i128::from(fraction.floor_mul(len)),
            ShiftBy::Fraction {
                negative: true,
                fraction,
            } => -i128::from(fraction.ceil_mul(len)),
        };
        if self.strand_aware && strand == Strand::Minus {
            -d
        } else {
            d
        }
    }

    /// `span` moved and clipped to `[0, limit]`.
    pub(crate) fn apply(&self, span: Interval, strand: Strand, limit: u64) -> Interval {
        let d = self.delta(span.len(), strand);
        let clip = |x: u64| (i128::from(x) + d).clamp(0, i128::from(limit)) as u64;
        Interval {
            start: clip(span.start),
            end: clip(span.end),
        }
    }
}

/// Write each record moved by `shift`, clipped at 0 and, with a genome, at
/// the chromosome length.
pub fn shift(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    shift: &Shift,
    genome: Option<&Genome>,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    run(Mode::Shift(shift), input, chroms, genome, w)
}
