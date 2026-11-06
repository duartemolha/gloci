//! On-the-fly sort-order checking (spec §4.5).

use crate::chrom::{ChromId, ChromTable};
use std::cmp::Ordering;

/// A record that breaks the expected `(chrom, key)` order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrderViolation {
    /// The chromosome sorts before the previous record's chromosome
    /// (or reappears after the input has moved on).
    ChromBackwards { prev: ChromId, next: ChromId },
    /// Same chromosome, but the key is smaller than the previous key.
    StartBackwards {
        chrom: ChromId,
        prev: u64,
        next: u64,
    },
}

/// Tracks the previous `(chrom, key)` and reports the first backwards step.
#[derive(Debug, Clone, Default)]
pub struct OrderChecker {
    prev: Option<(ChromId, u64)>,
}

impl OrderChecker {
    pub fn new() -> OrderChecker {
        OrderChecker { prev: None }
    }

    /// Checks the next record. After a violation the new record becomes the
    /// previous one, so scanning can continue.
    pub fn check(
        &mut self,
        table: &ChromTable,
        chrom: ChromId,
        key: u64,
    ) -> Result<(), OrderViolation> {
        let result = match self.prev {
            None => Ok(()),
            Some((prev_chrom, prev_key)) => match table.cmp(chrom, prev_chrom) {
                Ordering::Less => Err(OrderViolation::ChromBackwards {
                    prev: prev_chrom,
                    next: chrom,
                }),
                Ordering::Equal if key < prev_key => Err(OrderViolation::StartBackwards {
                    chrom,
                    prev: prev_key,
                    next: key,
                }),
                _ => Ok(()),
            },
        };
        self.prev = Some((chrom, key));
        result
    }
}
