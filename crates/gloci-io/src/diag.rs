//! Policy and bookkeeping for invalid records (`--skip-invalid`).

use crate::error::Loc;
use gloci_core::thousands;

/// What to do when a record cannot be parsed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InvalidPolicy {
    /// The first invalid record is an error.
    #[default]
    Fail,
    /// Count invalid records, keep a few examples, continue.
    Skip,
}

/// Maximum number of example messages kept.
const MAX_EXAMPLES: usize = 10;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Diagnostics {
    pub skipped: u64,
    /// At most 10 entries, each `"line N: msg"` (`"record N: msg"` for BCF).
    pub examples: Vec<String>,
}

impl Diagnostics {
    pub fn record(&mut self, line: u64, msg: &str) {
        self.record_at(Loc::Line(line), msg);
    }

    pub fn record_at(&mut self, loc: Loc, msg: &str) {
        self.record_place(&loc.to_string(), msg);
    }

    /// Like `record_at`, for a place that is not a line or record number.
    pub fn record_place(&mut self, place: &str, msg: &str) {
        self.skipped += 1;
        if self.examples.len() < MAX_EXAMPLES {
            self.examples.push(format!("{place}: {msg}"));
        }
    }

    /// A multi-line warning, or `None` when nothing was skipped.
    pub fn summary(&self, name: &str) -> Option<String> {
        if self.skipped == 0 {
            return None;
        }
        Some(format!(
            "skipped {} invalid records in {name}; first {}:\n  {}",
            thousands(self.skipped),
            self.examples.len(),
            self.examples.join("\n  ")
        ))
    }
}
