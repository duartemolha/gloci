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
    /// Blank lines a text reader dropped, and the first 10 line numbers.
    pub blank_lines: u64,
    pub blank_examples: Vec<u64>,
    /// Unmapped alignments an interval reader left out (D101).
    pub unmapped: u64,
    /// Alignments the read filters dropped (D102).
    pub filtered: u64,
    /// Mapped reads that start at or past the end of their reference, so cover no base of it.
    pub off_reference: u64,
    /// The one-time warnings of an interval reader that was told to defer them (a chunk of a
    /// per-chromosome run): the first of each kind, for the run to print once.
    pub deferred_off_reference: Option<String>,
    /// As `deferred_off_reference`, for a read cut at the end of its reference.
    pub deferred_clipped: Option<String>,
    /// `(example index, line, msg)` for the examples recorded by line, so
    /// `shift_lines` can renumber them.
    lines: Vec<(usize, u64, String)>,
}

impl Diagnostics {
    /// Note a blank line that was dropped.
    pub fn record_blank(&mut self, line: u64) {
        self.blank_lines += 1;
        if self.blank_examples.len() < MAX_EXAMPLES {
            self.blank_examples.push(line);
        }
    }

    pub fn record(&mut self, line: u64, msg: &str) {
        if self.examples.len() < MAX_EXAMPLES {
            self.lines
                .push((self.examples.len(), line, msg.to_string()));
        }
        self.record_at(Loc::Line(line), msg);
    }

    /// Add `by` to the line numbers of the examples recorded by line: for a
    /// reader that started `by` lines into its file.
    pub fn shift_lines(&mut self, by: u64) {
        for (i, line, msg) in &mut self.lines {
            *line += by;
            self.examples[*i] = format!("{}: {msg}", Loc::Line(*line));
        }
        for l in &mut self.blank_examples {
            *l += by;
        }
    }

    /// Add the counts and examples of `other` (keeping at most 10 examples).
    pub fn absorb(&mut self, other: Diagnostics) {
        self.skipped += other.skipped;
        self.blank_lines += other.blank_lines;
        self.unmapped += other.unmapped;
        self.filtered += other.filtered;
        self.off_reference += other.off_reference;
        // Chunks are absorbed in file order, so the first warning is the one streaming gives.
        if self.deferred_off_reference.is_none() {
            self.deferred_off_reference = other.deferred_off_reference;
        }
        if self.deferred_clipped.is_none() {
            self.deferred_clipped = other.deferred_clipped;
        }
        let room = MAX_EXAMPLES.saturating_sub(self.blank_examples.len());
        self.blank_examples
            .extend(other.blank_examples.into_iter().take(room));
        let base = self.examples.len();
        let room = MAX_EXAMPLES.saturating_sub(base);
        self.lines.extend(
            other
                .lines
                .into_iter()
                .filter(|(i, _, _)| *i < room)
                .map(|(i, line, msg)| (base + i, line, msg)),
        );
        self.examples.extend(other.examples.into_iter().take(room));
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

    /// Say at `-v` what an alignment reader left out of `name` (D101, D102): unmapped reads,
    /// reads the filters dropped, reads off the end of their reference and blank lines. Every
    /// command that reads alignments calls this when it is done with an input.
    pub fn log_alignment_counts(&self, name: &str) {
        self.log_alignment_counts_among(name, "");
    }

    /// As [`Diagnostics::log_alignment_counts`], with `note` after each line (for an input
    /// read only where it was needed: the counts are of the reads read).
    pub fn log_alignment_counts_among(&self, name: &str, note: &str) {
        let line = |n: u64, what: &str| {
            if n > 0 {
                log::info!("{name}: {} {what}{note}", thousands(n));
            }
        };
        line(self.unmapped, "unmapped reads left out");
        line(
            self.filtered,
            "reads dropped by --min-mapq, --require-flags or --exclude-flags",
        );
        line(
            self.off_reference,
            "reads left out for starting at or past the end of their reference",
        );
        line(self.blank_lines, "blank lines skipped");
    }

    /// Print the one-time warnings a deferring reader kept (see `deferred_off_reference`),
    /// the cut read's first. (Streaming prints them as met; a run that met them in the other
    /// order prints the same two lines, swapped.)
    pub fn log_deferred_warnings(&self) {
        for w in [&self.deferred_clipped, &self.deferred_off_reference]
            .into_iter()
            .flatten()
        {
            log::warn!("{w}");
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shift_and_absorb() {
        let mut a = Diagnostics::default();
        a.record(3, "bad a");
        let mut b = Diagnostics::default();
        b.record_place("record at BGZF offset 1:2", "bad b");
        b.record(1, "bad c");
        b.shift_lines(1_000);
        assert_eq!(
            b.examples,
            ["record at BGZF offset 1:2: bad b", "line 1,001: bad c"]
        );
        a.absorb(b);
        assert_eq!(a.skipped, 3);
        a.shift_lines(1);
        assert_eq!(
            a.examples,
            [
                "line 4: bad a",
                "record at BGZF offset 1:2: bad b",
                "line 1,002: bad c"
            ]
        );
        let mut many = Diagnostics::default();
        for i in 0..12 {
            many.record(i, "x");
        }
        let mut c = Diagnostics::default();
        c.record(1, "y");
        c.absorb(many);
        assert_eq!((c.skipped, c.examples.len()), (13, 10));
        c.shift_lines(10);
        assert_eq!(c.examples[0], "line 11: y");
        assert_eq!(c.examples[9], "line 18: x");
    }
}
