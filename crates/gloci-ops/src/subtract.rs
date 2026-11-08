//! `subtract`: remove from A's records the bases covered by B.
//!
//! A reporting layer over the relation scaffold (`relate`), like
//! `intersect`. Records are written in A's input order (D38), in A's
//! format. For each A the hits that count (`Acceptor`) decide what goes:
//!
//! - **Default.** A's span minus the union of the accepted hits' windowed
//!   spans; each piece left is A with rewritten coordinates (D44), in
//!   ascending order. A with nothing removed is written as read. A point A
//!   has no bases, so it is dropped when it has an accepted hit.
//! - **`whole`.** A is dropped when it has an accepted hit.
//! - **`sum_overlaps`.** Only the strand rule (and window) decide which
//!   hits take part; A is dropped when their overlap lengths add up to at
//!   least `-f` of A, or when any takes part and there is no `-f` (D60).
//! - **VCF/BCF A.** A record is a variant, not a span to cut: it is
//!   dropped whole, exactly when `intersect --mode none` would drop it (or,
//!   with `--sum-overlaps -f`, when the summed overlaps reach `-f` of its
//!   footprint).
//! - **`--split`.** B is its BED12 blocks, not its span (D60).

use crate::relate::{Acceptor, RelateInputs, RelateOptions, rewrite_coords, run_relation};
use gloci_core::{Fraction, Interval, OverlapCriteria, Window, overlap_len};
use gloci_engine::{ARec, Hit, Neighbours, SweepSummary, Visitor};
use gloci_io::{Format, IoError, RecordWriter, UsageError};

#[derive(Debug, Clone)]
pub struct SubtractOptions {
    /// Drop A when it has an accepted hit, instead of cutting it.
    pub whole: bool,
    /// Drop A when the summed overlaps reach `-f` (D60).
    pub sum_overlaps: bool,
    pub relate: RelateOptions,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// The option checks that need no input, so a command can run them before
/// it creates any output.
pub fn validate(o: &SubtractOptions) -> Result<(), IoError> {
    if !o.sum_overlaps {
        return Ok(());
    }
    if o.whole {
        return Err(usage(
            "--sum-overlaps and --whole both drop A as a whole; give one of them",
        ));
    }
    let c = &o.relate.criteria;
    let given = [
        (c.min_frac_b.is_some(), "-F"),
        (c.reciprocal, "-r"),
        (c.either, "-e"),
        (c.min_bases.is_some(), "--min-bases"),
    ];
    if let Some((_, flag)) = given.iter().find(|(on, _)| *on) {
        return Err(usage(format!(
            "{flag} cannot be combined with --sum-overlaps, which takes only -f: it adds up the \
bases several B share with A; drop {flag} or use --whole"
        )));
    }
    if o.relate.window.is_some() && c.min_frac_a.is_some() {
        return Err(usage(
            "-w, --left and --right cannot be combined with -f: a fraction of a padded interval \
has no clear meaning; drop -f to remove A when any padded B hits it",
        ));
    }
    Ok(())
}

struct Subtractor<'o, 'w> {
    whole: bool,
    sum_overlaps: bool,
    window: Option<Window>,
    split: bool,
    acceptor: Acceptor<'o>,
    w: &'w mut dyn RecordWriter,
    format: Format,
    /// The `-f` of `--sum-overlaps`, taken out of the criteria the
    /// acceptor sees so that hits are chosen by strand alone.
    sum_frac: Option<Fraction>,
    removed: Vec<Interval>,
    row: Vec<u8>,
}

impl Subtractor<'_, '_> {
    fn put(&mut self) -> Result<(), IoError> {
        self.w
            .write_record(&self.row)
            .map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })
    }

    fn put_raw(&mut self, raw: &[u8]) -> Result<(), IoError> {
        self.row.clear();
        self.row.extend_from_slice(raw);
        self.put()
    }

    /// Add what `h` removes from `a`: its windowed span, or windowed blocks.
    fn add_removed(&mut self, a: &ARec, h: &Hit) -> Result<(), IoError> {
        let pad = |iv: Interval| match self.window {
            Some(w) => w.apply(iv, a.meta.strand),
            None => iv,
        };
        if self.split {
            for b in self.acceptor.hit_blocks(h)? {
                self.removed.push(pad(b));
            }
        } else {
            self.removed.push(pad(h.span));
        }
        Ok(())
    }

    /// The pieces of `span` left after `removed`, ascending, written as `a`
    /// with new coordinates. Nothing removed writes `a` as read.
    fn write_pieces(&mut self, a: &ARec) -> Result<(), IoError> {
        let span = a.meta.span;
        self.removed.sort_by_key(|r| (r.start, r.end));
        let mut pieces: Vec<Interval> = Vec::new();
        let mut at = span.start;
        for r in &self.removed {
            let (start, end) = (r.start.max(span.start), r.end.min(span.end));
            if end <= start {
                continue;
            }
            if start > at {
                pieces.push(Interval {
                    start: at,
                    end: start,
                });
            }
            at = at.max(end);
        }
        if at < span.end {
            pieces.push(Interval {
                start: at,
                end: span.end,
            });
        }
        if pieces.len() == 1 && pieces[0] == span {
            return self.put_raw(a.raw);
        }
        for p in pieces {
            self.row.clear();
            rewrite_coords(self.format, a.raw, p, &mut self.row);
            self.put()?;
        }
        Ok(())
    }
}

impl Visitor for Subtractor<'_, '_> {
    fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], _: &Neighbours<'_>) -> Result<(), IoError> {
        let span = a.meta.span;
        let variant = matches!(self.format, Format::Vcf | Format::Bcf);
        // Only the existence of an accepted hit matters, except for the
        // `--sum-overlaps -f` test, which every record kind goes through.
        let summing = self.sum_overlaps && self.sum_frac.is_some();
        let any_enough = self.whole
            || span.is_point()
            || (self.sum_overlaps && !summing)
            || (variant && !summing);
        self.removed.clear();
        let mut any = false;
        for h in hits {
            if !self.acceptor.accepted(a, h)? {
                continue;
            }
            any = true;
            if any_enough {
                break;
            }
            self.add_removed(a, h)?;
        }
        if !any {
            return self.put_raw(a.raw);
        }
        if any_enough {
            return Ok(());
        }
        if let Some(f) = self.sum_frac {
            let sum: u64 = self.removed.iter().map(|r| overlap_len(span, *r)).sum();
            return if f.met(sum, span.len()) {
                Ok(())
            } else {
                self.put_raw(a.raw)
            };
        }
        self.write_pieces(a)
    }
}

/// Run `subtract` over `inputs`, writing A's records to `w`. The writer
/// must be in A's format (`open_report` with `ReportOutput::SameAsA`).
/// Nothing is finished here; the caller calls `finish`.
pub fn subtract(
    inputs: RelateInputs,
    o: &SubtractOptions,
    w: &mut dyn RecordWriter,
) -> Result<SweepSummary, IoError> {
    validate(o)?;
    // `--sum-overlaps` picks hits by strand alone; its `-f` is applied to
    // their total.
    let mut relate = o.relate.clone();
    let mut sum_frac = None;
    if o.sum_overlaps {
        sum_frac = relate.criteria.min_frac_a.take();
        relate.criteria = OverlapCriteria {
            min_frac_a: None,
            ..relate.criteria
        };
    }
    let relate = &relate;
    let format = inputs.a.format();
    let acceptor = Acceptor::new(relate, &inputs);
    let mut s = Subtractor {
        whole: o.whole,
        sum_overlaps: o.sum_overlaps,
        window: relate.window,
        split: relate.split,
        acceptor,
        w,
        format,
        sum_frac,
        removed: Vec::new(),
        row: Vec::new(),
    };
    run_relation(inputs, relate, None, &mut s)
}
