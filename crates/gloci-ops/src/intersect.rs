//! `intersect`: report A's records by how they relate to the B inputs.
//!
//! A reporting layer over the relation scaffold (`relate`): the engine
//! yields candidate hits (D36), the `Acceptor` decides which count, and
//! this module writes them per `--mode`. Records are written in A's input
//! order, hits by B file and then B input order (D38).

use crate::relate::{
    Acceptor, RelateInputs, RelateOptions, b_placeholder, rewrite_coords, run_relation,
};
use gloci_core::{ChromTable, Interval, StrandRule, Window, trimmed_alleles};
use gloci_engine::{ARec, Hit, Neighbours, SweepSummary, Visitor};
use gloci_io::{
    BcfText, Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordReader,
    RecordWriter, UsageError, column_count, field,
};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

/// What `intersect` writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// A's records with at least one accepted hit, as read.
    Any,
    /// A's records with none, as read.
    None,
    /// One TSV row per (A, accepted B) pair; with `keep_unmatched`, one row
    /// with a placeholder B for each A without a hit (D45).
    Pairs {
        overlap_bases: bool,
        keep_unmatched: bool,
    },
    /// A's record once per accepted hit, cut to the shared span (D44).
    Clip,
    /// A's record with its number of accepted hits (all B files).
    Count,
}

/// What makes a hit count besides the overlap criteria.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Match {
    Overlap,
    /// Some trimmed ALT of A equals some trimmed ALT of B, at the same
    /// trimmed position with the same trimmed REF (VCF/BCF inputs only).
    Alleles,
}

#[derive(Debug, Clone)]
pub struct IntersectOptions {
    pub mode: Mode,
    pub matching: Match,
    pub relate: RelateOptions,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn is_variant(f: Format) -> bool {
    matches!(f, Format::Vcf | Format::Bcf)
}

/// The option checks that need only the inputs' formats, so a command can
/// run them before it creates any output.
pub fn validate(a: Format, b: &[Format], o: &IntersectOptions) -> Result<(), IoError> {
    if o.mode == Mode::Clip {
        if is_variant(a) {
            return Err(usage(
                "a variant cannot be clipped; use --mode any to keep the whole record",
            ));
        }
        if o.relate.window.is_some() {
            return Err(usage(
                "--mode clip cannot be combined with -w, --left or --right: a window pads B, \
so there is no shared span to clip to; use --mode pairs",
            ));
        }
    }
    if o.matching == Match::Alleles {
        if !(is_variant(a) && b.iter().all(|f| is_variant(*f))) {
            return Err(usage(
                "--match alleles compares ALT alleles, so -a and every -b must be VCF or BCF files; \
use --match overlap for other formats",
            ));
        }
        let r = &o.relate;
        let c = &r.criteria;
        let given = [
            (c.min_frac_a.is_some(), "-f"),
            (c.min_frac_b.is_some(), "-F"),
            (c.reciprocal, "-r"),
            (c.either, "-e"),
            (c.min_bases.is_some(), "--min-bases"),
            (c.strand != StrandRule::Any, "-s/-S"),
            (c.inclusive_points, "--inclusive-points"),
            (r.split, "--split"),
            (r.window.is_some(), "-w/--left/--right"),
        ];
        if let Some((_, flag)) = given.iter().find(|(on, _)| *on) {
            return Err(usage(format!(
                "{flag} cannot be combined with --match alleles: allele matching compares the \
alleles of two records, not how much they overlap; drop {flag} or use --match overlap"
            )));
        }
    }
    Ok(())
}

/// Passes records through and notes the column count of the first one.
struct ColumnWatch {
    inner: Box<dyn RecordReader>,
    columns: Arc<AtomicUsize>,
}

impl RecordReader for ColumnWatch {
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
        let before = batch.len();
        let more = self.inner.read_batch(chroms, batch, max)?;
        if self.columns.load(Ordering::Relaxed) == 0 && batch.len() > before {
            self.columns
                .store(column_count(batch.raw(before)), Ordering::Relaxed);
        }
        Ok(more)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
}

/// How to print one input's records as text.
struct Text {
    format: Format,
    bcf: Option<BcfText>,
}

impl Text {
    fn new(r: &dyn RecordReader) -> Result<Text, IoError> {
        Ok(Text {
            format: r.format(),
            bcf: BcfText::for_reader(r)?,
        })
    }

    /// The record as a text line (a BCF record is rendered as VCF).
    fn of<'a>(&'a mut self, raw: &'a [u8], loc: u64) -> Result<&'a [u8], IoError> {
        match &mut self.bcf {
            Some(b) => b.render(raw, loc),
            None => Ok(raw),
        }
    }
}

/// Calls `f(pos, ref, alt)` with the trimmed (0-based position, REF, ALT)
/// of each usable ALT of a VCF line until it returns true. Symbolic,
/// breakend, `*` and missing ALTs have no sequence to compare and are
/// skipped. Allocates nothing; compare the slices case-insensitively.
fn any_allele(line: &[u8], mut f: impl FnMut(u64, &[u8], &[u8]) -> bool) -> bool {
    let pos = field(line, 1)
        .and_then(|p| std::str::from_utf8(p).ok())
        .and_then(|p| p.parse::<u64>().ok());
    let (Some(pos), Some(reference), Some(alt)) = (pos, field(line, 3), field(line, 4)) else {
        return false;
    };
    alt.split(|&c| c == b',')
        .filter(|a| {
            !a.is_empty()
                && !a
                    .iter()
                    .any(|c| matches!(c, b'<' | b'[' | b']' | b'*' | b'.'))
        })
        .any(|a| {
            let (p, r, a) = trimmed_alleles(pos.saturating_sub(1), reference, a);
            f(p, r, a)
        })
}

/// Do two VCF lines share a trimmed (position, REF, ALT)?
fn share_allele(a: &[u8], b: &[u8]) -> bool {
    any_allele(b, |bp, br, ba| {
        any_allele(a, |ap, ar, aa| {
            ap == bp && ar.eq_ignore_ascii_case(br) && aa.eq_ignore_ascii_case(ba)
        })
    })
}

struct Reporter<'o, 'w> {
    mode: Mode,
    matching: Match,
    acceptor: Acceptor<'o>,
    w: &'w mut dyn RecordWriter,
    a: Text,
    b: Vec<Text>,
    /// One label per B file when there are several, else empty.
    labels: &'o [String],
    /// Columns of a missing B (D45): known up front for VCF/GFF, else read
    /// from B's first record.
    b0_columns: Columns,
    row: Vec<u8>,
    /// A's VCF line, kept for the allele comparison.
    a_line: Vec<u8>,
}

enum Columns {
    Fixed(usize),
    FirstRecord(Arc<AtomicUsize>),
}

impl Reporter<'_, '_> {
    /// Is `h` a hit that counts for `a`?
    fn counts(&mut self, a: &ARec, h: &Hit) -> Result<bool, IoError> {
        if self.matching == Match::Overlap {
            return self.acceptor.accepted(a, h);
        }
        // Allele matching: the shared allele is the whole test (the
        // overlap criteria are refused by `validate`).
        let line = self.b[usize::from(h.file)].of(h.raw, h.loc)?;
        Ok(share_allele(&self.a_line, line))
    }

    fn put(&mut self) -> Result<(), IoError> {
        self.w
            .write_record(&self.row)
            .map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })
    }

    /// `row` = A as text.
    fn start_row(&mut self, a: &ARec) -> Result<(), IoError> {
        self.row.clear();
        let line = self.a.of(a.raw, a.meta.loc)?;
        self.row.extend_from_slice(line);
        Ok(())
    }
}

impl Visitor for Reporter<'_, '_> {
    fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], _: &Neighbours<'_>) -> Result<(), IoError> {
        if self.matching == Match::Alleles {
            let line = self.a.of(a.raw, a.meta.loc)?;
            self.a_line.clear();
            self.a_line.extend_from_slice(line);
        }
        match self.mode {
            Mode::Any | Mode::None => {
                let mut found = false;
                for h in hits {
                    if self.counts(a, h)? {
                        found = true;
                        break;
                    }
                }
                if found == (self.mode == Mode::Any) {
                    self.row.clear();
                    self.row.extend_from_slice(a.raw);
                    self.put()?;
                }
            }
            Mode::Count => {
                let mut n = 0u64;
                for h in hits {
                    n += u64::from(self.counts(a, h)?);
                }
                self.start_row(a)?;
                self.row.extend_from_slice(format!("\t{n}").as_bytes());
                self.put()?;
            }
            Mode::Clip => {
                for h in hits {
                    if !self.counts(a, h)? {
                        continue;
                    }
                    // Overlapping spans share [max start, min end), a point
                    // when one side is a point.
                    let shared = Interval {
                        start: a.meta.span.start.max(h.span.start),
                        end: a.meta.span.end.min(h.span.end),
                    };
                    self.row.clear();
                    rewrite_coords(self.a.format, a.raw, shared, &mut self.row);
                    self.put()?;
                }
            }
            Mode::Pairs {
                overlap_bases,
                keep_unmatched,
            } => {
                let mut any = false;
                for h in hits {
                    if !self.counts(a, h)? {
                        continue;
                    }
                    any = true;
                    self.start_row(a)?;
                    if !self.labels.is_empty() {
                        self.row.push(b'\t');
                        self.row
                            .extend_from_slice(self.labels[usize::from(h.file)].as_bytes());
                    }
                    self.row.push(b'\t');
                    let line = self.b[usize::from(h.file)].of(h.raw, h.loc)?;
                    self.row.extend_from_slice(line);
                    if overlap_bases {
                        let n = self.acceptor.overlap_len(a, h)?;
                        self.row.extend_from_slice(format!("\t{n}").as_bytes());
                    }
                    self.put()?;
                }
                if !any && keep_unmatched {
                    self.start_row(a)?;
                    if !self.labels.is_empty() {
                        self.row.extend_from_slice(b"\t.");
                    }
                    self.row.push(b'\t');
                    let columns = match &self.b0_columns {
                        Columns::Fixed(n) => *n,
                        // B is empty or unread: the shortest BED line
                        Columns::FirstRecord(c) => match c.load(Ordering::Relaxed) {
                            0 => 3,
                            n => n,
                        },
                    };
                    b_placeholder(self.b[0].format, columns, &mut self.row);
                    if overlap_bases {
                        self.row.extend_from_slice(b"\t0");
                    }
                    self.put()?;
                }
            }
        }
        Ok(())
    }
}

/// Run `intersect` over `inputs`, writing to `w`. The writer must already
/// match the mode: A's format for `any`, `none` and `clip`, TSV otherwise
/// (`open_report`). Nothing is finished here; the caller calls `finish`.
pub fn intersect(
    mut inputs: RelateInputs,
    o: &IntersectOptions,
    w: &mut dyn RecordWriter,
) -> Result<SweepSummary, IoError> {
    let b_formats: Vec<Format> = inputs.b.iter().map(|r| r.format()).collect();
    validate(inputs.a.format(), &b_formats, o)?;
    let multi = inputs.b.len() > 1;
    if multi && matches!(o.mode, Mode::Pairs { .. }) && o.relate.labels.len() != inputs.b.len() {
        return Err(usage(format!(
            "{} -b files but {} labels; give one label per -b",
            inputs.b.len(),
            o.relate.labels.len()
        )));
    }
    let a = Text::new(inputs.a.as_ref())?;
    let b = inputs
        .b
        .iter()
        .map(|r| Text::new(r.as_ref()))
        .collect::<Result<Vec<_>, _>>()?;
    let b0_columns = match inputs.b.first().map(|r| (r.format(), r.header())) {
        Some((Format::Bed | Format::BedGraph, _)) | None => {
            let seen = Arc::new(AtomicUsize::new(0));
            if matches!(
                o.mode,
                Mode::Pairs {
                    keep_unmatched: true,
                    ..
                }
            ) {
                let first = inputs.b.remove(0);
                inputs.b.insert(
                    0,
                    Box::new(ColumnWatch {
                        inner: first,
                        columns: seen.clone(),
                    }),
                );
            }
            Columns::FirstRecord(seen)
        }
        Some((Format::Gff3 | Format::Gtf, _)) => Columns::Fixed(9),
        Some((_, Header::Vcf(h))) => Columns::Fixed(if h.samples.is_empty() {
            8
        } else {
            9 + h.samples.len()
        }),
        Some(_) => Columns::Fixed(8),
    };
    // Candidates for allele matching: a shared trimmed allele lies at most
    // one base from each record's footprint (an insertion is a point just
    // past an anchor base, and a footprint can end where the other starts),
    // so pad B by one base on each side. The candidates are then accepted
    // by allele equality alone.
    let mut relate = o.relate.clone();
    if o.matching == Match::Alleles {
        relate.window = Some(Window {
            left: 1,
            right: 1,
            strand_aware: false,
        });
    }
    let relate = &relate;
    let acceptor = Acceptor::new(relate, &inputs);
    let mut reporter = Reporter {
        mode: o.mode,
        matching: o.matching,
        acceptor,
        w,
        a,
        b,
        labels: if multi && matches!(o.mode, Mode::Pairs { .. }) {
            &relate.labels
        } else {
            &[]
        },
        b0_columns,
        row: Vec::new(),
        a_line: Vec::new(),
    };
    run_relation(inputs, relate, None, &mut reporter)
}
