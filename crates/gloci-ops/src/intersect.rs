//! `intersect`: report A's records by how they relate to the B inputs.
//!
//! A reporting layer over the relation scaffold (`relate`): the engine
//! yields candidate hits (D36), the `Acceptor` decides which count, and
//! this module writes them per `--mode`. Records are written in A's input
//! order, hits by B file and then B input order (D38).

use crate::relate::{
    Acceptor, RelateInputs, RelateOptions, b_placeholder, rewrite_coords, run_relation,
};
use gloci_core::{ChromTable, Interval, trimmed_alleles};
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
    if o.matching == Match::Alleles && !(is_variant(a) && b.iter().all(|f| is_variant(*f))) {
        return Err(usage(
            "--match alleles compares ALT alleles, so -a and every -b must be VCF or BCF files; \
use --match overlap for other formats",
        ));
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

/// POS (0-based) and the trimmed (position, REF, ALT) of each usable ALT of
/// a VCF line. Symbolic, breakend, `*` and missing ALTs have no sequence to
/// compare and are skipped.
fn alleles(line: &[u8]) -> Vec<(u64, Vec<u8>, Vec<u8>)> {
    let pos = field(line, 1)
        .and_then(|p| std::str::from_utf8(p).ok())
        .and_then(|p| p.parse::<u64>().ok());
    let (Some(pos), Some(reference), Some(alt)) = (pos, field(line, 3), field(line, 4)) else {
        return Vec::new();
    };
    alt.split(|&c| c == b',')
        .filter(|a| {
            !a.is_empty()
                && !a
                    .iter()
                    .any(|c| matches!(c, b'<' | b'[' | b']' | b'*' | b'.'))
        })
        .map(|a| {
            let (p, r, a) = trimmed_alleles(pos.saturating_sub(1), reference, a);
            (p, r.to_ascii_uppercase(), a.to_ascii_uppercase())
        })
        .collect()
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
    a_alleles: Vec<(u64, Vec<u8>, Vec<u8>)>,
}

enum Columns {
    Fixed(usize),
    FirstRecord(Arc<AtomicUsize>),
}

impl Reporter<'_, '_> {
    /// Is `h` a hit that counts for `a`?
    fn counts(&mut self, a: &ARec, h: &Hit) -> Result<bool, IoError> {
        if !self.acceptor.accepted(a, h)? {
            return Ok(false);
        }
        if self.matching == Match::Overlap {
            return Ok(true);
        }
        let line = self.b[usize::from(h.file)].of(h.raw, h.loc)?;
        let b_alleles = alleles(line);
        Ok(self
            .a_alleles
            .iter()
            .any(|x| b_alleles.iter().any(|y| x == y)))
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
            self.a_alleles = alleles(line);
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
    // A shared allele is the test, not the overlap: an insertion sits as a
    // point on the edge of a multi-allelic neighbour's span, where the
    // exact rule says no overlap. Equal trimmed alleles at the same trimmed
    // position are the same variant, so candidates may touch at an edge.
    let mut relate = o.relate.clone();
    if o.matching == Match::Alleles {
        relate.criteria.inclusive_points = true;
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
        a_alleles: Vec::new(),
    };
    run_relation(inputs, relate, None, &mut reporter)
}
