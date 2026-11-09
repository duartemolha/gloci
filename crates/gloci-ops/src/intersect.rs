//! `intersect`: report A's records by how they relate to the B inputs.
//!
//! A reporting layer over the relation scaffold (`relate`): the engine
//! yields candidate hits (D36), the `Acceptor` decides which count, and
//! this module writes them per `--mode`. Records are written in A's input
//! order, hits by B file and then B input order (D38).

use crate::relate::{
    Acceptor, RelateInputs, RelateOptions, RelationFiles, ReportOutput, RunOptions, b_placeholder,
    b_width, relation_chroms, report_skipped, rewrite_coords, run_files, run_relation_in,
};
use gloci_core::{ChromTable, Interval, StrandRule, Window, trimmed_alleles};
use gloci_engine::{ARec, Hit, Neighbours, SweepSummary, Visitor};
use gloci_io::{BcfText, Format, IoError, RecordReader, RecordWriter, UsageError, field};

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

/// How to print one input's records as text.
pub(crate) struct Text {
    pub(crate) format: Format,
    bcf: Option<BcfText>,
}

impl Text {
    pub(crate) fn new(r: &dyn RecordReader) -> Result<Text, IoError> {
        Ok(Text {
            format: r.format(),
            bcf: BcfText::for_reader(r)?,
        })
    }

    /// The record as a text line (a BCF record is rendered as VCF).
    pub(crate) fn of<'a>(&'a mut self, raw: &'a [u8], loc: u64) -> Result<&'a [u8], IoError> {
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
    /// Columns of a missing B (D45): those of B file 1.
    b_width: usize,
    row: Vec<u8>,
    /// A's VCF line, kept for the allele comparison.
    a_line: Vec<u8>,
    /// Clip pieces of a GFF/GTF A that were points, not written.
    point_pieces: u64,
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
                    // GFF/GTF cannot hold an empty feature (start would
                    // exceed end); BED can.
                    if shared.is_point() && matches!(self.a.format, Format::Gff3 | Format::Gtf) {
                        self.point_pieces += 1;
                        continue;
                    }
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
                    b_placeholder(self.b[0].format, self.b_width, &mut self.row);
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

/// Whether the output can hold a missing-B placeholder (D45), whose width
/// comes from B file 1.
fn needs_b_width(o: &IntersectOptions) -> bool {
    matches!(
        o.mode,
        Mode::Pairs {
            keep_unmatched: true,
            ..
        }
    )
}

/// The checks that need only the inputs' formats and count.
fn check(a: Format, b: &[Format], o: &IntersectOptions) -> Result<(), IoError> {
    validate(a, b, o)?;
    if b.len() > 1 && matches!(o.mode, Mode::Pairs { .. }) && o.relate.labels.len() != b.len() {
        return Err(usage(format!(
            "{} -b files but {} labels; give one label per -b",
            b.len(),
            o.relate.labels.len()
        )));
    }
    Ok(())
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
    check(inputs.a.format(), &b_formats, o)?;
    let (chroms, width) = if needs_b_width(o) {
        let mut chroms = relation_chroms(&inputs, &o.relate)?;
        let width = b_width(&mut inputs, &mut chroms)?;
        (Some(chroms), width)
    } else {
        (None, 3)
    };
    let (summary, points) = intersect_in(inputs, chroms, width, o, w)?;
    report_skipped(&summary);
    warn_point_pieces(points);
    Ok(summary)
}

/// One warning for the clip pieces of a GFF/GTF A that were points.
fn warn_point_pieces(n: u64) {
    if n > 0 {
        log::warn!(
            "{} clip pieces were zero-length (a GFF/GTF feature overlapping a point) and were \
not written: GFF/GTF cannot hold an empty feature; use --mode any or pairs to keep these \
features",
            gloci_core::thousands(n)
        );
    }
}

/// Run `intersect` over the files a command opened: per chromosome when
/// every input is indexed (D41), else streamed. Writes as [`intersect`].
pub fn intersect_files(
    mut files: RelationFiles,
    o: &IntersectOptions,
    w: &mut dyn RecordWriter,
    run: &RunOptions,
) -> Result<SweepSummary, IoError> {
    check(files.a().format(), &files.b_formats(), o)?;
    let width = if needs_b_width(o) {
        files.b_width(&o.relate)?
    } else {
        3
    };
    let (summary, points) = run_files(files, &o.relate, run, w, &|inputs, chroms, w| {
        intersect_in(inputs, chroms, width, o, w)
    })?;
    warn_point_pieces(points.iter().sum());
    Ok(summary)
}

/// What `intersect` writes for `o.mode` over `files`: A's records, or TSV
/// with the appended column names (D42). B's column names are needed only
/// for a VCF/BCF A, which gets a column-name line.
pub fn intersect_output(
    files: &mut RelationFiles,
    o: &IntersectOptions,
) -> Result<ReportOutput, IoError> {
    Ok(match o.mode {
        Mode::Pairs { overlap_bases, .. } => {
            let mut names = Vec::new();
            if files.b_count() > 1 {
                names.push("b_name".to_string());
            }
            if is_variant(files.a().format()) {
                names.extend(files.b_column_names(&o.relate)?);
            }
            if overlap_bases {
                names.push("overlap_bases".to_string());
            }
            ReportOutput::Tsv { appended: names }
        }
        Mode::Count => ReportOutput::Tsv {
            appended: vec!["count".to_string()],
        },
        _ => ReportOutput::SameAsA,
    })
}

/// [`intersect`] once the options are checked, with B file 1's width and,
/// when something was read with it, the chromosome table. Also returns the
/// number of clip pieces not written because they were points of a GFF/GTF A.
fn intersect_in(
    inputs: RelateInputs,
    chroms: Option<ChromTable>,
    b_width: usize,
    o: &IntersectOptions,
    w: &mut dyn RecordWriter,
) -> Result<(SweepSummary, u64), IoError> {
    let multi = inputs.b.len() > 1;
    let a = Text::new(inputs.a.as_ref())?;
    let b = inputs
        .b
        .iter()
        .map(|r| Text::new(r.as_ref()))
        .collect::<Result<Vec<_>, _>>()?;
    // Candidates for allele matching: a shared trimmed allele lies at most
    // one base from each record's footprint (an insertion is a point just
    // past an anchor base, and a footprint can end where the other starts),
    // so pad B by one base on each side. The candidates are then accepted
    // by allele equality alone. Points count at the edge of the padded
    // interval (`inclusive_points`): the padding saturates at position 0, so
    // a point at 0 (an insertion at POS 1) would otherwise sit on the edge
    // of its own padded copy and never be a candidate.
    let mut relate = o.relate.clone();
    if o.matching == Match::Alleles {
        relate.criteria.inclusive_points = true;
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
        b_width,
        row: Vec::new(),
        a_line: Vec::new(),
        point_pieces: 0,
    };
    let summary = run_relation_in(inputs, chroms, relate, None, &mut reporter)?;
    Ok((summary, reporter.point_pieces))
}
