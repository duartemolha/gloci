//! `closest`: for each A, the nearest B (spec §4.4).
//!
//! A reporting layer over the relation scaffold (`relate`). The engine
//! yields, per A and B file, the accepted overlapping hits (distance 0) and
//! the `k` nearest non-overlapping B on each side plus ties
//! ([`Neighbours`]). This module pools them, drops what the ignore flags
//! exclude, selects the `k` nearest and writes TSV rows.
//!
//! - **Up and downstream** is the sign of `signed_distance` in the
//!   `--signed` frame (`ref` when absent). An overlap has sign 0: it is
//!   neither, so `--ignore-upstream` and `--ignore-downstream` keep it.
//! - **Selection** sorts by (distance, B file, B input order). The `k`
//!   nearest are kept, and the ties at the k-th distance: all of them with
//!   [`Ties::All`], else only as many of them as `k` still needs, the first
//!   (or the last) by (file, input order). So with `k = 1` there is exactly
//!   one row, and with `k > 1` a tie below the k-th distance is never cut.
//! - **`--per-file`** selects within each B file and writes each file's
//!   rows in file order; a file with no candidate writes nothing.
//! - **`--each-side`** writes one row per rank: A, the upstream B with its
//!   distance, then the downstream B with its distance. An overlapping B
//!   counts on both sides. The distance there is always signed.
//! - **No candidate:** one row with the D45 placeholder for B and `-1` as
//!   the distance.

use crate::intersect::Text;
use crate::relate::{
    Acceptor, RelateInputs, RelateOptions, RelationFiles, ReportOutput, RunOptions, b_placeholder,
    b_width, relation_chroms, run_files, run_relation_in,
};
use gloci_core::{ChromTable, SignFrame, signed_distance};
use gloci_engine::{ARec, Hit, NearestConfig, Neighbours, SweepSummary, Visitor};
use gloci_io::{Format, IoError, RecordWriter, UsageError};

/// Which of several equally near B to report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ties {
    All,
    First,
    Last,
}

#[derive(Debug, Clone)]
pub struct ClosestOptions {
    /// How many nearest B to report per A (at least 1); per side with
    /// `each_side`, per file with `per_file`.
    pub k: usize,
    pub ties: Ties,
    /// Write the distance column. `signed` implies it.
    pub distance: bool,
    /// The frame of signed distances, and of up/downstream.
    pub signed: Option<SignFrame>,
    pub ignore_overlaps: bool,
    pub ignore_upstream: bool,
    pub ignore_downstream: bool,
    pub max_distance: Option<u64>,
    pub each_side: bool,
    pub per_file: bool,
    /// Only the strand rule, `inclusive_points`, `unsorted`, labels and the
    /// genome options apply; windows and thresholds are not offered.
    pub relate: RelateOptions,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// The option checks that need no input.
pub fn validate(o: &ClosestOptions) -> Result<(), IoError> {
    if o.k == 0 {
        return Err(usage("-k must be at least 1; give -k 1 for the nearest B"));
    }
    if o.relate.window.is_some() {
        return Err(usage("closest takes no window (-w, --left, --right)"));
    }
    Ok(())
}

/// A candidate: a B with its distance from A.
#[derive(Clone, Copy)]
struct Cand<'a> {
    /// In the `--signed` frame; 0 for an overlap.
    signed: i64,
    abs: u64,
    hit: Hit<'a>,
}

impl Cand<'_> {
    /// Order by distance, then B file, then B input order.
    fn key(&self) -> (u64, u16, u64) {
        (self.abs, self.hit.file, self.hit.loc)
    }
}

/// The candidates of `c` (sorted by distance) to report: the `k` nearest and
/// their ties, per `ties`. Returned as two slices to chain.
fn select<'c, 'a>(c: &'c [Cand<'a>], k: usize, ties: Ties) -> (&'c [Cand<'a>], &'c [Cand<'a>]) {
    if c.len() <= k {
        return (c, &[]);
    }
    let d = c[k - 1].abs;
    let below = c.partition_point(|x| x.abs < d);
    let through = c.partition_point(|x| x.abs <= d);
    let group = &c[below..through];
    let need = k - below;
    let chosen = match ties {
        Ties::All => group,
        Ties::First => &group[..need],
        Ties::Last => &group[group.len() - need..],
    };
    (&c[..below], chosen)
}

struct Reporter<'o, 'w> {
    o: &'o ClosestOptions,
    frame: SignFrame,
    acceptor: Acceptor<'o>,
    w: &'w mut dyn RecordWriter,
    a: Text,
    b: Vec<Text>,
    /// One label per B file when there are several, else empty.
    labels: &'o [String],
    /// Columns of a missing B (D45): those of B file 1.
    b_width: usize,
    row: Vec<u8>,
    files: u16,
}

impl Reporter<'_, '_> {
    fn put(&mut self) -> Result<(), IoError> {
        self.w
            .write_record(&self.row)
            .map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })
    }

    /// Append `\t[label\t]B[\tdistance]`; `None` is the missing-B
    /// placeholder.
    fn push_b(&mut self, c: Option<&Cand>, with_distance: bool) -> Result<(), IoError> {
        match c {
            Some(c) => {
                if !self.labels.is_empty() {
                    self.row.push(b'\t');
                    self.row
                        .extend_from_slice(self.labels[usize::from(c.hit.file)].as_bytes());
                }
                self.row.push(b'\t');
                let line = self.b[usize::from(c.hit.file)].of(c.hit.raw, c.hit.loc)?;
                self.row.extend_from_slice(line);
                if with_distance {
                    // each-side and --signed print the signed value
                    let d = if self.shows_sign() {
                        c.signed.to_string()
                    } else {
                        c.abs.to_string()
                    };
                    self.row.push(b'\t');
                    self.row.extend_from_slice(d.as_bytes());
                }
            }
            None => {
                if !self.labels.is_empty() {
                    self.row.extend_from_slice(b"\t.");
                }
                self.row.push(b'\t');
                b_placeholder(self.b[0].format, self.b_width, &mut self.row);
                if with_distance {
                    self.row.extend_from_slice(b"\t-1");
                }
            }
        }
        Ok(())
    }

    fn shows_sign(&self) -> bool {
        self.o.signed.is_some() || self.o.each_side
    }

    fn shows_distance(&self) -> bool {
        self.o.distance || self.o.signed.is_some()
    }

    fn start_row(&mut self, a: &ARec) -> Result<(), IoError> {
        self.row.clear();
        let line = self.a.of(a.raw, a.meta.loc)?;
        self.row.extend_from_slice(line);
        Ok(())
    }

    /// One row per selected candidate of `c`, or nothing when it is empty.
    fn rows(&mut self, a: &ARec, c: &[Cand]) -> Result<usize, IoError> {
        let (head, tail) = select(c, self.o.k, self.o.ties);
        let with_distance = self.shows_distance();
        for x in head.iter().chain(tail) {
            self.start_row(a)?;
            self.push_b(Some(x), with_distance)?;
            self.put()?;
        }
        Ok(head.len() + tail.len())
    }

    /// The rows of `--each-side` for the candidates `c`: ranks of the
    /// upstream selection beside those of the downstream one.
    fn side_rows(&mut self, a: &ARec, c: &[Cand]) -> Result<usize, IoError> {
        let pick = |keep: fn(i64) -> bool| -> Vec<Cand> {
            let side: Vec<Cand> = c.iter().filter(|x| keep(x.signed)).copied().collect();
            let (head, tail) = select(&side, self.o.k, self.o.ties);
            head.iter().chain(tail).copied().collect()
        };
        let with_distance = self.shows_distance();
        let up = pick(|s| s <= 0);
        let down = pick(|s| s >= 0);
        for i in 0..up.len().max(down.len()) {
            self.start_row(a)?;
            self.push_b(up.get(i), with_distance)?;
            self.push_b(down.get(i), with_distance)?;
            self.put()?;
        }
        Ok(up.len().max(down.len()))
    }
}

impl Visitor for Reporter<'_, '_> {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        let o = self.o;
        let inclusive = o.relate.criteria.inclusive_points;
        let mut cands: Vec<Cand> = Vec::new();
        if !o.ignore_overlaps {
            for h in hits {
                if self.acceptor.accepted(a, h)? {
                    cands.push(Cand {
                        signed: 0,
                        abs: 0,
                        hit: *h,
                    });
                }
            }
        }
        for f in 0..self.files {
            for h in near.left(f).iter().chain(near.right(f)) {
                let signed = signed_distance(
                    a.meta.span,
                    a.meta.strand,
                    h.span,
                    h.strand,
                    self.frame,
                    inclusive,
                );
                if (signed < 0 && o.ignore_upstream) || (signed > 0 && o.ignore_downstream) {
                    continue;
                }
                cands.push(Cand {
                    signed,
                    abs: signed.unsigned_abs(),
                    hit: *h,
                });
            }
        }
        if o.per_file {
            cands.sort_by_key(|c| (c.hit.file, c.key()));
        } else {
            cands.sort_by_key(Cand::key);
        }
        let mut written = 0;
        if o.per_file {
            let mut rest = &cands[..];
            while let Some(first) = rest.first() {
                let n = rest.partition_point(|c| c.hit.file == first.hit.file);
                let (group, tail) = rest.split_at(n);
                written += if o.each_side {
                    self.side_rows(a, group)?
                } else {
                    self.rows(a, group)?
                };
                rest = tail;
            }
        } else {
            written += if o.each_side {
                self.side_rows(a, &cands)?
            } else {
                self.rows(a, &cands)?
            };
        }
        if written == 0 {
            self.start_row(a)?;
            let with_distance = self.shows_distance();
            self.push_b(None, with_distance)?;
            if o.each_side {
                self.push_b(None, with_distance)?;
            }
            self.put()?;
        }
        Ok(())
    }
}

/// The checks that need only the number of B files.
fn check(n_b: usize, o: &ClosestOptions) -> Result<(), IoError> {
    validate(o)?;
    if n_b > 1 && o.relate.labels.len() != n_b {
        return Err(usage(format!(
            "{} -b files but {} labels; give one label per -b",
            n_b,
            o.relate.labels.len()
        )));
    }
    u16::try_from(n_b).map_err(|_| usage("too many -b files; closest takes at most 65535"))?;
    Ok(())
}

/// Run `closest` over `inputs`, writing TSV rows to `w` (`open_report` with
/// `ReportOutput::Tsv`). Nothing is finished here; the caller calls
/// `finish`.
pub fn closest(
    mut inputs: RelateInputs,
    o: &ClosestOptions,
    w: &mut dyn RecordWriter,
) -> Result<SweepSummary, IoError> {
    check(inputs.b.len(), o)?;
    let mut chroms = relation_chroms(&inputs, &o.relate)?;
    let width = b_width(&mut inputs, &mut chroms)?;
    closest_in(inputs, Some(chroms), width, o, w)
}

/// Run `closest` over the files a command opened: per chromosome when
/// every input is indexed (D41), else streamed. Writes as [`closest`].
pub fn closest_files(
    mut files: RelationFiles,
    o: &ClosestOptions,
    w: &mut dyn RecordWriter,
    run: &RunOptions,
) -> Result<SweepSummary, IoError> {
    check(files.b_count(), o)?;
    let width = files.b_width(&o.relate)?;
    let (summary, _) = run_files(files, &o.relate, run, w, &|inputs, chroms, w| {
        Ok((closest_in(inputs, chroms, width, o, w)?, ()))
    })?;
    Ok(summary)
}

/// The TSV columns `closest` appends (D42). They are named only for a
/// VCF/BCF A, which gets a column-name line.
pub fn closest_output(
    files: &mut RelationFiles,
    o: &ClosestOptions,
) -> Result<ReportOutput, IoError> {
    let mut names = Vec::new();
    if matches!(files.a().format(), Format::Vcf | Format::Bcf) {
        let multi = files.b_count() > 1;
        let shows_distance = o.distance || o.signed.is_some();
        let sides: &[&str] = if o.each_side {
            &["up_", "down_"]
        } else {
            &[""]
        };
        let b = files.b_column_names(&o.relate)?;
        for side in sides {
            if multi {
                names.push(format!("{side}b_name"));
            }
            names.extend(b.iter().map(|n| format!("{side}{n}")));
            if shows_distance {
                names.push(format!("{side}distance"));
            }
        }
    }
    Ok(ReportOutput::Tsv { appended: names })
}

/// [`closest`] once the options are checked, with B file 1's width and,
/// when something was read with it, the chromosome table.
fn closest_in(
    inputs: RelateInputs,
    chroms: Option<ChromTable>,
    b_width: usize,
    o: &ClosestOptions,
    w: &mut dyn RecordWriter,
) -> Result<SweepSummary, IoError> {
    let multi = inputs.b.len() > 1;
    let files = u16::try_from(inputs.b.len())
        .map_err(|_| usage("too many -b files; closest takes at most 65535"))?;
    let a = Text::new(inputs.a.as_ref())?;
    let b = inputs
        .b
        .iter()
        .map(|r| Text::new(r.as_ref()))
        .collect::<Result<Vec<_>, _>>()?;
    let frame = o.signed.unwrap_or(SignFrame::Ref);
    // `--signed b` with a filter on direction: which side of A a B is on
    // does not tell its direction, so keep k per B strand.
    let split_by_b_strand =
        frame == SignFrame::B && (o.ignore_upstream || o.ignore_downstream || o.each_side);
    let nearest = NearestConfig {
        k: o.k,
        max_distance: o.max_distance,
        strand: o.relate.criteria.strand,
        split_by_b_strand,
    };
    let relate = &o.relate;
    let acceptor = Acceptor::new(relate, &inputs);
    let mut reporter = Reporter {
        o,
        frame,
        acceptor,
        w,
        a,
        b,
        labels: if multi { &relate.labels } else { &[] },
        b_width,
        row: Vec::new(),
        files,
    };
    run_relation_in(inputs, chroms, relate, Some(nearest), &mut reporter)
}
