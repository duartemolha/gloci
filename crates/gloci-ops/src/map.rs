//! `map`: aggregate columns of the B records over each A record
//! (bedtools map, bedops bedmap).
//!
//! A reporting layer over the relation scaffold (`relate`). For each A, in
//! A's input order, the accepted hits (D36) are folded in hit order (D38)
//! by the aggregation engine (D46), and one TSV row is written: A as text,
//! then one column per aggregation. Columns are resolved against the first
//! `-b` file's format and header; all `-b` files must be one format family
//! (D57).
//!
//! `bases` is the overlap of A's span with the hit's own (unpadded) span,
//! so with a window it is 0 for a hit that is within the window only.

use crate::agg::{Accumulators, AggSpec, HitInfo, resolve};
use crate::intersect::Text;
use crate::relate::{
    Acceptor, RelateInputs, RelateOptions, RelationFiles, RunOptions, report_skipped, run_files,
    run_relation_in,
};
use gloci_core::{ChromTable, Interval, overlap_len};
use gloci_engine::{ARec, Hit, Neighbours, SweepSummary, Visitor};
use gloci_io::{Format, IoError, RecordWriter, UsageError};

#[derive(Debug, Clone)]
pub struct MapOptions {
    pub aggs: Vec<AggSpec>,
    /// Drop an A with no accepted hit.
    pub skip_unmatched: bool,
    /// Value of an aggregation with no values.
    pub empty: String,
    /// Separator of `collapse` and `distinct` lists.
    pub delim: String,
    pub relate: RelateOptions,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn format_name(f: Format) -> &'static str {
    match f {
        Format::Bed => "BED",
        Format::BedGraph => "bedGraph",
        Format::Wig => "WIG",
        Format::Vcf => "VCF",
        Format::Bcf => "BCF",
        Format::Gff3 => "GFF3",
        Format::Gtf => "GTF",
        Format::Bam
        | Format::Cram
        | Format::Sam
        | Format::Fasta
        | Format::Psl
        | Format::Rmsk
        | Format::Bedpe => f.name(),
    }
}

fn family(f: Format) -> u8 {
    match f {
        Format::Vcf | Format::Bcf => 0,
        Format::Bed | Format::BedGraph | Format::Wig => 1,
        Format::Gff3 | Format::Gtf => 2,
        // Refused at input (the gate) before a map can compare families.
        Format::Bam
        | Format::Cram
        | Format::Sam
        | Format::Fasta
        | Format::Psl
        | Format::Rmsk
        | Format::Bedpe => 3,
    }
}

/// The option checks that need only the inputs' formats; `names` are the
/// `-b` files, for the message.
pub fn validate(o: &MapOptions, b: &[Format], names: &[String]) -> Result<(), IoError> {
    if o.aggs.is_empty() {
        return Err(usage(
            "map needs at least one -c COL:OP, for example -c 5:sum",
        ));
    }
    if let Some(first) = b.first()
        && let Some(i) = b.iter().position(|f| family(*f) != family(*first))
    {
        let name = |i: usize| names.get(i).map_or("?", String::as_str);
        return Err(usage(format!(
            "all -b files must be one format family so that -c columns mean the same thing, \
but {} is {} and {} is {}; convert them or run map once per format",
            name(0),
            format_name(*first),
            name(i),
            format_name(b[i])
        )));
    }
    Ok(())
}

struct Reporter<'o, 'w> {
    skip_unmatched: bool,
    wants_overlap: bool,
    acceptor: Acceptor<'o>,
    accs: Accumulators,
    w: &'w mut dyn RecordWriter,
    a: Text,
    b: Vec<Text>,
    b_names: Vec<String>,
    row: Vec<u8>,
    /// Scratch: the shared intervals of the current hit.
    shared: Vec<Interval>,
    split: bool,
    /// A's blocks and their total length under `--split`, kept per A record.
    a_blocks: Option<(u64, Vec<Interval>, u64)>,
}

impl Reporter<'_, '_> {
    /// The overlap facts of `h` with `a` (`shared` is filled; the returned
    /// slice is left empty for the caller to point at it). With `--split`
    /// lengths and shared bases are those of the blocks.
    fn overlap(&mut self, a: &ARec<'_>, h: &Hit<'_>) -> Result<HitInfo<'static>, IoError> {
        self.shared.clear();
        let bases = self.acceptor.overlap_len(a, h)?;
        if !self.split {
            let (x, y) = (a.meta.span, h.span);
            if overlap_len(x, y) > 0 {
                self.shared.push(Interval {
                    start: x.start.max(y.start),
                    end: x.end.min(y.end),
                });
            }
            return Ok(HitInfo {
                bases,
                a_len: x.len(),
                b_len: y.len(),
                shared: &[],
            });
        }
        let total = |v: &[Interval]| v.iter().map(Interval::len).sum::<u64>();
        if self
            .a_blocks
            .as_ref()
            .is_none_or(|(l, ..)| *l != a.meta.loc)
        {
            let ab = self.acceptor.a_blocks(a)?.to_vec();
            let len = total(&ab);
            self.a_blocks = Some((a.meta.loc, ab, len));
        }
        let (_, ab, a_len) = self.a_blocks.as_ref().expect("set above");
        let bb = self.acceptor.hit_blocks(h)?;
        for x in ab {
            for y in &bb {
                if overlap_len(*x, *y) > 0 {
                    self.shared.push(Interval {
                        start: x.start.max(y.start),
                        end: x.end.min(y.end),
                    });
                }
            }
        }
        Ok(HitInfo {
            bases,
            a_len: *a_len,
            b_len: total(&bb),
            shared: &[],
        })
    }
}

impl Visitor for Reporter<'_, '_> {
    fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], _: &Neighbours<'_>) -> Result<(), IoError> {
        let mut any = false;
        for h in hits {
            if !self.acceptor.accepted(a, h)? {
                continue;
            }
            any = true;
            let f = usize::from(h.file);
            let info = if self.wants_overlap {
                Some(self.overlap(a, h)?)
            } else {
                None
            };
            let line = self.b[f].of(h.raw, h.loc)?;
            if let Some(info) = info {
                let hit = HitInfo {
                    shared: &self.shared,
                    ..info
                };
                self.accs.add_hit(line, &hit, &self.b_names[f], h.loc)?;
            } else {
                self.accs.add(line, 0, &self.b_names[f], h.loc)?;
            }
        }
        if !any && self.skip_unmatched {
            // nothing was added; leave the accumulators clean
            self.row.clear();
            self.accs.finish_into(&mut self.row);
            return Ok(());
        }
        self.row.clear();
        let line = self.a.of(a.raw, a.meta.loc)?;
        self.row.extend_from_slice(line);
        self.accs.finish_into(&mut self.row);
        self.w
            .write_record(&self.row)
            .map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })
    }
}

/// Run `map` over `inputs`, writing TSV rows to `w` (`open_report` with
/// `ReportOutput::Tsv`). Nothing is finished here; the caller calls
/// `finish`.
pub fn map(
    inputs: RelateInputs,
    o: &MapOptions,
    w: &mut dyn RecordWriter,
) -> Result<SweepSummary, IoError> {
    let summary = map_in(inputs, None, o, w)?;
    report_skipped(&summary);
    Ok(summary)
}

/// Run `map` over the files a command opened: per chromosome when every
/// input is indexed (D41), else streamed. Writes as [`map`].
pub fn map_files(
    files: RelationFiles,
    o: &MapOptions,
    w: &mut dyn RecordWriter,
    run: &RunOptions,
) -> Result<SweepSummary, IoError> {
    let (summary, _) = run_files(files, &o.relate, run, w, &|inputs, chroms, w| {
        Ok((map_in(inputs, chroms, o, w)?, ()))
    })?;
    Ok(summary)
}

fn map_in(
    inputs: RelateInputs,
    chroms: Option<ChromTable>,
    o: &MapOptions,
    w: &mut dyn RecordWriter,
) -> Result<SweepSummary, IoError> {
    let b_formats: Vec<Format> = inputs.b.iter().map(|r| r.format()).collect();
    let b_names: Vec<String> = inputs.b.iter().map(|r| r.name().to_string()).collect();
    validate(o, &b_formats, &b_names)?;
    let resolved = match inputs.b.first() {
        Some(first) => resolve(&o.aggs, first.format(), first.header())?,
        None => return Err(usage("map needs at least one -b file")),
    };
    let a = Text::new(inputs.a.as_ref())?;
    let b = inputs
        .b
        .iter()
        .map(|r| Text::new(r.as_ref()))
        .collect::<Result<Vec<_>, _>>()?;
    let acceptor = Acceptor::new(&o.relate, &inputs);
    let mut reporter = Reporter {
        skip_unmatched: o.skip_unmatched,
        wants_overlap: o.aggs.iter().any(|a| a.op.needs_overlap()),
        acceptor,
        accs: Accumulators::new(resolved, o.empty.clone(), o.delim.clone()),
        w,
        a,
        b,
        b_names,
        row: Vec::new(),
        shared: Vec::new(),
        split: o.relate.split,
        a_blocks: None,
    };
    run_relation_in(inputs, chroms, &o.relate, None, &mut reporter)
}
