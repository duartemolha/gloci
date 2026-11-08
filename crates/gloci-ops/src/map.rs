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

use crate::agg::{Accumulators, AggSpec, Op, resolve};
use crate::intersect::Text;
use crate::relate::{Acceptor, RelateInputs, RelateOptions, run_relation};
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
        Format::Vcf => "VCF",
        Format::Bcf => "BCF",
        Format::Gff3 => "GFF3",
        Format::Gtf => "GTF",
    }
}

fn family(f: Format) -> u8 {
    match f {
        Format::Vcf | Format::Bcf => 0,
        Format::Bed | Format::BedGraph => 1,
        Format::Gff3 | Format::Gtf => 2,
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
    wants_bases: bool,
    acceptor: Acceptor<'o>,
    accs: Accumulators,
    w: &'w mut dyn RecordWriter,
    a: Text,
    b: Vec<Text>,
    b_names: Vec<String>,
    row: Vec<u8>,
}

impl Visitor for Reporter<'_, '_> {
    fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], _: &Neighbours<'_>) -> Result<(), IoError> {
        let mut any = false;
        for h in hits {
            if !self.acceptor.accepted(a, h)? {
                continue;
            }
            any = true;
            let bases = if self.wants_bases {
                self.acceptor.overlap_len(a, h)?
            } else {
                0
            };
            let f = usize::from(h.file);
            let line = self.b[f].of(h.raw, h.loc)?;
            self.accs.add(line, bases, &self.b_names[f], h.loc)?;
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
        wants_bases: o.aggs.iter().any(|a| a.op == Op::Bases),
        acceptor,
        accs: Accumulators::new(resolved, o.empty.clone(), o.delim.clone()),
        w,
        a,
        b,
        b_names,
        row: Vec::new(),
    };
    run_relation(inputs, &o.relate, None, &mut reporter)
}
