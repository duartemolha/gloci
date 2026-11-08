pub mod check;
pub mod closest;
pub mod complement;
pub mod convert;
pub mod coverage;
pub mod flank;
pub mod genomes;
pub mod index;
pub mod intersect;
pub mod map;
pub mod merge;
pub mod slop;
pub mod sort;
pub mod subtract;
pub mod view;
pub mod windows;

use crate::Ctx;
use crate::cli::{ChromOrderArg, SortArg, VcfSpan};
use gloci_io::{
    Format, InvalidPolicy, ReaderOptions, RecordReader, UsageError, chrom_table_for, open_reader,
};
use gloci_ops::{SortSettings, sort_input_options, sort_together};

/// Appended to the hint of an order error from a command with `--sort`.
pub const SORT_FLAG_HINT: &str = ", or pass --sort to sort the inputs first";

/// What the `merge`/`complement` commands need to open their inputs.
pub struct InputArgs<'a> {
    pub files: &'a [String],
    pub format: Option<&'a str>,
    pub vcf_span: VcfSpan,
    pub invalid: bool,
    pub genome: Option<&'a str>,
    pub chrom_order: &'a ChromOrderArg,
    pub sort: &'a SortArg,
}

/// Open the inputs of a streaming command and build their chromosome table.
/// With `--sort` the inputs are sorted first (ops: [`sort_together`]) and a single sorted reader comes back. `check` sees the
/// inputs' formats before anything is sorted.
pub fn open_inputs(
    ctx: &Ctx,
    a: &InputArgs,
    check: &dyn Fn(&[Format]) -> Result<(), UsageError>,
) -> anyhow::Result<(Vec<Box<dyn RecordReader>>, gloci_core::ChromTable)> {
    let ro = sort_aware(
        a.sort,
        ReaderOptions {
            threads: ctx.threads,
            format: match a.format {
                None => None,
                Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
            },
            vcf_span: user_span(a.vcf_span),
            invalid: if a.invalid {
                InvalidPolicy::Skip
            } else {
                InvalidPolicy::Fail
            },
        },
    );
    let readers = a
        .files
        .iter()
        .map(|f| open_reader(f, &ro))
        .collect::<Result<Vec<_>, _>>()?;
    let formats: Vec<Format> = readers.iter().map(|r| r.format()).collect();
    check(&formats)?;
    let refs: Vec<&dyn RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let mut chroms = chrom_table_for(a.genome, &refs, a.chrom_order.chrom_order)?;
    match sort_settings(ctx, a.sort, a.vcf_span) {
        None => Ok((readers, chroms)),
        Some(settings) => {
            let sorted = sort_together(readers, &mut chroms, &settings)?;
            Ok((sorted, chroms))
        }
    }
}

/// Appended to the hint of an order error from a relation command: both ways
/// out, sorting first (sorted output) or `--unsorted` (A's order kept, B in
/// memory).
pub const RELATION_FLAGS_HINT: &str = ", or pass --sort to sort the inputs first (output in sorted order), or --unsorted to keep A's order (every -b file is then held in memory)";

/// The hint extra for a relation command's order errors: none once `--sort`
/// has sorted the inputs (the sorted stream cannot be out of order).
pub fn relation_hint(sort: &SortArg) -> &'static str {
    if sort.sort { "" } else { RELATION_FLAGS_HINT }
}

pub fn user_span(span: VcfSpan) -> gloci_core::VcfSpan {
    match span {
        VcfSpan::Allele => gloci_core::VcfSpan::Allele,
        VcfSpan::Ref => gloci_core::VcfSpan::Ref,
    }
}

/// Reader options for the inputs of a command with `--sort`: with the flag
/// on they are opened the way the sorter needs (see [`sort_input_options`]).
pub fn sort_aware(sort: &SortArg, ro: ReaderOptions) -> ReaderOptions {
    if sort.sort {
        sort_input_options(&ro)
    } else {
        ro
    }
}

/// The settings of `--sort`, if given; `span` is the user's `--vcf-span`.
pub fn sort_settings(ctx: &Ctx, sort: &SortArg, span: VcfSpan) -> Option<SortSettings> {
    sort.sort.then(|| SortSettings {
        max_mem: sort.max_mem,
        tmp_dir: ctx.tmp_dir.clone(),
        threads: ctx.threads,
        vcf_span: user_span(span),
    })
}
