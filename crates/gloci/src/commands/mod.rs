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
pub mod vcf;
pub mod view;
pub mod windows;

use crate::Ctx;
use crate::cli::{ChromOrderArg, RelationSortArg, SortArg, VcfSpan};
use gloci_io::{
    Format, InvalidPolicy, ReaderOptions, RecordReader, UsageError, chrom_table_for, open_reader,
};
use gloci_ops::{
    InputPaths, RelationFiles, RunOptions, SortSettings, sort_input_options, sort_together,
};

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

pub trait SortArgLike {
    fn is_sort(&self) -> bool;
    fn max_mem(&self) -> u64;
}

impl SortArgLike for SortArg {
    fn is_sort(&self) -> bool {
        self.sort
    }
    fn max_mem(&self) -> u64 {
        self.max_mem
    }
}

impl SortArgLike for RelationSortArg {
    fn is_sort(&self) -> bool {
        self.sort
    }
    fn max_mem(&self) -> u64 {
        self.max_mem
    }
}

impl<T: SortArgLike> SortArgLike for &T {
    fn is_sort(&self) -> bool {
        (*self).is_sort()
    }
    fn max_mem(&self) -> u64 {
        (*self).max_mem()
    }
}

/// The hint extra for a relation command's order errors: none once `--sort`
/// has sorted the inputs (the sorted stream cannot be out of order).
pub fn relation_hint(sort: &RelationSortArg) -> &'static str {
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
pub fn sort_aware(sort: impl SortArgLike, ro: ReaderOptions) -> ReaderOptions {
    if sort.is_sort() {
        sort_input_options(&ro)
    } else {
        ro
    }
}

/// The settings of `--sort`, if given; `span` is the user's `--vcf-span`.
pub fn sort_settings(ctx: &Ctx, sort: impl SortArgLike, span: VcfSpan) -> Option<SortSettings> {
    sort.is_sort().then(|| SortSettings {
        max_mem: sort.max_mem(),
        tmp_dir: ctx.tmp_dir.clone(),
        threads: ctx.threads,
        vcf_span: user_span(span),
    })
}

/// Open a relation's A and B inputs, each once.
pub fn relation_files(
    a: &str,
    b: &[String],
    ro_a: ReaderOptions,
    ro_b: ReaderOptions,
) -> anyhow::Result<RelationFiles> {
    Ok(RelationFiles::open(InputPaths {
        a: a.to_string(),
        b: b.to_vec(),
        ro_a,
        ro_b,
    })?)
}

/// How a relation over files runs: threads and the spill directory.
pub fn run_options(ctx: &Ctx) -> RunOptions {
    RunOptions {
        threads: ctx.threads,
        tmp_dir: ctx.tmp_dir.clone(),
    }
}

/// The output's name for messages: the path, or `stdout`. Write errors from
/// the ops name the placeholder `output`; `IoError::with_output_name` puts
/// this in its place.
pub fn dest_name(spec: &gloci_io::OutputSpec) -> String {
    spec.path
        .as_ref()
        .map_or_else(|| "stdout".to_string(), |p| p.display().to_string())
}
