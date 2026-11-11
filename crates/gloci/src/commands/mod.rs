pub mod cat;
pub mod check;
pub mod closest;
pub mod cluster;
pub mod complement;
pub mod convert;
pub mod coverage;
pub mod dict;
pub mod faidx;
pub mod flags;
pub mod flank;
pub mod genomecov;
pub mod genomes;
pub mod getfasta;
pub mod groupby;
pub mod index;
pub mod intersect;
pub mod map;
pub mod maskfasta;
pub mod merge;
pub mod nuc;
pub mod partition;
pub mod shift;
pub mod slop;
pub mod sort;
pub mod spacing;
pub mod subtract;
pub mod vcf;
pub mod view;
pub mod windows;

use crate::Ctx;
use crate::cli::{ChromMapArg, ChromOrderArg, RegionArgs, RelationSortArg, SortArg, VcfSpan};
use gloci_io::{
    Format, InvalidPolicy, ReaderOptions, RecordReader, UsageError, chrom_table_for, open_reader,
};
use gloci_ops::{
    InputPaths, RelationFiles, RunOptions, SortSettings, ViewOptions, sort_input_options,
    sort_readers, sort_together,
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
    pub chrom_map: &'a ChromMapArg,
    /// `Some` when the command reads alignments: how (the read filters).
    pub alignments: Option<gloci_io::AlignOptions>,
}

/// Open the inputs of a streaming command and build their chromosome table.
/// With `--sort` the inputs are sorted first (ops: [`sort_together`]) and a single sorted reader comes back. `check` sees the
/// inputs' formats before anything is sorted.
pub fn open_inputs(
    ctx: &Ctx,
    a: &InputArgs,
    check: &dyn Fn(&[Format]) -> Result<(), UsageError>,
) -> anyhow::Result<(Vec<Box<dyn RecordReader>>, gloci_core::ChromTable)> {
    let (readers, mut chroms) = open_readers(ctx, a, check)?;
    match sort_settings(ctx, a.sort, a.vcf_span) {
        None => Ok((readers, chroms)),
        Some(settings) => {
            let sorted = sort_together(readers, &mut chroms, &settings)?;
            Ok((sorted, chroms))
        }
    }
}

/// [`open_inputs`] for commands that must tell their inputs apart: with
/// `--sort` each input is sorted by itself and one reader per input comes
/// back, in order.
pub fn open_inputs_each(
    ctx: &Ctx,
    a: &InputArgs,
    check: &dyn Fn(&[Format]) -> Result<(), UsageError>,
) -> anyhow::Result<(Vec<Box<dyn RecordReader>>, gloci_core::ChromTable)> {
    let (readers, mut chroms) = open_readers(ctx, a, check)?;
    match sort_settings(ctx, a.sort, a.vcf_span) {
        None => Ok((readers, chroms)),
        Some(settings) => {
            let sorted = readers
                .into_iter()
                .map(|r| sort_readers(vec![r], &mut chroms, &settings))
                .collect::<Result<Vec<_>, _>>()?;
            Ok((sorted, chroms))
        }
    }
}

/// The options the inputs of [`open_inputs`] are read with.
pub fn reader_options(ctx: &Ctx, a: &InputArgs) -> anyhow::Result<ReaderOptions> {
    Ok(sort_aware(
        a.sort,
        ReaderOptions {
            accept: gloci_io::Accept {
                alignments: a.alignments.clone(),
                ..gloci_io::Accept::default()
            },
            chrom_map: a.chrom_map.load()?,
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
    ))
}

fn open_readers(
    ctx: &Ctx,
    a: &InputArgs,
    check: &dyn Fn(&[Format]) -> Result<(), UsageError>,
) -> anyhow::Result<(Vec<Box<dyn RecordReader>>, gloci_core::ChromTable)> {
    let ro = reader_options(ctx, a)?;
    let readers = a
        .files
        .iter()
        .map(|f| open_reader(f, &ro))
        .collect::<Result<Vec<_>, _>>()?;
    let formats: Vec<Format> = readers.iter().map(|r| r.format()).collect();
    check(&formats)?;
    let refs: Vec<&dyn RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let chroms = chrom_table_for(a.genome, &refs, a.chrom_order.chrom_order)?;
    Ok((readers, chroms))
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

/// Restrict a relation's A to `--region`/`--regions-file` (D75); `unordered`
/// when `--unsorted` or `--sort` lets A be out of order.
pub fn restrict_a(
    files: RelationFiles,
    r: &RegionArgs,
    inclusive_points: bool,
    unordered: bool,
) -> anyhow::Result<RelationFiles> {
    let regions = ViewOptions {
        regions: r.region.clone(),
        regions_file: r
            .regions_file
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        inclusive_points,
        ..ViewOptions::default()
    };
    Ok(files.restrict_a(&regions, unordered)?)
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
