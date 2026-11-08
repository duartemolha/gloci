pub mod check;
pub mod closest;
pub mod complement;
pub mod coverage;
pub mod genomes;
pub mod index;
pub mod intersect;
pub mod map;
pub mod merge;
pub mod sort;
pub mod subtract;
pub mod view;

use crate::Ctx;
use crate::cli::{ChromOrderArg, SortArg, VcfSpan};
use gloci_io::IoError;
use gloci_io::{
    Format, InvalidPolicy, ReaderOptions, RecordReader, UsageError, chrom_table_for, open_reader,
};
use gloci_ops::{SortSettings, sort_readers};

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
/// With `--sort` the inputs are sorted first (one reusable step, see
/// [`sort_readers`]) and a single sorted reader comes back. `check` sees the
/// inputs' formats before anything is sorted.
pub fn open_inputs(
    ctx: &Ctx,
    a: &InputArgs,
    check: &dyn Fn(&[Format]) -> Result<(), UsageError>,
) -> anyhow::Result<(Vec<Box<dyn RecordReader>>, gloci_core::ChromTable)> {
    let ro = ReaderOptions {
        threads: ctx.threads,
        format: match a.format {
            None => None,
            Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
        },
        // `gloci sort` keys VCF records by their REF span, so sorting reads
        // them that way; the sorted stream is read with the user's span.
        vcf_span: if a.sort.sort {
            gloci_core::VcfSpan::Ref
        } else {
            match a.vcf_span {
                VcfSpan::Allele => gloci_core::VcfSpan::Allele,
                VcfSpan::Ref => gloci_core::VcfSpan::Ref,
            }
        },
        invalid: if a.invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let readers = a
        .files
        .iter()
        .map(|f| open_reader(f, &ro))
        .collect::<Result<Vec<_>, _>>()?;
    let formats: Vec<Format> = readers.iter().map(|r| r.format()).collect();
    check(&formats)?;
    let refs: Vec<&dyn RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let mut chroms = chrom_table_for(a.genome, &refs, a.chrom_order.chrom_order)?;
    if !a.sort.sort {
        return Ok((readers, chroms));
    }
    let reopen = ReaderOptions {
        vcf_span: match a.vcf_span {
            VcfSpan::Allele => gloci_core::VcfSpan::Allele,
            VcfSpan::Ref => gloci_core::VcfSpan::Ref,
        },
        ..ro
    };
    let settings = SortSettings {
        max_mem: a.sort.max_mem,
        tmp_dir: ctx.tmp_dir.clone(),
    };
    let sorted = sort_readers(readers, &mut chroms, &settings, &reopen)?;
    Ok((vec![sorted], chroms))
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

/// The VCF span to open a relation command's inputs with: `gloci sort` keys
/// VCF records by their REF span, so with `--sort` they are read that way and
/// the sorted stream is read back with the user's span.
pub fn open_span(sort: &SortArg, span: VcfSpan) -> gloci_core::VcfSpan {
    if sort.sort {
        return gloci_core::VcfSpan::Ref;
    }
    user_span(span)
}

fn user_span(span: VcfSpan) -> gloci_core::VcfSpan {
    match span {
        VcfSpan::Allele => gloci_core::VcfSpan::Allele,
        VcfSpan::Ref => gloci_core::VcfSpan::Ref,
    }
}

/// A's reader and the readers of the `-b` files.
pub type ARelation = (Box<dyn RecordReader>, Vec<Box<dyn RecordReader>>);

/// With `--sort`, sort A and then each B on its own (so the B files stay
/// distinct and hit order is file, then input order) and return the sorted
/// readers; without it, return them unchanged. The inputs must have been
/// opened with [`open_span`].
pub fn sort_relation_inputs(
    ctx: &Ctx,
    sort: &SortArg,
    genome: Option<&str>,
    chrom_order: &ChromOrderArg,
    span: VcfSpan,
    a: Box<dyn RecordReader>,
    b: Vec<Box<dyn RecordReader>>,
) -> Result<ARelation, IoError> {
    if !sort.sort {
        return Ok((a, b));
    }
    let mut chroms = {
        let all: Vec<&dyn RecordReader> = std::iter::once(a.as_ref())
            .chain(b.iter().map(|r| r.as_ref()))
            .collect();
        chrom_table_for(genome, &all, chrom_order.chrom_order)?
    };
    let reopen = ReaderOptions {
        threads: ctx.threads,
        format: None,
        vcf_span: user_span(span),
        invalid: InvalidPolicy::Fail,
    };
    let settings = SortSettings {
        max_mem: sort.max_mem,
        tmp_dir: ctx.tmp_dir.clone(),
    };
    let a = sort_readers(vec![a], &mut chroms, &settings, &reopen)?;
    let b = b
        .into_iter()
        .map(|r| sort_readers(vec![r], &mut chroms, &settings, &reopen))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((a, b))
}
