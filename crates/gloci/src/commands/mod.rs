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
