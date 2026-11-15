//! `gloci igv`: an IGV batch script with a snapshot of each feature. Also holds the plumbing
//! `gloci links` shares.
use crate::Ctx;
use crate::cli::IgvArgs;
use gloci_core::ChromTable;
use gloci_io::{
    Accept, Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, RecordReader,
    UsageError, open_output, open_reader, plain_text_output_spec,
};
use gloci_ops::{IgvOptions, igv};
use std::path::Path;

/// The input as BED, GFF or VCF records, read in order.
pub fn open_features(
    ctx: &Ctx,
    file: &str,
    format: Option<&str>,
) -> anyhow::Result<Box<dyn RecordReader>> {
    let ro = ReaderOptions {
        accept: Accept::default(),
        chrom_map: None,
        threads: ctx.threads,
        format: match format {
            None => None,
            Some(f) => Some(f.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
    };
    Ok(open_reader(file, &ro)?)
}

/// Write what `write` makes to `output` (stdout when none) as text; a name ending in `.gz` or
/// `.bgz` is BGZF-compressed. A failure removes a file; what was written before it stays on
/// stdout.
pub fn write_text(
    ctx: &Ctx,
    output: Option<&Path>,
    write: impl FnOnce(&mut dyn std::io::Write) -> Result<(), gloci_io::IoError>,
) -> anyhow::Result<()> {
    let mut spec = plain_text_output_spec(output);
    spec.bgzf = output.is_some_and(super::compressed_name);
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    write(&mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}

pub fn run(ctx: &Ctx, a: &IgvArgs) -> anyhow::Result<()> {
    let reader = open_features(ctx, &a.file, a.format.format.as_deref())?;
    let o = IgvOptions {
        path: a.path.clone(),
        session: a.session.clone(),
        sort: a.sort.clone(),
        collapse: a.collapse,
        name: a.name,
        slop: a.slop,
        image: a.image.clone(),
    };
    let mut chroms = ChromTable::first_appearance();
    write_text(ctx, a.output.as_deref(), |w| igv(reader, &mut chroms, &o, w))
}
