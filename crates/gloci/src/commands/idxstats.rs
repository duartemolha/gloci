//! `gloci idxstats`: reads per reference of a BAM, CRAM or SAM (D127).
use crate::Ctx;
use crate::cli::IdxstatsArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions,
    open_output,
};
use gloci_ops::idxstats;
use std::io::Write;

pub fn run(ctx: &Ctx, a: &IdxstatsArgs) -> anyhow::Result<()> {
    // TSV of idxstats' own layout, BGZF-compressed when -o ends in .gz or .bgz.
    let spec = super::own_tsv_output_spec(a.output.as_deref(), "idxstats")?;
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                reference: a.reference.clone(),
                ..AlignOptions::default()
            }),
            convert_only: false,
        },
        chrom_map: None,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    // Counted whole before anything is opened for writing, so a failing scan leaves no output.
    let mut text = Vec::new();
    idxstats(&a.file, &ro, &mut text)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    out.write_all(&text)
        .map_err(|source| IoError::Write { path: dest, source })?;
    out.finish()?;
    Ok(())
}
