//! `gloci flagstat`: samtools' flag counts of a BAM, CRAM or SAM (D128).
use crate::Ctx;
use crate::cli::FlagstatArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_aln_source, open_output, parse_output_spec,
};
use gloci_ops::{FlagstatFormat, flagstat, render_flagstat};
use std::io::Write;

/// `-O`, else the extension of `-o` (`.json`, `.tsv`), else text.
fn format_of(a: &FlagstatArgs) -> Result<FlagstatFormat, UsageError> {
    if let Some(s) = &a.output_format {
        return FlagstatFormat::from_name(s).ok_or_else(|| {
            UsageError(format!(
                "unknown output format '{s}' for -O; use one of: {}",
                FlagstatFormat::NAMES
            ))
        });
    }
    let ext = a
        .output
        .as_ref()
        .and_then(|p| p.extension())
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase);
    Ok(match ext.as_deref() {
        Some("json") => FlagstatFormat::Json,
        Some("tsv") => FlagstatFormat::Tsv,
        _ => FlagstatFormat::Text,
    })
}

pub fn run(ctx: &Ctx, a: &FlagstatArgs) -> anyhow::Result<()> {
    let format = format_of(a)?;
    // Plain text whatever -o is called, so the output spec is fixed (as for `dict`).
    let spec = parse_output_spec(a.output.as_deref(), Some("tsv"))?;
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
    let src = open_aln_source(&a.file, &ro)?;
    // Counted whole before anything is opened for writing, so a failing scan leaves no output.
    let stats = flagstat(src)?;
    let text = render_flagstat(&stats, format);
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    out.write_all(text.as_bytes())
        .map_err(|source| gloci_io::IoError::Write {
            path: dest.clone(),
            source,
        })?;
    out.finish()?;
    Ok(())
}
