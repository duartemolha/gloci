//! `gloci flagstat`: samtools' flag counts of a BAM, CRAM or SAM (D128).
use crate::Ctx;
use crate::cli::FlagstatArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_aln_source, open_output,
};
use gloci_ops::{FlagstatFormat, flagstat, render_flagstat};
use std::io::Write;

/// The extension of `-o` that names the layout, lower-cased, after any compression suffix
/// (`x.json.gz` is `json`).
fn layout_ext(a: &FlagstatArgs) -> Option<String> {
    let p = a.output.as_deref()?;
    let name = p.file_name()?.to_string_lossy().to_ascii_lowercase();
    let name = if super::compressed_name(p) {
        name.rsplit_once('.')
            .map_or(name.as_str(), |(stem, _)| stem)
    } else {
        name.as_str()
    };
    name.rsplit_once('.').map(|(_, ext)| ext.to_string())
}

/// `-O`, else the extension of `-o` (`.json`, `.tsv`, also compressed), else text.
fn format_of(a: &FlagstatArgs) -> Result<FlagstatFormat, UsageError> {
    if let Some(s) = &a.output_format {
        return FlagstatFormat::from_name(s).ok_or_else(|| {
            UsageError(format!(
                "unknown output format '{s}' for -O; use one of: {}",
                FlagstatFormat::NAMES
            ))
        });
    }
    let ext = layout_ext(a);
    Ok(match ext.as_deref() {
        Some("json") => FlagstatFormat::Json,
        Some("tsv") => FlagstatFormat::Tsv,
        _ => FlagstatFormat::Text,
    })
}

pub fn run(ctx: &Ctx, a: &FlagstatArgs) -> anyhow::Result<()> {
    let format = format_of(a)?;
    if let Some(ext) = layout_ext(a)
        && !matches!(ext.as_str(), "txt" | "tsv" | "json")
        && ext.parse::<gloci_io::Format>().is_ok()
    {
        return Err(UsageError(format!(
            "flagstat writes text, tsv or json; -o *.{ext} names another format \
             (use -o x.txt, x.tsv or x.json, or -O)"
        ))
        .into());
    }
    // Text of flagstat's own layout, BGZF-compressed when -o ends in .gz or .bgz.
    let spec = super::own_tsv_output_spec(a.output.as_deref())?;
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
