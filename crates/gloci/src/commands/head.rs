//! `gloci head`: the header and the first records of any input as text (D134).
use crate::Ctx;
use crate::cli::HeadArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, OutputSpec, ReaderOptions,
    UsageError, open_output, plain_text_output_spec,
};
use gloci_ops::{HeadOptions, head};
use std::path::Path;

/// head writes text whatever the input: a name ending in `.gz` or `.bgz` is BGZF-compressed,
/// a binary format's name (`.bam`, `.bcf`, `.cram`) is refused, and any other name gets the
/// text as it is (`.sam`, `.vcf`, `.txt`).
fn output_spec(output: Option<&Path>) -> Result<OutputSpec, UsageError> {
    if let Some(ext) = output
        .and_then(Path::extension)
        .map(|e| e.to_string_lossy().to_ascii_lowercase())
        && matches!(ext.as_str(), "bam" | "bcf" | "cram")
    {
        return Err(UsageError(format!(
            "head writes text; -o *.{ext} names a binary format (use .sam, .vcf, .txt or \
             another text name, or gloci view to write {ext})"
        )));
    }
    let mut spec = plain_text_output_spec(output);
    spec.bgzf = output.is_some_and(super::compressed_name);
    Ok(spec)
}

pub fn run(ctx: &Ctx, a: &HeadArgs) -> anyhow::Result<()> {
    let spec = output_spec(a.output.as_deref())?;
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
    let o = HeadOptions {
        headers: a.headers,
        records: a.records,
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    // On a failure, what was read before it stays written to stdout (dropping the output
    // flushes it), as `head` does with a cut file; a FILE is removed, as every output is.
    head(&a.file, &o, &ro, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
