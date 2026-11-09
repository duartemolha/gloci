//! `gloci groupby`: aggregate runs of lines with equal key columns.
use crate::Ctx;
use crate::cli::GroupbyArgs;
use gloci_io::{
    Format, InputOptions, OutFormat, OutputOptions, UsageError, detect_format, open_input,
    open_output, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::groupby::{GroupbyOptions, groupby, parse_keys};
use std::io::{BufRead, IsTerminal};

pub fn run(ctx: &Ctx, a: &GroupbyArgs) -> anyhow::Result<()> {
    let keys = parse_keys(&a.group)?;
    if a.out.write_index {
        return Err(UsageError(
            "--write-index does not apply to groupby: its output has no coordinates to index; \
             drop --write-index"
                .into(),
        )
        .into());
    }
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Tsv,
        spec.bgzf,
    )?;
    // Any tab-separated text is accepted: the format is fixed so that text that
    // does not look like BED passes detection, then peeked at to refuse VCF/BCF.
    let mut input = open_input(
        &a.file,
        &InputOptions {
            threads: ctx.threads,
            format: Some(Format::Bed),
        },
    )?;
    let head = input
        .reader
        .fill_buf()
        .map_err(|source| gloci_io::IoError::Read {
            path: a.file.clone(),
            source,
        })?;
    if let Ok(f @ (Format::Vcf | Format::Bcf)) = detect_format(head) {
        return Err(UsageError(format!(
            "{} is {} input, and groupby groups tab-separated columns of plain text; \
             pull the fields out with `gloci vcf extract` and group that output",
            a.file,
            f.name()
        ))
        .into());
    }
    let o = GroupbyOptions {
        keys,
        aggs: a.columns.clone(),
        header: a.header,
        empty: a.empty.clone(),
        delim: a.delim.clone(),
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: false,
    };
    let mut out = open_output(&spec, OutFormat::Tsv, &oo).map_err(|e| e.with_output_name(&dest))?;
    groupby(&mut input.reader, &input.name, &o, &mut out, &dest)?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
