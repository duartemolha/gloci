//! `gloci complement`: the bases no feature covers.
use super::{InputArgs, SORT_FLAG_HINT, open_inputs};
use crate::Ctx;
use crate::cli::ComplementArgs;
use gloci_io::{
    Format, Header, OutFormat, OutputOptions, TsvHeader, WriterOptions, open_output,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::complement::complement;
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &ComplementArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let genome = match a.genome.genome.as_deref() {
        Some(g) => Some(resolve_genome(g)?.0),
        None => None,
    };
    let (readers, mut chroms) = open_inputs(
        ctx,
        &InputArgs {
            files: &a.files,
            format: a.format.format.as_deref(),
            vcf_span: a.vcf_span.vcf_span,
            invalid: a.invalid.skip_invalid,
            genome: a.genome.genome.as_deref(),
            chrom_order: &a.chrom_order,
            sort: &a.sort,
        },
        &|_| Ok(()),
    )?;

    let to = spec.format.unwrap_or(OutFormat::Bed);
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let out = open_output(&spec, to, &oo)?;
    let mut w = if to == OutFormat::Tsv {
        open_tsv_writer(out, Format::Bed, &[], TsvHeader::None, oo.write_index)
    } else {
        let wo = WriterOptions {
            write_header: false,
            command_line: Some(ctx.command_line.clone()),
            version: ctx.version,
        };
        open_record_writer(out, Format::Bed, to, &Header::None, &wo)
    }
    .map_err(|e| e.with_output_name(&dest))?;
    complement(readers, &mut chroms, genome.as_ref(), w.as_mut()).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(if a.sort.sort { "" } else { SORT_FLAG_HINT })
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
