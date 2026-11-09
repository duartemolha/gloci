//! `gloci merge`: cluster the footprints of sorted inputs.
use super::{InputArgs, SORT_FLAG_HINT, open_inputs};
use crate::Ctx;
use crate::cli::MergeArgs;
use gloci_io::{
    Format, Header, OutFormat, OutputOptions, TsvHeader, WriterOptions, open_output,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::agg::parse_aggs;
use gloci_ops::merge::{MergeOptions, merge, validate};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &MergeArgs) -> anyhow::Result<()> {
    let aggs = parse_aggs(&a.columns)?;
    let o = MergeOptions {
        max_gap: a.max_gap,
        per_strand: a.strand,
        count: a.count,
        aggs,
        empty: a.empty.clone(),
        delim: a.delim.clone(),
    };
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
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
            chrom_map: &a.chrom_map,
        },
        &|formats| validate(formats, &o),
    )?;

    // The rows are new BED lines whatever the input format.
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
    merge(readers, &mut chroms, &o, w.as_mut()).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(if a.sort.sort { "" } else { SORT_FLAG_HINT })
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
