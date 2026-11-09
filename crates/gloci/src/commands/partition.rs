//! `gloci partition`: split sorted inputs at every boundary.
use super::{InputArgs, SORT_FLAG_HINT, open_inputs_each};
use crate::Ctx;
use crate::cli::PartitionArgs;
use gloci_io::{
    Format, Header, OutFormat, OutputOptions, TsvHeader, WriterOptions, open_output,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::partition::{PartitionOptions, partition, validate};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &PartitionArgs) -> anyhow::Result<()> {
    let names = if a.names.is_empty() {
        a.files.iter().map(|f| label_of(f)).collect()
    } else {
        a.names.clone()
    };
    let o = PartitionOptions {
        names,
        min_files: a.min_files,
        max_files: a.max_files,
        all: a.all,
        merge_adjacent: a.merge_adjacent,
        membership: !a.no_membership,
        bedgraph: a.bedgraph,
        empty: a.empty.clone(),
    };
    let genome = match a.genome.genome.as_deref() {
        Some(g) if o.empty.is_some() => Some(resolve_genome(g)?.0),
        _ => None,
    };
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let (readers, mut chroms) = open_inputs_each(
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
    partition(readers, &mut chroms, genome.as_ref(), &o, w.as_mut()).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(if a.sort.sort { "" } else { SORT_FLAG_HINT })
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

/// The default label of an input: its file name without directories (`stdin`
/// for `-`).
fn label_of(file: &str) -> String {
    if file == "-" {
        return "stdin".into();
    }
    std::path::Path::new(file)
        .file_name()
        .map_or_else(|| file.to_string(), |n| n.to_string_lossy().into_owned())
}
