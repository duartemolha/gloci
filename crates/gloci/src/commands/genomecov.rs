//! `gloci genomecov`: depth of coverage of a genome by BED features.
use super::{InputArgs, SORT_FLAG_HINT, open_inputs};
use crate::Ctx;
use crate::cli::GenomecovArgs;
use gloci_core::Strand;
use gloci_io::{
    Format, Header, OutFormat, OutputOptions, TsvHeader, WriterOptions, open_output,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::genomecov::{GenomecovMode, GenomecovOptions, Piece, genomecov, validate};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &GenomecovArgs) -> anyhow::Result<()> {
    let o = GenomecovOptions {
        mode: if a.bedgraph {
            GenomecovMode::BedGraph
        } else if a.bga {
            GenomecovMode::BedGraphAll
        } else if a.per_base {
            GenomecovMode::PerBase
        } else {
            GenomecovMode::Hist
        },
        scale: a.scale,
        strand: a.strand.as_deref().map(|s| {
            if s == "+" {
                Strand::Plus
            } else {
                Strand::Minus
            }
        }),
        split: a.split,
        piece: if a.five_prime {
            Piece::FivePrime
        } else if a.three_prime {
            Piece::ThreePrime
        } else {
            Piece::Whole
        },
        max_depth: a.max_depth,
    };
    validate(&o)?;
    let genome = resolve_genome(&a.genome)?.0;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let files = [a.file.clone()];
    let (mut readers, mut chroms) = open_inputs(
        ctx,
        &InputArgs {
            files: &files,
            format: a.format.format.as_deref(),
            vcf_span: a.vcf_span.vcf_span,
            invalid: a.invalid.skip_invalid,
            genome: Some(&a.genome),
            chrom_order: &a.chrom_order,
            sort: &a.sort,
            chrom_map: &a.chrom_map,
        },
        &|_| Ok(()),
    )?;
    anyhow::ensure!(readers.len() == 1, "genomecov takes exactly one input");
    let reader = readers.remove(0);

    // The rows are new TSV/BED lines whatever the input format.
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
    genomecov(reader, &mut chroms, &genome, &o, w.as_mut()).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(if a.sort.sort { "" } else { SORT_FLAG_HINT })
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
