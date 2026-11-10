//! `gloci genomecov`: depth of coverage of a genome by BED features.
use super::{InputArgs, SORT_FLAG_HINT, open_inputs, reader_options};
use crate::Ctx;
use crate::cli::GenomecovArgs;
use gloci_core::{Genome, Strand};
use gloci_io::{
    Format, Header, OutFormat, OutputOptions, TsvHeader, UsageError, WriterOptions, open_output,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::genomecov::{
    ChunkRun, GenomecovMode, GenomecovOptions, Piece, genomecov_files, validate,
};
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
        split_deletions: a.split_deletions,
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
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let files = [a.file.clone()];
    // Reads are slim: genomecov never writes one (D115).
    let input = InputArgs {
        files: &files,
        format: a.format.format.as_deref(),
        vcf_span: a.vcf_span.vcf_span,
        invalid: a.invalid.skip_invalid,
        genome: a.genome.as_deref(),
        chrom_order: &a.chrom_order,
        sort: &a.sort,
        chrom_map: &a.chrom_map,
        alignments: Some(a.align.options(false)?),
    };
    let (mut readers, mut chroms) = open_inputs(ctx, &input, &|formats| {
        if a.split_deletions && !formats.iter().any(|f| f.is_alignment()) {
            return Err(UsageError(
                "--split-deletions cuts alignments at deletions, and the input is not BAM \
or SAM; drop it, or use --split for BED12 blocks"
                    .into(),
            ));
        }
        Ok(())
    })?;
    anyhow::ensure!(readers.len() == 1, "genomecov takes exactly one input");
    let reader = readers.remove(0);
    let genome = genome_of(a.genome.as_deref(), reader.as_ref())?;
    // An indexed BAM runs per chromosome; `--sort` streams the sorted copy.
    let chunked = (!a.sort.sort && a.file != "-").then(|| -> anyhow::Result<ChunkRun> {
        Ok(ChunkRun {
            path: a.file.clone().into(),
            ro: reader_options(ctx, &input)?,
            threads: ctx.threads,
            tmp_dir: ctx.tmp_dir.clone(),
        })
    });
    let chunked = chunked.transpose()?;

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
    genomecov_files(
        reader,
        chunked.as_ref(),
        &mut chroms,
        &genome,
        &o,
        w.as_mut(),
    )
    .map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(if a.sort.sort { "" } else { SORT_FLAG_HINT })
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

/// The genome to cover: `-g`, else the lengths in an alignment's header.
fn genome_of(arg: Option<&str>, reader: &dyn gloci_io::RecordReader) -> anyhow::Result<Genome> {
    if let Some(g) = arg {
        return Ok(resolve_genome(g)?.0);
    }
    match reader.header() {
        Header::Alignment(h) if !h.refs.is_empty() => {
            Genome::new(h.refs.iter().map(|(n, l)| (n.clone(), Some(*l))).collect()).map_err(|e| {
                UsageError(format!(
                    "{}: the header's references do not make a genome: {e}",
                    reader.name()
                ))
                .into()
            })
        }
        Header::Alignment(_) => Err(UsageError(format!(
            "{} has no @SQ lines, so there are no chromosome lengths; \
pass -g with a built-in genome such as hg38 or a .fai file",
            reader.name()
        ))
        .into()),
        _ => Err(UsageError(
            "-g is required unless the input is BAM or SAM: give a built-in genome such as hg38 \
(see gloci genomes), or a .fai, .genome or .dict file"
                .into(),
        )
        .into()),
    }
}
