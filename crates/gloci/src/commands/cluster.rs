//! `gloci cluster` and the shared plumbing of the commands that append one
//! column to each record of a single sorted input (`spacing` too).
use super::{InputArgs, SORT_FLAG_HINT, dest_name, open_inputs};
use crate::Ctx;
use crate::cli::{
    ChromMapArg, ChromOrderArg, ClusterArgs, FormatArg, GenomeArg, HeaderArg, InvalidArg,
    OutputArgs, SortArg, VcfSpanArg,
};
use gloci_core::ChromTable;
use gloci_io::{
    IoError, OutputOptions, RecordReader, RecordWriter, WriterOptions, parse_output_spec,
};
use gloci_ops::cluster::{ClusterOptions, cluster, validate};
use gloci_ops::{ReportOutput, open_report};

/// The options every append-one-column command takes.
pub struct AppendArgs<'a> {
    pub file: &'a str,
    pub genome: &'a GenomeArg,
    pub chrom_order: &'a ChromOrderArg,
    pub vcf_span: &'a VcfSpanArg,
    pub chrom_map: &'a ChromMapArg,
    pub format: &'a FormatArg,
    pub header: &'a HeaderArg,
    pub invalid: &'a InvalidArg,
    pub sort: &'a SortArg,
    pub out: &'a OutputArgs,
}

/// Open the one input, then write its records with the column `name` added
/// by `op`.
pub fn run_appending(
    ctx: &Ctx,
    a: &AppendArgs,
    name: &str,
    op: impl FnOnce(
        Box<dyn RecordReader>,
        &mut ChromTable,
        &mut dyn RecordWriter,
    ) -> Result<(), IoError>,
) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let files = [a.file.to_string()];
    let (mut readers, mut chroms) = open_inputs(
        ctx,
        &InputArgs {
            files: &files,
            format: a.format.format.as_deref(),
            vcf_span: a.vcf_span.vcf_span,
            invalid: a.invalid.skip_invalid,
            genome: a.genome.genome.as_deref(),
            chrom_order: a.chrom_order,
            sort: a.sort,
            chrom_map: a.chrom_map,
            alignments_planned: false,
        },
        &|formats| validate(formats[0]),
    )?;
    let reader = readers.remove(0);
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let wo = WriterOptions {
        write_header: a.header.header,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let mut w = open_report(
        reader.as_ref(),
        &ReportOutput::Tsv {
            appended: vec![name.to_string()],
        },
        &spec,
        &oo,
        &wo,
        a.header.header,
    )?;
    let dest = dest_name(&spec);
    op(reader, &mut chroms, w.as_mut()).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(if a.sort.sort { "" } else { SORT_FLAG_HINT })
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

pub fn run(ctx: &Ctx, a: &ClusterArgs) -> anyhow::Result<()> {
    let o = ClusterOptions {
        max_gap: a.max_gap,
        per_strand: a.strand,
    };
    run_appending(
        ctx,
        &AppendArgs {
            file: &a.file,
            genome: &a.genome,
            chrom_order: &a.chrom_order,
            vcf_span: &a.vcf_span,
            chrom_map: &a.chrom_map,
            format: &a.format,
            header: &a.header,
            invalid: &a.invalid,
            sort: &a.sort,
            out: &a.out,
        },
        "cluster",
        |reader, chroms, w| cluster(reader, chroms, &o, w),
    )
}
