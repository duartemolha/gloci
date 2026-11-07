//! `gloci view`: region queries, streaming, VCF/BCF conversion and header options.
use crate::Ctx;
use crate::cli::{VcfSpan, ViewArgs};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions,
    parse_output_spec,
};
use gloci_ops::{ViewOptions, run_view};

pub fn run(ctx: &Ctx, a: &ViewArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let format = match a.format.format.as_deref() {
        None => None,
        Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
    };
    let ro = ReaderOptions {
        threads: ctx.threads,
        format,
        vcf_span: match a.vcf_span.vcf_span {
            VcfSpan::Allele => gloci_core::VcfSpan::Allele,
            VcfSpan::Ref => gloci_core::VcfSpan::Ref,
        },
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let wo = WriterOptions {
        write_header: !a.no_header,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let vo = ViewOptions {
        header_only: a.header_only,
        no_header: a.no_header,
        regions: a.regions.iter().chain(&a.region.region).cloned().collect(),
        regions_file: a
            .region
            .regions_file
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        inclusive_points: a.inclusive_points,
        genome: a.genome.genome.clone(),
        chrom_order: a.chrom_order.chrom_order,
        assume_sorted: a.assume_sorted,
    };
    run_view(&a.file, &vo, &ro, &spec, &oo, &wo)?;
    Ok(())
}
