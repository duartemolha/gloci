//! `gloci subtract`: remove from A the bases covered by B.
use crate::Ctx;
use crate::cli::SubtractArgs;
use crate::commands::{
    dest_name, relation_files, relation_hint, restrict_a, run_options, sort_aware, sort_settings,
    user_span,
};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions,
    parse_output_spec,
};
use gloci_ops::subtract::{SubtractOptions, subtract_files, validate};
use gloci_ops::{RelateOptions, ReportOutput, open_report};

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

pub fn run(ctx: &Ctx, a: &SubtractArgs) -> anyhow::Result<()> {
    let labels = a.b.labels().map_err(UsageError)?;
    let o = SubtractOptions {
        whole: a.whole,
        sum_overlaps: a.sum_overlaps,
        relate: RelateOptions {
            criteria: a.o.criteria(),
            window: a.w.window(),
            split: a.o.split,
            unsorted: a.u.unsorted,
            sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
            labels,
            genome: a.genome.genome.clone(),
            chrom_order: a.chrom_order.chrom_order,
            header: a.header.header,
        },
    };
    validate(&o)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let chrom_map = a.chrom_map.load()?;
    let reader_opts = |format| {
        sort_aware(
            &a.sort,
            ReaderOptions {
                chrom_map: chrom_map.clone(),
                threads: ctx.threads,
                format,
                vcf_span: user_span(a.vcf_span.vcf_span),
                invalid: if a.invalid.skip_invalid {
                    InvalidPolicy::Skip
                } else {
                    InvalidPolicy::Fail
                },
            },
        )
    };
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?);
    let ro_b = reader_opts(parse_format(a.b.b_format.as_deref())?);
    let files = relation_files(&a.a.a, &a.b.b, ro_a, ro_b)?;
    let files = restrict_a(
        files,
        &a.region,
        a.o.inclusive_points,
        a.u.unsorted || a.sort.sort,
    )?;
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
        files.a(),
        &ReportOutput::SameAsA,
        &spec,
        &oo,
        &wo,
        a.header.header,
    )?;
    let dest = dest_name(&spec);
    subtract_files(files, &o, w.as_mut(), &run_options(ctx)).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(relation_hint(&a.sort))
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
