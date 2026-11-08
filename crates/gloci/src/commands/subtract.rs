//! `gloci subtract`: remove from A the bases covered by B.
use crate::Ctx;
use crate::cli::SubtractArgs;
use crate::commands::{open_span, relation_hint, sort_relation_inputs};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::subtract::{SubtractOptions, subtract, validate};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, open_report};

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
            labels,
            genome: a.genome.genome.clone(),
            chrom_order: a.chrom_order.chrom_order,
            header: a.header.header,
        },
    };
    validate(&o)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let reader_opts = |format| ReaderOptions {
        threads: ctx.threads,
        format,
        vcf_span: open_span(&a.sort, a.vcf_span.vcf_span),
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?);
    let ro_b = reader_opts(parse_format(a.b.b_format.as_deref())?);
    let reader_a = open_reader(&a.a.a, &ro_a)?;
    let readers_b =
        a.b.b
            .iter()
            .map(|f| open_reader(f, &ro_b))
            .collect::<Result<Vec<_>, _>>()?;
    let (reader_a, readers_b) = sort_relation_inputs(
        ctx,
        &a.sort,
        a.genome.genome.as_deref(),
        &a.chrom_order,
        a.vcf_span.vcf_span,
        reader_a,
        readers_b,
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
        reader_a.as_ref(),
        &ReportOutput::SameAsA,
        &spec,
        &oo,
        &wo,
        a.header.header,
    )?;
    let inputs = RelateInputs {
        a: reader_a,
        b: readers_b,
    };
    subtract(inputs, &o, w.as_mut()).map_err(|e| {
        e.with_order_hint(relation_hint(&a.sort))
            .with_sorted_input_loc()
    })?;
    w.finish()?;
    Ok(())
}
