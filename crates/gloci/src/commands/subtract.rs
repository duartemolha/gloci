//! `gloci subtract`: remove from A the bases covered by B.
use crate::Ctx;
use crate::cli::SubtractArgs;
use crate::commands::{relation_hint, sort_aware, sort_settings, user_span};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::subtract::{SubtractOptions, subtract, validate};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, all_indexed, open_report, run_indexed};

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
    let reader_opts = |format| {
        sort_aware(
            &a.sort,
            ReaderOptions {
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
    let reader_a = open_reader(&a.a.a, &ro_a)?;
    let readers_b =
        a.b.b
            .iter()
            .map(|f| open_reader(f, &ro_b))
            .collect::<Result<Vec<_>, _>>()?;
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
    let all_paths: Vec<&str> = std::iter::once(a.a.a.as_str())
        .chain(a.b.b.iter().map(|s| s.as_str()))
        .collect();
    let use_indexed = all_indexed(&all_paths) && !a.u.unsorted && o.relate.sort.is_none();

    if a.sort.tmp_dir.is_some() && !use_indexed && o.relate.sort.is_none() {
        return Err(UsageError("--tmp-dir requires --sort or indexed inputs".into()).into());
    }

    if use_indexed {
        run_indexed(
            &a.a.a,
            &a.b.b,
            &ro_a,
            &ro_b,
            &o.relate,
            a.o.criteria().inclusive_points,
            ctx.threads,
            &ctx.tmp_dir,
            w.as_mut(),
            &|chunk_inputs, chunk_sink| {
                let summary = subtract(chunk_inputs, &o, chunk_sink)?;
                Ok((summary, ()))
            },
        )
        .map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    } else {
        let inputs = RelateInputs {
            a: reader_a,
            b: readers_b,
        };
        subtract(inputs, &o, w.as_mut()).map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    }
    w.finish()?;
    Ok(())
}
