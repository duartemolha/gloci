//! `gloci map`: aggregate B columns over each A.
use crate::Ctx;
use crate::cli::MapArgs;
use crate::commands::{relation_hint, sort_aware, sort_settings, user_span};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::agg::parse_aggs;
use gloci_ops::map::{MapOptions, map, validate};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, all_indexed, open_report, run_indexed};

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

pub fn run(ctx: &Ctx, a: &MapArgs) -> anyhow::Result<()> {
    let aggs = parse_aggs(&a.columns)?;
    let labels = a.b.labels().map_err(UsageError)?;
    let o = MapOptions {
        aggs,
        skip_unmatched: a.skip_unmatched,
        empty: a.empty.clone(),
        delim: a.delim.clone(),
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
    let b_formats: Vec<Format> = readers_b.iter().map(|r| r.format()).collect();
    validate(&o, &b_formats, &a.b.b)?;

    let names = o.aggs.iter().map(|s| s.label.clone()).collect();
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
        &ReportOutput::Tsv { appended: names },
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
                let summary = map(chunk_inputs, &o, chunk_sink)?;
                Ok((summary, ()))
            },
        )
        .map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    } else {
        let inputs = RelateInputs {
            a: reader_a,
            b: readers_b,
        };
        map(inputs, &o, w.as_mut()).map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    }
    w.finish()?;
    Ok(())
}
