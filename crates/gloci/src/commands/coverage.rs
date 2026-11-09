//! `gloci coverage`: how much of each A the B records cover.
use crate::Ctx;
use crate::cli::CoverageArgs;
use crate::commands::{relation_hint, sort_aware, sort_settings, user_span};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::coverage::{
    CoverageOptions, CoverageReport, CoverageTotals, column_names, coverage, validate,
    write_hist_totals,
};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, all_indexed, open_report, run_indexed};

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

pub fn run(ctx: &Ctx, a: &CoverageArgs) -> anyhow::Result<()> {
    let labels = a.b.labels().map_err(UsageError)?;
    let report = match (a.hist, a.per_base, a.mean) {
        (true, _, _) => CoverageReport::Hist,
        (_, true, _) => CoverageReport::PerBase,
        (_, _, true) => CoverageReport::Mean,
        _ => CoverageReport::Summary,
    };
    let o = CoverageOptions {
        report,
        relate: RelateOptions {
            criteria: a.o.criteria(),
            window: None,
            split: a.o.split,
            unsorted: a.u.unsorted,
            sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
            labels,
            genome: a.genome.genome.clone(),
            chrom_order: a.chrom_order.chrom_order,
            header: a.header.header,
        },
    };
    validate(&o, a.b.b.len())?;
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
        &ReportOutput::Tsv {
            appended: column_names(&o, a.b.b.len()),
        },
        &spec,
        &oo,
        &wo,
        a.header.header,
    )?;
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "stdout".to_string(), |p| p.display().to_string());

    let all_paths: Vec<&str> = std::iter::once(a.a.a.as_str())
        .chain(a.b.b.iter().map(|s| s.as_str()))
        .collect();
    let use_indexed = all_indexed(&all_paths) && !a.u.unsorted && o.relate.sort.is_none();

    let totals = if use_indexed {
        let chunk_results = run_indexed(
            &a.a.a,
            &a.b.b,
            &ro_a,
            &ro_b,
            &o.relate,
            a.o.criteria().inclusive_points,
            ctx.threads,
            &ctx.tmp_dir,
            w.as_mut(),
            &|chunk_inputs, chunk_sink| coverage(chunk_inputs, &o, chunk_sink),
        )
        .map_err(|e| {
            e.with_output_name(&dest)
                .with_order_hint(relation_hint(&a.sort))
        })?;
        let mut combined_totals = CoverageTotals::default();
        for (_, t) in chunk_results {
            combined_totals.merge(&t);
        }
        combined_totals
    } else {
        let inputs = RelateInputs {
            a: reader_a,
            b: readers_b,
        };
        let (_, totals) = coverage(inputs, &o, w.as_mut()).map_err(|e| {
            e.with_output_name(&dest)
                .with_order_hint(relation_hint(&a.sort))
        })?;
        totals
    };

    if report == CoverageReport::Hist {
        write_hist_totals(&totals, w.as_mut()).map_err(|e| e.with_output_name(&dest))?;
    }
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
