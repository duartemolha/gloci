//! `gloci closest`: the nearest B to each A.
use crate::Ctx;
use crate::cli::{ClosestArgs, ClosestTies, SignedFrame};
use crate::commands::{
    relation_files, relation_hint, run_options, sort_aware, sort_settings, user_span,
};
use gloci_core::{OverlapCriteria, SignFrame, StrandRule};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions,
    parse_output_spec,
};
use gloci_ops::closest::{ClosestOptions, Ties, closest_files, closest_output, validate};
use gloci_ops::{RelateOptions, open_report};

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

pub fn run(ctx: &Ctx, a: &ClosestArgs) -> anyhow::Result<()> {
    let labels = a.b.labels().map_err(UsageError)?;
    let o = ClosestOptions {
        k: usize::try_from(a.k).unwrap_or(usize::MAX),
        ties: match a.ties {
            ClosestTies::All => Ties::All,
            ClosestTies::First => Ties::First,
            ClosestTies::Last => Ties::Last,
        },
        distance: a.distance,
        signed: a.signed.map(|f| match f {
            SignedFrame::Ref => SignFrame::Ref,
            SignedFrame::A => SignFrame::A,
            SignedFrame::B => SignFrame::B,
        }),
        ignore_overlaps: a.ignore_overlaps,
        ignore_upstream: a.ignore_upstream,
        ignore_downstream: a.ignore_downstream,
        max_distance: a.max_distance,
        each_side: a.each_side,
        per_file: a.per_file,
        relate: RelateOptions {
            criteria: OverlapCriteria {
                min_frac_a: None,
                min_frac_b: None,
                reciprocal: false,
                either: false,
                min_bases: None,
                strand: match (a.same_strand, a.opposite_strand) {
                    (true, _) => StrandRule::Same,
                    (_, true) => StrandRule::Opposite,
                    _ => StrandRule::Any,
                },
                inclusive_points: a.inclusive_points,
            },
            window: None,
            split: false,
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
    let mut files = relation_files(&a.a.a, &a.b.b, ro_a, ro_b)?;

    let kind = closest_output(&mut files, &o)?;
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
    let mut w = open_report(files.a(), &kind, &spec, &oo, &wo, a.header.header)?;
    closest_files(files, &o, w.as_mut(), &run_options(ctx))
        .map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    w.finish()?;
    Ok(())
}
