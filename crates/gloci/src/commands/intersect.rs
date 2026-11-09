//! `gloci intersect`: A's records by their overlaps with B.
use crate::Ctx;
use crate::cli::{IntersectArgs, IntersectMatch, IntersectMode};
use crate::commands::{
    dest_name, relation_files, relation_hint, run_options, sort_aware, sort_settings, user_span,
};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions,
    parse_output_spec,
};
use gloci_ops::intersect::{
    IntersectOptions, Match, Mode, intersect_files, intersect_output, validate,
};
use gloci_ops::{RelateOptions, open_report};

fn usage(msg: impl Into<String>) -> anyhow::Error {
    UsageError(msg.into()).into()
}

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

pub fn run(ctx: &Ctx, a: &IntersectArgs) -> anyhow::Result<()> {
    if a.mode != IntersectMode::Pairs && (a.overlap_bases || a.keep_unmatched) {
        let flag = if a.overlap_bases {
            "--overlap-bases"
        } else {
            "--keep-unmatched"
        };
        return Err(usage(format!(
            "{flag} only applies to --mode pairs; add --mode pairs"
        )));
    }
    let labels = a.b.labels().map_err(usage)?;
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

    let mode = match a.mode {
        IntersectMode::Any => Mode::Any,
        IntersectMode::NoHit => Mode::None,
        IntersectMode::Pairs => Mode::Pairs {
            overlap_bases: a.overlap_bases,
            keep_unmatched: a.keep_unmatched,
        },
        IntersectMode::Clip => Mode::Clip,
        IntersectMode::Count => Mode::Count,
    };
    let o = IntersectOptions {
        mode,
        matching: match a.matching {
            IntersectMatch::Overlap => Match::Overlap,
            IntersectMatch::Alleles => Match::Alleles,
        },
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
    validate(files.a().format(), &files.b_formats(), &o)?;
    let kind = intersect_output(&mut files, &o)?;
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
    let dest = dest_name(&spec);
    intersect_files(files, &o, w.as_mut(), &run_options(ctx)).map_err(|e| {
        e.with_output_name(&dest)
            .with_order_hint(relation_hint(&a.sort))
    })?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
