//! `gloci intersect`: A's records by their overlaps with B.
use crate::Ctx;
use crate::cli::{IntersectArgs, IntersectMatch, IntersectMode};
use crate::commands::{
    dest_name, relation_files, relation_hint, restrict_a, run_options, sort_aware, sort_settings,
    user_span,
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
    let chrom_map = a.chrom_map.load()?;
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
    // Reads are kept whole only when they are written: A in `any` and `none` (D115).
    let write_a = matches!(mode, Mode::Any | Mode::None);
    let accept = |whole: bool| -> anyhow::Result<gloci_io::Accept> {
        Ok(gloci_io::Accept {
            alignments: Some(a.align.options(whole)?),
            ..gloci_io::Accept::default()
        })
    };
    let reader_opts = |format, accept| {
        sort_aware(
            &a.sort,
            ReaderOptions {
                accept,
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
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?, accept(write_a)?);
    let ro_b = reader_opts(parse_format(a.b.b_format.as_deref())?, accept(false)?);
    let mut files = relation_files(&a.a.a, &a.b.b, ro_a, ro_b)?;
    files = restrict_a(
        files,
        &a.region,
        a.o.inclusive_points,
        a.u.unsorted || a.sort.sort,
    )?;

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
            split_deletions: a.split_deletions,
            unsorted: a.u.unsorted,
            sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
            labels,
            genome: a.genome.genome.clone(),
            chrom_order: a.chrom_order.chrom_order,
            header: a.header.header,
        },
    };
    validate(files.a().format(), &files.b_formats(), &o)?;
    files = files.with_sorted_a_header(&o.relate)?;
    // SAM lines are written as read; --chrom-map renames the header, so the two would
    // disagree. (BAM records name references by number, so BAM is fine.)
    if write_a && chrom_map.is_some() && files.a().format() == Format::Sam {
        return Err(usage(
            "--chrom-map cannot be used when writing a SAM A as alignments (--mode any or none): \
the reads keep the chromosome names they have, so they would not match the renamed header; \
convert A to BAM first (gloci view a.sam -O bam -o a.bam) and pass that, or use --mode pairs",
        ));
    }
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
