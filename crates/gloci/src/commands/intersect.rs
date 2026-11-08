//! `gloci intersect`: A's records by their overlaps with B.
use crate::Ctx;
use crate::cli::{IntersectArgs, IntersectMatch, IntersectMode, VcfSpan};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::intersect::{IntersectOptions, Match, Mode, intersect, validate};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, open_report};

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
    let reader_opts = |format| ReaderOptions {
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
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?);
    let ro_b = reader_opts(parse_format(a.b.b_format.as_deref())?);
    let reader_a = open_reader(&a.a.a, &ro_a)?;
    let readers_b =
        a.b.b
            .iter()
            .map(|f| open_reader(f, &ro_b))
            .collect::<Result<Vec<_>, _>>()?;

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
            labels,
            genome: a.genome.genome.clone(),
            chrom_order: a.chrom_order.chrom_order,
            header: a.header.header,
        },
    };
    let b_formats: Vec<Format> = readers_b.iter().map(|r| r.format()).collect();
    validate(reader_a.format(), &b_formats, &o)?;

    let kind = match mode {
        Mode::Pairs { overlap_bases, .. } => {
            let mut names = Vec::new();
            if readers_b.len() > 1 {
                names.push("b_name".to_string());
            }
            names.push("b_record".to_string());
            if overlap_bases {
                names.push("overlap_bases".to_string());
            }
            ReportOutput::Tsv { appended: names }
        }
        Mode::Count => ReportOutput::Tsv {
            appended: vec!["count".to_string()],
        },
        _ => ReportOutput::SameAsA,
    };
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
    let mut w = open_report(reader_a.as_ref(), &kind, &spec, &oo, &wo, a.header.header)?;
    let inputs = RelateInputs {
        a: reader_a,
        b: readers_b,
    };
    intersect(inputs, &o, w.as_mut())?;
    w.finish()?;
    Ok(())
}
