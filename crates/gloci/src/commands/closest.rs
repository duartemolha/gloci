//! `gloci closest`: the nearest B to each A.
use crate::Ctx;
use crate::cli::{ClosestArgs, ClosestTies, SignedFrame, VcfSpan};
use crate::commands::intersect::b_column_names;
use gloci_core::{OverlapCriteria, SignFrame, StrandRule};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::closest::{ClosestOptions, Ties, closest, validate};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, open_report};

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

    // Column names matter only for the column-name line of a VCF/BCF A
    // (D42), so only then is B's width needed.
    let mut names = Vec::new();
    let multi = readers_b.len() > 1;
    let shows_distance = o.distance || o.signed.is_some();
    let sides: &[&str] = if o.each_side {
        &["up_", "down_"]
    } else {
        &[""]
    };
    if matches!(reader_a.format(), Format::Vcf | Format::Bcf) {
        let b = b_column_names(readers_b[0].as_ref(), &a.b.b[0], &ro_b);
        for side in sides {
            if multi {
                names.push(format!("{side}b_name"));
            }
            names.extend(b.iter().map(|n| format!("{side}{n}")));
            if shows_distance {
                names.push(format!("{side}distance"));
            }
        }
    }
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
    let inputs = RelateInputs {
        a: reader_a,
        b: readers_b,
    };
    closest(inputs, &o, w.as_mut())?;
    w.finish()?;
    Ok(())
}
