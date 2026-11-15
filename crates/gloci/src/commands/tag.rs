//! `gloci tag`: the alignments of a file, each tagged with the features it overlaps (D168).
use crate::Ctx;
use crate::cli::TagArgs;
use gloci_io::{
    Accept, AlignOptions, Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions,
    UsageError, WriterOptions, open_aln_source, open_output, open_reader, parse_output_spec,
    refuse_binary_to_tty,
};
use gloci_ops::tag::{
    DEFAULT_FRACTION, TagOptions, TagStrand, TagValue, load_features, parse_tag_name,
    tag_alignments, validate,
};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &TagArgs) -> anyhow::Result<()> {
    let labels = || a.labels.clone().unwrap_or_default();
    let value = if a.names {
        TagValue::Names
    } else if a.scores {
        TagValue::Scores
    } else if a.intervals {
        TagValue::Intervals(labels())
    } else {
        TagValue::Labels(labels())
    };
    let o = TagOptions {
        tag: parse_tag_name(&a.tag)?,
        value,
        strand: match (a.same_strand, a.opposite_strand) {
            (true, _) => TagStrand::Same,
            (_, true) => TagStrand::Opposite,
            _ => TagStrand::Any,
        },
        fraction: a.fraction.unwrap_or(DEFAULT_FRACTION),
    };
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let base = ReaderOptions {
        accept: Accept::default(),
        chrom_map: None,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::default(),
        invalid: InvalidPolicy::Fail,
    };
    let src = open_aln_source(
        &a.input,
        &ReaderOptions {
            accept: Accept {
                alignments: Some(AlignOptions {
                    full_records: true,
                    reference: a.reference.clone(),
                    ..AlignOptions::default()
                }),
                ..Accept::default()
            },
            ..base.clone()
        },
    )?;
    gloci_io::warn_unused_reference(a.reference.as_deref(), src.format() == Format::Cram);
    // The input's format: SAM for SAM, else BAM (CRAM is never written); -o and -O choose.
    let to = spec.format.unwrap_or(if src.format() == Format::Sam {
        OutFormat::Sam
    } else {
        OutFormat::Bam
    });
    if !matches!(to, OutFormat::Bam | OutFormat::Sam) {
        return Err(UsageError(format!(
            "{} output is not possible for tagged alignments; write BAM or SAM (-o x.bam, \
-o x.sam); CRAM input is written as BAM",
            to.name()
        ))
        .into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    // The features are read before the output is created.
    let readers =
        a.b.iter()
            .map(|p| open_reader(p, &base))
            .collect::<Result<Vec<_>, _>>()?;
    let formats: Vec<(&str, Format)> = readers.iter().map(|r| (r.name(), r.format())).collect();
    validate(&o, &formats)?;
    let features = load_features(readers, &o)?;
    log::info!("holding {} features", features.len());
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let wo = WriterOptions {
        write_header: true,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let dest = super::dest_name(&spec);
    let out = open_output(&spec, to, &oo).map_err(|e| e.with_output_name(&dest))?;
    tag_alignments(src, &features, &o, out, to, &wo).map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
