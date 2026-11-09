//! `gloci intersect`: A's records by their overlaps with B.
use crate::Ctx;
use crate::cli::{IntersectArgs, IntersectMatch, IntersectMode};
use crate::commands::{relation_hint, sort_aware, sort_settings, user_span};
use gloci_core::ChromTable;
use gloci_io::{
    Format, Header, InvalidPolicy, OutputOptions, ReaderOptions, RecordBatch, RecordReader,
    UsageError, WriterOptions, column_count, open_reader, parse_output_spec,
};
use gloci_ops::intersect::{IntersectOptions, Match, Mode, intersect, validate};
use gloci_ops::{RelateInputs, RelateOptions, ReportOutput, all_indexed, open_report, run_indexed};

fn usage(msg: impl Into<String>) -> anyhow::Error {
    UsageError(msg.into()).into()
}

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

/// One name per column of B's rows, for the column-name line of a VCF/BCF
/// A's TSV: `b_` plus each field of B's `#CHROM` line for a VCF/BCF B, else
/// `b_1`..`b_n` with n the column count of B's first record (3 for an empty
/// B). With several `-b` files the widths come from the first one, and so
/// does the width of the placeholder row of `--keep-unmatched`. The first
/// record is peeked through a second reader; stdin cannot be reopened, so a
/// BED `-b -` gets 3 names.
pub(crate) fn b_column_names(b: &dyn RecordReader, path: &str, ro: &ReaderOptions) -> Vec<String> {
    let numbered = |n: usize| (1..=n).map(|i| format!("b_{i}")).collect::<Vec<_>>();
    match (b.format(), b.header()) {
        (Format::Vcf | Format::Bcf, Header::Vcf(h)) => {
            let mut names: Vec<String> =
                ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                    .iter()
                    .map(|n| format!("b_{n}"))
                    .collect();
            if !h.samples.is_empty() {
                names.push("b_FORMAT".to_string());
                names.extend(
                    h.samples
                        .iter()
                        .map(|s| format!("b_{}", String::from_utf8_lossy(s))),
                );
            }
            names
        }
        (Format::Gff3 | Format::Gtf, _) => numbered(9),
        _ => {
            let first = (path != "-")
                .then(|| open_reader(path, ro).ok())
                .flatten()
                .and_then(|mut r| {
                    let mut batch = RecordBatch::default();
                    r.read_batch(&mut ChromTable::karyotype(), &mut batch, 1)
                        .ok()?;
                    (!batch.is_empty()).then(|| column_count(batch.raw(0)))
                });
            numbered(first.unwrap_or(3))
        }
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
            sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
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
            // Only a VCF/BCF A gets a column-name line (D42), so only then
            // is B's width needed.
            if matches!(reader_a.format(), Format::Vcf | Format::Bcf) {
                names.extend(b_column_names(readers_b[0].as_ref(), &a.b.b[0], &ro_b));
            }
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
    let all_paths: Vec<&str> = std::iter::once(a.a.a.as_str())
        .chain(a.b.b.iter().map(|s| s.as_str()))
        .collect();
    let use_indexed = all_indexed(&all_paths) && !a.u.unsorted && o.relate.sort.is_none();

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
            &|chunk_inputs, chunk_sink| intersect(chunk_inputs, &o, chunk_sink),
        )
        .map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    } else {
        let inputs = RelateInputs {
            a: reader_a,
            b: readers_b,
        };
        intersect(inputs, &o, w.as_mut()).map_err(|e| e.with_order_hint(relation_hint(&a.sort)))?;
    }
    w.finish()?;
    Ok(())
}
