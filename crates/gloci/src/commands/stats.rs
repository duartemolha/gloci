//! `gloci stats`: statistics of interval sets.
use crate::Ctx;
use crate::cli::{FisherArgs, JaccardArgs, ReldistArgs};
use crate::commands::{dest_name, relation_hint, sort_aware, sort_settings, user_span};
use gloci_io::{
    Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError, WriterOptions,
    open_output, parse_output_spec, resolve_genome,
};
use gloci_ops::stats::{
    FisherOptions, JaccardOptions, ReldistOptions, fisher_counts, fisher_table, jaccard, reldist,
    write_fisher, write_jaccard, write_reldist,
};
use gloci_ops::{RelateInputs, ReportOutput, open_report};
use std::io::Write;

fn parse_format(s: Option<&str>) -> anyhow::Result<Option<Format>> {
    match s {
        None => Ok(None),
        Some(s) => Ok(Some(s.parse::<Format>().map_err(UsageError)?)),
    }
}

pub fn run_jaccard(ctx: &Ctx, a: &JaccardArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let chrom_map = a.chrom_map.load()?;
    // Alignments are only read for their spans, never written, so they are slim (D115).
    let accept = gloci_io::Accept {
        alignments: Some(a.align.options(false)?),
        ..gloci_io::Accept::default()
    };
    let reader_opts = |format| {
        sort_aware(
            &a.sort,
            ReaderOptions {
                accept: accept.clone(),
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
    let open = |path: &str, ro: &ReaderOptions| -> anyhow::Result<_> {
        let r = gloci_io::open_reader(path, ro)?;
        gloci_io::check_map_against_index(path, r.as_ref(), ro)?;
        Ok(r)
    };
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?);
    let ro_b = reader_opts(parse_format(a.b_format.as_deref())?);
    let inputs = RelateInputs {
        a: open(&a.a.a, &ro_a)?,
        b: vec![open(&a.b, &ro_b)?],
    };
    let o = JaccardOptions {
        criteria: a.criteria(),
        sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
        genome: a.genome.genome.clone(),
        chrom_order: a.chrom_order.chrom_order,
    };
    let dest = dest_name(&spec);
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let mut out = open_output(&spec, OutFormat::Tsv, &oo).map_err(|e| e.with_output_name(&dest))?;
    let result = jaccard(inputs, &o).map_err(|e| {
        e.with_order_hint(if a.sort.sort {
            ""
        } else {
            ", or pass --sort to sort A and B first (output is the same as for sorted inputs)"
        })
    })?;
    let mut text = Vec::new();
    write_jaccard(&result, &mut text);
    out.write_all(&text)
        .map_err(|source| gloci_io::IoError::Write {
            path: dest.clone(),
            source,
        })?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

pub fn run_fisher(ctx: &Ctx, a: &FisherArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let chrom_map = a.chrom_map.load()?;
    let genome = resolve_genome(&a.genome)?.0;
    let accept = gloci_io::Accept {
        alignments: Some(a.align.options(false)?),
        ..gloci_io::Accept::default()
    };
    let reader_opts = |format| {
        sort_aware(
            &a.sort,
            ReaderOptions {
                accept: accept.clone(),
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
    let open = |path: &str, ro: &ReaderOptions| -> anyhow::Result<_> {
        let r = gloci_io::open_reader(path, ro)?;
        gloci_io::check_map_against_index(path, r.as_ref(), ro)?;
        Ok(r)
    };
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?);
    let ro_b = reader_opts(parse_format(a.b_format.as_deref())?);
    let inputs = RelateInputs {
        a: open(&a.a.a, &ro_a)?,
        b: vec![open(&a.b, &ro_b)?],
    };
    let o = FisherOptions {
        criteria: a.criteria(),
        merge: a.merge_overlaps,
        exclude: a.exclude.clone(),
        sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
        genome: Some(a.genome.clone()),
    };
    let dest = dest_name(&spec);
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let mut out = open_output(&spec, OutFormat::Tsv, &oo).map_err(|e| e.with_output_name(&dest))?;
    let counts = fisher_counts(inputs, &genome, &o).map_err(|e| {
        e.with_order_hint(if a.sort.sort {
            ""
        } else {
            ", or pass --sort to sort A, B and the exclude file first (output is the same as for sorted inputs)"
        })
    })?;
    let table = fisher_table(&counts);
    let mut text = Vec::new();
    write_fisher(&counts, &table, &mut text);
    out.write_all(&text)
        .map_err(|source| gloci_io::IoError::Write {
            path: dest.clone(),
            source,
        })?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

pub fn run_reldist(ctx: &Ctx, a: &ReldistArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let chrom_map = a.chrom_map.load()?;
    // Alignments are only read for their spans (or printed as BED6 rows), D115.
    let accept = gloci_io::Accept {
        alignments: Some(a.align.options(a.detail)?),
        ..gloci_io::Accept::default()
    };
    let reader_opts = |format| {
        sort_aware(
            &a.sort,
            ReaderOptions {
                accept: accept.clone(),
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
    let open = |path: &str, ro: &ReaderOptions| -> anyhow::Result<_> {
        let r = gloci_io::open_reader(path, ro)?;
        gloci_io::check_map_against_index(path, r.as_ref(), ro)?;
        Ok(r)
    };
    let ro_a = reader_opts(parse_format(a.a.a_format.as_deref())?);
    let ro_b = reader_opts(parse_format(a.b_format.as_deref())?);
    let reader_a = open(&a.a.a, &ro_a)?;
    let reader_b = open(&a.b, &ro_b)?;
    let o = ReldistOptions {
        detail: a.detail,
        unsorted: a.u.unsorted,
        sort: sort_settings(ctx, &a.sort, a.vcf_span.vcf_span),
        genome: a.genome.genome.clone(),
        chrom_order: a.chrom_order.chrom_order,
    };
    let dest = dest_name(&spec);
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let hint = relation_hint(&a.sort);
    let order = |e: gloci_io::IoError| e.with_order_hint(hint);
    if a.detail {
        let wo = WriterOptions {
            write_header: false,
            command_line: Some(ctx.command_line.clone()),
            version: ctx.version,
        };
        let mut w = open_report(
            reader_a.as_ref(),
            &ReportOutput::Tsv {
                appended: vec!["reldist".to_string()],
            },
            &spec,
            &oo,
            &wo,
            false,
        )?;
        let inputs = RelateInputs {
            a: reader_a,
            b: vec![reader_b],
        };
        reldist(inputs, &o, Some(w.as_mut())).map_err(|e| order(e.with_output_name(&dest)))?;
        w.finish().map_err(|e| e.with_output_name(&dest))?;
        return Ok(());
    }
    let mut out = open_output(&spec, OutFormat::Tsv, &oo).map_err(|e| e.with_output_name(&dest))?;
    let inputs = RelateInputs {
        a: reader_a,
        b: vec![reader_b],
    };
    let summary = reldist(inputs, &o, None).map_err(order)?;
    let mut text = Vec::new();
    write_reldist(&summary, &mut text);
    out.write_all(&text)
        .map_err(|source| gloci_io::IoError::Write {
            path: dest.clone(),
            source,
        })?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
