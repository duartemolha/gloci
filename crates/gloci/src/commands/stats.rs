//! `gloci stats`: statistics of interval sets.
use crate::Ctx;
use crate::cli::JaccardArgs;
use crate::commands::{dest_name, sort_aware, sort_settings, user_span};
use gloci_io::{
    Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError, open_output,
    parse_output_spec,
};
use gloci_ops::RelateInputs;
use gloci_ops::stats::{JaccardOptions, jaccard, write_jaccard};
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
