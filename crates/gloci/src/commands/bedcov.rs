//! `gloci bedcov`: depth of coverage summed over each BED region (D140).
use crate::Ctx;
use crate::cli::BedcovArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_output, open_reader, parse_flags, parse_output_spec,
};
use gloci_ops::bedcov::{BedcovOptions, bedcov};
use gloci_ops::depth::CountRules;

pub fn run(ctx: &Ctx, a: &BedcovArgs) -> anyhow::Result<()> {
    let flags =
        |name: &str, v: &str| parse_flags(v).map_err(|e| UsageError(format!("{name}: {e}")));
    let align = AlignOptions {
        min_mapq: a.min_mapq,
        require_flags: match &a.require_flags {
            Some(s) => flags("--require-flags", s)?,
            None => 0,
        },
        exclude_flags: flags("--exclude-flags", &a.exclude_flags)?,
        full_records: false,
        reference: a.reference.clone(),
        score_tag: None,
    };
    let spec = parse_output_spec(a.output.as_deref(), Some("tsv"))?;
    let chrom_map = a.chrom_map.load()?;
    let invalid = if a.invalid.skip_invalid {
        InvalidPolicy::Skip
    } else {
        InvalidPolicy::Fail
    };
    let base = ReaderOptions {
        accept: Accept::default(),
        chrom_map,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::default(),
        invalid,
    };
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(align),
            ..Accept::default()
        },
        ..base.clone()
    };
    let o = BedcovOptions {
        rules: CountRules {
            deletions: !a.split_deletions,
            ref_skips: !a.split,
            min_baseq: 0,
        },
        bases_at_depth: a.bases_at_depth,
        read_count: a.read_count,
        column_names: a.columns.column_names,
    };
    // Open the BED before the output is created.
    let bed = open_reader(&a.bed, &base)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    bedcov(bed, &a.files, &o, &ro, ctx.threads, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
