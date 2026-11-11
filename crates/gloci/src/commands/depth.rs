//! `gloci depth`: depth of coverage of alignments at each position (D137, D138).
use crate::Ctx;
use crate::cli::DepthArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_output, parse_flags, parse_output_spec,
};
use gloci_ops::depth::{
    AllPositions, CountRules, DepthOptions, depth_regions, depth_sources, open_depth_inputs,
};

pub fn run(ctx: &Ctx, a: &DepthArgs) -> anyhow::Result<()> {
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
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(align),
            ..Accept::default()
        },
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::default(),
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    // Each input is opened once: the regions are checked against the first one's header.
    let inputs = open_depth_inputs(&a.files, &ro, ctx.threads)?;
    let regions = depth_regions(
        inputs[0].header(),
        inputs[0].name(),
        &a.region.region,
        // Read as view reads it: a path that is not UTF-8 fails to open, by name, and is never
        // dropped.
        a.region
            .regions_file
            .as_ref()
            .map(|p| p.to_string_lossy())
            .as_deref(),
        &ro,
    )?;
    let o = DepthOptions {
        rules: CountRules {
            deletions: a.count_deletions,
            ref_skips: false,
            min_baseq: a.min_baseq,
        },
        min_read_length: a.min_read_length,
        all: if a.all_references {
            AllPositions::AllRefs
        } else if a.all_positions {
            AllPositions::RefsWithReads
        } else {
            AllPositions::Covered
        },
        regions,
        skip_mate_overlap: a.skip_mate_overlap,
        column_names: a.columns.column_names,
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    depth_sources(inputs, &a.files, &o, &ro, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
