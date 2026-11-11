//! `gloci depth`: depth of coverage of alignments at each position (D137, D138).
use crate::Ctx;
use crate::cli::DepthArgs;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_output, parse_flags, parse_output_spec,
};
use gloci_ops::depth::{AllPositions, CountRules, DepthOptions, depth_in, depth_regions};

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
    let regions = depth_regions(
        &a.files[0],
        &a.region.region,
        a.region.regions_file.as_ref().and_then(|p| p.to_str()),
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
    depth_in(&a.files, &o, &ro, ctx.threads, &ctx.tmp_dir, &mut out)
        .map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
