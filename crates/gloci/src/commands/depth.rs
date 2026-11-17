//! `gloci depth`: depth of coverage of alignments at each position (D137, D138).
use crate::Ctx;
use crate::cli::DepthArgs;
use gloci_core::{Interval, Region};
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_output, parse_flags,
};
use gloci_ops::depth::{
    AllPositions, By, CountRules, CoverageOptions, DepthOptions, MosdepthOptions, coverage_summary,
    depth_regions, depth_sources, mosdepth_files, open_depth_inputs,
};

pub fn run(ctx: &Ctx, a: &DepthArgs) -> anyhow::Result<()> {
    let flags =
        |name: &str, v: &str| parse_flags(v).map_err(|e| UsageError(format!("{name}: {e}")));
    if let Some(prefix) = &a.out_prefix {
        return run_files(ctx, a, prefix, flags);
    }
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
    let spec = super::own_tsv_output_spec(a.output.as_deref(), "depth")?;
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
    if a.summary {
        return run_summary(ctx, a, &spec, &ro, inputs);
    }
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

/// `--summary`: samtools coverage's table, for the one `--region` if there is one.
fn run_summary(
    ctx: &Ctx,
    a: &DepthArgs,
    spec: &gloci_io::OutputSpec,
    ro: &ReaderOptions,
    inputs: Vec<Box<dyn gloci_io::AlnSource>>,
) -> anyhow::Result<()> {
    if a.region.region.len() > 1 {
        return Err(
            UsageError("--summary takes one --region, as samtools coverage -r".into()).into(),
        );
    }
    let region = match a.region.region.first() {
        None => None,
        Some(s) => {
            let hdr = inputs[0].header();
            let known: std::collections::HashSet<&[u8]> =
                hdr.refs.iter().map(|(n, _)| &n[..]).collect();
            let r =
                Region::parse(s, |n| known.contains(n)).map_err(|e| UsageError(e.to_string()))?;
            let Some((_, ln)) = hdr.refs.iter().find(|(n, _)| *n == r.chrom) else {
                return Err(UsageError(format!(
                    "--region {s}: {} has no reference {}",
                    inputs[0].name(),
                    String::from_utf8_lossy(&r.chrom)
                ))
                .into());
            };
            let window = Interval {
                start: r.start,
                end: r.end.unwrap_or(*ln),
            };
            Some((r.chrom, window))
        }
    };
    let o = CoverageOptions {
        min_baseq: a.min_baseq,
        min_read_length: a.min_read_length,
        region,
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(spec);
    let mut out = open_output(spec, OutFormat::Tsv, &oo)?;
    coverage_summary(inputs, &a.files, &o, ro, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}

/// `--out-prefix`: mosdepth's files. The reader filters nothing: the depth core applies
/// `--exclude-flags` and `--min-mapq` after it notes each record's reference (D173).
fn run_files(
    ctx: &Ctx,
    a: &DepthArgs,
    prefix: &std::path::Path,
    flags: impl Fn(&str, &str) -> Result<u16, UsageError>,
) -> anyhow::Result<()> {
    if a.files.len() != 1 {
        return Err(
            UsageError("--out-prefix writes the files of one FILE, as mosdepth".into()).into(),
        );
    }
    let by = match a.by.as_deref() {
        None => By::Nothing,
        Some("") => {
            return Err(UsageError("--by needs a BED file or a window size".into()).into());
        }
        Some(s) if !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()) => {
            By::Windows(s.parse().map_err(|_| {
                UsageError(format!("--by {s}: a window size that large is not a size"))
            })?)
        }
        Some(s) => By::Bed(s.to_string()),
    };
    if a.thresholds.contains(&0) {
        return Err(UsageError(
            "--thresholds 0: every base has depth at least 0; use 1 or more".into(),
        )
        .into());
    }
    let mut thresholds = a.thresholds.clone();
    thresholds.sort_unstable();
    let o = MosdepthOptions {
        prefix: prefix.to_path_buf(),
        by,
        thresholds,
        exclude_flags: flags("--exclude-flags", &a.exclude_flags)?,
        min_mapq: a.min_mapq,
    };
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                min_mapq: 0,
                require_flags: 0,
                exclude_flags: 0,
                full_records: false,
                reference: a.reference.clone(),
                score_tag: None,
            }),
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
    let mut inputs = open_depth_inputs(&a.files, &ro, ctx.threads)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = format!("{}.mosdepth.*", prefix.display());
    mosdepth_files(inputs.remove(0), &a.files[0], &o, &ro, &oo)
        .map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
