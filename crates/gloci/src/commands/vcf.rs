//! `gloci vcf`: selection pipeline and other VCF commands (spec §7.8).

use crate::Ctx;
use crate::cli::{
    ChromMapArg, ExtractLayout, InvalidArg, VcfExtractArgs, VcfFilterArgs, VcfSelectArgs, VcfSpan,
    VcfSpanArg,
};
use gloci_core::{Interval, RegionSet};
use gloci_io::{
    Format, InvalidPolicy, OutFormat, OutputOptions, OutputSpec, ReaderOptions, UsageError,
    WriterOptions, open_output, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{
    ExtractMode, ExtractOptions, FilterExpr, VcfFilterOptions, VcfSelect, check_keys,
    extract_matrix, extract_table, open_selected, open_vcf_filter, read_regions_file,
    validate_thresholds, vcf_filter_file, vcf_filter_output, vcf_header,
};
use gloci_vcf::{
    MaskRules, RewritePlan, SiteFilters, VcfType, read_ids, read_positions, read_samples,
};
use std::collections::{HashMap, HashSet};
use std::io::IsTerminal;

/// Build the selection of a `vcf` command from its flags, reading the list
/// files the flags name. `make_ro` builds the reader options and is called
/// after the lists are read (so a list-file error is reported first); it is
/// also what `--exclude-regions-file` is read with. The pipeline stage that
/// applies the result is `gloci_ops::SelectStage`.
pub fn select_from_args(
    a: &VcfSelectArgs,
    make_ro: impl FnOnce() -> anyhow::Result<ReaderOptions>,
) -> anyhow::Result<(VcfSelect, ReaderOptions)> {
    let (min_alleles, max_alleles) = if a.biallelic {
        (Some(2), Some(2))
    } else {
        (a.min_alleles, a.max_alleles)
    };
    // Usage errors first: before any list file is read.
    validate_thresholds(
        a.min_maf,
        a.max_maf,
        a.max_missing,
        min_alleles,
        max_alleles,
    )?;
    let include_samples = match (&a.samples, &a.samples_file) {
        (Some(list), None) => Some(list.clone()),
        (None, Some(path)) => Some(read_samples(&path.to_string_lossy())?),
        _ => None,
    };
    let exclude_samples = match (&a.exclude_samples, &a.exclude_samples_file) {
        (Some(list), None) => list.clone(),
        (None, Some(path)) => read_samples(&path.to_string_lossy())?,
        _ => Vec::new(),
    };

    let positions = match &a.positions {
        Some(p) => Some(read_positions(&p.to_string_lossy())?),
        None => None,
    };
    let exclude_positions = match &a.exclude_positions {
        Some(p) => Some(read_positions(&p.to_string_lossy())?),
        None => None,
    };
    let ids = match &a.ids {
        Some(p) => Some(read_ids(&p.to_string_lossy())?),
        None => None,
    };
    let exclude_ids = match &a.exclude_ids {
        Some(p) => Some(read_ids(&p.to_string_lossy())?),
        None => None,
    };

    let types = if a.types.is_empty() {
        None
    } else {
        Some(a.types.iter().map(|&t| VcfType::from(t)).collect())
    };
    let exclude_types = if a.exclude_types.is_empty() {
        None
    } else {
        Some(a.exclude_types.iter().map(|&t| VcfType::from(t)).collect())
    };

    let keep_filter = if a.keep_filter.is_empty() {
        None
    } else {
        Some(
            a.keep_filter
                .iter()
                .map(|s| s.as_bytes().to_vec())
                .collect(),
        )
    };
    let remove_filter = if a.remove_filter.is_empty() {
        None
    } else {
        Some(
            a.remove_filter
                .iter()
                .map(|s| s.as_bytes().to_vec())
                .collect(),
        )
    };

    let require_info: Vec<Vec<u8>> = a
        .require_info
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let exclude_info: Vec<Vec<u8>> = a
        .exclude_info
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let chroms = if a.chrom.is_empty() {
        None
    } else {
        Some(a.chrom.iter().map(|s| s.as_bytes().to_vec()).collect())
    };
    let exclude_chroms: HashSet<Vec<u8>> = a
        .exclude_chrom
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();

    let ro = make_ro()?;

    let exclude_regions = match &a.exclude_regions_file {
        Some(p) => {
            let ivs = read_regions_file(&p.to_string_lossy(), &ro)?;
            let mut by_chrom: HashMap<Vec<u8>, Vec<Interval>> = HashMap::new();
            for (c, iv) in ivs {
                by_chrom.entry(c).or_default().push(iv);
            }
            Some(
                by_chrom
                    .into_iter()
                    .map(|(c, list)| (c, RegionSet::new(list)))
                    .collect(),
            )
        }
        None => None,
    };

    let mask_rules = MaskRules {
        min_gq: a.min_gq,
        min_dp: a.min_dp,
        max_dp: a.max_dp,
        mask_filtered: a.mask_filtered_genotypes,
    };
    let select = VcfSelect {
        include_samples,
        exclude_samples,
        mask: mask_rules,
        sites: SiteFilters {
            chroms,
            exclude_chroms,
            positions,
            exclude_positions,
            ids,
            exclude_ids,
            types,
            exclude_types,
            min_alleles,
            max_alleles,
            min_qual: a.min_qual,
            pass_only: a.pass_only,
            keep_filter,
            remove_filter,
            require_info,
            exclude_info,
            min_maf: a.min_maf,
            max_maf: a.max_maf,
            min_mac: a.min_mac,
            max_mac: a.max_mac,
            max_missing: a.max_missing,
            max_missing_count: a.max_missing_count,
            min_mean_dp: a.min_mean_dp,
            max_mean_dp: a.max_mean_dp,
            hwe: a.hwe,
            thin: a.thin,
            phased_only: a.phased_only,
            exclude_regions,
            inclusive_points: a.inclusive_points,
        },
        expr: match (&a.include, &a.exclude) {
            (Some(src), _) => Some(FilterExpr {
                src: src.clone(),
                exclude: false,
            }),
            (_, Some(src)) => Some(FilterExpr {
                src: src.clone(),
                exclude: true,
            }),
            _ => None,
        },
        regions: a.region.region.clone(),
        regions_file: a
            .region
            .regions_file
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
    };
    select.validate()?;
    Ok((select, ro))
}

fn reader_options(
    ctx: &Ctx,
    chrom_map: &ChromMapArg,
    format: Option<&str>,
    span: &VcfSpanArg,
    invalid: &InvalidArg,
) -> anyhow::Result<ReaderOptions> {
    Ok(ReaderOptions {
        chrom_map: chrom_map.load()?,
        threads: ctx.threads,
        format: match format {
            None => None,
            Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: match span.vcf_span {
            VcfSpan::Allele => gloci_core::VcfSpan::Allele,
            VcfSpan::Ref => gloci_core::VcfSpan::Ref,
        },
        invalid: if invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    })
}

pub fn run_filter(ctx: &Ctx, a: &VcfFilterArgs) -> anyhow::Result<()> {
    let (select, ro) = select_from_args(&a.select, || {
        reader_options(
            ctx,
            &a.chrom_map,
            a.format.format.as_deref(),
            &a.vcf_span,
            &a.invalid,
        )
    })?;
    let options = VcfFilterOptions {
        select,
        plan: RewritePlan {
            selection: None,
            mask: MaskRules::default(),
            keep_info: a
                .keep_info
                .as_ref()
                .map(|v| v.iter().map(|s| s.as_bytes().to_vec()).collect()),
            drop_info: a.drop_info.iter().map(|s| s.as_bytes().to_vec()).collect(),
            recalc_ac: a.recalc_ac,
        },
    };

    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let input = open_vcf_filter(&a.file, &ro, &options)?;
    let (to, header) = vcf_filter_output(input.reader(), spec.format, &options)?;
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
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
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "stdout".to_string(), |p| p.display().to_string());
    let out = open_output(&spec, to, &oo)?;
    let record_format = if to == OutFormat::Bcf {
        Format::Bcf
    } else {
        Format::Vcf
    };
    let mut writer = open_record_writer(out, record_format, to, &header, &wo)
        .map_err(|e| e.with_output_name(&dest))?;
    vcf_filter_file(input, &options, writer.as_mut()).map_err(|e| e.with_output_name(&dest))?;
    writer.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

pub fn run_extract(ctx: &Ctx, a: &VcfExtractArgs) -> anyhow::Result<()> {
    if a.layout.is_some() && a.format.is_none() {
        return Err(UsageError(
            "--layout applies to --format only; add --format KEYS or drop --layout".into(),
        )
        .into());
    }
    if a.out_prefix.is_some() && a.genotype_matrix.is_none() {
        return Err(UsageError(
            "--out-prefix names the --genotype-matrix files; add --genotype-matrix 012 \
or use -o FILE for a table"
                .into(),
        )
        .into());
    }
    let mode = match (&a.info, &a.format) {
        (Some(keys), _) => {
            check_keys("--info", keys)?;
            ExtractMode::Info(keys.clone())
        }
        (_, Some(keys)) => {
            check_keys("--format", keys)?;
            ExtractMode::Format {
                keys: keys.clone(),
                long: a.layout == Some(ExtractLayout::Long),
            }
        }
        _ => ExtractMode::Matrix012,
    };
    let (select, ro) = select_from_args(&a.select, || {
        reader_options(ctx, &a.chrom_map, None, &a.vcf_span, &a.invalid)
    })?;
    let options = ExtractOptions {
        select,
        mode,
        skip_invalid: a.invalid.skip_invalid,
        max_mem: a.max_mem,
        tmp_dir: ctx.tmp_dir.clone(),
    };
    let input = open_selected(&a.file, &ro, &options.select)?;
    vcf_header("vcf extract", input.reader())?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.compression_level,
        write_index: false,
    };

    if options.mode == ExtractMode::Matrix012 {
        let prefix = a
            .out_prefix
            .as_ref()
            .ok_or_else(|| UsageError("--genotype-matrix needs --out-prefix PREFIX".into()))?;
        let open = |suffix: &str| {
            let mut path = prefix.clone().into_os_string();
            path.push(suffix);
            let spec = OutputSpec {
                path: Some(path.into()),
                format: None,
                bgzf: false,
                direct: false,
            };
            open_output(&spec, OutFormat::Tsv, &oo)
        };
        let mut matrix = open(".012")?;
        let mut indv = open(".012.indv")?;
        let mut pos = open(".012.pos")?;
        let dest = prefix.display().to_string();
        extract_matrix(input, &options, &mut matrix, &mut indv, &mut pos)
            .map_err(|e| e.with_output_name(&dest))?;
        indv.finish()?;
        pos.finish()?;
        matrix.finish()?;
        return Ok(());
    }

    let spec = parse_output_spec(a.output.as_deref(), None)?;
    if let Some(f) = spec.format
        && f != OutFormat::Tsv
    {
        return Err(UsageError(format!(
            "vcf extract writes a TSV, not {}; name the output .tsv or .tsv.gz",
            f.name()
        ))
        .into());
    }
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "stdout".to_string(), |p| p.display().to_string());
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    extract_table(input, &options, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
