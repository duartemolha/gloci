//! `gloci vcf`: selection pipeline and other VCF commands (spec §7.8).

use crate::Ctx;
use crate::cli::{VcfFilterArgs, VcfSpan};
use gloci_core::{Interval, RegionSet};
use gloci_io::{
    Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError, WriterOptions,
    open_output, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{
    VcfFilterOptions, open_vcf_filter, read_regions_file, validate_thresholds, vcf_filter_file,
    vcf_filter_output,
};
use gloci_vcf::{
    MaskRules, RewritePlan, SiteFilters, VcfType, read_ids, read_positions, read_samples,
};
use std::collections::{HashMap, HashSet};
use std::io::IsTerminal;

pub fn run_filter(ctx: &Ctx, a: &VcfFilterArgs) -> anyhow::Result<()> {
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
    let keep_info: Option<Vec<Vec<u8>>> = a
        .keep_info
        .as_ref()
        .map(|v| v.iter().map(|s| s.as_bytes().to_vec()).collect());
    let drop_info: Vec<Vec<u8>> = a.drop_info.iter().map(|s| s.as_bytes().to_vec()).collect();

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

    let ro = ReaderOptions {
        threads: ctx.threads,
        format: match a.format.format.as_deref() {
            None => None,
            Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
        },
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
    let options = VcfFilterOptions {
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
        plan: RewritePlan {
            selection: None,
            mask: mask_rules,
            keep_info,
            drop_info,
            recalc_ac: a.recalc_ac,
        },
        regions: a.region.region.clone(),
        regions_file: a
            .region
            .regions_file
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
    };
    options.validate()?;

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
