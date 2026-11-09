//! `gloci vcf`: selection pipeline and other VCF commands (spec §7.8).

use crate::Ctx;
use crate::cli::{VcfFilterArgs, VcfSpan};
use gloci_core::{Interval, RegionSet};
use gloci_io::{
    Format, Header, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions, UsageError,
    WriterOptions, open_output, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{VcfFilterOptions, ViewOptions, open_view, read_regions_file, vcf_filter};
use gloci_vcf::{
    MaskRules, RewritePlan, SiteFilters, VcfType, read_ids, read_positions, read_samples,
    rewrite_header, select_samples,
};
use std::collections::{HashMap, HashSet};
use std::io::IsTerminal;

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

pub fn run_filter(ctx: &Ctx, a: &VcfFilterArgs) -> anyhow::Result<()> {
    if a.biallelic && (a.min_alleles.is_some() || a.max_alleles.is_some()) {
        return Err(
            usage("--biallelic cannot be combined with --min-alleles or --max-alleles").into(),
        );
    }
    if let Some(m) = a.min_maf
        && !(0.0..=1.0).contains(&m)
    {
        return Err(usage(format!("--min-maf {m} is not in [0, 1]")).into());
    }
    if let Some(m) = a.max_maf
        && !(0.0..=1.0).contains(&m)
    {
        return Err(usage(format!("--max-maf {m} is not in [0, 1]")).into());
    }
    if let Some(m) = a.max_missing
        && !(0.0..=1.0).contains(&m)
    {
        return Err(usage(format!("--max-missing {m} is not in [0, 1]")).into());
    }
    if a.samples.is_some() && a.samples_file.is_some() {
        return Err(usage("--samples and --samples-file cannot be combined").into());
    }
    if a.exclude_samples.is_some() && a.exclude_samples_file.is_some() {
        return Err(
            usage("--exclude-samples and --exclude-samples-file cannot be combined").into(),
        );
    }

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

    let (min_alleles, max_alleles) = if a.biallelic {
        (Some(2), Some(2))
    } else {
        (a.min_alleles, a.max_alleles)
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

    let vo = ViewOptions {
        header_only: false,
        no_header: false,
        regions: a.region.region.clone(),
        regions_file: a
            .region
            .regions_file
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        inclusive_points: a.inclusive_points,
        genome: None,
        chrom_order: gloci_core::ChromOrder::Karyotype,
        assume_sorted: false,
    };

    let (input, mut chrom_table) = open_view(&a.file, &vo, &ro)?;
    let in_format = input.format();
    if in_format != Format::Vcf && in_format != Format::Bcf {
        return Err(usage(format!(
            "vcf filter works on VCF and BCF input, but {} is {:?}",
            input.name(),
            in_format
        ))
        .into());
    }

    let Header::Vcf(h) = input.header().clone() else {
        return Err(IoError::Detect {
            path: input.name().to_string(),
            msg: "a VCF/BCF input has no VCF header".into(),
        }
        .into());
    };

    let selection = if include_samples.is_some() || !exclude_samples.is_empty() {
        Some(select_samples(
            &h.samples,
            include_samples.as_deref(),
            &exclude_samples,
        )?)
    } else {
        None
    };

    let mask_rules = MaskRules {
        min_gq: a.min_gq,
        min_dp: a.min_dp,
        max_dp: a.max_dp,
        mask_filtered: a.mask_filtered_genotypes,
    };

    let plan = RewritePlan {
        selection,
        mask: mask_rules,
        keep_info,
        drop_info,
        recalc_ac: a.recalc_ac,
    };

    let rewritten_header_raw = rewrite_header(&h, &plan);
    let out_header =
        gloci_io::VcfHeader::parse(rewritten_header_raw).map_err(|e| IoError::Detect {
            path: input.name().to_string(),
            msg: e,
        })?;

    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let family = OutFormat::family(in_format);
    let to = spec.format.unwrap_or(family);
    if to != OutFormat::Vcf && to != OutFormat::Bcf {
        return Err(usage(format!(
            "output format for vcf filter must be VCF or BCF, not {:?}",
            to
        ))
        .into());
    }

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

    let out = open_output(&spec, to, &oo)?;
    let mut writer = open_record_writer(
        out,
        if to == OutFormat::Bcf {
            Format::Bcf
        } else {
            Format::Vcf
        },
        to,
        &Header::Vcf(out_header),
        &wo,
    )?;

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
        plan,
        // open_view wraps the reader with in-stream or indexed region selection.
        regions: Vec::new(),
        regions_file: None,
    };

    vcf_filter(input, &mut chrom_table, &options, writer.as_mut())?;
    writer.finish()?;

    Ok(())
}
