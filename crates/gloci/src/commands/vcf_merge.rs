//! `gloci vcf merge` (P3b Task 23).

use crate::Ctx;
use crate::cli::{FilterLogicArg, MergeModeArg, VcfMergeArgs};
use gloci_io::{
    Format, Header, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions, UsageError,
    VcfHeader, WriterOptions, chrom_table_lenient, open_output, open_reader, open_record_writer,
    parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{FilterList, vcf_header, vcf_merge, vcf_merge_output};
use gloci_vcf::merge::{MergeMode, MergeOptions, Merger, header::HeaderIn, header::merge_headers};
use std::io::IsTerminal;

/// The bcftools merge options gloci does not offer, each a usage error.
fn refuse_unsupported(a: &VcfMergeArgs) -> anyhow::Result<()> {
    let no = |opt: &str, why: &str| -> anyhow::Result<()> {
        Err(UsageError(format!(
            "bcftools merge {opt} is not available in gloci vcf merge: {why}"
        ))
        .into())
    };
    if a.gvcf.is_some() {
        return no(
            "-g/--gvcf",
            "gVCF blocks are not merged; merge the called VCFs, or use bcftools merge -g",
        );
    }
    if a.local_alleles.is_some() {
        return no(
            "-L/--local-alleles",
            "the merged records keep their full allele lists",
        );
    }
    if a.print_header {
        return no(
            "--print-header",
            "write the merged file and read its header with gloci view --header-only",
        );
    }
    if a.use_header.is_some() {
        return no(
            "--use-header",
            "merge, then replace the header with bcftools reheader",
        );
    }
    if a.regions.is_some() || a.regions_file.is_some() || a.targets_file.is_some() {
        return no(
            "-r/--regions, -R/--regions-file and -T/--targets-file",
            "gloci reads the inputs as streams without an index; select regions after merging, e.g. gloci vcf merge ... | gloci view - REGION",
        );
    }
    if a.missing_rules.is_some() {
        return no(
            "-M/--missing-rules",
            "missing values of a re-indexed vector are always .",
        );
    }
    if a.regions_overlap.is_some() {
        return no(
            "--regions-overlap",
            "regions are not available either; select them after merging",
        );
    }
    if a.merge_logic.is_some() {
        return no(
            "--merge-logic",
            "it belongs to the gVCF mode (-g), which is not available",
        );
    }
    Ok(())
}

fn input_files(a: &VcfMergeArgs) -> anyhow::Result<Vec<String>> {
    let mut files = a.files.clone();
    if let Some(list) = &a.file_list {
        if !files.is_empty() {
            return Err(UsageError(
                "give the input files either as arguments or with --file-list, not both"
                    .to_string(),
            )
            .into());
        }
        let text = std::fs::read(list).map_err(|source| IoError::Open {
            path: list.display().to_string(),
            source,
        })?;
        files = text
            .split(|&b| b == b'\n')
            .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
            .filter(|l| !l.is_empty())
            .map(|l| String::from_utf8_lossy(l).into_owned())
            .collect();
    }
    if files.len() < 2 && !(a.force_single && files.len() == 1) {
        return Err(UsageError(format!(
            "vcf merge needs at least two input files, got {}; pass them as arguments or with --file-list FILE",
            files.len()
        ))
        .into());
    }
    if files.iter().any(|f| f == "-") {
        return Err(UsageError(
            "vcf merge cannot read - (stdin): the inputs are read together, so each needs a path (process substitution works)"
                .to_string(),
        )
        .into());
    }
    Ok(files)
}

pub fn run(ctx: &Ctx, a: &VcfMergeArgs) -> anyhow::Result<()> {
    refuse_unsupported(a)?;
    let files = input_files(a)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    // Many inputs are open at once: each reads on its own thread only when there are few.
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: if files.len() <= 8 {
            ctx.threads
        } else {
            std::num::NonZeroUsize::MIN
        },
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
    };
    let mut readers = Vec::with_capacity(files.len());
    for f in &files {
        readers.push(open_reader(f, &ro).map_err(|e| match e {
            IoError::Open { path, source } if source.raw_os_error() == Some(24) => {
                IoError::Detect {
                    path: path.clone(),
                    msg: format!(
                        "cannot open {path}: too many open files (merging {} inputs needs that many file handles at once); raise the limit with ulimit -n {} and run again",
                        files.len(),
                        (files.len() + 64).next_power_of_two()
                    ),
                }
            }
            other => other,
        })?);
    }
    let headers: Vec<VcfHeader> = readers
        .iter()
        .map(|r| vcf_header("vcf merge", r.as_ref()))
        .collect::<Result<_, _>>()?;
    let inputs: Vec<HeaderIn> = headers
        .iter()
        .zip(&files)
        .map(|(h, n)| HeaderIn {
            name: n,
            raw: &h.raw,
        })
        .collect();
    let mut merged = merge_headers(&inputs, a.force_samples).map_err(|msg| IoError::Detect {
        path: "vcf merge".to_string(),
        msg,
    })?;
    let mode = match a.merge {
        MergeModeArg::None => MergeMode::None,
        MergeModeArg::Snps => MergeMode::Snps,
        MergeModeArg::Indels => MergeMode::Indels,
        MergeModeArg::Both => MergeMode::Both,
        MergeModeArg::All => MergeMode::All,
        MergeModeArg::Id => MergeMode::Id,
    };
    let options = MergeOptions {
        mode,
        missing_to_ref: a.missing_to_ref,
        info_rules: a.info_rules.clone(),
        filter_x: a.filter_logic == FilterLogicArg::Pass,
    };
    let mut merger = Merger::new(options, &mut merged, &files).map_err(UsageError)?;
    let (to, header) = vcf_merge_output(readers[0].as_ref(), spec.format, &merged.text)?;
    let Header::Vcf(out_header) = &header else {
        unreachable!("vcf_merge_output returns a VCF header")
    };
    let out_header: VcfHeader = out_header.clone();
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let refs: Vec<&dyn gloci_io::RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let mut chroms = chrom_table_lenient(None, &refs, a.chrom_order.chrom_order)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    // As bcftools merge --no-version: no command line is added to the header.
    let wo = WriterOptions {
        write_header: true,
        command_line: None,
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
    let filter = a.apply_filters.as_deref().map(FilterList::parse);
    let sum = vcf_merge(
        readers,
        &mut chroms,
        &mut merger,
        filter,
        &out_header,
        writer.as_mut(),
    )
    .map_err(|e| e.with_output_name(&dest))?;
    writer.finish().map_err(|e| e.with_output_name(&dest))?;
    log::info!(
        "vcf merge: {} records read from {} files, {} written ({} merged from several files, {} skipped by --apply-filters); at most {} records held at once",
        sum.records_in,
        files.len(),
        sum.records_out,
        sum.records_merged,
        sum.records_filtered,
        sum.peak_buffered
    );
    Ok(())
}
