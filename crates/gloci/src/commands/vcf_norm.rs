//! `gloci vcf norm` (P3b Task 21).

use crate::Ctx;
use crate::cli::{CheckRefArg, RmDupArg, SplitTypesArg, VcfNormArgs};
use gloci_core::ChromTable;
use gloci_io::{
    Format, IndexedFasta, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    WriterOptions, open_output, open_reader, open_record_writer, parse_output_spec,
    refuse_binary_to_tty,
};
use gloci_ops::{NormOptions, vcf_norm, vcf_norm_output};
use gloci_vcf::{CheckRef, RmDup, SplitTypes};
use std::io::IsTerminal;

/// The bcftools norm options gloci does not offer, each a usage error.
fn refuse_unsupported(a: &VcfNormArgs) -> anyhow::Result<()> {
    let no = |opt: &str, why: &str| -> anyhow::Result<()> {
        Err(UsageError(format!(
            "bcftools norm {opt} is not available in gloci vcf norm: {why}"
        ))
        .into())
    };
    if a.atomize {
        return no(
            "-a/--atomize",
            "gloci does not decompose complex variants; use gloci vcf norm for left-alignment and trimming only",
        );
    }
    if let Some(m) = &a.multiallelics {
        return if m.starts_with('+') {
            no(
                &format!("-m{m}"),
                "gloci does not join records into multi-allelic ones",
            )
        } else {
            no(
                &format!("-m{m}"),
                "use --split[=snps|indels|both|any] to split multi-allelic records",
            )
        };
    }
    if a.old_rec_tag.is_some() {
        return no(
            "--old-rec-tag",
            "gloci does not annotate modified records with their original variant",
        );
    }
    if a.do_not_normalize {
        return no(
            "-N/--do-not-normalize",
            "for splitting without normalising use gloci vcf split; for REF checks the records are always normalised",
        );
    }
    Ok(())
}

pub fn run(ctx: &Ctx, a: &VcfNormArgs) -> anyhow::Result<()> {
    refuse_unsupported(a)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let options = NormOptions {
        split: a.split.map(|t| match t {
            SplitTypesArg::Snps => SplitTypes::Snps,
            SplitTypesArg::Indels => SplitTypes::Indels,
            SplitTypesArg::Both => SplitTypes::Both,
            SplitTypesArg::Any => SplitTypes::Any,
        }),
        check_ref: match a.check_ref {
            CheckRefArg::E => CheckRef::Exit,
            CheckRefArg::W => CheckRef::Warn,
            CheckRefArg::X => CheckRef::Exclude,
            CheckRefArg::S => CheckRef::Set,
        },
        rm_dup: a.rm_dup.map(|d| match d {
            RmDupArg::Snps => RmDup::Snps,
            RmDupArg::Indels => RmDup::Indels,
            RmDupArg::Both => RmDup::Both,
            RmDupArg::All => RmDup::All,
            RmDupArg::Exact => RmDup::Exact,
        }),
        window: a.window,
        skip_invalid: a.invalid.skip_invalid,
    };
    // A missing or stale index is reported here, before any output is made.
    let fasta = IndexedFasta::open(&a.fasta)?;
    let input = open_reader(&a.file, &ro)?;
    let (to, header) = vcf_norm_output("vcf norm", input.as_ref(), spec.format)?;
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
    // As bcftools norm --no-version: no command line is added to the header.
    let wo = WriterOptions {
        write_header: true,
        command_line: None,
        version: ctx.version,
    };
    let dest = spec.path.as_ref().map_or_else(
        || "stdout".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
    let out = open_output(&spec, to, &oo)?;
    let record_format = if to == OutFormat::Bcf {
        Format::Bcf
    } else {
        Format::Vcf
    };
    let mut writer = open_record_writer(out, record_format, to, &header, &wo)
        .map_err(|e| e.with_output_name(&dest))?;
    let mut chroms = ChromTable::first_appearance();
    let sum = vcf_norm(input, &mut chroms, fasta, &options, writer.as_mut())
        .map_err(|e| e.with_output_name(&dest))?;
    writer.finish().map_err(|e| e.with_output_name(&dest))?;
    log::info!(
        "vcf norm: {} records read, {} split, {} realigned, {} excluded, {} duplicates removed, {} written, {} skipped; peak reference windows {} bases",
        sum.records_in,
        sum.records_split,
        sum.realigned,
        sum.excluded,
        sum.duplicates,
        sum.records_out,
        sum.records_skipped,
        sum.peak_reference_bases
    );
    Ok(())
}
