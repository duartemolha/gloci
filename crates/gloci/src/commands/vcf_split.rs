//! `gloci vcf split` (P3b Task 20).

use crate::Ctx;
use crate::cli::{SplitTypesArg, VcfSplitArgs};
use gloci_core::ChromTable;
use gloci_io::{
    Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, WriterOptions, open_output,
    open_reader, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{SplitOptions, vcf_split, vcf_split_output};
use gloci_vcf::SplitTypes;
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &VcfSplitArgs) -> anyhow::Result<()> {
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
    let options = SplitOptions {
        types: match a.types {
            SplitTypesArg::Snps => SplitTypes::Snps,
            SplitTypesArg::Indels => SplitTypes::Indels,
            SplitTypesArg::Both => SplitTypes::Both,
            SplitTypesArg::Any => SplitTypes::Any,
        },
        skip_invalid: a.invalid.skip_invalid,
    };
    let input = open_reader(&a.file, &ro)?;
    let (to, header) = vcf_split_output(input.as_ref(), spec.format)?;
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
    let sum = vcf_split(input, &mut chroms, &options, writer.as_mut())
        .map_err(|e| e.with_output_name(&dest))?;
    writer.finish().map_err(|e| e.with_output_name(&dest))?;
    log::info!(
        "vcf split: {} records read, {} split, {} written, {} skipped",
        sum.records_in,
        sum.records_split,
        sum.records_out,
        sum.records_skipped
    );
    Ok(())
}
