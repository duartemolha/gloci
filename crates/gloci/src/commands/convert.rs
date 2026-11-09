//! `gloci convert`: VCF/BCF, GFF3/GTF and BED12 to BED.
use super::user_span;
use crate::Ctx;
use crate::cli::ConvertArgs;
use gloci_core::{ChromTable, VariantKind};
use gloci_io::{
    Format, Header, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions, TsvHeader,
    UsageError, WriterOptions, open_output, open_reader, open_record_writer, open_tsv_writer,
    parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::convert::{ConvertOptions, check_convert, convert};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &ConvertArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let o = ConvertOptions {
        vcf_span: user_span(a.vcf_span.vcf_span),
        split_alleles: a.split_alleles,
        types: (!a.types.is_empty())
            .then(|| a.types.iter().map(|&k| VariantKind::from(k)).collect()),
        features: a.feature.clone(),
        name_attr: a.name_attr.clone(),
        blocks: a.blocks,
    };
    let to = spec.format.unwrap_or(OutFormat::Bed);
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let ro = ReaderOptions {
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: match a.format.format.as_deref() {
            None => None,
            Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: user_span(a.vcf_span.vcf_span),
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    // Open and check the input before the output is created.
    let input = open_reader(&a.file, &ro)?;
    check_convert(input.format(), &o)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let out = open_output(&spec, to, &oo)?;
    let mut w = if to == OutFormat::Tsv {
        open_tsv_writer(out, Format::Bed, &[], TsvHeader::None, oo.write_index)
    } else {
        let wo = WriterOptions {
            write_header: false,
            command_line: Some(ctx.command_line.clone()),
            version: ctx.version,
        };
        open_record_writer(out, Format::Bed, to, &Header::None, &wo)
    }
    .map_err(|e| e.with_output_name(&dest))?;
    // Records are converted in input order: the table only names chromosomes.
    let mut chroms = ChromTable::first_appearance();
    let named = |e: IoError| unsortable(e.with_output_name(&dest), a.out.output.is_some());
    convert(input, &mut chroms, &o, w.as_mut()).map_err(named)?;
    w.finish().map_err(named)?;
    Ok(())
}

/// The index of an out-of-order conversion cannot be built, and the file it
/// names is gone: say how to get a sorted, indexed result instead.
fn unsortable(e: IoError, to_file: bool) -> IoError {
    match e {
        IoError::Order(mut o) if to_file => {
            o.hint = "the output is out of order (an unsorted input, or allele footprints that \
start after the next record); write to stdout and pipe through gloci sort, then gloci index, \
or drop --write-index"
                .into();
            IoError::Order(o)
        }
        e => e,
    }
}
