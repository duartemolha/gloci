//! `gloci windows`: cut chromosomes or intervals into windows.
use super::user_span;
use crate::Ctx;
use crate::cli::WindowsArgs;
use gloci_core::ChromTable;
use gloci_io::{
    Format, Header, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, TsvHeader, UsageError,
    WriterOptions, open_output, open_reader, open_record_writer, open_tsv_writer,
    parse_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::windows::{Split, WindowsOptions, windows_genome, windows_intervals};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &WindowsArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let split = match (a.size, a.count) {
        (Some(size), None) => Split::Size {
            size,
            step: a.step.unwrap_or(size),
        },
        (None, Some(n)) => Split::Count(n),
        _ => unreachable!("clap requires exactly one of -w and -n"),
    };
    let o = WindowsOptions {
        split,
        reverse: a.reverse,
        id_column: a.id_column,
    };
    let to = spec.format.unwrap_or(OutFormat::Bed);
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    // Resolve inputs before the output is created.
    let genome = match a.genome.as_deref() {
        Some(g) => Some(resolve_genome(g)?.0),
        None => None,
    };
    let input = match a.intervals.as_deref() {
        Some(f) => {
            let ro = ReaderOptions {
                accept: gloci_io::Accept::default(),
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
            Some(open_reader(f, &ro)?)
        }
        None => None,
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = spec.path.as_ref().map_or_else(
        || "<stdout>".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
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
    match (genome, input) {
        (Some(g), None) => windows_genome(&g, &o, w.as_mut()),
        (None, Some(r)) => {
            // Records are cut in input order: the table only names chromosomes.
            let mut chroms = ChromTable::first_appearance();
            windows_intervals(r, &mut chroms, &o, w.as_mut())
        }
        _ => unreachable!("clap requires exactly one of -g and -i"),
    }
    .map_err(|e| e.with_output_name(&dest))?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
