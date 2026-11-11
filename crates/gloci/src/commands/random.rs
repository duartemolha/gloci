//! `gloci random`: random intervals over a genome.
use crate::Ctx;
use crate::cli::RandomArgs;
use gloci_io::{
    Format, Header, OutFormat, OutputOptions, TsvHeader, WriterOptions, open_output,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::random::{RandomOptions, random_intervals};
use gloci_ops::rng::seed_or_entropy;
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &RandomArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let to = spec.format.unwrap_or(OutFormat::Bed);
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    // Resolve the genome before the output is created.
    let (genome, _) = resolve_genome(&a.genome)?;
    let o = RandomOptions {
        length: a.length,
        count: a.count,
        seed: seed_or_entropy(a.seed.seed),
    };
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
    random_intervals(&genome, &o, ctx.threads, w.as_mut())
        .map_err(|e| e.with_output_name(&dest))?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
