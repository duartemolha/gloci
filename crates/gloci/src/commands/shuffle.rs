//! `gloci shuffle`: move every record to a random place in the genome.
use crate::Ctx;
use crate::cli::ShuffleArgs;
use gloci_core::ChromTable;
use gloci_io::{
    BedpeReader, Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    WriterOptions, open_output, open_reader, parse_output_spec, refuse_binary_to_tty,
    resolve_genome,
};
use gloci_ops::rng::seed_or_entropy;
use gloci_ops::shuffle::{ShuffleOptions, shuffle, shuffle_bedpe};
use gloci_ops::{ReportOutput, open_report};
use std::io::{IsTerminal, Write};

pub fn run(ctx: &Ctx, a: &ShuffleArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let format = match a.format.format.as_deref() {
        None => None,
        Some(f) => Some(f.parse::<Format>().map_err(UsageError)?),
    };
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let (genome, source) = resolve_genome(&a.genome)?;
    let mut chroms = ChromTable::declared(&genome, source);
    let o = ShuffleOptions {
        seed: seed_or_entropy(a.seed.seed),
        include: a.include.clone(),
        exclude: a.exclude.clone(),
        same_chrom: a.same_chrom,
        chrom_first: a.chrom_first,
        no_overlapping: a.no_overlapping,
        max_tries: a.max_tries,
        skip_unplaceable: a.skip_unplaceable,
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: false,
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    if format == Some(Format::Bedpe) {
        let to = spec.format.unwrap_or(OutFormat::Bedpe);
        if !matches!(to, OutFormat::Bedpe | OutFormat::Tsv) {
            return Err(UsageError(format!(
                "shuffle --format bedpe writes BEDPE: use -o FILE.bedpe or FILE.bedpe.gz (or \
-O bedpe), not {}",
                to.name()
            ))
            .into());
        }
        refuse_binary_to_tty(
            spec.path.is_none() && std::io::stdout().is_terminal(),
            to,
            spec.bgzf,
        )?;
        let reader = BedpeReader::open(&a.file, &ro)?;
        let mut out = open_output(&spec, to, &oo)?;
        if a.header.header {
            for line in reader.header() {
                out.write_all(line)
                    .and_then(|()| out.write_all(b"\n"))
                    .map_err(|source| gloci_io::IoError::Write {
                        path: dest.clone(),
                        source,
                    })?;
            }
        }
        shuffle_bedpe(reader, &mut chroms, &genome, &o, &mut out)
            .map_err(|e| e.with_output_name(&dest))?;
        out.finish().map_err(|e| e.with_output_name(&dest))?;
        return Ok(());
    }
    let reader = open_reader(&a.file, &ro)?;
    let wo = WriterOptions {
        write_header: a.header.header,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let mut w = open_report(
        reader.as_ref(),
        &ReportOutput::SameAsA,
        &spec,
        &oo,
        &wo,
        a.header.header,
    )?;
    shuffle(reader, &mut chroms, &genome, &o, ctx.threads, w.as_mut())
        .map_err(|e| e.with_output_name(&dest))?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
