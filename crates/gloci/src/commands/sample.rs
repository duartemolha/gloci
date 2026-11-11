//! `gloci sample`: a random subset of the records of a file.
use crate::Ctx;
use crate::cli::SampleArgs;
use gloci_core::{ChromTable, Strand};
use gloci_io::{
    Format, InvalidPolicy, OutputOptions, ReaderOptions, UsageError, WriterOptions, open_reader,
    parse_output_spec,
};
use gloci_ops::rng::seed_or_entropy;
use gloci_ops::sample::{SampleOptions, sample};
use gloci_ops::{ReportOutput, open_report};

pub fn run(ctx: &Ctx, a: &SampleArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: match a.format.format.as_deref() {
            None => None,
            Some(f) => Some(f.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let reader = open_reader(&a.file, &ro)?;
    let o = SampleOptions {
        count: a.count,
        seed: seed_or_entropy(a.seed.seed),
        strand: a.strand.as_deref().map(|s| {
            if s == "+" {
                Strand::Plus
            } else {
                Strand::Minus
            }
        }),
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
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
    // Records are written in input order, so the table only gives names ids.
    let mut chroms = ChromTable::first_appearance();
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    sample(reader, &mut chroms, &o, w.as_mut()).map_err(|e| e.with_output_name(&dest))?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
