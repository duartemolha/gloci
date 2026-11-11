//! `gloci nuc`: base composition under each interval.
use super::user_span;
use crate::Ctx;
use crate::cli::NucArgs;
use gloci_core::ChromTable;
use gloci_io::{
    Format, IndexedFasta, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions,
    UsageError, open_output, open_reader, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::nuc::{NucOptions, nuc};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &NucArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref()).map_err(
        |e| match (
            e.0.starts_with("cannot tell the output format"),
            &a.out.output,
        ) {
            (true, Some(p)) => UsageError(format!(
                "cannot tell the output format from '{}'; nuc writes .tsv or .txt (optionally \
with .gz for BGZF), or pass -O",
                p.display()
            )),
            _ => e,
        },
    )?;
    if let Some(f) = spec.format
        && !matches!(f, OutFormat::Tsv | OutFormat::Bed)
    {
        return Err(UsageError(format!(
            "nuc writes tab-separated text, not {}; name the output like -o nuc.tsv or \
-o nuc.tsv.gz",
            f.name()
        ))
        .into());
    }
    if a.out.write_index {
        return Err(UsageError(
            "--write-index does not apply to nuc, which writes a table, not intervals".into(),
        )
        .into());
    }
    let pattern = match a.pattern.as_deref() {
        Some("") => {
            return Err(UsageError(
                "--pattern is empty; give the bases to count, e.g. --pattern CG".into(),
            )
            .into());
        }
        Some(p) => Some(p.as_bytes().to_vec()),
        None => None,
    };
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Tsv,
        spec.bgzf,
    )?;
    let o = NucOptions {
        strand: a.strand,
        seq: a.seq,
        pattern,
        ignore_case: a.ignore_case,
        skip_outside: a.invalid.skip_invalid,
    };
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
    // Open the inputs before the output is created.
    let input = open_reader(&a.file, &ro)?;
    let mut fasta = IndexedFasta::open(&a.fasta)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: false,
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    let mut chroms = ChromTable::first_appearance();
    let named = |e: IoError| e.with_output_name(&dest);
    nuc(input, &mut fasta, &mut chroms, &o, &mut out).map_err(named)?;
    out.finish().map_err(named)?;
    Ok(())
}
