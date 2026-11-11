//! `gloci maskfasta`: mask the bases under intervals (D113).
use crate::Ctx;
use crate::cli::MaskfastaArgs;
use gloci_core::ChromTable;
use gloci_io::{
    FastaStream, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError, open_output,
    open_reader, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{Mask, MaskSet, PastEnd, maskfasta_with};
use std::io::IsTerminal;

/// The most chromosome names a warning lists.
const SHOWN: usize = 5;

pub fn run(ctx: &Ctx, a: &MaskfastaArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    if let Some(f) = spec.format
        && f != OutFormat::Fasta
    {
        return Err(UsageError(format!(
            "maskfasta writes FASTA but the output format is {}; use -o FILE.fa (or .fa.gz), \
or -O fa for standard output",
            f.name()
        ))
        .into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Fasta,
        spec.bgzf,
    )?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: None,
        vcf_span: super::user_span(a.vcf_span.vcf_span),
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let readers = a
        .masks
        .iter()
        .map(|f| open_reader(f, &ro))
        .collect::<Result<Vec<_>, _>>()?;
    let mut chroms = ChromTable::first_appearance();
    let masks = MaskSet::load(readers, &mut chroms)?;
    let fasta = FastaStream::open(&a.fasta, ctx.threads)?;
    let mask = if a.soft {
        Mask::Soft
    } else {
        Mask::Hard(a.mask_char)
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Fasta, &oo)?;
    let past_end = if a.invalid.skip_invalid {
        PastEnd::Clip
    } else {
        PastEnd::Error
    };
    let summary = maskfasta_with(fasta, &masks, mask, past_end, &mut out)
        .map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    if !summary.missing_chroms.is_empty() {
        let n = summary.missing_chroms.len();
        let mut names: Vec<String> = summary
            .missing_chroms
            .iter()
            .take(SHOWN)
            .map(|c| String::from_utf8_lossy(c).into_owned())
            .collect();
        if n > SHOWN {
            names.push(format!("and {} more", n - SHOWN));
        }
        log::warn!(
            "{n} chromosome{} with masks not in {}: {}; those masks were ignored (check the \
chromosome names, or rename them with --chrom-map)",
            if n == 1 { "" } else { "s" },
            a.fasta,
            names.join(", ")
        );
    }
    if let Some(first) = &summary.first_past_end {
        let n = summary.past_end;
        log::warn!(
            "{n} mask{} reached past the end of {} and {} clipped to it (a mask that starts at or \
past the end masks nothing); first: {first}",
            if n == 1 { "" } else { "s" },
            if n == 1 {
                "its sequence"
            } else {
                "their sequences"
            },
            if n == 1 { "was" } else { "were" },
        );
    }
    log::info!("masked {} bases", summary.masked_bases);
    Ok(())
}
