//! `gloci getfasta`: the sequence under each interval.
use super::user_span;
use crate::Ctx;
use crate::cli::GetfastaArgs;
use gloci_core::ChromTable;
use gloci_io::{
    Format, IndexedFasta, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions,
    UsageError, open_output, open_reader, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::getfasta::{FastaName, GetfastaOptions, SeqLayout, getfasta};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &GetfastaArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref()).map_err(
        |e| match (
            e.0.starts_with("cannot tell the output format"),
            &a.out.output,
        ) {
            (true, Some(p)) => UsageError(format!(
                "cannot tell the output format from '{}'; getfasta writes .fa, .fasta, .fna, \
.tsv or .bed (each optionally with .gz for BGZF), or pass -O",
                p.display()
            )),
            _ => e,
        },
    )?;
    // The layout flags choose what is written; a name only chooses BGZF.
    if let Some(f) = spec.format
        && !matches!(f, OutFormat::Fasta | OutFormat::Tsv | OutFormat::Bed)
    {
        return Err(UsageError(format!(
            "getfasta writes FASTA, or with --tab / --bed-out tab-separated text, not {}; \
name the output like -o seqs.fa, -o seqs.fa.gz or -o seqs.tsv",
            f.name()
        ))
        .into());
    }
    if a.out.write_index {
        return Err(UsageError(
            "--write-index does not apply to getfasta; index the FASTA it wrote with: \
gloci index OUT.fa"
                .into(),
        )
        .into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Fasta,
        spec.bgzf,
    )?;
    let o = GetfastaOptions {
        strand: a.strand,
        name: if a.name {
            FastaName::NameCoords
        } else if a.name_only {
            FastaName::NameOnly
        } else {
            FastaName::Coords
        },
        split: a.split,
        layout: if a.tab {
            SeqLayout::Tab
        } else if a.bed_out {
            SeqLayout::BedOut
        } else {
            SeqLayout::Fasta
        },
        rna: a.rna,
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
    let dest = spec.path.as_ref().map_or_else(
        || "<stdout>".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
    let mut out = open_output(&spec, OutFormat::Fasta, &oo)?;
    let mut chroms = ChromTable::first_appearance();
    let named = |e: IoError| e.with_output_name(&dest);
    getfasta(input, &mut fasta, &mut chroms, &o, &mut out).map_err(named)?;
    out.finish().map_err(named)?;
    Ok(())
}
