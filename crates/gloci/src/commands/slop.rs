//! `gloci slop`: grow each feature on its sides. Also holds the plumbing
//! `gloci flank` shares.
use crate::Ctx;
use crate::cli::{PadArgs, SlopArgs};
use gloci_core::Genome;
use gloci_core::{ChromTable, Fraction};
use gloci_io::{
    Format, InvalidPolicy, IoError, OutputOptions, ReaderOptions, RecordReader, RecordWriter,
    UsageError, WriterOptions, open_reader, parse_output_spec, resolve_genome,
};
use gloci_ops::transform::{Amount, Pad, check_format};
use gloci_ops::{ReportOutput, open_report};

/// `slop` and `flank` differ only in the operation.
pub type Op = fn(
    Box<dyn RecordReader>,
    &mut ChromTable,
    &Pad,
    Option<&Genome>,
    &mut dyn RecordWriter,
) -> Result<(), IoError>;

/// The arguments the two commands have in common.
pub struct Shared<'a> {
    pub name: &'static str,
    pub file: &'a str,
    pub pad: &'a PadArgs,
    pub genome: Option<&'a str>,
    pub format: Option<&'a str>,
    pub header: bool,
    pub skip_invalid: bool,
    pub chrom_map: &'a crate::cli::ChromMapArg,
    pub out: &'a crate::cli::OutputArgs,
}

fn amount(flag: &str, v: Option<&str>, pct: bool) -> Result<Amount, UsageError> {
    let Some(v) = v else {
        return Ok(Amount::Bases(0));
    };
    if pct {
        Fraction::parse_ratio(v)
            .map(Amount::Fraction)
            .map_err(|e| UsageError(format!("{flag}: {e}")))
    } else {
        v.parse::<u64>().map(Amount::Bases).map_err(|_| {
            UsageError(format!(
                "{flag}: \"{v}\" is not a whole number of bases; use --pct for a fraction"
            ))
        })
    }
}

pub fn pad(a: &PadArgs) -> Result<Pad, UsageError> {
    let (l, r) = match a.both.as_deref() {
        Some(b) => (Some(b), Some(b)),
        None => (a.left.as_deref(), a.right.as_deref()),
    };
    let flag = |own: &'static str| if a.both.is_some() { "-b" } else { own };
    Ok(Pad {
        left: amount(flag("--left"), l, a.pct)?,
        right: amount(flag("--right"), r, a.pct)?,
        strand_aware: a.strand,
    })
}

pub fn run_shared(ctx: &Ctx, s: &Shared, op: Op) -> anyhow::Result<()> {
    let pad = pad(s.pad)?;
    let spec = parse_output_spec(s.out.output.as_deref(), s.out.output_format.as_deref())?;
    let genome = match s.genome {
        Some(g) => Some(resolve_genome(g)?),
        None => None,
    };
    let ro = ReaderOptions {
        chrom_map: s.chrom_map.load()?,
        threads: ctx.threads,
        format: match s.format {
            None => None,
            Some(f) => Some(f.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: if s.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let reader = open_reader(s.file, &ro)?;
    check_format(s.name, reader.format())?;
    // Records are transformed in input order, so the table only gives names
    // ids; with -g it is the genome's (a chromosome outside it is an error).
    let mut chroms = match &genome {
        Some((g, source)) => ChromTable::closed(g, source.clone()),
        None => ChromTable::first_appearance(),
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: s.out.compression_level,
        write_index: s.out.write_index,
    };
    let wo = WriterOptions {
        write_header: s.header,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let mut w = open_report(
        reader.as_ref(),
        &ReportOutput::SameAsA,
        &spec,
        &oo,
        &wo,
        s.header,
    )?;
    op(
        reader,
        &mut chroms,
        &pad,
        genome.as_ref().map(|g| &g.0),
        w.as_mut(),
    )?;
    w.finish()?;
    Ok(())
}

pub fn run(ctx: &Ctx, a: &SlopArgs) -> anyhow::Result<()> {
    run_shared(
        ctx,
        &Shared {
            name: "slop",
            file: &a.file,
            pad: &a.pad,
            genome: a.genome.genome.as_deref(),
            format: a.format.format.as_deref(),
            header: a.header.header,
            skip_invalid: a.invalid.skip_invalid,
            chrom_map: &a.chrom_map,
            out: &a.out,
        },
        gloci_ops::transform::slop,
    )
}
