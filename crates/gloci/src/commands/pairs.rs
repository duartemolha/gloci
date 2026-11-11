//! `gloci pair-intersect` and `gloci pair-pair`: BEDPE pairs by how their ends overlap B.
use crate::Ctx;
use crate::cli::{OutputArgs, PairIntersectArgs, PairPairArgs, PairPairTypeArg, PairTypeArg};
use crate::commands::{dest_name, user_span};
use gloci_core::ChromTable;
use gloci_io::{
    BedpeReader, Format, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, UsageError,
    open_output, open_reader, parse_output_spec, tsv_output_ok,
};
use gloci_ops::pairs::{
    PairIntersectOptions, PairPairOptions, PairPairType, PairType, pair_intersect, pair_pair,
    validate,
};

fn usage(msg: impl Into<String>) -> anyhow::Error {
    UsageError(msg.into()).into()
}

fn kind(t: PairTypeArg) -> PairType {
    match t {
        PairTypeArg::Either => PairType::Either,
        PairTypeArg::Neither => PairType::Neither,
        PairTypeArg::Both => PairType::Both,
        PairTypeArg::Xor => PairType::Xor,
        PairTypeArg::Notboth => PairType::NotBoth,
        PairTypeArg::Ispan => PairType::InnerSpan,
        PairTypeArg::Ospan => PairType::OuterSpan,
        PairTypeArg::Notispan => PairType::NotInnerSpan,
        PairTypeArg::Notospan => PairType::NotOuterSpan,
    }
}

pub fn run_intersect(ctx: &Ctx, a: &PairIntersectArgs) -> anyhow::Result<()> {
    if a.o.split {
        return Err(usage(
            "--split has no meaning for pair-intersect: ends are whole intervals; drop it",
        ));
    }
    let o = PairIntersectOptions {
        kind: kind(a.kind),
        criteria: a.o.criteria(),
    };
    validate(&o)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    check_no_index(&a.out, "pair-intersect")?;
    let to = spec.format.unwrap_or(if o.kind.writes_a_alone() {
        OutFormat::Bedpe
    } else {
        OutFormat::Tsv
    });
    tsv_output_ok(Format::Bedpe, to)?;
    let chrom_map = a.chrom_map.load()?;
    let ro = |format| ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: chrom_map.clone(),
        threads: ctx.threads,
        format,
        vcf_span: user_span(a.vcf_span.vcf_span),
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let b_format = match a.b_format.as_deref() {
        None => None,
        Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
    };
    // BEDPE is never detected: A is read as BEDPE whatever it looks like.
    let reader_a = BedpeReader::open(&a.a, &ro(None))?;
    let reader_b = open_reader(&a.b, &ro(b_format))?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: false,
    };
    let dest = dest_name(&spec);
    let mut out = open_output(&spec, to, &oo).map_err(|e| e.with_output_name(&dest))?;
    let mut chroms = ChromTable::first_appearance();
    pair_intersect(reader_a, reader_b, &mut chroms, &o, ctx.threads, &mut out)
        .map_err(|e| e.with_output_name(&dest))?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

fn pair_kind(t: PairPairTypeArg) -> PairPairType {
    match t {
        PairPairTypeArg::Neither => PairPairType::Neither,
        PairPairTypeArg::Either => PairPairType::Either,
        PairPairTypeArg::Both => PairPairType::Both,
        PairPairTypeArg::Notboth => PairPairType::NotBoth,
    }
}

pub fn run_pair(ctx: &Ctx, a: &PairPairArgs) -> anyhow::Result<()> {
    if a.o.split {
        return Err(usage(
            "--split has no meaning for pair-pair: ends are whole intervals; drop it",
        ));
    }
    let o = PairPairOptions {
        kind: pair_kind(a.kind),
        criteria: a.o.criteria(),
        window: a.w.window(),
        different_names: a.require_different_names,
    };
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    check_no_index(&a.out, "pair-pair")?;
    let to = spec.format.unwrap_or(if o.kind.writes_a_alone() {
        OutFormat::Bedpe
    } else {
        OutFormat::Tsv
    });
    tsv_output_ok(Format::Bedpe, to)?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let reader_a = BedpeReader::open(&a.a, &ro)?;
    let reader_b = BedpeReader::open(&a.b, &ro)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: false,
    };
    let dest = dest_name(&spec);
    let mut out = open_output(&spec, to, &oo).map_err(|e| e.with_output_name(&dest))?;
    let mut chroms = ChromTable::first_appearance();
    pair_pair(reader_a, reader_b, &mut chroms, &o, ctx.threads, &mut out)
        .map_err(|e| e.with_output_name(&dest))?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

fn check_no_index(out: &OutputArgs, cmd: &str) -> anyhow::Result<()> {
    if out.write_index {
        return Err(usage(format!(
            "--write-index cannot index {cmd} output, because its columns are not the BED ones; \
write -o FILE.bedpe.gz and then run: gloci index -s 1 -b 2 -e 3 --zero-based FILE.bedpe.gz"
        )));
    }
    Ok(())
}
