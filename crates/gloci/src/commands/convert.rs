//! `gloci convert`: VCF/BCF, GFF3/GTF and BED12 to BED, and alignments to BED, BEDPE or FASTQ.
use super::user_span;
use crate::Ctx;
use crate::cli::ConvertArgs;
use gloci_core::{ChromTable, VariantKind};
use gloci_io::{
    Format, Header, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions, TsvHeader,
    UsageError, WriterOptions, open_aln_source, open_output, open_reader_noting_empty,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::convert::{
    BedpeOptions, ConvertOptions, FastqOptions, FastqSink, alignments_to_bedpe,
    alignments_to_fastq, check_convert, convert,
};
use std::io::{IsTerminal, Write};

pub fn run(ctx: &Ctx, a: &ConvertArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let score_tag = match a.score_tag.as_deref() {
        None => None,
        Some(t) => match t.as_bytes() {
            &[x, y] if x.is_ascii_alphabetic() && y.is_ascii_alphanumeric() => Some([x, y]),
            _ => {
                return Err(UsageError(format!(
                    "--score-tag {t} is not a tag name; use two characters such as NM or AS"
                ))
                .into());
            }
        },
    };
    let o = ConvertOptions {
        vcf_span: user_span(a.vcf_span.vcf_span),
        split_alleles: a.split_alleles,
        types: (!a.types.is_empty())
            .then(|| a.types.iter().map(|&k| VariantKind::from(k)).collect()),
        features: a.feature.clone(),
        name_attr: a.name_attr.clone(),
        blocks: a.blocks,
        bed12: a.bed12,
        color: a.color.clone(),
        score_tag,
        cigar: a.cigar,
        split_deletions: a.split_deletions,
        bedpe: a.bedpe,
        mate1_first: a.mate1_first,
    };
    // One output mode per run: each takes the input as alignments and writes its own format.
    if a.bedpe {
        return run_bedpe(ctx, a, &o, &spec, score_tag);
    }
    if a.fastq {
        return run_fastq(ctx, a, &spec);
    }
    let to = spec.format.unwrap_or(OutFormat::Bed);
    if !matches!(to, OutFormat::Bed | OutFormat::BedGraph | OutFormat::Tsv) {
        let view = if matches!(to, OutFormat::Bam | OutFormat::Sam) {
            "; gloci view converts between alignment formats (gloci view x.cram -o x.bam)"
        } else {
            ""
        };
        return Err(UsageError(format!(
            "convert writes BED, not {}: pick an output ending in .bed or .bed.gz (or -O bed){view}",
            to.name()
        ))
        .into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let ro = reader_options(ctx, a, score_tag, false)?;
    // Open and check the input before the output is created.
    let (input, empty) = open_reader_noting_empty(&a.file, &ro)?;
    gloci_io::warn_unused_reference(a.align.reference.as_deref(), input.format() == Format::Cram);
    // An empty file is detected as an empty BED: fine for --blocks (no rows), but an empty BAM
    // or SAM would be told "a BED input is already BED".
    check_convert(input.format(), &o).map_err(|e| {
        if empty {
            IoError::empty_input(input.name())
        } else {
            e.into()
        }
    })?;
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
    // Records are converted in input order: the table only names chromosomes.
    let mut chroms = ChromTable::first_appearance();
    let named = |e: IoError| unsortable(e.with_output_name(&dest), a.out.output.is_some());
    convert(input, &mut chroms, &o, w.as_mut()).map_err(named)?;
    w.finish().map_err(named)?;
    Ok(())
}

/// The reader options of `convert`.
fn reader_options(
    ctx: &Ctx,
    a: &ConvertArgs,
    score_tag: Option<[u8; 2]>,
    full_records: bool,
) -> anyhow::Result<ReaderOptions> {
    Ok(ReaderOptions {
        // `convert` reads alignments, PSL and RepeatMasker too (D100).
        accept: gloci_io::Accept {
            alignments: Some(gloci_io::AlignOptions {
                score_tag,
                ..a.align.options(full_records)?
            }),
            convert_only: true,
        },
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
    })
}

/// `convert --bedpe` (D156): BEDPE pairs from alignments sorted by name.
fn run_bedpe(
    ctx: &Ctx,
    a: &ConvertArgs,
    o: &ConvertOptions,
    spec: &gloci_io::OutputSpec,
    score_tag: Option<[u8; 2]>,
) -> anyhow::Result<()> {
    let to = spec.format.unwrap_or(OutFormat::Bedpe);
    if !matches!(to, OutFormat::Bedpe | OutFormat::Tsv) {
        return Err(UsageError(format!(
            "convert --bedpe writes BEDPE, not {}: pick an output ending in .bedpe or \
.bedpe.gz (or -O bedpe)",
            to.name()
        ))
        .into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let ro = reader_options(ctx, a, score_tag, false)?;
    // Open and check the input before the output is created.
    let src = open_aln_source(&a.file, &ro).map_err(|e| match e {
        IoError::Detect { path, msg } if msg.ends_with("is empty: it has no records") => {
            IoError::empty_input(&path)
        }
        e => e,
    })?;
    gloci_io::warn_unused_reference(a.align.reference.as_deref(), src.format() == Format::Cram);
    check_convert(src.format(), o)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = spec.path.as_ref().map_or_else(
        || "<stdout>".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
    let mut out = open_output(spec, to, &oo).map_err(|e| e.with_output_name(&dest))?;
    let bo = BedpeOptions {
        mate1_first: a.mate1_first,
        score_tag,
    };
    alignments_to_bedpe(src, &bo, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish().map_err(|e| e.with_output_name(&dest))?;
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

/// Where `convert --fastq` writes: one output (the reads, or both ends interleaved), or the
/// first ends in `first` and the second ends in `second`.
struct Fastq {
    first: gloci_io::Output,
    second: Option<gloci_io::Output>,
}

impl FastqSink for Fastq {
    fn first(&mut self, record: &[u8]) -> std::io::Result<()> {
        self.first.write_all(record)
    }
    fn second(&mut self, record: &[u8]) -> std::io::Result<()> {
        match &mut self.second {
            Some(b) => b.write_all(record),
            None => self.first.write_all(record),
        }
    }
}

/// `convert --fastq` (D166): the reads of alignments as FASTQ, in pairs with `--fastq2`.
fn run_fastq(ctx: &Ctx, a: &ConvertArgs, spec: &gloci_io::OutputSpec) -> anyhow::Result<()> {
    let spec2 = match &a.fastq2 {
        None => None,
        Some(p) => Some(parse_output_spec(Some(p), None)?),
    };
    for (s, flag) in [(Some(spec), "-o"), (spec2.as_ref(), "--fastq2")] {
        if let Some(to) = s.and_then(|s| s.format)
            && to != OutFormat::Fastq
        {
            return Err(UsageError(format!(
                "convert --fastq writes FASTQ, not {} (given with {flag}): pick an output ending \
in .fq, .fq.gz, .fastq or .fastq.gz (or -O fq)",
                to.name()
            ))
            .into());
        }
    }
    // Both ends go to one output when both name it: stdout, or the same file.
    let shared = spec2
        .as_ref()
        .is_some_and(|s2| s2.path == spec.path && s2.bgzf == spec.bgzf);
    if let Some(s2) = &spec2
        && s2.path == spec.path
        && !shared
    {
        return Err(UsageError(
            "-o and --fastq2 name the same file with different compression; use the same ending \
for both, or two files"
                .into(),
        )
        .into());
    }
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    for s in std::iter::once(spec).chain(spec2.as_ref()) {
        refuse_binary_to_tty(
            s.path.is_none() && std::io::stdout().is_terminal(),
            OutFormat::Fastq,
            s.bgzf,
        )?;
    }
    let ro = reader_options(ctx, a, None, true)?;
    // Open and check the input before the output is created.
    let src = open_aln_source(&a.file, &ro).map_err(|e| match e {
        IoError::Detect { path, msg } if msg.ends_with("is empty: it has no records") => {
            IoError::empty_input(&path)
        }
        e => e,
    })?;
    gloci_io::warn_unused_reference(a.align.reference.as_deref(), src.format() == Format::Cram);
    let name_of = |s: &gloci_io::OutputSpec| {
        s.path.as_ref().map_or_else(
            || "<stdout>".to_string(),
            |p| gloci_io::memio::display_name(&p.to_string_lossy()),
        )
    };
    let open = |s: &gloci_io::OutputSpec| {
        open_output(s, OutFormat::Fastq, &oo).map_err(|e| e.with_output_name(&name_of(s)))
    };
    let dest = name_of(spec);
    let first = open(spec)?;
    let second = match &spec2 {
        Some(s2) if !shared => Some(open(s2)?),
        _ => None,
    };
    let mut sink = Fastq { first, second };
    let o = FastqOptions {
        paired: spec2.is_some(),
        tags: a.tags,
    };
    alignments_to_fastq(src, &o, &mut sink).map_err(|e| e.with_output_name(&dest))?;
    sink.first.finish().map_err(|e| e.with_output_name(&dest))?;
    if let (Some(second), Some(s2)) = (sink.second, &spec2) {
        second
            .finish()
            .map_err(|e| e.with_output_name(&name_of(s2)))?;
    }
    Ok(())
}
