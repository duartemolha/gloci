//! `gloci convert`: VCF/BCF, GFF3/GTF and BED12 to BED, alignments to BED, BEDPE or FASTQ,
//! features or BEDPE pairs to BAM, and bedGraph, WIG, BED and bigWig to bigWig and bigBed.
use super::user_span;
use crate::Ctx;
use crate::cli::ConvertArgs;
use gloci_core::{ChromTable, Genome, VariantKind};
use gloci_io::{
    BbiKind, BbiSizes, BbiWriteOptions, BbiWriter, BedpeReader, Format, GateKind, Header,
    InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions, RecordReader, TsvHeader,
    UsageError, WriterOptions, open_alignment_writer, open_aln_source, open_output, open_reader,
    open_reader_noting_empty, open_record_writer, open_tsv_writer, parse_output_spec,
    parse_track_output_spec, refuse_binary_to_tty, resolve_genome,
};
use gloci_ops::convert::{
    BamOptions, BedpeOptions, ConvertOptions, FastqOptions, FastqSink, alignments_to_bedpe,
    alignments_to_fastq, bam_header, bedpe_to_bam, check_convert, check_track_input, convert,
    convert_to_track, features_to_bam,
};
use std::io::{IsTerminal, Write};

pub fn run(ctx: &Ctx, a: &ConvertArgs) -> anyhow::Result<()> {
    let spec =
        match parse_track_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref()) {
            Ok(spec) => spec,
            // gloci writes no CRAM; for features, the alignments convert can write are --bam's.
            Err(e) if !(a.bam || a.fastq || a.bedpe) && names_cram(a) => {
                return Err(match input_format(ctx, a) {
                    Some(f) if is_feature(f) => UsageError(format!(
                        "gloci does not write CRAM; to write the features of {} as alignments, add \
--bam -g GENOME and name the output .bam (bedtools bedtobam), then convert it with `samtools \
view -C -T ref.fa`",
                        if a.file == "-" { "stdin" } else { &a.file }
                    )),
                    _ => e,
                }
                .into());
            }
            Err(e) => return Err(e.into()),
        };
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
    let track = match spec.format {
        Some(OutFormat::BigWig) => Some(BbiKind::BigWig),
        Some(OutFormat::BigBed) => Some(BbiKind::BigBed),
        _ => None,
    };
    if a.autosql.is_some() && track != Some(BbiKind::BigBed) {
        return Err(UsageError(
            "--autosql describes the columns of a bigBed: name the output .bb or .bigbed (or -O bigbed)"
                .into(),
        )
        .into());
    }
    if a.genome.is_some() && !a.bam && track.is_none() {
        return Err(UsageError(
            "-g applies with --bam, or for a bigWig or bigBed output (-o x.bw or -o x.bb); drop it, or add one of them"
                .into(),
        )
        .into());
    }
    if let Some(kind) = track {
        return run_track(ctx, a, &spec, kind, score_tag);
    }
    // One output mode per run: each takes the input as alignments and writes its own format.
    if a.bedpe {
        return run_bedpe(ctx, a, &o, &spec, score_tag);
    }
    if a.fastq {
        return run_fastq(ctx, a, &spec);
    }
    if a.bam {
        return run_bam(ctx, a, &spec);
    }
    let to = spec.format.unwrap_or(OutFormat::Bed);
    if !matches!(to, OutFormat::Bed | OutFormat::BedGraph | OutFormat::Tsv) {
        return Err(UsageError(wrong_output(to, &a.file, input_format(ctx, a))).into());
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

/// `convert` to a bigWig or bigBed (D171): bedGraph, a 4-column BED, WIG or bigWig, or BED or
/// bigBed, with the chromosome sizes of `-g` or of the bigWig or bigBed being converted.
fn run_track(
    ctx: &Ctx,
    a: &ConvertArgs,
    spec: &gloci_io::OutputSpec,
    kind: BbiKind,
    score_tag: Option<[u8; 2]>,
) -> anyhow::Result<()> {
    let out = match kind {
        BbiKind::BigWig => "bigWig",
        BbiKind::BigBed => "bigBed",
    };
    if a.out.write_index {
        return Err(UsageError(format!(
            "--write-index does not apply to a {out} output: a bigWig or bigBed carries its own \
index; drop --write-index"
        ))
        .into());
    }
    for (set, flag) in [
        (a.bam, "--bam"),
        (a.fastq, "--fastq"),
        (a.bedpe, "--bedpe"),
        (a.fastq2.is_some(), "--fastq2"),
    ] {
        if set {
            return Err(UsageError(format!(
                "{flag} chooses another output format than the {out} that -o names; drop {flag}, \
or name an output it writes"
            ))
            .into());
        }
    }
    for (set, flag) in [
        (a.split_alleles, "--split-alleles"),
        (!a.types.is_empty(), "--types"),
        (!a.feature.is_empty(), "--feature"),
        (a.name_attr.is_some(), "--name-attr"),
        (a.blocks, "--blocks"),
        (a.bed12, "--bed12"),
        (a.color.is_some(), "--color"),
        (score_tag.is_some(), "--score-tag"),
        (a.cigar, "--cigar"),
        (a.split_deletions, "--split-deletions"),
        (a.mate1_first, "--mate1-first"),
        (a.tags, "--tags"),
        (a.mapq.is_some(), "--mapq"),
    ] {
        if set {
            return Err(UsageError(format!(
                "{flag} does not apply to a {out} output: a {out} holds the records as they are; \
drop {flag}"
            ))
            .into());
        }
    }
    let autosql = match &a.autosql {
        None => None,
        Some(p) => Some(std::fs::read_to_string(p).map_err(|e| {
            UsageError(format!(
                "cannot read the --autosql file {}: {e}; check the --autosql path",
                gloci_io::memio::display_name(&p.to_string_lossy())
            ))
        })?),
    };
    let ro = reader_options(ctx, a, score_tag, false)?;
    // Open and check the input before the output is created.
    let (input, empty) = open_reader_noting_empty(&a.file, &ro)?;
    gloci_io::warn_unused_reference(a.align.reference.as_deref(), input.format() == Format::Cram);
    check_track_input(input.format(), kind).map_err(|e| {
        if empty {
            IoError::empty_input(input.name())
        } else {
            e.into()
        }
    })?;
    let sizes = match a.genome.as_deref() {
        Some(g) => BbiSizes {
            genome: resolve_genome(g)?.0,
            source: format!("-g {g}"),
        },
        None if matches!(input.format(), Format::BigWig | Format::BigBed) => {
            let map = ro.chrom_map.as_deref();
            let entries = gloci_io::bbi_chrom_sizes(&a.file)?
                .into_iter()
                .map(|(name, len)| {
                    let name = map.map_or(name.clone(), |m| m.rename(&name).to_vec());
                    (name, Some(len))
                })
                .collect();
            BbiSizes {
                genome: Genome::new(entries).map_err(|e| {
                    UsageError(format!(
                        "the chromosome list of {} cannot be used as sizes: {e}",
                        a.file
                    ))
                })?,
                source: format!("the header of {}", a.file),
            }
        }
        None => {
            return Err(UsageError(format!(
                "a {out} needs the chromosome sizes: add -g GENOME (a built-in genome such as \
hg38, or a .genome, .fai or .dict file); only a bigWig or bigBed input brings its own"
            ))
            .into());
        }
    };
    // A bigBed converted to a bigBed keeps its autoSql and its definedFieldCount, unless
    // --autosql replaces them (D171).
    let (autosql, autosql_name, defined_field_count, autosql_inherited) =
        match (autosql, input.format()) {
            (None, Format::BigBed) if kind == BbiKind::BigBed => {
                let (text, defined) = gloci_io::bbi_bed_columns(&a.file)?;
                (
                    text,
                    format!("the autoSql of {}", a.file),
                    (defined >= 3).then_some(usize::from(defined)),
                    true,
                )
            }
            (text, _) => (
                text,
                a.autosql.as_ref().map_or_else(String::new, |p| {
                    gloci_io::memio::display_name(&p.to_string_lossy())
                }),
                None,
                false,
            ),
        };
    let dest = spec.path.as_ref().map_or_else(
        || "<stdout>".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
    let mut w = BbiWriter::create(
        spec,
        BbiWriteOptions {
            kind,
            sizes,
            threads: ctx.threads,
            autosql,
            autosql_name,
            defined_field_count,
            autosql_inherited,
        },
    )
    .map_err(|e| e.with_output_name(&dest))?;
    let mut chroms = ChromTable::first_appearance();
    convert_to_track(input, &mut chroms, &mut w).map_err(|e| e.with_output_name(&dest))?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}

/// Whether `-o` or `-O` names CRAM.
fn names_cram(a: &ConvertArgs) -> bool {
    let ext = a
        .out
        .output
        .as_deref()
        .and_then(|p| p.extension())
        .is_some_and(|e| e.eq_ignore_ascii_case("cram"));
    let flag = a
        .out
        .output_format
        .as_deref()
        .is_some_and(|s| s.trim_start_matches('.').eq_ignore_ascii_case("cram"));
    ext || flag
}

/// The formats `convert --bam` reads as features.
fn is_feature(f: Format) -> bool {
    matches!(
        f,
        Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf | Format::Vcf
    )
}

/// The detected format of the input, for a message about the output that depends on it;
/// `None` when the input cannot be opened (the message then does not depend on it). Only
/// called on the way to an error, since it reads the start of the input.
fn input_format(ctx: &Ctx, a: &ConvertArgs) -> Option<Format> {
    let ro = reader_options(ctx, a, None, false).ok()?;
    open_reader_noting_empty(&a.file, &ro)
        .ok()
        .map(|(r, _)| r.format())
}

/// The message for an output format plain `convert` (no --bam, --fastq or --bedpe) does not
/// write: it writes BED, and the option that writes `to` depends on what the input is.
fn wrong_output(to: OutFormat, arg: &str, input: Option<Format>) -> String {
    let pick = "pick an output ending in .bed or .bed.gz (or -O bed)";
    let alignments = input.is_some_and(Format::is_alignment);
    let features = input.is_some_and(is_feature);
    // `arg` as the text names it; the command line names standard input `-`.
    let name = if arg == "-" { "stdin" } else { arg };
    let what = |f: Format| format!("{name} is {}", f.name());
    match (to, input) {
        (OutFormat::Bam, _) if features => format!(
            "convert writes BED unless --bam asks for BAM: add --bam -g GENOME to write the \
features of {name} as BAM records (bedtools bedtobam)"
        ),
        (OutFormat::Sam, _) if features => format!(
            "convert writes BED, or BAM with --bam -g GENOME, not SAM: add --bam -g GENOME and \
name the output .bam to write the features of {name} as alignments (gloci view x.bam -O sam \
shows it as SAM)"
        ),
        (OutFormat::Bam | OutFormat::Sam, _) if alignments => format!(
            "convert writes BED, not {}: {pick}; gloci view converts between alignment formats \
(gloci view {arg} -o x.bam)",
            to.name()
        ),
        (OutFormat::Fastq, _) if alignments => format!(
            "convert writes BED unless --fastq asks for FASTQ: add --fastq to write the reads of \
{name} as FASTQ (bedtools bamtofastq)"
        ),
        (OutFormat::Bedpe, _) if alignments => format!(
            "convert writes BED unless --bedpe asks for BEDPE: add --bedpe to write the read \
pairs of {name} as BEDPE (bedtools bamtobed -bedpe)"
        ),
        (OutFormat::Bam | OutFormat::Sam, Some(f)) => format!(
            "convert writes BED, not {}: {pick}; --bam writes BAM from BED, GFF or VCF features, \
and {}",
            to.name(),
            what(f)
        ),
        (OutFormat::Fastq | OutFormat::Bedpe, Some(f)) => format!(
            "convert writes {} only from alignments (BAM, CRAM or SAM), and {}: {pick}",
            to.name(),
            what(f)
        ),
        (OutFormat::Bam | OutFormat::Sam, None) => format!(
            "convert writes BED, not {}: {pick}; --bam -g GENOME writes BAM from BED, GFF or VCF \
features, and gloci view converts between alignment formats",
            to.name()
        ),
        (OutFormat::Fastq, None) => {
            format!("convert writes BED, not FASTQ: {pick}; --fastq writes FASTQ from alignments")
        }
        (OutFormat::Bedpe, None) => {
            format!("convert writes BED, not BEDPE: {pick}; --bedpe writes BEDPE from alignments")
        }
        _ => format!("convert writes BED, not {}: {pick}", to.name()),
    }
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
    /// A write to `second` failed, so the error names the `--fastq2` file.
    second_failed: bool,
}

impl FastqSink for Fastq {
    fn first(&mut self, record: &[u8]) -> std::io::Result<()> {
        self.first.write_all(record)
    }
    fn second(&mut self, record: &[u8]) -> std::io::Result<()> {
        match &mut self.second {
            Some(b) => b
                .write_all(record)
                .inspect_err(|_| self.second_failed = true),
            None => self.first.write_all(record),
        }
    }
}

/// Whether the outputs `a` and `b` are one output, whatever their spelling (`same.fq`,
/// `./same.fq`, `d/../same.fq`): both stdout; for a file (written to a temp file and renamed
/// into place), the same name in the same directory, resolved, since the file may not exist
/// yet and two renames onto one name keep only the second; for a device or pipe written
/// directly, the same file. A second hard link or a symlink to a file is another name, which
/// its own rename replaces, so it is another output.
fn same_output(a: &gloci_io::OutputSpec, b: &gloci_io::OutputSpec) -> bool {
    use std::path::{Path, PathBuf};
    let (pa, pb) = match (&a.path, &b.path) {
        (None, None) => return true,
        (Some(pa), Some(pb)) => (pa, pb),
        _ => return false,
    };
    if pa == pb {
        return true;
    }
    if a.direct || b.direct {
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            if let (Ok(ma), Ok(mb)) = (std::fs::metadata(pa), std::fs::metadata(pb)) {
                return ma.dev() == mb.dev() && ma.ino() == mb.ino();
            }
        }
        return false;
    }
    let landing = |p: &Path| -> Option<PathBuf> {
        let name = p.file_name()?;
        let dir = match p.parent() {
            Some(d) if !d.as_os_str().is_empty() => d,
            _ => Path::new("."),
        };
        Some(std::fs::canonicalize(dir).ok()?.join(name))
    };
    matches!((landing(pa), landing(pb)), (Some(x), Some(y)) if x == y)
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
    // Both ends go to one output when both name it: stdout, or the same file by any spelling.
    let one_output = spec2.as_ref().is_some_and(|s2| same_output(spec, s2));
    let shared = one_output && spec2.as_ref().is_some_and(|s2| s2.bgzf == spec.bgzf);
    if one_output && !shared {
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
    let dest2 = spec2.as_ref().map(name_of);
    let mut sink = Fastq {
        first,
        second,
        second_failed: false,
    };
    let o = FastqOptions {
        paired: spec2.is_some(),
        tags: a.tags,
    };
    alignments_to_fastq(src, &o, &mut sink).map_err(|e| match &dest2 {
        Some(d2) if sink.second_failed => e.with_output_name(d2),
        _ => e.with_output_name(&dest),
    })?;
    // The two files are published together: a failure in finishing either leaves neither in
    // place (and the previous ones as they were).
    let mut staged = gloci_io::StagedSet::default();
    sink.first
        .finish_staged(&mut staged)
        .map_err(|e| e.with_output_name(&dest))?;
    if let (Some(second), Some(d2)) = (sink.second, &dest2) {
        second
            .finish_staged(&mut staged)
            .map_err(|e| e.with_output_name(d2))?;
    }
    staged.publish()?;
    Ok(())
}

/// What `convert --bam` reads: features, or BEDPE pairs.
enum BamInput {
    Features(Box<dyn RecordReader>),
    Pairs(Box<BedpeReader>),
}

/// `convert --bam` (D167): BED, GFF or VCF features, or BEDPE pairs, as BAM records.
fn run_bam(ctx: &Ctx, a: &ConvertArgs, spec: &gloci_io::OutputSpec) -> anyhow::Result<()> {
    if let Some(to) = spec.format
        && to != OutFormat::Bam
    {
        return Err(UsageError(format!(
            "convert --bam writes BAM, not {}: pick an output ending in .bam (or -O bam); \
gloci view x.bam -O sam shows it as SAM",
            to.name()
        ))
        .into());
    }
    let format = match a.format.format.as_deref() {
        None => None,
        Some(f) => Some(f.parse::<Format>().map_err(UsageError)?),
    };
    let pairs = format == Some(Format::Bedpe);
    let lower = a.file.to_ascii_lowercase();
    if format.is_none() && (lower.ends_with(".bedpe") || lower.ends_with(".bedpe.gz")) {
        return Err(UsageError(format!(
            "{} is named as BEDPE, but BEDPE is read only when asked for: add --format bedpe to \
write two records per pair (bedtools bedpetobam), or --format bed to read it as BED",
            a.file
        ))
        .into());
    }
    if pairs && a.bed12 {
        return Err(UsageError(
            "--bed12 reads the blocks of BED12 records, and BEDPE has none; drop --bed12".into(),
        )
        .into());
    }
    let g = a
        .genome
        .as_deref()
        .ok_or_else(|| UsageError("--bam needs -g GENOME for the BAM header".into()))?;
    let (genome, source) = resolve_genome(g)?;
    let header = bam_header(&genome, g).map_err(UsageError)?;
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Bam,
        true,
    )?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: format.filter(|_| !pairs),
        // bedtools' footprint of a VCF record is its REF bases.
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    // Open and check the input before the output is created.
    let input = if pairs {
        BamInput::Pairs(Box::new(BedpeReader::open(&a.file, &ro)?))
    } else {
        let r = open_reader(&a.file, &ro).map_err(|e| match e {
            // `--bam` writes alignments from features; from alignments, that is view's work.
            IoError::Gate {
                kind: GateKind::NoAlignments,
                ..
            } => UsageError(format!(
                "convert --bam writes BAM from BED, GFF or VCF features, and {} holds alignments \
(BAM, CRAM or SAM); gloci view converts between alignment formats: gloci view {} -o x.bam",
                if a.file == "-" { "stdin" } else { &a.file },
                a.file
            ))
            .into(),
            e => e,
        })?;
        match r.format() {
            Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf | Format::Vcf => {}
            f => {
                return Err(UsageError(format!(
                    "convert --bam reads BED, GFF or VCF, or BEDPE with --format bedpe; {} is {}",
                    r.name(),
                    f.name()
                ))
                .into());
            }
        }
        BamInput::Features(r)
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
    let named = |e: IoError| e.with_output_name(&dest);
    let out = open_output(spec, OutFormat::Bam, &oo).map_err(named)?;
    let wo = WriterOptions {
        write_header: true,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    // The records are SAM lines, which the writer's encoder turns into BAM.
    let mut w =
        open_alignment_writer(out, Format::Sam, OutFormat::Bam, &header, &wo).map_err(named)?;
    let mut chroms = ChromTable::first_appearance();
    let o = BamOptions {
        mapq: a.mapq.unwrap_or(255),
        bed12: a.bed12,
    };
    match input {
        BamInput::Features(r) => {
            features_to_bam(r, &mut chroms, &header, &source, &o, w.as_mut()).map_err(named)?
        }
        BamInput::Pairs(r) => {
            bedpe_to_bam(*r, &mut chroms, &header, &source, o.mapq, w.as_mut()).map_err(named)?
        }
    };
    w.finish().map_err(named)?;
    Ok(())
}
