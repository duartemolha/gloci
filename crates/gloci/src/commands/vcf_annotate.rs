//! `gloci vcf annotate` (P3b Task 22).

use crate::Ctx;
use crate::cli::VcfAnnotateArgs;
use gloci_io::{
    Format, Header, InputOptions, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions,
    UsageError, VcfHeader, WriterOptions, chrom_table_lenient, open_input, open_output,
    open_reader, open_reader_from, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{AnnotationSource, is_bed_name, vcf_annotate, vcf_annotate_output, vcf_header};
use gloci_vcf::annotate::{AnnotateOptions, Annotator, SourceKind};
use std::io::{BufRead, IsTerminal};

/// The bcftools annotate options gloci does not offer, each a usage error.
fn refuse_unsupported(a: &VcfAnnotateArgs) -> anyhow::Result<()> {
    let no = |opt: &str, why: &str| -> anyhow::Result<()> {
        Err(UsageError(format!(
            "bcftools annotate {opt} is not available in gloci vcf annotate: {why}"
        ))
        .into())
    };
    if a.rename_chrs.is_some() {
        return no(
            "--rename-chrs",
            "rename the chromosomes first with --chrom-map, e.g. gloci view FILE --chrom-map MAP | gloci vcf annotate - ...",
        );
    }
    if a.samples.is_some() || a.samples_file.is_some() {
        return no(
            "-s/--samples and -S/--samples-file",
            "annotate acts on sites only; to drop samples first use gloci vcf filter --samples or --exclude-samples",
        );
    }
    if a.include.is_some() || a.exclude.is_some() {
        return no(
            "-i/--include and -e/--exclude",
            "select sites first by piping through gloci vcf filter, e.g. gloci vcf filter FILE --include 'QUAL>20' | gloci vcf annotate - ...",
        );
    }
    if a.merge_logic.is_some() {
        return no(
            "-l/--merge-logic",
            "the first matching annotation is used for each column",
        );
    }
    if a.single_overlaps {
        return no(
            "--single-overlaps",
            "gloci always reads the annotation file as a stream with the rows that can still overlap",
        );
    }
    if a.rename_annots.is_some() {
        return no(
            "--rename-annots",
            "transfer a renamed tag with --columns NEW:=OLD from a VCF source, or remove one with --remove",
        );
    }
    if a.pair_logic.is_some() {
        return no(
            "--pair-logic",
            "gloci pairs records the way bcftools does by default (--pair-logic some: they must share an allele)",
        );
    }
    if a.columns_file.is_some() {
        return no("-C/--columns-file", "give the columns with --columns");
    }
    if a.header_line.is_some() {
        return no(
            "-H/--header-line",
            "put the header lines in a file and pass --header-lines FILE",
        );
    }
    if a.keep_sites {
        return no(
            "-k/--keep-sites",
            "it goes with -i/-e, which are not available either; pipe through gloci vcf filter",
        );
    }
    if a.regions.is_some() || a.regions_file.is_some() {
        return no(
            "-r/--regions and -R/--regions-file",
            "select regions first, e.g. gloci view FILE REGION | gloci vcf annotate - ... (or gloci view --regions-file BED)",
        );
    }
    Ok(())
}

/// The lines of `--header-lines FILE`, blank lines skipped.
fn read_header_lines(path: &std::path::Path) -> anyhow::Result<Vec<Vec<u8>>> {
    let text = std::fs::read(path).map_err(|source| IoError::Open {
        path: path.display().to_string(),
        source,
    })?;
    Ok(text
        .split(|&b| b == b'\n')
        .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
        .filter(|l| !l.is_empty())
        .map(<[u8]>::to_vec)
        .collect())
}

pub fn run(ctx: &Ctx, a: &VcfAnnotateArgs) -> anyhow::Result<()> {
    refuse_unsupported(a)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let header_lines = match &a.header_lines {
        Some(p) => read_header_lines(p)?,
        None => Vec::new(),
    };
    let input = open_reader(&a.file, &ro)?;
    let file_header = vcf_header("vcf annotate", input.as_ref())?;

    // What kind of file --annotations names: a VCF/BCF by its content, anything
    // else is a tab or BED file.
    let mut src_reader = None;
    let mut tab = None;
    if let Some(path) = &a.annotations {
        if path == "-" {
            return Err(UsageError(
                "--annotations - is not supported: the annotation file is read alongside FILE, so it needs a path (FILE itself can be -)"
                    .to_string(),
            )
            .into());
        }
        // Open once and look at the first bytes without consuming them, so a
        // pipe or process substitution keeps all its data.
        let input = open_input(
            path,
            &InputOptions {
                threads: ctx.threads,
                format: Some(Format::Bed),
            },
        )?;
        let mut reader = input.reader;
        let head = reader.fill_buf().map_err(|source| IoError::Read {
            path: path.clone(),
            source,
        })?;
        let kind = if head.starts_with(b"BCF\x02") {
            Some(Format::Bcf)
        } else if head.starts_with(b"##fileformat=VCF") || head.starts_with(b"#CHROM") {
            Some(Format::Vcf)
        } else {
            None
        };
        match kind {
            Some(format) => {
                let o = ReaderOptions {
                    format: Some(format),
                    ..ro.clone()
                };
                src_reader = Some(open_reader_from(path, Box::new(reader), &o)?);
            }
            None => tab = Some((reader, path.clone(), is_bed_name(path))),
        }
    }
    let src_header = match &src_reader {
        Some(r) => Some(vcf_header("vcf annotate", r.as_ref())?),
        None => None,
    };
    let kind = match (&src_header, &tab) {
        (Some(h), _) => SourceKind::Vcf(&h.raw),
        (None, Some((_, _, bed))) => SourceKind::Tab { bed: *bed },
        (None, None) => SourceKind::None,
    };
    let options = AnnotateOptions {
        columns: a.columns.clone(),
        remove: a.remove.clone(),
        header_lines,
        set_id: a.set_id.clone(),
        mark_sites: a.mark_sites.clone(),
        min_overlap: a.min_overlap.clone(),
    };
    let mut annotator = Annotator::new(&options, &file_header.raw, kind).map_err(UsageError)?;
    let out_text = annotator.out_header();
    let (to, header) = vcf_annotate_output(input.as_ref(), spec.format, &out_text)?;
    let Header::Vcf(out_header) = &header else {
        unreachable!("vcf_annotate_output returns a VCF header")
    };
    let out_header: VcfHeader = out_header.clone();
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let mut refs: Vec<&dyn gloci_io::RecordReader> = vec![input.as_ref()];
    if let Some(r) = &src_reader {
        refs.push(r.as_ref());
    }
    let mut chroms = chrom_table_lenient(None, &refs, a.chrom_order.chrom_order)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    // As bcftools annotate --no-version: no command line is added to the header.
    let wo = WriterOptions {
        write_header: true,
        command_line: None,
        version: ctx.version,
    };
    let dest = spec.path.as_ref().map_or_else(
        || "stdout".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
    let out = open_output(&spec, to, &oo)?;
    let record_format = if to == OutFormat::Bcf {
        Format::Bcf
    } else {
        Format::Vcf
    };
    let mut writer = open_record_writer(out, record_format, to, &header, &wo)
        .map_err(|e| e.with_output_name(&dest))?;
    let source = match (src_reader, tab) {
        (Some(r), _) => AnnotationSource::Vcf(r),
        (None, Some((t, name, _))) => AnnotationSource::Tab(t, name),
        (None, None) => AnnotationSource::None,
    };
    let sum = vcf_annotate(
        input,
        source,
        &mut chroms,
        &mut annotator,
        &out_header,
        writer.as_mut(),
    )
    .map_err(|e| e.with_output_name(&dest))?;
    writer.finish().map_err(|e| e.with_output_name(&dest))?;
    log::info!(
        "vcf annotate: {} records read, {} annotated, {} written; at most {} records held at once",
        sum.records_in,
        sum.matched,
        sum.records_out,
        sum.peak_buffered
    );
    Ok(())
}
