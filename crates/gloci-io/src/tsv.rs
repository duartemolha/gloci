//! The TSV results writer for the commands that append columns (`count`,
//! `pairs`, `map`, `coverage`, `closest`): an optional header, then rows the
//! caller has already joined (A's columns plus the appended ones).

use crate::detect::Format;
use crate::error::{IoError, UsageError};
use crate::index::TsvColumns;
use crate::output::{OutFormat, Output};
use crate::writer::{RecordWriter, TextWriter};
use std::io::Write;

/// What comes before the rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TsvHeader {
    None,
    /// A's own header lines (BED/GFF with `--header`), copied as they are;
    /// each line ends with a newline.
    Lines(Vec<u8>),
    /// A VCF A: the column-name line, with the FORMAT column and the sample
    /// names when there are samples, then the appended column names.
    VcfColumns {
        samples: Vec<Vec<u8>>,
    },
}

fn is_bed(f: Format) -> bool {
    matches!(f, Format::Bed | Format::BedGraph | Format::Wig)
}

/// D42: output that appends columns is TSV. `.tsv`/`.tsv.gz` always work;
/// a BED-family A also takes the BED extensions, as the result is BED-like.
pub fn tsv_output_ok(a_format: Format, to: OutFormat) -> Result<(), UsageError> {
    let bed_like = is_bed(a_format) && matches!(to, OutFormat::Bed | OutFormat::BedGraph);
    if to == OutFormat::Tsv || bed_like {
        return Ok(());
    }
    let works = if is_bed(a_format) {
        ".tsv, .tsv.gz, .bed, .bed.gz, .bedgraph or .bedgraph.gz"
    } else {
        ".tsv or .tsv.gz (.bed only works when the input is BED)"
    };
    Err(UsageError(format!(
        "this command appends columns to each record, so the output is TSV, not a copy of the \
input format; use {works} for -o, or pass -O tsv"
    )))
}

/// The `gloci index` columns that match A's format, for the usage hint.
fn index_hint(a: Format) -> &'static str {
    match a {
        Format::Vcf | Format::Bcf => "-s 1 -b 2 -e 2",
        Format::Gff3 | Format::Gtf => "-s 1 -b 4 -e 5",
        Format::Bed | Format::BedGraph | Format::Wig => "-s 1 -b 2 -e 3",
        // Alignments and sequence files are never TSV inputs; BED columns are the safe hint.
        Format::Bam | Format::Cram | Format::Sam | Format::Fasta | Format::Psl | Format::Rmsk => {
            "-s 1 -b 2 -e 3"
        }
    }
}

/// Write `header` and return a writer for the rows. With `write_index`
/// (D43), a BED-family A is indexed on columns 1, 2, 3 (zero-based) when the
/// output is finished; any other A is a usage error that names the
/// `gloci index` command to run on the result.
pub fn open_tsv_writer(
    mut out: Output,
    a_format: Format,
    appended: &[String],
    header: TsvHeader,
    write_index: bool,
) -> Result<Box<dyn RecordWriter>, IoError> {
    if write_index || out.wants_index() {
        if !is_bed(a_format) {
            return Err(UsageError(format!(
                "--write-index cannot index TSV output of a {} input, because its columns are not \
the BED ones; write -o FILE.tsv.gz and then run: gloci index {} FILE.tsv.gz",
                match a_format {
                    Format::Vcf | Format::Bcf => "VCF/BCF",
                    _ => "GFF/GTF",
                },
                index_hint(a_format)
            ))
            .into());
        }
        let skip_lines = match &header {
            TsvHeader::Lines(l) => l.iter().filter(|&&b| b == b'\n').count(),
            _ => 0,
        };
        out.set_tsv_index(TsvColumns {
            seq: 1,
            begin: 2,
            end: Some(3),
            zero_based: true,
            skip_lines,
            comment: b'#',
        });
    }
    let werr = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    match &header {
        TsvHeader::None => {}
        TsvHeader::Lines(lines) => out.write_all(lines).map_err(werr)?,
        TsvHeader::VcfColumns { samples } => {
            let mut line = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_vec();
            if !samples.is_empty() {
                line.extend_from_slice(b"\tFORMAT");
                for s in samples {
                    line.push(b'\t');
                    line.extend_from_slice(s);
                }
            }
            for name in appended {
                line.push(b'\t');
                line.extend_from_slice(name.as_bytes());
            }
            line.push(b'\n');
            out.write_all(&line).map_err(werr)?;
        }
    }
    Ok(Box::new(TextWriter { out }))
}
