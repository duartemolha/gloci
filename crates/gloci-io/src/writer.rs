//! Record output shared by the commands that write records (`sort`, `view`).

use crate::bcf::write_bcf_header;
use crate::detect::Format;
use crate::error::{IoError, UsageError};
use crate::output::{OutFormat, Output};
use crate::record::Header;
use std::io::{self, Write};

/// Writes records whose raw bytes come from a reader of the same family.
pub trait RecordWriter {
    /// Write one record: text formats get `raw` plus a newline, BCF gets the
    /// framed bytes as they are.
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()>;
    /// Flush and publish the output (BGZF EOF block, index, atomic rename).
    fn finish(self: Box<Self>) -> Result<(), IoError>;
}

#[derive(Debug, Clone)]
pub struct WriterOptions {
    /// Write the header lines of formats whose header is optional (BED, GFF).
    /// VCF and BCF always get their header.
    pub write_header: bool,
    /// The command line recorded in the `##glociCommand` line of VCF/BCF.
    pub command_line: Option<String>,
    pub version: &'static str,
}

/// Do the records of `records` fit the output format `to`? BED and bedGraph
/// are one family: the bytes pass through unchanged, and most bedGraph
/// files have no track line, so detection calls them BED. Format conversion
/// is not a writer's job (D18).
pub fn check_family(records: Format, to: OutFormat) -> Result<(), UsageError> {
    let is_bed = |f: OutFormat| matches!(f, OutFormat::Bed | OutFormat::BedGraph);
    let from = OutFormat::family(records);
    if from == to || (is_bed(from) && is_bed(to)) {
        return Ok(());
    }
    Err(family_mismatch(records, to))
}

/// The one message for output that does not match the records' format,
/// used by `sort` and `view`.
pub fn family_mismatch(records: Format, to: OutFormat) -> UsageError {
    let from = OutFormat::family(records);
    let fits = match from {
        OutFormat::Bed | OutFormat::BedGraph => ".bed, .bed.gz, .bedgraph or .bedgraph.gz",
        OutFormat::Vcf | OutFormat::Bcf => ".vcf, .vcf.gz or .bcf",
        OutFormat::Gff3 => ".gff, .gff.gz, .gff3 or .gff3.gz",
        OutFormat::Gtf => ".gtf or .gtf.gz",
        OutFormat::Tsv => ".tsv or .tsv.gz",
    };
    UsageError(format!(
        "the input is {} but the output format is {}; sort and view keep the format, so pick an \
output extension that matches ({fits}) or pass -O with a matching format (other conversions are \
not available yet)",
        family_name(from),
        family_name(to)
    ))
}

fn family_name(f: OutFormat) -> &'static str {
    match f {
        OutFormat::Bed => "BED",
        OutFormat::BedGraph => "bedGraph",
        OutFormat::Vcf => "VCF",
        OutFormat::Bcf => "BCF",
        OutFormat::Gff3 => "GFF3",
        OutFormat::Gtf => "GTF",
        OutFormat::Tsv => "TSV",
    }
}

struct TextWriter {
    out: Output,
}

impl RecordWriter for TextWriter {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        self.out.write_all(raw)?;
        self.out.write_all(b"\n")
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        self.out.finish()
    }
}

struct BcfWriter {
    out: Output,
}

impl RecordWriter for BcfWriter {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        self.out.write_all(raw)
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        self.out.finish()
    }
}

/// Write the header and return a writer for records of format `records`.
/// `Lines` headers are written only when `o.write_header`; a VCF/BCF header
/// always is, with the `##glociCommand` line. `to` must be `records`' family.
pub fn open_record_writer(
    mut out: Output,
    records: Format,
    to: OutFormat,
    header: &Header,
    o: &WriterOptions,
) -> Result<Box<dyn RecordWriter>, IoError> {
    check_family(records, to)?;
    let werr = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    match (records, header) {
        (Format::Vcf | Format::Bcf, Header::Vcf(h)) => {
            let text = match &o.command_line {
                Some(cmd) => h.with_command_line(cmd, o.version),
                None => h.raw.clone(),
            };
            if records == Format::Bcf {
                write_bcf_header(&mut out, &text).map_err(werr)?;
            } else {
                out.write_all(&text).map_err(werr)?;
            }
        }
        (Format::Vcf | Format::Bcf, _) => {
            return Err(IoError::Detect {
                path: "output".into(),
                msg: "internal error: a VCF/BCF input has no VCF header".into(),
            });
        }
        (_, Header::Lines(lines)) if o.write_header => out.write_all(lines).map_err(werr)?,
        _ => {}
    }
    Ok(if records == Format::Bcf {
        Box::new(BcfWriter { out })
    } else {
        Box::new(TextWriter { out })
    })
}

/// A writer for text records with no header at all (`view --no-header`).
/// `to` must be `records`' family, and a text one: BCF cannot do without its
/// header.
pub fn open_headerless_writer(
    out: Output,
    records: Format,
    to: OutFormat,
) -> Result<Box<dyn RecordWriter>, IoError> {
    check_family(records, to)?;
    if records == Format::Bcf {
        return Err(UsageError(
            "--no-header cannot be used with BCF output; BCF records need their header".into(),
        )
        .into());
    }
    Ok(Box::new(TextWriter { out }))
}
