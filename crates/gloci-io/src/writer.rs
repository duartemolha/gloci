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

/// Is `to` the output family of `records`? Format conversion is not a
/// writer's job (D18).
pub fn check_family(records: Format, to: OutFormat) -> Result<(), UsageError> {
    if OutFormat::family(records) == to {
        return Ok(());
    }
    let (from, into) = (family_name(OutFormat::family(records)), family_name(to));
    Err(UsageError(format!(
        "the input is {from} but the output format is {into}; sort and view keep the format \
(convert VCF and BCF with gloci view; other conversions are not available yet)"
    )))
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
