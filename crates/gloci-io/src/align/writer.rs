//! BAM and SAM output (D104): the alignment writer, and the SAM line to BAM record encoder.
//!
//! A record's raw form is the reader's: BAM bytes (without `block_size`) for BAM and CRAM,
//! the line (without the newline) for SAM. The writer converts only when the layouts differ:
//! BAM to BAM and SAM to SAM copy the bytes.

use super::Aln;
use super::encode::sam_to_bam;
use super::header::AlnHeader;
use crate::detect::Format;
use crate::error::{IoError, UsageError};
use crate::output::{OutFormat, Output};
use crate::writer::{RecordWriter, WriterOptions};
use std::io::{self, Write};

/// Encodes SAM record lines as BAM records against one header (kept by value; the name
/// lookup is the header's own).
pub struct SamToBam {
    header: AlnHeader,
}

impl SamToBam {
    pub fn new(h: &AlnHeader) -> SamToBam {
        SamToBam { header: h.clone() }
    }

    /// Append the BAM record (without `block_size`) for the SAM line `line` to `out`, as
    /// `samtools view -b` writes it. On an error `out` is as it was.
    pub fn encode(&self, line: &[u8], out: &mut Vec<u8>) -> Result<(), String> {
        sam_line_to_bam(line, &self.header, out)
    }
}

/// Append the BAM record (without `block_size`) for the SAM line `line` to `out`, with the
/// reference ids of `h`. On an error `out` is as it was.
pub fn sam_line_to_bam(line: &[u8], h: &AlnHeader, out: &mut Vec<u8>) -> Result<(), String> {
    sam_to_bam(line, &|n| h.ref_id(n), out)
}

/// Renders BAM records as SAM lines.
struct BamToSam {
    header: AlnHeader,
    line: Vec<u8>,
}

impl BamToSam {
    /// The SAM line (with its newline) for the BAM record `raw`, in a buffer valid until the
    /// next call.
    fn line(&mut self, raw: &[u8]) -> Result<&[u8], String> {
        self.line.clear();
        Aln::parse(raw, &self.header)?.write_sam(&self.header, &mut self.line)?;
        self.line.push(b'\n');
        Ok(&self.line)
    }
}

enum Convert {
    /// The raw bytes are the output's.
    Same,
    /// SAM lines to BAM records.
    ToBam(Box<SamToBam>, Vec<u8>),
    /// BAM records to SAM lines.
    ToSam(Box<BamToSam>),
}

struct AlignmentWriter {
    out: Output,
    /// BAM output frames each record with its `block_size`; SAM ends each line.
    bam: bool,
    convert: Convert,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl RecordWriter for AlignmentWriter {
    fn write_alignment(&mut self, raw: &[u8], bam: Option<&[u8]>) -> io::Result<()> {
        match (&self.convert, bam) {
            (Convert::ToBam(..), Some(b)) => {
                self.out.write_all(&(b.len() as u32).to_le_bytes())?;
                self.out.write_all(b)
            }
            _ => self.write_record(raw),
        }
    }

    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        match &mut self.convert {
            Convert::Same if self.bam => {
                self.out.write_all(&(raw.len() as u32).to_le_bytes())?;
                self.out.write_all(raw)
            }
            Convert::Same => {
                self.out.write_all(raw)?;
                self.out.write_all(b"\n")
            }
            Convert::ToBam(enc, buf) => {
                buf.clear();
                enc.encode(raw, buf).map_err(invalid)?;
                self.out.write_all(&(buf.len() as u32).to_le_bytes())?;
                self.out.write_all(buf)
            }
            Convert::ToSam(dec) => {
                let line = dec
                    .line(raw)
                    .map_err(|e| invalid(format!("cannot write a BAM record as SAM: {e}")))?;
                self.out.write_all(line)
            }
        }
    }

    fn finish(self: Box<Self>) -> Result<(), IoError> {
        self.out.finish()
    }
}

/// Write the header (unless `o.write_header` is off) and return a writer for records in the
/// layout of `records` (`Bam` and `Cram` give BAM bytes, `Sam` gives lines) as `to` (BAM or
/// SAM). With `o.command_line` the header gets gloci's `@PG` line (D104).
pub fn open_alignment_writer(
    mut out: Output,
    records: Format,
    to: OutFormat,
    h: &AlnHeader,
    o: &WriterOptions,
) -> Result<Box<dyn RecordWriter>, IoError> {
    let bam = match to {
        OutFormat::Bam => true,
        OutFormat::Sam => false,
        other => {
            return Err(UsageError(format!(
                "{} output is not possible for alignments; write BAM or SAM (-o x.bam, -o x.sam), \
or pipe BAM through: gloci view in.bam -O bam | gloci convert - -o x.bed",
                other.name()
            ))
            .into());
        }
    };
    if bam && !o.write_header {
        return Err(UsageError(
            "--no-header cannot be used with BAM output; BAM records need their header; \
write SAM with -O sam"
                .into(),
        )
        .into());
    }
    let from_sam = match records {
        Format::Bam | Format::Cram => false,
        Format::Sam => true,
        other => {
            return Err(IoError::Detect {
                path: "output".into(),
                msg: format!(
                    "internal error: {} records are not alignments",
                    other.name()
                ),
            });
        }
    };
    let werr = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    // SAM text shows the references of a header that kept them only in its binary list, and
    // so does BAM output that gets a @PG line: samtools rebuilds the text then. Without a
    // @PG line BAM output keeps the text as stored.
    let shown;
    let base = if bam && o.command_line.is_none() {
        h
    } else {
        shown = h.with_sq_lines();
        &shown
    };
    let with_pg;
    let header = match &o.command_line {
        Some(cmd) => {
            with_pg = base.with_program(cmd, o.version);
            &with_pg
        }
        None => base,
    };
    if o.write_header {
        if bam {
            header.write_bam(&mut out).map_err(werr)?;
        } else {
            out.write_all(&header.text).map_err(werr)?;
            if header.text.last().is_some_and(|&b| b != b'\n') {
                out.write_all(b"\n").map_err(werr)?;
            }
        }
    }
    let convert = match (from_sam, bam) {
        (true, true) => Convert::ToBam(Box::new(SamToBam::new(h)), Vec::new()),
        (false, false) => Convert::ToSam(Box::new(BamToSam {
            header: h.clone(),
            line: Vec::new(),
        })),
        _ => Convert::Same,
    };
    Ok(Box::new(AlignmentWriter { out, bam, convert }))
}
