//! The every-record alignment source (D120): `view`, `sort`, `cat` and P3c's per-record
//! commands read BAM and SAM through it, so unmapped reads are kept, where the interval
//! readers skip them (D101). Records come in the reader's raw layout: BAM bytes (without
//! `block_size`) for BAM and CRAM, the line (without the newline) for SAM.

use super::bam::BamFrames;
use super::sam::SamLines;
use super::{AlignOptions, Aln, AlnHeader};
use crate::detect::Format;
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::IoError;
use crate::input::{Input, InputOptions, open_input, open_input_from};
use crate::record::ReaderOptions;
use std::io::Read;
use std::sync::Arc;

/// Every record of an alignment file, in file order.
pub trait AlnSource: Send {
    /// The input's name, for messages.
    fn name(&self) -> &str;
    /// The layout of the raw records: `Format::Bam` for BAM and CRAM, `Format::Sam` for SAM.
    fn layout(&self) -> Format;
    fn header(&self) -> &Arc<AlnHeader>;
    /// Append the next record (no `block_size`) to `buf` and return its record number,
    /// counted from 1 over the file. `None` at the end. Records the read filters drop are
    /// counted in the diagnostics and not returned.
    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError>;
    fn diagnostics(&self) -> &Diagnostics;
}

struct BamSource {
    frames: BamFrames,
    diag: Diagnostics,
    policy: InvalidPolicy,
    opts: AlignOptions,
}

impl AlnSource for BamSource {
    fn name(&self) -> &str {
        &self.frames.name
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn header(&self) -> &Arc<AlnHeader> {
        &self.frames.aln
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }

    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        loop {
            let start = buf.len();
            if !self.frames.next(buf)? {
                return Ok(None);
            }
            let size = buf.len() - start;
            let problem = if size < 32 {
                format!(
                    "block_size {size} is below the 32 bytes of fixed fields of a BAM record; \
the file is corrupt"
                )
            } else {
                match Aln::parse(&buf[start..], &self.frames.aln) {
                    Ok(a) if self.opts.drops(a.flag(), a.mapq()) => {
                        self.diag.filtered += 1;
                        buf.truncate(start);
                        continue;
                    }
                    Ok(_) => return Ok(Some(self.frames.rec)),
                    Err(msg) => msg,
                }
            };
            buf.truncate(start);
            self.frames.invalid(self.policy, &mut self.diag, &problem)?;
        }
    }
}

struct SamSource {
    lines: SamLines,
    diag: Diagnostics,
    policy: InvalidPolicy,
    opts: AlignOptions,
}

impl AlnSource for SamSource {
    fn name(&self) -> &str {
        &self.lines.name
    }
    fn layout(&self) -> Format {
        Format::Sam
    }
    fn header(&self) -> &Arc<AlnHeader> {
        &self.lines.aln
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }

    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        loop {
            let start = buf.len();
            if !self.lines.next(buf, &mut self.diag)? {
                return Ok(None);
            }
            let problem = match self.lines.parse(&buf[start..]) {
                Ok(a) if self.opts.drops(a.flag(), a.mapq()) => {
                    self.diag.filtered += 1;
                    buf.truncate(start);
                    continue;
                }
                Ok(_) => return Ok(Some(self.lines.rec)),
                Err(msg) => msg,
            };
            buf.truncate(start);
            self.lines.invalid(self.policy, &mut self.diag, &problem)?;
        }
    }
}

pub(crate) fn source_from_input(
    input: Input,
    o: &ReaderOptions,
) -> Result<Box<dyn AlnSource>, IoError> {
    let opts = o.accept.alignments.clone();
    let not_alignments = |input: &Input| IoError::Detect {
        path: input.name.clone(),
        msg: format!(
            "{} is {}, but this command reads alignments (BAM, CRAM or SAM)",
            input.name,
            input.format.name()
        ),
    };
    if !input.format.is_alignment() {
        return Err(not_alignments(&input));
    }
    crate::align::gate(&input.name, input.format, &o.accept)?;
    let opts = opts.unwrap_or_default();
    let map = o.chrom_map.as_deref();
    match input.format {
        Format::Bam => Ok(Box::new(BamSource {
            frames: BamFrames::open(input, map)?,
            diag: Diagnostics::default(),
            policy: o.invalid,
            opts,
        })),
        Format::Sam => Ok(Box::new(SamSource {
            lines: SamLines::open(input, map)?,
            diag: Diagnostics::default(),
            policy: o.invalid,
            opts,
        })),
        other => Err(IoError::Detect {
            path: input.name.clone(),
            msg: format!(
                "{} is {}, and reading it is not available in this build yet",
                input.name,
                other.name()
            ),
        }),
    }
}

/// Open the BAM or SAM file `path` (`-` is stdin) as every-record source. The reader options
/// say which command asks (`accept.alignments` must be `Some`), the read filters, `--skip-invalid`
/// and the threads. CRAM arrives with its reader.
pub fn open_aln_source(path: &str, o: &ReaderOptions) -> Result<Box<dyn AlnSource>, IoError> {
    let input = open_input(
        path,
        &InputOptions {
            threads: o.threads,
            format: o.format,
        },
    )?;
    source_from_input(input, o)
}

/// As [`open_aln_source`], for bytes already in hand, named `name`.
pub fn open_aln_source_from(
    name: &str,
    src: Box<dyn Read + Send>,
    o: &ReaderOptions,
) -> Result<Box<dyn AlnSource>, IoError> {
    let input = open_input_from(
        name,
        src,
        &InputOptions {
            threads: o.threads,
            format: o.format,
        },
    )?;
    source_from_input(input, o)
}
