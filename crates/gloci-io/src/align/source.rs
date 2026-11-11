//! The every-record alignment source (D120): `view`, `sort`, `cat` and P3c's per-record
//! commands read BAM, CRAM and SAM through it, so unmapped reads are kept, where the interval
//! readers skip them (D101). Records come in the reader's raw layout: BAM bytes (without
//! `block_size`) for BAM and CRAM, the line (without the newline) for SAM.

use super::bam::BamFrames;
use super::sam::SamLines;
use super::{AlignOptions, Aln, AlnHeader, aln_ref, alignment_query_span};
use crate::detect::Format;
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, UsageError};
use crate::input::{Input, InputOptions, open_input, open_input_from};
use crate::index::{IndexedFile, find_index};
use crate::record::{Header, ReaderOptions};
use gloci_core::{Interval, RegionSet};
use std::collections::HashMap;
use std::path::Path;
use std::io::Read;
use std::sync::Arc;

/// Every record of an alignment file, in file order.
pub trait AlnSource: Send {
    /// The input's name, for messages.
    fn name(&self) -> &str;
    /// The layout of the raw records: `Format::Bam` for BAM and CRAM, `Format::Sam` for SAM.
    fn layout(&self) -> Format;
    fn header(&self) -> &Arc<AlnHeader>;
    /// The input's format: `Format::Cram` for a CRAM, whose records come as BAM bytes.
    fn format(&self) -> Format {
        self.layout()
    }
    /// Append the next record (no `block_size`) to `buf` and return its record number,
    /// counted from 1 over the file. `None` at the end. Records the read filters drop are
    /// counted in the diagnostics and not returned.
    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError>;
    fn diagnostics(&self) -> &Diagnostics;
    /// The BAM record (without `block_size`) of the record last returned, when the source has
    /// it already (SAM encodes each record to check it) and BAM can hold it. A writer that
    /// writes BAM can use it instead of encoding the line again.
    fn last_bam(&self) -> Option<&[u8]> {
        None
    }
}

pub(crate) struct BamSource {
    frames: BamFrames,
    diag: Diagnostics,
    policy: InvalidPolicy,
    opts: AlignOptions,
}

impl BamSource {
    /// The source over `frames` (a BAM's, or a CRAM's transcoded records).
    pub(crate) fn new(frames: BamFrames, policy: InvalidPolicy, opts: AlignOptions) -> BamSource {
        BamSource {
            frames,
            diag: Diagnostics::default(),
            policy,
            opts,
        }
    }
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

    fn last_bam(&self) -> Option<&[u8]> {
        self.lines.bam()
    }

    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        loop {
            let start = buf.len();
            if !self.lines.next(buf, &mut self.diag)? {
                return Ok(None);
            }
            let problem = match self.lines.check(&buf[start..]) {
                Ok(c) if self.opts.drops(c.flag, c.mapq) => {
                    self.diag.filtered += 1;
                    buf.truncate(start);
                    continue;
                }
                Ok(c) => {
                    self.lines.rewrite_normalised(buf, start, &c);
                    return Ok(Some(self.lines.rec));
                }
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
        Format::Bam => Ok(Box::new(BamSource::new(
            BamFrames::open(input, map)?,
            o.invalid,
            opts,
        ))),
        Format::Cram => Ok(Box::new(super::cram::CramSource::open(input, o)?)),
        Format::Sam if map.is_some() => Err(UsageError(format!(
            "{} is SAM, and --chrom-map cannot rename the reference names in SAM text that is \
written as it is; convert it to BAM first (gloci view {0} -O bam -o x.bam) and read that",
            input.name
        ))
        .into()),
        // SAM: only BAM, CRAM and SAM get past the check above.
        _ => Ok(Box::new(SamSource {
            lines: SamLines::open(input, map)?,
            diag: Diagnostics::default(),
            policy: o.invalid,
            opts,
        })),
    }
}

/// Open the BAM, CRAM or SAM file `path` (`-` is stdin) as every-record source. The reader
/// options say which command asks (`accept.alignments` must be `Some`), the read filters,
/// `--reference` (for CRAM), `--skip-invalid` and the threads.
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

/// The records of `path` (BAM, CRAM or SAM) whose D116 span overlaps one of `regions` under
/// §4.3, each once, in file order. `regions` name references as the header lists them (after
/// `--chrom-map`); a name the header lacks matches nothing. The read filters of `o` apply.
///
/// A BAM or CRAM with an index (`.bai`, `.csi`, `.crai`) is read through it. Anything else is
/// scanned whole, with one warning, and filtered by the same span, so the two give the same
/// records. The span of a placed unmapped read is `[pos, pos+1)`.
pub fn open_aln_source_regions(
    path: &str,
    regions: &[(Vec<u8>, Interval)],
    o: &ReaderOptions,
) -> Result<Box<dyn AlnSource>, IoError> {
    let src = open_aln_source(path, o)?;
    let format = src.format();
    let name = src.name().to_string();
    let header = src.header().clone();
    let file = Some(Path::new(path)).filter(|p| path != "-" && p.is_file());
    let index = file
        .filter(|_| matches!(format, Format::Bam | Format::Cram))
        .and_then(|p| find_index(p).map(|i| (p, i)));
    if let Some((p, info)) = index {
        info.warn_if_stale(&name);
        drop(src);
        let ix = IndexedFile::open(p, format, Header::Alignment(header), o)?;
        return ix.aln_regions_source(regions, false);
    }
    if path != "-" {
        let hint = match format {
            Format::Bam | Format::Cram => format!("; create an index with: gloci index {name}"),
            _ => String::new(),
        };
        log::warn!("scanning the whole file {name} to find the regions{hint}");
    }
    let mut by_name: HashMap<&[u8], Vec<Interval>> = HashMap::new();
    for (chrom, iv) in regions {
        by_name.entry(&chrom[..]).or_default().push(*iv);
    }
    let by_ref = header
        .refs
        .iter()
        .map(|(n, _)| by_name.remove(&n[..]).map(RegionSet::new))
        .collect();
    Ok(Box::new(ScanRegions {
        inner: src,
        by_ref,
        yielded: 0,
    }))
}

/// A whole-file scan that yields the records that overlap the regions.
struct ScanRegions {
    inner: Box<dyn AlnSource>,
    /// The region set of each reference, by reference id.
    by_ref: Vec<Option<RegionSet>>,
    /// Records yielded, which numbers them as the index reader does.
    yielded: u64,
}

impl AlnSource for ScanRegions {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn layout(&self) -> Format {
        self.inner.layout()
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn header(&self) -> &Arc<AlnHeader> {
        self.inner.header()
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn last_bam(&self) -> Option<&[u8]> {
        self.inner.last_bam()
    }

    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        loop {
            let start = buf.len();
            let Some(_) = self.inner.next_record(buf)? else {
                return Ok(None);
            };
            let layout = self.inner.layout();
            let hit = aln_ref(layout, &buf[start..], self.inner.header())
                .ok()
                .and_then(|a| alignment_query_span(&a))
                .is_some_and(|(rid, iv)| {
                    self.by_ref
                        .get(rid)
                        .and_then(Option::as_ref)
                        .is_some_and(|set| set.matches(iv, false))
                });
            if hit {
                self.yielded += 1;
                return Ok(Some(self.yielded));
            }
            buf.truncate(start);
        }
    }
}
