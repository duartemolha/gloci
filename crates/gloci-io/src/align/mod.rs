//! Alignment input (BAM, CRAM, SAM): options, flag parsing and the input gate (D100), and the
//! BAM, CRAM and SAM readers and the BAM and SAM writers.

mod bam;
mod bed;
pub(crate) mod cram;
mod encode;
mod header;
mod place;
mod record;
mod region;
mod sam;
mod source;
mod writer;

pub(crate) use bam::BamReader;
pub use bed::{AlnBed, AlnBedLayout, AlnBedOptions, AlnBlocks};
pub use cram::{CramReader, read_cram_header, warn_unused_reference};
pub use header::{AlnHeader, SortOrder};
pub use record::{Aln, AlnRef, aln_ref, blocks, blocks_within};
pub use region::alignment_query_span;
pub(crate) use region::{BamRegionParse, index_span};
pub use sam::SamAln;
pub(crate) use sam::SamReader;
pub(crate) use source::source_from_input;
pub use source::{AlnSource, open_aln_source, open_aln_source_from};
pub use writer::{SamToBam, open_alignment_writer, sam_line_to_bam};

use crate::detect::Format;
use crate::error::{GateKind, IoError};
use std::path::PathBuf;

/// How alignments are read: the filters of the `AlignArgs` flags (D101).
#[derive(Debug, Clone, Default)]
pub struct AlignOptions {
    pub min_mapq: u8,
    pub require_flags: u16,
    pub exclude_flags: u16,
    /// Keep whole records (sequence, qualities, tags), for `view` and `sort`.
    pub full_records: bool,
    /// The FASTA reference CRAM input needs.
    pub reference: Option<PathBuf>,
    /// An aux tag a slim record keeps (`--score-tag`, D115).
    pub score_tag: Option<[u8; 2]>,
}

impl AlignOptions {
    /// Do the read filters (D102) drop a record with this flag and MAPQ?
    pub fn drops(&self, flag: u16, mapq: u8) -> bool {
        mapq < self.min_mapq
            || flag & self.require_flags != self.require_flags
            || flag & self.exclude_flags != 0
    }
}

/// What a command accepts as input beyond interval formats. The default accepts
/// none of them, so a command that does not read alignments refuses them.
#[derive(Debug, Clone, Default)]
pub struct Accept {
    /// `Some` when the command reads alignments.
    pub alignments: Option<AlignOptions>,
    /// The command is `convert`, which also reads PSL and RepeatMasker.
    pub convert_only: bool,
    /// The command will read alignments once their reader lands (intersect, coverage,
    /// genomecov, view, sort, index): until then the gate says so instead of claiming
    /// that the command never reads them.
    pub alignments_planned: bool,
    /// The command reads CRAM as well (`view`); elsewhere a CRAM is "not available in this
    /// build yet" until the command opts in.
    pub cram: bool,
}

impl Accept {
    /// For a command that will read alignments but whose reader is not built yet.
    pub fn planned() -> Accept {
        Accept {
            alignments_planned: true,
            ..Accept::default()
        }
    }
}

const FLAG_NAMES: [(&str, u16); 12] = [
    ("PAIRED", 0x1),
    ("PROPER_PAIR", 0x2),
    ("UNMAP", 0x4),
    ("MUNMAP", 0x8),
    ("REVERSE", 0x10),
    ("MREVERSE", 0x20),
    ("READ1", 0x40),
    ("READ2", 0x80),
    ("SECONDARY", 0x100),
    ("QCFAIL", 0x200),
    ("DUP", 0x400),
    ("SUPPLEMENTARY", 0x800),
];

/// Parse a SAM flag set: a decimal or `0x` hexadecimal number, or comma-separated
/// names (`PAIRED`, `UNMAP`, `DUP`, ...), as `samtools view -f/-F` take.
pub fn parse_flags(s: &str) -> Result<u16, String> {
    let t = s.trim();
    let names = || {
        FLAG_NAMES
            .iter()
            .map(|(n, _)| *n)
            .collect::<Vec<_>>()
            .join(", ")
    };
    if t.is_empty() {
        return Err(format!(
            "empty flag value; give a number such as 0x704 or names such as secondary,dup \
(names: {})",
            names()
        ));
    }
    let number = if let Some(h) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        Some(u32::from_str_radix(h, 16))
    } else if t.bytes().all(|b| b.is_ascii_digit()) {
        Some(t.parse::<u32>())
    } else {
        None
    };
    if let Some(n) = number {
        return n
            .ok()
            .and_then(|n| u16::try_from(n).ok())
            .ok_or_else(|| format!("flag value '{t}' is not a number from 0 to 65535"));
    }
    let mut out = 0u16;
    for part in t.split(',') {
        let p = part.trim();
        match FLAG_NAMES.iter().find(|(n, _)| n.eq_ignore_ascii_case(p)) {
            Some((_, v)) => out |= v,
            None => {
                return Err(format!(
                    "unknown flag name '{p}'; use a number or names from: {}",
                    names()
                ));
            }
        }
    }
    Ok(out)
}

fn gate_error(kind: GateKind, msg: String) -> IoError {
    IoError::Gate { kind, msg }
}

/// Refuse a format the command does not read (D100), or say that its reader is not built yet.
pub fn gate(name: &str, format: Format, accept: &Accept) -> Result<(), IoError> {
    let not_yet = || {
        gate_error(
            GateKind::NotYet,
            format!(
                "{name} is {}, and reading it is not available in this build yet",
                format.name()
            ),
        )
    };
    if format.is_alignment() {
        if matches!(format, Format::Bam | Format::Sam) && accept.alignments.is_some() {
            return Ok(());
        }
        if format == Format::Cram && accept.alignments.is_some() {
            // CRAM is read by the commands that opt in (view); the interval commands follow.
            return if accept.cram { Ok(()) } else { Err(not_yet()) };
        }
        if accept.alignments.is_some() {
            return Err(IoError::Detect {
                path: name.to_string(),
                msg: not_yet().to_string(),
            });
        }
        if accept.alignments_planned {
            return Err(not_yet());
        }
        return Err(gate_error(
            GateKind::NoAlignments,
            format!(
                "{name} is {}; this command does not read alignments (intersect, coverage, \
genomecov, convert, view, index, sort and cat do); convert it to BED first: \
gloci convert {name} -o x.bed",
                format.name()
            ),
        ));
    }
    if format == Format::Fasta {
        return Err(gate_error(
            GateKind::Sequence,
            format!(
                "{name} is FASTA, which is sequence, not intervals; only getfasta, maskfasta or \
nuc read it, as the reference for the intervals in a BED, GFF or VCF file, and -g ref.fa \
reads its .fai for the chromosome sizes (create it with: gloci index {name})"
            ),
        ));
    }
    if format.is_convert_only() {
        if accept.convert_only {
            // `convert` reads PSL and RepeatMasker (D114).
            return Ok(());
        }
        return Err(gate_error(
            GateKind::ConvertOnly,
            format!(
                "{name} is {}; only gloci convert reads it; convert it to BED first: \
gloci convert {name} -o x.bed (if {name} is not {} at all, \
pass --format bed or the format it is)",
                format.name(),
                format.name()
            ),
        ));
    }
    Ok(())
}
