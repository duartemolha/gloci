//! Alignment input (BAM, CRAM, SAM): options, flag parsing and the input gate (D100).
//! The readers themselves arrive in later P3a tasks.

use crate::detect::Format;
use crate::error::{IoError, UsageError};
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
}

/// What a command accepts as input beyond interval formats. The default accepts
/// none of them, so a command that does not read alignments refuses them.
#[derive(Debug, Clone, Default)]
pub struct Accept {
    /// `Some` when the command reads alignments.
    pub alignments: Option<AlignOptions>,
    /// The command is `convert`, which also reads PSL and RepeatMasker.
    pub convert_only: bool,
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

/// Refuse a format the command does not read (D100), or say that its reader is not built yet.
pub fn gate(name: &str, format: Format, accept: &Accept) -> Result<(), IoError> {
    if format.is_alignment() {
        if accept.alignments.is_some() {
            return Err(IoError::Detect {
                path: name.to_string(),
                msg: format!(
                    "{name} is {}, and reading it is not available in this build yet",
                    format.name()
                ),
            });
        }
        return Err(UsageError(format!(
            "{name} is {}; this command does not read alignments (intersect, coverage, \
genomecov, convert, view, index, sort and cat do); convert it to BED first: \
gloci convert {name} -o x.bed",
            format.name()
        ))
        .into());
    }
    if format == Format::Fasta {
        return Err(UsageError(format!(
            "{name} is FASTA, which is sequence, not intervals; only getfasta, maskfasta or nuc \
read it, as the reference for the intervals in a BED, GFF or VCF file"
        ))
        .into());
    }
    if format.is_convert_only() {
        if accept.convert_only {
            return Err(IoError::Detect {
                path: name.to_string(),
                msg: format!(
                    "{name} is {}, and reading it is not available in this build yet",
                    format.name()
                ),
            });
        }
        return Err(UsageError(format!(
            "{name} is {}; only gloci convert reads it; convert it to BED first: \
gloci convert {name} -o x.bed",
            format.name()
        ))
        .into());
    }
    Ok(())
}
