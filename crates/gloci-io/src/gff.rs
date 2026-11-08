//! GFF3 and GTF reader (1-based closed coordinates, 9 tab-separated columns).

use crate::diag::InvalidPolicy;
use crate::error::IoError;
use crate::input::Input;
use crate::record::{Class, Parsed, TextKind, TextReader, coord, intern_chrom};
use gloci_core::{ChromTable, Strand, one_based_closed};

#[derive(Default)]
pub(crate) struct Gff;
pub(crate) type GffReader = TextReader<Gff>;

impl TextReader<Gff> {
    pub(crate) fn open_gff(input: Input, policy: InvalidPolicy) -> Result<Self, IoError> {
        Self::open(input, policy, Gff)
    }
}

impl TextKind for Gff {
    fn classify(line: &[u8], _seen: bool) -> Class {
        if line.is_empty() || line.starts_with(b"###") {
            Class::Skip
        } else if line.starts_with(b"##FASTA") {
            Class::Fatal(
                "##FASTA sections are not supported yet; remove the FASTA section first".into(),
            )
        } else if line[0] == b'#' {
            Class::Header
        } else {
            Class::Record
        }
    }

    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let mut f = line.split(|&b| b == b'\t');
        let c = f.next().unwrap_or_default();
        let (Some(s), Some(e)) = (f.nth(2), f.next()) else {
            return Err(too_few(line));
        };
        // Remaining: score, strand, phase, attributes.
        let (Some(_), Some(st), Some(_), Some(_)) = (f.next(), f.next(), f.next(), f.next()) else {
            return Err(too_few(line));
        };
        let strand = Strand::parse(st).unwrap_or(Strand::Unknown);
        let start = coord("start", s)?;
        let end = coord("end", e)?;
        let span = one_based_closed(start, end).map_err(|e| e.to_string())?;
        let chrom = intern_chrom(chroms, c)?;
        Ok((chrom, span.start, span, strand))
    }
}

fn too_few(line: &[u8]) -> String {
    let n = line.split(|&b| b == b'\t').count();
    format!("expected 9 tab-separated columns, found {n}")
}
