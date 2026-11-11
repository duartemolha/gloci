//! RepeatMasker `.out` reader (D114): whitespace-separated columns, the
//! query (`sequence`, `begin`, `end`) is the footprint, 1-based closed.
//! `convert` is the only command that reads it.

use crate::chrom_map::ChromMap;
use crate::diag::InvalidPolicy;
use crate::error::IoError;
use crate::input::Input;
use crate::record::{Class, Parsed, TextKind, TextReader, coord, intern_chrom};
use gloci_core::{ChromTable, Strand, one_based_closed};
use std::sync::Arc;

pub(crate) struct Rmsk {
    map: Option<Arc<ChromMap>>,
}
pub(crate) type RmskReader = TextReader<Rmsk>;

impl TextReader<Rmsk> {
    pub(crate) fn open_rmsk(
        input: Input,
        policy: InvalidPolicy,
        map: Option<Arc<ChromMap>>,
    ) -> Result<Self, IoError> {
        Self::open(input, policy, Rmsk { map })
    }
}

/// The whitespace-separated fields of a RepeatMasker line.
pub fn rmsk_fields(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|b| b.is_ascii_whitespace())
        .filter(|f| !f.is_empty())
}

/// A data line has at least score, three percentages, sequence, begin, end,
/// `(left)`, strand, repeat, class, three repeat positions; the ID and the
/// `*` flag are optional.
const MIN_FIELDS: usize = 14;

impl TextKind for Rmsk {
    /// The 3-line header (two titles and a blank line) comes first.
    fn classify(line: &[u8], seen_record: bool) -> Class {
        let t = line.trim_ascii();
        if t.is_empty() {
            return Class::Skip;
        }
        if !seen_record && (t.starts_with(b"SW ") || t.starts_with(b"score ")) {
            return Class::Header;
        }
        Class::Record
    }

    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let f: Vec<&[u8]> = rmsk_fields(line).collect();
        if f.len() < MIN_FIELDS {
            return Err(format!(
                "expected at least {MIN_FIELDS} whitespace-separated RepeatMasker columns, \
found {}",
                f.len()
            ));
        }
        let strand = match f[8] {
            b"+" => Strand::Plus,
            b"C" => Strand::Minus,
            other => {
                return Err(format!(
                    "strand \"{}\" is neither + nor C (complement)",
                    String::from_utf8_lossy(other)
                ));
            }
        };
        let begin = coord("begin", f[5])?;
        let end = coord("end", f[6])?;
        if begin < 1 {
            return Err("begin must be >= 1 (RepeatMasker positions are 1-based)".into());
        }
        let span = one_based_closed(begin, end)
            .map_err(|_| format!("begin {begin} is greater than end {end}"))?;
        let name = match &self.map {
            Some(m) => m.rename(f[4]),
            None => f[4],
        };
        let chrom = intern_chrom(chroms, name)?;
        Ok((chrom, span.start, span, strand))
    }
}
