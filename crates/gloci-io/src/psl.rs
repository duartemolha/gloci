//! PSL reader (BLAT output; D114): 21 tab-separated columns, the target
//! (`tName`, `tStart`, `tEnd`) is the footprint. `convert` is the only
//! command that reads it.

use crate::chrom_map::ChromMap;
use crate::diag::InvalidPolicy;
use crate::error::IoError;
use crate::fields::split_fields;
use crate::input::Input;
use crate::record::{Class, Parsed, TextKind, TextReader, coord, intern_chrom};
use gloci_core::{ChromTable, Interval, Strand};
use std::sync::Arc;

pub(crate) struct Psl {
    map: Option<Arc<ChromMap>>,
}
pub(crate) type PslReader = TextReader<Psl>;

impl TextReader<Psl> {
    pub(crate) fn open_psl(
        input: Input,
        policy: InvalidPolicy,
        map: Option<Arc<ChromMap>>,
    ) -> Result<Self, IoError> {
        Self::open(input, policy, Psl { map })
    }
}

/// The columns of a PSL row.
const COLUMNS: usize = 21;
const STRAND: usize = 8;
const T_NAME: usize = 13;
const T_START: usize = 15;
const T_END: usize = 16;

impl TextKind for Psl {
    /// A BLAT file starts with a 5-line header (`psLayout version 3`, a blank
    /// line, two column-title lines and a dashed line); a header-less file is
    /// as common.
    fn classify(line: &[u8], seen_record: bool) -> Class {
        if line.is_empty() {
            return Class::Skip;
        }
        if !seen_record
            && (line.starts_with(b"psLayout")
                || line.starts_with(b"match")
                || line.starts_with(b"---"))
        {
            return Class::Header;
        }
        Class::Record
    }

    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let f: Vec<&[u8]> = split_fields(line).collect();
        if f.len() < COLUMNS {
            return Err(format!(
                "expected {COLUMNS} tab-separated PSL columns, found {}; \
is this a PSL file (fields must be separated by tabs)?",
                f.len()
            ));
        }
        // A translated search writes two strand characters: query then target.
        let strand = match f[STRAND] {
            b"+" => Strand::Plus,
            b"-" => Strand::Minus,
            b"++" | b"-+" => Strand::Plus,
            b"+-" | b"--" => Strand::Minus,
            other => {
                return Err(format!(
                    "strand \"{}\" is not +, -, or a two-character pair such as +-",
                    String::from_utf8_lossy(other)
                ));
            }
        };
        let start = coord("tStart", f[T_START])?;
        let end = coord("tEnd", f[T_END])?;
        let span = Interval::new(start, end).map_err(|e| format!("tStart/tEnd: {e}"))?;
        let name = match &self.map {
            Some(m) => m.rename(f[T_NAME]),
            None => f[T_NAME],
        };
        let chrom = intern_chrom(chroms, name)?;
        Ok((chrom, start, span, strand))
    }
}
