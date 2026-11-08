//! Tabix-style indexing of any BGZF-compressed TSV by user-given columns
//! (`gloci index -s/-b/-e`), on the same builder as the BED/GFF/VCF/BCF
//! indexes: BGZF virtual offsets, the contiguity-and-order check, TBI or CSI
//! by the largest end, and the staged publish.

use super::*;
use crate::record::{coord, intern_chrom};
use noodles_csi::binning_index::index::header::Format as CsiFormat;
use noodles_csi::binning_index::index::header::format::CoordinateSystem;

/// Which columns of a TSV hold the interval. Column numbers are 1-based, as
/// the user gives them to `tabix`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TsvColumns {
    pub seq: usize,
    pub begin: usize,
    /// The end column; without one the interval is a single base.
    pub end: Option<usize>,
    /// The begin column is 0-based (`[begin, end)`); otherwise it is 1-based
    /// and the interval is `[begin-1, end)`.
    pub zero_based: bool,
    /// Leading lines to skip.
    pub skip_lines: usize,
    /// Lines starting with this byte are comments.
    pub comment: u8,
}

impl TsvColumns {
    /// The generic tabix header for these columns: format 0, with the UCSC
    /// zero-based flag (0x10000) when `zero_based`. Reference names are
    /// filled in by the caller.
    pub(super) fn header_builder(&self) -> HeaderBuilder {
        let system = if self.zero_based {
            CoordinateSystem::Bed
        } else {
            CoordinateSystem::Gff
        };
        HeaderBuilder::bed()
            .set_format(CsiFormat::Generic(system))
            .set_reference_sequence_name_index(self.seq - 1)
            .set_start_position_index(self.begin - 1)
            .set_end_position_index(self.end.map(|e| e - 1))
            .set_line_comment_prefix(self.comment)
            .set_line_skip_count(u32::try_from(self.skip_lines).unwrap_or(u32::MAX))
    }

    fn check(&self) -> Result<(), UsageError> {
        for (flag, n) in [
            ("sequence", Some(self.seq)),
            ("begin", Some(self.begin)),
            ("end", self.end),
        ] {
            if n == Some(0) {
                return Err(UsageError(format!(
                    "the {flag} column is numbered from 1; column 0 does not exist"
                )));
            }
        }
        Ok(())
    }
}

/// Parses the interval of one TSV line.
pub(super) struct TsvKind {
    seq: usize,
    begin: usize,
    end: Option<usize>,
    zero_based: bool,
    comment: u8,
}

impl TsvKind {
    pub(super) fn new(c: &TsvColumns) -> TsvKind {
        TsvKind {
            seq: c.seq - 1,
            begin: c.begin - 1,
            end: c.end.map(|e| e - 1),
            zero_based: c.zero_based,
            comment: c.comment,
        }
    }

    pub(super) fn classify(&self, line: &[u8]) -> Class {
        if line.is_empty() || line[0] == self.comment {
            Class::Skip
        } else {
            Class::Record
        }
    }

    /// `(chrom, begin, [begin, end))` with `end > begin` (an empty interval
    /// is indexed as one base, as for BED).
    pub(super) fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let last = [Some(self.seq), Some(self.begin), self.end]
            .into_iter()
            .flatten()
            .max()
            .expect("seq is present");
        let mut fields = [&[][..]; 3];
        let mut found = 0;
        let want = [Some(self.seq), Some(self.begin), self.end];
        for (i, f) in line.split(|&b| b == b'\t').enumerate().take(last + 1) {
            for (slot, w) in want.iter().enumerate() {
                if *w == Some(i) {
                    fields[slot] = f;
                }
            }
            found = i + 1;
        }
        if found <= last {
            return Err(format!(
                "expected at least {} tab-separated columns (the index uses column {}), found {found}",
                last + 1,
                last + 1
            ));
        }
        let b = coord(&format!("begin (column {})", self.begin + 1), fields[1])?;
        let begin = if self.zero_based {
            b
        } else {
            b.checked_sub(1).ok_or_else(|| {
                format!(
                    "begin (column {}) is 0 but 1-based coordinates start at 1; \
add --zero-based if the file is 0-based",
                    self.begin + 1
                )
            })?
        };
        let end = match self.end {
            // With no end column, or the begin column as its own end, a
            // 1-based position covers one base and so does a 0-based one.
            None => begin + 1,
            Some(e) if e == self.begin => begin + 1,
            Some(e) => coord(&format!("end (column {})", e + 1), fields[2])?,
        };
        if end < begin {
            return Err(format!("end {end} is before begin {begin}"));
        }
        let end = end.max(begin + 1);
        let chrom = intern_chrom(chroms, fields[0])?;
        Ok((
            chrom,
            begin,
            Interval { start: begin, end },
            Strand::Unknown,
        ))
    }
}

/// Build a TBI or CSI index for the BGZF file `data` using `cols`, write it
/// atomically next to `final_data_path` and return its path (an index of the
/// other kind is removed). Messages name `final_data_path`.
///
/// Chromosomes must be contiguous and begins sorted within each one. The
/// index is a `.tbi` unless an end exceeds 2^29, then a `.csi`.
pub fn write_tsv_index(
    data: &Path,
    final_data_path: &Path,
    cols: &TsvColumns,
    threads: NonZeroUsize,
) -> Result<PathBuf, IoError> {
    cols.check()?;
    let (kind, bytes) =
        build::build_index_with(data, final_data_path, Format::Bed, Some(cols), threads)?;
    build::publish_index(final_data_path, kind, &bytes)
}

/// `write_tsv_index` split for the output path: build (and so validate) the
/// index bytes of `data` without writing anything, so the caller can publish
/// them in its own order (data renamed, then index).
pub(crate) fn build_tsv_index(
    data: &Path,
    final_data_path: &Path,
    cols: &TsvColumns,
    threads: NonZeroUsize,
) -> Result<(IndexKind, Vec<u8>), IoError> {
    cols.check()?;
    build::build_index_with(data, final_data_path, Format::Bed, Some(cols), threads)
}
