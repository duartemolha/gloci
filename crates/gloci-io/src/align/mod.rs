//! Alignment input (BAM, CRAM, SAM): options, flag parsing and the input gate (D100), and the
//! BAM, CRAM and SAM readers and the BAM and SAM writers.

mod bam;
mod bed;
pub(crate) mod cram;
mod cram_prelude;
mod encode;
mod flags;
mod header;
mod place;
mod record;
mod region;
mod reheader;
mod sam;
mod source;
mod writer;

pub(crate) use bam::{BamReader, bam_record_place};
pub use bed::{AlnBed, AlnBedLayout, AlnBedOptions, AlnBlocks};
pub(crate) use cram::CramFile;
pub use cram::{CramReader, ReferenceCounts, read_cram_header, warn_unused_reference};
pub use cram_prelude::{
    CramPrelude, check_cram_version, cram_eof_present, cram_first_container, read_cram_prelude,
};
pub use flags::{SAM_FLAGS, flag_names, parse_flags};
pub use header::{AlnHeader, AlnSortOrder, SortOrder};
pub use record::{Aln, AlnRef, aln_ref, aln_view, blocks, blocks_within};
pub use region::alignment_query_span;
pub(crate) use region::{BamRegionParse, index_span};
pub use reheader::{
    bam_header_blocks, check_sq_compatible, check_sq_compatible_named, read_bam_file_header,
    read_bam_header_with_extent, read_new_header, write_cram_header_container,
};
pub use sam::SamAln;
pub(crate) use sam::SamReader;
pub(crate) use source::source_from_input;
pub use source::{
    AlnSource, FlagFields, FlagFieldsVisitor, RecordVisitor, open_aln_source, open_aln_source_from,
    open_aln_source_regions, regions_source_from,
};
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

/// Warn, once per call, that the read filters given (`--min-mapq`, `--require-flags`,
/// `--exclude-flags`) are not used because no input is BAM, CRAM or SAM, as
/// [`warn_unused_reference`] does for `--reference`. `opts` is `None` for a command that
/// reads no alignments.
pub fn warn_unused_read_filters(opts: Option<&AlignOptions>, any_alignment: bool) {
    let Some(o) = opts else { return };
    if any_alignment {
        return;
    }
    let given: Vec<&str> = [
        (o.min_mapq > 0, "--min-mapq"),
        (o.require_flags != 0, "--require-flags"),
        (o.exclude_flags != 0, "--exclude-flags"),
    ]
    .into_iter()
    .filter_map(|(set, name)| set.then_some(name))
    .collect();
    let (list, verb) = match given.as_slice() {
        [] => return,
        [one] => ((*one).to_string(), "is"),
        [rest @ .., last] => (format!("{} and {last}", rest.join(", ")), "are"),
    };
    log::warn!(
        "{list} {verb} not used: no input is BAM, CRAM or SAM (the read filters apply to \
alignments only)"
    );
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

fn gate_error(kind: GateKind, msg: String) -> IoError {
    IoError::Gate { kind, msg }
}

/// `name` as a command line names it: standard input is `-`.
fn as_argument(name: &str) -> &str {
    if name == "stdin" { "-" } else { name }
}

/// What follows `gloci convert NAME -o x.bed` in a message that sends `name`
/// (of `format`) to `convert`: the BED keeps the chromosome order of the
/// alignment header, or the record order of a PSL or RepeatMasker file, which
/// a command that checks the order may not expect. Empty for other formats.
pub(crate) fn converted_order_note(name: &str, format: Format) -> String {
    if format.is_alignment() {
        let genome = if name == "stdin" || name == "-" {
            "-g with the alignment file".to_string()
        } else {
            format!("-g {name}")
        };
        format!(
            " (the BED keeps the header's chromosome order: pass {genome}, or --sort, to a \
command that needs sorted input)"
        )
    } else if format.is_convert_only() {
        " (the BED keeps the file's record order: sort it with gloci sort, or pass --sort, for a \
command that needs sorted input)"
            .to_string()
    } else {
        String::new()
    }
}

/// Refuse a format the command does not read (D100).
pub fn gate(name: &str, format: Format, accept: &Accept) -> Result<(), IoError> {
    if format == Format::Bedpe {
        return Err(gate_error(
            GateKind::PairedOnly,
            format!(
                "{name} is read as BEDPE, which only pair-intersect, pair-pair and \
shuffle --format bedpe read; to treat its first end as a BED interval, pass --format bed"
            ),
        ));
    }
    if format.is_alignment() {
        // BAM, CRAM and SAM are read by every command that reads alignments (D100).
        if accept.alignments.is_some() {
            return Ok(());
        }
        return Err(gate_error(
            GateKind::NoAlignments,
            format!(
                "{name} is {}; this command does not read alignments (intersect, coverage, \
genomecov, convert, view, index, sort, cat, idxstats, flagstat, check, reheader, head, \
samples, depth, bedcov, tag and stats jaccard, fisher and reldist do); convert it to BED first: \
gloci convert {} -o x.bed{}",
                format.name(),
                as_argument(name),
                converted_order_note(name, format)
            ),
        ));
    }
    if format == Format::Fasta {
        return Err(gate_error(
            GateKind::Sequence,
            format!(
                "{name} is FASTA, which is sequence, not intervals; only getfasta, maskfasta, \
nuc, faidx or dict read it, as the reference for the intervals in a BED, GFF or VCF file, and \
-g ref.fa reads its .fai for the chromosome sizes (create it with: gloci index {name})"
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
gloci convert {} -o x.bed{}; if {name} is not {} at all, \
pass --format bed or the format it is",
                format.name(),
                as_argument(name),
                converted_order_note(name, format),
                format.name()
            ),
        ));
    }
    Ok(())
}
