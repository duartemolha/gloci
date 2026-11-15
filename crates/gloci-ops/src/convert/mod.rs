//! `convert`: VCF/BCF, GFF3/GTF, WIG, PSL, RepeatMasker and BED12 to BED (replaces bedops
//! `convert2bed`, `vcf2bed`, `gff2bed`, `gtf2bed` and bedtools
//! `bed12tobed6`).
//!
//! Records are converted in input order and the input needs no sorting. Only
//! the conversion itself can break the order of the output (an allele
//! footprint can start after its record's key, and a BED12 block after the
//! next record's start); that prints one warning (D49).

mod align;
mod align_bedpe;
mod align_fastq;
mod bed12;
mod gff;
mod psl;
mod rmsk;
mod vcf;
mod wig;

pub use align_bedpe::{BedpeOptions, BedpeSummary, alignments_to_bedpe};
pub use align_fastq::{FastqOptions, FastqSink, FastqSummary, alignments_to_fastq};

use gloci_core::{ChromId, ChromTable, Strand, VariantKind, VcfSpan};
use gloci_io::{Format, IoError, RecordReader, RecordWriter, UsageError};

#[derive(Debug, Clone, Default)]
pub struct ConvertOptions {
    /// VCF/BCF: the footprint rule. Must match the reader's span: `Ref` gives
    /// every `--split-alleles` row the record's REF footprint.
    pub vcf_span: VcfSpan,
    /// VCF/BCF: one row per ALT with that ALT's footprint.
    pub split_alleles: bool,
    /// VCF/BCF: keep records (with `split_alleles`, ALTs) of these kinds.
    pub types: Option<Vec<VariantKind>>,
    /// GFF3/GTF: keep these column-3 types.
    pub features: Vec<String>,
    /// GFF3/GTF: the attribute for the name column.
    pub name_attr: Option<String>,
    /// BED12 and alignments: one row per block.
    pub blocks: bool,
    /// Alignments: one BED12 row per read.
    pub bed12: bool,
    /// Alignments: the BED12 `itemRgb` (`R,G,B`, default `255,0,0`).
    pub color: Option<String>,
    /// Alignments: the integer tag the score column takes (default MAPQ).
    pub score_tag: Option<[u8; 2]>,
    /// Alignments: append the CIGAR (BED6 only).
    pub cigar: bool,
    /// Alignments: cut blocks at `D` as well as `N`.
    pub split_deletions: bool,
    /// Alignments: write BEDPE pairs (D156, `alignments_to_bedpe`) instead of BED.
    pub bedpe: bool,
    /// With `bedpe`: keep READ1 as the first end.
    pub mate1_first: bool,
}

/// The output rows, with the order warning.
pub(crate) struct Sink<'a> {
    w: &'a mut dyn RecordWriter,
    /// Chromosome, input key and output start of the previous row.
    last: Option<(ChromId, u64, u64)>,
    rows: u64,
    warned: bool,
}

impl Sink<'_> {
    /// Write one row. `key` is the position the input is sorted by.
    pub(crate) fn row(
        &mut self,
        chrom: ChromId,
        key: u64,
        start: u64,
        row: &[u8],
    ) -> Result<(), IoError> {
        self.rows += 1;
        // Warn only when the conversion itself breaks the order: same
        // chromosome, input in order, output start going down.
        if !self.warned
            && self
                .last
                .is_some_and(|(c, k, s)| c == chrom && k <= key && s > start)
        {
            self.warned = true;
            log::warn!(
                "output row {} is out of order; the output of this command is not sorted \
(pipe the output through gloci sort)",
                self.rows
            );
        }
        self.last = Some((chrom, key, start));
        self.w.write_record(row).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })
    }
}

/// The BED strand character of a record whose reader found `+` or `-`.
pub(crate) fn strand_byte(s: Strand) -> u8 {
    match s {
        Strand::Plus => b'+',
        Strand::Minus => b'-',
        Strand::Unknown => b'.',
    }
}

/// Append `n` as decimal digits.
pub(crate) fn push_num(row: &mut Vec<u8>, n: u64) {
    row.extend_from_slice(n.to_string().as_bytes());
}

fn what(format: Format) -> &'static str {
    match format {
        Format::Vcf | Format::Bcf => "VCF/BCF",
        Format::Gff3 | Format::Gtf => "GFF3/GTF",
        Format::Bed | Format::BedGraph => "BED",
        Format::Wig => "WIG",
        Format::Bam
        | Format::Cram
        | Format::Sam
        | Format::Fasta
        | Format::Psl
        | Format::Rmsk
        | Format::Bedpe => format.name(),
    }
}

/// What `convert` does for each input, for usage errors.
const SUMMARY: &str = "convert turns VCF/BCF, GFF3/GTF, WIG, PSL, RepeatMasker .out, BAM/CRAM/SAM or (with --blocks) BED12 into BED";

impl ConvertOptions {
    fn check(&self, format: Format) -> Result<(), UsageError> {
        let usage = |flag: &str, only: &str| {
            Err(UsageError(format!(
                "{flag} applies only to {only} input, not {}; {SUMMARY}",
                what(format)
            )))
        };
        let vcf = matches!(format, Format::Vcf | Format::Bcf);
        let gff = matches!(format, Format::Gff3 | Format::Gtf);
        if !vcf {
            if self.split_alleles {
                return usage("--split-alleles", "VCF/BCF");
            }
            if self.types.is_some() {
                return usage("--types", "VCF/BCF");
            }
        }
        if !gff {
            if !self.features.is_empty() {
                return usage("--feature", "GFF3/GTF");
            }
            if self.name_attr.is_some() {
                return usage("--name-attr", "GFF3/GTF");
            }
        }
        let aln = matches!(format, Format::Bam | Format::Cram | Format::Sam);
        if self.mate1_first && !self.bedpe {
            return Err(UsageError(format!(
                "--mate1-first orders the ends of --bedpe rows and needs --bedpe; {SUMMARY}"
            )));
        }
        if self.bedpe {
            if !aln {
                return usage("--bedpe", "BAM/CRAM/SAM");
            }
            for (set, flag) in [
                (self.bed12, "--bed12"),
                (self.blocks, "--blocks"),
                (self.cigar, "--cigar"),
                (self.color.is_some(), "--color"),
                (self.split_deletions, "--split-deletions"),
            ] {
                if set {
                    return Err(UsageError(format!(
                        "--bedpe writes one row per pair of reads and cannot be used with {flag}; \
drop one of them"
                    )));
                }
            }
        }
        if !aln {
            for (set, flag) in [
                (self.bed12, "--bed12"),
                (self.color.is_some(), "--color"),
                (self.score_tag.is_some(), "--score-tag"),
                (self.cigar, "--cigar"),
                (self.split_deletions, "--split-deletions"),
            ] {
                if set {
                    return usage(flag, "BAM/CRAM/SAM");
                }
            }
        }
        if self.blocks && format != Format::Bed && !aln {
            return Err(UsageError(format!(
                "--blocks needs BED12, BAM, CRAM or SAM input, not {}; {SUMMARY}",
                what(format)
            )));
        }
        if self.color.is_some() && !self.bed12 {
            return Err(UsageError(format!(
                "--color is the itemRgb of BED12 rows and needs --bed12; {SUMMARY}"
            )));
        }
        if self.cigar && (self.bed12 || self.blocks || self.split_deletions) {
            return Err(UsageError(
                "--cigar adds a column to BED6 rows and cannot be used with --bed12, --blocks or --split-deletions; drop one of them (bedtools -cigar -split silently drops the CIGAR)"
                    .into(),
            ));
        }
        if let Some(c) = &self.color
            && (c.split(',').count() != 3 || c.split(',').any(|v| v.parse::<u8>().is_err()))
        {
            return Err(UsageError(format!(
                "--color {c} is not R,G,B with each value from 0 to 255 (for example 255,0,0)"
            )));
        }
        if format == Format::Bed && !self.blocks {
            return Err(UsageError(format!(
                "a BED input is already BED; pass --blocks to write one row per block of \
BED12 records. {SUMMARY}"
            )));
        }
        if format == Format::BedGraph {
            return Err(UsageError(format!(
                "a bedGraph input is already BED and has no blocks; {SUMMARY}"
            )));
        }
        Ok(())
    }
}

/// Check that `o` suits `format` (a usage error otherwise), before any
/// output is created.
pub fn check_convert(format: Format, o: &ConvertOptions) -> Result<(), UsageError> {
    o.check(format)
}

/// Convert every record of `input` to BED rows on `w`.
pub fn convert(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &ConvertOptions,
    w: &mut dyn RecordWriter,
) -> Result<(), IoError> {
    let format = input.format();
    o.check(format)?;
    let mut sink = Sink {
        w,
        last: None,
        rows: 0,
        warned: false,
    };
    match format {
        Format::Bam | Format::Cram | Format::Sam => align::convert(input, chroms, o, &mut sink),
        Format::Fasta => Err(UsageError(format!(
            "convert does not read {}, which is sequence and not intervals; {SUMMARY}",
            format.name()
        ))
        .into()),
        Format::Bedpe => Err(UsageError(format!(
            "convert does not read BEDPE; it writes BEDPE from alignments with --bedpe; {SUMMARY}"
        ))
        .into()),
        Format::Psl => psl::convert(input, chroms, &mut sink),
        Format::Rmsk => rmsk::convert(input, chroms, &mut sink),
        Format::Vcf | Format::Bcf => vcf::convert(input, chroms, o, &mut sink),
        Format::Gff3 | Format::Gtf => gff::convert(input, chroms, o, &mut sink),
        Format::Bed => bed12::convert(input, chroms, &mut sink),
        Format::Wig => wig::convert(input, chroms, &mut sink),
        Format::BedGraph => unreachable!("rejected by check"),
    }
}
