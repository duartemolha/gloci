//! `gloci getfasta`: write the sequence under each interval (replaces
//! `bedtools getfasta`, D111).
//!
//! Intervals are read in input order and each sequence is streamed from the
//! indexed FASTA in chunks of at most [`gloci_io::FETCH_CHUNK`], so memory
//! does not depend on the interval length. A `-` strand record under `-s` is
//! read last chunk first and each chunk is reverse-complemented as it goes.
//!
//! Where bedtools is matched, the rules are its own: `--rna` only changes the
//! complement of `A` (so it has an effect under `-s` only), a VCF record is
//! its REF span with the name `REF/ALT_ID`, a GFF record is named by its
//! feature type, and under `-s` the header carries the record's own strand
//! column (empty for a BED record without one).

use crate::seq::{SeqSummary, revcomp_in_place};
use gloci_core::{ChromTable, Interval, Strand};
use gloci_io::{
    BcfText, Format, IndexedFasta, IoError, Loc, RecordBatch, RecordReader, bed12_blocks,
    column_count, field, next_batch_size,
};
use std::io::Write;

/// How the header (or the first column) names a sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastaName {
    /// `chrom:start-end`.
    Coords,
    /// `name::chrom:start-end`.
    NameCoords,
    /// `name`.
    NameOnly,
}

/// How a record and its sequence are written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqLayout {
    /// `>header` and the sequence on the next line.
    Fasta,
    /// `header TAB sequence`.
    Tab,
    /// The input record, TAB, the sequence.
    BedOut,
}

#[derive(Debug, Clone, Copy)]
pub struct GetfastaOptions {
    /// Reverse-complement `-` strand records and add `(+)`/`(-)` to the header.
    pub strand: bool,
    pub name: FastaName,
    /// BED12: join the blocks instead of reading the whole span.
    pub split: bool,
    pub layout: SeqLayout,
    /// Write `U` for the complement of `A` (only reverse-complemented
    /// sequences change).
    pub rna: bool,
    /// `--skip-invalid`: skip (and count) an interval outside the FASTA
    /// instead of failing.
    pub skip_outside: bool,
}

pub(crate) fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// The name of a record, as bedtools reads it from each format.
fn record_name(format: Format, line: &[u8], out: &mut Vec<u8>) {
    match format {
        Format::Gff3 | Format::Gtf => out.extend_from_slice(field(line, 2).unwrap_or(b"")),
        Format::Vcf | Format::Bcf => {
            let col = |n| field(line, n).unwrap_or(b".");
            out.extend_from_slice(col(3));
            out.push(b'/');
            out.extend_from_slice(col(4));
            if col(2) != b"." {
                out.push(b'_');
                out.extend_from_slice(col(2));
            }
        }
        _ => out.extend_from_slice(field(line, 3).unwrap_or(b"")),
    }
}

/// The text between the parentheses of `-s` headers.
fn strand_text(format: Format, strand: Strand, line: &[u8]) -> &'static [u8] {
    match (strand, format) {
        (Strand::Plus, _) | (_, Format::Vcf | Format::Bcf) => b"+",
        (Strand::Minus, _) => b"-",
        // A BED record without a strand column has an empty strand.
        (_, Format::Bed | Format::BedGraph) if column_count(line) < 6 => b"",
        _ => b".",
    }
}

/// Why an interval cannot be read from the FASTA, if it cannot.
pub(crate) fn outside(fasta: &IndexedFasta, chrom: &[u8], span: Interval) -> Option<String> {
    let shown = String::from_utf8_lossy(chrom);
    match fasta.len_of(chrom) {
        None => Some(format!(
            "{shown} is not in the FASTA index {}.fai",
            fasta.name()
        )),
        Some(len) if span.end > len => Some(format!(
            "{shown}:{}-{} ends past the {len} bases of {shown}",
            span.start, span.end
        )),
        Some(_) => None,
    }
}

/// Write the sequence of every interval of `input`.
pub fn getfasta(
    mut input: Box<dyn RecordReader>,
    fasta: &mut IndexedFasta,
    chroms: &mut ChromTable,
    o: &GetfastaOptions,
    out: &mut dyn Write,
) -> Result<SeqSummary, IoError> {
    let format = input.format();
    let name = input.name().to_string();
    let mut bcf = BcfText::for_reader(input.as_ref())?;
    let binary = bcf.is_some();
    let mut sum = SeqSummary::default();
    let mut first_outside: Option<String> = None;
    let (mut n_outside, mut bad_blocks) = (0u64, 0u64);
    let mut first_bad: Option<String> = None;
    let mut batch = RecordBatch::default();
    let mut head: Vec<u8> = Vec::new();
    let mut label: Vec<u8> = Vec::new();
    let mut scratch: Vec<u8> = Vec::new();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let at = |msg: String| IoError::Parse {
                path: name.clone(),
                loc: if binary {
                    Loc::Record(m.loc)
                } else {
                    Loc::Line(m.loc)
                },
                msg,
            };
            if m.span.is_point() {
                sum.points += 1;
                continue;
            }
            let chrom = chroms.name(m.chrom);
            if let Some(why) = outside(fasta, chrom, m.span) {
                first_outside.get_or_insert_with(|| why.clone());
                if !o.skip_outside {
                    // Default: a data error (exit 1); --skip-invalid skips.
                    return Err(at(format!(
                        "{why}; fix the interval or the FASTA, or pass --skip-invalid to skip \
intervals outside the FASTA"
                    )));
                }
                n_outside += 1;
                sum.skipped += 1;
                continue;
            }
            let line = match &mut bcf {
                Some(b) => b.render(batch.raw(i), m.loc)?,
                None => batch.raw(i),
            };
            let blocks = if o.split && format == Format::Bed && column_count(line) >= 12 {
                match bed12_blocks(line) {
                    Ok(b) => b.unwrap_or_else(|| vec![m.span]),
                    Err(e) => {
                        if !o.skip_outside {
                            // bedtools writes such a record's blocks unchecked.
                            return Err(at(format!(
                                "invalid BED12: {e}; fix the blocks, or pass --skip-invalid to \
skip records whose blocks are inconsistent"
                            )));
                        }
                        bad_blocks += 1;
                        first_bad.get_or_insert_with(|| format!("{name} line {}: {e}", m.loc));
                        sum.skipped += 1;
                        continue;
                    }
                }
            } else {
                vec![m.span]
            };
            // The header, or the leading columns.
            label.clear();
            let cs = String::from_utf8_lossy(chrom);
            let coords = format!("{cs}:{}-{}", m.span.start, m.span.end);
            match o.name {
                FastaName::Coords => label.extend_from_slice(coords.as_bytes()),
                FastaName::NameCoords => {
                    record_name(format, line, &mut label);
                    label.extend_from_slice(b"::");
                    label.extend_from_slice(coords.as_bytes());
                }
                FastaName::NameOnly => record_name(format, line, &mut label),
            }
            if o.strand {
                label.push(b'(');
                label.extend_from_slice(strand_text(format, m.strand, line));
                label.push(b')');
            }
            head.clear();
            match o.layout {
                SeqLayout::Fasta => {
                    head.push(b'>');
                    head.extend_from_slice(&label);
                    head.push(b'\n');
                }
                SeqLayout::Tab => {
                    head.extend_from_slice(&label);
                    head.push(b'\t');
                }
                SeqLayout::BedOut => {
                    head.extend_from_slice(line);
                    head.push(b'\t');
                }
            }
            out.write_all(&head).map_err(write_err)?;
            let minus = o.strand && m.strand == Strand::Minus;
            let mut emit = |chunk: &[u8]| -> Result<(), IoError> {
                if minus {
                    scratch.clear();
                    scratch.extend_from_slice(chunk);
                    revcomp_in_place(&mut scratch, o.rna);
                    out.write_all(&scratch)
                } else {
                    out.write_all(chunk)
                }
                .map_err(write_err)
            };
            if minus {
                for b in blocks.iter().rev() {
                    fasta.fetch_rev(chrom, *b, &mut emit)?;
                }
            } else {
                for b in &blocks {
                    fasta.fetch(chrom, *b, &mut emit)?;
                }
            }
            out.write_all(b"\n").map_err(write_err)?;
            sum.written += 1;
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    if let Some(s) = input.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    if sum.points > 0 {
        log::warn!(
            "skipped {} zero-length interval{} (they have no sequence to write)",
            sum.points,
            plural(sum.points)
        );
    }
    let mut parts = Vec::new();
    if let Some(first) = first_outside {
        parts.push(format!(
            "{n_outside} interval{} outside the FASTA (first: {first})",
            plural(n_outside)
        ));
    }
    if let Some(first) = first_bad {
        parts.push(format!(
            "{bad_blocks} BED12 record{} with inconsistent blocks (first: {first})",
            plural(bad_blocks)
        ));
    }
    match parts.len() {
        0 => {}
        1 => log::warn!("skipped {}", parts[0]),
        _ => log::warn!("skipped {} records: {}", sum.skipped, parts.join("; ")),
    }
    Ok(sum)
}

pub(crate) fn plural(n: u64) -> &'static str {
    if n == 1 { "" } else { "s" }
}
