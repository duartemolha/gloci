//! `gloci tag` (D168, bedtools tag): every record of a BAM, CRAM or SAM file, in input order,
//! with one aux tag added to each mapped read that overlaps features of the `-b` files.
//!
//! - **Features.** Every `-b` file (BED, bedGraph, GFF or GTF) is read into one
//!   [`IntervalIndex`], the in-memory index `--unsorted` and the `pair-*` commands hold B in
//!   (D61, D154). Each feature keeps its file, its coordinates as bedtools reads them (a
//!   zero-length BED record is `[start - 1, end + 1)`), its strand, its place in bedtools' bins
//!   and the text the tag takes from it. Memory grows with the features alone: about 140 bytes
//!   each while the index is built, plus that text.
//! - **Reads.** A record is tagged only when it is mapped (FLAGS bit 0x4 clear). Its span runs
//!   from POS over the reference bases of its CIGAR (`M`, `D`, `N`, `=`, `X`), as BamTools'
//!   `GetEndPosition` counts them: clips and padding add nothing, introns and deletions count,
//!   and a read with no such base spans nothing and is never tagged. Every other record is
//!   written as it is.
//! - **Hits.** A feature hits a read when the bases they share, divided by the read's span, reach
//!   `fraction`; the division and the comparison are in 32-bit floating point, as bedtools
//!   makes them. Under [`TagStrand::Same`] the feature's strand must be the read's, under
//!   [`TagStrand::Opposite`] the other one; a feature of unknown strand (`.`, no strand column, any
//!   text but `+` or `-`) matches neither, as in `intersect` (spec section 4.3), and its file gets
//!   one warning. bedtools tag counts such a feature as on the opposite strand
//!   (`#tag-unknown-strand`).
//! - **Order.** bedtools finds the features of one file by walking its bins (16 kb bins, then
//!   bins 8 times larger, 7 levels), so a file's hits come in bin order: finer levels first, then
//!   by bin, then in the file's order. [`bin_key`] is that place.
//! - **Value.** For each file with a hit, in `-b` order: its label, or the text of each hit joined
//!   by `,`; the files' parts are joined by `;`. An empty part (a name the file lacks) still
//!   counts, so the value can be empty or end in `;`, as bedtools writes it. The tag is a `Z`
//!   string added after the record's other tags, even when the record already has one of the
//!   same name (bedtools adds a second one).

use crate::skipped::log_reader;
use gloci_core::{ChromId, ChromTable, Interval, Strand, thousands};
use gloci_engine::{IndexEntry, IntervalIndex};
use gloci_io::{
    AlnSource, Format, IoError, Loc, OutFormat, Output, RecordBatch, RecordReader, UsageError,
    WriterOptions, aln_view, field, next_batch_size, open_alignment_writer,
};

/// What the tag holds for each file with a hit (bedtools `-labels`, `-names`, `-scores`,
/// `-intervals`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagValue {
    /// The file's label, one per `-b` file.
    Labels(Vec<String>),
    /// The name of each feature hit: BED column 4, GFF column 3.
    Names,
    /// The score of each feature hit, as written: BED column 5, GFF column 6.
    Scores,
    /// `LABEL:CHROM:START-END,NAME,SCORE,STRAND` for each feature hit, with the file's label.
    Intervals(Vec<String>),
}

/// The strand rule (bedtools `-s`, `-S`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TagStrand {
    #[default]
    Any,
    /// The feature's strand is the read's; a feature of unknown strand never matches.
    Same,
    /// The feature's strand is the other one; a feature of unknown strand never matches.
    Opposite,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TagOptions {
    /// The aux tag's name, two characters.
    pub tag: [u8; 2],
    pub value: TagValue,
    pub strand: TagStrand,
    /// The least share of the read's span a feature must cover (bedtools `-f`, default 1e-9).
    pub fraction: f32,
}

/// bedtools' default `-f`: any shared base.
pub const DEFAULT_FRACTION: f32 = 1e-9;

/// What a run did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TagSummary {
    /// Records written.
    pub records: u64,
    /// Records given the tag.
    pub tagged: u64,
}

/// A two-character tag name as the SAM specification allows it: a letter, then a letter or a
/// digit.
pub fn parse_tag_name(s: &str) -> Result<[u8; 2], UsageError> {
    match s.as_bytes() {
        [a, b] if a.is_ascii_alphabetic() && b.is_ascii_alphanumeric() => Ok([*a, *b]),
        _ => Err(UsageError(format!(
            "--tag {s:?} is not a tag name: give two characters, a letter and then a letter or a \
digit, such as YB or XT"
        ))),
    }
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Check the options against the number of `-b` files and their formats, before anything is
/// read: one label per file where labels are used, and features in BED or GFF.
pub fn validate(o: &TagOptions, files: &[(&str, Format)]) -> Result<(), IoError> {
    if let TagValue::Labels(l) | TagValue::Intervals(l) = &o.value
        && l.len() != files.len()
    {
        return Err(usage(format!(
            "--labels gives {} label(s) but there are {} -b file(s); give one label per -b, in \
the same order",
            l.len(),
            files.len()
        )));
    }
    if files.len() > usize::from(u16::MAX) {
        return Err(usage(format!(
            "{} -b files are more than tag takes ({}); join some of them first",
            files.len(),
            u16::MAX
        )));
    }
    for (name, f) in files {
        if !matches!(
            f,
            Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf
        ) {
            return Err(usage(format!(
                "{name} is {}, and tag reads its features from BED or GFF files; convert it to BED \
first (gloci convert {name} -o features.bed) and pass that",
                f.name()
            )));
        }
    }
    Ok(())
}

/// bedtools' bin of a feature `[start, end)`, as a sort key: the finest level whose one bin
/// holds the feature (level 0 has 16 kb bins, each next level bins 8 times larger, 7 levels), in
/// the top 8 bits, and the bin's number within its level below. A feature no level holds (one
/// that starts before 0) is in bin 0 of the last level, where bedtools' getBin puts it.
pub fn bin_key(start: i64, end: i64) -> u64 {
    let (mut s, mut e) = (start >> 14, (end - 1) >> 14);
    for level in 0..7u64 {
        if s == e {
            return (level << 56) | (s as u64 & ((1 << 56) - 1));
        }
        s >>= 3;
        e >>= 3;
    }
    6 << 56
}

/// One feature: bedtools' coordinates, its bin, its file and strand, and where its text is.
#[derive(Debug, Clone, Copy)]
struct Feature {
    start: i64,
    end: i64,
    key: u64,
    text: usize,
    text_len: u32,
    file: u16,
    strand: Strand,
}

/// The features of the `-b` files, held in memory.
pub struct Features {
    chroms: ChromTable,
    index: IntervalIndex,
    feats: Vec<Feature>,
    text: Vec<u8>,
}

impl Features {
    /// The number of features held.
    pub fn len(&self) -> usize {
        self.feats.len()
    }

    pub fn is_empty(&self) -> bool {
        self.feats.is_empty()
    }
}

/// Append the text of a feature for `value` to `out`: its name, its score, or its whole interval
/// (`label:chrom:start-end,name,score,strand`); nothing for labels.
fn feature_text(
    value: &TagValue,
    file: usize,
    gff: bool,
    raw: &[u8],
    chrom: &[u8],
    (start, end): (i64, i64),
    out: &mut Vec<u8>,
) {
    // bedtools reads a BED's name, score and strand from columns 4 to 6 and a GFF's from columns
    // 3, 6 and 7; a column the line lacks is empty.
    let col =
        |bed: usize, gff_col: usize| field(raw, if gff { gff_col } else { bed }).unwrap_or(b"");
    let (name, score, strand) = (col(3, 2), col(4, 5), col(5, 6));
    match value {
        TagValue::Labels(_) => {}
        TagValue::Names => out.extend_from_slice(name),
        TagValue::Scores => out.extend_from_slice(score),
        TagValue::Intervals(labels) => {
            out.extend_from_slice(labels[file].as_bytes());
            out.push(b':');
            out.extend_from_slice(chrom);
            out.extend_from_slice(format!(":{start}-{end},").as_bytes());
            out.extend_from_slice(name);
            out.push(b',');
            out.extend_from_slice(score);
            out.push(b',');
            out.extend_from_slice(strand);
        }
    }
}

/// Read every feature of `readers` (the `-b` files, in order, already checked by [`validate`])
/// into memory.
pub fn load_features(
    readers: Vec<Box<dyn RecordReader>>,
    o: &TagOptions,
) -> Result<Features, IoError> {
    let mut chroms = ChromTable::first_appearance();
    let mut feats: Vec<Feature> = Vec::new();
    let mut entries: Vec<IndexEntry> = Vec::new();
    let mut text: Vec<u8> = Vec::new();
    for (file, mut r) in readers.into_iter().enumerate() {
        let gff = matches!(r.format(), Format::Gff3 | Format::Gtf);
        let mut unknown_strand = false;
        let (mut batch, mut max) = (RecordBatch::default(), 1);
        loop {
            batch.clear();
            let more = r.read_batch(&mut chroms, &mut batch, max)?;
            max = next_batch_size(&batch, max);
            for (i, m) in batch.recs.iter().enumerate() {
                let id = u32::try_from(feats.len())
                    .ok()
                    .filter(|&id| id < u32::MAX)
                    .ok_or_else(|| {
                        usage(format!(
                            "the -b files hold more than {} features, more than tag holds; split \
them, or tag in several passes",
                            thousands(u64::from(u32::MAX) - 1)
                        ))
                    })?;
                let (mut start, mut end) = (m.span.start as i64, m.span.end as i64);
                if !gff && start == end {
                    // bedtools widens a zero-length BED record by a base on each side.
                    start -= 1;
                    end += 1;
                }
                let at = text.len();
                feature_text(
                    &o.value,
                    file,
                    gff,
                    batch.raw(i),
                    chroms.name(m.chrom),
                    (start, end),
                    &mut text,
                );
                feats.push(Feature {
                    start,
                    end,
                    key: bin_key(start, end),
                    text: at,
                    text_len: (text.len() - at) as u32,
                    file: file as u16,
                    strand: m.strand,
                });
                unknown_strand |= m.strand == Strand::Unknown;
                entries.push(IndexEntry {
                    chrom: m.chrom,
                    span: Interval {
                        start: start.max(0) as u64,
                        end: end.max(0) as u64,
                    },
                    id,
                });
            }
            if !more {
                break;
            }
        }
        log_reader(r.as_ref());
        if unknown_strand && o.strand != TagStrand::Any {
            log::warn!(
                "{} has records with unknown strand; they never match under -s or -S",
                r.name()
            );
        }
    }
    text.shrink_to_fit();
    feats.shrink_to_fit();
    Ok(Features {
        chroms,
        index: IntervalIndex::build(entries),
        feats,
        text,
    })
}

/// Scratch space for [`Features::value`].
#[derive(Default)]
struct Scratch {
    cand: Vec<u32>,
    hits: Vec<u32>,
}

impl Features {
    /// Write to `out` the tag value of a mapped read on `chrom` spanning `[start, end)` on
    /// `strand`, and say whether it gets one.
    fn value(
        &self,
        chrom: ChromId,
        (start, end): (u64, u64),
        strand: Strand,
        o: &TagOptions,
        s: &mut Scratch,
        out: &mut Vec<u8>,
    ) -> bool {
        out.clear();
        if end <= start {
            // bedtools divides the shared bases by a length of 0, which never reaches -f.
            return false;
        }
        self.index.query(
            chrom,
            Interval { start, end },
            None,
            Strand::Unknown,
            false,
            &mut s.cand,
        );
        let len = (end - start) as f32;
        let (start, end) = (start as i64, end as i64);
        s.hits.clear();
        s.hits.extend(s.cand.iter().copied().filter(|&id| {
            let f = &self.feats[id as usize];
            let shared = end.min(f.end) - start.max(f.start);
            let strand_ok = match o.strand {
                TagStrand::Any => true,
                TagStrand::Same => f.strand == strand,
                TagStrand::Opposite => f.strand != Strand::Unknown && f.strand != strand,
            };
            strand_ok && shared as f32 / len >= o.fraction
        }));
        // Ids ascend with the files, so each file's hits are one run.
        let hits = &mut s.hits[..];
        let mut i = 0;
        while i < hits.len() {
            let file = self.feats[hits[i] as usize].file;
            let n = hits[i..]
                .iter()
                .take_while(|&&id| self.feats[id as usize].file == file)
                .count();
            let run = &mut hits[i..i + n];
            i += n;
            match &o.value {
                TagValue::Labels(labels) => out.extend_from_slice(labels[file as usize].as_bytes()),
                _ => {
                    run.sort_unstable_by_key(|&id| (self.feats[id as usize].key, id));
                    for (k, &id) in run.iter().enumerate() {
                        if k > 0 {
                            out.push(b',');
                        }
                        let f = &self.feats[id as usize];
                        out.extend_from_slice(&self.text[f.text..f.text + f.text_len as usize]);
                    }
                }
            }
            out.push(b';');
        }
        // bedtools drops the last `;`.
        out.pop().is_some()
    }
}

/// Write every record of `src`, in input order, to `out` as `to` (BAM or SAM), each mapped read
/// that overlaps `features` with the tag added (see the module docs). The header is `src`'s, with
/// gloci's `@PG` line when `w` asks for one.
pub fn tag_alignments(
    mut src: Box<dyn AlnSource>,
    features: &Features,
    o: &TagOptions,
    mut out: Output,
    to: OutFormat,
    w: &WriterOptions,
) -> Result<TagSummary, IoError> {
    let header = src.header().clone();
    let layout = src.layout();
    let name = src.name().to_string();
    // An order error found while indexing names the input.
    out.set_index_input(&name);
    let mut writer = open_alignment_writer(out, layout, to, &header, w)?;
    // The features' chromosome of each reference, by reference id.
    let by_ref: Vec<Option<ChromId>> = header
        .refs
        .iter()
        .map(|(n, _)| features.chroms.get(n))
        .collect();
    let mut summary = TagSummary::default();
    let (mut buf, mut value, mut line, mut bam) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut scratch = Scratch::default();
    let write_err = |e: std::io::Error| {
        if e.kind() == std::io::ErrorKind::InvalidData {
            // The record cannot be written in this format; not a write failure.
            IoError::Detect {
                path: name.clone(),
                msg: format!("{name}: {e}"),
            }
        } else {
            IoError::Write {
                path: "output".to_string(),
                source: e,
            }
        }
    };
    while let Some(rec) = src.next_record(&mut buf)? {
        summary.records += 1;
        let r = aln_view(layout, &buf, &header).map_err(|msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(rec),
            msg,
        })?;
        let mut tagged = false;
        if r.flag() & 0x4 == 0
            && let (Some(rid), Some(pos)) = (r.ref_id(), r.pos())
            && let Some(Some(chrom)) = by_ref.get(rid)
        {
            let mut ref_len = 0u64;
            r.for_each_op(&mut |op, n| {
                if matches!(op, b'M' | b'D' | b'N' | b'=' | b'X') {
                    ref_len += u64::from(n);
                }
            });
            let span = (pos, pos + ref_len);
            tagged = features.value(*chrom, span, r.strand(), o, &mut scratch, &mut value);
        }
        if !tagged {
            writer
                .write_alignment(&buf, src.last_bam())
                .map_err(write_err)?;
        } else {
            summary.tagged += 1;
            line.clear();
            line.extend_from_slice(&buf);
            if layout == Format::Sam {
                line.push(b'\t');
                line.extend_from_slice(&o.tag);
                line.extend_from_slice(b":Z:");
                line.extend_from_slice(&value);
                // The BAM bytes a SAM source made of the line, with the tag, spare the writer
                // encoding the line again.
                let with_bam = src.last_bam().map(|b| {
                    bam.clear();
                    bam.extend_from_slice(b);
                    push_bam_tag(&mut bam, o.tag, &value);
                    &bam[..]
                });
                writer.write_alignment(&line, with_bam).map_err(write_err)?;
            } else {
                push_bam_tag(&mut line, o.tag, &value);
                writer.write_record(&line).map_err(write_err)?;
            }
        }
        buf.clear();
    }
    writer.finish()?;
    src.diagnostics().log_alignment_counts(&name);
    if let Some(s) = src.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    log::info!(
        "tagged {} of {} records of {name}",
        thousands(summary.tagged),
        thousands(summary.records)
    );
    Ok(summary)
}

/// Append a `Z` aux field to the BAM record `rec`.
fn push_bam_tag(rec: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) {
    rec.extend_from_slice(&tag);
    rec.push(b'Z');
    rec.extend_from_slice(value);
    rec.push(0);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_follow_bedtools_levels() {
        // In one 16 kb bin: level 0, by bin number.
        assert_eq!(bin_key(100, 200), 0);
        assert_eq!(bin_key(16384, 16390), 1);
        // Across a 16 kb boundary but in one 128 kb bin: level 1.
        assert_eq!(bin_key(16380, 16390), 1 << 56);
        assert_eq!(bin_key(15000, 17000), 1 << 56);
        // A zero-length record at 0, widened to [-1, 1): no level holds it.
        assert_eq!(bin_key(-1, 1), 6 << 56);
        // Level 0 always sorts before level 1, whatever the bin numbers.
        assert!(bin_key(16_000_000, 16_000_010) < bin_key(16380, 16390));
    }

    #[test]
    fn tag_names_follow_the_sam_specification() {
        assert_eq!(parse_tag_name("YB").unwrap(), *b"YB");
        assert_eq!(parse_tag_name("X1").unwrap(), *b"X1");
        for bad in ["X", "XYZ", "1X", "X-", "", "é1"] {
            assert!(parse_tag_name(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn labels_must_match_the_files() {
        let o = TagOptions {
            tag: *b"YB",
            value: TagValue::Labels(vec!["a".into()]),
            strand: TagStrand::Any,
            fraction: DEFAULT_FRACTION,
        };
        let bed = [("a.bed", Format::Bed), ("b.bed", Format::Bed)];
        assert!(validate(&o, &bed).is_err());
        assert!(validate(&o, &bed[..1]).is_ok());
        assert!(validate(&o, &[("x.vcf", Format::Vcf)]).is_err());
        assert!(validate(&o, &[("x.gff", Format::Gff3)]).is_ok());
        let names = TagOptions {
            value: TagValue::Names,
            ..o
        };
        assert!(validate(&names, &bed).is_ok());
    }
}
