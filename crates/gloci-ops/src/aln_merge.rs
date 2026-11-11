//! `gloci sort -m` of alignments: merge BAM, SAM and CRAM files that are already in order, as
//! `samtools merge -c -p` does (D122).
//!
//! The keys and ties are those of `gloci sort` (D121); ties go to the earlier input, then to
//! record order. The header is built by [`merge_alignment_headers`]; the records by the k-way
//! merge in `sort::merge`, which cascades above 128 inputs.

use crate::aln_cat::sq_difference;
use crate::sort::{MAX_FAN_IN, merge_alignment_sources};
use gloci_io::{
    AlnHeader, AlnSortOrder, AlnSource, Format, IoError, OutFormat, Output, SortOrder, UsageError,
    WriterOptions, open_alignment_writer,
};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Arc;

fn field<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    line.split(|&b| b == b'\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key))
}

fn show(line: &[u8]) -> String {
    String::from_utf8_lossy(line).into_owned()
}

fn detect(path: &str, msg: String) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg,
    }
}

/// Refuse an input whose `@HD SO`/`SS` says it is in another order than `order`. An absent
/// or unknown order, `unsorted` and an absent `SS` are accepted: the records are checked as
/// they are read.
fn check_declared_order(name: &str, h: &AlnHeader, order: AlnSortOrder) -> Result<(), IoError> {
    let ss = h.sub_sort.as_deref();
    let (declared, hint) = match (order, h.sort_order) {
        (AlnSortOrder::Coordinate, SortOrder::Queryname) => (
            "SO:queryname",
            format!("pass -n/-N, or sort it: gloci sort {name} -o sorted.bam"),
        ),
        (AlnSortOrder::Name | AlnSortOrder::NameBytes, SortOrder::Coordinate) => (
            "SO:coordinate",
            format!(
                "drop -n/-N to merge by coordinate, or sort it: gloci sort {name}{} -o sorted.bam",
                if order == AlnSortOrder::Name {
                    " -n"
                } else {
                    " -N"
                }
            ),
        ),
        (AlnSortOrder::Name, SortOrder::Queryname)
            if ss.is_some_and(|s| s.starts_with(b"queryname:") && s != b"queryname:natural") =>
        {
            (
                "SS:queryname:lexicographical or another sub-order",
                format!("pass -N for byte order, or sort it: gloci sort {name} -n -o sorted.bam"),
            )
        }
        (AlnSortOrder::NameBytes, SortOrder::Queryname)
            if ss.is_some_and(|s| {
                s.starts_with(b"queryname:") && s != b"queryname:lexicographical"
            }) =>
        {
            (
                "SS:queryname:natural or another sub-order",
                format!(
                    "pass -n for natural order, or sort it: gloci sort {name} -N -o sorted.bam"
                ),
            )
        }
        _ => return Ok(()),
    };
    let wanted = match order {
        AlnSortOrder::Coordinate => "coordinate",
        AlnSortOrder::Name => "name (-n)",
        AlnSortOrder::NameBytes => "name (-N)",
    };
    Err(detect(
        name,
        format!(
            "{name}: its @HD line says {declared}, but this merge is in {wanted} order; {hint}"
        ),
    ))
}

/// The header lines of `text` that start with `tag` (`@RG`, ...), without their line ends.
fn lines<'a>(text: &'a [u8], tag: &'a [u8]) -> impl Iterator<Item = &'a [u8]> {
    text.split(|&b| b == b'\n').filter_map(move |l| {
        let l = l.strip_suffix(b"\r").unwrap_or(l);
        (l.starts_with(tag) && matches!(l.get(tag.len()), Some(b'\t') | None)).then_some(l)
    })
}

/// The header of a merge of `inputs` (`(name, header)` in order) in `order` (D122):
/// - the `@SQ` lists must be identical, else an error naming the first difference
///   (`#alignment-merge-sq`);
/// - `@RG` lines with one ID must be identical (kept once), else an error naming both lines
///   (`#alignment-merge-rg`); `@PG` lines with one ID keep the first;
/// - `@CO` lines come from every input, in order;
/// - `@HD` is the first input's, with `SO`/`SS` set for `order` as `gloci sort` does.
///
/// An input whose `@HD SO`/`SS` contradicts `order` is an error.
pub fn merge_alignment_headers(
    inputs: &[(&str, &AlnHeader)],
    order: AlnSortOrder,
) -> Result<AlnHeader, IoError> {
    let Some(&(first_name, first)) = inputs.first() else {
        return Err(UsageError("sort -m needs at least one alignment file".into()).into());
    };
    for &(name, h) in inputs {
        check_declared_order(name, h, order)?;
    }
    for &(name, h) in &inputs[1..] {
        if let Some(diff) = sq_difference(first_name, first, name, h) {
            return Err(detect(
                name,
                format!(
                    "{name}: the @SQ list differs from that of {first_name}: {diff}. Records name \
their reference by number, so merging these files would attach reads to the wrong contigs \
(samtools merge only warns). Merge files aligned to one reference, in one @SQ order \
(see docs/differences.md#alignment-merge-sq)"
                ),
            ));
        }
    }
    let first = first.with_sq_lines();
    let mut text: Vec<u8> = Vec::new();
    let mut put = |line: &[u8]| {
        text.extend_from_slice(line);
        text.push(b'\n');
    };
    // The first input's @HD, then its @SQ lines, then its lines of other kinds.
    for l in lines(&first.text, b"@HD") {
        put(l);
    }
    for l in lines(&first.text, b"@SQ") {
        put(l);
    }
    for l in first.text.split(|&b| b == b'\n') {
        let l = l.strip_suffix(b"\r").unwrap_or(l);
        if l.starts_with(b"@")
            && ![&b"@HD"[..], b"@SQ", b"@RG", b"@PG", b"@CO"]
                .iter()
                .any(|t| l.starts_with(t) && matches!(l.get(3), Some(b'\t') | None))
        {
            put(l);
        }
    }
    // @RG: one line per ID; the same ID with another line is a conflict.
    let mut rgs: Vec<(&[u8], &[u8], &str)> = Vec::new();
    for &(name, h) in inputs {
        for l in lines(&h.text, b"@RG") {
            let id = field(l, b"ID:");
            match rgs.iter().find(|(i, ..)| id.is_some() && Some(*i) == id) {
                Some(&(_, kept, kept_name)) => {
                    if kept != l {
                        return Err(detect(
                            name,
                            format!(
                                "{name}: the @RG line with ID {} differs from that in \
{kept_name}:\n  {kept_name}: {}\n  {name}: {}\nReads would be attributed to the wrong sample or \
library. Give the files' read groups distinct IDs (samtools merge -c keeps the first line \
without telling you; see docs/differences.md#alignment-merge-rg)",
                                show(id.unwrap_or_default()),
                                show(kept),
                                show(l)
                            ),
                        ));
                    }
                }
                None if id.is_none() && rgs.iter().any(|(_, k, _)| *k == l) => {}
                None => rgs.push((id.unwrap_or_default(), l, name)),
            }
        }
    }
    for (_, l, _) in &rgs {
        put(l);
    }
    // @PG: the first line of each ID.
    let mut pgs: Vec<(&[u8], &[u8])> = Vec::new();
    for &(_, h) in inputs {
        for l in lines(&h.text, b"@PG") {
            let id = field(l, b"ID:");
            let seen = match id {
                Some(id) => pgs.iter().any(|(i, _)| *i == id),
                None => pgs.iter().any(|(_, k)| *k == l),
            };
            if !seen {
                pgs.push((id.unwrap_or_default(), l));
            }
        }
    }
    for (_, l) in &pgs {
        put(l);
    }
    for &(_, h) in inputs {
        for l in lines(&h.text, b"@CO") {
            put(l);
        }
    }
    let merged = AlnHeader::parse_sam(text).map_err(IoError::Internal)?;
    Ok(merged.with_sort_order(order))
}

/// Merge `srcs`, each in `order`, and write the records to `out` as `to` (BAM or SAM) under
/// the merged header and gloci's `@PG` line (when `w.command_line` is set). More than 128
/// inputs are merged in rounds through spill files in `tmp_dir`. Returns the records written.
pub fn merge_alignments(
    srcs: Vec<Box<dyn AlnSource>>,
    order: AlnSortOrder,
    tmp_dir: &Path,
    out: Output,
    to: OutFormat,
    w: &WriterOptions,
) -> Result<u64, IoError> {
    let headers: Vec<(String, Arc<AlnHeader>)> = srcs
        .iter()
        .map(|s| (s.name().to_string(), s.header().clone()))
        .collect();
    let mut srcs: Vec<Option<Box<dyn AlnSource>>> = srcs.into_iter().map(Some).collect();
    let mut open = |i: usize| {
        srcs[i]
            .take()
            .ok_or_else(|| IoError::Internal(format!("alignment input {i} was opened twice")))
    };
    merge_alignments_lazy(&headers, &mut open, order, MAX_FAN_IN, tmp_dir, out, to, w)
}

/// As [`merge_alignments`], for inputs opened as they are needed: `headers` are the inputs'
/// names and headers, and `open(i)` opens input `i` for its merge round. At most `fan_in`
/// inputs are open at once, so memory and open files stay bounded however many there are.
#[allow(clippy::too_many_arguments)]
pub fn merge_alignments_lazy(
    headers: &[(String, Arc<AlnHeader>)],
    open: &mut dyn FnMut(usize) -> Result<Box<dyn AlnSource>, IoError>,
    order: AlnSortOrder,
    fan_in: usize,
    tmp_dir: &Path,
    out: Output,
    to: OutFormat,
    w: &WriterOptions,
) -> Result<u64, IoError> {
    let header = {
        let hs: Vec<(&str, &AlnHeader)> = headers
            .iter()
            .map(|(n, h)| (n.as_str(), h.as_ref()))
            .collect();
        merge_alignment_headers(&hs, order)?
    };
    // The records of every source reach the writer as BAM bytes.
    let mut writer = open_alignment_writer(out, Format::Bam, to, &header, w)?;
    let n = merge_alignment_sources(headers.len(), open, order, fan_in, tmp_dir, &mut |raw| {
        writer.write_record(raw)
    })?;
    writer.finish()?;
    Ok(n)
}

/// How many threads each of `open` simultaneously open readers gets, so that together they
/// stay within `threads`: at least 1 each.
pub fn threads_per_reader(threads: NonZeroUsize, open: usize) -> NonZeroUsize {
    NonZeroUsize::new(threads.get() / open.max(1)).unwrap_or(NonZeroUsize::MIN)
}
