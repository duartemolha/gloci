//! `gloci sort` of alignments (D121): BAM, SAM and CRAM records in samtools 1.21's coordinate,
//! `-n` and `-N` orders, through the external sorter. Records travel as BAM bytes (a SAM line is
//! encoded once, when it is read), so the spill runs and the output writer see one layout.

use crate::sort::{Order, SortConfig, SortStats, Sorter, coordinate_key, name_cmp_raw};
use gloci_io::{
    AlnHeader, AlnSortOrder, AlnSource, Format, IoError, Loc, OutFormat, Output, WriterOptions,
    open_alignment_writer, sam_line_to_bam,
};
use std::cmp::Ordering;
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// How to sort alignments.
#[derive(Debug, Clone)]
pub struct AlnSortOptions {
    pub order: AlnSortOrder,
    /// Memory budget in bytes for the records held in memory.
    pub max_mem: u64,
    pub tmp_dir: PathBuf,
    /// More than 1 sorts chunks in parallel; never changes the output.
    pub threads: NonZeroUsize,
}

/// The records of a source as BAM bytes, each SAM line encoded once.
struct BamRecords {
    src: Box<dyn AlnSource>,
    header: std::sync::Arc<AlnHeader>,
    buf: Vec<u8>,
    enc: Vec<u8>,
    /// Where the current record's BAM bytes are.
    at: At,
}

enum At {
    Buf,
    Source,
    Enc,
}

impl BamRecords {
    fn new(src: Box<dyn AlnSource>) -> BamRecords {
        BamRecords {
            header: src.header().clone(),
            src,
            buf: Vec::new(),
            enc: Vec::new(),
            at: At::Buf,
        }
    }

    /// Advance to the next record and return its number in the file.
    fn next(&mut self) -> Result<Option<u64>, IoError> {
        self.buf.clear();
        let Some(n) = self.src.next_record(&mut self.buf)? else {
            return Ok(None);
        };
        if self.src.layout() != Format::Sam {
            self.at = At::Buf;
        } else if self.src.last_bam().is_some() {
            self.at = At::Source;
        } else {
            self.enc.clear();
            sam_line_to_bam(&self.buf, &self.header, &mut self.enc).map_err(|msg| {
                IoError::Parse {
                    path: self.src.name().to_string(),
                    loc: Loc::Record(n),
                    msg: format!("{msg}; correct the record or drop it with --skip-invalid"),
                }
            })?;
            self.at = At::Enc;
        }
        Ok(Some(n))
    }

    /// The current record, BAM bytes without `block_size`.
    fn bam(&self) -> &[u8] {
        match self.at {
            At::Buf => &self.buf,
            At::Source => self.src.last_bam().unwrap_or(&self.buf),
            At::Enc => &self.enc,
        }
    }

    fn finish(self) {
        let name = self.src.name();
        self.src.diagnostics().log_alignment_counts(name);
        if let Some(s) = self.src.diagnostics().summary(name) {
            log::warn!("{s}");
        }
    }
}

/// Sort every record of `src` and write them to `out` as `to` (BAM or SAM), with the header's
/// `@HD` set for the order (D121) and gloci's `@PG` line when `w.command_line` is set.
pub fn sort_alignments(
    src: Box<dyn AlnSource>,
    o: &AlnSortOptions,
    out: Output,
    to: OutFormat,
    w: &WriterOptions,
) -> Result<SortStats, IoError> {
    // The spill directory comes first, so a bad --tmp-dir fails before any input is read.
    let mut sorter = Sorter::new(SortConfig {
        max_mem: o.max_mem,
        tmp_dir: o.tmp_dir.clone(),
        threads: o.threads,
        order: Order::Aln(o.order),
        dedup: crate::sort::Dedup::All,
    })?;
    let mut recs = BamRecords::new(src);
    let header = recs.header.with_sort_order(o.order);
    while recs.next()?.is_some() {
        sorter.push_raw(recs.bam())?;
    }
    recs.finish();
    let mut writer = open_alignment_writer(out, Format::Bam, to, &header, w)?;
    let stats = sorter.finish_raw(&mut |raw| writer.write_record(raw))?;
    writer.finish()?;
    Ok(stats)
}

fn order_name(order: AlnSortOrder) -> &'static str {
    match order {
        AlnSortOrder::Coordinate => "coordinate",
        AlnSortOrder::Name => "name (-n)",
        AlnSortOrder::NameBytes => "name (-N)",
    }
}

fn cmp(order: AlnSortOrder, a: &[u8], b: &[u8]) -> Ordering {
    match order {
        AlnSortOrder::Coordinate => coordinate_key(a, u32::MAX).cmp(&coordinate_key(b, u32::MAX)),
        AlnSortOrder::Name => name_cmp_raw(a, b, true),
        AlnSortOrder::NameBytes => name_cmp_raw(a, b, false),
    }
}

/// Whether the records of `src` are in `order`. At the first record out of order it logs an
/// error naming both records and returns `false`.
pub fn check_alignments(src: Box<dyn AlnSource>, order: AlnSortOrder) -> Result<bool, IoError> {
    let mut recs = BamRecords::new(src);
    let mut prev: Vec<u8> = Vec::new();
    let mut prev_n = 0;
    let mut sorted = true;
    while let Some(n) = recs.next()? {
        if prev_n != 0 && cmp(order, &prev, recs.bam()) == Ordering::Greater {
            let name = recs.src.name().to_string();
            let label =
                |raw: &[u8]| String::from_utf8_lossy(crate::sort::qname_of(raw)).into_owned();
            log::error!(
                "{name} is not in {} order: record {n} ({}) comes after record {prev_n} ({}); \
sort it with: gloci sort {name} -o sorted.bam{}",
                order_name(order),
                label(recs.bam()),
                label(&prev),
                match order {
                    AlnSortOrder::Coordinate => "",
                    AlnSortOrder::Name => " -n",
                    AlnSortOrder::NameBytes => " -N",
                }
            );
            sorted = false;
            break;
        }
        prev.clear();
        prev.extend_from_slice(recs.bam());
        prev_n = n;
    }
    recs.finish();
    Ok(sorted)
}
