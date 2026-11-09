//! Tabix (TBI) and CSI indexes: building one from a finished BGZF file
//! (D20), finding the index next to a data file (D26), and reading the
//! records that overlap a set of regions through it (§4.3, D9, D10).
//!
//! Index spans are a superset of htslib's, so htslib queries against gloci
//! indexes stay correct:
//! - BED `[start, end)`, with a zero-length record indexed as `[s, s+1)`;
//! - GFF/GTF `[start-1, end)`, a zero-length feature as for BED;
//! - VCF/BCF `[POS-1, max(ref end, INFO/END))`. The ref end is gloci's `ref`
//!   footprint (so SVLEN-only SVs are covered), INFO/END counts whatever the
//!   ALT alleles are, and BCF also counts `rlen`.

use crate::atomic::AtomicFile;
use crate::bcf::{self, BcfDicts, Decoder, read_bcf_header};
use crate::bed::Bed;
use crate::detect::{Compression, Format, detect_compression};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc, UsageError};
use crate::gff::Gff;
use crate::input::{BGZF_EOF, read_up_to};
use crate::record::{
    Class, Header, Parsed, ReaderOptions, RecordBatch, RecordMeta, RecordReader, TextKind,
    open_reader_from, warn_undeclared,
};
use crate::sorted::OrderGuard;
use crate::vcf::{Vcf, VcfHeader};
use gloci_core::{
    ChromId, ChromTable, Interval, RegionSet, Strand, VcfSpan, parse_coord, thousands,
};
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::MultithreadedReader;
use noodles_core::Position;
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::ReferenceSequence;
use noodles_csi::binning_index::index::header::{
    Builder as HeaderBuilder, Header as CsiHeader, ReferenceSequenceNames,
};
use noodles_csi::binning_index::index::reference_sequence::Index as BinIndexKind;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex};
use noodles_csi::binning_index::indexer::AddRecordError;
use noodles_csi::binning_index::{
    Index as BinningIndexOf, Indexer, ReferenceSequence as _, optimize_chunks,
};
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

/// Width of the smallest bin and of a linear-index window: 2^14 = 16 kb.
const MIN_SHIFT: u8 = 14;
/// TBI depth; 2^(14 + 3*5) = 2^29 is the largest end a TBI index can hold.
const TBI_DEPTH: u8 = 5;
/// The deepest CSI index noodles reads and writes covers 2^(14 + 3*9).
const MAX_INDEX_END: u64 = 1 << 41;
/// Records checked for order and indexed at a time ...
const BATCH: usize = 4096;
/// ... or fewer, once their raw bytes reach this, so memory stays flat with
/// multi-megabyte VCF/BCF records (100k+ samples).
const BATCH_BYTES: usize = 16 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    Tbi,
    Csi,
}

impl IndexKind {
    pub fn extension(self) -> &'static str {
        match self {
            IndexKind::Tbi => "tbi",
            IndexKind::Csi => "csi",
        }
    }

    fn other(self) -> IndexKind {
        match self {
            IndexKind::Tbi => IndexKind::Csi,
            IndexKind::Csi => IndexKind::Tbi,
        }
    }
}

/// An index found next to a data file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexInfo {
    pub kind: IndexKind,
    pub path: PathBuf,
    /// The index is older than the data file (D26).
    pub stale: bool,
}

/// `DATA.tbi` or `DATA.csi`.
pub(crate) fn index_path(data: &Path, kind: IndexKind) -> PathBuf {
    let mut s: OsString = data.as_os_str().to_owned();
    s.push(".");
    s.push(kind.extension());
    PathBuf::from(s)
}

/// Look for `DATA.tbi` and `DATA.csi`, preferring the newer one (the CSI
/// when both are equally new, since it can describe any file).
pub fn find_index(data: &Path) -> Option<IndexInfo> {
    let data_mtime = fs::metadata(data).and_then(|m| m.modified()).ok();
    let mut best: Option<(Option<SystemTime>, IndexInfo)> = None;
    for kind in [IndexKind::Tbi, IndexKind::Csi] {
        let path = index_path(data, kind);
        let Ok(meta) = fs::metadata(&path) else {
            continue;
        };
        if !meta.is_file() {
            continue;
        }
        let mtime = meta.modified().ok();
        if best.as_ref().is_some_and(|(t, _)| *t > mtime) {
            continue;
        }
        let stale = matches!((mtime, data_mtime), (Some(i), Some(d)) if i < d);
        best = Some((mtime, IndexInfo { kind, path, stale }));
    }
    best.map(|(_, info)| info)
}

fn loc(binary: bool, n: u64) -> Loc {
    if binary { Loc::Record(n) } else { Loc::Line(n) }
}

/// A noodles position; index positions are within `1..=2^41`.
fn position(v: u64) -> Position {
    usize::try_from(v)
        .ok()
        .and_then(Position::new)
        .expect("index positions are within 1..=2^41")
}

/// Remove a trailing `\n` (and a `\r` before it) from the line at `buf[off..]`.
fn strip_eol(buf: &mut Vec<u8>, off: usize) {
    if buf.last() == Some(&b'\n') {
        buf.pop();
        if buf.len() > off && buf.last() == Some(&b'\r') {
            buf.pop();
        }
    }
}

/// INFO/END of a VCF line, if present and an integer.
fn vcf_info_end(line: &[u8]) -> Option<u64> {
    let info = line.split(|&b| b == b'\t').nth(7)?;
    info.split(|&b| b == b';')
        .find_map(|kv| kv.strip_prefix(b"END="))
        .and_then(|v| parse_coord(v).ok())
}

fn le_u32(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])
}

// ---- per-format record parsing ----

/// Record parsing shared by index building and region reading.
enum Parser {
    Bed(Bed),
    Gff(Gff),
    Vcf(Vcf),
    Bcf(Box<(Decoder, BcfDicts)>),
    /// A generic TSV indexed by user-given columns (build only).
    Tsv(TsvKind),
}

impl Parser {
    fn text(format: Format, span: VcfSpan) -> Parser {
        match format {
            // WIG is never indexed; its lines would fail to parse as BED.
            Format::Bed | Format::BedGraph | Format::Wig => Parser::Bed(Bed),
            Format::Gff3 | Format::Gtf => Parser::Gff(Gff),
            Format::Vcf => Parser::Vcf(Vcf::new(span)),
            Format::Bcf => unreachable!("BCF is not a text format"),
        }
    }

    fn bcf(span: VcfSpan, dicts: BcfDicts) -> Parser {
        Parser::Bcf(Box::new((Decoder::new(span, &dicts), dicts)))
    }

    fn classify(&self, line: &[u8], seen_record: bool) -> Class {
        match self {
            Parser::Bed(_) => <Bed as TextKind>::classify(line, seen_record),
            Parser::Gff(_) => <Gff as TextKind>::classify(line, seen_record),
            Parser::Vcf(_) => <Vcf as TextKind>::classify(line, seen_record),
            Parser::Bcf(_) => Class::Record,
            Parser::Tsv(k) => k.classify(line),
        }
    }

    fn parse(&mut self, raw: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        match self {
            Parser::Bed(k) => k.parse(raw, chroms),
            Parser::Gff(k) => k.parse(raw, chroms),
            Parser::Vcf(k) => k.parse(raw, chroms),
            Parser::Bcf(b) => {
                let (dec, dicts) = &mut **b;
                let (chrom, key, span) = dec.decode(raw, dicts, chroms)?;
                Ok((chrom, key, span, Strand::Unknown))
            }
            Parser::Tsv(k) => k.parse(raw, chroms),
        }
    }

    /// The record's index span `(chrom, begin, end)` (D20). VCF/BCF
    /// parsers must have been made with the `ref` span.
    fn index_span(
        &mut self,
        raw: &[u8],
        chroms: &mut ChromTable,
    ) -> Result<(ChromId, u64, u64), String> {
        let (chrom, key, span, _) = self.parse(raw, chroms)?;
        let end = match self {
            Parser::Bed(_) | Parser::Gff(_) if span.is_point() => span.start.saturating_add(1),
            Parser::Bed(_) | Parser::Gff(_) => span.end,
            // parse() already widened an empty TSV interval to one base
            Parser::Tsv(_) => span.end,
            Parser::Vcf(_) => span.end.max(vcf_info_end(raw).unwrap_or(0)),
            Parser::Bcf(b) => {
                // shared section: CHROM, POS, rlen, ...
                let rlen = le_u32(raw, 16) as i32;
                let by_rlen = if rlen > 0 { key + rlen as u64 } else { 0 };
                let info_end = b.0.info_end(raw)?.unwrap_or(0);
                span.end.max(info_end).max(by_rlen)
            }
        };
        Ok((chrom, key, end))
    }

    /// Is the raw record on reference `rid` (BCF) or chromosome `name`?
    /// A record too short to tell is let through, so parsing reports it.
    fn on_reference(&self, raw: &[u8], rid: usize, name: &[u8]) -> bool {
        match self {
            Parser::Bcf(_) => {
                raw.len() < 12
                    || le_u32(raw, 0) < 24
                    || usize::try_from(le_u32(raw, 8) as i32).ok() == Some(rid)
            }
            _ => {
                let end = memchr::memchr(b'\t', raw).unwrap_or(raw.len());
                &raw[..end] == name
            }
        }
    }
}

// ---- index building ----

mod build;
mod query;
mod tsv;

use build::check_bgzf;
pub use build::write_index;
#[cfg(test)]
use build::{Entry, Linear, bin_start, csi_depth, scan};
pub(crate) use build::{build_index, remove_indexes, stage_index};
#[cfg(test)]
use query::Plan;
pub use query::{
    IndexedFile, file_chrom_name, index_non_empty_reference_names, index_reference_names,
    open_region_reader,
};
use tsv::TsvKind;
pub(crate) use tsv::build_tsv_index;
pub use tsv::{TsvColumns, write_tsv_index};

#[cfg(test)]
mod tests {
    use super::*;

    /// Records over 1 MB each (a 300k-sample VCF): the scan holds at most
    /// `BATCH_BYTES` plus one record, not a whole 4096-record batch, and
    /// queries still find the right records.
    #[test]
    fn scan_memory_is_bounded_by_bytes() {
        let samples = "\t0|1".repeat(300_000);
        let mut vcf = b"##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n"
            .to_vec();
        let mut longest = 0;
        let mut records = 0;
        for i in 0..24u64 {
            let line = format!(
                "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT{samples}\n",
                1000 + i * 100
            );
            longest = longest.max(line.len());
            records += line.len();
            vcf.extend_from_slice(line.as_bytes());
        }
        assert!(records > BATCH_BYTES + 8 * longest);
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(&vcf).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("big.vcf.gz");
        fs::write(&p, w.finish().unwrap()).unwrap();

        let mut n = 0;
        let mut sink = |_: &Entry| -> Result<(), IoError> {
            n += 1;
            Ok(())
        };
        let threads = NonZeroUsize::new(2).unwrap();
        let s = scan(&p, "big.vcf.gz", Format::Vcf, None, threads, &mut sink).unwrap();
        assert_eq!(n, 24);
        assert!(
            s.peak_bytes <= BATCH_BYTES + longest,
            "held {} bytes, budget {} + one record of {longest}",
            s.peak_bytes,
            BATCH_BYTES
        );

        write_index(&p, &p, Format::Vcf, threads).unwrap();
        let o = ReaderOptions {
            chrom_map: None,
            threads,
            format: None,
            vcf_span: VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
        };
        let mut chroms = ChromTable::first_appearance();
        let mut batch = RecordBatch::default();
        let mut query = |iv: Interval| -> Vec<(u64, u32)> {
            let regions = [(b"chr1".to_vec(), iv)];
            let mut r = open_region_reader(&p, &regions, &o, false).unwrap();
            let mut out = Vec::new();
            loop {
                batch.clear();
                let more = r.read_batch(&mut chroms, &mut batch, 4).unwrap();
                out.extend(batch.recs.iter().map(|m| (m.key, m.len)));
                if !more {
                    break;
                }
            }
            out
        };
        let len = |i: u64| {
            format!("chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT{samples}", 1000 + i * 100).len() as u32
        };
        assert_eq!(
            query(Interval {
                start: 2099,
                end: 2100
            }),
            [(2099, len(11))]
        );
        let all = query(Interval {
            start: 0,
            end: 1 << 20,
        });
        assert_eq!(all.len(), 24);
        assert!(
            all.iter()
                .enumerate()
                .all(|(i, &(k, l))| k == 999 + 100 * i as u64 && l == len(i as u64))
        );
    }

    #[test]
    fn bin_starts_follow_the_csi_layout() {
        assert_eq!(bin_start(0, 14, 5), 0);
        assert_eq!(bin_start(1, 14, 5), 0);
        assert_eq!(bin_start(2, 14, 5), 1 << 26);
        assert_eq!(bin_start(4681, 14, 5), 0); // first leaf
        assert_eq!(bin_start(4682, 14, 5), 1 << 14);
        assert_eq!(bin_start(585 + 3, 14, 5), 3 << 17);
    }

    #[test]
    fn csi_depth_covers_the_largest_end() {
        assert_eq!(csi_depth(1), 5);
        assert_eq!(csi_depth(1 << 29), 5);
        assert_eq!(csi_depth((1 << 29) + 1), 6);
        assert_eq!(csi_depth(600_000_010), 6);
        assert_eq!(csi_depth(1 << 40), 9);
        assert_eq!(csi_depth(1 << 41), 9);
    }

    #[test]
    fn sparse_linear_index() {
        let vp = |n: u64| VirtualPosition::from(n);
        let mut l = Linear::default();
        l.add(10, vp(1)); // window 0
        l.add(100_000, vp(2)); // windows 1..=6
        l.add(50_000, vp(3)); // already covered
        l.add(200_000, vp(4)); // windows 7..=12
        assert_eq!(l.at(0), vp(1));
        assert_eq!(l.at(1), vp(2));
        assert_eq!(l.at(6), vp(2));
        assert_eq!(l.at(7), vp(4));
        assert_eq!(l.at(99), vp(4));
    }

    #[test]
    fn plan_windows_widen_merge_and_clamp() {
        let iv = |s, e| Interval { start: s, end: e };
        let p = Plan::new(
            0,
            b"c".to_vec(),
            vec![iv(10, 20), iv(0, 0), iv(21, 30), iv(100, 100), iv(10, 20)],
            1 << 29,
        );
        assert_eq!(p.set.regions().len(), 4);
        // [0,0) -> (0,1); [10,20) and [21,30) widen into touching windows
        assert_eq!(p.windows, [(0, 1), (9, 31), (99, 101)]);
        let p = Plan::new(
            0,
            b"c".to_vec(),
            vec![iv(5, u64::MAX), iv(1 << 30, 1 << 31)],
            1 << 29,
        );
        assert_eq!(p.windows, [(4, 1 << 29)]);
        assert!(p.set.matches(iv(7, 7), false));
        assert!(!p.set.matches(iv(5, 5), false));
        assert!(p.set.matches(iv(5, 5), true));
    }
}
