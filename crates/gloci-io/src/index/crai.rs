//! CRAM indexes (`.crai`, D107, D109), built by gloci itself: noodles' `fs::index` decodes
//! slices of several references with an empty repository, and can panic there.
//!
//! A `.crai` is gzip-compressed text, one line per slice and reference: reference id, start,
//! span, container offset, slice offset (in the container's data) and slice size. A slice on
//! one reference gives the line its header describes; a slice of several references is decoded
//! (positions only: no reference bases are read, so indexing never needs `--reference`, where
//! `samtools index` does) and gives one line per run of reads on one reference, as htslib's
//! `cram_index_build_multiref` writes them.

use crate::align::cram::{CramInfo, ExtentDecoder, RawContainer, open_for_index, slice_heads};
use crate::atomic::AtomicFile;
use crate::error::IoError;
use crate::pool::Ordered;
use gloci_core::thousands;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// One line of a `.crai`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CraiEntry {
    /// The reference id (-1 for unplaced reads).
    pub ref_id: i32,
    /// The first reference position covered (1-based; 0 for unplaced reads).
    pub start: i64,
    pub span: i64,
    /// Where the container starts in the CRAM file.
    pub container_offset: u64,
    /// Where the slice starts in the container's data (after its header).
    pub slice_offset: u64,
    /// The slice's bytes.
    pub slice_size: u64,
}

/// The entries of a `.crai`, in file order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CramIndex {
    pub entries: Vec<CraiEntry>,
}

/// `CRAM.crai`.
pub fn crai_path(cram: &Path) -> PathBuf {
    let mut s: OsString = cram.as_os_str().to_owned();
    s.push(".crai");
    PathBuf::from(s)
}

/// The `.crai` of `cram` (`x.cram.crai`, else `x.crai`), and whether it is older than the CRAM
/// (D26, D107).
pub fn find_crai(cram: &Path) -> Option<(PathBuf, bool)> {
    let data = std::fs::metadata(cram).and_then(|m| m.modified()).ok();
    let mut candidates = vec![crai_path(cram)];
    if cram.extension().is_some_and(|e| e == "cram") {
        candidates.push(cram.with_extension("crai"));
    }
    candidates.into_iter().find_map(|p| {
        let meta = std::fs::metadata(&p).ok().filter(|m| m.is_file())?;
        let stale = matches!((meta.modified().ok(), data), (Some(i), Some(d)) if i < d);
        Some((p, stale))
    })
}

impl CramIndex {
    /// Read the `.crai` at `path`.
    pub fn read(path: &Path) -> Result<CramIndex, IoError> {
        let name = path.display().to_string();
        let f = File::open(path).map_err(|source| IoError::Open {
            path: name.clone(),
            source,
        })?;
        let mut text = Vec::new();
        flate2::read::MultiGzDecoder::new(BufReader::new(f))
            .read_to_end(&mut text)
            .map_err(|source| IoError::Read {
                path: name.clone(),
                source,
            })?;
        let mut entries = Vec::new();
        for (i, line) in text.split(|&b| b == b'\n').enumerate() {
            if line.is_empty() {
                continue;
            }
            let bad = || IoError::Detect {
                path: name.clone(),
                msg: format!(
                    "{name} line {}: a .crai line has six whole numbers (reference, start, span, \
container offset, slice offset, slice size); rebuild the index: gloci index CRAM",
                    thousands(i as u64 + 1)
                ),
            };
            let f: Vec<i64> = line
                .split(|&b| b == b'\t')
                .map(|v| std::str::from_utf8(v).ok().and_then(|v| v.parse().ok()))
                .collect::<Option<_>>()
                .ok_or_else(bad)?;
            let [r, start, span, c, s, z] = f[..] else {
                return Err(bad());
            };
            let u = |v: i64| u64::try_from(v).map_err(|_| bad());
            entries.push(CraiEntry {
                ref_id: i32::try_from(r).map_err(|_| bad())?,
                start,
                span,
                container_offset: u(c)?,
                slice_offset: u(s)?,
                slice_size: u(z)?,
            });
        }
        Ok(CramIndex { entries })
    }

    /// The `.crai` text (before compression).
    pub fn text(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for e in &self.entries {
            out.extend_from_slice(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    e.ref_id, e.start, e.span, e.container_offset, e.slice_offset, e.slice_size
                )
                .as_bytes(),
            );
        }
        out
    }

    /// Index the CRAM at `path`, decoding slices of several references on `threads` workers.
    /// With `reference`, its lengths are checked against the CRAM's header (D109); indexing
    /// itself reads no reference bases.
    pub fn build(
        path: &Path,
        reference: Option<&Path>,
        threads: NonZeroUsize,
    ) -> Result<CramIndex, IoError> {
        let name = path.display().to_string();
        let f = File::open(path).map_err(|source| IoError::Open {
            path: name.clone(),
            source,
        })?;
        let src: Box<dyn BufRead + Send> = Box::new(BufReader::with_capacity(1 << 16, f));
        let (info, mut frames) = open_for_index(&name, src, reference)?;
        let mut entries = Vec::new();
        let n = threads.get();
        if n == 1 {
            let mut dec = ExtentDecoder::new(Arc::clone(&info));
            while let Some(c) = frames.next()? {
                entries.extend(index_container(&mut dec, c)?);
            }
        } else {
            let decoders = (0..n)
                .map(|_| ExtentDecoder::new(Arc::clone(&info)))
                .collect();
            let mut pool = Ordered::new(decoders, index_container)?;
            let cap = 2 * n as u64;
            let mut reading = true;
            loop {
                while reading && pool.in_flight() < cap {
                    match frames.next() {
                        Ok(Some(c)) => pool.submit(c),
                        Ok(None) => reading = false,
                        Err(e) => {
                            pool.submit_result(Err(e));
                            reading = false;
                        }
                    }
                }
                match pool.next() {
                    Some(r) => entries.extend(r?),
                    None => break,
                }
            }
        }
        check_order(&name, &info, &entries)?;
        Ok(CramIndex { entries })
    }
}

/// The usage error for a CRAM that is not sorted (htslib's `cram_index_build` refuses it too).
fn unsorted(name: &str, offset: u64, chrom: &str, at: i64, after: i64) -> IoError {
    IoError::Detect {
        path: name.to_string(),
        msg: format!(
            "{name} is not sorted by position: at byte {} a read at {chrom}:{} follows one at \
{chrom}:{}; sort it first (gloci sort {name} -o sorted.bam, or samtools sort -O cram \
--reference ref.fa -o sorted.cram {name}), then index it",
            thousands(offset),
            thousands(at.max(0) as u64),
            thousands(after.max(0) as u64)
        ),
    }
}

/// The `.crai` lines of container `c`; the container is parsed once.
fn index_container(dec: &mut ExtentDecoder, c: RawContainer) -> Result<Vec<CraiEntry>, IoError> {
    let name = dec.info().name().to_string();
    let heads = slice_heads(&c, &name, dec.info().ref_count())?;
    let extents = dec.extents(&c, &heads)?;
    let mut out = Vec::new();
    for (h, ext) in heads.iter().zip(extents) {
        let entry = |ref_id: i32, start: i64, span: i64| CraiEntry {
            ref_id,
            start,
            span,
            container_offset: c.offset,
            slice_offset: h.landmark as u64,
            slice_size: h.size as u64,
        };
        let Some(ext) = ext else {
            out.push(entry(h.ref_id, i64::from(h.start), i64::from(h.span)));
            continue;
        };
        // htslib's runs: a line per run of reads on one reference, from the first read's start
        // to the furthest end; a read before the previous one on its reference is an error.
        let mut run: Option<(i32, i64, i64)> = None;
        let mut last: Option<(i32, i64)> = None;
        for (rid, start, end) in ext {
            if let Some((lr, lp)) = last
                && lr == rid
                && start < lp
            {
                return Err(unsorted(
                    &name,
                    c.offset,
                    &dec.info().chrom_or_unplaced(rid),
                    start,
                    lp,
                ));
            }
            last = Some((rid, start));
            match &mut run {
                Some((r, _, e)) if *r == rid => *e = (*e).max(end),
                _ => {
                    if let Some((r, s, e)) = run {
                        out.push(entry(r, s, e - s + 1));
                    }
                    run = Some((rid, start, end));
                }
            }
        }
        if let Some((r, s, e)) = run {
            out.push(entry(r, s, e - s + 1));
        }
    }
    Ok(out)
}

/// Refuse entries (in file order) whose start goes back on the reference of the one before:
/// htslib's `cram_index_build` refuses such a file, and region reads trust the index.
fn check_order(name: &str, info: &CramInfo, entries: &[CraiEntry]) -> Result<(), IoError> {
    for w in entries.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        if a.ref_id == b.ref_id && b.start < a.start {
            return Err(unsorted(
                name,
                b.container_offset,
                &info.chrom_or_unplaced(b.ref_id),
                b.start,
                a.start,
            ));
        }
    }
    Ok(())
}

/// Write `PATH.crai` for the CRAM at `path` atomically and return its path. `reference`, when
/// given, is checked against the CRAM's header; it is not needed (see the module notes).
pub fn write_crai(
    path: &Path,
    reference: Option<&Path>,
    threads: NonZeroUsize,
) -> Result<PathBuf, IoError> {
    let index = CramIndex::build(path, reference, threads)?;
    let dest = crai_path(path);
    let werr = |source| IoError::Write {
        path: dest.display().to_string(),
        source,
    };
    let mut af = AtomicFile::create(&dest).map_err(werr)?;
    {
        let mut gz = flate2::write::GzEncoder::new(&mut af, flate2::Compression::default());
        gz.write_all(&index.text()).map_err(werr)?;
        gz.finish().map_err(werr)?;
    }
    af.commit().map_err(werr)?;
    Ok(dest)
}
