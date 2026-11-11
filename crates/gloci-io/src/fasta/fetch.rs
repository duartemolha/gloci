//! Random access to a FASTA through its `.fai` (and `.gzi` for BGZF).

use super::index::{FaiEntry, parse_fai};
use super::{find_fai, gzi_path, gzip_error};
use crate::detect::{Compression, detect_compression};
use crate::error::{IoError, UsageError};
use gloci_core::Interval;
use memchr::memchr2;
use noodles_bgzf::{self as bgzf, gzi};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Most bases handed to the callback of one `fetch` call (1 MiB).
pub const FETCH_CHUNK: u64 = 1 << 20;

struct Shared {
    entries: Vec<FaiEntry>,
    by_name: HashMap<Vec<u8>, usize>,
    fai: PathBuf,
    gzi: Option<gzi::Index>,
}

enum Source {
    Plain(File),
    Bgzf(Box<bgzf::io::Reader<File>>),
}

/// A FASTA opened for random access.
///
/// - [`IndexedFasta::open`] needs `PATH.fai`, plus `PATH.gzi` for a BGZF file.
///   A missing one is a usage error naming `gloci index PATH` (exit 2); a
///   plain-gzip file is a usage error with a `bgzip` recipe.
/// - [`IndexedFasta::fetch`] passes the bases of a half-open interval to a
///   callback in chunks of at most [`FETCH_CHUNK`], in order, free of line
///   breaks, with the case as stored. [`IndexedFasta::fetch_rev`] passes the
///   same chunks, last first. Every line break met is checked against the
///   index: an edited FASTA is an error saying `does not match`. Stale-index
///   detection covers only the bytes a fetch reads: a same-length edit, or a
///   sequence appended after the last entry, is not detected (the trust level
///   of samtools/htslib). The "index older than FASTA" mtime warning covers
///   most real edits.
/// - For BGZF, the `.gzi` maps an uncompressed offset to the block holding
///   it; that block is decoded and the reader then reads on from there.
pub struct IndexedFasta {
    name: String,
    path: PathBuf,
    shared: Arc<Shared>,
    src: Source,
    buf: Vec<u8>,
}

impl IndexedFasta {
    pub fn open(path: &Path) -> Result<Self, IoError> {
        let name = path.display().to_string();
        let mut f = File::open(path).map_err(|source| IoError::Open {
            path: name.clone(),
            source,
        })?;
        let mut head = [0u8; 18];
        let mut got = 0;
        while got < head.len() {
            match f.read(&mut head[got..]) {
                Ok(0) => break,
                Ok(k) => got += k,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(source) => {
                    return Err(IoError::Read { path: name, source });
                }
            }
        }
        let compression = detect_compression(&head[..got]);
        if compression == Compression::Gzip {
            return Err(gzip_error(&name));
        }
        let (fai, stale) = find_fai(path)?;
        if stale {
            crate::index::warn_stale_index(&fai, &name);
        }
        let text = std::fs::read(&fai).map_err(|source| IoError::Read {
            path: fai.display().to_string(),
            source,
        })?;
        let entries = parse_fai(&text, &fai.display().to_string())?;
        let gzi = if compression == Compression::Bgzf {
            let g = gzi_path(path);
            if !g.exists() {
                return Err(UsageError(format!(
                    "{name} is BGZF-compressed and has no index {}; create it with: gloci index {name}",
                    g.display()
                ))
                .into());
            }
            Some(gzi::fs::read(&g).map_err(|source| IoError::Read {
                path: g.display().to_string(),
                source,
            })?)
        } else {
            None
        };
        let by_name = entries
            .iter()
            .enumerate()
            .map(|(i, e)| (e.name.clone(), i))
            .collect();
        let shared = Arc::new(Shared {
            entries,
            by_name,
            fai,
            gzi,
        });
        Ok(IndexedFasta {
            src: Self::source(f, compression),
            name,
            path: path.to_path_buf(),
            shared,
            buf: Vec::new(),
        })
    }

    fn source(f: File, compression: Compression) -> Source {
        match compression {
            Compression::Bgzf => Source::Bgzf(Box::new(bgzf::io::Reader::new(f))),
            _ => Source::Plain(f),
        }
    }

    /// Another handle on the same file that shares the parsed index.
    pub fn reopen(&self) -> Result<Self, IoError> {
        let f = File::open(&self.path).map_err(|source| IoError::Open {
            path: self.name.clone(),
            source,
        })?;
        let compression = if self.shared.gzi.is_some() {
            Compression::Bgzf
        } else {
            Compression::None
        };
        Ok(IndexedFasta {
            name: self.name.clone(),
            path: self.path.clone(),
            shared: Arc::clone(&self.shared),
            src: Self::source(f, compression),
            buf: Vec::new(),
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// The `.fai` rows, in file order.
    pub fn entries(&self) -> &[FaiEntry] {
        &self.shared.entries
    }

    pub fn len_of(&self, chrom: &[u8]) -> Option<u64> {
        self.shared
            .by_name
            .get(chrom)
            .map(|&i| self.shared.entries[i].len)
    }

    /// Pass the bases of `iv` on `chrom` to `f`, forwards.
    pub fn fetch(
        &mut self,
        chrom: &[u8],
        iv: Interval,
        f: &mut dyn FnMut(&[u8]) -> Result<(), IoError>,
    ) -> Result<(), IoError> {
        self.fetch_dir(chrom, iv, false, f)
    }

    /// The chunks of [`IndexedFasta::fetch`], last first (each chunk is still
    /// in sequence order).
    pub fn fetch_rev(
        &mut self,
        chrom: &[u8],
        iv: Interval,
        f: &mut dyn FnMut(&[u8]) -> Result<(), IoError>,
    ) -> Result<(), IoError> {
        self.fetch_dir(chrom, iv, true, f)
    }

    fn fetch_dir(
        &mut self,
        chrom: &[u8],
        iv: Interval,
        rev: bool,
        f: &mut dyn FnMut(&[u8]) -> Result<(), IoError>,
    ) -> Result<(), IoError> {
        let shown = String::from_utf8_lossy(chrom).into_owned();
        let Some(&idx) = self.shared.by_name.get(chrom) else {
            return Err(IoError::Detect {
                path: self.name.clone(),
                msg: format!(
                    "{}: sequence {shown} is not in the index {}; check the chromosome name",
                    self.name,
                    self.shared.fai.display()
                ),
            });
        };
        let shared = Arc::clone(&self.shared);
        let e = &shared.entries[idx];
        if iv.end > e.len {
            return Err(IoError::Detect {
                path: self.name.clone(),
                msg: format!(
                    "{}: {shown}:{}-{} is past the end of {shown}, which has {} bases",
                    self.name, iv.start, iv.end, e.len
                ),
            });
        }
        if iv.start == iv.end {
            return Ok(());
        }
        let n = (iv.end - iv.start).div_ceil(FETCH_CHUNK);
        for k in 0..n {
            let k = if rev { n - 1 - k } else { k };
            let s = iv.start + k * FETCH_CHUNK;
            let end = (s + FETCH_CHUNK).min(iv.end);
            let mut buf = std::mem::take(&mut self.buf);
            let r = self.read_range(e, &shown, s, end, &mut buf);
            let r = r.and_then(|w| f(&buf[..w]));
            self.buf = buf;
            r?;
        }
        Ok(())
    }

    fn mismatch(&self, shown: &str, what: &str) -> IoError {
        IoError::Detect {
            path: self.name.clone(),
            msg: format!(
                "{}: sequence {shown} does not match the index {} ({what}); rebuild it: gloci index {}",
                self.name,
                self.shared.fai.display(),
                self.name
            ),
        }
    }

    /// Read bases `[s, end)` (at most one chunk) into `buf`, line breaks
    /// removed in place; returns the number of bases.
    fn read_range(
        &mut self,
        e: &FaiEntry,
        shown: &str,
        s: u64,
        end: u64,
        buf: &mut Vec<u8>,
    ) -> Result<usize, IoError> {
        let term_len = e.line_width.wrapping_sub(e.line_bases);
        if e.line_bases == 0 || !(1..=2).contains(&term_len) {
            return Err(self.mismatch(
                shown,
                "its line width is not the line bases plus a line break",
            ));
        }
        let byte_of = |p: u64| -> Option<u64> {
            (p / e.line_bases)
                .checked_mul(e.line_width)?
                .checked_add(p % e.line_bases)?
                .checked_add(e.offset)
        };
        let overflow = || self.mismatch(shown, "an offset overflows");
        let b0 = byte_of(s).ok_or_else(overflow)?;
        let b1 = byte_of(end - 1)
            .and_then(|b| b.checked_add(1))
            .ok_or_else(overflow)?;
        let raw_len = usize::try_from(b1 - b0).map_err(|_| overflow())?;
        buf.clear();
        buf.resize(raw_len, 0);
        let read = match &mut self.src {
            Source::Plain(f) => f.seek(SeekFrom::Start(b0)).and_then(|_| f.read_exact(buf)),
            Source::Bgzf(r) => match self.shared.gzi.as_ref() {
                Some(g) => r
                    .seek_by_uncompressed_position(g, b0)
                    .and_then(|_| r.read_exact(buf)),
                None => Err(std::io::Error::other("no .gzi")),
            },
        };
        if let Err(source) = read {
            return Err(match source.kind() {
                std::io::ErrorKind::UnexpectedEof | std::io::ErrorKind::InvalidData => self
                    .mismatch(
                        shown,
                        "the file is shorter than the index says, or the .gzi is wrong",
                    ),
                _ => IoError::Read {
                    path: self.name.clone(),
                    source,
                },
            });
        }
        let term: &[u8] = if term_len == 2 { b"\r\n" } else { b"\n" };
        let (mut r, mut w) = (0usize, 0usize);
        let mut col = s % e.line_bases;
        let mut line = s / e.line_bases;
        let mut need = end - s;
        while need > 0 {
            let cols = e.line_bases.min(e.len - line * e.line_bases);
            let run = (cols - col).min(need) as usize;
            let Some(bases) = buf.get(r..r + run) else {
                return Err(self.mismatch(shown, "the lines are shorter than the index says"));
            };
            if memchr2(b'\n', b'\r', bases).is_some() {
                return Err(self.mismatch(shown, "a line break is inside a line"));
            }
            buf.copy_within(r..r + run, w);
            r += run;
            w += run;
            need -= run as u64;
            if need > 0 {
                if buf.get(r..r + term.len()) != Some(term) {
                    return Err(self.mismatch(shown, "a line break is missing or different"));
                }
                r += term.len();
                col = 0;
                line += 1;
            }
        }
        Ok(w)
    }
}

impl std::fmt::Debug for IndexedFasta {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedFasta")
            .field("name", &self.name)
            .field("sequences", &self.shared.entries.len())
            .finish_non_exhaustive()
    }
}
