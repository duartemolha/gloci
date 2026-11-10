//! Building and reading `.fai` and `.gzi` files. Output equals samtools
//! 1.21 (`samtools faidx`) and htslib's `bgzip -r` index byte for byte.

use super::{fai_path, gzi_path, gzip_error};
use crate::atomic::AtomicFile;
use crate::detect::{Compression, Format};
use crate::error::{IoError, Loc};
use crate::input::{InputOptions, open_input};
use gloci_core::thousands;
use memchr::memchr;
use std::collections::HashSet;
use std::io::{BufRead, Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// A header line longer than this is not a FASTA header.
const MAX_HEADER: usize = 1 << 20;

/// One `.fai` row: `name len offset line_bases line_width`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
    /// The first word of the header line.
    pub name: Vec<u8>,
    /// Bases in the sequence.
    pub len: u64,
    /// Byte offset (in the uncompressed file) of the first base.
    pub offset: u64,
    /// Bases per full line.
    pub line_bases: u64,
    /// Bytes per full line, line terminator included.
    pub line_width: u64,
}

/// The `.fai` text of `entries`.
pub fn fai_bytes(entries: &[FaiEntry]) -> Vec<u8> {
    let mut out = Vec::new();
    for e in entries {
        out.extend_from_slice(&e.name);
        out.extend_from_slice(
            format!(
                "\t{}\t{}\t{}\t{}\n",
                e.len, e.offset, e.line_bases, e.line_width
            )
            .as_bytes(),
        );
    }
    out
}

#[derive(PartialEq, Eq, Clone, Copy)]
enum Kind {
    /// No byte of the line seen yet (blank so far).
    Blank,
    Header,
    Seq,
}

struct Cur {
    name: Vec<u8>,
    offset: u64,
    header_line: u64,
    len: u64,
    line_bases: u64,
    line_width: u64,
    term_len: u64,
    short_seen: bool,
    blank_seen: bool,
}

struct Builder<'a> {
    name: &'a str,
    /// Bytes consumed before the current chunk.
    pos: u64,
    line_no: u64,
    kind: Kind,
    header: Vec<u8>,
    // tallies of the sequence line being read
    bytes: u64,
    graph: u64,
    bad: bool,
    last_cr: bool,
    cur: Option<Cur>,
    entries: Vec<FaiEntry>,
    seen: HashSet<Vec<u8>>,
    empties: u64,
}

impl<'a> Builder<'a> {
    fn new(name: &'a str) -> Self {
        Builder {
            name,
            pos: 0,
            line_no: 1,
            kind: Kind::Blank,
            header: Vec::new(),
            bytes: 0,
            graph: 0,
            bad: false,
            last_cr: false,
            cur: None,
            entries: Vec::new(),
            seen: HashSet::new(),
            empties: 0,
        }
    }

    fn err(&self, line: u64, msg: String) -> IoError {
        IoError::Parse {
            path: self.name.to_string(),
            loc: Loc::Line(line),
            msg,
        }
    }

    fn feed(&mut self, chunk: &[u8]) -> Result<(), IoError> {
        let mut i = 0;
        while i < chunk.len() {
            match memchr(b'\n', &chunk[i..]) {
                Some(nl) => {
                    self.piece(&chunk[i..i + nl])?;
                    let next = self.pos + (i + nl + 1) as u64;
                    self.end_line(next, true)?;
                    i += nl + 1;
                }
                None => {
                    self.piece(&chunk[i..])?;
                    i = chunk.len();
                }
            }
        }
        self.pos += chunk.len() as u64;
        Ok(())
    }

    /// Part of the current line (no `\n` in it).
    fn piece(&mut self, seg: &[u8]) -> Result<(), IoError> {
        if seg.is_empty() {
            return Ok(());
        }
        let mut seg = seg;
        if self.kind == Kind::Blank {
            if seg[0] == b'>' {
                self.kind = Kind::Header;
                seg = &seg[1..];
            } else {
                self.kind = Kind::Seq;
            }
            self.bytes = 0;
            self.graph = 0;
            self.bad = false;
            self.last_cr = false;
            self.header.clear();
            if self.kind == Kind::Header {
                self.bytes = 1;
            }
        }
        match self.kind {
            Kind::Header => {
                if self.header.len() + seg.len() > MAX_HEADER {
                    return Err(self.err(
                        self.line_no,
                        "header line is longer than 1 MiB; this is not a FASTA header, check that the file is FASTA".into(),
                    ));
                }
                self.header.extend_from_slice(seg);
            }
            Kind::Seq => {
                if self.last_cr && !seg.is_empty() {
                    self.bad = true;
                }
                self.bytes += seg.len() as u64;
                let (body, cr) = match seg.split_last() {
                    Some((b'\r', rest)) => (rest, true),
                    _ => (seg, false),
                };
                self.last_cr = cr;
                let graph = body.iter().filter(|b| b.is_ascii_graphic()).count();
                if graph != body.len() {
                    self.bad = true;
                }
                self.graph += graph as u64;
            }
            Kind::Blank => {}
        }
        Ok(())
    }

    /// The line ends (`terminated`: at a `\n`; else at the end of the file).
    /// `next` is the byte offset after the line.
    fn end_line(&mut self, next: u64, terminated: bool) -> Result<(), IoError> {
        let line = self.line_no;
        match self.kind {
            Kind::Blank => {
                if let Some(c) = &mut self.cur {
                    c.blank_seen = true;
                }
            }
            Kind::Header => {
                self.finish_entry(false)?;
                self.start_entry(next, line)?;
                if !terminated {
                    // A header with nothing after it: the entry has no sequence.
                    return self.finish_entry(true);
                }
            }
            Kind::Seq => self.seq_line(line, terminated)?,
        }
        self.kind = Kind::Blank;
        self.line_no += 1;
        Ok(())
    }

    fn start_entry(&mut self, offset: u64, line: u64) -> Result<(), IoError> {
        let end = self
            .header
            .iter()
            .position(|b| matches!(b, b' ' | b'\t' | b'\r' | 0x0b | 0x0c))
            .unwrap_or(self.header.len());
        let name = self.header[..end].to_vec();
        if name.is_empty() {
            return Err(self.err(
                line,
                "a header line has no name after '>'; name every sequence".into(),
            ));
        }
        if !self.seen.insert(name.clone()) {
            return Err(self.err(
                line,
                format!(
                    "sequence name {} is used more than once; names must be unique to index the file (samtools keeps the first and warns); rename the duplicate",
                    String::from_utf8_lossy(&name)
                ),
            ));
        }
        self.cur = Some(Cur {
            name,
            offset,
            header_line: line,
            len: 0,
            line_bases: 0,
            line_width: 0,
            term_len: 0,
            short_seen: false,
            blank_seen: false,
        });
        Ok(())
    }

    fn seq_line(&mut self, line: u64, terminated: bool) -> Result<(), IoError> {
        let blank = self.bytes == 0 || (self.bytes == 1 && self.last_cr && self.graph == 0);
        if blank {
            if let Some(c) = &mut self.cur {
                c.blank_seen = true;
            }
            return Ok(());
        }
        let Some(c) = self.cur.as_mut() else {
            return Err(IoError::Parse {
                path: self.name.to_string(),
                loc: Loc::Line(line),
                msg: "expected a header line starting with '>' before the sequence; this does not look like FASTA".into(),
            });
        };
        let err = |msg: String| IoError::Parse {
            path: self.name.to_string(),
            loc: Loc::Line(line),
            msg,
        };
        if self.bad {
            return Err(err(
                "a sequence line contains a space, a control character or a non-ASCII byte; remove it (only a final CR before the line break is allowed)".into(),
            ));
        }
        if c.blank_seen {
            return Err(err(format!(
                "sequence {} continues after a blank line; remove the blank lines inside the sequence",
                String::from_utf8_lossy(&c.name)
            )));
        }
        let (bases, width) = (self.graph, self.bytes + 1);
        if c.len == 0 {
            c.line_bases = bases;
            c.line_width = width;
            c.term_len = width - bases;
        } else {
            let sname = String::from_utf8_lossy(&c.name).into_owned();
            if c.short_seen {
                return Err(err(format!(
                    "a line of {sname} follows a line shorter than {} bases; all lines but the last must have equal length (rewrap the file, for example with: seqkit seq -w {})",
                    c.line_bases, c.line_bases
                )));
            }
            if bases > c.line_bases {
                return Err(err(format!(
                    "this line of {sname} has {} bases but the first has {}; all lines but the last must have equal length (rewrap the file, for example with: seqkit seq -w {})",
                    bases, c.line_bases, c.line_bases
                )));
            }
            if terminated && width - bases != c.term_len {
                return Err(err(format!(
                    "line ends of {sname} are mixed (CRLF and LF); use one line ending in each sequence"
                )));
            }
            if bases < c.line_bases {
                c.short_seen = true;
            }
        }
        c.len += bases;
        Ok(())
    }

    fn finish_entry(&mut self, at_eof: bool) -> Result<(), IoError> {
        let Some(c) = self.cur.take() else {
            return Ok(());
        };
        if c.len == 0 {
            if at_eof {
                return Err(self.err(
                    c.header_line,
                    format!(
                        "the last sequence, {}, has no bases; remove its header or add the sequence",
                        String::from_utf8_lossy(&c.name)
                    ),
                ));
            }
            // htslib leaves empty sequences out of the index.
            self.empties += 1;
            return Ok(());
        }
        self.entries.push(FaiEntry {
            name: c.name,
            len: c.len,
            offset: c.offset,
            line_bases: c.line_bases,
            line_width: c.line_width,
        });
        Ok(())
    }

    fn finish(mut self) -> Result<Vec<FaiEntry>, IoError> {
        if self.kind != Kind::Blank {
            let next = self.pos;
            self.end_line(next, false)?;
        }
        self.finish_entry(true)?;
        if self.entries.is_empty() {
            return Err(IoError::Parse {
                path: self.name.to_string(),
                loc: Loc::Line(1),
                msg: "no sequences; the file is empty or is not FASTA".into(),
            });
        }
        if self.empties > 0 {
            log::warn!(
                "{}: {} empty sequence(s) are left out of the index, as samtools faidx does",
                self.name,
                thousands(self.empties)
            );
        }
        Ok(self.entries)
    }
}

/// Index the (uncompressed) FASTA text in `src`; `name` is used in errors.
///
/// Matches `samtools faidx`: names are the first word of the header, every
/// line but the last of a sequence has the same width, CRLF is allowed
/// (the CR is not a base), and empty sequences are left out. gloci is
/// stricter in two ways, each an error rather than a wrong index: a duplicate
/// or empty name, and a sequence line holding a space or other
/// non-graphic byte.
pub fn build_fai(src: &mut dyn BufRead, name: &str) -> Result<Vec<FaiEntry>, IoError> {
    let mut b = Builder::new(name);
    loop {
        let chunk = match src.fill_buf() {
            Ok(c) => c,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(source) => {
                return Err(IoError::Read {
                    path: name.to_string(),
                    source,
                });
            }
        };
        if chunk.is_empty() {
            break;
        }
        let n = chunk.len();
        b.feed(chunk)?;
        src.consume(n);
    }
    b.finish()
}

/// Parse `.fai` text. Rows have five columns; a row that cannot be used for
/// random access (a zero line width, duplicate names) is an error naming its line.
pub fn parse_fai(text: &[u8], name: &str) -> Result<Vec<FaiEntry>, IoError> {
    let mut out: Vec<FaiEntry> = Vec::new();
    let mut seen = HashSet::new();
    for (i, line) in text.split(|&b| b == b'\n').enumerate() {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        let err = |msg: String| IoError::Parse {
            path: name.to_string(),
            loc: Loc::Line(i as u64 + 1),
            msg: format!("{msg}; rebuild the index with: gloci index <the FASTA file>"),
        };
        let f: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        if f.len() != 5 {
            return Err(err(format!(
                "expected 5 tab-separated columns (name, length, offset, line bases, line width), found {}",
                f.len()
            )));
        }
        let num = |k: usize, what: &str| -> Result<u64, IoError> {
            std::str::from_utf8(f[k])
                .ok()
                .and_then(|s| s.parse::<u64>().ok())
                .ok_or_else(|| err(format!("{what} is not a number")))
        };
        let e = FaiEntry {
            name: f[0].to_vec(),
            len: num(1, "the length")?,
            offset: num(2, "the offset")?,
            line_bases: num(3, "the line-bases column")?,
            line_width: num(4, "the line-width column")?,
        };
        if e.name.is_empty() {
            return Err(err("the name is empty".into()));
        }
        if e.len > 0 && (e.line_bases == 0 || e.line_width < e.line_bases) {
            return Err(err(
                "line bases is 0 or larger than line width for a non-empty sequence".into(),
            ));
        }
        if !seen.insert(e.name.clone()) {
            return Err(err(format!(
                "the name {} appears twice",
                String::from_utf8_lossy(&e.name)
            )));
        }
        out.push(e);
    }
    Ok(out)
}

/// The `.gzi` of a BGZF file: for each non-empty block but the first, its
/// compressed offset and the uncompressed bytes before it (as `bgzip -r`).
fn gzi_entries(path: &Path, name: &str) -> Result<Vec<(u64, u64)>, IoError> {
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let mut f = std::fs::File::open(path).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let total = f.metadata().map_err(read_err)?.len();
    let bad = |msg: &str| IoError::Detect {
        path: name.to_string(),
        msg: format!("{name}: {msg}; the file is probably corrupt"),
    };
    let mut entries = Vec::new();
    let (mut coff, mut uoff, mut first) = (0u64, 0u64, true);
    let mut buf = Vec::new();
    f.seek(SeekFrom::Start(0)).map_err(read_err)?;
    while coff < total {
        let mut head = [0u8; 12];
        f.read_exact(&mut head)
            .map_err(|_| bad("a BGZF block header is cut off"))?;
        if head[0] != 0x1f || head[1] != 0x8b || head[3] & 4 == 0 {
            return Err(bad("a BGZF block header is invalid"));
        }
        let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
        let mut extra = vec![0u8; xlen];
        f.read_exact(&mut extra)
            .map_err(|_| bad("a BGZF block header is cut off"))?;
        let mut bsize = None;
        let mut p = 0;
        while p + 4 <= extra.len() {
            let slen = u16::from_le_bytes([extra[p + 2], extra[p + 3]]) as usize;
            if extra[p] == b'B' && extra[p + 1] == b'C' && slen == 2 && p + 6 <= extra.len() {
                bsize = Some(u16::from_le_bytes([extra[p + 4], extra[p + 5]]) as usize + 1);
            }
            p += 4 + slen;
        }
        let bsize = bsize.ok_or_else(|| bad("a block has no BGZF size field"))?;
        let head_len = 12 + xlen;
        if bsize < head_len + 8 {
            return Err(bad("a BGZF block size is invalid"));
        }
        buf.resize(bsize - head_len, 0);
        f.read_exact(&mut buf)
            .map_err(|_| bad("a BGZF block is cut off"))?;
        let isize = u32::from_le_bytes([
            buf[buf.len() - 4],
            buf[buf.len() - 3],
            buf[buf.len() - 2],
            buf[buf.len() - 1],
        ]) as u64;
        if isize > 0 {
            if !first {
                entries.push((coff, uoff));
            }
            first = false;
            uoff += isize;
        }
        coff += bsize as u64;
    }
    Ok(entries)
}

fn gzi_bytes(entries: &[(u64, u64)]) -> Vec<u8> {
    let mut out = Vec::with_capacity(8 + 16 * entries.len());
    out.extend_from_slice(&(entries.len() as u64).to_le_bytes());
    for (c, u) in entries {
        out.extend_from_slice(&c.to_le_bytes());
        out.extend_from_slice(&u.to_le_bytes());
    }
    out
}

fn publish(dest: &Path, bytes: &[u8]) -> Result<(), IoError> {
    let werr = |source| IoError::Write {
        path: dest.display().to_string(),
        source,
    };
    let mut af = AtomicFile::create(dest).map_err(werr)?;
    af.write_all(bytes).map_err(werr)?;
    af.commit().map_err(werr)
}

/// Write `PATH.fai` (and `PATH.gzi` when `path` is BGZF) atomically; returns
/// the paths written. A plain-gzip FASTA is a usage error with a `bgzip` recipe.
pub fn write_fasta_index(path: &Path, threads: NonZeroUsize) -> Result<Vec<PathBuf>, IoError> {
    let name = path.display().to_string();
    let mut input = open_input(
        &name,
        &InputOptions {
            threads,
            format: Some(Format::Fasta),
        },
    )?;
    if input.compression == Compression::Gzip {
        return Err(gzip_error(&name));
    }
    let entries = build_fai(&mut input.reader, &name)?;
    let gzi = if input.compression == Compression::Bgzf {
        Some(gzi_bytes(&gzi_entries(path, &name)?))
    } else {
        None
    };
    drop(input);
    let mut written = Vec::new();
    if let Some(g) = gzi {
        let p = gzi_path(path);
        publish(&p, &g)?;
        written.push(p);
    }
    let p = fai_path(path);
    publish(&p, &fai_bytes(&entries))?;
    written.insert(0, p);
    Ok(written)
}
