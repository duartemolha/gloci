//! Copying BGZF blocks without recompressing them: the loop behind `vcf concat --naive`
//! (D74) and `gloci cat` for BAM (D123).
//!
//! A BGZF file is a series of independent gzip members ("blocks"), so the records of
//! several files can be joined by copying their compressed bytes. [`copy_from`] reads the
//! stream one block at a time and never seeks, so it works on stdin and holds at most one
//! block (stored and inflated) in memory. Of the data, everything from an uncompressed
//! offset `start` on goes to the output: the block that holds `start` is compressed again
//! from there (the only recompression) and every later block is copied byte for byte.
//! Empty blocks, which includes each file's EOF marker, are dropped; the caller writes one
//! [`BGZF_EOF`] at the very end.

use crate::input::BGZF_EOF;
use std::fmt;
use std::io::{self, Read, Write};

/// The largest block BGZF allows, stored or inflated.
const MAX_BLOCK: usize = 1 << 16;

/// One BGZF block as it is stored.
pub struct Block {
    /// The block's bytes, header and trailer included.
    pub raw: Vec<u8>,
    /// Size of the data when inflated.
    pub isize: usize,
}

/// The reason inside an `InvalidData` error for a stream that is not well-formed BGZF.
#[derive(Debug)]
pub struct MalformedBgzf(pub String);

impl fmt::Display for MalformedBgzf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for MalformedBgzf {}

fn malformed(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, MalformedBgzf(msg.into()))
}

/// What is wrong with the BGZF stream, if `e` says the stream is malformed (as opposed to a
/// read or write that failed).
pub fn malformed_reason(e: &io::Error) -> Option<&str> {
    e.get_ref()
        .and_then(|inner| inner.downcast_ref::<MalformedBgzf>())
        .map(|m| m.0.as_str())
}

/// A writer that remembers whether a write to it failed, so a caller can tell a write
/// error from a read error in the `io::Error` that [`copy_from`] returns.
pub struct TrackedWrite<W: Write> {
    pub inner: W,
    pub failed: bool,
}

impl<W: Write> TrackedWrite<W> {
    pub fn new(inner: W) -> Self {
        TrackedWrite {
            inner,
            failed: false,
        }
    }
}

impl<W: Write> Write for TrackedWrite<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner.write(buf).inspect_err(|_| self.failed = true)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush().inspect_err(|_| self.failed = true)
    }
}

/// Where one block of the input is and how big it is, as [`copy_from_observed`] reports it.
#[derive(Debug, Clone, Copy)]
pub struct BlockInfo {
    /// Offset of the block in the input.
    pub offset: u64,
    /// Stored size.
    pub len: usize,
    /// Size when inflated; 0 for an empty block.
    pub isize: usize,
}

/// What a copy did.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct CopyStats {
    /// The input ended with the BGZF EOF marker.
    pub eof_marker: bool,
    /// Size of the whole input when inflated.
    pub uncompressed: u64,
    /// Blocks written as they were stored.
    pub blocks_copied: u64,
    /// Blocks compressed again (0 or 1).
    pub blocks_recompressed: u64,
    /// The most bytes held at once for one block: stored, inflated and compressed again.
    pub peak_block_bytes: usize,
}

fn le32(b: &[u8]) -> u32 {
    u32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

fn read_full(src: &mut dyn Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut got = 0;
    while got < buf.len() {
        match src.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(got)
}

/// The next block of `src`; `None` at a clean end of the stream.
pub fn read_block(src: &mut dyn Read) -> io::Result<Option<Block>> {
    let truncated = || malformed("the file ends inside a BGZF block");
    let mut head = [0u8; 12];
    let got = read_full(src, &mut head)?;
    if got == 0 {
        return Ok(None);
    }
    if got < head.len() {
        return Err(truncated());
    }
    if head[..4] != [0x1f, 0x8b, 8, 4] {
        return Err(malformed(
            "a block does not start with a BGZF header (is the file BGZF-compressed?)",
        ));
    }
    let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
    let mut raw = head.to_vec();
    raw.resize(12 + xlen, 0);
    if read_full(src, &mut raw[12..])? < xlen {
        return Err(truncated());
    }
    // The BC subfield holds the block size minus one.
    let mut bsize = None;
    let mut extra = &raw[12..];
    while extra.len() >= 4 {
        let len = u16::from_le_bytes([extra[2], extra[3]]) as usize;
        if extra[..2] == *b"BC" && len == 2 && extra.len() >= 6 {
            bsize = Some(u16::from_le_bytes([extra[4], extra[5]]) as usize + 1);
        }
        extra = &extra[(4 + len).min(extra.len())..];
    }
    let total = bsize.ok_or_else(|| {
        malformed("a gzip block has no BGZF size field (is the file plain gzip?)")
    })?;
    if total < 12 + xlen + 8 {
        return Err(malformed("a BGZF block is too small"));
    }
    let start = raw.len();
    raw.resize(total, 0);
    if read_full(src, &mut raw[start..])? < total - start {
        return Err(truncated());
    }
    let isize = le32(&raw[total - 4..]) as usize;
    if isize > MAX_BLOCK {
        return Err(malformed("a BGZF block claims more than 64 KiB of data"));
    }
    Ok(Some(Block { raw, isize }))
}

/// The data of `block`.
pub fn inflate(block: &Block) -> io::Result<Vec<u8>> {
    let xlen = u16::from_le_bytes([block.raw[10], block.raw[11]]) as usize;
    let body = &block.raw[12 + xlen..block.raw.len() - 8];
    let mut out = Vec::with_capacity(block.isize);
    flate2::read::DeflateDecoder::new(body)
        .take(block.isize as u64 + 1)
        .read_to_end(&mut out)
        .map_err(|e| malformed(format!("a BGZF block does not inflate: {e}")))?;
    if out.len() != block.isize {
        return Err(malformed("a BGZF block's size does not match its content"));
    }
    Ok(out)
}

/// `data` as BGZF blocks, without an EOF marker.
pub fn deflate_blocks(data: &[u8]) -> io::Result<Vec<u8>> {
    let mut w = noodles_bgzf::io::Writer::new(Vec::new());
    w.write_all(data)?;
    let mut out = w.finish()?;
    debug_assert!(out.ends_with(&BGZF_EOF));
    out.truncate(out.len() - BGZF_EOF.len());
    Ok(out)
}

/// Copy the BGZF stream `src` to `out` from uncompressed offset `start` on, as the module
/// comment describes. With `start` 0 every non-empty block is copied as it is.
pub fn copy_from(src: &mut dyn Read, out: &mut dyn Write, start: u64) -> io::Result<CopyStats> {
    copy_from_observed(src, out, start, &mut |_| {})
}

/// [`copy_from`] that also calls `observe` for every block read, empty ones included.
pub fn copy_from_observed(
    src: &mut dyn Read,
    out: &mut dyn Write,
    start: u64,
    observe: &mut dyn FnMut(BlockInfo),
) -> io::Result<CopyStats> {
    let mut stats = CopyStats::default();
    let mut offset = 0u64;
    let mut last_was_eof = false;
    while let Some(block) = read_block(src)? {
        last_was_eof = block.raw == BGZF_EOF;
        observe(BlockInfo {
            offset,
            len: block.raw.len(),
            isize: block.isize,
        });
        offset += block.raw.len() as u64;
        if block.isize == 0 {
            continue;
        }
        let begin = stats.uncompressed;
        stats.uncompressed += block.isize as u64;
        let mut held = block.raw.len();
        if stats.uncompressed <= start {
            // Wholly before `start`.
        } else if begin >= start {
            out.write_all(&block.raw)?;
            stats.blocks_copied += 1;
        } else {
            // `start` is inside this block: only what follows it is output, as a new block.
            let data = inflate(&block)?;
            let packed = deflate_blocks(&data[(start - begin) as usize..])?;
            held += data.len() + packed.len();
            out.write_all(&packed)?;
            stats.blocks_recompressed += 1;
        }
        stats.peak_block_bytes = stats.peak_block_bytes.max(held);
    }
    stats.eof_marker = last_was_eof;
    Ok(stats)
}

/// A reader of the data of a BGZF stream that keeps the stored bytes of every block it has
/// read, so that a caller can read a header from a stream (stdin, say), learn how long the
/// header is, and then hand the same stream, from its first byte, to [`copy_from`]: see
/// [`HeadReader::finish`].
pub struct HeadReader<R: Read> {
    src: R,
    stored: Vec<u8>,
    data: Vec<u8>,
    at: usize,
    produced: u64,
}

impl<R: Read> HeadReader<R> {
    pub fn new(src: R) -> Self {
        HeadReader {
            src,
            stored: Vec::new(),
            data: Vec::new(),
            at: 0,
            produced: 0,
        }
    }

    /// The bytes of the stream read so far (stored form) and what remains of it, plus how
    /// many bytes of data were handed out: the header's length once it is read. Chain the
    /// first two (`Cursor::new(stored).chain(rest)`) to get the whole stream again.
    pub fn finish(self) -> (Vec<u8>, R, u64) {
        (self.stored, self.src, self.produced)
    }
}

impl<R: Read> Read for HeadReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.at == self.data.len() {
            let Some(block) = read_block(&mut self.src)? else {
                return Ok(0);
            };
            self.data = inflate(&block)?;
            self.stored.extend_from_slice(&block.raw);
            self.at = 0;
        }
        let n = buf.len().min(self.data.len() - self.at);
        buf[..n].copy_from_slice(&self.data[self.at..self.at + n]);
        self.at += n;
        self.produced += n as u64;
        Ok(n)
    }
}
