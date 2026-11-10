//! CRAM input (D109, D110): the file's framing, the reference, the decode pool, and the records
//! transcoded to BAM bytes.
//!
//! - **Framing** is gloci's own: the file definition, the header container and each data
//!   container are read in order and CRC-checked, and every slice header is checked (landmarks,
//!   block counts, reference ids, spans) before noodles-cram decodes anything, so its `expect`s on
//!   a bad reference id or an empty repository stay unreachable. A file that ends without the
//!   end-of-file container is an error naming it (`#cram-eof-container`).
//! - **The reference** is `--reference FASTA` only (bedtools `--cram-ref`): no `REF_PATH`, no
//!   cache directory, no network. Each `@SQ` the `.fai` names must have its `LN`, checked before
//!   any record. [`RefAdapter`] reads whole sequences through [`IndexedFasta`], uppercased as
//!   htslib reads them, and never answers `None`: a missing name is `Some(Err)` with the D109
//!   message.
//! - **Each decoder** owns a noodles repository and clears it when the slice's reference id
//!   changes, so it holds one sequence (a slice holding several references holds those).
//! - **The pool** (D110): containers are read in order, decoded on up to `threads` workers with
//!   at most `2 × threads` in flight, and released in order, so the output does not depend on
//!   the thread count.
//! - **Transcoding**: each record is encoded as BAM bytes once, with htslib's `bin` and MAPQ 0
//!   for unmapped reads (CRAM does not store their MAPQ), and the bytes go through the BAM
//!   reader's framing ([`BamFrames`]), so interval commands and the every-record source treat a
//!   CRAM exactly as the BAM it encodes. MD and NM are not regenerated (`#cram-md-nm`).

use super::AlignOptions;
use super::bam::{BamFrames, BamReader};
use super::encode::reg2bin;
use super::header::AlnHeader;
use super::source::{AlnSource, BamSource};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::Diagnostics;
use crate::error::IoError;
use crate::fasta::IndexedFasta;
use crate::input::{Input, is_regular_file};
use crate::record::{Header, ReaderOptions, RecordBatch, RecordReader};
use gloci_core::{ChromTable, Interval, thousands};
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_sam::alignment::Record as _;
use noodles_sam::alignment::io::Write as _;
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::io::{self, BufRead, Read};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, mpsc};
use std::thread::JoinHandle;

/// The 38 bytes of a CRAM 3 end-of-file container.
pub(crate) const CRAM_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// Most records a slice may claim. noodles allocates a record for each before decoding, so a
/// corrupt count must not ask for gigabytes; samtools writes 10,000 per slice by default.
const MAX_SLICE_RECORDS: usize = 1 << 22;
/// Most slices a container may claim.
const MAX_LANDMARKS: usize = 1 << 20;

// ---- numbers ----

/// An ITF8 integer (1 to 5 bytes).
fn itf8(src: &mut &[u8]) -> Option<i32> {
    let (&b0, rest) = src.split_first()?;
    let extra = match b0 {
        0x00..=0x7f => 0,
        0x80..=0xbf => 1,
        0xc0..=0xdf => 2,
        0xe0..=0xef => 3,
        _ => 4,
    };
    let b = rest.get(..extra)?;
    let v = match extra {
        0 => u32::from(b0),
        1 => (u32::from(b0 & 0x7f) << 8) | u32::from(b[0]),
        2 => (u32::from(b0 & 0x3f) << 16) | (u32::from(b[0]) << 8) | u32::from(b[1]),
        3 => {
            (u32::from(b0 & 0x1f) << 24)
                | (u32::from(b[0]) << 16)
                | (u32::from(b[1]) << 8)
                | u32::from(b[2])
        }
        _ => {
            (u32::from(b0 & 0x0f) << 28)
                | (u32::from(b[0]) << 20)
                | (u32::from(b[1]) << 12)
                | (u32::from(b[2]) << 4)
                | u32::from(b[3] & 0x0f)
        }
    };
    *src = &rest[extra..];
    Some(v as i32)
}

/// An LTF8 integer (1 to 9 bytes).
fn ltf8(src: &mut &[u8]) -> Option<i64> {
    let (&b0, rest) = src.split_first()?;
    let extra = b0.leading_ones() as usize;
    let b = rest.get(..extra)?;
    let mut v = if extra >= 7 {
        0
    } else {
        u64::from(b0 & (0x7f >> extra))
    };
    for &x in b {
        v = (v << 8) | u64::from(x);
    }
    *src = &rest[extra..];
    Some(v as i64)
}

fn u32_le(src: &mut &[u8]) -> Option<u32> {
    let (b, rest) = src.split_first_chunk::<4>()?;
    *src = rest;
    Some(u32::from_le_bytes(*b))
}

fn crc32(b: &[u8]) -> u32 {
    let mut c = flate2::Crc::new();
    c.update(b);
    c.sum()
}

/// Append up to `n` bytes of `r` to `buf`; returns how many came (fewer only at the end).
fn take_into(r: &mut dyn Read, buf: &mut Vec<u8>, n: usize) -> io::Result<usize> {
    let before = buf.len();
    r.take(n as u64).read_to_end(buf)?;
    Ok(buf.len() - before)
}

// ---- framing ----

/// A container header (CRAM 3).
#[derive(Debug, Clone)]
pub(crate) struct ContainerHeader {
    /// Bytes of data after the header.
    pub len: usize,
    pub ref_id: i32,
    pub start: i32,
    pub blocks: i32,
    pub landmarks: Vec<usize>,
    /// Bytes of the header itself.
    pub header_len: usize,
    pub crc: u32,
}

impl ContainerHeader {
    /// The end-of-file container (as noodles and htslib recognise it).
    fn is_eof(&self) -> bool {
        self.len == 15
            && self.ref_id == -1
            && self.start == 4_542_278
            && self.blocks == 1
            && self.crc == 0x4fd9_bd05
    }
}

/// One data container: its header, and its bytes (header and data) as read.
pub(crate) struct RawContainer {
    /// Counted from 1 over the data containers of the file.
    pub number: u64,
    /// Where the container starts in the file.
    pub offset: u64,
    pub header: ContainerHeader,
    pub bytes: Vec<u8>,
}

impl RawContainer {
    pub fn data(&self) -> &[u8] {
        &self.bytes[self.header.header_len..]
    }

    fn place(&self) -> String {
        format!(
            "container {} (at byte {})",
            thousands(self.number),
            thousands(self.offset)
        )
    }
}

/// Why a container header could not be read.
enum HeaderProblem {
    /// The stream ended before the first byte.
    End,
    /// The stream ended inside the header.
    Short,
    Io(io::Error),
    Bad(String),
}

/// Read a container header from `r`, appending its bytes to `raw`.
fn read_container_header(
    r: &mut dyn Read,
    raw: &mut Vec<u8>,
) -> Result<ContainerHeader, HeaderProblem> {
    let start = raw.len();
    let need = |r: &mut dyn Read, raw: &mut Vec<u8>, n: usize| -> Result<(), HeaderProblem> {
        match take_into(r, raw, n) {
            Ok(k) if k == n => Ok(()),
            Ok(_) if raw.len() == start => Err(HeaderProblem::End),
            Ok(_) => Err(HeaderProblem::Short),
            Err(e) => Err(HeaderProblem::Io(e)),
        }
    };
    // An ITF8 or LTF8 integer: its first byte says how many follow.
    let int = |r: &mut dyn Read, raw: &mut Vec<u8>, long: bool| -> Result<i64, HeaderProblem> {
        let at = raw.len();
        need(r, raw, 1)?;
        let b0 = raw[at];
        let extra = if long {
            b0.leading_ones() as usize
        } else {
            match b0 {
                0x00..=0x7f => 0,
                0x80..=0xbf => 1,
                0xc0..=0xdf => 2,
                0xe0..=0xef => 3,
                _ => 4,
            }
        };
        need(r, raw, extra)?;
        let mut s = &raw[at..];
        let v = if long {
            ltf8(&mut s)
        } else {
            itf8(&mut s).map(i64::from)
        };
        v.ok_or(HeaderProblem::Short)
    };
    need(r, raw, 4)?;
    let len = i32::from_le_bytes([raw[start], raw[start + 1], raw[start + 2], raw[start + 3]]);
    let ref_id = int(r, raw, false)? as i32;
    let pos = int(r, raw, false)? as i32;
    let _span = int(r, raw, false)?;
    let _records = int(r, raw, false)?;
    let _counter = int(r, raw, true)?;
    let _bases = int(r, raw, true)?;
    let blocks = int(r, raw, false)? as i32;
    let n = int(r, raw, false)?;
    let n = usize::try_from(n)
        .ok()
        .filter(|&n| n <= MAX_LANDMARKS)
        .ok_or_else(|| HeaderProblem::Bad(format!("a container claims {n} slices")))?;
    let mut landmarks = Vec::with_capacity(n.min(1024));
    for _ in 0..n {
        let l = int(r, raw, false)?;
        landmarks.push(
            usize::try_from(l).map_err(|_| {
                HeaderProblem::Bad(format!("a container has a slice at offset {l}"))
            })?,
        );
    }
    let crc = crc32(&raw[start..]);
    need(r, raw, 4)?;
    let at = raw.len() - 4;
    let stored = u32::from_le_bytes([raw[at], raw[at + 1], raw[at + 2], raw[at + 3]]);
    if stored != crc {
        return Err(HeaderProblem::Bad(
            "a container header's checksum does not match its bytes".into(),
        ));
    }
    let len = usize::try_from(len)
        .map_err(|_| HeaderProblem::Bad(format!("a container claims a length of {len} bytes")))?;
    let h = ContainerHeader {
        len,
        ref_id,
        start: pos,
        blocks,
        landmarks,
        header_len: raw.len() - start,
        crc: stored,
    };
    if h.is_eof() {
        return Ok(h);
    }
    if h.landmarks.iter().any(|&l| l >= h.len) || h.landmarks.windows(2).any(|w| w[0] >= w[1]) {
        return Err(HeaderProblem::Bad(
            "a container's slice offsets are not inside it in order".into(),
        ));
    }
    Ok(h)
}

/// One block of a container or slice, as stored.
struct Block<'a> {
    method: u8,
    content_type: u8,
    raw_size: usize,
    data: &'a [u8],
}

/// Read the block at the front of `src` (CRAM 3: with its CRC32).
fn read_block<'a>(src: &mut &'a [u8]) -> Result<Block<'a>, String> {
    let all = *src;
    let short = || "a block is cut short".to_string();
    let (&method, rest) = src.split_first().ok_or_else(short)?;
    let (&content_type, mut rest) = rest.split_first().ok_or_else(short)?;
    let _content_id = itf8(&mut rest).ok_or_else(short)?;
    let size = itf8(&mut rest).ok_or_else(short)?;
    let raw_size = itf8(&mut rest).ok_or_else(short)?;
    let size = usize::try_from(size).map_err(|_| format!("a block claims {size} bytes"))?;
    let raw_size =
        usize::try_from(raw_size).map_err(|_| format!("a block claims {raw_size} bytes"))?;
    let data = rest.get(..size).ok_or_else(short)?;
    let mut rest = &rest[size..];
    let used = all.len() - rest.len();
    let stored = u32_le(&mut rest).ok_or_else(short)?;
    if crc32(&all[..used]) != stored {
        return Err("a block's checksum does not match its bytes".into());
    }
    *src = rest;
    Ok(Block {
        method,
        content_type,
        raw_size,
        data,
    })
}

/// The contents of a raw or gzip block (the only methods header blocks use).
fn block_contents<'a>(b: &Block<'a>, what: &str) -> Result<Cow<'a, [u8]>, String> {
    match b.method {
        0 => Ok(Cow::Borrowed(b.data)),
        1 => {
            let mut out = Vec::new();
            flate2::read::MultiGzDecoder::new(b.data)
                .take(b.raw_size as u64)
                .read_to_end(&mut out)
                .map_err(|e| format!("the {what} block cannot be decompressed: {e}"))?;
            Ok(Cow::Owned(out))
        }
        m => Err(format!(
            "the {what} block uses compression method {m}, which gloci does not read there"
        )),
    }
}

/// The start of a CRAM file: its version and SAM header text.
pub(crate) struct CramHead {
    pub text: Vec<u8>,
    /// Bytes of the file definition and header container.
    pub len: u64,
}

fn detect_err(name: &str, msg: String) -> IoError {
    IoError::Detect {
        path: name.to_string(),
        msg,
    }
}

fn corrupt(name: &str, what: &str) -> IoError {
    detect_err(
        name,
        format!("{name}: {what}; the file is not valid CRAM, or it is corrupt"),
    )
}

/// A container whose encodings gloci (noodles) cannot decode.
fn unsupported_codec(name: &str, c: &RawContainer, what: &str) -> IoError {
    detect_err(
        name,
        format!(
            "{name}: {} uses the CRAM encoding {what}, which gloci does not decode (samtools \
does not write it); rewrite the file with samtools: samtools view -C -T ref.fa -o new.cram {name}",
            c.place()
        ),
    )
}

fn truncated(name: &str, after: &str) -> IoError {
    detect_err(
        name,
        format!(
            "{name}: the CRAM end-of-file container is missing ({after}); the file is probably \
truncated"
        ),
    )
}

/// Read the file definition and the header container of `name` from `r`.
pub(crate) fn read_head(r: &mut dyn Read, name: &str) -> Result<CramHead, IoError> {
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let mut def = Vec::new();
    take_into(r, &mut def, 26).map_err(read_err)?;
    if def.len() < 26 || !def.starts_with(b"CRAM") {
        return Err(corrupt(
            name,
            "the CRAM file definition is missing or cut short",
        ));
    }
    let version = (def[4], def[5]);
    if version.0 != 3 || version.1 > 1 {
        return Err(detect_err(
            name,
            format!(
                "{name} is CRAM {}.{}, which gloci does not read; it reads CRAM 3.0 and 3.1; \
rewrite it with: samtools view -C -T ref.fa --output-fmt-option version=3.0 -o x3.cram {name}",
                version.0, version.1
            ),
        ));
    }
    let mut raw = Vec::new();
    let h = match read_container_header(r, &mut raw) {
        Ok(h) => h,
        Err(HeaderProblem::End | HeaderProblem::Short) => {
            return Err(truncated(name, "the file ends inside its header"));
        }
        Err(HeaderProblem::Io(e)) => return Err(read_err(e)),
        Err(HeaderProblem::Bad(m)) => return Err(corrupt(name, &format!("in the header, {m}"))),
    };
    let mut data = Vec::new();
    if take_into(r, &mut data, h.len).map_err(read_err)? < h.len {
        return Err(truncated(name, "the file ends inside its header"));
    }
    let mut src = &data[..];
    let block = read_block(&mut src).map_err(|m| corrupt(name, &format!("in the header, {m}")))?;
    if block.content_type != 0 {
        return Err(corrupt(
            name,
            "the CRAM header container does not hold the SAM header",
        ));
    }
    let contents = block_contents(&block, "header").map_err(|m| corrupt(name, &m))?;
    let mut c = &contents[..];
    let text = u32_le(&mut c)
        .and_then(|n| c.get(..n as usize))
        .ok_or_else(|| {
            corrupt(
                name,
                "the SAM header in the CRAM header container is cut short",
            )
        })?;
    let mut text = text.to_vec();
    if let Some(nul) = memchr::memchr(0, &text) {
        text.truncate(nul);
    }
    Ok(CramHead {
        text,
        len: 26 + (raw.len() + data.len()) as u64,
    })
}

/// The data containers of a CRAM file, in order, after its header.
pub(crate) struct CramFrames {
    name: String,
    src: Box<dyn BufRead + Send>,
    /// Bytes read so far.
    pos: u64,
    containers: u64,
    done: bool,
}

impl CramFrames {
    /// Frames of `src`, whose header (`head.len` bytes) has been read.
    pub fn new(name: &str, src: Box<dyn BufRead + Send>, head: &CramHead) -> CramFrames {
        CramFrames {
            name: name.to_string(),
            src,
            pos: head.len,
            containers: 0,
            done: false,
        }
    }

    /// The next data container, or `None` after the end-of-file container. A stream that ends
    /// without one is an error.
    pub fn next(&mut self) -> Result<Option<RawContainer>, IoError> {
        if self.done {
            return Ok(None);
        }
        let name = self.name.clone();
        let read_err = |source| IoError::Read {
            path: name.clone(),
            source,
        };
        let after = |n: u64| {
            if n == 0 {
                "the file ends after its header".to_string()
            } else {
                format!("the file ends after container {}", thousands(n))
            }
        };
        let offset = self.pos;
        let mut bytes = Vec::new();
        let h = match read_container_header(&mut self.src, &mut bytes) {
            Ok(h) => h,
            Err(HeaderProblem::End) => {
                self.done = true;
                return Err(truncated(&self.name, &after(self.containers)));
            }
            Err(HeaderProblem::Short) => {
                self.done = true;
                return Err(truncated(
                    &self.name,
                    &format!(
                        "it ends inside container {} at byte {}",
                        thousands(self.containers + 1),
                        thousands(offset)
                    ),
                ));
            }
            Err(HeaderProblem::Io(e)) => return Err(read_err(e)),
            Err(HeaderProblem::Bad(m)) => {
                self.done = true;
                return Err(corrupt(
                    &self.name,
                    &format!(
                        "at container {} (byte {}), {m}",
                        thousands(self.containers + 1),
                        thousands(offset)
                    ),
                ));
            }
        };
        let got = take_into(&mut self.src, &mut bytes, h.len).map_err(read_err)?;
        self.pos += bytes.len() as u64;
        if got < h.len {
            self.done = true;
            return Err(truncated(
                &self.name,
                &format!(
                    "it ends inside container {} at byte {}",
                    thousands(self.containers + 1),
                    thousands(offset)
                ),
            ));
        }
        if h.is_eof() {
            self.done = true;
            return Ok(None);
        }
        self.containers += 1;
        Ok(Some(RawContainer {
            number: self.containers,
            offset,
            header: h,
            bytes,
        }))
    }
}

/// A slice header, checked against the CRAM header.
#[derive(Debug, Clone)]
pub(crate) struct SliceHead {
    /// The reference id: -1 for unplaced reads, -2 for several references.
    pub ref_id: i32,
    pub start: i32,
    pub span: i32,
    /// The external block holding the slice's reference bases, or -1.
    pub embedded: i32,
    /// Where the slice starts in the container's data, and its bytes.
    pub landmark: usize,
    pub size: usize,
}

impl SliceHead {
    /// The last reference base the slice covers (1-based), for a slice on one reference.
    fn end(&self) -> i64 {
        i64::from(self.start) + i64::from(self.span) - 1
    }
}

/// The slice headers of container `c` of `name`, whose header lists `n_refs` references.
pub(crate) fn slice_heads(
    c: &RawContainer,
    name: &str,
    n_refs: usize,
) -> Result<Vec<SliceHead>, IoError> {
    let data = c.data();
    let bad = |m: String| corrupt(name, &format!("in {}, {m}", c.place()));
    let landmarks = &c.header.landmarks;
    let mut out = Vec::with_capacity(landmarks.len());
    for (i, &at) in landmarks.iter().enumerate() {
        let end = landmarks.get(i + 1).copied().unwrap_or(data.len());
        let mut src = &data[at..end];
        let block = read_block(&mut src).map_err(&bad)?;
        if block.content_type != 2 {
            return Err(bad(format!(
                "slice {} does not start with a slice header",
                i + 1
            )));
        }
        let contents = block_contents(&block, "slice header").map_err(&bad)?;
        let mut s = &contents[..];
        let short = || bad(format!("the header of slice {} is cut short", i + 1));
        let ref_id = itf8(&mut s).ok_or_else(short)?;
        let start = itf8(&mut s).ok_or_else(short)?;
        let span = itf8(&mut s).ok_or_else(short)?;
        let records = itf8(&mut s).ok_or_else(short)?;
        let _counter = ltf8(&mut s).ok_or_else(short)?;
        let blocks = itf8(&mut s).ok_or_else(short)?;
        let ids = itf8(&mut s).ok_or_else(short)?;
        let ids = usize::try_from(ids).map_err(|_| short())?;
        for _ in 0..ids {
            itf8(&mut s).ok_or_else(short)?;
        }
        let embedded = itf8(&mut s).ok_or_else(short)?;
        if s.len() < 16 {
            return Err(short());
        }
        usize::try_from(records)
            .ok()
            .filter(|&n| n <= MAX_SLICE_RECORDS)
            .ok_or_else(|| {
                bad(format!(
                    "slice {} claims {records} records; gloci reads at most {} in one slice",
                    i + 1,
                    thousands(MAX_SLICE_RECORDS as u64)
                ))
            })?;
        if blocks < 1 {
            return Err(bad(format!("slice {} has no data blocks", i + 1)));
        }
        if ref_id < -2 || ref_id >= 0 && ref_id as usize >= n_refs {
            return Err(bad(format!(
                "slice {} names reference id {ref_id}, but the header has {} references",
                i + 1,
                n_refs
            )));
        }
        if ref_id >= 0 && (start < 1 || span < 1) {
            return Err(bad(format!(
                "slice {} has start {start} and span {span} on its reference",
                i + 1
            )));
        }
        out.push(SliceHead {
            ref_id,
            start,
            span,
            embedded,
            landmark: at,
            size: end - at,
        });
    }
    Ok(out)
}

/// What gloci reads of a container's compression header before noodles does.
struct CompressionInfo {
    /// The preservation map's `RR` (true unless the map says otherwise): the container needs
    /// reference bases. A header gloci cannot read is taken to need them.
    required: bool,
    /// An encoding noodles cannot decode (it would panic), named.
    unsupported: Option<String>,
}

/// The preservation map's `RR`, from its contents `s`.
fn preservation_rr(mut s: &[u8]) -> bool {
    let Some(n) = itf8(&mut s) else {
        return true;
    };
    for _ in 0..n.max(0) {
        let Some((key, rest)) = s.split_first_chunk::<2>() else {
            return true;
        };
        s = rest;
        match key {
            b"RN" | b"AP" => match s.split_first() {
                Some((_, r)) => s = r,
                None => return true,
            },
            b"RR" => return s.first().is_none_or(|&b| b != 0),
            b"SM" => match s.get(5..) {
                Some(r) => s = r,
                None => return true,
            },
            b"TD" => {
                let Some(len) = itf8(&mut s).and_then(|l| usize::try_from(l).ok()) else {
                    return true;
                };
                match s.get(len..) {
                    Some(r) => s = r,
                    None => return true,
                }
            }
            _ => return true,
        }
    }
    true
}

/// Check the encoding at the front of `s` (and those nested in it): noodles 0.102 panics on the
/// Golomb, Subexp and Golomb-Rice codecs (`todo!`) and on an empty Huffman table. `None` when
/// the encoding cannot be read (noodles reports that itself).
fn encoding_problem(s: &mut &[u8], depth: u8) -> Option<Option<String>> {
    let id = itf8(s)?;
    let len = usize::try_from(itf8(s)?).ok()?;
    let mut params = s.get(..len)?;
    *s = &s[len..];
    Some(match id {
        2 => Some("Golomb".to_string()),
        7 => Some("Subexp".to_string()),
        8 => Some("Golomb-Rice".to_string()),
        3 => {
            let n = itf8(&mut params)?;
            for _ in 0..n.max(0) {
                itf8(&mut params)?;
            }
            let m = itf8(&mut params)?;
            (n < 1 || n != m).then(|| format!("Huffman with {n} symbols and {m} lengths"))
        }
        4 if depth < 8 => {
            let first = encoding_problem(&mut params, depth + 1)?;
            match first {
                Some(p) => Some(p),
                None => encoding_problem(&mut params, depth + 1)?,
            }
        }
        _ => None,
    })
}

/// The first unsupported encoding of a data series or tag encoding map (contents `s`); `wide`
/// keys are ITF8 (tags), else two bytes (data series).
fn map_problem(mut s: &[u8], wide: bool) -> Option<String> {
    let n = itf8(&mut s)?;
    for _ in 0..n.max(0) {
        if wide {
            itf8(&mut s)?;
        } else {
            s = s.get(2..)?;
        }
        if let Some(p) = encoding_problem(&mut s, 0)? {
            return Some(p);
        }
    }
    None
}

/// Read the compression header of container `c` (see [`CompressionInfo`]).
fn compression(c: &RawContainer) -> CompressionInfo {
    let unknown = CompressionInfo {
        required: true,
        unsupported: None,
    };
    let data = c.data();
    let end = c.header.landmarks.first().copied().unwrap_or(data.len());
    let mut src = &data[..end];
    let Ok(block) = read_block(&mut src) else {
        return unknown;
    };
    let Ok(contents) = block_contents(&block, "compression header") else {
        return unknown;
    };
    // Three maps, each its byte size, then its contents.
    let mut s = &contents[..];
    let mut maps = Vec::with_capacity(3);
    for _ in 0..3 {
        let Some(size) = itf8(&mut s).and_then(|n| usize::try_from(n).ok()) else {
            break;
        };
        let Some(m) = s.get(..size) else {
            break;
        };
        maps.push(m);
        s = &s[size..];
    }
    let Some(preservation) = maps.first() else {
        return unknown;
    };
    CompressionInfo {
        required: preservation_rr(preservation),
        unsupported: maps
            .get(1)
            .and_then(|m| map_problem(m, false))
            .or_else(|| maps.get(2).and_then(|m| map_problem(m, true))),
    }
}

// ---- the CRAM header and the reference ----

/// What every decoder of one CRAM shares.
pub(crate) struct CramInfo {
    /// The input's name, for messages.
    name: String,
    /// The header for noodles: the references and read groups, in file order.
    sam: sam::Header,
    /// The header as read (references named as in the file).
    aln: AlnHeader,
    /// Each reference's `M5` and `UR`, for messages.
    m5: Vec<Option<String>>,
    ur: Vec<Option<String>>,
    /// The `--reference` FASTA: a handle each decoder reopens, and its name.
    fasta: Option<(Mutex<IndexedFasta>, String)>,
    /// Whether the FASTA has each reference.
    in_fasta: Vec<bool>,
}

fn field<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    line.split(|&b| b == b'\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key))
}

impl CramInfo {
    /// Build from the header `h` of the CRAM `name`, opening and checking `reference` (D109).
    fn new(name: &str, h: AlnHeader, reference: Option<&Path>) -> Result<CramInfo, IoError> {
        let n = h.refs.len();
        let mut m5 = Vec::with_capacity(n);
        let mut ur = Vec::with_capacity(n);
        let mut groups = Vec::new();
        for line in h.text.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            let text = |v: &[u8]| String::from_utf8_lossy(v).into_owned();
            if line.starts_with(b"@SQ\t") {
                m5.push(field(line, b"M5:").map(text));
                ur.push(field(line, b"UR:").map(text));
            } else if line.starts_with(b"@RG\t")
                && let Some(id) = field(line, b"ID:")
            {
                groups.push(id.to_vec());
            }
        }
        m5.resize(n, None);
        ur.resize(n, None);
        let mut b = sam::Header::builder();
        for (r, len) in &h.refs {
            let len = usize::try_from(*len)
                .ok()
                .and_then(NonZeroUsize::new)
                .ok_or_else(|| {
                    detect_err(
                        name,
                        format!(
                            "{name}: the @SQ line of {} has LN:{len}; a CRAM reference needs a \
length of at least 1",
                            String::from_utf8_lossy(r)
                        ),
                    )
                })?;
            b =
                b.add_reference_sequence(
                    r.clone(),
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(len),
                );
        }
        for g in groups {
            b = b.add_read_group(g, sam::header::record::value::Map::default());
        }
        let mut in_fasta = vec![false; n];
        let fasta = match reference {
            None => None,
            Some(p) => {
                let f = IndexedFasta::open(p)?;
                let shown = p.display().to_string();
                for (i, (r, len)) in h.refs.iter().enumerate() {
                    let Some(flen) = f.len_of(r) else {
                        continue;
                    };
                    if flen != *len {
                        let made = ur[i]
                            .as_deref()
                            .map(|u| format!(" (its header names {u})"))
                            .unwrap_or_default();
                        return Err(detect_err(
                            name,
                            format!(
                                "{name} was not made with {shown}: {} is {} bases long in its \
header and {} in {shown}; pass the FASTA it was made with{made}: --reference FASTA",
                                String::from_utf8_lossy(r),
                                thousands(*len),
                                thousands(flen)
                            ),
                        ));
                    }
                    in_fasta[i] = true;
                }
                Some((Mutex::new(f), shown))
            }
        };
        Ok(CramInfo {
            name: name.to_string(),
            sam: b.build(),
            aln: h,
            m5,
            ur,
            fasta,
            in_fasta,
        })
    }

    fn chrom(&self, rid: usize) -> String {
        String::from_utf8_lossy(&self.aln.refs[rid].0).into_owned()
    }

    /// The D109 error for a reference the decoder needs and does not have.
    fn needs_reference(&self, name: &[u8]) -> IoError {
        let chrom = String::from_utf8_lossy(name);
        let rid = self.aln.ref_id(name);
        let mut ids = Vec::new();
        if let Some(m) = rid.and_then(|i| self.m5[i].as_deref()) {
            ids.push(format!("M5 {m}"));
        }
        if let Some(u) = rid.and_then(|i| self.ur[i].as_deref()) {
            ids.push(format!("UR {u}"));
        }
        let ids = if ids.is_empty() {
            String::new()
        } else {
            format!(" ({})", ids.join(", "))
        };
        let cram = &self.name;
        let msg = match &self.fasta {
            None => format!(
                "{cram} needs the reference sequence of {chrom}{ids} to decode its reads; pass \
it with --reference FASTA"
            ),
            Some((_, fa)) => format!(
                "{cram} needs the reference sequence of {chrom}{ids} to decode its reads, and \
{fa} has no sequence {chrom}; pass the FASTA it was made with: --reference FASTA"
            ),
        };
        detect_err(cram, msg)
    }

    /// The D109 error for a slice whose reference bases do not have the checksum it stores.
    fn mismatch(&self, rid: usize, start: i64, end: i64) -> IoError {
        let fa = self
            .fasta
            .as_ref()
            .map_or("the reference", |(_, n)| n.as_str());
        let chrom = self.chrom(rid);
        detect_err(
            &self.name,
            format!(
                "the sequence {chrom} in {fa} does not match the one {} was made with \
({chrom}:{}-{}); pass the FASTA it was made with: --reference FASTA",
                self.name,
                thousands(start.max(0) as u64),
                thousands(end.max(0) as u64)
            ),
        )
    }

    pub(crate) fn ref_count(&self) -> usize {
        self.aln.refs.len()
    }

    pub(crate) fn name(&self) -> &str {
        &self.name
    }
}

/// What a decoder and its [`RefAdapter`] share.
#[derive(Default)]
pub(crate) struct RefState {
    /// `N` bases added after each sequence, for reads that run past the end of their reference
    /// (htslib reads them against `N`).
    pad: AtomicUsize,
    /// The `RR` of the container being decoded: without it no bases are fetched.
    required: AtomicBool,
    /// Without `RR`, the length of the `N` placeholder given for each name. noodles still asks
    /// for the sequence of each mapped read of a slice of several references, and slices it at
    /// the read's positions though it reads no base, so the placeholder must reach them.
    placeholder: AtomicUsize,
    /// The adapter's last failure, with gloci's message (noodles passes on only an `io::Error`).
    failure: Mutex<Option<IoError>>,
    /// The names fetched since the repository was last cleared.
    cached: Mutex<Vec<Vec<u8>>>,
    /// The slice being decoded may add sequences to the repository. Without it, a fetch while
    /// the repository holds another sequence fails with `evict` set, and the decoder clears the
    /// repository and decodes the slice again: a decoder holds only what one slice needs.
    grow: AtomicBool,
    evict: AtomicBool,
    /// Sequences read from the FASTA.
    fetches: AtomicUsize,
    /// The most sequences the decoder's repository held at once.
    max_live: AtomicUsize,
}

/// The repository adapter over the `--reference` FASTA (D109). It never answers `None`.
pub(crate) struct RefAdapter {
    info: Arc<CramInfo>,
    fasta: Option<Mutex<IndexedFasta>>,
    state: Arc<RefState>,
}

const ADAPTER_FAILED: &str = "the reference adapter failed";
const EVICT: &str = "the repository holds another sequence";

impl RefAdapter {
    fn fetch(&mut self, name: &[u8]) -> Result<fasta::Record, IoError> {
        let record = |seq: Vec<u8>| {
            fasta::Record::new(
                fasta::record::Definition::new(name.to_vec(), None),
                fasta::record::Sequence::from(seq),
            )
        };
        if !self.state.required.load(Ordering::Relaxed) {
            // A container stored without the reference reads no base of it.
            if let Ok(mut c) = self.state.cached.lock() {
                c.push(name.to_vec());
            }
            return Ok(record(vec![
                b'N';
                self.state.placeholder.load(Ordering::Relaxed)
            ]));
        }
        let Some(f) = self.fasta.as_mut().and_then(|m| m.get_mut().ok()) else {
            return Err(self.info.needs_reference(name));
        };
        let Some(len) = f.len_of(name) else {
            return Err(self.info.needs_reference(name));
        };
        let mut cached = self
            .state
            .cached
            .lock()
            .map_err(|_| IoError::Internal("the reference cache lock is poisoned".into()))?;
        if !cached.is_empty() && !self.state.grow.load(Ordering::Relaxed) {
            self.state.evict.store(true, Ordering::Relaxed);
            return Err(IoError::Internal(EVICT.into()));
        }
        cached.push(name.to_vec());
        drop(cached);
        let pad = self.state.pad.load(Ordering::Relaxed);
        let mut seq = Vec::with_capacity(usize::try_from(len).unwrap_or(0).saturating_add(pad));
        f.fetch(
            name,
            Interval { start: 0, end: len },
            &mut |chunk: &[u8]| {
                seq.extend(chunk.iter().map(u8::to_ascii_uppercase));
                Ok(())
            },
        )?;
        seq.resize(seq.len() + pad, b'N');
        self.state.fetches.fetch_add(1, Ordering::Relaxed);
        Ok(record(seq))
    }
}

impl fasta::repository::Adapter for RefAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        Some(self.fetch(name).map_err(|e| {
            if !self.state.evict.load(Ordering::Relaxed)
                && let Ok(mut slot) = self.state.failure.lock()
                && slot.is_none()
            {
                *slot = Some(e);
            }
            io::Error::other(ADAPTER_FAILED)
        }))
    }
}

// ---- decoding ----

/// Where a slice's reads find their reference bases.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Bases {
    /// The repository (`--reference`).
    Repository,
    /// The slice's embedded reference block.
    Embedded,
    /// Every base is stored, but noodles slices a placeholder sequence at each mapped read (a
    /// slice of several references without `RR`).
    Placeholder,
    /// Nowhere: every base is stored.
    Stored,
}

/// What the reads of a decoded slice need beyond what they were decoded with.
enum Extent {
    Fits,
    /// The repository's sequences need this many `N` after their end.
    Pad(usize),
    /// The embedded bases need to run to this reference position (1-based).
    EmbeddedTo(i64),
    /// The placeholder sequences need this length.
    Placeholder(usize),
}

/// Decodes containers to BAM records. Each worker of the pool has one.
pub(crate) struct Decoder {
    info: Arc<CramInfo>,
    repo: fasta::Repository,
    state: Arc<RefState>,
    /// The reference id of the last slice decoded.
    last_ref: Option<i32>,
    /// The `RR` of the last container decoded.
    last_required: Option<bool>,
    bam: noodles_bam::io::Writer<Vec<u8>>,
}

impl Decoder {
    fn new(info: Arc<CramInfo>) -> Result<Decoder, IoError> {
        let fasta = match &info.fasta {
            None => None,
            Some((f, _)) => {
                let f = f
                    .lock()
                    .map_err(|_| IoError::Internal("the reference lock is poisoned".into()))?;
                Some(Mutex::new(f.reopen()?))
            }
        };
        let state = Arc::new(RefState::default());
        let adapter = RefAdapter {
            info: info.clone(),
            fasta,
            state: state.clone(),
        };
        Ok(Decoder {
            info,
            repo: fasta::Repository::new(adapter),
            state,
            last_ref: None,
            last_required: None,
            bam: noodles_bam::io::Writer::from(Vec::new()),
        })
    }

    fn corrupt(&self, c: &RawContainer, e: &dyn std::fmt::Display) -> IoError {
        let name = &self.info.name;
        detect_err(
            name,
            format!(
                "{name}: {} cannot be decoded: {e}; the file is corrupt",
                c.place()
            ),
        )
    }

    /// The gloci error for a noodles failure while decoding `head` of `c`.
    fn decode_error(&self, c: &RawContainer, head: &SliceHead, e: io::Error) -> IoError {
        if let Some(f) = self.state.failure.lock().ok().and_then(|mut s| s.take()) {
            return f;
        }
        if e.kind() == io::ErrorKind::InvalidData
            && e.to_string()
                .starts_with("reference sequence checksum mismatch")
            && head.ref_id >= 0
        {
            return self
                .info
                .mismatch(head.ref_id as usize, i64::from(head.start), head.end());
        }
        self.corrupt(c, &e)
    }

    /// Empty the repository (and with `reset_pad`, forget the padding).
    fn clear(&mut self, reset_pad: bool) {
        self.repo.clear();
        if let Ok(mut c) = self.state.cached.lock() {
            c.clear();
        }
        if reset_pad {
            self.state.pad.store(0, Ordering::Relaxed);
            self.state.placeholder.store(0, Ordering::Relaxed);
        }
    }

    /// Where the reads of `head` find their bases, after checking that they can (D109).
    fn bases(
        &mut self,
        c: &RawContainer,
        head: &SliceHead,
        required: bool,
    ) -> Result<Bases, IoError> {
        // One sequence at a time (D109): a slice on another reference starts afresh. A slice of
        // several references keeps what the repository holds; the adapter asks for a clear if
        // the slice needs anything else.
        if head.ref_id >= 0 && self.last_ref != Some(head.ref_id) {
            let name = &self.info.aln.refs[head.ref_id as usize].0;
            let keep = self
                .state
                .cached
                .lock()
                .is_ok_and(|c| c.iter().all(|n| n == name));
            if !keep {
                self.clear(true);
            }
        }
        self.last_ref = Some(head.ref_id);
        let bases = match head.ref_id {
            -1 => Bases::Stored,
            _ if required => Bases::Repository,
            -2 => Bases::Placeholder,
            _ if head.embedded != -1 => Bases::Embedded,
            _ => Bases::Stored,
        };
        if head.ref_id >= 0 && required {
            let rid = head.ref_id as usize;
            let (name, len) = &self.info.aln.refs[rid];
            if !self.info.in_fasta[rid] {
                return Err(self.info.needs_reference(name));
            }
            if head.end() > *len as i64 {
                return Err(self.corrupt(
                    c,
                    &format!(
                        "a slice of {} runs to {}, past its end at {}",
                        self.info.chrom(rid),
                        thousands(head.end() as u64),
                        thousands(*len)
                    ),
                ));
            }
        }
        Ok(bases)
    }

    /// Append the BAM records (each with its `block_size`) of container `c` to `out`.
    pub fn decode(&mut self, c: &RawContainer, out: &mut Vec<u8>) -> Result<(), IoError> {
        let heads = slice_heads(c, &self.info.name, self.info.aln.refs.len())?;
        let CompressionInfo {
            required,
            unsupported,
        } = compression(c);
        if let Some(what) = unsupported {
            return Err(unsupported_codec(&self.info.name, c, &what));
        }
        self.state.required.store(required, Ordering::Relaxed);
        if self.last_required != Some(required) {
            // Placeholders never serve a container that reads the reference, nor the reverse.
            self.clear(true);
            self.last_ref = None;
            self.last_required = Some(required);
        }
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&c.bytes[..])
            .read_container(&mut container)
            .map_err(|e| self.corrupt(c, &e))?;
        let ch = container
            .compression_header()
            .map_err(|e| self.corrupt(c, &e))?;
        // The records borrow the header; the decoder itself stays free to transcode them.
        let info = self.info.clone();
        for (head, slice) in heads.iter().zip(container.slices()) {
            let slice = slice.map_err(|e| self.corrupt(c, &e))?;
            let bases = self.bases(c, head, required)?;
            let (core, ext) = slice.decode_blocks().map_err(|e| self.corrupt(c, &e))?;
            if bases == Bases::Embedded && !ext.iter().any(|(id, _)| *id == head.embedded) {
                // noodles finds the embedded block by its id, and `expect`s it.
                return Err(self.corrupt(
                    c,
                    &"a slice names an embedded reference block it does not have",
                ));
            }
            // A slice that needs a sequence the repository lacks, while it holds another, is
            // decoded again after a clear; reads past the end of their reference bases are
            // decoded again against `N`.
            let mut padded: Option<Vec<(i32, Cow<'_, [u8]>)>> = None;
            let mut done = false;
            self.state.grow.store(false, Ordering::Relaxed);
            self.state.evict.store(false, Ordering::Relaxed);
            for _ in 0..4 {
                let srcs = padded.as_deref().unwrap_or(&ext);
                let records = match slice.records(self.repo.clone(), &info.sam, &ch, &core, srcs) {
                    Ok(r) => r,
                    Err(_) if self.state.evict.swap(false, Ordering::Relaxed) => {
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        continue;
                    }
                    Err(e) => return Err(self.decode_error(c, head, e)),
                };
                self.state
                    .max_live
                    .fetch_max(self.repo.len(), Ordering::Relaxed);
                match self.extent(c, head, bases, &records, srcs)? {
                    Extent::Fits => {
                        self.transcode(c, &records, out)?;
                        done = true;
                        break;
                    }
                    Extent::Pad(pad) => {
                        drop(records);
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        self.state.pad.store(pad, Ordering::Relaxed);
                    }
                    Extent::Placeholder(len) => {
                        drop(records);
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        self.state.placeholder.store(len, Ordering::Relaxed);
                    }
                    Extent::EmbeddedTo(end) => {
                        drop(records);
                        let want = usize::try_from(end - i64::from(head.start) + 1).unwrap_or(0);
                        let mut v = ext.to_vec();
                        for (id, src) in &mut v {
                            if *id == head.embedded && src.len() < want {
                                let mut b = src.to_vec();
                                b.resize(want, b'N');
                                *src = Cow::Owned(b);
                            }
                        }
                        padded = Some(v);
                    }
                }
            }
            if !done {
                return Err(IoError::Internal(format!(
                    "{}: the reads of {} still run past their reference bases after padding",
                    self.info.name,
                    c.place()
                )));
            }
        }
        Ok(())
    }

    /// Whether every mapped read of `records` lies within the reference bases it decodes with.
    fn extent(
        &self,
        c: &RawContainer,
        head: &SliceHead,
        bases: Bases,
        records: &[cram::Record<'_>],
        srcs: &[(i32, Cow<'_, [u8]>)],
    ) -> Result<Extent, IoError> {
        if bases == Bases::Stored {
            return Ok(Extent::Fits);
        }
        let pad = self.state.pad.load(Ordering::Relaxed) as i64;
        let embedded_len = srcs
            .iter()
            .find(|(id, _)| *id == head.embedded)
            .map_or(0, |(_, s)| s.len() as i64);
        let placeholder = self.state.placeholder.load(Ordering::Relaxed) as i64;
        let mut need_pad = 0i64;
        let mut need_end = 0i64;
        let mut need_len = 0i64;
        for r in records {
            let flags = r.flags().map_err(|e| self.corrupt(c, &e))?;
            if flags.is_unmapped() {
                continue;
            }
            let Some(Ok(start)) = r.alignment_start() else {
                return Err(self.corrupt(c, &"a mapped read has no position"));
            };
            let start = usize::from(start) as i64;
            let end = match r.alignment_end() {
                Some(Ok(e)) => usize::from(e) as i64,
                _ => start - 1,
            };
            match bases {
                Bases::Repository => {
                    let len = match r.reference_sequence_id(&self.info.sam) {
                        Some(Ok(id)) if id < self.info.aln.refs.len() => {
                            self.info.aln.refs[id].1 as i64
                        }
                        _ => return Err(self.corrupt(c, &"a mapped read has no reference")),
                    };
                    if end > len + pad {
                        need_pad = need_pad.max(end - len);
                    }
                }
                Bases::Embedded => {
                    let first = i64::from(head.start);
                    if start < first {
                        return Err(self.corrupt(
                            c,
                            &"a read starts before the reference bases embedded in its slice",
                        ));
                    }
                    if end > first + embedded_len - 1 {
                        need_end = need_end.max(end);
                    }
                }
                Bases::Placeholder => {
                    if start.max(end) > placeholder {
                        need_len = need_len.max(start.max(end));
                    }
                }
                Bases::Stored => {}
            }
        }
        Ok(if need_pad > 0 {
            Extent::Pad(need_pad as usize)
        } else if need_end > 0 {
            Extent::EmbeddedTo(need_end)
        } else if need_len > 0 {
            Extent::Placeholder(need_len as usize)
        } else {
            Extent::Fits
        })
    }

    fn transcode(
        &mut self,
        c: &RawContainer,
        records: &[cram::Record<'_>],
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        for r in records {
            self.bam.get_mut().clear();
            if let Err(e) = self.bam.write_alignment_record(&self.info.sam, r) {
                return Err(self.corrupt(c, &format!("a record cannot be written as BAM: {e}")));
            }
            let rec = self.bam.get_mut();
            htslib_fields(&mut rec[4..]);
            out.extend_from_slice(rec);
        }
        Ok(())
    }
}

/// Set the fields htslib derives when it decodes a CRAM record: MAPQ 0 for an unmapped read
/// (CRAM does not store it) and `bin` from the footprint, as the SAM encoder computes it.
fn htslib_fields(r: &mut [u8]) {
    if r.len() < 32 {
        return;
    }
    let pos = i64::from(i32::from_le_bytes([r[4], r[5], r[6], r[7]]));
    let l_name = usize::from(r[8]);
    let n_cigar = usize::from(u16::from_le_bytes([r[12], r[13]]));
    let flag = u16::from_le_bytes([r[14], r[15]]);
    let unmapped = flag & 0x4 != 0;
    if unmapped && r[9] == 255 {
        r[9] = 0;
    }
    let at = 32 + l_name;
    let rlen: i64 = r
        .get(at..at + 4 * n_cigar)
        .unwrap_or(&[])
        .chunks_exact(4)
        .map(|w| u32::from_le_bytes([w[0], w[1], w[2], w[3]]))
        .filter(|v| matches!(v & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|v| i64::from(v >> 4))
        .sum();
    let end = if unmapped { pos + 1 } else { pos + rlen.max(1) };
    r[10..12].copy_from_slice(&reg2bin(pos, end).to_le_bytes());
}

// ---- the pool (D110) ----

type Job<J> = (u64, J);

/// Work done on up to `threads` workers, results released in submission order.
pub(crate) struct Ordered<J: Send + 'static, R: Send + 'static> {
    jobs: Option<mpsc::Sender<Job<J>>>,
    results: mpsc::Receiver<(u64, Result<R, IoError>)>,
    workers: Vec<JoinHandle<()>>,
    pending: BTreeMap<u64, Result<R, IoError>>,
    submitted: u64,
    released: u64,
}

impl<J: Send + 'static, R: Send + 'static> Ordered<J, R> {
    /// One worker thread per state in `states`, each running `work` on the jobs it takes.
    pub fn new<W: Send + 'static>(
        states: Vec<W>,
        work: fn(&mut W, J) -> Result<R, IoError>,
    ) -> Result<Self, IoError> {
        let (job_tx, job_rx) = mpsc::channel::<Job<J>>();
        let job_rx = Arc::new(Mutex::new(job_rx));
        let (res_tx, results) = mpsc::channel();
        let mut workers = Vec::with_capacity(states.len());
        for (i, mut state) in states.into_iter().enumerate() {
            let rx = job_rx.clone();
            let tx = res_tx.clone();
            let h = std::thread::Builder::new()
                .name(format!("gloci-cram-{i}"))
                .spawn(move || {
                    loop {
                        let job = match rx.lock() {
                            Ok(g) => g.recv(),
                            Err(_) => break,
                        };
                        let Ok((seq, j)) = job else { break };
                        if tx.send((seq, work(&mut state, j))).is_err() {
                            break;
                        }
                    }
                })
                .map_err(|e| IoError::Internal(format!("cannot start a decoder thread: {e}")))?;
            workers.push(h);
        }
        Ok(Ordered {
            jobs: Some(job_tx),
            results,
            workers,
            pending: BTreeMap::new(),
            submitted: 0,
            released: 0,
        })
    }

    /// Jobs submitted and not yet released.
    pub fn in_flight(&self) -> u64 {
        self.submitted - self.released
    }

    pub fn submit(&mut self, j: J) {
        let seq = self.submitted;
        self.submitted += 1;
        let sent = self
            .jobs
            .as_ref()
            .is_some_and(|tx| tx.send((seq, j)).is_ok());
        if !sent {
            self.pending.insert(
                seq,
                Err(IoError::Internal("the decoder threads stopped".into())),
            );
        }
    }

    /// Take a slot in the order for a result known already (a read error).
    pub fn submit_result(&mut self, r: Result<R, IoError>) {
        let seq = self.submitted;
        self.submitted += 1;
        self.pending.insert(seq, r);
    }

    /// The next result in order, waiting for it; `None` when nothing is in flight.
    pub fn next(&mut self) -> Option<Result<R, IoError>> {
        if self.released == self.submitted {
            return None;
        }
        loop {
            if let Some(r) = self.pending.remove(&self.released) {
                self.released += 1;
                return Some(r);
            }
            match self.results.recv() {
                Ok((seq, r)) => {
                    self.pending.insert(seq, r);
                }
                Err(_) => {
                    self.released += 1;
                    return Some(Err(IoError::Internal(
                        "a decoder thread stopped unexpectedly".into(),
                    )));
                }
            }
        }
    }
}

impl<J: Send + 'static, R: Send + 'static> Drop for Ordered<J, R> {
    fn drop(&mut self) {
        self.jobs = None;
        for h in self.workers.drain(..) {
            let _ = h.join();
        }
    }
}

fn decode_job(d: &mut Decoder, c: RawContainer) -> Result<Vec<u8>, IoError> {
    let mut out = Vec::new();
    d.decode(&c, &mut out)?;
    Ok(out)
}

enum Mode {
    Inline(Box<Decoder>),
    Pool(Ordered<RawContainer, Vec<u8>>, u64),
}

/// The decoded records of a CRAM as the byte stream of a BAM body: each record with its
/// `block_size`, in file order. A decoding failure is an `io::Error` holding the gloci error
/// ([`unwrap_failure`] takes it out).
struct CramStream {
    frames: CramFrames,
    mode: Mode,
    cur: Vec<u8>,
    at: usize,
    reading_done: bool,
    failed: bool,
    /// Each decoder's reference counters (the bounded-cache test reads them).
    states: Vec<Arc<RefState>>,
}

impl CramStream {
    fn new(
        frames: CramFrames,
        info: Arc<CramInfo>,
        threads: NonZeroUsize,
    ) -> Result<Self, IoError> {
        let n = threads.get();
        let (mode, states) = if n == 1 {
            let d = Decoder::new(info)?;
            let states = vec![d.state.clone()];
            (Mode::Inline(Box::new(d)), states)
        } else {
            let decoders = (0..n)
                .map(|_| Decoder::new(info.clone()))
                .collect::<Result<Vec<_>, _>>()?;
            let states = decoders.iter().map(|d| d.state.clone()).collect();
            (
                Mode::Pool(Ordered::new(decoders, decode_job)?, 2 * n as u64),
                states,
            )
        };
        Ok(CramStream {
            frames,
            mode,
            cur: Vec::new(),
            at: 0,
            reading_done: false,
            failed: false,
            states,
        })
    }

    /// The next decoded container, in file order.
    fn next_container(&mut self) -> Result<Option<Vec<u8>>, IoError> {
        match &mut self.mode {
            Mode::Inline(d) => match self.frames.next()? {
                None => Ok(None),
                Some(c) => {
                    let mut out = std::mem::take(&mut self.cur);
                    out.clear();
                    d.decode(&c, &mut out)?;
                    Ok(Some(out))
                }
            },
            Mode::Pool(pool, cap) => {
                while !self.reading_done && pool.in_flight() < *cap {
                    match self.frames.next() {
                        Ok(Some(c)) => pool.submit(c),
                        Ok(None) => self.reading_done = true,
                        Err(e) => {
                            pool.submit_result(Err(e));
                            self.reading_done = true;
                        }
                    }
                }
                pool.next().transpose()
            }
        }
    }
}

/// A CRAM decoding failure carried through `io::Error`.
#[derive(Debug)]
struct CramFailure(Box<IoError>);

impl std::fmt::Display for CramFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for CramFailure {}

impl Read for CramStream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let avail = self.fill_buf()?;
        let n = avail.len().min(buf.len());
        buf[..n].copy_from_slice(&avail[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for CramStream {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.at >= self.cur.len() {
            if self.failed {
                return Err(io::Error::other("the CRAM input failed earlier"));
            }
            match self.next_container() {
                Ok(Some(b)) => {
                    self.cur = b;
                    self.at = 0;
                }
                Ok(None) => return Ok(&[]),
                Err(e) => {
                    self.failed = true;
                    return Err(io::Error::other(CramFailure(Box::new(e))));
                }
            }
        }
        Ok(&self.cur[self.at..])
    }

    fn consume(&mut self, n: usize) {
        self.at = (self.at + n).min(self.cur.len());
    }
}

/// The gloci error inside a read error from a CRAM stream, or `e` itself.
pub(crate) fn unwrap_failure(e: IoError) -> IoError {
    match e {
        IoError::Read { path, source }
            if source
                .get_ref()
                .is_some_and(|inner| inner.is::<CramFailure>()) =>
        {
            match source.into_inner().map(|i| i.downcast::<CramFailure>()) {
                Some(Ok(f)) => *f.0,
                _ => IoError::Internal(format!("{path}: a CRAM failure was lost")),
            }
        }
        e => e,
    }
}

// ---- opening ----

/// Check that a seekable CRAM ends with the end-of-file container, before any record is read
/// (as gloci checks a BGZF file's end-of-file block).
fn check_eof_container(name: &str) -> Result<(), IoError> {
    use std::io::{Seek, SeekFrom};
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let mut f = std::fs::File::open(name).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let len = f.metadata().map_err(read_err)?.len();
    let mut tail = Vec::new();
    if len >= CRAM_EOF.len() as u64 {
        f.seek(SeekFrom::Start(len - CRAM_EOF.len() as u64))
            .map_err(read_err)?;
        take_into(&mut f, &mut tail, CRAM_EOF.len()).map_err(read_err)?;
    }
    if tail != CRAM_EOF {
        return Err(truncated(name, "the file does not end with one"));
    }
    Ok(())
}

/// A CRAM input's header, its decoded records as a BAM body, and the counters of its decoders.
struct Opened {
    frames: BamFrames,
    states: Vec<Arc<RefState>>,
}

/// Open the CRAM `input` for the BAM reader's framing: the header (renamed by `map`), and the
/// records decoded on `threads` workers with the reference `reference`.
fn open(
    input: Input,
    map: Option<&ChromMap>,
    reference: Option<&Path>,
    threads: NonZeroUsize,
) -> Result<Opened, IoError> {
    let Input {
        name, mut reader, ..
    } = input;
    let head = read_head(&mut reader, &name)?;
    if is_regular_file(&name) {
        check_eof_container(&name)?;
    }
    let h = parse_header(&name, &head)?;
    let info = Arc::new(CramInfo::new(&name, h.clone(), reference)?);
    let frames = CramFrames::new(&name, reader, &head);
    let stream = CramStream::new(frames, info, threads)?;
    let states = stream.states.clone();
    let mut aln = h;
    if let Some(m) = map {
        m.check_file_names(&name, aln.refs.iter().map(|r| &r.0[..]))?;
        aln = aln.renamed(m);
    }
    Ok(Opened {
        frames: BamFrames::from_parts(name, Compression::None, aln, Box::new(stream)),
        states,
    })
}

fn reference_of(o: &ReaderOptions) -> Option<&Path> {
    o.accept
        .alignments
        .as_ref()
        .and_then(|a| a.reference.as_deref())
}

/// Reads CRAM records for interval commands: the BAM reader over the transcoded records, so
/// spans, skips, filters and warnings are the BAM's (D101, D102).
pub struct CramReader {
    inner: BamReader,
    states: Vec<Arc<RefState>>,
}

impl CramReader {
    pub(crate) fn open_cram(
        input: Input,
        o: &ReaderOptions,
        opts: AlignOptions,
    ) -> Result<CramReader, IoError> {
        let Opened { frames, states } =
            open(input, o.chrom_map.as_deref(), reference_of(o), o.threads)?;
        Ok(CramReader {
            inner: BamReader::from_frames(frames, o.invalid, opts),
            states,
        })
    }

    /// Sequences each decoder read from the FASTA, and the most it held at once.
    pub fn reference_counts(&self) -> Vec<(usize, usize)> {
        counts(&self.states)
    }
}

fn counts(states: &[Arc<RefState>]) -> Vec<(usize, usize)> {
    states
        .iter()
        .map(|s| {
            (
                s.fetches.load(Ordering::Relaxed),
                s.max_live.load(Ordering::Relaxed),
            )
        })
        .collect()
}

impl RecordReader for CramReader {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        Format::Cram
    }
    fn compression(&self) -> Compression {
        Compression::None
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        self.inner
            .read_batch(chroms, batch, max)
            .map_err(unwrap_failure)
    }
}

/// Every record of a CRAM (D120): the BAM source over the transcoded records.
pub(crate) struct CramSource {
    inner: BamSource,
}

impl CramSource {
    pub(crate) fn open(input: Input, o: &ReaderOptions) -> Result<CramSource, IoError> {
        let Opened { frames, .. } =
            open(input, o.chrom_map.as_deref(), reference_of(o), o.threads)?;
        let opts = o.accept.alignments.clone().unwrap_or_default();
        Ok(CramSource {
            inner: BamSource::new(frames, o.invalid, opts),
        })
    }
}

impl AlnSource for CramSource {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn format(&self) -> Format {
        Format::Cram
    }
    fn header(&self) -> &Arc<AlnHeader> {
        self.inner.header()
    }
    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        self.inner.next_record(buf).map_err(unwrap_failure)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
}

fn parse_header(name: &str, head: &CramHead) -> Result<AlnHeader, IoError> {
    AlnHeader::parse_sam(head.text.clone()).map_err(|m| {
        detect_err(
            name,
            format!("{name}: the SAM header of the CRAM file: {m}"),
        )
    })
}

/// The alignment header of the CRAM file `name`, read from `r` (for `-g x.cram`); no reference
/// is needed.
pub fn read_cram_header(r: &mut dyn Read, name: &str) -> Result<AlnHeader, IoError> {
    parse_header(name, &read_head(r, name)?)
}

/// Warn that `--reference` is not used when no input is CRAM (D109).
pub fn warn_unused_reference(reference: Option<&Path>, any_cram: bool) {
    if let Some(r) = reference
        && !any_cram
    {
        log::warn!(
            "--reference {} is not used: no input is CRAM (only CRAM input needs a reference)",
            r.display()
        );
    }
}

// ---- for the index ----

/// The records of a slice that holds several references, for the `.crai`: each record's
/// reference id (-1 for none), start and end, as htslib's index builder sees them. No reference
/// bases are read: the adapter is told the container needs none (so it gives empty
/// placeholders), and the reads' bases are never put together.
pub(crate) struct ExtentDecoder {
    info: Arc<CramInfo>,
    repo: fasta::Repository,
}

impl ExtentDecoder {
    pub fn new(info: Arc<CramInfo>) -> ExtentDecoder {
        let state = Arc::new(RefState::default());
        state.required.store(false, Ordering::Relaxed);
        let adapter = RefAdapter {
            info: info.clone(),
            fasta: None,
            state,
        };
        ExtentDecoder {
            info,
            repo: fasta::Repository::new(adapter),
        }
    }

    pub fn info(&self) -> &CramInfo {
        &self.info
    }

    /// `(ref id, start, end)` of each record of slice `index` (0-based) of `c`, a slice of
    /// several references. A read covers at least its first base, and an unmapped read only
    /// that; an unplaced read is at 0.
    pub fn extents(
        &mut self,
        c: &RawContainer,
        index: usize,
    ) -> Result<Vec<(i32, i64, i64)>, IoError> {
        let name = &self.info.name;
        let bad = |e: &dyn std::fmt::Display| {
            detect_err(
                name,
                format!(
                    "{name}: {} cannot be decoded: {e}; the file is corrupt",
                    c.place()
                ),
            )
        };
        if let Some(what) = compression(c).unsupported {
            return Err(unsupported_codec(name, c, &what));
        }
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&c.bytes[..])
            .read_container(&mut container)
            .map_err(|e| bad(&e))?;
        let ch = container.compression_header().map_err(|e| bad(&e))?;
        let Some(slice) = container.slices().nth(index) else {
            return Err(bad(&"a slice is missing"));
        };
        let slice = slice.map_err(|e| bad(&e))?;
        let (core, ext) = slice.decode_blocks().map_err(|e| bad(&e))?;
        self.repo.clear();
        let records = slice
            .records(self.repo.clone(), &self.info.sam, &ch, &core, &ext)
            .map_err(|e| bad(&e))?;
        let mut out = Vec::with_capacity(records.len());
        for r in &records {
            let flags = r.flags().map_err(|e| bad(&e))?;
            let rid = match r.reference_sequence_id(&self.info.sam) {
                Some(Ok(id)) => id as i32,
                Some(Err(e)) => return Err(bad(&e)),
                None => -1,
            };
            let start = match r.alignment_start() {
                Some(Ok(p)) => usize::from(p) as i64,
                Some(Err(e)) => return Err(bad(&e)),
                None => 0,
            };
            // htslib: a read covers at least its first base, and an unmapped read only that.
            let span = if flags.is_unmapped() {
                1
            } else {
                match r.alignment_span() {
                    Some(Ok(s)) => (s as i64).max(1),
                    _ => 1,
                }
            };
            out.push((rid, start, start + span - 1));
        }
        Ok(out)
    }
}

/// The header of the CRAM `name` as decoders share it, and its frames.
pub(crate) fn open_for_index(
    name: &str,
    src: Box<dyn BufRead + Send>,
    reference: Option<&Path>,
) -> Result<(Arc<CramInfo>, CramFrames), IoError> {
    let mut src = src;
    let head = read_head(&mut src, name)?;
    if is_regular_file(name) {
        check_eof_container(name)?;
    }
    let h = parse_header(name, &head)?;
    let info = Arc::new(CramInfo::new(name, h, reference)?);
    Ok((info, CramFrames::new(name, src, &head)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diag::InvalidPolicy;
    use crate::input::{InputOptions, open_input};
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../xtask/diff/inputs")
            .join(name)
    }

    fn options(reference: Option<PathBuf>, threads: usize) -> ReaderOptions {
        ReaderOptions {
            threads: NonZeroUsize::new(threads).unwrap(),
            format: None,
            vcf_span: gloci_core::VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
            chrom_map: None,
            accept: crate::align::Accept {
                alignments: Some(AlignOptions {
                    full_records: true,
                    reference,
                    ..AlignOptions::default()
                }),
                ..crate::align::Accept::default()
            },
        }
    }

    fn drain(path: &Path, threads: usize) -> (Vec<Vec<u8>>, Vec<(usize, usize)>) {
        let o = options(Some(fixture("fa/ref.fa")), threads);
        let input = open_input(
            &path.to_string_lossy(),
            &InputOptions {
                threads: o.threads,
                format: None,
            },
        )
        .unwrap();
        let opts = o.accept.alignments.clone().unwrap();
        let mut r = CramReader::open_cram(input, &o, opts).unwrap();
        let Header::Alignment(h) = r.header().clone() else {
            panic!("no alignment header");
        };
        let mut chroms = ChromTable::first_appearance();
        let mut recs = Vec::new();
        let mut batch = RecordBatch::default();
        loop {
            batch.clear();
            let more = r.read_batch(&mut chroms, &mut batch, 3).unwrap();
            for i in 0..batch.len() {
                assert!(crate::align::Aln::parse(batch.raw(i), &h).is_ok());
                recs.push(batch.raw(i).to_vec());
            }
            if !more {
                break;
            }
        }
        (recs, r.reference_counts())
    }

    /// A counting adapter: each decoder holds at most one sequence of a single-reference run, and
    /// reads it once per run (D109).
    #[test]
    fn cram_reference_cache_bounded() {
        // flags.multi.cram: slices of several references (as multi_seq_per_slice=1 writes every
        // slice), five on chr1, then chr1 (an unmapped read) + chr2, then chr2 + unplaced reads:
        // chr1 and chr2 are each read once. tags.cram: chrM, chr1, then chr2 + chr10 (an
        // unmapped read).
        for (file, fetches) in [("aln/flags.multi.cram", 2), ("aln/tags.cram", 3)] {
            let (recs, counts) = drain(&fixture(file), 1);
            assert!(!recs.is_empty());
            assert_eq!(counts.len(), 1);
            let (got, live) = counts[0];
            assert_eq!(got, fetches, "{file}: sequences read");
            assert!(live <= 1, "{file}: {live} sequences held at once");
            // On every worker of a pool: one sequence each, at most.
            for t in [2, 8] {
                let (again, counts) = drain(&fixture(file), t);
                assert_eq!(again, recs, "{file} at -t {t}");
                assert_eq!(counts.len(), t);
                for (_, live) in counts {
                    assert!(live <= 1, "{file} at -t {t}: {live} sequences held at once");
                }
            }
        }
    }

    /// Encodings noodles 0.102 would panic on are named before it sees them.
    #[test]
    fn unsupported_encodings_are_found() {
        let check = |b: &[u8]| encoding_problem(&mut &b[..], 0).unwrap();
        // EXTERNAL(block 5), BETA(offset 0, 8 bits), a one-symbol HUFFMAN: fine.
        assert_eq!(check(&[1, 1, 5]), None);
        assert_eq!(check(&[6, 2, 0, 8]), None);
        assert_eq!(check(&[3, 4, 1, 0x41, 1, 0]), None);
        // GOLOMB, SUBEXP, GOLOMB_RICE, an empty HUFFMAN, and GOLOMB inside BYTE_ARRAY_LEN.
        assert_eq!(check(&[2, 2, 1, 3]).as_deref(), Some("Golomb"));
        assert_eq!(check(&[7, 2, 0, 7]).as_deref(), Some("Subexp"));
        assert_eq!(check(&[8, 2, 0, 2]).as_deref(), Some("Golomb-Rice"));
        assert!(check(&[3, 2, 0, 0]).is_some());
        assert_eq!(
            check(&[4, 7, 2, 2, 1, 3, 1, 1, 5]).as_deref(),
            Some("Golomb")
        );
        // A data series map of two entries, the second one GOLOMB.
        let map = [2, b'B', b'F', 1, 1, 5, b'R', b'L', 2, 2, 1, 3];
        assert_eq!(map_problem(&map, false).as_deref(), Some("Golomb"));
        assert_eq!(encoding_problem(&mut &[2, 9, 1][..], 0), None);
    }

    #[test]
    fn numbers_decode_as_htslib_writes_them() {
        let mut s: &[u8] = &[0xe0, 0x45, 0x4f, 0x46];
        assert_eq!(itf8(&mut s), Some(4_542_278));
        let mut s: &[u8] = &[0xff, 0xff, 0xff, 0xff, 0x0f];
        assert_eq!(itf8(&mut s), Some(-1));
        let mut s: &[u8] = &[0x80, 0xff];
        assert_eq!(itf8(&mut s), Some(255));
        let mut s: &[u8] = &[0x7f];
        assert_eq!(ltf8(&mut s), Some(127));
        let mut s: &[u8] = &[0x81, 0x00];
        assert_eq!(ltf8(&mut s), Some(256));
        let mut s: &[u8] = &[0xff, 0, 0, 0, 0, 0, 0, 1, 0];
        assert_eq!(ltf8(&mut s), Some(256));
        let mut s: &[u8] = &[0x80];
        assert_eq!(itf8(&mut s), None);
    }
}
