//! CRAM framing, read by gloci before noodles-cram sees a byte (D109): ITF8/LTF8 numbers, the
//! file definition and header container, data containers (CRC-checked, read through `Take`),
//! blocks, slice headers checked against the CRAM header, and the compression header's `RR` and
//! encodings (those noodles 0.102 would panic on are named).

use crate::error::IoError;
use gloci_core::thousands;
use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{self, BufRead, Read};

/// The 38 bytes of a CRAM 3 end-of-file container.
pub(crate) const CRAM_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// Most bytes a header block (file header, compression header, slice header) may expand to.
pub(crate) const MAX_HEADER_BLOCK: usize = 256 << 20;
/// Most bytes any block may expand to (samtools' slices hold a few megabytes of each series).
pub(crate) const MAX_BLOCK: usize = 1 << 30;
/// Most slices a container may claim.
pub(crate) const MAX_LANDMARKS: usize = 1 << 20;

// ---- numbers ----

/// An ITF8 integer (1 to 5 bytes).
pub(crate) fn itf8(src: &mut &[u8]) -> Option<i32> {
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
pub(crate) fn ltf8(src: &mut &[u8]) -> Option<i64> {
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

pub(crate) fn u32_le(src: &mut &[u8]) -> Option<u32> {
    let (b, rest) = src.split_first_chunk::<4>()?;
    *src = rest;
    Some(u32::from_le_bytes(*b))
}

pub(crate) fn crc32(b: &[u8]) -> u32 {
    let mut c = flate2::Crc::new();
    c.update(b);
    c.sum()
}

/// Append up to `n` bytes of `r` to `buf`; returns how many came (fewer only at the end).
pub(crate) fn take_into(r: &mut dyn Read, buf: &mut Vec<u8>, n: usize) -> io::Result<usize> {
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
    pub(crate) fn is_eof(&self) -> bool {
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

    pub(crate) fn place(&self) -> String {
        format!(
            "container {} (at byte {})",
            thousands(self.number),
            thousands(self.offset)
        )
    }
}

/// Why a container header could not be read.
pub(crate) enum HeaderProblem {
    /// The stream ended before the first byte.
    End,
    /// The stream ended inside the header.
    Short,
    Io(io::Error),
    Bad(String),
}

/// Read a container header from `r`, appending its bytes to `raw`.
pub(crate) fn read_container_header(
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
pub(crate) struct Block<'a> {
    pub(crate) method: u8,
    pub(crate) content_type: u8,
    pub(crate) content_id: i32,
    pub(crate) raw_size: usize,
    pub(crate) data: &'a [u8],
}

/// Read the block at the front of `src` (CRAM 3: with its CRC32).
pub(crate) fn read_block<'a>(src: &mut &'a [u8]) -> Result<Block<'a>, String> {
    let all = *src;
    let short = || "a block is cut short".to_string();
    let (&method, rest) = src.split_first().ok_or_else(short)?;
    let (&content_type, mut rest) = rest.split_first().ok_or_else(short)?;
    let content_id = itf8(&mut rest).ok_or_else(short)?;
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
        content_id,
        raw_size,
        data,
    })
}

/// The contents of a raw or gzip block (the only methods header blocks use).
pub(crate) fn block_contents<'a>(b: &Block<'a>, what: &str) -> Result<Cow<'a, [u8]>, String> {
    match b.method {
        0 => Ok(Cow::Borrowed(b.data)),
        1 => {
            if b.raw_size > MAX_HEADER_BLOCK {
                return Err(format!(
                    "the {what} block claims {} bytes; gloci reads at most {} in one",
                    thousands(b.raw_size as u64),
                    thousands(MAX_HEADER_BLOCK as u64)
                ));
            }
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

pub(crate) fn detect_err(name: &str, msg: String) -> IoError {
    IoError::Detect {
        path: name.to_string(),
        msg,
    }
}

pub(crate) fn corrupt(name: &str, what: &str) -> IoError {
    detect_err(
        name,
        format!("{name}: {what}; the file is not valid CRAM, or it is corrupt"),
    )
}

/// A container whose encodings gloci (noodles) cannot decode.
pub(crate) fn unsupported_codec(name: &str, c: &RawContainer, what: &str) -> IoError {
    detect_err(
        name,
        format!(
            "{name}: {} uses the CRAM encoding {what}, which gloci does not decode (samtools \
does not write it); rewrite the file with samtools: samtools view -C -T ref.fa -o new.cram {name}",
            c.place()
        ),
    )
}

pub(crate) fn truncated(name: &str, after: &str) -> IoError {
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

    /// Frames of `src`, which stands at the data container at byte `offset` of the file (a
    /// seek made through the `.crai`); containers are numbered from `first_number`.
    pub fn at(
        name: &str,
        src: Box<dyn BufRead + Send>,
        offset: u64,
        first_number: u64,
    ) -> CramFrames {
        CramFrames {
            name: name.to_string(),
            src,
            pos: offset,
            containers: first_number.saturating_sub(1),
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
    /// The records the slice claims (checked against its bytes before decoding).
    pub records: usize,
    /// The external block holding the slice's reference bases, or -1.
    pub embedded: i32,
    /// Where the slice starts in the container's data, and its bytes.
    pub landmark: usize,
    pub size: usize,
}

impl SliceHead {
    /// The last reference base the slice covers (1-based), for a slice on one reference.
    pub(crate) fn end(&self) -> i64 {
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
        let records = usize::try_from(records)
            .map_err(|_| bad(format!("slice {} claims {records} records", i + 1)))?;
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
            records,
            embedded,
            landmark: at,
            size: end - at,
        });
    }
    Ok(out)
}

/// What gloci reads of a container's compression header before noodles does.
#[derive(Debug, Default)]
pub(crate) struct CompressionInfo {
    /// The preservation map's `RR` (true unless the map says otherwise): the container needs
    /// reference bases. A header gloci cannot read is taken to need them.
    pub(crate) required: bool,
    /// An encoding noodles cannot decode (it would panic), named.
    pub(crate) unsupported: Option<String>,
    /// The external block of each data series stored in one (`EXTERNAL`).
    pub(crate) external: HashMap<[u8; 2], i32>,
    /// How many encodings (data series, tags, nested ones) read each external block.
    pub(crate) uses: HashMap<i32, usize>,
    /// The value of each data series coded as a one-symbol Huffman table: every record has it.
    pub(crate) constant: HashMap<[u8; 2], i32>,
}

impl CompressionInfo {
    /// The external block of data series `key` when no other encoding reads it, so its
    /// contents are that series' values alone.
    pub(crate) fn own_block(&self, key: &[u8; 2]) -> Option<i32> {
        let id = *self.external.get(key)?;
        (self.uses.get(&id) == Some(&1)).then_some(id)
    }
}

/// The preservation map's `RR`, from its contents `s`.
pub(crate) fn preservation_rr(mut s: &[u8]) -> bool {
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

/// One encoding read by [`scan_encoding`].
pub(crate) struct Scanned {
    /// The external block an `EXTERNAL` encoding reads.
    pub(crate) external: Option<i32>,
    /// The one symbol of a one-symbol Huffman table.
    pub(crate) constant: Option<i32>,
    /// A codec noodles 0.102 would panic on: Golomb, Subexp and Golomb-Rice (`todo!`), and a
    /// Huffman table that is empty or whose symbol and length counts differ.
    pub(crate) problem: Option<String>,
}

/// Read the encoding at the front of `s` (and those nested in it), counting in `uses` the
/// external blocks each reads. `None` when it cannot be read (noodles reports that itself).
pub(crate) fn scan_encoding(
    s: &mut &[u8],
    depth: u8,
    uses: &mut HashMap<i32, usize>,
) -> Option<Scanned> {
    let id = itf8(s)?;
    let len = usize::try_from(itf8(s)?).ok()?;
    let mut params = s.get(..len)?;
    *s = &s[len..];
    let mut out = Scanned {
        external: None,
        constant: None,
        problem: None,
    };
    match id {
        1 => {
            let block = itf8(&mut params)?;
            *uses.entry(block).or_default() += 1;
            out.external = Some(block);
        }
        2 => out.problem = Some("Golomb".to_string()),
        7 => out.problem = Some("Subexp".to_string()),
        8 => out.problem = Some("Golomb-Rice".to_string()),
        3 => {
            let n = itf8(&mut params)?;
            for k in 0..n.max(0) {
                let symbol = itf8(&mut params)?;
                if n == 1 && k == 0 {
                    out.constant = Some(symbol);
                }
            }
            let m = itf8(&mut params)?;
            out.problem =
                (n < 1 || n != m).then(|| format!("Huffman with {n} symbols and {m} lengths"));
        }
        4 if depth < 8 => {
            let first = scan_encoding(&mut params, depth + 1, uses)?;
            let second = scan_encoding(&mut params, depth + 1, uses)?;
            out.problem = first.problem.or(second.problem);
        }
        5 => {
            // BYTE_ARRAY_STOP: the stop byte, then the external block.
            params = params.get(1..)?;
            let block = itf8(&mut params)?;
            *uses.entry(block).or_default() += 1;
        }
        _ => {}
    }
    Some(out)
}

/// Scan a data series (`wide` false: two-byte keys) or tag encoding map (ITF8 keys) into
/// `info`; returns the first unsupported encoding.
pub(crate) fn scan_map(mut s: &[u8], wide: bool, info: &mut CompressionInfo) -> Option<String> {
    let n = itf8(&mut s)?;
    let mut problem = None;
    for _ in 0..n.max(0) {
        let key = if wide {
            itf8(&mut s)?;
            None
        } else {
            let (k, rest) = s.split_first_chunk::<2>()?;
            s = rest;
            Some(*k)
        };
        let e = scan_encoding(&mut s, 0, &mut info.uses)?;
        if let Some(k) = key {
            if let Some(b) = e.external {
                info.external.insert(k, b);
            }
            if let Some(v) = e.constant {
                info.constant.insert(k, v);
            }
        }
        problem = problem.or(e.problem);
    }
    problem
}

/// Read the compression header of container `c` (see [`CompressionInfo`]).
pub(crate) fn compression(c: &RawContainer) -> CompressionInfo {
    let mut info = CompressionInfo {
        required: true,
        ..CompressionInfo::default()
    };
    let data = c.data();
    let end = c.header.landmarks.first().copied().unwrap_or(data.len());
    let mut src = &data[..end];
    let Ok(block) = read_block(&mut src) else {
        return info;
    };
    let Ok(contents) = block_contents(&block, "compression header") else {
        return info;
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
        return info;
    };
    info.required = preservation_rr(preservation);
    let series = maps.get(1).and_then(|m| scan_map(m, false, &mut info));
    let tags = maps.get(2).and_then(|m| scan_map(m, true, &mut info));
    info.unsupported = series.or(tags);
    info
}

/// A uint7 (7 bits a byte, most significant first, up to 5 bytes) at the front of `src`.
fn uint7(src: &[u8]) -> Option<u64> {
    let mut n = 0u64;
    for (i, &b) in src.iter().take(5).enumerate() {
        n = (n << 7) | u64::from(b & 0x7f);
        if b & 0x80 == 0 {
            return Some(n);
        }
        if i == 4 {
            return None;
        }
    }
    None
}

/// Check block `b` before noodles decodes it: noodles allocates the decoded size a block
/// states, in the block header (gzip, bzip2, lzma) or in the codec's clear-text header (rANS
/// 4x8, rANS Nx16, the arithmetic coder, fqzcomp, the name tokeniser), so a corrupt one asks
/// for gigabytes. The block header's size (checksummed) must be at most [`MAX_BLOCK`], and the
/// codec's must equal it. A name-tokeniser block must hold no more names than its bytes, and,
/// as the block of the read names, no more than the slice's `records`.
pub(crate) fn check_block(b: &Block<'_>, names_of: Option<usize>) -> Result<(), String> {
    if b.raw_size == 0 || b.method == 0 {
        // noodles reads an empty or raw block as it is.
        return Ok(());
    }
    if b.raw_size > MAX_BLOCK {
        return Err(format!(
            "a block says it decodes to {} bytes; gloci decodes at most {} in one block",
            thousands(b.raw_size as u64),
            thousands(MAX_BLOCK as u64)
        ));
    }
    let u32_at = |at: usize| {
        b.data
            .get(at..at + 4)
            .map(|v| u64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])))
    };
    let (codec, stated) = match b.method {
        4 => ("rANS 4x8", u32_at(5)),
        5 | 6 => {
            let codec = if b.method == 5 {
                "rANS Nx16"
            } else {
                "arithmetic coder"
            };
            match b.data.first() {
                // NO_SIZE: noodles takes the block header's size.
                Some(&flags) if flags & 0x10 != 0 => return Ok(()),
                Some(_) => (codec, uint7(&b.data[1..])),
                None => (codec, None),
            }
        }
        7 => ("fqzcomp", uint7(b.data)),
        8 => ("name-tokeniser", u32_at(0)),
        _ => return Ok(()),
    };
    let Some(stated) = stated else {
        return Err(format!("a {codec} block is cut short"));
    };
    if stated != b.raw_size as u64 {
        return Err(format!(
            "a {codec} block says it decodes to {} bytes, and its header {}",
            thousands(stated),
            thousands(b.raw_size as u64)
        ));
    }
    if b.method == 8 {
        let names = u32_at(4).ok_or_else(|| "a name-tokeniser block is cut short".to_string())?;
        if names > stated || names_of.is_some_and(|n| names > n as u64) {
            return Err(format!(
                "a name-tokeniser block claims {} names in {} bytes, for a slice of {} records",
                thousands(names),
                thousands(stated),
                names_of.map_or_else(|| "?".to_string(), |n| thousands(n as u64))
            ));
        }
    }
    Ok(())
}

/// Check every block of slice `head` of `c` ([`check_block`]); `rn` is the external block of
/// the read names.
pub(crate) fn check_slice_blocks(
    c: &RawContainer,
    head: &SliceHead,
    rn: Option<i32>,
) -> Result<(), String> {
    let data = c.data();
    let Some(mut src) = data.get(head.landmark..head.landmark + head.size) else {
        return Ok(());
    };
    while !src.is_empty() {
        let Ok(b) = read_block(&mut src) else {
            // noodles reports a broken block itself.
            return Ok(());
        };
        let names_of = (rn == Some(b.content_id)).then_some(head.records);
        check_block(&b, names_of)?;
    }
    Ok(())
}

/// Check the compression header block of `c` ([`check_block`]).
pub(crate) fn check_compression_block(c: &RawContainer) -> Result<(), String> {
    let data = c.data();
    let end = c.header.landmarks.first().copied().unwrap_or(data.len());
    let mut src = &data[..end];
    match read_block(&mut src) {
        Ok(b) => check_block(&b, None),
        Err(_) => Ok(()),
    }
}

/// The ITF8 values of `block`, all of it; `None` if one is cut short.
pub(crate) fn itf8_values(block: &[u8]) -> Option<Vec<i32>> {
    let mut s = block;
    let mut out = Vec::new();
    while !s.is_empty() {
        out.push(itf8(&mut s)?);
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encodings noodles 0.102 would panic on are named before it sees them.
    #[test]
    fn unsupported_encodings_are_found() {
        let check = |b: &[u8]| {
            scan_encoding(&mut &b[..], 0, &mut HashMap::new())
                .unwrap()
                .problem
        };
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
        let mut info = CompressionInfo::default();
        assert_eq!(scan_map(&map, false, &mut info).as_deref(), Some("Golomb"));
        assert_eq!(info.own_block(b"BF"), Some(5));
        assert!(scan_encoding(&mut &[2, 9, 1][..], 0, &mut HashMap::new()).is_none());
        // A block two series read is no series' own.
        let map = [2, b'B', b'F', 1, 1, 5, b'R', b'L', 1, 1, 5];
        let mut info = CompressionInfo::default();
        assert_eq!(scan_map(&map, false, &mut info), None);
        assert_eq!(info.own_block(b"BF"), None);
    }

    /// Each codec's clear-text decoded size must be the block header's, and no block may say
    /// it decodes to more than 1 GiB: noodles allocates by them.
    #[test]
    fn block_sizes_are_checked() {
        let block = |method: u8, raw_size: usize, data: &'static [u8]| Block {
            method,
            content_type: 4,
            content_id: 1,
            raw_size,
            data,
        };
        // rANS 4x8: order, compressed size, then the decoded size (bytes 5-8).
        let ok = block(4, 100, &[0, 9, 0, 0, 0, 100, 0, 0, 0, 1, 2]);
        assert!(check_block(&ok, None).is_ok());
        let bad = block(4, 100, &[0, 9, 0, 0, 0, 0x64, 0, 0, 0xfe, 1, 2]);
        let e = check_block(&bad, None).unwrap_err();
        assert!(
            e.contains(
                "a rANS 4x8 block says it decodes to 4,261,412,964 bytes, and its header 100"
            ),
            "{e}"
        );
        // rANS Nx16 and the arithmetic coder: flags, then a uint7 size unless NO_SIZE (0x10).
        assert!(check_block(&block(5, 300, &[0x01, 0x82, 0x2c, 7]), None).is_ok());
        assert!(check_block(&block(6, 300, &[0x00, 0x82, 0x2d, 7]), None).is_err());
        assert!(check_block(&block(5, 300, &[0x10, 0xff, 0xff, 0x7f]), None).is_ok());
        // fqzcomp: a uint7 size first.
        assert!(check_block(&block(7, 5, &[5, 1, 2]), None).is_ok());
        assert!(check_block(&block(7, 5, &[0x86, 0, 1]), None).is_err());
        // gzip (and bzip2, lzma) take the block header's size: at most 1 GiB.
        let e = check_block(&block(1, (1 << 30) + 1, &[0x1f, 0x8b]), None).unwrap_err();
        assert!(
            e.contains("gloci decodes at most 1,073,741,824 in one block"),
            "{e}"
        );
        // The name tokeniser: decoded size, then names (no more than the slice's records).
        let tok: &'static [u8] = &[10, 0, 0, 0, 3, 0, 0, 0, 0];
        assert!(check_block(&block(8, 10, tok), Some(3)).is_ok());
        assert!(check_block(&block(8, 10, tok), Some(2)).is_err());
    }

    /// cat-big.v31.cram reaches the codecs whose tables noodles keeps on the stack (the pool's
    /// workers get 16 MiB stacks for them): order-1 rANS Nx16 and the name tokeniser.
    #[test]
    fn cat_big_31_reaches_order_1_rans_and_the_tokeniser() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../xtask/diff/inputs/aln/cat-big.v31.cram");
        let mut f = std::io::BufReader::new(std::fs::File::open(&path).unwrap());
        let head = read_head(&mut f, "x").unwrap();
        let mut frames = CramFrames::new("x", Box::new(f), &head);
        let (mut order_1, mut tokeniser) = (false, false);
        while let Some(c) = frames.next().unwrap() {
            let mut data = c.data();
            while !data.is_empty() {
                let b = read_block(&mut data).unwrap();
                // rANS Nx16 (method 5): the first byte holds the flags, ORDER = 1.
                order_1 |= b.method == 5 && b.data.first().is_some_and(|f| f & 1 != 0);
                tokeniser |= b.method == 8;
            }
        }
        assert!(order_1, "no order-1 rANS block");
        assert!(tokeniser, "no tokeniser block");
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
