//! BCF: header framing, dictionaries, and a record reader that passes framed
//! records through untouched (D22).
//!
//! Only CHROM, POS, the alleles and, for records with a symbolic ALT,
//! INFO/END and INFO/SVLEN are decoded. END and SVLEN are rendered back to
//! their VCF text form and go through the same parsing as the text reader, so
//! a BCF record gets exactly the footprint of the equivalent VCF line.

use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::Input;
use crate::record::{Header, RecordBatch, RecordMeta, RecordReader, coord, intern_chrom};
use crate::vcf::{VcfHeader, header_fields};
use gloci_core::{
    ChromId, ChromTable, FootprintError, Interval, Strand, Variant, VcfSpan, footprint,
};
use std::collections::HashMap;
use std::io::{self, BufRead, Read, Write};

/// The BCF dictionaries of a header; the index is the dictionary id.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BcfDicts {
    /// `##contig` IDs in header order; `IDX=` overrides the position.
    pub contigs: Vec<Vec<u8>>,
    /// `PASS`, then the FILTER/INFO/FORMAT IDs in header order without
    /// duplicates; `IDX=` overrides the position. Ids that no line claims
    /// (gaps left by `IDX=`) are empty.
    pub strings: Vec<Vec<u8>>,
}

impl BcfDicts {
    /// Build the dictionaries from header text (BCF2.2 §6.2.1, §6.2.2).
    pub(crate) fn from_header(raw: &[u8]) -> Result<BcfDicts, String> {
        // A valid IDX can't be much larger than the number of header lines;
        // the cap keeps a corrupt IDX from allocating gigabytes.
        let cap = raw.len() + 65_536;
        let mut strings = Dict::default();
        strings.insert(b"PASS", None, b"")?;
        let mut contigs = Dict::default();
        for line in raw.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            let (dict, body) = if let Some(b) = line.strip_prefix(b"##contig=<") {
                (&mut contigs, b)
            } else if let Some(b) = [&b"##INFO=<"[..], b"##FILTER=<", b"##FORMAT=<"]
                .iter()
                .find_map(|p| line.strip_prefix(*p))
            {
                (&mut strings, b)
            } else {
                continue;
            };
            let fields = header_fields(body);
            let field = |key: &[u8]| fields.iter().find(|f| f.0 == key).map(|f| f.1);
            let Some(id) = field(b"ID").filter(|id| !id.is_empty()) else {
                continue;
            };
            let idx = match field(b"IDX") {
                None => None,
                Some(v) => Some(
                    std::str::from_utf8(v)
                        .ok()
                        .and_then(|s| s.parse::<usize>().ok())
                        .filter(|&i| i <= cap)
                        .ok_or_else(|| format!("header line {} has an invalid IDX", lossy(line)))?,
                ),
            };
            dict.insert(id, idx, line)?;
        }
        Ok(BcfDicts {
            contigs: contigs.names,
            strings: strings.names,
        })
    }
}

#[derive(Default)]
struct Dict {
    names: Vec<Vec<u8>>,
    ids: HashMap<Vec<u8>, usize>,
}

impl Dict {
    fn insert(&mut self, id: &[u8], idx: Option<usize>, line: &[u8]) -> Result<(), String> {
        match (self.ids.get(id).copied(), idx) {
            (Some(have), Some(want)) if have != want => Err(format!(
                "header line {} has IDX={want}, but {} already has IDX={have}",
                lossy(line),
                lossy(id)
            )),
            (Some(_), _) => Ok(()),
            (None, idx) => {
                let i = idx.unwrap_or(self.names.len());
                if i >= self.names.len() {
                    self.names.resize(i + 1, Vec::new());
                } else if !self.names[i].is_empty() {
                    return Err(format!(
                        "header line {} has IDX={i}, which {} already uses",
                        lossy(line),
                        lossy(&self.names[i])
                    ));
                }
                self.names[i] = id.to_vec();
                self.ids.insert(id.to_vec(), i);
                Ok(())
            }
        }
    }
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// Read until `buf` is full or input ends; returns the bytes read.
fn read_full(r: &mut dyn Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut got = 0;
    while got < buf.len() {
        match r.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(k) => got += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(got)
}

/// Read the BCF magic, version and header text (BCF 2.1 or 2.2), and build
/// the dictionaries. Errors do not name the file; the caller adds it.
pub fn read_bcf_header(r: &mut dyn Read) -> Result<VcfHeader, String> {
    const HEADER_TRUNCATED: &str = "BCF header is truncated";
    let read_err = |e: io::Error| format!("cannot read the BCF header: {e}");
    let mut head = [0u8; 9];
    let got = read_full(r, &mut head).map_err(read_err)?;
    let m = got.min(3);
    if head[..m] != b"BCF"[..m] {
        return Err(
            "not a BCF file (it does not start with the BCF magic bytes); check --format".into(),
        );
    }
    if got < 5 {
        return Err(HEADER_TRUNCATED.into());
    }
    let (major, minor) = (head[3], head[4]);
    if major != 2 || !(1..=2).contains(&minor) {
        return Err(format!(
            "unsupported BCF version {major}.{minor}; gloci reads BCF 2.1 and 2.2"
        ));
    }
    if got < 9 {
        return Err(HEADER_TRUNCATED.into());
    }
    let l_text = u64::from(u32::from_le_bytes([head[5], head[6], head[7], head[8]]));
    let mut text = Vec::new();
    Read::take(&mut *r, l_text)
        .read_to_end(&mut text)
        .map_err(read_err)?;
    if (text.len() as u64) < l_text {
        return Err(HEADER_TRUNCATED.into());
    }
    if let Some(nul) = memchr::memchr(0, &text) {
        text.truncate(nul);
    }
    if !text.is_empty() && !text.ends_with(b"\n") {
        text.push(b'\n');
    }
    let mut h = VcfHeader::parse(text)?;
    h.bcf = Some(BcfDicts::from_header(&h.raw)?);
    Ok(h)
}

/// Write the BCF 2.2 magic, `l_text` and the NUL-terminated header text.
pub fn write_bcf_header(w: &mut dyn Write, header_text: &[u8]) -> io::Result<()> {
    let l_text = u32::try_from(header_text.len() + 1).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "the BCF header text is larger than 4 GiB",
        )
    })?;
    w.write_all(b"BCF\x02\x02")?;
    w.write_all(&l_text.to_le_bytes())?;
    w.write_all(header_text)?;
    w.write_all(&[0])
}

/// BCF typed values (BCF2.2 §6.3.3).
pub mod typed {
    use std::io::Write;

    const PAST_END: &str = "a typed value runs past the end of the record's shared section";
    const BAD_LEN: &str = "a typed value has an invalid length";

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Ty {
        Missing,
        Int8,
        Int16,
        Int32,
        Float,
        Char,
    }

    impl Ty {
        fn from_code(c: u8) -> Result<Ty, String> {
            Ok(match c {
                0 => Ty::Missing,
                1 => Ty::Int8,
                2 => Ty::Int16,
                3 => Ty::Int32,
                5 => Ty::Float,
                7 => Ty::Char,
                c => return Err(format!("unknown BCF value type {c}")),
            })
        }

        pub fn width(self) -> usize {
            match self {
                Ty::Missing => 0,
                Ty::Int8 | Ty::Char => 1,
                Ty::Int16 => 2,
                Ty::Int32 | Ty::Float => 4,
            }
        }
    }

    /// One element of an integer vector.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Int {
        Value(i32),
        Missing,
        /// End-of-vector padding: no further elements.
        End,
    }

    /// One element of a float vector.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum Float {
        Value(f32),
        Missing,
        End,
    }

    /// A typed value: `n` elements of `ty` in `data`.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct Value<'a> {
        pub ty: Ty,
        pub n: usize,
        pub data: &'a [u8],
    }

    fn take_u8(src: &mut &[u8]) -> Result<u8, String> {
        let (&b, rest) = src.split_first().ok_or(PAST_END)?;
        *src = rest;
        Ok(b)
    }

    /// Read one typed value from the front of `src` and advance past it.
    /// The type is in the low 4 bits of the descriptor and the count in the
    /// high 4 bits; a count of 15 means the count follows as a typed integer.
    pub fn read<'a>(src: &mut &'a [u8]) -> Result<Value<'a>, String> {
        let (ty, n) = read_descriptor(src)?;
        let len = n
            .checked_mul(ty.width())
            .filter(|&l| l <= src.len())
            .ok_or(PAST_END)?;
        let (data, rest) = src.split_at(len);
        *src = rest;
        Ok(Value { ty, n, data })
    }

    /// Read only a type descriptor (type and element count) and advance past
    /// it. FORMAT fields use this: the data that follows holds the count for
    /// every sample, not just one.
    pub fn read_descriptor(src: &mut &[u8]) -> Result<(Ty, usize), String> {
        let d = take_u8(src)?;
        let ty = Ty::from_code(d & 0x0f)?;
        let mut n = usize::from(d >> 4);
        if n == 15 {
            n = overflow_len(src)?;
        }
        Ok((ty, n))
    }

    /// The count of an overflowed typed value: a single int8/16/32 (never
    /// itself overflowed, so a corrupt record can't recurse).
    fn overflow_len(src: &mut &[u8]) -> Result<usize, String> {
        let d = take_u8(src)?;
        let ty = Ty::from_code(d & 0x0f)?;
        if d >> 4 != 1 || !matches!(ty, Ty::Int8 | Ty::Int16 | Ty::Int32) {
            return Err(BAD_LEN.into());
        }
        let w = ty.width();
        if src.len() < w {
            return Err(PAST_END.into());
        }
        let v = Value {
            ty,
            n: 1,
            data: &src[..w],
        };
        *src = &src[w..];
        match v.int(0) {
            Some(Int::Value(n)) if n >= 0 => Ok(n as usize),
            _ => Err(BAD_LEN.into()),
        }
    }

    impl<'a> Value<'a> {
        /// Element `i` of an integer vector; `None` for other types or past `n`.
        pub fn int(&self, i: usize) -> Option<Int> {
            if i >= self.n {
                return None;
            }
            let d = self.data;
            let (v, missing) = match self.ty {
                Ty::Int8 => (i32::from(d[i] as i8), i32::from(i8::MIN)),
                Ty::Int16 => (
                    i32::from(i16::from_le_bytes([d[2 * i], d[2 * i + 1]])),
                    i32::from(i16::MIN),
                ),
                Ty::Int32 => (
                    i32::from_le_bytes([d[4 * i], d[4 * i + 1], d[4 * i + 2], d[4 * i + 3]]),
                    i32::MIN,
                ),
                _ => return None,
            };
            Some(if v == missing {
                Int::Missing
            } else if v == missing + 1 {
                Int::End
            } else {
                Int::Value(v)
            })
        }

        /// Element `i` of a float vector; `None` for other types or past `n`.
        pub fn float(&self, i: usize) -> Option<Float> {
            if self.ty != Ty::Float || i >= self.n {
                return None;
            }
            let d = &self.data[4 * i..4 * i + 4];
            Some(match u32::from_le_bytes([d[0], d[1], d[2], d[3]]) {
                0x7f80_0001 => Float::Missing,
                0x7f80_0002 => Float::End,
                bits => Float::Value(f32::from_bits(bits)),
            })
        }

        /// A character string without NUL padding (empty for a missing
        /// value); `None` for numeric types.
        pub fn chars(&self) -> Option<&'a [u8]> {
            match self.ty {
                Ty::Char => {
                    let end = self.data.iter().rposition(|&b| b != 0).map_or(0, |p| p + 1);
                    Some(&self.data[..end])
                }
                Ty::Missing => Some(&[]),
                _ => None,
            }
        }

        /// Append the value as VCF text (`-50,.` for an integer vector with a
        /// missing element). Returns false for a value-less entry (a flag).
        pub fn render(&self, out: &mut Vec<u8>) -> bool {
            fn list<T>(
                out: &mut Vec<u8>,
                items: impl Iterator<Item = Option<T>>,
                put: impl Fn(&mut Vec<u8>, T),
            ) {
                for (k, item) in items.enumerate() {
                    if k > 0 {
                        out.push(b',');
                    }
                    match item {
                        Some(v) => put(out, v),
                        None => out.push(b'.'),
                    }
                }
            }
            match self.ty {
                Ty::Missing => return false,
                Ty::Char => out.extend_from_slice(self.chars().unwrap_or_default()),
                Ty::Int8 | Ty::Int16 | Ty::Int32 => {
                    let items = (0..self.n).map_while(|i| match self.int(i)? {
                        Int::Value(v) => Some(Some(v)),
                        Int::Missing => Some(None),
                        Int::End => None,
                    });
                    list(out, items, |o, v| {
                        write!(o, "{v}").expect("write to Vec");
                    });
                }
                Ty::Float => {
                    let items = (0..self.n).map_while(|i| match self.float(i)? {
                        Float::Value(v) => Some(Some(v)),
                        Float::Missing => Some(None),
                        Float::End => None,
                    });
                    list(out, items, |o, v| {
                        write!(o, "{v}").expect("write to Vec");
                    });
                }
            }
            true
        }
    }
}

use typed::Int;

/// Per-reader decoding state; the buffers are reused for every record.
pub(crate) struct Decoder {
    span: VcfSpan,
    end_key: Option<i32>,
    svlen_key: Option<i32>,
    /// Contig dictionary index to `ChromTable` id.
    chrom_ids: Vec<Option<ChromId>>,
    alts: Vec<u8>,
    end: Vec<u8>,
    svlen: Vec<u8>,
}

impl Decoder {
    pub(crate) fn new(span: VcfSpan, dicts: &BcfDicts) -> Decoder {
        let key = |id: &[u8]| {
            dicts
                .strings
                .iter()
                .position(|s| s == id)
                .and_then(|i| i32::try_from(i).ok())
        };
        Decoder {
            span,
            end_key: key(b"END"),
            svlen_key: key(b"SVLEN"),
            chrom_ids: vec![None; dicts.contigs.len()],
            alts: Vec::new(),
            end: Vec::new(),
            svlen: Vec::new(),
        }
    }

    /// Decode a whole framed record (`l_shared`, `l_indiv`, shared, indiv).
    pub(crate) fn decode(
        &mut self,
        rec: &[u8],
        dicts: &BcfDicts,
        chroms: &mut ChromTable,
    ) -> Result<(ChromId, u64, Interval), String> {
        let l_shared = u32::from_le_bytes([rec[0], rec[1], rec[2], rec[3]]) as usize;
        let shared = &rec[8..8 + l_shared];
        if shared.len() < 24 {
            return Err(format!(
                "the shared section is {} bytes; a BCF record needs at least 24",
                shared.len()
            ));
        }
        let word = |o: usize| [shared[o], shared[o + 1], shared[o + 2], shared[o + 3]];
        let contig = i32::from_le_bytes(word(0));
        let pos = match i32::from_le_bytes(word(4)) {
            -1 => return Err(FootprintError::PosZero.to_string()),
            p if p < -1 => return Err(format!("POS {} is negative", i64::from(p) + 1)),
            p => p as u64 + 1,
        };
        let n_allele_info = u32::from_le_bytes(word(16));
        let (n_allele, n_info) = (n_allele_info >> 16, n_allele_info & 0xffff);

        let mut cur = &shared[24..];
        typed::read(&mut cur)?; // ID
        let mut reference: &[u8] = b"";
        self.alts.clear();
        for i in 0..n_allele {
            let a = typed::read(&mut cur)?
                .chars()
                .ok_or("an allele is not a character string")?;
            match i {
                0 => reference = a,
                1 => self.alts.extend_from_slice(a),
                _ => {
                    self.alts.push(b',');
                    self.alts.extend_from_slice(a);
                }
            }
        }
        if n_allele <= 1 {
            self.alts.push(b'.');
        }

        // Like the text reader, END and SVLEN only matter for symbolic ALTs.
        let (mut has_end, mut has_svlen) = (false, false);
        if memchr::memchr(b'<', &self.alts).is_some()
            && (self.end_key.is_some() || self.svlen_key.is_some())
        {
            self.end.clear();
            self.svlen.clear();
            typed::read(&mut cur)?; // FILTER
            for _ in 0..n_info {
                let k = typed::read(&mut cur)?;
                let key = match k.int(0) {
                    Some(Int::Value(key)) if k.n == 1 => key,
                    _ => return Err("an INFO key is not a single integer".into()),
                };
                let v = typed::read(&mut cur)?;
                if Some(key) == self.end_key && !has_end {
                    has_end = v.render(&mut self.end);
                } else if Some(key) == self.svlen_key && !has_svlen {
                    has_svlen = v.render(&mut self.svlen);
                }
            }
        }
        let end = if has_end {
            Some(coord("INFO/END", &self.end)?)
        } else {
            None
        };
        let span = footprint(
            &Variant {
                pos,
                reference,
                alts: &self.alts,
                end,
                svlen: has_svlen.then_some(&self.svlen[..]),
            },
            self.span,
        )
        .map_err(|e| e.to_string())?;
        let chrom = self.chrom(contig, dicts, chroms)?;
        Ok((chrom, pos - 1, span))
    }

    /// INFO/END of a framed record, whatever its ALT alleles (index spans
    /// use it for every record, D20). Call after a successful `decode`.
    pub(crate) fn info_end(&self, rec: &[u8]) -> Result<Option<u64>, String> {
        let Some(end_key) = self.end_key else {
            return Ok(None);
        };
        let l_shared = u32::from_le_bytes([rec[0], rec[1], rec[2], rec[3]]) as usize;
        let shared = &rec[8..8 + l_shared];
        let n_allele_info = u32::from_le_bytes([shared[16], shared[17], shared[18], shared[19]]);
        let (n_allele, n_info) = (n_allele_info >> 16, n_allele_info & 0xffff);
        let mut cur = &shared[24..];
        typed::read(&mut cur)?; // ID
        for _ in 0..n_allele {
            typed::read(&mut cur)?;
        }
        typed::read(&mut cur)?; // FILTER
        for _ in 0..n_info {
            let k = typed::read(&mut cur)?;
            let key = match k.int(0) {
                Some(Int::Value(key)) if k.n == 1 => key,
                _ => return Err("an INFO key is not a single integer".into()),
            };
            let v = typed::read(&mut cur)?;
            if key == end_key {
                return Ok(match v.int(0) {
                    Some(Int::Value(e)) if e >= 0 => Some(e as u64),
                    _ => None,
                });
            }
        }
        Ok(None)
    }

    fn chrom(
        &mut self,
        contig: i32,
        dicts: &BcfDicts,
        chroms: &mut ChromTable,
    ) -> Result<ChromId, String> {
        let i = usize::try_from(contig)
            .ok()
            .filter(|&i| dicts.contigs.get(i).is_some_and(|n| !n.is_empty()))
            .ok_or_else(|| {
                format!(
                    "CHROM index {contig} is not in the contig dictionary (the header declares {} contigs)",
                    dicts.contigs.iter().filter(|n| !n.is_empty()).count()
                )
            })?;
        let name = &dicts.contigs[i];
        // The cached id is checked against the table, so a reader can be
        // used with more than one table.
        if let Some(id) = self.chrom_ids[i]
            && (id.0 as usize) < chroms.len()
            && chroms.name(id) == &name[..]
        {
            return Ok(id);
        }
        let id = intern_chrom(chroms, name)?;
        self.chrom_ids[i] = Some(id);
        Ok(id)
    }
}

/// Append up to `n` bytes of `src` to `out`; returns how many were read
/// (fewer only at end of input). `out` grows with the data actually read, so
/// a corrupt length can't cause a huge allocation.
pub(crate) fn append(src: &mut dyn BufRead, out: &mut Vec<u8>, n: usize) -> io::Result<usize> {
    let mut left = n;
    while left > 0 {
        let avail = match src.fill_buf() {
            Ok(b) => b,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        };
        if avail.is_empty() {
            break;
        }
        let k = avail.len().min(left);
        out.extend_from_slice(&avail[..k]);
        src.consume(k);
        left -= k;
    }
    Ok(n - left)
}

fn dicts(h: &Header) -> &BcfDicts {
    match h {
        Header::Vcf(VcfHeader { bcf: Some(d), .. }) => d,
        _ => unreachable!("a BCF reader always has BCF dictionaries"),
    }
}

pub(crate) struct BcfReader {
    name: String,
    compression: Compression,
    header: Header,
    diag: Diagnostics,
    policy: InvalidPolicy,
    src: Box<dyn BufRead + Send>,
    /// Records started so far (the 1-based number of the current one).
    rec: u64,
    eof: bool,
    dec: Decoder,
}

impl BcfReader {
    pub(crate) fn open_bcf(
        input: Input,
        policy: InvalidPolicy,
        span: VcfSpan,
    ) -> Result<Self, IoError> {
        let Input {
            name,
            compression,
            mut reader,
            ..
        } = input;
        let h = read_bcf_header(&mut reader).map_err(|msg| IoError::Detect {
            msg: format!("{name}: {msg}"),
            path: name.clone(),
        })?;
        let header = Header::Vcf(h);
        let dec = Decoder::new(span, dicts(&header));
        Ok(BcfReader {
            name,
            compression,
            header,
            diag: Diagnostics::default(),
            policy,
            src: reader,
            rec: 0,
            eof: false,
            dec,
        })
    }

    fn err(&self, msg: impl Into<String>) -> IoError {
        IoError::Parse {
            path: self.name.clone(),
            loc: Loc::Record(self.rec),
            msg: msg.into(),
        }
    }

    fn read_err(&self, source: io::Error) -> IoError {
        IoError::Read {
            path: self.name.clone(),
            source,
        }
    }
}

pub(crate) const RECORD_TRUNCATED: &str = "record is truncated; the file was probably cut off";

impl RecordReader for BcfReader {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        Format::Bcf
    }
    fn compression(&self) -> Compression {
        self.compression
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        let mut n = 0;
        while n < max {
            if self.eof {
                return Ok(false);
            }
            let start = batch.buf.len();
            let got = append(&mut self.src, &mut batch.buf, 8).map_err(|e| self.read_err(e))?;
            if got == 0 {
                self.eof = true;
                return Ok(false);
            }
            self.rec += 1;
            if got < 8 {
                batch.buf.truncate(start);
                return Err(self.err(RECORD_TRUNCATED));
            }
            let b = &batch.buf[start..];
            let l_shared = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
            let l_indiv = u32::from_le_bytes([b[4], b[5], b[6], b[7]]);
            let body = u64::from(l_shared) + u64::from(l_indiv);
            if body + 8 > u64::from(u32::MAX) {
                batch.buf.truncate(start);
                return Err(self.err("record is larger than 4 GiB"));
            }
            let body = body as usize;
            let got = append(&mut self.src, &mut batch.buf, body).map_err(|e| self.read_err(e))?;
            if got < body {
                batch.buf.truncate(start);
                return Err(self.err(RECORD_TRUNCATED));
            }
            let len = batch.buf.len() - start;
            let decoded = self
                .dec
                .decode(&batch.buf[start..], dicts(&self.header), chroms);
            match decoded {
                Ok((chrom, key, span)) => {
                    batch.recs.push(RecordMeta {
                        chrom,
                        key,
                        span,
                        strand: Strand::Unknown,
                        loc: self.rec,
                        off: start,
                        len: len as u32,
                    });
                    n += 1;
                }
                Err(msg) => {
                    batch.buf.truncate(start);
                    match self.policy {
                        InvalidPolicy::Fail => return Err(self.err(msg)),
                        InvalidPolicy::Skip => self.diag.record_at(Loc::Record(self.rec), &msg),
                    }
                }
            }
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::typed::{Float, Int, Ty, read};
    use super::*;

    fn one(bytes: &[u8]) -> typed::Value<'_> {
        let mut src = bytes;
        let v = read(&mut src).unwrap();
        assert!(src.is_empty(), "consumed exactly one value");
        v
    }
    fn ints(v: typed::Value) -> Vec<Int> {
        (0..v.n).map(|i| v.int(i).unwrap()).collect()
    }
    fn rendered(v: typed::Value) -> Option<String> {
        let mut out = Vec::new();
        v.render(&mut out).then(|| String::from_utf8(out).unwrap())
    }

    #[test]
    fn int8_vector_with_missing() {
        let v = one(&[0x21, 0x01, 0x80]);
        assert_eq!((v.ty, v.n), (Ty::Int8, 2));
        assert_eq!(ints(v), [Int::Value(1), Int::Missing]);
        assert_eq!(rendered(v).unwrap(), "1,.");
        assert_eq!(ints(one(&[0x11, 0x81])), [Int::End]);
        assert_eq!(ints(one(&[0x11, 0xce])), [Int::Value(-50)]);
    }

    #[test]
    fn int16_end() {
        // END=10000 as htslib/noodles write it: the smallest int type that fits
        let v = one(&[0x12, 0x10, 0x27]);
        assert_eq!((v.ty, ints(v)), (Ty::Int16, vec![Int::Value(10000)]));
        assert_eq!(rendered(v).unwrap(), "10000");
        assert_eq!(ints(one(&[0x12, 0x18, 0xfc])), [Int::Value(-1000)]);
        assert_eq!(
            ints(one(&[0x22, 0x00, 0x80, 0x01, 0x80])),
            [Int::Missing, Int::End]
        );
    }

    #[test]
    fn int32_svlen_with_end_of_vector() {
        let mut b = vec![0x33];
        for x in [-50i32, i32::MIN, i32::MIN + 1] {
            b.extend_from_slice(&x.to_le_bytes());
        }
        let v = one(&b);
        assert_eq!(ints(v), [Int::Value(-50), Int::Missing, Int::End]);
        assert_eq!(rendered(v).unwrap(), "-50,.");
        let mut b = vec![0x33];
        for x in [100_000i32, i32::MIN + 1, i32::MIN + 1] {
            b.extend_from_slice(&x.to_le_bytes());
        }
        assert_eq!(rendered(one(&b)).unwrap(), "100000");
    }

    #[test]
    fn floats_chars_and_missing() {
        let mut b = vec![0x45];
        for bits in [
            1.5f32.to_bits(),
            (-50.0f32).to_bits(),
            0x7f80_0001,
            0x7f80_0002,
        ] {
            b.extend_from_slice(&bits.to_le_bytes());
        }
        let v = one(&b);
        assert_eq!(
            (0..4).map(|i| v.float(i).unwrap()).collect::<Vec<_>>(),
            [
                Float::Value(1.5),
                Float::Value(-50.0),
                Float::Missing,
                Float::End
            ]
        );
        assert_eq!(rendered(v).unwrap(), "1.5,-50,.");
        assert_eq!(v.int(0), None);
        let v = one(b"\x47ACG\x00");
        assert_eq!((v.ty, v.chars()), (Ty::Char, Some(&b"ACG"[..])));
        assert_eq!(rendered(v).unwrap(), "ACG");
        let v = one(&[0x07]);
        assert_eq!((v.n, v.chars()), (0, Some(&b""[..])));
        let v = one(&[0x00]);
        assert_eq!(
            (v.ty, v.chars(), rendered(v)),
            (Ty::Missing, Some(&b""[..]), None)
        );
    }

    #[test]
    fn overflow_length() {
        let mut b = vec![0xf7, 0x11, 16];
        b.extend_from_slice(b"0123456789abcdef");
        let v = one(&b);
        assert_eq!(
            (v.ty, v.n, v.chars()),
            (Ty::Char, 16, Some(&b"0123456789abcdef"[..]))
        );
        let mut b = vec![0xf1, 0x12, 0x2c, 0x01];
        b.extend(std::iter::repeat_n(7u8, 300));
        let v = one(&b);
        assert_eq!(
            (v.ty, v.n, v.int(299)),
            (Ty::Int8, 300, Some(Int::Value(7)))
        );
        // the next value starts right after
        let mut src = &[
            0xf1, 0x11, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0x11, 9,
        ][..];
        assert_eq!(read(&mut src).unwrap().n, 15);
        assert_eq!(ints(read(&mut src).unwrap()), [Int::Value(9)]);
    }

    #[test]
    fn malformed_values_are_errors() {
        let err = |b: &[u8]| read(&mut &b[..]).unwrap_err();
        assert!(err(&[]).contains("past the end"));
        assert!(err(&[0x31, 1, 2]).contains("past the end"));
        assert!(err(&[0x14]).contains("unknown BCF value type 4"));
        assert!(err(&[0xf1, 0x21, 1, 1]).contains("invalid length")); // count 2
        assert!(err(&[0xf1, 0xf1, 0x11, 1]).contains("invalid length")); // nested overflow
        assert!(err(&[0xf1, 0x11, 0xff]).contains("invalid length")); // negative
        assert!(err(&[0xf1, 0x17, b'a']).contains("invalid length")); // not an int
        let mut huge = vec![0xf3, 0x13];
        huge.extend_from_slice(&i32::MAX.to_le_bytes());
        assert!(err(&huge).contains("past the end"));
    }

    #[test]
    fn header_fields_respect_quotes() {
        let f = header_fields(br#"ID=DP,Description="a, b \"c,d\" IDX=9",IDX=3>"#);
        assert_eq!(
            f,
            [
                (&b"ID"[..], &b"DP"[..]),
                (b"Description", br#"a, b \"c,d\" IDX=9"#),
                (b"IDX", b"3"),
            ]
        );
    }

    #[test]
    fn dictionaries_follow_idx_and_dedupe() {
        let d = BcfDicts::from_header(
            b"##fileformat=VCFv4.3\n\
##FILTER=<ID=q10,Description=\"x\">\n\
##INFO=<ID=DP,Description=\"x, IDX=7\">\n\
##FORMAT=<ID=DP,Description=\"x\">\n\
##FILTER=<ID=PASS,Description=\"All\">\n\
##contig=<ID=b,IDX=2>\n\
##contig=<ID=a,IDX=0>\n\
#CHROM\n",
        )
        .unwrap();
        assert_eq!(d.strings, [&b"PASS"[..], b"q10", b"DP"].map(<[u8]>::to_vec));
        assert_eq!(d.contigs, [&b"a"[..], b"", b"b"].map(<[u8]>::to_vec));
        let e =
            BcfDicts::from_header(b"##INFO=<ID=DP,IDX=1>\n##FORMAT=<ID=DP,IDX=2>\n").unwrap_err();
        assert!(e.contains("IDX=2, but DP already has IDX=1"), "{e}");
        let e = BcfDicts::from_header(b"##contig=<ID=a,IDX=x>\n").unwrap_err();
        assert!(e.contains("invalid IDX"), "{e}");
        let e = BcfDicts::from_header(b"##contig=<ID=a,IDX=99999999999>\n").unwrap_err();
        assert!(e.contains("invalid IDX"), "{e}");
    }
}
