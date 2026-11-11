//! VCF <-> BCF conversion (`view -O bcf` / `-O vcf`).
//!
//! The output header is the input's own header text with the
//! `##glociCommand` line inserted, never a header serialised by noodles:
//! noodles drops `IDX=` attributes when it writes a header, although BCF
//! records are encoded against the IDX-aware dictionaries. noodles is only
//! used to encode and decode records, against a header parsed from the same
//! text whose dictionary positions are set from gloci's own dictionaries
//! (`BcfDicts`), so the header text and the records always agree.

use crate::bcf::{BcfDicts, typed, write_bcf_header};
use crate::bcf_fields::{Def, Defs, FieldEncoder, Kind};
use crate::detect::Format;
use crate::error::{IoError, Loc};
use crate::output::{OutFormat, Output};
use crate::record::{Header, RecordBatch, RecordReader, next_batch_size};
use crate::sorted::OrderGuard;
use crate::vcf::VcfHeader;
use crate::writer::WriterOptions;
use gloci_core::ChromTable;
use noodles_vcf::header::StringMaps;
use noodles_vcf::variant::io::Write as _;
use std::collections::HashMap;
use std::io::Write;

/// Which parts of the output `convert_variants` writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeaderMode {
    /// Header, then records.
    Full,
    /// Records only (text output).
    None,
    /// The header only.
    Only,
}

fn internal(path: &str, msg: String) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg,
    }
}

/// Parse `h.raw` with noodles and give every contig, INFO, FILTER and FORMAT
/// the dictionary position gloci's dictionaries assign it.
fn noodles_header(name: &str, h: &VcfHeader) -> Result<noodles_vcf::Header, IoError> {
    let bad = |what: &str, e: &dyn std::fmt::Display| {
        internal(
            name,
            format!(
                "{name}: the VCF header cannot be used for conversion ({what}: {e}); fix the header lines it names"
            ),
        )
    };
    let text = std::str::from_utf8(&h.raw).map_err(|e| bad("not UTF-8", &e))?;
    let mut header: noodles_vcf::Header = text.parse().map_err(|e| bad("parse", &e))?;
    let dicts = match &h.bcf {
        Some(d) => d.clone(),
        None => {
            BcfDicts::from_header(&h.raw).map_err(|e| internal(name, format!("{name}: {e}")))?
        }
    };
    let positions = |names: &[Vec<u8>]| -> HashMap<String, usize> {
        names
            .iter()
            .enumerate()
            .filter(|(_, n)| !n.is_empty())
            .map(|(i, n)| (String::from_utf8_lossy(n).into_owned(), i))
            .collect()
    };
    let (contigs, strings) = (positions(&dicts.contigs), positions(&dicts.strings));
    for (id, c) in header.contigs_mut() {
        *c.idx_mut() = contigs.get(id.as_str()).copied();
    }
    for (id, i) in header.infos_mut() {
        *i.idx_mut() = strings.get(id.as_str()).copied();
    }
    for (id, f) in header.filters_mut() {
        *f.idx_mut() = strings.get(id.as_str()).copied();
    }
    for (id, f) in header.formats_mut() {
        *f.idx_mut() = strings.get(id.as_str()).copied();
    }
    *header.string_maps_mut() =
        StringMaps::try_from(&header).map_err(|e| bad("dictionaries", &e))?;
    Ok(header)
}

enum Encoder {
    /// VCF text records to BCF: noodles encodes into `buf`, with the string
    /// maps of the header.
    ToBcf(noodles_bcf::io::Writer<Vec<u8>>),
    /// BCF records to VCF text lines.
    ToVcf(noodles_vcf::io::Writer<Vec<u8>>),
}

/// Converts the records of batches between VCF and BCF.
struct Converter {
    name: String,
    from: Format,
    header: noodles_vcf::Header,
    enc: Encoder,
    /// The INFO and FORMAT fields of the header, for the native encoding
    /// of the record's INFO and sample columns.
    defs: Defs,
    n_samples: usize,
    /// String-dictionary index of FORMAT/GT, when declared.
    gt_key: Option<usize>,
    fields: FieldEncoder,
    /// The record's shared part (without INFO) and its INFO and sample
    /// parts, reused between records.
    shared: Vec<u8>,
    info: Vec<u8>,
    indiv: Vec<u8>,
    framed: Vec<u8>,
    widened: Vec<u8>,
    /// The record with GT moved first, when it was not (noodles' VCF writer needs it there).
    gt_moved: Vec<u8>,
}

impl Converter {
    fn new(name: &str, h: &VcfHeader, from: Format) -> Result<Converter, IoError> {
        let header = noodles_header(name, h)?;
        let enc = match from {
            Format::Vcf => {
                let mut w = noodles_bcf::io::Writer::from(Vec::new());
                // Sets the writer's string maps; the header bytes it writes
                // are discarded (our own header text is written instead).
                w.write_header(&header)
                    .map_err(|e| internal(name, format!("{name}: {e}")))?;
                w.get_mut().clear();
                Encoder::ToBcf(w)
            }
            _ => Encoder::ToVcf(noodles_vcf::io::Writer::new(Vec::new())),
        };
        let dicts = match &h.bcf {
            Some(d) => d.clone(),
            None => {
                BcfDicts::from_header(&h.raw).map_err(|e| internal(name, format!("{name}: {e}")))?
            }
        };
        let idx: HashMap<&[u8], usize> = dicts
            .strings
            .iter()
            .enumerate()
            .filter(|(_, n)| !n.is_empty())
            .map(|(i, n)| (n.as_slice(), i))
            .collect();
        let mut defs = Defs::default();
        {
            use noodles_vcf::header::record::value::map::{format, info};
            for (id, f) in header.formats() {
                let kind = match (id.as_str(), f.ty()) {
                    ("GT", _) => Kind::Gt,
                    (_, format::Type::Integer) => Kind::Int,
                    (_, format::Type::Float) => Kind::Float,
                    _ => Kind::Str,
                };
                if let Some(&i) = idx.get(id.as_bytes()) {
                    defs.format
                        .insert(id.as_bytes().to_vec(), Def { idx: i, kind });
                }
            }
            for (id, f) in header.infos() {
                let kind = match f.ty() {
                    info::Type::Flag => Kind::Flag,
                    info::Type::Integer => Kind::Int,
                    info::Type::Float => Kind::Float,
                    _ => Kind::Str,
                };
                if let Some(&i) = idx.get(id.as_bytes()) {
                    defs.info
                        .insert(id.as_bytes().to_vec(), Def { idx: i, kind });
                }
            }
        }
        let n_samples = header.sample_names().len();
        let defs_gt = defs.format.get(&b"GT"[..]).map(|d| d.idx);
        Ok(Converter {
            name: name.to_string(),
            from,
            header,
            enc,
            defs,
            n_samples,
            gt_key: defs_gt,
            fields: FieldEncoder::new(name),
            shared: Vec::new(),
            info: Vec::new(),
            indiv: Vec::new(),
            framed: Vec::new(),
            widened: Vec::new(),
            gt_moved: Vec::new(),
        })
    }

    /// Encode one VCF text record (no line terminator) as BCF into `out`.
    fn encode_text(&mut self, raw: &[u8], loc: Loc, out: &mut dyn Write) -> Result<(), IoError> {
        let fail = |e: std::io::Error| IoError::Parse {
            path: self.name.clone(),
            loc,
            msg: convert_failure(&e),
        };
        let Encoder::ToBcf(w) = &mut self.enc else {
            return Err(internal(
                &self.name,
                "internal error: not a BCF encoder".into(),
            ));
        };
        // BCF needs every contig in the header's dictionary.
        let chrom = raw.split(|&b| b == b'\t').next().unwrap_or_default();
        let declared =
            std::str::from_utf8(chrom).is_ok_and(|c| self.header.contigs().contains_key(c));
        if !declared {
            let c = String::from_utf8_lossy(chrom);
            return Err(IoError::Parse {
                path: self.name.clone(),
                loc,
                msg: format!(
                    "chromosome \"{c}\" is not in the ##contig lines, and BCF needs every contig in its header; \
add ##contig=<ID={c}> to the header"
                ),
            });
        }
        // bcftools refuses a FILTER missing from the header (BCF has no
        // place for it); so does gloci, and says what to add.
        if let Some(id) = undeclared_filter(raw, &self.header) {
            return Err(IoError::Parse {
                path: self.name.clone(),
                loc,
                msg: format!(
                    "FILTER \"{id}\" is not declared in the header, and BCF needs every FILTER in its header; \
add ##FILTER=<ID={id},Description=\"...\"> to the header"
                ),
            });
        }
        // noodles encodes the fixed part of the record (CHROM to FILTER); the
        // INFO and sample columns are encoded here (see `bcf_fields`).
        let mut cols = raw.split(|&b| b == b'\t');
        let fixed_len: usize = cols.by_ref().take(7).map(|c| c.len() + 1).sum();
        let info_col = cols.next().unwrap_or(b".");
        let mut fixed = Vec::with_capacity(fixed_len + 2);
        fixed.extend_from_slice(&raw[..fixed_len.min(raw.len())]);
        fixed.extend_from_slice(b".");
        let rec = noodles_vcf::Record::try_from(&fixed[..]).map_err(fail)?;
        if let Err(e) = w.write_variant_record(&self.header, &rec) {
            w.get_mut().clear();
            return Err(fail(e));
        }
        let parse = |msg: String| IoError::Parse {
            path: self.name.clone(),
            loc,
            msg,
        };
        let framed = w.get_ref();
        let l_shared = framed
            .get(..4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize);
        let (Some(l_shared), true) = (l_shared, framed.len() >= 8 + 24) else {
            w.get_mut().clear();
            return Err(internal(
                &self.name,
                "internal error: noodles wrote a malformed BCF record".into(),
            ));
        };
        self.shared.clear();
        self.shared.extend_from_slice(&framed[8..8 + l_shared]);
        w.get_mut().clear();
        self.info.clear();
        self.indiv.clear();
        self.fields.at(loc);
        let (n_info, end) = self
            .fields
            .info(&self.defs, info_col, &mut self.info)
            .map_err(parse)?;
        let n_fmt = self
            .fields
            .samples(&self.defs, raw, self.n_samples, &mut self.indiv)
            .map_err(|msg| IoError::Parse {
                path: self.name.clone(),
                loc,
                msg,
            })?;
        if n_info > 0xffff || n_fmt > 0xff {
            return Err(IoError::Parse {
                path: self.name.clone(),
                loc,
                msg:
                    "the record has too many INFO or FORMAT fields for BCF (65535 and 255 at most)"
                        .into(),
            });
        }
        // rlen is END - POS + 1 when INFO has a single END of at least POS
        // (bcftools 1.21 ignores SVLEN and any END below POS).
        if let Some(end) = end {
            let pos = raw
                .split(|&b| b == b'\t')
                .nth(1)
                .and_then(|p| std::str::from_utf8(p).ok())
                .and_then(|p| p.parse::<i64>().ok());
            match pos {
                Some(pos) if end >= pos => {
                    let rlen = i32::try_from(end - pos + 1).unwrap_or(i32::MAX);
                    self.shared[8..12].copy_from_slice(&rlen.to_le_bytes());
                }
                _ => self.fields.end_before_pos(end),
            }
        }
        let allele_info = u32::from_le_bytes(self.shared[16..20].try_into().unwrap());
        self.shared[16..20]
            .copy_from_slice(&((allele_info & 0xffff_0000) | n_info as u32).to_le_bytes());
        let fmt_sample = ((n_fmt as u32) << 24) | self.n_samples as u32;
        self.shared[20..24].copy_from_slice(&fmt_sample.to_le_bytes());
        self.framed.clear();
        let l_shared = self.shared.len() + self.info.len();
        self.framed
            .extend_from_slice(&(l_shared as u32).to_le_bytes());
        self.framed
            .extend_from_slice(&(self.indiv.len() as u32).to_le_bytes());
        self.framed.extend_from_slice(&self.shared);
        self.framed.extend_from_slice(&self.info);
        self.framed.extend_from_slice(&self.indiv);
        out.write_all(&self.framed)
            .map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })
    }

    fn convert(&mut self, batch: &RecordBatch, out: &mut dyn Write) -> Result<(), IoError> {
        for (i, m) in batch.recs.iter().enumerate() {
            let loc = if self.from == Format::Bcf {
                Loc::Record(m.loc)
            } else {
                Loc::Line(m.loc)
            };
            let fail = |e: std::io::Error| IoError::Parse {
                path: self.name.clone(),
                loc,
                msg: convert_failure(&e),
            };
            let raw = batch.raw(i);
            let werr = |source| IoError::Write {
                path: "output".to_string(),
                source,
            };
            match &mut self.enc {
                Encoder::ToBcf(_) => self.encode_text(raw, loc, out)?,
                Encoder::ToVcf(w) => {
                    let mut rec = noodles_bcf::Record::default();
                    let moved = gt_first(raw, self.gt_key, &mut self.gt_moved).map_err(fail)?;
                    let raw: &[u8] = if moved.is_some() { &self.gt_moved } else { raw };
                    let raw =
                        fill_empty_fields(raw, self.gt_key, &mut self.widened).map_err(fail)?;
                    noodles_bcf::io::Reader::from(raw)
                        .read_record(&mut rec)
                        .map_err(fail)?;
                    w.write_variant_record(&self.header, &rec).map_err(fail)?;
                    let line = w.get_mut();
                    if let Some(p) = moved {
                        restore_gt_position(line, p);
                    }
                    if line.contains(&b'%') {
                        let had_eol = line.last() == Some(&b'\n');
                        if had_eol {
                            line.pop();
                        }
                        tidy_line(line);
                        if had_eol {
                            line.push(b'\n');
                        }
                    }
                    out.write_all(w.get_ref()).map_err(werr)?;
                    w.get_mut().clear();
                }
            }
        }
        Ok(())
    }
}

/// Renders the records of a BCF input as VCF text lines, as `view -O vcf`
/// writes them, for commands that print records in a TSV.
pub struct BcfText {
    name: String,
    conv: Converter,
    line: Vec<u8>,
}

impl BcfText {
    /// A renderer for `reader` when it is BCF; `None` for any other
    /// format, whose raw records are already text.
    pub fn for_reader(reader: &dyn RecordReader) -> Result<Option<BcfText>, IoError> {
        let (Format::Bcf, Header::Vcf(h)) = (reader.format(), reader.header()) else {
            return Ok(None);
        };
        Ok(Some(BcfText {
            name: reader.name().to_string(),
            conv: Converter::new(reader.name(), h, Format::Bcf)?,
            line: Vec::new(),
        }))
    }

    /// A renderer for BCF records encoded against the header `h`.
    pub fn for_header(name: &str, h: &VcfHeader) -> Result<BcfText, IoError> {
        Ok(BcfText {
            name: name.to_string(),
            conv: Converter::new(name, h, Format::Bcf)?,
            line: Vec::new(),
        })
    }

    /// The VCF line of the BCF record `raw` (record number `loc`), without
    /// a newline. Valid until the next call.
    pub fn render(&mut self, raw: &[u8], loc: u64) -> Result<&[u8], IoError> {
        let fail = |e: std::io::Error| IoError::Parse {
            path: self.name.clone(),
            loc: Loc::Record(loc),
            msg: format!("cannot render the record as VCF text: {e}"),
        };
        let Encoder::ToVcf(w) = &mut self.conv.enc else {
            return Err(internal(
                &self.name,
                "internal error: not a BCF renderer".into(),
            ));
        };
        let mut rec = noodles_bcf::Record::default();
        let moved = gt_first(raw, self.conv.gt_key, &mut self.conv.gt_moved).map_err(fail)?;
        let raw: &[u8] = if moved.is_some() {
            &self.conv.gt_moved
        } else {
            raw
        };
        let raw = fill_empty_fields(raw, self.conv.gt_key, &mut self.conv.widened).map_err(fail)?;
        noodles_bcf::io::Reader::from(raw)
            .read_record(&mut rec)
            .map_err(fail)?;
        w.write_variant_record(&self.conv.header, &rec)
            .map_err(fail)?;
        self.line.clear();
        self.line.extend_from_slice(w.get_ref());
        w.get_mut().clear();
        while self.line.last() == Some(&b'\n') {
            self.line.pop();
        }
        if let Some(p) = moved {
            restore_gt_position(&mut self.line, p);
        }
        tidy_line(&mut self.line);
        Ok(&self.line)
    }
}

/// The message for a record that failed to convert: the declared-in-the-header
/// hint fits every cause but an unsupported one (the 16-bit GT).
fn convert_failure(e: &std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::Unsupported {
        format!("cannot convert the record: {e}")
    } else {
        format!(
            "cannot convert the record: {e}; check that every contig, INFO, FILTER and FORMAT field it uses is declared in the header"
        )
    }
}

/// Stands for an empty String value while noodles decodes the record.
const EMPTY_CELL: u8 = 0x01;

/// Undoes what noodles' VCF writer does to two values: the INFO string `.`
/// (written `%2E`, though it is the missing value) and the empty-value
/// placeholder of `fill_empty_fields` (written `%01`), a whole sample cell
/// that bcftools prints as nothing.
fn tidy_line(line: &mut Vec<u8>) {
    if !line.contains(&b'%') {
        return;
    }
    let mut out = Vec::with_capacity(line.len());
    for (i, col) in line.split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        if i == 7 {
            for (j, item) in col.split(|&b| b == b';').enumerate() {
                if j > 0 {
                    out.push(b';');
                }
                match item.strip_suffix(b"=%2E") {
                    Some(key) => {
                        out.extend_from_slice(key);
                        out.extend_from_slice(b"=.");
                    }
                    None => out.extend_from_slice(item),
                }
            }
        } else if i >= 9 {
            for (j, cell) in col.split(|&b| b == b':').enumerate() {
                if j > 0 {
                    out.push(b':');
                }
                if cell != b"%01" {
                    out.extend_from_slice(cell);
                }
            }
        } else {
            out.extend_from_slice(col);
        }
    }
    *line = out;
}

/// htslib writes a FORMAT key that no sample has a value for (every sample
/// cut short before it) as a zero-width field, which noodles cannot decode, and an empty String value
/// as NULs, which it refuses. Returns `raw` with each zero-width field
/// replaced by one missing value per sample and each empty value by `.`
/// (written to `buf`), or `raw` itself when there is none.
fn fill_empty_fields<'a>(
    raw: &'a [u8],
    gt_key: Option<usize>,
    buf: &'a mut Vec<u8>,
) -> std::io::Result<&'a [u8]> {
    let bad = || std::io::Error::new(std::io::ErrorKind::InvalidData, "malformed BCF record");
    let u32_at = |at: usize| -> std::io::Result<usize> {
        raw.get(at..at + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
            .ok_or_else(bad)
    };
    let (l_shared, l_indiv) = (u32_at(0)?, u32_at(4)?);
    let indiv_at = 8 + l_shared;
    let n_fmt = u32_at(8 + 20)? >> 24;
    let n_samples = u32_at(8 + 20)? & 0x00ff_ffff;
    let indiv = raw.get(indiv_at..indiv_at + l_indiv).ok_or_else(bad)?;
    let mut cur = indiv;
    let mut body: Option<Vec<u8>> = None;
    for _ in 0..n_fmt {
        let start = indiv.len() - cur.len();
        let key = typed::read(&mut cur).map_err(|_| bad())?;
        let is_gt = matches!(key.int(0), Some(typed::Int::Value(v)) if Some(v as usize) == gt_key);
        let desc = indiv.len() - cur.len();
        let (ty, n) = typed::read_descriptor(&mut cur).map_err(|_| bad())?;
        if is_gt && matches!(ty, typed::Ty::Int16 | typed::Ty::Int32) {
            // noodles-bcf panics on a GT wider than 8 bits.
            return Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "the genotype uses an allele index of 63 or more (a 16-bit GT), which gloci \
cannot render as VCF text; `gloci vcf stats` reads such a record",
            ));
        }
        let len = n
            .checked_mul(ty.width())
            .and_then(|l| l.checked_mul(n_samples))
            .filter(|&l| l <= cur.len())
            .ok_or_else(bad)?;
        cur = &cur[len..];
        let end = indiv.len() - cur.len();
        if n == 0 && n_samples > 0 {
            let b = body.get_or_insert_with(|| indiv[..start].to_vec());
            b.extend_from_slice(&indiv[start..desc]);
            // One missing value per sample, of the type the field has (GT and
            // the others are decoded by their declared type).
            match ty {
                typed::Ty::Char => {
                    b.push(0x17);
                    b.extend(std::iter::repeat_n(b'.', n_samples));
                }
                typed::Ty::Float => {
                    b.push(0x15);
                    for _ in 0..n_samples {
                        b.extend_from_slice(&0x7F80_0001u32.to_le_bytes());
                    }
                }
                // GT: the missing allele (a GT field is never zero-width
                // for any other type).
                typed::Ty::Missing => {
                    b.push(0x11);
                    b.extend(std::iter::repeat_n(0, n_samples));
                }
                _ => {
                    b.push(0x11);
                    b.extend(std::iter::repeat_n(0x80, n_samples));
                }
            }
        } else if is_gt
            && ty == typed::Ty::Int8
            && n > 0
            && indiv[end - n * n_samples..end]
                .chunks(n)
                .any(|c| c[0] == 0x80)
        {
            // A sample cut short before its GT has the integer-missing value
            // there; noodles reads it as an allele, so it is made the
            // missing allele.
            let b = body.get_or_insert_with(|| indiv[..start].to_vec());
            b.extend_from_slice(&indiv[start..end]);
            let at = b.len() - n * n_samples;
            for cell in b[at..].chunks_mut(n) {
                if cell[0] == 0x80 {
                    cell[0] = 0;
                }
            }
        } else if ty == typed::Ty::Char
            && n > 0
            && indiv[end - n * n_samples..end]
                .chunks(n)
                .any(|c| c.iter().all(|&b| b == 0))
        {
            // An empty string is all NULs, which noodles refuses; it decodes
            // a placeholder that `tidy_line` turns back into an empty value.
            let b = body.get_or_insert_with(|| indiv[..start].to_vec());
            b.extend_from_slice(&indiv[start..end]);
            let at = b.len() - n * n_samples;
            for cell in b[at..].chunks_mut(n) {
                if cell.iter().all(|&x| x == 0) {
                    cell[0] = EMPTY_CELL;
                }
            }
        } else if let Some(b) = body.as_mut() {
            b.extend_from_slice(&indiv[start..end]);
        }
    }
    let Some(body) = body else {
        return Ok(raw);
    };
    buf.clear();
    buf.extend_from_slice(&raw[..4]);
    buf.extend_from_slice(&(body.len() as u32).to_le_bytes());
    buf.extend_from_slice(&raw[8..indiv_at]);
    buf.extend_from_slice(&body);
    Ok(buf)
}

/// When GT is a FORMAT field of the BCF record `raw` but not the first, write
/// the record to `buf` with GT first, which noodles' VCF writer insists on,
/// and return GT's original position (see [`restore_gt_position`]).
fn gt_first(
    raw: &[u8],
    gt_key: Option<usize>,
    buf: &mut Vec<u8>,
) -> std::io::Result<Option<usize>> {
    let bad = || std::io::Error::new(std::io::ErrorKind::InvalidData, "malformed BCF record");
    let Some(gt_key) = gt_key else {
        return Ok(None);
    };
    let u32_at = |at: usize| -> std::io::Result<usize> {
        raw.get(at..at + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
            .ok_or_else(bad)
    };
    let (l_shared, l_indiv) = (u32_at(0)?, u32_at(4)?);
    let indiv_at = 8 + l_shared;
    let n_fmt = u32_at(8 + 20)? >> 24;
    let n_samples = u32_at(8 + 20)? & 0x00ff_ffff;
    let indiv = raw.get(indiv_at..indiv_at + l_indiv).ok_or_else(bad)?;
    let mut cur = indiv;
    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let mut gt_at = None;
    for k in 0..n_fmt {
        let start = indiv.len() - cur.len();
        let key = typed::read(&mut cur).map_err(|_| bad())?;
        if matches!(key.int(0), Some(typed::Int::Value(v)) if v as usize == gt_key) {
            gt_at = Some(k);
        }
        let (ty, n) = typed::read_descriptor(&mut cur).map_err(|_| bad())?;
        let len = n
            .checked_mul(ty.width())
            .and_then(|l| l.checked_mul(n_samples))
            .filter(|&l| l <= cur.len())
            .ok_or_else(bad)?;
        cur = &cur[len..];
        blocks.push((start, indiv.len() - cur.len()));
    }
    let Some(p) = gt_at.filter(|&p| p > 0) else {
        return Ok(None);
    };
    buf.clear();
    buf.extend_from_slice(&raw[..indiv_at]);
    buf.extend_from_slice(&indiv[blocks[p].0..blocks[p].1]);
    for (k, &(a, b)) in blocks.iter().enumerate() {
        if k != p {
            buf.extend_from_slice(&indiv[a..b]);
        }
    }
    buf.extend_from_slice(&indiv[blocks.last().map_or(0, |b| b.1)..]);
    Ok(Some(p))
}

/// Put the first `:` field of the FORMAT column and of every sample back at
/// index `p` (GT was moved first by [`gt_first`]) in the VCF line `line`.
fn restore_gt_position(line: &mut Vec<u8>, p: usize) {
    let had_eol = line.last() == Some(&b'\n');
    let body_end = line.len() - usize::from(had_eol);
    let mut out = Vec::with_capacity(line.len());
    for (i, col) in line[..body_end].split(|&b| b == b'\t').enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        if i < 8 {
            out.extend_from_slice(col);
            continue;
        }
        let mut fields: Vec<&[u8]> = col.split(|&b| b == b':').collect();
        let first = fields.remove(0);
        fields.insert(p.min(fields.len()), first);
        for (k, f) in fields.iter().enumerate() {
            if k > 0 {
                out.push(b':');
            }
            out.extend_from_slice(f);
        }
    }
    if had_eol {
        out.push(b'\n');
    }
    *line = out;
}

/// The first FILTER of the text record `raw` that the header does not
/// declare (`PASS` and `.` never need a declaration).
fn undeclared_filter(raw: &[u8], header: &noodles_vcf::Header) -> Option<String> {
    let col = raw.split(|&b| b == b'\t').nth(6)?;
    if col == b"." {
        return None;
    }
    col.split(|&b| b == b';')
        .filter(|id| !id.is_empty() && *id != b"PASS")
        .filter_map(|id| std::str::from_utf8(id).ok())
        .find(|id| !header.filters().contains_key(*id))
        .map(str::to_string)
}

/// Encodes VCF text records as BCF records against a header, for commands
/// that rewrite records as text and write BCF.
pub struct TextToBcf {
    conv: Converter,
}

impl TextToBcf {
    /// An encoder for records that use the dictionaries of `h`. The BCF
    /// header text to write is `h.raw`, via [`write_bcf_header`].
    pub fn new(name: &str, h: &VcfHeader) -> Result<TextToBcf, IoError> {
        Ok(TextToBcf {
            conv: Converter::new(name, h, Format::Vcf)?,
        })
    }

    /// Append the BCF record for the VCF line `line` (no terminator; record
    /// number `loc`, for messages) to `out`.
    pub fn encode(&mut self, line: &[u8], loc: u64, out: &mut Vec<u8>) -> Result<(), IoError> {
        self.conv.encode_text(line, Loc::Record(loc), out)
    }
}

/// Convert `reader`'s records (VCF to BCF or BCF to VCF) into `out`, in the
/// format `to`, and publish the output. Returns the number of records
/// written. `guard` checks the input's order batch by batch. Write errors
/// name the placeholder path `output`; the caller names the destination.
pub fn convert_variants(
    reader: &mut dyn RecordReader,
    chroms: &mut ChromTable,
    mut out: Output,
    to: OutFormat,
    o: &WriterOptions,
    mode: HeaderMode,
    mut guard: Option<&mut OrderGuard>,
) -> Result<u64, IoError> {
    let from = reader.format();
    let name = reader.name().to_string();
    match (from, to) {
        (Format::Vcf, OutFormat::Bcf) | (Format::Bcf, OutFormat::Vcf) => {}
        _ => {
            return Err(internal(
                &name,
                format!("internal error: cannot convert {from:?} to {to:?}"),
            ));
        }
    }
    let Header::Vcf(h) = reader.header().clone() else {
        return Err(internal(
            &name,
            "internal error: a VCF/BCF input has no VCF header".into(),
        ));
    };
    let werr = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    let mut conv = if mode == HeaderMode::Only {
        None
    } else {
        Some(Converter::new(&name, &h, from)?)
    };
    if mode != HeaderMode::None {
        let text = match &o.command_line {
            Some(cmd) => h.with_command_line(cmd, o.version),
            None => h.raw.clone(),
        };
        if to == OutFormat::Bcf {
            write_bcf_header(&mut out, &text).map_err(werr)?;
        } else {
            out.write_all(&text).map_err(werr)?;
        }
    }
    let mut n = 0u64;
    if let Some(conv) = conv.as_mut() {
        let mut batch = RecordBatch::default();
        let mut max = 1;
        loop {
            batch.clear();
            let more = reader.read_batch(chroms, &mut batch, max)?;
            if let Some(g) = guard.as_deref_mut() {
                g.check(chroms, &batch)?;
            }
            conv.convert(&batch, &mut out)?;
            n += batch.len() as u64;
            max = next_batch_size(&batch, max);
            if !more {
                break;
            }
        }
    }
    out.finish()?;
    Ok(n)
}
