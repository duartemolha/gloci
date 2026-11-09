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
use crate::bcf_encode::encode_gt;
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
    /// Dictionary index of the FORMAT key `GT`, when the header declares it.
    gt_key: Option<usize>,
    /// The record with its GT values blanked, reused between records.
    scratch: Vec<u8>,
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
        let gt_key = match &h.bcf {
            Some(d) => d.strings.iter().position(|n| n == b"GT"),
            None => BcfDicts::from_header(&h.raw)
                .map_err(|e| internal(name, format!("{name}: {e}")))?
                .strings
                .iter()
                .position(|n| n == b"GT"),
        };
        Ok(Converter {
            name: name.to_string(),
            from,
            header,
            enc,
            gt_key,
            scratch: Vec::new(),
        })
    }

    /// Encode one VCF text record (no line terminator) as BCF into `out`.
    fn encode_text(&mut self, raw: &[u8], loc: Loc, out: &mut dyn Write) -> Result<(), IoError> {
        let fail = |e: std::io::Error| IoError::Parse {
            path: self.name.clone(),
            loc,
            msg: format!(
                "cannot convert the record: {e}; check that every contig, INFO, FILTER and FORMAT field it uses is declared in the header"
            ),
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
        // noodles cannot encode GT vectors of different lengths, so it
        // encodes the record with every GT blanked to `./.` and the real GT
        // field is spliced in afterwards (see `splice_gt`).
        let gt = match self.gt_key {
            Some(_) => blank_gt(raw, &mut self.scratch),
            None => None,
        };
        let src = if gt.is_some() { &self.scratch[..] } else { raw };
        let rec = noodles_vcf::Record::try_from(src).map_err(fail)?;
        if let Err(e) = w.write_variant_record(&self.header, &rec) {
            w.get_mut().clear();
            return Err(fail(e));
        }
        match (gt, self.gt_key) {
            (Some(gt), Some(key)) => {
                let spliced = splice_gt(w.get_ref(), key, raw, gt);
                w.get_mut().clear();
                let bytes = spliced.map_err(|msg| IoError::Parse {
                    path: self.name.clone(),
                    loc,
                    msg,
                })?;
                out.write_all(&bytes).map_err(|source| IoError::Write {
                    path: "output".to_string(),
                    source,
                })?;
            }
            _ => {
                out.write_all(w.get_ref())
                    .map_err(|source| IoError::Write {
                        path: "output".to_string(),
                        source,
                    })?;
                w.get_mut().clear();
            }
        }
        Ok(())
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
                msg: format!(
                    "cannot convert the record: {e}; check that every contig, INFO, FILTER and FORMAT field it uses is declared in the header"
                ),
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
                    noodles_bcf::io::Reader::from(raw)
                        .read_record(&mut rec)
                        .map_err(fail)?;
                    w.write_variant_record(&self.header, &rec).map_err(fail)?;
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
        Ok(&self.line)
    }
}

/// Where the GT values of a text record are: the index of `GT` among the
/// FORMAT subfields.
#[derive(Clone, Copy)]
struct GtPos(usize);

/// When `raw` has samples and a GT FORMAT subfield, writes a copy of it to
/// `scratch` with every sample's GT replaced by `./.` and returns the GT
/// position; otherwise `None`.
fn blank_gt(raw: &[u8], scratch: &mut Vec<u8>) -> Option<GtPos> {
    let mut cols = raw.split(|&b| b == b'\t');
    let format = cols.nth(8)?;
    let gt = format.split(|&b| b == b':').position(|k| k == b"GT")?;
    // The samples start after the FORMAT column.
    let samples_at = format.as_ptr() as usize - raw.as_ptr() as usize + format.len();
    if samples_at >= raw.len() {
        return None;
    }
    scratch.clear();
    scratch.extend_from_slice(&raw[..samples_at]);
    for sample in raw[samples_at + 1..].split(|&b| b == b'\t') {
        scratch.push(b'\t');
        for (k, f) in sample.split(|&b| b == b':').enumerate() {
            if k > 0 {
                scratch.push(b':');
            }
            if k == gt {
                scratch.extend_from_slice(b"./.");
            } else {
                scratch.extend_from_slice(f);
            }
        }
        // A sample with fewer subfields than the GT position has no GT.
    }
    Some(GtPos(gt))
}

/// The GT value text of every sample of `raw`.
fn gt_values(raw: &[u8], gt: GtPos) -> Vec<&[u8]> {
    raw.split(|&b| b == b'\t')
        .skip(9)
        .map(|s| s.split(|&b| b == b':').nth(gt.0).unwrap_or_default())
        .collect()
}

/// Replaces the GT field of the framed BCF record `rec` (encoded with GT
/// blanked) by the field encoded from the GT values of the text record
/// `raw`, and returns the new framed record. `key` is the string
/// dictionary index of `GT`.
fn splice_gt(rec: &[u8], key: usize, raw: &[u8], gt: GtPos) -> Result<Vec<u8>, String> {
    const BAD: &str = "internal error: noodles wrote a malformed BCF record";
    let u32_at = |at: usize| -> Result<usize, String> {
        rec.get(at..at + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
            .ok_or_else(|| BAD.to_string())
    };
    let (l_shared, l_indiv) = (u32_at(0)?, u32_at(4)?);
    let indiv_at = 8 + l_shared;
    let n_fmt = u32_at(8 + 20)? >> 24;
    let n_samples = u32_at(8 + 20)? & 0x00ff_ffff;
    let indiv = rec.get(indiv_at..indiv_at + l_indiv).ok_or(BAD)?;
    let mut cur = indiv;
    let mut span = None;
    for _ in 0..n_fmt {
        let start = indiv.len() - cur.len();
        let k = typed::read(&mut cur).map_err(|_| BAD.to_string())?;
        let body = indiv.len() - cur.len();
        let is_gt = matches!(k.int(0), Some(typed::Int::Value(v)) if v as usize == key);
        let (ty, n) = typed::read_descriptor(&mut cur).map_err(|_| BAD.to_string())?;
        let len = n * ty.width() * n_samples;
        if len > cur.len() {
            return Err(BAD.into());
        }
        cur = &cur[len..];
        if is_gt {
            span = Some((start, body, indiv.len() - cur.len()));
            break;
        }
    }
    let Some((_, body, end)) = span else {
        return Err(BAD.into());
    };
    let mut block = Vec::new();
    encode_gt(&gt_values(raw, gt), &mut block)
        .map_err(|e| format!("cannot write the GT field as BCF: {e}"))?;
    let new_indiv_len = l_indiv - (end - body) + block.len();
    let mut out = Vec::with_capacity(rec.len() + block.len());
    out.extend_from_slice(&rec[..4]);
    out.extend_from_slice(&(new_indiv_len as u32).to_le_bytes());
    out.extend_from_slice(&rec[8..indiv_at + body]);
    out.extend_from_slice(&block);
    out.extend_from_slice(&rec[indiv_at + end..]);
    Ok(out)
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
