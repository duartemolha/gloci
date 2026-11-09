//! VCF <-> BCF conversion (`view -O bcf` / `-O vcf`).
//!
//! The output header is the input's own header text with the
//! `##glociCommand` line inserted, never a header serialised by noodles:
//! noodles drops `IDX=` attributes when it writes a header, although BCF
//! records are encoded against the IDX-aware dictionaries. noodles is only
//! used to encode and decode records, against a header parsed from the same
//! text whose dictionary positions are set from gloci's own dictionaries
//! (`BcfDicts`), so the header text and the records always agree.

use crate::bcf::{BcfDicts, write_bcf_header};
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
        Ok(Converter {
            name: name.to_string(),
            from,
            header,
            enc,
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
        let rec = noodles_vcf::Record::try_from(raw).map_err(fail)?;
        if let Err(e) = w.write_variant_record(&self.header, &rec) {
            w.get_mut().clear();
            if mixes_missing_haploid(raw) {
                let mut f = raw.split(|&b| b == b'\t');
                let chrom = String::from_utf8_lossy(f.next().unwrap_or_default());
                let pos = String::from_utf8_lossy(f.next().unwrap_or_default());
                return Err(IoError::Parse {
                    path: self.name.clone(),
                    loc,
                    msg: format!(
                        "{chrom}:{pos}: BCF output of a record that mixes a missing haploid \
genotype (GT \".\") with diploid ones is not supported yet; write VCF instead, for example \
-o out.vcf.gz (see docs/differences.md#bcf-mixed-ploidy-missing)"
                    ),
                });
            }
            return Err(fail(e));
        }
        out.write_all(w.get_ref())
            .map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })?;
        w.get_mut().clear();
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

/// Whether the VCF line `raw` has a sample whose GT is a bare `.` (a
/// missing haploid call) next to a sample with two or more alleles, which
/// the BCF encoder cannot write yet.
fn mixes_missing_haploid(raw: &[u8]) -> bool {
    let mut fields = raw.split(|&b| b == b'\t');
    let Some(format) = fields.nth(8) else {
        return false;
    };
    let Some(gt) = format.split(|&b| b == b':').position(|k| k == b"GT") else {
        return false;
    };
    let (mut haploid_missing, mut multi) = (false, false);
    for sample in fields {
        let Some(v) = sample.split(|&b| b == b':').nth(gt) else {
            continue;
        };
        if v == b"." {
            haploid_missing = true;
        } else if v.iter().any(|&b| b == b'/' || b == b'|') {
            multi = true;
        }
    }
    haploid_missing && multi
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

#[cfg(test)]
mod tests {
    use super::mixes_missing_haploid;

    #[test]
    fn detects_missing_haploid_next_to_diploid() {
        let line = |gts: &str| format!("chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP\t{gts}");
        assert!(mixes_missing_haploid(line("0/1:3\t.:4").as_bytes()));
        assert!(mixes_missing_haploid(line(".:1\t1|1:2").as_bytes()));
        assert!(!mixes_missing_haploid(line("0/1:3\t./.:4").as_bytes()));
        assert!(!mixes_missing_haploid(line("0:3\t.:4").as_bytes()));
        assert!(!mixes_missing_haploid(b"chr1\t10\t.\tA\tG\t.\tPASS\t."));
    }
}
