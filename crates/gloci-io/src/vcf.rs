//! VCF text reader, header parsing and provenance line.

use crate::bcf::BcfDicts;
use crate::diag::InvalidPolicy;
use crate::error::IoError;
use crate::input::Input;
use crate::record::{Class, Header, Parsed, TextKind, TextReader, coord, intern_chrom};
use gloci_core::{ChromTable, Strand, Variant, VcfSpan, footprint};
use std::ffi::OsString;

/// A parsed VCF header; `raw` holds every header line (including `#CHROM`),
/// each terminated by `\n`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfHeader {
    pub raw: Vec<u8>,
    /// `##contig` IDs with their optional lengths, in header order.
    pub contigs: Vec<(Vec<u8>, Option<u64>)>,
    /// `#CHROM` columns from the tenth onward.
    pub samples: Vec<Vec<u8>>,
    /// BCF dictionaries; `None` for text VCF.
    pub bcf: Option<BcfDicts>,
}

impl VcfHeader {
    /// Parse `##contig` and `#CHROM` lines by hand (D22).
    pub fn parse(raw: Vec<u8>) -> Result<VcfHeader, String> {
        let mut contigs = Vec::new();
        let mut samples = None;
        for line in raw.split(|&b| b == b'\n').filter(|l| !l.is_empty()) {
            // BCF header text does not go through the text reader's CR stripping.
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if line == b"#CHROM" || line.starts_with(b"#CHROM\t") {
                if samples.is_some() {
                    return Err("VCF header has more than one #CHROM line".into());
                }
                samples = Some(
                    line.split(|&b| b == b'\t')
                        .skip(9)
                        .map(<[u8]>::to_vec)
                        .collect(),
                );
            } else if let Some(body) = line.strip_prefix(b"##contig=<") {
                contigs.push(parse_contig(line, body)?);
            }
        }
        let samples = samples.ok_or("VCF header has no #CHROM line")?;
        Ok(VcfHeader {
            raw,
            contigs,
            samples,
            bcf: None,
        })
    }

    /// The structured definition lines (`##INFO`, `##FILTER`, `##FORMAT`,
    /// `##ALT`, `##contig`) by kind and ID, each as its whole line without
    /// the line terminator. A kind and ID defined twice keeps the first.
    pub fn definitions(&self) -> std::collections::BTreeMap<(&'static str, Vec<u8>), Vec<u8>> {
        const KINDS: [(&str, &[u8]); 5] = [
            ("INFO", b"##INFO=<"),
            ("FILTER", b"##FILTER=<"),
            ("FORMAT", b"##FORMAT=<"),
            ("ALT", b"##ALT=<"),
            ("contig", b"##contig=<"),
        ];
        let mut out = std::collections::BTreeMap::new();
        for line in self.raw.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            for (kind, prefix) in KINDS {
                if let Some(body) = line.strip_prefix(prefix) {
                    if let Some((_, id)) =
                        header_fields(body).into_iter().find(|(k, _)| *k == b"ID")
                    {
                        out.entry((kind, id.to_vec()))
                            .or_insert_with(|| line.to_vec());
                    }
                    break;
                }
            }
        }
        out
    }

    /// `raw` with a `##glociCommand` line inserted immediately before `#CHROM`.
    pub fn with_command_line(&self, cmd: &str, version: &str) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.raw.len() + cmd.len() + 40);
        let mut done = false;
        for line in self.raw.split_inclusive(|&b| b == b'\n') {
            let clean = line.strip_suffix(b"\n").unwrap_or(line);
            let clean = clean.strip_suffix(b"\r").unwrap_or(clean);
            if !done && (clean == b"#CHROM" || clean.starts_with(b"#CHROM\t")) {
                out.extend_from_slice(
                    format!("##glociCommand={cmd}; Version={version}\n").as_bytes(),
                );
                done = true;
            }
            out.extend_from_slice(line);
        }
        out
    }
}

fn parse_contig(line: &[u8], body: &[u8]) -> Result<(Vec<u8>, Option<u64>), String> {
    let shown = String::from_utf8_lossy(line);
    let (mut id, mut len) = (None, None);
    for (k, v) in header_fields(body) {
        if k == b"ID" {
            id.get_or_insert(v.to_vec());
        } else if k == b"length" {
            let n = std::str::from_utf8(v)
                .ok()
                .and_then(|s| s.parse::<u64>().ok())
                .ok_or_else(|| format!("header line {shown} has a non-integer length"))?;
            len = Some(n);
        }
    }
    let id = id
        .filter(|i| !i.is_empty())
        .ok_or_else(|| format!("header line {shown} has no ID"))?;
    Ok((id, len))
}

/// Split the body of a structured header line (after `##KEY=<`) into
/// `key=value` pairs. Commas inside double quotes do not split, quotes are
/// removed from the value, and a trailing `>` is ignored.
pub(crate) fn header_fields(body: &[u8]) -> Vec<(&[u8], &[u8])> {
    let body = body.strip_suffix(b">").unwrap_or(body);
    let mut out = Vec::new();
    let mut i = 0;
    while i < body.len() {
        let start = i;
        let (mut quoted, mut escaped) = (false, false);
        while i < body.len() {
            let b = body[i];
            if escaped {
                escaped = false;
            } else if quoted && b == b'\\' {
                escaped = true;
            } else if b == b'"' {
                quoted = !quoted;
            } else if b == b',' && !quoted {
                break;
            }
            i += 1;
        }
        let field = &body[start..i];
        i += 1; // the comma
        let (k, v) = match memchr::memchr(b'=', field) {
            Some(eq) => (&field[..eq], &field[eq + 1..]),
            None => (field, &field[field.len()..]),
        };
        let v = match v {
            [b'"', inner @ .., b'"'] => inner,
            v => v,
        };
        out.push((k, v));
    }
    out
}

/// Quote a command line for the `##glociCommand` header line: arguments with
/// whitespace or shell metacharacters are single-quoted, and tabs and
/// newlines are escaped so the result stays on one line.
pub fn shell_quote_args(args: &[OsString]) -> String {
    let safe = |c: char| c.is_ascii_alphanumeric() || "_-./=:,@%+".contains(c);
    let quoted: Vec<String> = args
        .iter()
        .map(|a| {
            let a = a.to_string_lossy();
            if !a.is_empty() && a.chars().all(safe) {
                return a.into_owned();
            }
            let mut s = String::from("'");
            for c in a.chars() {
                match c {
                    '\'' => s.push_str("'\\''"),
                    '\t' => s.push_str("\\t"),
                    '\n' => s.push_str("\\n"),
                    '\r' => s.push_str("\\r"),
                    c => s.push(c),
                }
            }
            s.push('\'');
            s
        })
        .collect();
    quoted.join(" ")
}

pub(crate) struct Vcf {
    span: VcfSpan,
}
pub(crate) type VcfReader = TextReader<Vcf>;

impl Vcf {
    pub(crate) fn new(span: VcfSpan) -> Vcf {
        Vcf { span }
    }
}

impl TextReader<Vcf> {
    pub(crate) fn open_vcf(
        input: Input,
        policy: InvalidPolicy,
        span: VcfSpan,
    ) -> Result<Self, IoError> {
        let name = input.name.clone();
        let mut r = Self::open(input, policy, Vcf { span })?;
        let raw = match std::mem::replace(&mut r.header, Header::None) {
            Header::Lines(l) => l,
            _ => Vec::new(),
        };
        let h = VcfHeader::parse(raw).map_err(|msg| IoError::Detect {
            msg: format!("{name}: {msg}"),
            path: name,
        })?;
        r.header = Header::Vcf(h);
        Ok(r)
    }
}

impl TextKind for Vcf {
    fn classify(line: &[u8], _seen: bool) -> Class {
        if line.is_empty() {
            Class::Skip
        } else if line[0] == b'#' {
            Class::Header
        } else {
            Class::Record
        }
    }

    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let mut f = line.splitn(9, |&b| b == b'\t');
        let cols: Vec<&[u8]> = f.by_ref().take(8).collect();
        if cols.len() < 8 {
            return Err(format!(
                "expected at least 8 tab-separated columns, found {}",
                cols.len()
            ));
        }
        let pos = coord("POS", cols[1])?;
        if pos == 0 {
            return Err(gloci_core::FootprintError::PosZero.to_string());
        }
        let (reference, alts, info) = (cols[3], cols[4], cols[7]);
        let (mut end, mut svlen) = (None, None);
        if memchr::memchr(b'<', alts).is_some() {
            for kv in info.split(|&b| b == b';') {
                if end.is_none()
                    && let Some(v) = kv.strip_prefix(b"END=")
                {
                    end = Some(coord("INFO/END", v)?);
                }
                if svlen.is_none()
                    && let Some(v) = kv.strip_prefix(b"SVLEN=")
                {
                    svlen = Some(v);
                }
            }
        }
        let span = footprint(
            &Variant {
                pos,
                reference,
                alts,
                end,
                svlen,
            },
            self.span,
        )
        .map_err(|e| e.to_string())?;
        let chrom = intern_chrom(chroms, cols[0])?;
        Ok((chrom, pos - 1, span, Strand::Unknown))
    }
}
