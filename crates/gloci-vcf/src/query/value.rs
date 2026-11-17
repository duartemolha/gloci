//! Typed values, `%TYPE`, `%END` and `%TGT` (D172).

use gloci_io::kputd;
use std::collections::HashSet;
use std::sync::{Arc, Mutex};

use super::TagType;
use crate::srsort::{VCF_BND, VCF_INDEL, VCF_MNP, VCF_OTHER, VCF_OVERLAP, VCF_SNP, alt_type};

/// The Integer range VCF and BCF hold: htslib keeps the 8 values above `i32::MIN` for
/// missing and vector end and others it reserves; bcftools reads a value outside as missing.
const INT_MIN: i64 = -2_147_483_640;
const INT_MAX: i64 = 2_147_483_647;

/// The tags already warned about in one run, shared by the clones of a query (the workers),
/// so an out-of-range Integer warns once per tag (D172). Queries compare equal whatever
/// they have warned about.
#[derive(Clone, Debug, Default)]
pub struct Warned(Arc<Mutex<HashSet<String>>>);

impl PartialEq for Warned {
    fn eq(&self, _: &Warned) -> bool {
        true
    }
}

impl Warned {
    /// Warn that `what` (`INFO/I`) holds `value`, the first time for that tag.
    pub(super) fn out_of_range(&self, what: &str, value: &str) {
        let first = self
            .0
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .insert(what.to_string());
        if first {
            log::warn!(
                "vcf query: {what} value {value} is outside the 32-bit range of VCF integers; \
                 printed as missing (further values of this tag are not reported)"
            );
        }
    }
}

/// A comma-separated value as the typed tag `ty` prints it, or the subscript `sub`
/// of it (`.` past the end). Err is the element that does not read as the type; Ok is the
/// first Integer element outside the range VCF holds (printed as missing), if any.
pub(super) fn write_typed(
    value: &[u8],
    ty: TagType,
    sub: Option<usize>,
    out: &mut Vec<u8>,
) -> Result<Option<String>, String> {
    let mut range = None;
    match sub {
        Some(n) => match value.split(|&b| b == b',').nth(n) {
            Some(e) => write_element(e, ty, out, &mut range)?,
            None => out.push(b'.'),
        },
        None => {
            for (i, e) in value.split(|&b| b == b',').enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                write_element(e, ty, out, &mut range)?;
            }
        }
    }
    Ok(range)
}

fn write_element(
    e: &[u8],
    ty: TagType,
    out: &mut Vec<u8>,
    range: &mut Option<String>,
) -> Result<(), String> {
    if e.is_empty() || e == b"." {
        out.push(b'.');
        return Ok(());
    }
    let bad = || String::from_utf8_lossy(e).into_owned();
    match ty {
        TagType::Integer => {
            let text = std::str::from_utf8(e).map_err(|_| bad())?;
            match text.parse::<i64>() {
                Ok(n) if (INT_MIN..=INT_MAX).contains(&n) => {
                    out.extend_from_slice(itoa::Buffer::new().format(n).as_bytes());
                }
                // Beyond i64 is out of range too, when it is digits (bcftools: missing).
                Ok(_) => {
                    range.get_or_insert_with(|| text.to_string());
                    out.push(b'.');
                }
                Err(_) if is_integer_text(e) => {
                    range.get_or_insert_with(|| text.to_string());
                    out.push(b'.');
                }
                Err(_) => return Err(bad()),
            }
        }
        TagType::Float => write_float(e, out).ok_or_else(bad)?,
        TagType::String | TagType::Flag => out.extend_from_slice(e),
    }
    Ok(())
}

/// An optional sign and at least one digit, nothing else.
fn is_integer_text(e: &[u8]) -> bool {
    let digits = e.strip_prefix(b"-").or(e.strip_prefix(b"+")).unwrap_or(e);
    !digits.is_empty() && digits.iter().all(u8::is_ascii_digit)
}

/// A float as `TextToBcf` stores it and bcftools prints it: read as 32 bits, written by `kputd`.
pub(super) fn write_float(e: &[u8], out: &mut Vec<u8>) -> Option<()> {
    let x: f64 = std::str::from_utf8(e).ok()?.parse().ok()?;
    kputd(f64::from(x as f32), out);
    Some(())
}

/// htslib's variant classes of the ALTs, in the order bcftools lists them, or `REF`.
pub(super) fn write_types(reference: &[u8], alt: &[u8], out: &mut Vec<u8>) {
    let mut bits = 0;
    if alt != b"." {
        for a in alt.split(|&b| b == b',') {
            if a != b"." {
                bits |= alt_type(reference, a);
            }
        }
    }
    let names = [
        (VCF_SNP, "SNP"),
        (VCF_MNP, "MNP"),
        (VCF_INDEL, "INDEL"),
        (VCF_OTHER, "OTHER"),
        (VCF_BND, "BND"),
        (VCF_OVERLAP, "OVERLAP"),
    ];
    let mut any = false;
    for (bit, name) in names {
        if bits & bit != 0 {
            if any {
                out.push(b',');
            }
            out.extend_from_slice(name.as_bytes());
            any = true;
        }
    }
    if !any {
        out.extend_from_slice(b"REF");
    }
}

/// INFO/END when it is at least POS, else POS + len(REF) - 1.
pub(super) fn end_of(pos: i64, reference: &[u8], info: &[u8]) -> i64 {
    if info != b"." {
        for item in info.split(|&b| b == b';') {
            if let Some(v) = item.strip_prefix(b"END=")
                && let Some(end) = std::str::from_utf8(v)
                    .ok()
                    .and_then(|t| t.parse::<i64>().ok())
                && end >= pos
            {
                return end;
            }
        }
    }
    pos + reference.len() as i64 - 1
}

/// A genotype with each allele number replaced by its bases.
pub(super) fn write_tgt(
    gt: &[u8],
    reference: &[u8],
    alt: &[u8],
    out: &mut Vec<u8>,
) -> Result<(), String> {
    if gt.is_empty() {
        out.push(b'.');
        return Ok(());
    }
    let mut start = 0;
    for i in 0..=gt.len() {
        if i < gt.len() && gt[i] != b'/' && gt[i] != b'|' {
            continue;
        }
        let allele = &gt[start..i];
        if allele == b"." {
            out.push(b'.');
        } else {
            let k: usize = std::str::from_utf8(allele)
                .ok()
                .and_then(|t| t.parse().ok())
                .ok_or_else(|| {
                    format!("genotype \"{}\" is not valid", String::from_utf8_lossy(gt))
                })?;
            let bases = if k == 0 {
                Some(reference)
            } else if alt == b"." {
                None
            } else {
                alt.split(|&b| b == b',').nth(k - 1)
            };
            match bases {
                Some(b) => out.extend_from_slice(b),
                None => return Err(format!("allele {k} is not in the record's alleles")),
            }
        }
        if i < gt.len() {
            out.push(gt[i]);
        }
        start = i + 1;
    }
    Ok(())
}
