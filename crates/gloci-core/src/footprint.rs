//! Variant footprints: VCF/BCF record to interval (spec 4.2, decision D28).

use crate::error::FootprintError;
use crate::interval::Interval;

/// Which footprint rule to apply to VCF/BCF records.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum VcfSpan {
    /// The reference bases each ALT actually changes (default).
    #[default]
    Allele,
    /// The bases REF covers (extended by END for symbolic ALTs).
    Ref,
}

/// The parts of a VCF record needed to compute its footprint.
#[derive(Clone, Copy, Debug)]
pub struct Variant<'a> {
    /// 1-based POS.
    pub pos: u64,
    pub reference: &'a [u8],
    /// Raw ALT column (comma-separated).
    pub alts: &'a [u8],
    /// INFO/END, if present.
    pub end: Option<u64>,
    /// Raw INFO/SVLEN value (comma-separated, possibly negative or `.`).
    pub svlen: Option<&'a [u8]>,
}

pub(crate) enum Alt<'a> {
    /// `.`, empty, `*`: no footprint.
    None,
    /// `<*>` / `<NON_REF>`: gVCF block.
    Gvcf,
    /// `<INS...>`
    Ins,
    /// Other symbolic ALTs.
    Sv(&'a [u8]),
    Breakend,
    Seq(&'a [u8]),
}

fn symbolic_name(alt: &[u8]) -> Option<&[u8]> {
    if alt.len() >= 2 && alt[0] == b'<' && alt[alt.len() - 1] == b'>' {
        let inner = &alt[1..alt.len() - 1];
        Some(match memchr::memchr(b':', inner) {
            Some(i) => &inner[..i],
            None => inner,
        })
    } else {
        None
    }
}

pub(crate) fn classify(alt: &[u8]) -> Alt<'_> {
    if alt.is_empty() || alt == b"." || alt == b"*" {
        return Alt::None;
    }
    if let Some(name) = symbolic_name(alt) {
        return if name == b"*" || name.eq_ignore_ascii_case(b"NON_REF") {
            Alt::Gvcf
        } else if name.eq_ignore_ascii_case(b"INS") {
            Alt::Ins
        } else {
            Alt::Sv(alt)
        };
    }
    if alt.iter().any(|&b| b == b'[' || b == b']')
        || (alt.len() > 1 && (alt[0] == b'.' || alt[alt.len() - 1] == b'.'))
    {
        return Alt::Breakend;
    }
    Alt::Seq(alt)
}

pub(crate) fn split_alts(alts: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut rest = alts;
    while let Some(i) = memchr::memchr(b',', rest) {
        out.push(&rest[..i]);
        rest = &rest[i + 1..];
    }
    out.push(rest);
    out
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

fn parse_svlen(raw: &[u8]) -> Result<Vec<Option<u64>>, FootprintError> {
    split_alts(raw)
        .into_iter()
        .map(|t| {
            if t == b"." {
                return Ok(None);
            }
            let bad = || FootprintError::BadSvlen { value: lossy(t) };
            let digits = t.strip_prefix(b"-").or(t.strip_prefix(b"+")).unwrap_or(t);
            if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
                return Err(bad());
            }
            std::str::from_utf8(digits)
                .ok()
                .and_then(|s| s.parse::<u64>().ok())
                .map(Some)
                .ok_or_else(bad)
        })
        .collect()
}

fn iv(start: u64, end: u64) -> Interval {
    Interval { start, end }
}

/// Exclusive end for an SV ALT (index `i` of `n`): END, else POS + |SVLEN|.
fn sv_end(
    v: &Variant,
    alt: &[u8],
    i: usize,
    n: usize,
    svlen: &mut Option<Vec<Option<u64>>>,
) -> Result<u64, FootprintError> {
    if let Some(e) = v.end {
        return Ok(e);
    }
    let missing = || FootprintError::MissingSvEnd { alt: lossy(alt) };
    let raw = v.svlen.ok_or_else(missing)?;
    if svlen.is_none() {
        *svlen = Some(parse_svlen(raw)?);
    }
    let vals = svlen.as_ref().expect("just set");
    let val = match vals.len() {
        1 => vals[0],
        k if k == n => vals[i],
        k => return Err(FootprintError::SvlenCount { svlen: k, alts: n }),
    };
    let len = val.ok_or_else(missing)?;
    v.pos.checked_add(len).ok_or(FootprintError::BadSvlen {
        value: len.to_string(),
    })
}

fn validate(v: &Variant, alts: &[Alt]) -> Result<(), FootprintError> {
    if v.pos == 0 {
        return Err(FootprintError::PosZero);
    }
    if v.reference.is_empty() {
        return Err(FootprintError::EmptyRef);
    }
    // `ref_end` and the SEQ footprint add the REF length to POS.
    if (v.pos - 1).checked_add(v.reference.len() as u64).is_none() {
        return Err(FootprintError::PosOverflow { pos: v.pos });
    }
    let symbolic = alts
        .iter()
        .any(|a| matches!(a, Alt::Gvcf | Alt::Ins | Alt::Sv(_)));
    if symbolic
        && let Some(end) = v.end
        && end < v.pos
    {
        return Err(FootprintError::EndBeforePos { pos: v.pos, end });
    }
    Ok(())
}

fn ref_end(v: &Variant) -> u64 {
    v.pos - 1 + v.reference.len() as u64
}

fn allele_footprint(
    v: &Variant,
    alt: &Alt,
    i: usize,
    n: usize,
    svlen: &mut Option<Vec<Option<u64>>>,
) -> Result<Option<Interval>, FootprintError> {
    Ok(match alt {
        Alt::None | Alt::Gvcf => None,
        Alt::Ins => Some(Interval::point(v.pos)),
        Alt::Sv(a) => Some(iv(v.pos, sv_end(v, a, i, n, svlen)?)),
        Alt::Breakend => Some(iv(v.pos - 1, ref_end(v))),
        Alt::Seq(a) => {
            let r = v.reference;
            let prefix = r
                .iter()
                .zip(a.iter())
                .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
                .count();
            let room = r.len().min(a.len()) - prefix;
            let suffix = r[prefix..]
                .iter()
                .rev()
                .zip(a[prefix..].iter().rev())
                .take(room)
                .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
                .count();
            let base = v.pos - 1;
            Some(iv(base + prefix as u64, base + (r.len() - suffix) as u64))
        }
    })
}

/// Allele footprint per ALT, `None` for `*`, `.`, `<*>` and `<NON_REF>`.
pub fn alt_footprints(v: &Variant) -> Result<Vec<Option<Interval>>, FootprintError> {
    let raw = split_alts(v.alts);
    let alts: Vec<Alt> = raw.iter().map(|a| classify(a)).collect();
    validate(v, &alts)?;
    let mut svlen = None;
    alts.iter()
        .enumerate()
        .map(|(i, a)| allele_footprint(v, a, i, alts.len(), &mut svlen))
        .collect()
}

/// Footprint of the whole record under the chosen span rule.
pub fn footprint(v: &Variant, span: VcfSpan) -> Result<Interval, FootprintError> {
    let raw = split_alts(v.alts);
    let alts: Vec<Alt> = raw.iter().map(|a| classify(a)).collect();
    validate(v, &alts)?;
    let n = alts.len();
    let mut svlen = None;
    let mut sv_max_end: Option<u64> = None;
    let mut hull: Option<Interval> = None;
    for (i, a) in alts.iter().enumerate() {
        // Evaluate every SV ALT in both modes so a missing END/SVLEN is always reported.
        if let Alt::Sv(s) = a {
            let e = sv_end(v, s, i, n, &mut svlen)?;
            sv_max_end = Some(sv_max_end.map_or(e, |m| m.max(e)));
        }
        if span == VcfSpan::Allele
            && let Some(f) = allele_footprint(v, a, i, n, &mut svlen)?
        {
            hull = Some(match hull {
                None => f,
                Some(h) => iv(h.start.min(f.start), h.end.max(f.end)),
            });
        }
    }
    if let Some(h) = hull {
        return Ok(h);
    }
    let symbolic = alts
        .iter()
        .any(|a| matches!(a, Alt::Gvcf | Alt::Ins | Alt::Sv(_)));
    let mut end = ref_end(v);
    if symbolic && let Some(e) = v.end.or(sv_max_end) {
        end = end.max(e);
    }
    Ok(iv(v.pos - 1, end))
}
