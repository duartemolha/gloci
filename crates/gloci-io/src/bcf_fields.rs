//! Native BCF encoding of the INFO and FORMAT fields of a VCF text record,
//! byte for byte what bcftools 1.21 (htslib `vcf_parse_format`) writes.
//!
//! noodles encodes the fixed part of the record (CHROM to FILTER); this module
//! encodes the rest. noodles cannot: it writes a count-0 descriptor for an
//! integer field that is missing everywhere yet still writes data, it fails on
//! a missing float vector and panics on a missing scalar INFO value.
//!
//! The rules, as measured against bcftools 1.21 (`bcf-conformance` golden):
//! * a numeric field is `n` elements wide, `n` being the longest list of any
//!   sample (1 when every sample is `.`); a shorter list, a `.` and a sample
//!   cut short before the field are padded with the end-of-vector value after
//!   the first element, which is the missing value;
//! * integers use the narrowest of int8, int16 and int32 that holds every
//!   value (int8 is -120..=127, the codes above it being reserved);
//! * a String or Character field is as wide as its longest value, plus one
//!   byte (a NUL after every value) unless it is the first FORMAT key;
//! * a FORMAT key that every sample is cut short before is a zero-width field.

use crate::bcf_encode::GtEncoder;
use std::collections::HashMap;

/// How the value of a field is encoded.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Kind {
    Int,
    Float,
    Str,
    Flag,
    /// FORMAT/GT.
    Gt,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Def {
    /// Position in the BCF string dictionary.
    pub idx: usize,
    pub kind: Kind,
}

/// The INFO and FORMAT fields a header declares.
#[derive(Default)]
pub(crate) struct Defs {
    pub info: HashMap<Vec<u8>, Def>,
    pub format: HashMap<Vec<u8>, Def>,
}

const MISSING: i64 = i64::MIN;
const FLOAT_MISSING: u32 = 0x7F80_0001;
const FLOAT_END: u32 = 0x7F80_0002;

#[derive(Default)]
pub(crate) struct FieldEncoder {
    vals: Vec<i64>,
    lens: Vec<usize>,
    /// Per sample: where its cells start in `cells`, and how many there are.
    rows: Vec<(usize, usize)>,
    cells: Vec<(usize, usize)>,
    gt: GtEncoder,
}

fn put_desc(out: &mut Vec<u8>, n: usize, code: u8) {
    if n < 15 {
        out.push((n as u8) << 4 | code);
    } else {
        out.push(15 << 4 | code);
        put_int(out, n as i64);
    }
}

/// A typed integer scalar (its own descriptor, then the narrowest width).
fn put_int(out: &mut Vec<u8>, v: i64) {
    if let Ok(v) = i8::try_from(v) {
        out.extend_from_slice(&[1 << 4 | 1, v as u8]);
    } else if let Ok(v) = i16::try_from(v) {
        out.push(1 << 4 | 2);
        out.extend_from_slice(&v.to_le_bytes());
    } else {
        out.push(1 << 4 | 3);
        out.extend_from_slice(&(v as i32).to_le_bytes());
    }
}

fn parse_int(tok: &[u8], what: &str) -> Result<i64, String> {
    if tok == b"." {
        return Ok(MISSING);
    }
    let bad = || {
        format!(
            "{what} value {:?} is not an integer; write a whole number or . for missing",
            String::from_utf8_lossy(tok)
        )
    };
    let digits = tok.strip_prefix(b"-").or_else(|| tok.strip_prefix(b"+"));
    let body = digits.unwrap_or(tok);
    if body.is_empty() {
        // htslib reads an empty value as an out-of-range integer and sets it
        // missing, with a warning.
        log::warn!("{what}: an empty integer value is read as missing");
        return Ok(MISSING);
    }
    if !body.iter().all(u8::is_ascii_digit) {
        return Err(bad());
    }
    let mut n: i64 = 0;
    for &b in body {
        n = n.saturating_mul(10).saturating_add(i64::from(b - b'0'));
    }
    if tok[0] == b'-' {
        n = -n;
    }
    if !(-2_147_483_640..=i64::from(i32::MAX)).contains(&n) {
        log::warn!("{what}: the value {n} is out of range for BCF and is read as missing");
        return Ok(MISSING);
    }
    Ok(n)
}

/// `at_end`: the token ends the record's line. htslib reads an empty float
/// token as 0 (strtod finds no number), except at the very end of the line,
/// where it is missing.
fn parse_float(tok: &[u8], what: &str, at_end: bool) -> Result<i64, String> {
    if tok == b"." {
        return Ok(MISSING);
    }
    if tok.is_empty() {
        return Ok(if at_end {
            MISSING
        } else {
            i64::from(0f32.to_bits())
        });
    }
    let s = std::str::from_utf8(tok).ok();
    // htslib parses a double and narrows it to a float.
    match s.and_then(|s| s.parse::<f64>().ok()) {
        Some(v) => Ok(i64::from((v as f32).to_bits())),
        None => Err(format!(
            "{what} value {:?} is not a number; write a decimal number or . for missing",
            String::from_utf8_lossy(tok)
        )),
    }
}

impl FieldEncoder {
    /// Appends the encoded INFO entries of the INFO column `col` to `out` and
    /// returns how many there are.
    pub(crate) fn info(
        &mut self,
        defs: &Defs,
        col: &[u8],
        out: &mut Vec<u8>,
    ) -> Result<usize, String> {
        if col == b"." || col.is_empty() {
            return Ok(0);
        }
        let mut n_info = 0;
        for item in col.split(|&b| b == b';') {
            if item.is_empty() {
                continue;
            }
            let (key, val) = match item.iter().position(|&b| b == b'=') {
                Some(i) => (&item[..i], Some(&item[i + 1..])),
                None => (item, None),
            };
            let name = String::from_utf8_lossy(key);
            let def = defs.info.get(key).ok_or_else(|| {
                format!(
                    "INFO field {name} is not declared in the header, and BCF needs every INFO field in its header; \
add ##INFO=<ID={name},Number=1,Type=String,Description=\"...\"> to the header"
                )
            })?;
            put_int(out, def.idx as i64);
            let what = format!("INFO/{name}");
            self.vals.clear();
            self.lens.clear();
            match def.kind {
                Kind::Flag => out.push(0),
                Kind::Int | Kind::Float | Kind::Gt => {
                    let val = val.unwrap_or(b".");
                    for tok in val.split(|&b| b == b',') {
                        let v = if def.kind == Kind::Float {
                            parse_float(tok, &what, false)?
                        } else {
                            parse_int(tok, &what)?
                        };
                        self.vals.push(v);
                    }
                    self.lens.push(self.vals.len());
                    let n = self.vals.len();
                    self.write_numeric(def.kind == Kind::Float, n, out);
                }
                Kind::Str => {
                    let val = val.unwrap_or(b".");
                    put_desc(out, val.len(), 7);
                    out.extend_from_slice(val);
                }
            }
            n_info += 1;
        }
        Ok(n_info)
    }

    /// Appends the FORMAT fields of the record `raw` to `out` (its columns
    /// from FORMAT on) and returns how many there are.
    pub(crate) fn samples(
        &mut self,
        defs: &Defs,
        raw: &[u8],
        n_samples: usize,
        out: &mut Vec<u8>,
    ) -> Result<usize, String> {
        let mut cols = raw.split(|&b| b == b'\t').skip(8);
        let Some(format) = cols.next() else {
            if n_samples == 0 {
                return Ok(0);
            }
            return Err(format!(
                "the record has no FORMAT and sample columns but the header declares {n_samples} sample(s) \
(bcftools refuses it too); add a FORMAT column and a column per sample"
            ));
        };
        if n_samples == 0 {
            return Ok(0);
        }
        // Where every sample's cells are, as offsets into `raw`.
        self.rows.clear();
        self.cells.clear();
        let base = raw.as_ptr() as usize;
        let mut n_cols = 0;
        for s in cols {
            let start = self.cells.len();
            let mut at = s.as_ptr() as usize - base;
            for cell in s.split(|&b| b == b':') {
                self.cells.push((at, cell.len()));
                at += cell.len() + 1;
            }
            self.rows.push((start, self.cells.len() - start));
            n_cols += 1;
        }
        if n_cols != n_samples {
            return Err(format!(
                "the record has {n_cols} sample column(s) but the header declares {n_samples}; \
fix the line so every sample has a column"
            ));
        }
        let mut n_fmt = 0;
        for (k, key) in format.split(|&b| b == b':').enumerate() {
            let name = String::from_utf8_lossy(key);
            let def = defs.format.get(key).ok_or_else(|| {
                format!(
                    "FORMAT field {name} is not declared in the header, and BCF needs every FORMAT field in its header; \
add ##FORMAT=<ID={name},Number=1,Type=String,Description=\"...\"> to the header"
                )
            })?;
            put_int(out, def.idx as i64);
            let what = format!("FORMAT/{name}");
            let cell = |s: usize| -> Option<&[u8]> {
                let (start, n) = self.rows[s];
                (k < n).then(|| {
                    let (at, len) = self.cells[start + k];
                    &raw[at..at + len]
                })
            };
            let cell_end = |c: &[u8]| c.as_ptr() as usize - raw.as_ptr() as usize + c.len();
            let absent_everywhere = (0..n_samples).all(|s| cell(s).is_none());
            match def.kind {
                Kind::Gt => {
                    if absent_everywhere {
                        out.push(0);
                    } else {
                        let gt = &mut self.gt;
                        gt.encode((0..n_samples).map(|s| cell(s).unwrap_or(b"")), out)
                            .map_err(|e| format!("cannot write the GT field as BCF: {e}"))?;
                    }
                }
                Kind::Int | Kind::Float => {
                    let float = def.kind == Kind::Float;
                    {
                        self.vals.clear();
                        self.lens.clear();
                        for s in 0..n_samples {
                            let before = self.vals.len();
                            match cell(s) {
                                None => self.vals.push(MISSING),
                                Some(c) => {
                                    for tok in c.split(|&b| b == b',') {
                                        let v = if float {
                                            parse_float(tok, &what, cell_end(tok) == raw.len())?
                                        } else {
                                            parse_int(tok, &what)?
                                        };
                                        self.vals.push(v);
                                    }
                                }
                            }
                            self.lens.push(self.vals.len() - before);
                        }
                        let n = self.lens.iter().copied().max().unwrap_or(1);
                        self.write_numeric(float, n, out);
                    }
                }
                Kind::Str | Kind::Flag => {
                    if absent_everywhere {
                        put_desc(out, 0, 7);
                    } else {
                        // An empty value is a zero-length string, but at the
                        // very end of the line htslib reads it as missing.
                        let value = |s: usize| match cell(s) {
                            None => &b"."[..],
                            Some(c) if c.is_empty() && cell_end(c) == raw.len() => b".",
                            Some(c) => c,
                        };
                        let longest = (0..n_samples).map(|s| value(s).len()).max().unwrap_or(1);
                        let width = longest + usize::from(k > 0);
                        put_desc(out, width, 7);
                        for s in 0..n_samples {
                            let v = value(s);
                            out.extend_from_slice(v);
                            out.resize(out.len() + (width - v.len()), 0);
                        }
                    }
                }
            }
            n_fmt += 1;
        }
        Ok(n_fmt)
    }

    /// Writes the typed block of `self.vals` (a list per entry of `self.lens`)
    /// `n` elements wide.
    fn write_numeric(&self, float: bool, n: usize, out: &mut Vec<u8>) {
        if float {
            put_desc(out, n, 5);
            let mut at = 0;
            for &len in &self.lens {
                for j in 0..n {
                    let bits = if j >= len {
                        FLOAT_END
                    } else {
                        match self.vals[at + j] {
                            MISSING => FLOAT_MISSING,
                            v => v as u32,
                        }
                    };
                    out.extend_from_slice(&bits.to_le_bytes());
                }
                at += len;
            }
            return;
        }
        let (mut lo, mut hi) = (0i64, 0i64);
        for &v in self.vals.iter().filter(|&&v| v != MISSING) {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        let code = if lo >= -120 && hi <= 127 {
            1u8
        } else if lo >= -32760 && hi <= 32767 {
            2
        } else {
            3
        };
        put_desc(out, n, code);
        let mut at = 0;
        for &len in &self.lens {
            for j in 0..n {
                let v = if j >= len {
                    None
                } else {
                    Some(self.vals[at + j])
                };
                match (code, v) {
                    (1, None) => out.push(0x81),
                    (2, None) => out.extend_from_slice(&0x8001u16.to_le_bytes()),
                    (_, None) => out.extend_from_slice(&0x8000_0001u32.to_le_bytes()),
                    (1, Some(MISSING)) => out.push(0x80),
                    (2, Some(MISSING)) => out.extend_from_slice(&0x8000u16.to_le_bytes()),
                    (_, Some(MISSING)) => out.extend_from_slice(&0x8000_0000u32.to_le_bytes()),
                    (1, Some(v)) => out.push(v as i8 as u8),
                    (2, Some(v)) => out.extend_from_slice(&(v as i16).to_le_bytes()),
                    (_, Some(v)) => out.extend_from_slice(&(v as i32).to_le_bytes()),
                }
            }
            at += len;
        }
    }
}
