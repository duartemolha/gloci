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
//! * a String or Character field is as wide as the longest text a sample has
//!   for the key, plus one byte (a NUL after every value) unless it is the
//!   first FORMAT key; a String or GT key that every sample is cut short
//!   before is a zero-width field (an integer or float one is a missing value);
//! * INFO/END (one value, at least POS) is the record length's end; INFO
//!   numbers are read by their leading digits.

use crate::bcf_encode::GtEncoder;
use crate::error::Loc;
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

/// Values htslib sets to missing or ignores with a warning. They are counted
/// and reported once, with the first place, when the encoder is dropped.
#[derive(Default)]
struct Notes {
    counts: [u64; 4],
    first: [Option<String>; 4],
    name: String,
    at: Option<Loc>,
}

const NOTE_TEXT: [&str; 4] = [
    "an empty integer value is read as missing",
    "an integer value out of the range BCF holds is read as missing",
    "an INFO number that is not a whole or decimal number is read as missing, or by its leading digits",
    "INFO/END smaller than POS is ignored for the record's length",
];

impl Notes {
    fn note(&mut self, kind: usize, what: &str) {
        self.counts[kind] += 1;
        if self.first[kind].is_none() {
            let place = match self.at {
                Some(l) => format!("{} {l}", self.name),
                None => self.name.clone(),
            };
            self.first[kind] = Some(format!("{place}: {what}"));
        }
    }
}

impl Drop for Notes {
    fn drop(&mut self) {
        for (k, text) in NOTE_TEXT.iter().enumerate() {
            if let Some(first) = &self.first[k] {
                log::warn!(
                    "{first}; {} value(s) in all were affected ({text})",
                    self.counts[k]
                );
            }
        }
    }
}

#[derive(Default)]
pub(crate) struct FieldEncoder {
    notes: Notes,
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

/// A FORMAT or INFO integer. `info`: htslib reads an INFO number by its
/// leading digits (`12abc` is 12, `abc` missing) where a FORMAT one must be
/// a whole number.
fn parse_int(tok: &[u8], what: &str, info: bool, notes: &mut Notes) -> Result<i64, String> {
    if tok == b"." {
        return Ok(MISSING);
    }
    let sign = tok.first().is_some_and(|b| matches!(b, b'-' | b'+'));
    let body = if sign { &tok[1..] } else { tok };
    if body.is_empty() {
        // htslib reads an empty value as an out-of-range integer and sets it
        // missing, with a warning.
        notes.note(0, what);
        return Ok(MISSING);
    }
    let digits = body.iter().take_while(|b| b.is_ascii_digit()).count();
    if digits < body.len() {
        if !info {
            return Err(format!(
                "{what} value {:?} is not an integer; write a whole number or . for missing",
                String::from_utf8_lossy(tok)
            ));
        }
        notes.note(
            2,
            &format!("{what} value {:?}", String::from_utf8_lossy(tok)),
        );
        if digits == 0 {
            return Ok(MISSING);
        }
    }
    let mut n: i64 = 0;
    for &b in &body[..digits] {
        n = n.saturating_mul(10).saturating_add(i64::from(b - b'0'));
    }
    if tok[0] == b'-' {
        n = -n;
    }
    if !(-2_147_483_640..=i64::from(i32::MAX)).contains(&n) {
        notes.note(1, &format!("{what} value {n}"));
        return Ok(MISSING);
    }
    Ok(n)
}

/// A FORMAT or INFO float. `at_end`: the token ends the record's line. An
/// empty FORMAT token is 0 for htslib (strtod finds no number), except at
/// the very end of the line, where it is missing; an empty INFO token is
/// missing, and an INFO number is read by its leading characters.
fn parse_float(
    tok: &[u8],
    what: &str,
    info: bool,
    at_end: bool,
    notes: &mut Notes,
) -> Result<i64, String> {
    if tok == b"." {
        return Ok(MISSING);
    }
    if tok.is_empty() {
        return Ok(if info || at_end {
            MISSING
        } else {
            i64::from(0f32.to_bits())
        });
    }
    // htslib parses a double and narrows it to a float.
    let parse = |t: &[u8]| {
        std::str::from_utf8(t)
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
    };
    if let Some(v) = parse(tok) {
        return Ok(i64::from((v as f32).to_bits()));
    }
    if info {
        notes.note(
            2,
            &format!("{what} value {:?}", String::from_utf8_lossy(tok)),
        );
        let ok = |c: &u8| matches!(c, b'+' | b'-' | b'.' | b'0'..=b'9' | b'e' | b'E');
        let mut n = tok.iter().take_while(|c| ok(c)).count();
        while n > 0 {
            if let Some(v) = parse(&tok[..n]) {
                return Ok(i64::from((v as f32).to_bits()));
            }
            n -= 1;
        }
        return Ok(MISSING);
    }
    Err(format!(
        "{what} value {:?} is not a number; write a decimal number or . for missing",
        String::from_utf8_lossy(tok)
    ))
}

impl FieldEncoder {
    /// An encoder whose warnings name `name` (the input).
    pub(crate) fn new(name: &str) -> FieldEncoder {
        let mut e = FieldEncoder::default();
        e.notes.name = name.to_string();
        e
    }

    /// The record the next calls encode, for warnings.
    pub(crate) fn at(&mut self, loc: Loc) {
        self.notes.at = Some(loc);
    }

    /// Notes that `END` is below `POS`, which htslib ignores for the length.
    pub(crate) fn end_before_pos(&mut self, end: i64) {
        self.notes.note(3, &format!("INFO/END={end}"));
    }

    /// Appends the encoded INFO entries of the INFO column `col` to `out` and
    /// returns how many there are, and the value of `END` when it has one.
    pub(crate) fn info(
        &mut self,
        defs: &Defs,
        col: &[u8],
        out: &mut Vec<u8>,
    ) -> Result<(usize, Option<i64>), String> {
        if col == b"." || col.is_empty() {
            return Ok((0, None));
        }
        let mut end = None;
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
            match (def.kind, val) {
                // A flag, or a key with no `=` of any type: no value.
                (Kind::Flag, _) | (_, None) => out.push(0),
                (Kind::Int | Kind::Float | Kind::Gt, Some(val)) => {
                    for tok in val.split(|&b| b == b',') {
                        let v = if def.kind == Kind::Float {
                            parse_float(tok, &what, true, false, &mut self.notes)?
                        } else {
                            parse_int(tok, &what, true, &mut self.notes)?
                        };
                        self.vals.push(v);
                    }
                    self.lens.push(self.vals.len());
                    let n = self.vals.len();
                    if key == b"END" && def.kind == Kind::Int && n == 1 && self.vals[0] != MISSING {
                        end = Some(self.vals[0]);
                    }
                    self.write_numeric(def.kind == Kind::Float, n, out);
                }
                (Kind::Str, Some(val)) => {
                    put_desc(out, val.len(), 7);
                    out.extend_from_slice(val);
                }
            }
            n_info += 1;
        }
        Ok((n_info, end))
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
        // No FORMAT (`.`) means no sample data, whatever the columns hold.
        if n_samples == 0 || format == b"." {
            return Ok(0);
        }
        // Where every sample's cells are, as offsets into `raw`.
        self.rows.clear();
        self.cells.clear();
        let base = raw.as_ptr() as usize;
        let mut n_cols = 0;
        // A column beyond the header's samples is ignored, as htslib does.
        for s in cols.take(n_samples) {
            let start = self.cells.len();
            let mut at = s.as_ptr() as usize - base;
            for cell in s.split(|&b| b == b':') {
                self.cells.push((at, cell.len()));
                at += cell.len() + 1;
            }
            self.rows.push((start, self.cells.len() - start));
            n_cols += 1;
        }
        if n_cols < n_samples {
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
                                            parse_float(
                                                tok,
                                                &what,
                                                false,
                                                cell_end(tok) == raw.len(),
                                                &mut self.notes,
                                            )?
                                        } else {
                                            parse_int(tok, &what, false, &mut self.notes)?
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
                        // The width is that of the longest text a sample has for
                        // the key (a `.` standing in for a missing or cut-short
                        // value does not count), plus one unless first.
                        let longest = (0..n_samples)
                            .filter_map(|s| cell(s).map(<[u8]>::len))
                            .max()
                            .unwrap_or(0);
                        let width = longest + usize::from(k > 0);
                        put_desc(out, width, 7);
                        for s in 0..n_samples {
                            let v = value(s);
                            let v = &v[..v.len().min(width)];
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
