//! The VCF line of a BCF record, as htslib writes it.
//!
//! bcftools prints a BCF record through htslib's `vcf_format`, and gloci means to write the
//! same bytes. This module is that function for a framed record (`l_shared`, `l_indiv`, the
//! shared part, the individual part), with nothing in between: every value is formatted from
//! its BCF bytes, so a string is written as it is stored (a `,`, a `;` or a `%` in it is not
//! percent-encoded), a number is in htslib's notation (`kputd`: six significant digits,
//! `8.08575e-05`), and the quirks are htslib's too (a GT of the integer-missing value, which
//! a sample cut short before its GT holds, prints as `-65`).
//!
//! Ported from htslib 1.21: `vcf_format`, `bcf_fmt_array` and `bcf_fmt_array1` (vcf.c),
//! `bcf_format_gt` (htslib/vcf.h) and `kputd` (kstring.c), which are distributed under the
//! following licence.
//!
//! ```text
//! The MIT/Expat License
//!
//! Copyright (C) 2012-2024 Genome Research Ltd.
//!
//! Permission is hereby granted, free of charge, to any person obtaining a copy
//! of this software and associated documentation files (the "Software"), to deal
//! in the Software without restriction, including without limitation the rights
//! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! copies of the Software, and to permit persons to whom the Software is
//! furnished to do so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included in
//! all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//! DEALINGS IN THE SOFTWARE.
//! ```

use crate::bcf::BcfDicts;
use crate::bcf::typed::{self, Ty};
use gloci_core::cfmt::{write_g6, write_i64};

const MALFORMED: &str = "malformed BCF record";

/// htslib's `kputd`: a number as C's `%g` writes it (six significant digits, no trailing
/// zeros), computed as htslib computes it so that the last digit agrees.
pub(crate) fn kputd(d: f64, out: &mut Vec<u8>) {
    if d == 0.0 {
        out.extend_from_slice(if d.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    let mut d = d;
    if d < 0.0 {
        out.push(b'-');
        d = -d;
    }
    if !(0.0001..=999999.0).contains(&d) {
        // htslib leaves the exponent cases (and nan and inf) to printf's %g.
        write_g6(out, d);
        return;
    }
    // Six digits (seven when the scaled value rounds up to 1000000), written at the end of
    // `buf`; `cp` is where the text starts and `first` where the digits were to start.
    let mut buf = [0u8; 21];
    let (scale, back) = if d < 0.001 {
        (1_000_000_000.0, 1)
    } else if d < 0.01 {
        (100_000_000.0, 2)
    } else if d < 0.1 {
        (10_000_000.0, 3)
    } else if d < 1.0 {
        (1_000_000.0, 4)
    } else if d < 10.0 {
        (100_000.0, 5)
    } else if d < 100.0 {
        (10_000.0, 6)
    } else if d < 1000.0 {
        (1000.0, 7)
    } else if d < 10_000.0 {
        (100.0, 8)
    } else if d < 100_000.0 {
        (10.0, 9)
    } else {
        (1.0, 10)
    };
    let mut i = (d * scale).round_ties_even() as u32;
    let mut cp = 20 - back;
    for k in 0..3 {
        cp -= 2;
        let two = i % 100;
        buf[cp] = b'0' + (two / 10) as u8;
        buf[cp + 1] = b'0' + (two % 10) as u8;
        if k < 2 {
            i /= 100;
        }
    }
    if i >= 100 {
        cp -= 1;
        buf[cp] = b'0' + (i / 100) as u8;
    }
    let mut p = 20 - cp;
    let mut ep;
    if p <= 10 {
        // Below 1: 0.00123 is 123, so zeros and "0." go in front.
        ep = cp + 5;
        while p < 10 {
            cp -= 1;
            buf[cp] = b'0';
            p += 1;
        }
        cp -= 1;
        buf[cp] = b'.';
        cp -= 1;
        buf[cp] = b'0';
    } else {
        // 123.001 is 123001: the integer part moves down one place and "." goes in.
        cp -= 1;
        let mut xp = cp;
        ep = cp + 6;
        while p > 10 {
            buf[xp] = buf[xp + 1];
            xp += 1;
            p -= 1;
        }
        buf[xp] = b'.';
    }
    while buf[ep] == b'0' && ep > cp {
        ep -= 1;
    }
    if buf[ep] != b'.' {
        ep += 1;
    }
    out.extend_from_slice(&buf[cp..ep]);
}

/// The integer at `at` of a vector of `ty`, sign-extended, with its type's missing value.
#[inline]
fn int_at(ty: Ty, data: &[u8], at: usize) -> (i32, i32) {
    match ty {
        Ty::Int8 => (i32::from(data[at] as i8), i32::from(i8::MIN)),
        Ty::Int16 => (
            i32::from(i16::from_le_bytes([data[2 * at], data[2 * at + 1]])),
            i32::from(i16::MIN),
        ),
        _ => (
            i32::from_le_bytes([
                data[4 * at],
                data[4 * at + 1],
                data[4 * at + 2],
                data[4 * at + 3],
            ]),
            i32::MIN,
        ),
    }
}

const FLOAT_MISSING: u32 = 0x7f80_0001;
const FLOAT_END: u32 = 0x7f80_0002;
/// htslib's `bcf_str_missing`.
const STR_MISSING: u8 = 0x07;

#[inline]
fn float_bits(data: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([
        data[4 * at],
        data[4 * at + 1],
        data[4 * at + 2],
        data[4 * at + 3],
    ])
}

/// htslib's `bcf_fmt_array`: `n` values of `ty` in `data` (which holds at least that many).
fn fmt_array(out: &mut Vec<u8>, ty: Ty, n: usize, data: &[u8]) -> Result<(), String> {
    if n == 0 {
        out.push(b'.');
        return Ok(());
    }
    match ty {
        Ty::Char => {
            let s = &data[..n];
            let end = memchr::memchr(0, s).unwrap_or(n);
            out.extend_from_slice(&s[..end]);
        }
        Ty::Int8 | Ty::Int16 | Ty::Int32 => {
            for j in 0..n {
                let (v, missing) = int_at(ty, data, j);
                if v == missing + 1 {
                    break;
                }
                if j > 0 {
                    out.push(b',');
                }
                if v == missing {
                    out.push(b'.');
                } else {
                    write_i64(out, i64::from(v));
                }
            }
        }
        Ty::Float => {
            for j in 0..n {
                let bits = float_bits(data, j);
                if bits == FLOAT_END {
                    break;
                }
                if j > 0 {
                    out.push(b',');
                }
                if bits == FLOAT_MISSING {
                    out.push(b'.');
                } else {
                    kputd(f64::from(f32::from_bits(bits)), out);
                }
            }
        }
        Ty::Missing => return Err("a value has no type but a length".into()),
    }
    Ok(())
}

/// htslib's `bcf_fmt_array1`: the one value of a FORMAT field of length 1. Unlike
/// `bcf_fmt_array` it writes a character whatever it is, the missing string as `.`.
#[inline]
fn fmt_array1(out: &mut Vec<u8>, ty: Ty, data: &[u8]) -> Result<(), String> {
    match ty {
        Ty::Char => out.push(if data[0] == STR_MISSING {
            b'.'
        } else {
            data[0]
        }),
        Ty::Int8 | Ty::Int16 | Ty::Int32 => {
            let (v, missing) = int_at(ty, data, 0);
            if v != missing + 1 {
                if v == missing {
                    out.push(b'.');
                } else {
                    write_i64(out, i64::from(v));
                }
            }
        }
        Ty::Float => {
            let bits = float_bits(data, 0);
            if bits != FLOAT_END {
                if bits == FLOAT_MISSING {
                    out.push(b'.');
                } else {
                    kputd(f64::from(f32::from_bits(bits)), out);
                }
            }
        }
        Ty::Missing => return Err("a FORMAT value has no type but a length".into()),
    }
    Ok(())
}

/// htslib's `bcf_format_gt`: the genotype of one sample. An allele is `(value >> 1) - 1`, 0
/// the missing allele, and the low bit says the allele is phased with the one before.
#[inline]
fn format_gt(out: &mut Vec<u8>, ty: Ty, n: usize, data: &[u8]) -> Result<(), String> {
    match ty {
        Ty::Missing => out.push(b'.'),
        Ty::Int8 if n == 2 => {
            // The common case, diploid calls of small alleles, without the general loop.
            let (a, b) = (data[0] as i8, data[1] as i8);
            if a == i8::MIN + 1 {
                out.push(b'.');
                return Ok(());
            }
            put_allele(out, i32::from(a));
            if b != i8::MIN + 1 {
                out.push(if b & 1 != 0 { b'|' } else { b'/' });
                put_allele(out, i32::from(b));
            }
        }
        Ty::Int8 | Ty::Int16 | Ty::Int32 => {
            let mut i = 0;
            while i < n {
                let (v, missing) = int_at(ty, data, i);
                if v == missing + 1 {
                    break;
                }
                if i > 0 {
                    out.push(if v & 1 != 0 { b'|' } else { b'/' });
                }
                put_allele(out, v);
                i += 1;
            }
            if i == 0 {
                out.push(b'.');
            }
        }
        Ty::Float | Ty::Char => return Err("the GT field is not stored as integers".into()),
    }
    Ok(())
}

#[inline]
fn put_allele(out: &mut Vec<u8>, v: i32) {
    let a = v >> 1;
    if a == 0 {
        out.push(b'.');
    } else if (1..=10).contains(&a) {
        out.push(b'0' + (a - 1) as u8);
    } else {
        write_i64(out, i64::from(a) - 1);
    }
}

/// One FORMAT field of a record: its type, the values per sample and the bytes of all
/// samples.
struct Field<'a> {
    ty: Ty,
    n: usize,
    size: usize,
    data: &'a [u8],
    gt: bool,
}

fn u32_at(b: &[u8], at: usize) -> Result<u32, String> {
    b.get(at..at + 4)
        .map(|w| u32::from_le_bytes([w[0], w[1], w[2], w[3]]))
        .ok_or_else(|| MALFORMED.to_string())
}

/// The dictionary id a typed integer holds (the key of an INFO or FORMAT field).
fn key_of(cur: &mut &[u8]) -> Result<i64, String> {
    let v = typed::read(cur)?;
    match v.int(0) {
        Some(typed::Int::Value(k)) => Ok(i64::from(k)),
        _ => Err(MALFORMED.into()),
    }
}

fn name_of<'d>(strings: &'d [Vec<u8>], id: i64, what: &str) -> Result<&'d [u8], String> {
    usize::try_from(id)
        .ok()
        .and_then(|i| strings.get(i))
        .filter(|s| !s.is_empty())
        .map(Vec::as_slice)
        .ok_or_else(|| format!("the {what} id {id} is not in the header"))
}

/// Append the VCF line of the framed BCF record `raw` to `out`, without a line terminator,
/// as htslib's `vcf_format` writes it against the dictionaries `dicts`.
pub(crate) fn write_vcf_line(
    raw: &[u8],
    dicts: &BcfDicts,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let (l_shared, l_indiv) = (u32_at(raw, 0)? as usize, u32_at(raw, 4)? as usize);
    let shared = raw.get(8..8 + l_shared).ok_or(MALFORMED)?;
    let indiv = raw
        .get(8 + l_shared..)
        .and_then(|r| r.get(..l_indiv))
        .ok_or(MALFORMED)?;
    if shared.len() < 24 {
        return Err(MALFORMED.into());
    }
    let rid = u32_at(shared, 0)? as i32;
    let pos = u32_at(shared, 4)? as i32;
    let qual = u32_at(shared, 12)?;
    let n_allele_info = u32_at(shared, 16)?;
    let (n_allele, n_info) = (n_allele_info >> 16, n_allele_info & 0xffff);
    let n_fmt_sample = u32_at(shared, 20)?;
    let (n_fmt, n_sample) = (
        (n_fmt_sample >> 24) as usize,
        (n_fmt_sample & 0x00ff_ffff) as usize,
    );
    let mut cur = &shared[24..];

    // CHROM, POS
    let chrom = usize::try_from(rid)
        .ok()
        .and_then(|i| dicts.contigs.get(i))
        .filter(|c| !c.is_empty())
        .ok_or_else(|| format!("the contig id {rid} is not in the header"))?;
    out.extend_from_slice(chrom);
    out.push(b'\t');
    write_i64(out, i64::from(pos) + 1);

    // ID
    out.push(b'\t');
    let id = typed::read(&mut cur)?;
    fmt_array(out, id.ty, id.n, id.data)?;

    // REF, ALT
    out.push(b'\t');
    if n_allele == 0 {
        out.extend_from_slice(b".\t.");
    } else {
        for i in 0..n_allele {
            let a = typed::read(&mut cur)?;
            match i {
                0 => {}
                1 => out.push(b'\t'),
                _ => out.push(b','),
            }
            fmt_array(out, a.ty, a.n, a.data)?;
        }
        if n_allele == 1 {
            out.extend_from_slice(b"\t.");
        }
    }

    // QUAL
    out.push(b'\t');
    if qual == FLOAT_MISSING {
        out.push(b'.');
    } else {
        kputd(f64::from(f32::from_bits(qual)), out);
    }

    // FILTER
    out.push(b'\t');
    let flt = typed::read(&mut cur)?;
    if flt.n == 0 {
        out.push(b'.');
    } else {
        if !matches!(flt.ty, Ty::Int8 | Ty::Int16 | Ty::Int32) {
            return Err(MALFORMED.into());
        }
        for i in 0..flt.n {
            let (v, _) = int_at(flt.ty, flt.data, i);
            if i > 0 {
                out.push(b';');
            }
            out.extend_from_slice(name_of(&dicts.strings, i64::from(v), "FILTER")?);
        }
    }

    // INFO
    out.push(b'\t');
    if n_info == 0 {
        out.push(b'.');
    }
    for i in 0..n_info {
        let key = key_of(&mut cur)?;
        let v = typed::read(&mut cur)?;
        if i > 0 {
            out.push(b';');
        }
        out.extend_from_slice(name_of(&dicts.strings, key, "INFO")?);
        if v.n == 0 {
            continue;
        }
        out.push(b'=');
        fmt_array(out, v.ty, v.n, v.data)?;
    }

    // FORMAT and the samples
    if n_sample == 0 {
        return Ok(());
    }
    if n_fmt == 0 {
        for _ in 0..=n_sample {
            out.extend_from_slice(b"\t.");
        }
        return Ok(());
    }
    let mut fields: Vec<Field<'_>> = Vec::with_capacity(n_fmt);
    let mut cur = indiv;
    for i in 0..n_fmt {
        let key = key_of(&mut cur)?;
        let (ty, n) = typed::read_descriptor(&mut cur)?;
        let size = n.checked_mul(ty.width()).ok_or(MALFORMED)?;
        let len = size
            .checked_mul(n_sample)
            .filter(|&l| l <= cur.len())
            .ok_or(MALFORMED)?;
        let (data, rest) = cur.split_at(len);
        cur = rest;
        let name = name_of(&dicts.strings, key, "FORMAT")?;
        out.push(if i == 0 { b'\t' } else { b':' });
        out.extend_from_slice(name);
        fields.push(Field {
            ty,
            n,
            size,
            data,
            gt: name == b"GT",
        });
    }
    // htslib takes the last key named GT as the genotype.
    if let Some(last) = fields.iter().rposition(|f| f.gt) {
        for (i, f) in fields.iter_mut().enumerate() {
            f.gt = i == last;
        }
    }
    for j in 0..n_sample {
        for (i, f) in fields.iter().enumerate() {
            out.push(if i == 0 { b'\t' } else { b':' });
            let cell = &f.data[j * f.size..(j + 1) * f.size];
            if f.gt {
                format_gt(out, f.ty, f.n, cell)?;
            } else if f.n == 1 {
                fmt_array1(out, f.ty, cell)?;
            } else {
                fmt_array(out, f.ty, f.n, cell)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(x: f64) -> String {
        let mut o = Vec::new();
        kputd(x, &mut o);
        String::from_utf8(o).unwrap()
    }

    fn f(x: f32) -> String {
        d(f64::from(x))
    }

    /// What bcftools 1.21 prints for these values of a VCF once they have been through BCF
    /// (32-bit floats).
    #[test]
    fn kputd_is_bcftools_notation() {
        for (x, want) in [
            ("0.0000808575", "8.08575e-05"),
            ("0.000012345", "1.2345e-05"),
            ("1e-10", "1e-10"),
            ("123456789.0", "1.23457e+08"),
            ("0.333333333", "0.333333"),
            ("3.14159265", "3.14159"),
            ("1e6", "1e+06"),
            ("999999.4", "999999"),
            ("999999.5", "1e+06"),
            ("0.0001", "0.0001"),
            ("0.00009999995", "9.99999e-05"),
            ("0.009999999", "0.01"),
            ("0.001", "0.001"),
            ("0.01", "0.01"),
            ("0.1", "0.1"),
            ("1.0", "1"),
            ("10.0", "10"),
            ("100.5", "100.5"),
            ("12345.678", "12345.7"),
            ("123456.7", "123457"),
            ("-0.5", "-0.5"),
            ("-1e-8", "-1e-08"),
            ("-123.456", "-123.456"),
            ("1e10", "1e+10"),
            ("0.30000001", "0.3"),
            ("16777217.0", "1.67772e+07"),
            ("0.00012345678", "0.000123457"),
            ("0.0009999999", "0.001"),
            ("99999.95", "100000"),
            ("1.999999", "2"),
            ("0.123456789", "0.123457"),
            ("0.5", "0.5"),
            ("2.0", "2"),
            ("50.0", "50"),
            ("3.0", "3"),
            ("-20.25", "-20.25"),
            ("-1e-7", "-1e-07"),
        ] {
            assert_eq!(f(x.parse().unwrap()), want, "{x}");
        }
        assert_eq!(d(0.0), "0");
        assert_eq!(d(-0.0), "-0");
        assert_eq!(d(f64::NAN), "nan");
        assert_eq!(d(f64::INFINITY), "inf");
        assert_eq!(d(f64::NEG_INFINITY), "-inf");
    }

    /// A record of the smallest valid shape, with `body` after the fixed part and `indiv`.
    fn framed(
        n_allele: u16,
        n_info: u16,
        n_fmt: u8,
        n_sample: u32,
        body: &[u8],
        indiv: &[u8],
    ) -> Vec<u8> {
        let mut shared = Vec::new();
        shared.extend_from_slice(&0i32.to_le_bytes()); // CHROM
        shared.extend_from_slice(&9i32.to_le_bytes()); // POS - 1
        shared.extend_from_slice(&1i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&FLOAT_MISSING.to_le_bytes()); // QUAL
        shared.extend_from_slice(&((u32::from(n_allele) << 16) | u32::from(n_info)).to_le_bytes());
        shared.extend_from_slice(&((u32::from(n_fmt) << 24) | n_sample).to_le_bytes());
        shared.extend_from_slice(body);
        let mut r = Vec::new();
        r.extend_from_slice(&(shared.len() as u32).to_le_bytes());
        r.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
        r.extend_from_slice(&shared);
        r.extend_from_slice(indiv);
        r
    }

    fn dicts() -> BcfDicts {
        BcfDicts {
            contigs: vec![b"chr1".to_vec()],
            strings: vec![b"PASS".to_vec(), b"DP".to_vec(), Vec::new(), b"GT".to_vec()],
        }
    }

    fn line(raw: &[u8]) -> Result<String, String> {
        let mut o = Vec::new();
        write_vcf_line(raw, &dicts(), &mut o)?;
        Ok(String::from_utf8(o).unwrap())
    }

    // ID ".", REF "A", no FILTER.
    const SITE: &[u8] = &[0x07, 0x17, b'A', 0x00];

    #[test]
    fn a_site_without_samples() {
        assert_eq!(
            line(&framed(1, 0, 0, 0, SITE, &[])).unwrap(),
            "chr1\t10\t.\tA\t.\t.\t.\t."
        );
        // No allele at all: REF and ALT are both `.`.
        assert_eq!(
            line(&framed(0, 0, 0, 0, &[0x07, 0x00], &[])).unwrap(),
            "chr1\t10\t.\t.\t.\t.\t.\t."
        );
    }

    /// Samples without any FORMAT field: htslib writes `.` for FORMAT and for each sample.
    #[test]
    fn samples_without_format_fields() {
        assert_eq!(
            line(&framed(1, 0, 0, 2, SITE, &[])).unwrap(),
            "chr1\t10\t.\tA\t.\t.\t.\t.\t.\t.\t."
        );
    }

    #[test]
    fn info_flag_value_and_missing() {
        // DP as a flag (no value), then DP=5, then DP=. (int8 missing).
        let mut body = SITE.to_vec();
        body.extend_from_slice(&[0x11, 1, 0x00, 0x11, 1, 0x11, 5, 0x11, 1, 0x11, 0x80]);
        assert_eq!(
            line(&framed(1, 3, 0, 0, &body, &[])).unwrap(),
            "chr1\t10\t.\tA\t.\t.\t.\tDP;DP=5;DP=."
        );
    }

    /// A GT of the integer-missing value prints as `-65`, as in htslib; a vector that ends at
    /// once is `.`; 16-bit alleles are written like any other.
    #[test]
    fn genotypes_as_htslib_writes_them() {
        let gt = |ty: u8, cells: &[u8]| {
            let mut indiv = vec![0x11, 3, ty];
            indiv.extend_from_slice(cells);
            let n = if ty & 0x0f == 2 {
                cells.len() / 4
            } else {
                cells.len() / 2
            };
            line(&framed(1, 0, 1, n as u32, SITE, &indiv)).unwrap()
        };
        // 0/1, 1|1, ./., a haploid 0 (the vector ends), nothing at all, the missing value.
        assert_eq!(
            gt(0x21, &[2, 4, 4, 5, 0, 0, 2, 0x81, 0x81, 0x81, 0x80, 0x81]),
            "chr1\t10\t.\tA\t.\t.\t.\t.\tGT\t0/1\t1/1\t./.\t0\t.\t-65".replace("1/1", "1|1")
        );
        // 1/65 and 70|0 need 16 bits.
        let cells: Vec<u8> = [4i16, 132, 142, 3]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        assert_eq!(
            gt(0x22, &cells),
            "chr1\t10\t.\tA\t.\t.\t.\t.\tGT\t1/65\t70|0"
        );
    }

    #[test]
    fn a_key_the_header_lacks_is_an_error() {
        let mut body = SITE.to_vec();
        body.extend_from_slice(&[0x11, 2, 0x00]);
        let e = line(&framed(1, 1, 0, 0, &body, &[])).unwrap_err();
        assert_eq!(e, "the INFO id 2 is not in the header");
        let e = line(&framed(1, 0, 1, 1, SITE, &[0x11, 9, 0x11, 1])).unwrap_err();
        assert_eq!(e, "the FORMAT id 9 is not in the header");
        // A record cut short is an error, not a panic.
        let whole = framed(1, 0, 1, 2, SITE, &[0x11, 3, 0x21, 2, 4, 2, 4]);
        for cut in 0..whole.len() {
            let _ = line(&whole[..cut]);
        }
        assert!(line(&framed(1, 0, 1, 3, SITE, &[0x11, 3, 0x21, 2, 4])).is_err());
    }
}
