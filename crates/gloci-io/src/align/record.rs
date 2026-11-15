//! One alignment record as an interval (D101): a borrowed view of a BAM record,
//! and the blocks of its CIGAR (D105).

use super::encode::CIGAR_LETTERS;
use super::header::AlnHeader;
use super::sam::SamAln;
use crate::detect::Format;
use gloci_core::{Interval, Strand};
use std::borrow::Cow;

/// Bytes before the read name in a BAM record, after `block_size`.
const FIXED: usize = 32;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;

/// A view of one BAM record (the bytes after `block_size`), checked by [`Aln::parse`].
#[derive(Debug, Clone, Copy)]
pub struct Aln<'a> {
    raw: &'a [u8],
    ref_id: i32,
    pos: i32,
    flag: u16,
    mapq: u8,
    name: &'a [u8],
    /// The CIGAR words as stored in the record.
    inline_cigar: &'a [u8],
    /// The CIGAR words in use: `inline_cigar`, or the `CG:B:I` array that replaces
    /// a placeholder (htslib's convention for more than 65,535 operations).
    cigar: &'a [u8],
    /// The aux data after the sequence and qualities.
    tags: &'a [u8],
}

fn u32_at(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])
}

fn i32_at(b: &[u8], at: usize) -> i32 {
    u32_at(b, at) as i32
}

fn words(c: &[u8]) -> impl Iterator<Item = (u8, u32)> + '_ {
    c.as_chunks::<4>().0.iter().map(|w| {
        let v = u32::from_le_bytes(*w);
        ((v & 0xf) as u8, v >> 4)
    })
}

/// Check the operation codes of the CIGAR words `c` and return the bases of the read they
/// cover (`M I S = X`), in one pass.
fn scan_cigar(c: &[u8]) -> Result<u64, String> {
    // Bit `op` is set for the operations that consume the read: 0 M, 1 I, 4 S, 7 =, 8 X.
    const QUERY: u32 = 1 << 0 | 1 << 1 | 1 << 4 | 1 << 7 | 1 << 8;
    let mut q = 0u64;
    for w in c.as_chunks::<4>().0 {
        let v = u32::from_le_bytes(*w);
        let op = v & 0xf;
        if op > 9 {
            return Err(format!(
                "the CIGAR has the operation code {op}; BAM defines 0 to 9"
            ));
        }
        if QUERY >> op & 1 != 0 {
            q += u64::from(v >> 4);
        }
    }
    Ok(q)
}

/// One aux field: `(tag, type, value bytes)`. For `B` arrays the value starts with
/// the element type and the count.
type Tag<'a> = ([u8; 2], u8, &'a [u8]);

/// The aux field at the front of `data` and its encoded length, or `None` if malformed.
fn one_tag(data: &[u8]) -> Option<(usize, Tag<'_>)> {
    let head = data.get(..3)?;
    let (tag, ty) = ([head[0], head[1]], head[2]);
    let body = &data[3..];
    let size = match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => memchr::memchr(0, body)? + 1,
        b'B' => {
            let w = match *body.first()? {
                b'c' | b'C' => 1usize,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            let n = u32::from_le_bytes(body.get(1..5)?.try_into().ok()?) as usize;
            5usize.checked_add(n.checked_mul(w)?)?
        }
        _ => return None,
    };
    let value = body.get(..size)?;
    Some((3 + size, (tag, ty, value)))
}

/// The aux fields of `data`, stopping (with `Err`) at one that is malformed.
fn tags(mut data: &[u8]) -> impl Iterator<Item = Result<Tag<'_>, ()>> {
    std::iter::from_fn(move || {
        if data.is_empty() {
            return None;
        }
        match one_tag(data) {
            Some((used, t)) => {
                data = &data[used..];
                Some(Ok(t))
            }
            None => {
                data = &[];
                Some(Err(()))
            }
        }
    })
}

impl<'a> Aln<'a> {
    /// [`Aln::parse`] of a record a BAM stream framed (`raw` is its `block_size` bytes): one
    /// shorter than the fixed fields gets the message every BAM reader gives for it.
    pub(crate) fn parse_framed(raw: &'a [u8], h: &AlnHeader) -> Result<Aln<'a>, String> {
        if raw.len() < FIXED {
            return Err(format!(
                "block_size {} is below the {FIXED} bytes of fixed fields of a BAM record; the \
file is corrupt",
                raw.len()
            ));
        }
        Aln::parse(raw, h)
    }

    /// Check and view the record `raw` (without its `block_size`).
    #[inline]
    pub fn parse(raw: &'a [u8], h: &AlnHeader) -> Result<Aln<'a>, String> {
        match Aln::parse_common(raw, h) {
            Some(a) => Ok(a),
            None => Aln::parse_in_full(raw, h),
        }
    }

    /// [`Aln::parse`] of a record that passes every check and has no `CG` placeholder, as
    /// nearly every record does, without building a message; `None` for any other, which
    /// [`Aln::parse_in_full`] then checks in order.
    #[inline]
    fn parse_common(raw: &'a [u8], h: &AlnHeader) -> Option<Aln<'a>> {
        // Bit `op` is set for the operations that consume the read: 0 M, 1 I, 4 S, 7 =, 8 X.
        const QUERY: u32 = 1 << 0 | 1 << 1 | 1 << 4 | 1 << 7 | 1 << 8;
        let fixed = raw.get(..FIXED)?;
        let ref_id = i32_at(fixed, 0);
        let pos = i32_at(fixed, 4);
        let l_name = usize::from(fixed[8]);
        let n_cigar = usize::from(u16::from_le_bytes([fixed[12], fixed[13]]));
        let flag = u16::from_le_bytes([fixed[14], fixed[15]]);
        let l_seq = u32_at(fixed, 16) as usize;
        let name_end = FIXED + l_name;
        let cigar_end = name_end + 4 * n_cigar;
        let tags_at = cigar_end.checked_add(l_seq.checked_add(l_seq.div_ceil(2))?)?;
        if ref_id < -1
            || ref_id as i64 >= h.refs.len() as i64
            || pos < -1
            || l_name < 2
            || tags_at > raw.len()
            || raw[name_end - 1] != 0
        {
            return None;
        }
        let inline_cigar = &raw[name_end..cigar_end];
        // A placeholder CIGAR (`kS nN`) may have its real one in a CG tag.
        if n_cigar == 2 && inline_cigar[0] & 0xf == 4 && inline_cigar[4] & 0xf == 3 {
            return None;
        }
        let mut q = 0u64;
        let mut bad = false;
        for w in inline_cigar.as_chunks::<4>().0 {
            let v = u32::from_le_bytes(*w);
            let op = v & 0xf;
            bad |= op > 9;
            if QUERY >> op & 1 != 0 {
                q += u64::from(v >> 4);
            }
        }
        let mapped = flag & FLAG_UNMAPPED == 0;
        if bad
            || mapped && (ref_id < 0 || pos < 0 || (l_seq > 0 && n_cigar > 0 && q != l_seq as u64))
        {
            return None;
        }
        Some(Aln {
            raw,
            ref_id,
            pos,
            flag,
            mapq: fixed[9],
            name: &raw[FIXED..name_end - 1],
            inline_cigar,
            cigar: inline_cigar,
            tags: &raw[tags_at..],
        })
    }

    /// [`Aln::parse`] with every check in order, each failing with its message, and a `CG`
    /// placeholder resolved.
    #[cold]
    #[inline(never)]
    fn parse_in_full(raw: &'a [u8], h: &AlnHeader) -> Result<Aln<'a>, String> {
        if raw.len() < FIXED {
            return Err(format!(
                "the record is {} bytes, shorter than the {FIXED} fixed bytes of a BAM record",
                raw.len()
            ));
        }
        let ref_id = i32_at(raw, 0);
        let pos = i32_at(raw, 4);
        if ref_id < -1 || ref_id as i64 >= h.refs.len() as i64 {
            return Err(format!(
                "refID {ref_id} is outside the {} references of the header",
                h.refs.len()
            ));
        }
        if pos < -1 {
            return Err(format!("pos {pos} is below -1"));
        }
        let l_name = usize::from(raw[8]);
        let mapq = raw[9];
        let n_cigar = usize::from(u16::from_le_bytes([raw[12], raw[13]]));
        let flag = u16::from_le_bytes([raw[14], raw[15]]);
        let l_seq = u32_at(raw, 16) as usize;
        let name_end = FIXED + l_name;
        if l_name < 2 {
            return Err("the read name is empty; a BAM read name has 1 to 254 characters".into());
        }
        if name_end > raw.len() || raw[name_end - 1] != 0 {
            return Err("the read name is not terminated by a NUL byte".into());
        }
        let name = &raw[FIXED..name_end - 1];
        let cigar_end = name_end + 4 * n_cigar;
        if cigar_end > raw.len() {
            return Err(format!(
                "the CIGAR of {n_cigar} operations runs past the end of the record"
            ));
        }
        let seq_bytes = l_seq / 2 + (l_seq & 1);
        let tags_at = l_seq
            .checked_add(seq_bytes)
            .and_then(|s| cigar_end.checked_add(s))
            .filter(|&t| t <= raw.len())
            .ok_or_else(|| {
                format!("the sequence of {l_seq} bases runs past the end of the record")
            })?;
        let inline_cigar = &raw[name_end..cigar_end];
        let tags = &raw[tags_at..];
        let mut a = Aln {
            raw,
            ref_id,
            pos,
            flag,
            mapq,
            name,
            inline_cigar,
            cigar: inline_cigar,
            tags,
        };
        let mut q = scan_cigar(inline_cigar)?;
        // A placeholder CIGAR (`kS nN`) plus a CG tag holds the real one.
        let mut it = words(inline_cigar);
        if n_cigar == 2
            && it.next().is_some_and(|(op, _)| op == 4)
            && it.next().is_some_and(|(op, _)| op == 3)
            && ref_id >= 0
            && let Some(cg) = a.cg_array()?
        {
            q = scan_cigar(cg)?;
            a.cigar = cg;
        }
        if flag & FLAG_UNMAPPED == 0 {
            if ref_id < 0 || pos < 0 {
                return Err(format!(
                    "a mapped read (flag {flag}) has refID {ref_id} and pos {pos}; \
the flag should include 4 (unmapped), or the record needs a reference and position"
                ));
            }
            if l_seq > 0 && !a.cigar.is_empty() && q != l_seq as u64 {
                return Err(format!(
                    "the CIGAR covers {q} bases of the read but the sequence has {l_seq}"
                ));
            }
        }
        Ok(a)
    }

    /// View the record `raw` that an [`AlnSource`](super::AlnSource) returned, and so checked:
    /// the fields [`Aln::parse`] gives, without checking the CIGAR again. A record that does
    /// not frame (or holds a `CG` placeholder) is parsed in full.
    pub fn view(raw: &'a [u8], h: &AlnHeader) -> Result<Aln<'a>, String> {
        let Some(fixed) = raw.get(..FIXED) else {
            return Aln::parse(raw, h);
        };
        let l_name = usize::from(fixed[8]);
        let n_cigar = usize::from(u16::from_le_bytes([fixed[12], fixed[13]]));
        let l_seq = u32_at(fixed, 16) as usize;
        let name_end = FIXED + l_name;
        let cigar_end = name_end + 4 * n_cigar;
        let tags_at = l_seq
            .checked_add(l_seq.div_ceil(2))
            .and_then(|s| cigar_end.checked_add(s))
            .filter(|&t| t <= raw.len());
        let Some(tags_at) = tags_at.filter(|_| l_name >= 2) else {
            return Aln::parse(raw, h);
        };
        let inline_cigar = &raw[name_end..cigar_end];
        // The placeholder of a CIGAR held in a `CG` tag (`kS nN`).
        if n_cigar == 2 && inline_cigar[0] & 0xf == 4 && inline_cigar[4] & 0xf == 3 {
            return Aln::parse(raw, h);
        }
        Ok(Aln {
            raw,
            ref_id: i32_at(fixed, 0),
            pos: i32_at(fixed, 4),
            flag: u16::from_le_bytes([fixed[14], fixed[15]]),
            mapq: fixed[9],
            name: &raw[FIXED..name_end - 1],
            inline_cigar,
            cigar: inline_cigar,
            tags: &raw[tags_at..],
        })
    }

    /// The words of a `CG:B:I` tag.
    fn cg_array(&self) -> Result<Option<&'a [u8]>, String> {
        for t in tags(self.tags) {
            let (tag, ty, v) =
                t.map_err(|()| "the aux data is malformed where the CG tag should be".to_string())?;
            if tag == *b"CG" {
                if ty != b'B' || v[0] != b'I' {
                    return Err("the CG tag is not a B:I array".into());
                }
                return Ok(Some(&v[5..]));
            }
        }
        Ok(None)
    }

    pub fn ref_id(&self) -> Option<usize> {
        usize::try_from(self.ref_id).ok()
    }

    /// The 0-based position, or `None` when the record has none.
    pub fn pos(&self) -> Option<u64> {
        u64::try_from(self.pos).ok()
    }

    pub fn flag(&self) -> u16 {
        self.flag
    }

    pub fn mapq(&self) -> u8 {
        self.mapq
    }

    /// The reference of the mate, or `None` when the record has none.
    pub fn mate_ref_id(&self) -> Option<usize> {
        usize::try_from(i32_at(self.raw, 20)).ok()
    }

    /// The 0-based position of the mate, or `None` when the record has none.
    pub fn mate_pos(&self) -> Option<u64> {
        u64::try_from(i32_at(self.raw, 24)).ok()
    }

    /// The number of bases in the read (0 for SEQ `*`).
    pub fn seq_len(&self) -> usize {
        u32_at(self.raw, 16) as usize
    }

    /// The raw Phred qualities: empty for SEQ `*`, and `0xFF` for every base of a record
    /// whose QUAL is missing.
    pub fn quals(&self) -> Cow<'a, [u8]> {
        let l_seq = self.seq_len();
        let tags_at = self.raw.len() - self.tags.len();
        Cow::Borrowed(&self.raw[tags_at - l_seq..tags_at])
    }

    /// Append the bases of the read as letters (nothing for SEQ `*`).
    pub fn seq_into(&self, out: &mut Vec<u8>) {
        let l_seq = self.seq_len();
        let tags_at = self.raw.len() - self.tags.len();
        let seq = &self.raw[tags_at - l_seq - l_seq.div_ceil(2)..tags_at - l_seq];
        out.extend((0..l_seq).map(|i| {
            let code = (seq[i / 2] >> if i % 2 == 0 { 4 } else { 0 }) & 0xf;
            b"=ACMGRSVTWYHKDBN"[usize::from(code)]
        }));
    }

    /// The value of a `Z` string aux field, without its NUL; `None` when the record has no
    /// such field or it has another type.
    pub fn str_tag(&self, tag: [u8; 2]) -> Option<&'a [u8]> {
        let (_, ty, v) = tags(self.tags).map_while(Result::ok).find(|t| t.0 == tag)?;
        (ty == b'Z').then(|| &v[..v.len() - 1])
    }

    pub fn name(&self) -> &'a [u8] {
        self.name
    }

    pub fn raw(&self) -> &'a [u8] {
        self.raw
    }

    /// Call `f(letter, length)` for each CIGAR operation, `CG`-resolved, with SAM letters.
    pub fn for_each_op(&self, f: &mut dyn FnMut(u8, u32)) {
        for (op, n) in words(self.cigar) {
            f(CIGAR_LETTERS[usize::from(op)], n);
        }
    }

    /// An integer aux field of any integer type.
    pub fn int_tag(&self, tag: [u8; 2]) -> Option<i64> {
        let (_, ty, v) = tags(self.tags).map_while(Result::ok).find(|t| t.0 == tag)?;
        Some(match ty {
            b'c' => i64::from(v[0] as i8),
            b'C' => i64::from(v[0]),
            b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
            b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
            b'i' => i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]])),
            b'I' => i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])),
            _ => return None,
        })
    }

    /// The type of an aux field, with the BAM integer types folded into `i`
    /// (`i`, `f`, `A`, `Z`, `H` or `B`); `None` when the record has no such tag.
    pub fn tag_type(&self, tag: [u8; 2]) -> Option<u8> {
        let (_, ty, _) = tags(self.tags).map_while(Result::ok).find(|t| t.0 == tag)?;
        Some(match ty {
            b'c' | b'C' | b's' | b'S' | b'I' => b'i',
            t => t,
        })
    }

    /// The reference length the CIGAR covers (`M D N = X`), at least 1 for a read with none
    /// (htslib's `bam_endpos`).
    fn ref_len(&self) -> u64 {
        let n: u64 = words(self.cigar)
            .filter(|(op, _)| matches!(op, 0 | 2 | 3 | 7 | 8))
            .map(|(_, n)| u64::from(n))
            .sum();
        n.max(1)
    }

    /// The footprint of a mapped read, cut at `ref_len` (the length of its reference),
    /// and whether it was cut. `None` for an unmapped read (D101).
    pub fn span(&self, ref_len: u64) -> Option<(Interval, bool)> {
        footprint(self.flag, self.pos(), self.ref_len(), ref_len)
    }

    pub fn strand(&self) -> Strand {
        if self.flag & FLAG_REVERSE != 0 {
            Strand::Minus
        } else {
            Strand::Plus
        }
    }

    /// Append the SAM text of this record (no newline) to `out`, as `samtools view` writes it:
    /// the resolved CIGAR (a `CG` array that replaced a placeholder is not repeated as a tag),
    /// `=` for a mate on the same reference, `*` for absent fields, and floats in `%g` form.
    pub fn write_sam(&self, h: &AlnHeader, out: &mut Vec<u8>) -> Result<(), String> {
        let raw = self.raw;
        let ref_name = |id: i32| -> Option<&[u8]> {
            usize::try_from(id)
                .ok()
                .and_then(|i| h.refs.get(i))
                .map(|(n, _)| &n[..])
        };
        let num = |out: &mut Vec<u8>, n: i64| out.extend_from_slice(n.to_string().as_bytes());
        out.extend_from_slice(self.name);
        out.push(b'\t');
        num(out, i64::from(self.flag));
        out.push(b'\t');
        out.extend_from_slice(ref_name(self.ref_id).unwrap_or(b"*"));
        out.push(b'\t');
        num(out, i64::from(self.pos) + 1);
        out.push(b'\t');
        num(out, i64::from(self.mapq));
        out.push(b'\t');
        if self.cigar.is_empty() {
            out.push(b'*');
        }
        for (op, n) in words(self.cigar) {
            num(out, i64::from(n));
            out.push(CIGAR_LETTERS[usize::from(op)]);
        }
        out.push(b'\t');
        let mate = i32_at(raw, 20);
        match ref_name(mate) {
            _ if mate < 0 => out.push(b'*'),
            _ if mate == self.ref_id => out.push(b'='),
            Some(n) => out.extend_from_slice(n),
            None => out.push(b'*'),
        }
        out.push(b'\t');
        num(out, i64::from(i32_at(raw, 24)) + 1);
        out.push(b'\t');
        num(out, i64::from(i32_at(raw, 28)));
        out.push(b'\t');
        let l_seq = u32_at(raw, 16) as usize;
        let tags_at = raw.len() - self.tags.len();
        let qual = &raw[tags_at - l_seq..tags_at];
        let seq = &raw[tags_at - l_seq - l_seq.div_ceil(2)..tags_at - l_seq];
        if l_seq == 0 {
            out.extend_from_slice(b"*\t*");
        } else {
            for i in 0..l_seq {
                let code = (seq[i / 2] >> if i % 2 == 0 { 4 } else { 0 }) & 0xf;
                out.push(b"=ACMGRSVTWYHKDBN"[usize::from(code)]);
            }
            out.push(b'\t');
            if qual[0] == 0xff {
                out.push(b'*');
            } else {
                out.extend(qual.iter().map(|q| q.wrapping_add(33)));
            }
        }
        let resolved_cg = self.cigar.as_ptr() != self.inline_cigar.as_ptr();
        for t in tags(self.tags) {
            let (tag, ty, v) = t.map_err(|()| "the aux data is malformed".to_string())?;
            if resolved_cg && tag == *b"CG" {
                continue;
            }
            out.push(b'\t');
            out.extend_from_slice(&tag);
            out.push(b':');
            let int = |v: &[u8]| -> i64 {
                match ty {
                    b'c' => i64::from(v[0] as i8),
                    b'C' => i64::from(v[0]),
                    b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
                    b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
                    b'i' => i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]])),
                    _ => i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])),
                }
            };
            match ty {
                b'A' => {
                    out.extend_from_slice(b"A:");
                    out.push(v[0]);
                }
                b'c' | b'C' | b's' | b'S' | b'i' | b'I' => {
                    out.extend_from_slice(b"i:");
                    num(out, int(v));
                }
                b'f' => {
                    out.extend_from_slice(b"f:");
                    push_float(out, f32::from_le_bytes([v[0], v[1], v[2], v[3]]));
                }
                b'Z' | b'H' => {
                    out.push(ty);
                    out.push(b':');
                    out.extend_from_slice(&v[..v.len() - 1]);
                }
                _ => {
                    // `B`: the element type, the count, then the elements.
                    let elem = v[0];
                    out.extend_from_slice(b"B:");
                    out.push(elem);
                    for e in v[5..].chunks_exact(match elem {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        _ => 4,
                    }) {
                        out.push(b',');
                        match elem {
                            b'f' => push_float(out, f32::from_le_bytes([e[0], e[1], e[2], e[3]])),
                            b'c' => num(out, i64::from(e[0] as i8)),
                            b'C' => num(out, i64::from(e[0])),
                            b's' => num(out, i64::from(i16::from_le_bytes([e[0], e[1]]))),
                            b'S' => num(out, i64::from(u16::from_le_bytes([e[0], e[1]]))),
                            b'i' => {
                                num(out, i64::from(i32::from_le_bytes([e[0], e[1], e[2], e[3]])))
                            }
                            _ => num(out, i64::from(u32::from_le_bytes([e[0], e[1], e[2], e[3]]))),
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Append the slim form of this record (D115) to `out`: the fixed fields (with no
    /// sequence), the name, the CIGAR words as stored, and the `CG` tag and `keep` tag.
    /// It parses with [`Aln::parse`] and gives the same name, flag, MAPQ, operations and span.
    pub fn write_slim(&self, keep: Option<[u8; 2]>, out: &mut Vec<u8>) {
        let at = out.len();
        out.extend_from_slice(&self.raw[..FIXED]);
        out[at + 16..at + 20].fill(0); // l_seq = 0
        out.extend_from_slice(&self.raw[FIXED..FIXED + usize::from(self.raw[8])]);
        out.extend_from_slice(self.inline_cigar);
        let mut rest = self.tags;
        for t in tags(self.tags) {
            let Ok((tag, _, v)) = t else { break };
            let total = 3 + v.len();
            if tag == *b"CG" || Some(tag) == keep {
                out.extend_from_slice(&rest[..total]);
            }
            rest = &rest[total..];
        }
    }
}

/// Append `f` as C's `%g` does (six significant digits, no trailing zeros), which is how
/// htslib writes the float of an aux field.
fn push_float(out: &mut Vec<u8>, f: f32) {
    if f.is_nan() {
        out.extend_from_slice(if f.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
        return;
    }
    if f.is_infinite() {
        out.extend_from_slice(if f < 0.0 { b"-inf" } else { b"inf" });
        return;
    }
    let d = f64::from(f);
    if d == 0.0 {
        out.extend_from_slice(if d.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    // `{:.5e}` is d.dddddE<exp>, rounded as printf rounds.
    let sci = format!("{d:.5e}");
    let (mant, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let strip = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if !(-4..6).contains(&exp) {
        let sign = if exp < 0 { '-' } else { '+' };
        out.extend_from_slice(format!("{}e{sign}{:02}", strip(mant), exp.abs()).as_bytes());
    } else {
        let decimals = usize::try_from(5 - exp).unwrap_or(0);
        out.extend_from_slice(strip(&format!("{d:.decimals$}")).as_bytes());
    }
}

/// The footprint of a read that covers `covered` reference bases from `pos`, cut at `ref_len`
/// (the length of its reference), and whether it was cut. `None` for an unmapped read (D101).
pub(super) fn footprint(
    flag: u16,
    pos: Option<u64>,
    covered: u64,
    ref_len: u64,
) -> Option<(Interval, bool)> {
    if flag & FLAG_UNMAPPED != 0 {
        return None;
    }
    let start = pos?;
    let end = start.saturating_add(covered);
    if end > ref_len {
        let cut = Interval {
            start: start.min(ref_len),
            end: ref_len,
        };
        return Some((cut, true));
    }
    Some((Interval { start, end }, false))
}

/// An alignment record of any format, viewed the same way.
#[derive(Debug, Clone, Copy)]
pub enum AlnRef<'a> {
    Bam(Aln<'a>),
    Sam(SamAln<'a>),
}

macro_rules! delegate {
    ($(#[$m:meta])* $name:ident($($arg:ident: $t:ty),*) -> $ret:ty) => {
        $(#[$m])*
        pub fn $name(&self $(, $arg: $t)*) -> $ret {
            match self {
                AlnRef::Bam(a) => a.$name($($arg),*),
                AlnRef::Sam(a) => a.$name($($arg),*),
            }
        }
    };
}

impl<'a> AlnRef<'a> {
    delegate!(ref_id() -> Option<usize>);
    delegate!(pos() -> Option<u64>);
    delegate!(flag() -> u16);
    delegate!(mapq() -> u8);
    delegate!(mate_ref_id() -> Option<usize>);
    delegate!(mate_pos() -> Option<u64>);
    delegate!(seq_len() -> usize);
    delegate!(quals() -> Cow<'a, [u8]>);
    delegate!(seq_into(out: &mut Vec<u8>) -> ());
    delegate!(str_tag(tag: [u8; 2]) -> Option<&'a [u8]>);
    delegate!(strand() -> Strand);
    delegate!(int_tag(tag: [u8; 2]) -> Option<i64>);
    delegate!(tag_type(tag: [u8; 2]) -> Option<u8>);
    delegate!(span(ref_len: u64) -> Option<(Interval, bool)>);
    delegate!(for_each_op(f: &mut dyn FnMut(u8, u32)) -> ());

    pub fn name(&self) -> &'a [u8] {
        match self {
            AlnRef::Bam(a) => a.name(),
            AlnRef::Sam(a) => a.name(),
        }
    }
}

/// View the record `raw` of a file in `f`'s layout.
pub fn aln_ref<'a>(f: Format, raw: &'a [u8], h: &AlnHeader) -> Result<AlnRef<'a>, String> {
    match f {
        Format::Bam | Format::Cram => Aln::parse(raw, h).map(AlnRef::Bam),
        Format::Sam => SamAln::parse(raw, h).map(AlnRef::Sam),
        other => Err(format!(
            "{} records cannot be read as alignments here",
            other.name()
        )),
    }
}

/// [`aln_ref`] of a record that an [`AlnSource`](super::AlnSource) returned, and so checked:
/// BAM bytes are viewed with [`Aln::view`], without checking them again.
pub fn aln_view<'a>(f: Format, raw: &'a [u8], h: &AlnHeader) -> Result<AlnRef<'a>, String> {
    match f {
        Format::Bam | Format::Cram => Aln::view(raw, h).map(AlnRef::Bam),
        _ => aln_ref(f, raw, h),
    }
}

/// The reference blocks of a mapped read (D105), appended to `out`: the segments between
/// `N` operations, also cut at `D` with `split_deletions`. Empty blocks are dropped; with
/// none left, the span of the read is the block. Blocks are not cut at the reference length.
pub fn blocks(r: &AlnRef<'_>, split_deletions: bool, out: &mut Vec<Interval>) {
    let Some(pos) = r.pos() else { return };
    let first = out.len();
    let mut start = pos;
    let mut cur = pos;
    let flush = |start: u64, end: u64, out: &mut Vec<Interval>| {
        if end > start {
            out.push(Interval { start, end });
        }
    };
    r.for_each_op(&mut |op, n| {
        let n = u64::from(n);
        match op {
            b'M' | b'=' | b'X' => cur = cur.saturating_add(n),
            b'D' if !split_deletions => cur = cur.saturating_add(n),
            b'D' | b'N' => {
                flush(start, cur, out);
                cur = cur.saturating_add(n);
                start = cur;
            }
            _ => {}
        }
    });
    flush(start, cur, out);
    if out.len() == first {
        // The read has no aligned bases: its span is one base (D101).
        let end = r.span(u64::MAX).map_or(pos + 1, |(iv, _)| iv.end);
        out.push(Interval { start: pos, end });
    }
}

/// [`blocks`] with D101's cut at the end of the reference (`ref_len`, its `LN`): blocks that
/// start at or past it are dropped and the end of each other block is cut there, so the blocks
/// agree with the span of the read ([`AlnRef::span`]). With no block left, the span is the one
/// block (an empty `out` when the read has no span at all). Every command that compares or
/// prints blocks uses this one rule.
pub fn blocks_within(r: &AlnRef<'_>, split_deletions: bool, ref_len: u64, out: &mut Vec<Interval>) {
    let first = out.len();
    blocks(r, split_deletions, out);
    let mut keep = first;
    for i in first..out.len() {
        let b = Interval {
            start: out[i].start,
            end: out[i].end.min(ref_len),
        };
        if b.start < b.end {
            out[keep] = b;
            keep += 1;
        }
    }
    out.truncate(keep);
    if out.len() == first
        && let Some((span, _)) = r.span(ref_len)
        && span.start < span.end
    {
        out.push(span);
    }
}
