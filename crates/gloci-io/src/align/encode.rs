//! SAM text to a BAM record, with the choices samtools makes (D104): the smallest integer
//! type that holds each `i` aux value, `bin` from the footprint, and a `CG` array when the
//! CIGAR has more than 65,535 operations. The output is byte for byte what `samtools view -b`
//! writes for the same line.

use memchr::memchr;

const FLAG_UNMAPPED: u16 = 0x4;
/// CIGAR letters by BAM operation code (9 is `B`, which consumes neither query nor reference).
pub(crate) const CIGAR_LETTERS: &[u8; 10] = b"MIDNSHP=XB";

/// htslib's `seq_nt16_table`: the 4-bit code of a base letter. `=` is 0, the IUPAC letters
/// (either case) are 1 to 14, the digits 0 to 3 are A, C, G and T, and everything else is N.
const SEQ_CODE: [u8; 256] = {
    let mut t = [15u8; 256];
    let letters = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;
    while i < 16 {
        t[letters[i] as usize] = i as u8;
        if letters[i].is_ascii_uppercase() {
            t[letters[i].to_ascii_lowercase() as usize] = i as u8;
        }
        i += 1;
    }
    t[b'0' as usize] = 1;
    t[b'1' as usize] = 2;
    t[b'2' as usize] = 4;
    t[b'3' as usize] = 8;
    t
};

/// The letter BAM stores for the base letter `b` (`seq_nt16_str` of `seq_nt16_table`): the
/// IUPAC letter in upper case, `N` for anything else.
pub(crate) fn bam_base(b: u8) -> u8 {
    b"=ACMGRSVTWYHKDBN"[usize::from(SEQ_CODE[usize::from(b)])]
}

/// What htslib does when it reads a SAM record whose place or CIGAR make no sense as mapped.
pub(crate) const ZERO_POS: u8 = 1;
pub(crate) const NO_CIGAR: u8 = 2;

/// htslib's reading rules for a record's flag and reference: a reference with no position is no
/// reference (`ZERO_POS` when it was named), a read with no reference is unmapped, and so is a
/// mapped read with no CIGAR (`NO_CIGAR`). Returns the flag, the reference id and the notes.
pub(crate) fn normalise(flag: u16, ref_id: i32, pos: i64, has_cigar: bool) -> (u16, i32, u8) {
    let (mut flag, mut ref_id, mut notes) = (flag, ref_id, 0u8);
    if pos < 0 {
        if ref_id >= 0 {
            notes |= ZERO_POS;
        }
        ref_id = -1;
    }
    if ref_id < 0 {
        flag |= FLAG_UNMAPPED;
    }
    if flag & FLAG_UNMAPPED == 0 && !has_cigar {
        flag |= FLAG_UNMAPPED;
        notes |= NO_CIGAR;
    }
    (flag, ref_id, notes)
}

/// A FLAG as `strtol` base 0 reads it: decimal, `0x` hexadecimal or leading-`0` octal.
pub(crate) fn parse_flag(field: &[u8]) -> Result<u16, String> {
    let bad = || {
        format!(
            "FLAG '{}' is not a number from 0 to 65535 (decimal, 0x hexadecimal or 0 octal)",
            text(field)
        )
    };
    let t = std::str::from_utf8(field).map_err(|_| bad())?;
    // `strtol` takes a sign, samtools does not take a FLAG with one.
    if !t.bytes().all(|b| b.is_ascii_alphanumeric()) {
        return Err(bad());
    }
    let n = if let Some(h) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        u32::from_str_radix(h, 16)
    } else if t.len() > 1 && t.starts_with('0') {
        u32::from_str_radix(&t[1..], 8)
    } else {
        t.parse::<u32>()
    };
    n.ok().and_then(|n| u16::try_from(n).ok()).ok_or_else(bad)
}

/// htslib's `hts_reg2bin(beg, end, 14, 5)`.
pub(crate) fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let bin = if beg >> 14 == end >> 14 {
        4681 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        585 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        73 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        9 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        1 + (beg >> 26)
    } else {
        0
    };
    bin as u16
}

fn text(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

fn num<T: std::str::FromStr>(field: &[u8], what: &str) -> Result<T, String> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<T>().ok())
        .ok_or_else(|| format!("{what} '{}' is not a valid number for SAM", text(field)))
}

/// A float as htslib reads it: text to `double`, then to `float`.
fn float(field: &[u8], what: &str) -> Result<f32, String> {
    num::<f64>(field, what).map(|d| d as f32)
}

/// The type letter and bytes of an integer aux value, in the smallest type that holds it.
fn push_int(out: &mut Vec<u8>, v: i64, field: &[u8]) -> Result<(), String> {
    if v < 0 {
        if v >= i64::from(i8::MIN) {
            out.extend_from_slice(&[b'c', v as i8 as u8]);
        } else if v >= i64::from(i16::MIN) {
            out.push(b's');
            out.extend_from_slice(&(v as i16).to_le_bytes());
        } else if v >= i64::from(i32::MIN) {
            out.push(b'i');
            out.extend_from_slice(&(v as i32).to_le_bytes());
        } else {
            return Err(format!(
                "the integer {} is below -2147483648, which BAM cannot hold",
                text(field)
            ));
        }
    } else if v <= i64::from(u8::MAX) {
        out.extend_from_slice(&[b'C', v as u8]);
    } else if v <= i64::from(u16::MAX) {
        out.push(b'S');
        out.extend_from_slice(&(v as u16).to_le_bytes());
    } else if v <= i64::from(u32::MAX) {
        out.push(b'I');
        out.extend_from_slice(&(v as u32).to_le_bytes());
    } else {
        return Err(format!(
            "the integer {} is above 4294967295, which BAM cannot hold",
            text(field)
        ));
    }
    Ok(())
}

fn push_aux(out: &mut Vec<u8>, field: &[u8]) -> Result<(), String> {
    if field.len() < 5 || field[2] != b':' || field[4] != b':' {
        return Err(format!(
            "the optional field '{}' is not TAG:TYPE:VALUE",
            text(field)
        ));
    }
    let (tag, ty, value) = (&field[..2], field[3], &field[5..]);
    out.extend_from_slice(tag);
    match ty {
        // htslib reads `c` as a character, like `A`, and `I` as an integer, like `i`.
        b'A' | b'c' => {
            if value.len() != 1 {
                return Err(format!(
                    "the character field '{}' must hold one character",
                    text(field)
                ));
            }
            out.extend_from_slice(&[b'A', value[0]]);
        }
        b'i' | b'I' => push_int(out, num::<i64>(value, "integer value")?, value)?,
        b'f' => {
            out.push(b'f');
            out.extend_from_slice(&float(value, "float value")?.to_le_bytes());
        }
        b'Z' | b'H' => {
            if ty == b'H' && value.len() % 2 != 0 {
                return Err(format!(
                    "the hex field '{}' has an odd number of digits; it needs whole bytes",
                    text(field)
                ));
            }
            out.push(ty);
            out.extend_from_slice(value);
            out.push(0);
        }
        b'B' => {
            let mut parts = value.split(|&b| b == b',');
            let elem = match parts.next() {
                Some([c]) if b"cCsSiIf".contains(c) => *c,
                _ => {
                    return Err(format!(
                        "the B field '{}' needs an element type c, C, s, S, i, I or f",
                        text(field)
                    ));
                }
            };
            out.extend_from_slice(&[b'B', elem]);
            let count_at = out.len();
            out.extend_from_slice(&[0; 4]);
            let mut n = 0u32;
            for p in parts {
                n += 1;
                let bad = |t: &str| format!("'{}' is not a {t} in the B field", text(p));
                match elem {
                    b'c' => out.push(num::<i8>(p, "c element")? as u8),
                    b'C' => out.push(num::<u8>(p, "C element").map_err(|_| bad("C value"))?),
                    b's' => out.extend_from_slice(&num::<i16>(p, "s element")?.to_le_bytes()),
                    b'S' => out.extend_from_slice(&num::<u16>(p, "S element")?.to_le_bytes()),
                    b'i' => out.extend_from_slice(&num::<i32>(p, "i element")?.to_le_bytes()),
                    b'I' => out.extend_from_slice(&num::<u32>(p, "I element")?.to_le_bytes()),
                    _ => out.extend_from_slice(&float(p, "f element")?.to_le_bytes()),
                }
            }
            out[count_at..count_at + 4].copy_from_slice(&n.to_le_bytes());
        }
        other => {
            return Err(format!(
                "the optional field '{}' has the type '{}'; SAM types are A, i, f, Z, H and B (c and I are read as A and i)",
                text(field),
                other as char
            ));
        }
    }
    Ok(())
}

/// What encoding a SAM line found out about it, besides its bytes.
#[derive(Debug, Default)]
pub(crate) struct Encoded {
    /// `ZERO_POS` and `NO_CIGAR`: what htslib would warn about.
    pub notes: u8,
    /// htslib reads the record as a different one: its flag, reference or mate changed, so the
    /// text of the line is not the record's.
    pub normalised: bool,
    /// Why the record cannot be written as BAM (a 32-bit field overflows), when it cannot. Only
    /// set by a lenient encoding, whose bytes are then not a valid BAM record.
    pub not_bam: Option<String>,
}

fn not_bam(reason: &str) -> String {
    format!("this record cannot be written as BAM ({reason}); write SAM with -O sam to keep it")
}

/// A 32-bit field overflows: an error for BAM, a note for SAM.
fn too_big(found: &mut Encoded, strict: bool, reason: String) -> Result<(), String> {
    if strict {
        Err(not_bam(&reason))
    } else {
        found.not_bam.get_or_insert(reason);
        Ok(())
    }
}

/// Append the BAM record (without `block_size`) for the SAM line `line` to `out`, or say why it
/// cannot be written as BAM. `rid` finds a reference by name. On an error `out` is as it was.
pub(crate) fn sam_to_bam(
    line: &[u8],
    rid: &dyn Fn(&[u8]) -> Option<usize>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    encode_sam(line, rid, out, true).map(|_| ())
}

/// As [`sam_to_bam`], with the findings. With `strict` off a record that BAM cannot hold (see
/// [`Encoded::not_bam`]) is not an error, because SAM can hold it.
pub(crate) fn encode_sam(
    line: &[u8],
    rid: &dyn Fn(&[u8]) -> Option<usize>,
    out: &mut Vec<u8>,
    strict: bool,
) -> Result<Encoded, String> {
    let start = out.len();
    let r = encode(line, rid, out, strict);
    if r.is_err() {
        out.truncate(start);
    }
    r
}

thread_local! {
    /// The CIGAR words of the record being encoded, kept between records.
    static CIGAR: std::cell::RefCell<Vec<u32>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn encode(
    line: &[u8],
    rid: &dyn Fn(&[u8]) -> Option<usize>,
    out: &mut Vec<u8>,
    strict: bool,
) -> Result<Encoded, String> {
    let mut cigar = CIGAR.with(|c| std::mem::take(&mut *c.borrow_mut()));
    cigar.clear();
    let r = encode_with(line, rid, out, strict, &mut cigar);
    CIGAR.with(|c| *c.borrow_mut() = cigar);
    r
}

fn encode_with(
    line: &[u8],
    rid: &dyn Fn(&[u8]) -> Option<usize>,
    out: &mut Vec<u8>,
    strict: bool,
    cigar: &mut Vec<u32>,
) -> Result<Encoded, String> {
    let mut found = Encoded::default();
    let mut f: [&[u8]; 11] = [&[]; 11];
    let mut n = 0;
    let mut rest = line;
    while n < 11 {
        match memchr(b'\t', rest) {
            Some(i) => {
                f[n] = &rest[..i];
                rest = &rest[i + 1..];
                n += 1;
            }
            None => {
                f[n] = rest;
                rest = &[];
                n += 1;
                break;
            }
        }
    }
    if n < 11 {
        return Err(format!(
            "the line has {n} tab-separated fields; a SAM record has at least 11"
        ));
    }
    let name = f[0];
    if name.is_empty() || name.len() > 254 {
        return Err(format!(
            "the read name is {} characters long; a SAM read name has 1 to 254",
            name.len()
        ));
    }
    let flag0 = parse_flag(f[1])?;
    let named_ref = if f[2] == b"*" {
        -1i32
    } else {
        rid(f[2]).ok_or_else(|| {
            format!(
                "RNAME '{}' is not among the @SQ lines of the header; add an @SQ line for it, \
fix the name, or pass --skip-invalid to skip such reads",
                text(f[2])
            )
        })? as i32
    };
    let pos1: u32 = num(f[3], "POS")?;
    if pos1 > i32::MAX as u32 + 1 {
        return Err(format!(
            "POS {pos1} is above 2147483648; BAM positions are 32-bit"
        ));
    }
    let pos = i64::from(pos1) - 1;
    let mapq: u8 = num(f[4], "MAPQ")?;
    let pnext: u32 = num(f[7], "PNEXT")?;
    if pnext > i32::MAX as u32 + 1 {
        too_big(
            &mut found,
            strict,
            format!("PNEXT {pnext} is above 2147483648"),
        )?;
    }
    let mate_pos = i64::from(pnext) - 1;
    let tlen: i64 = num(f[8], "TLEN")?;
    if i32::try_from(tlen).is_err() {
        too_big(
            &mut found,
            strict,
            format!("TLEN {tlen} does not fit in 32 bits"),
        )?;
    }

    // The CIGAR, and what it covers.
    let (mut rlen, mut qlen) = (0u64, 0u64);
    if f[5] != b"*" {
        let mut c = f[5];
        while !c.is_empty() {
            let digits = c.iter().take_while(|b| b.is_ascii_digit()).count();
            let op = c
                .get(digits)
                .and_then(|l| CIGAR_LETTERS.iter().position(|x| x == l));
            let (Some(op), true) = (op, (1..=9).contains(&digits)) else {
                return Err(format!(
                    "CIGAR '{}' is not a list of operations such as 10M2D5M",
                    text(f[5])
                ));
            };
            let n: u32 = num(&c[..digits], "CIGAR length")?;
            if n >= 1 << 28 {
                return Err(format!("the CIGAR length {n} is too large for BAM"));
            }
            if matches!(op, 0 | 2 | 3 | 7 | 8) {
                rlen += u64::from(n);
            }
            if matches!(op, 0 | 1 | 4 | 7 | 8) {
                qlen += u64::from(n);
            }
            cigar.push(n << 4 | op as u32);
            c = &c[digits + 1..];
        }
    }
    let seq = f[9];
    let l_seq = if seq == b"*" { 0 } else { seq.len() };
    if l_seq > 0 && !cigar.is_empty() && qlen != l_seq as u64 {
        return Err(format!(
            "the CIGAR covers {qlen} bases of the read but the sequence has {l_seq}"
        ));
    }
    // As htslib reads SAM (see `normalise`); the mate's `=` is the read's reference as normalised.
    let (flag, ref_id, notes) = normalise(flag0, named_ref, pos, !cigar.is_empty());
    found.notes = notes;
    let mut mate = match f[6] {
        b"*" => -1i32,
        b"=" => ref_id,
        n => rid(n).ok_or_else(|| {
            format!(
                "RNEXT '{}' is not among the @SQ lines of the header; add an @SQ line for it, \
fix the name, or pass --skip-invalid to skip such reads",
                text(n)
            )
        })? as i32,
    };
    if mate_pos < 0 {
        mate = -1;
    }
    // What the mate would be without htslib's rules: the line's own reading of RNEXT.
    let written_mate = match f[6] {
        b"*" => -1,
        b"=" => named_ref,
        n => rid(n).map_or(-1, |i| i as i32),
    };
    found.normalised = flag != flag0 || ref_id != named_ref || mate != written_mate;
    let end = if flag & FLAG_UNMAPPED != 0 {
        pos + 1
    } else {
        pos + (rlen.max(1) as i64)
    };
    // More operations than the 16-bit count holds: a placeholder, and the array in `CG`.
    let long = cigar.len() > usize::from(u16::MAX);
    if long && (rlen >= 1 << 28 || l_seq >= 1 << 28) {
        let what = if rlen >= 1 << 28 {
            "reference span"
        } else {
            "sequence"
        };
        too_big(
            &mut found,
            strict,
            format!(
                "the CIGAR has {} operations, which BAM keeps in a CG tag, and the read's {what} \
is 2^28 or more",
                cigar.len()
            ),
        )?;
    }
    let placeholder = [(l_seq as u32) << 4 | 4, (rlen as u32) << 4 | 3];
    let inline: &[u32] = if long { &placeholder } else { cigar };

    out.extend_from_slice(&ref_id.to_le_bytes());
    out.extend_from_slice(&(pos as i32).to_le_bytes());
    out.push(name.len() as u8 + 1);
    out.push(mapq);
    out.extend_from_slice(&reg2bin(pos, end).to_le_bytes());
    out.extend_from_slice(&(inline.len() as u16).to_le_bytes());
    out.extend_from_slice(&flag.to_le_bytes());
    out.extend_from_slice(&(l_seq as u32).to_le_bytes());
    out.extend_from_slice(&mate.to_le_bytes());
    out.extend_from_slice(&(mate_pos as i32).to_le_bytes());
    out.extend_from_slice(&(tlen as i32).to_le_bytes());
    out.extend_from_slice(name);
    out.push(0);
    for w in inline {
        out.extend_from_slice(&w.to_le_bytes());
    }
    for pair in seq[..l_seq].chunks(2) {
        out.push(
            SEQ_CODE[usize::from(pair[0])] << 4
                | pair.get(1).map_or(0, |&b| SEQ_CODE[usize::from(b)]),
        );
    }
    let qual = f[10];
    if qual == b"*" {
        out.extend(std::iter::repeat_n(0xffu8, l_seq));
    } else {
        if qual.len() != l_seq {
            return Err(format!(
                "the quality string has {} characters but the sequence has {l_seq}",
                qual.len()
            ));
        }
        for &q in qual {
            if !(33..=126).contains(&q) {
                return Err(format!(
                    "the quality character {:?} is outside ! to ~",
                    q as char
                ));
            }
            out.push(q - 33);
        }
    }
    // The optional fields (`rest` is what follows the eleventh field).
    if !rest.is_empty() {
        for field in rest.split(|&b| b == b'\t') {
            push_aux(out, field)?;
        }
    }
    if long {
        out.extend_from_slice(b"CGBI");
        out.extend_from_slice(&(cigar.len() as u32).to_le_bytes());
        for w in cigar.iter() {
            out.extend_from_slice(&w.to_le_bytes());
        }
    }
    Ok(found)
}
