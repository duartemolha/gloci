//! SAM text to a BAM record, with the choices samtools makes (D104): the smallest integer
//! type that holds each `i` aux value, `bin` from the footprint, and a `CG` array when the
//! CIGAR has more than 65,535 operations. The output is byte for byte what `samtools view -b`
//! writes for the same line.

use memchr::memchr;

const FLAG_UNMAPPED: u16 = 0x4;
const SEQ_CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// htslib's `hts_reg2bin(beg, end, 14, 5)`.
fn reg2bin(beg: i64, end: i64) -> u16 {
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
        b'A' => {
            if value.len() != 1 {
                return Err(format!(
                    "the A field '{}' must hold one character",
                    text(field)
                ));
            }
            out.extend_from_slice(&[b'A', value[0]]);
        }
        b'i' => push_int(out, num::<i64>(value, "integer value")?, value)?,
        b'f' => {
            out.push(b'f');
            out.extend_from_slice(&float(value, "float value")?.to_le_bytes());
        }
        b'Z' | b'H' => {
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
                "the optional field '{}' has the type '{}'; SAM types are A, i, f, Z, H and B",
                text(field),
                other as char
            ));
        }
    }
    Ok(())
}

/// Append the BAM record (without `block_size`) for the SAM line `line` to `out`. `rid` finds
/// a reference by name. On an error `out` is as it was.
pub(crate) fn sam_to_bam(
    line: &[u8],
    rid: &dyn Fn(&[u8]) -> Option<usize>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let start = out.len();
    let r = encode(line, rid, out);
    if r.is_err() {
        out.truncate(start);
    }
    r
}

fn encode(
    line: &[u8],
    rid: &dyn Fn(&[u8]) -> Option<usize>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let mut f: Vec<&[u8]> = Vec::with_capacity(16);
    let mut rest = line;
    while f.len() < 11 {
        match memchr(b'\t', rest) {
            Some(i) => {
                f.push(&rest[..i]);
                rest = &rest[i + 1..];
            }
            None => {
                f.push(rest);
                rest = &[];
                break;
            }
        }
    }
    if f.len() < 11 {
        return Err(format!(
            "the line has {} tab-separated fields; a SAM record has at least 11",
            f.len()
        ));
    }
    let name = f[0];
    if name.is_empty() || name.len() > 254 {
        return Err(format!(
            "the read name is {} characters long; a SAM read name has 1 to 254",
            name.len()
        ));
    }
    let mut flag: u16 = num(f[1], "FLAG")?;
    let mut ref_id = if f[2] == b"*" {
        -1i32
    } else {
        rid(f[2]).ok_or_else(|| {
            format!(
                "RNAME '{}' is not among the @SQ lines of the header",
                text(f[2])
            )
        })? as i32
    };
    let pos1: u32 = num(f[3], "POS")?;
    if pos1 > i32::MAX as u32 {
        return Err(format!("POS {pos1} is above 2147483647"));
    }
    let pos = i64::from(pos1) - 1;
    let mapq: u8 = num(f[4], "MAPQ")?;
    let mut mate = match f[6] {
        b"*" => -1i32,
        b"=" => ref_id,
        n => rid(n).ok_or_else(|| {
            format!(
                "RNEXT '{}' is not among the @SQ lines of the header",
                text(n)
            )
        })? as i32,
    };
    let mate_pos = i64::from(num::<u32>(f[7], "PNEXT")?) - 1;
    // As htslib reads SAM: a reference without a position is no reference (for the read and
    // for its mate), and a read with no reference is unmapped.
    if pos < 0 {
        ref_id = -1;
    }
    if mate_pos < 0 {
        mate = -1;
    }
    if ref_id < 0 {
        flag |= FLAG_UNMAPPED;
    }
    let tlen: i32 = num(f[8], "TLEN")?;

    // The CIGAR, and what it covers.
    let mut cigar: Vec<u32> = Vec::new();
    let (mut rlen, mut qlen) = (0u64, 0u64);
    if f[5] != b"*" {
        let mut c = f[5];
        while !c.is_empty() {
            let digits = c.iter().take_while(|b| b.is_ascii_digit()).count();
            let op = c
                .get(digits)
                .and_then(|l| b"MIDNSHP=X".iter().position(|x| x == l));
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
    if l_seq > 0 && !cigar.is_empty() && flag & FLAG_UNMAPPED == 0 && qlen != l_seq as u64 {
        return Err(format!(
            "the CIGAR covers {qlen} bases of the read but the sequence has {l_seq}"
        ));
    }
    let end = if flag & FLAG_UNMAPPED != 0 {
        pos + 1
    } else {
        pos + (rlen.max(1) as i64)
    };
    // More operations than the 16-bit count holds: a placeholder, and the array in `CG`.
    let long = cigar.len() > usize::from(u16::MAX);
    let inline: Vec<u32> = if long {
        vec![(l_seq as u32) << 4 | 4, (rlen as u32) << 4 | 3]
    } else {
        cigar.clone()
    };

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
    out.extend_from_slice(&tlen.to_le_bytes());
    out.extend_from_slice(name);
    out.push(0);
    for w in &inline {
        out.extend_from_slice(&w.to_le_bytes());
    }
    for pair in seq[..l_seq].chunks(2) {
        let code = |b: u8| {
            SEQ_CODES
                .iter()
                .position(|c| c.eq_ignore_ascii_case(&b))
                .unwrap_or(15) as u8
        };
        out.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |&b| code(b)));
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
        for w in &cigar {
            out.extend_from_slice(&w.to_le_bytes());
        }
    }
    Ok(())
}
