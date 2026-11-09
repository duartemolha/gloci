//! Native BCF encoding of the FORMAT `GT` field (BCF 2.2 section 6.3.3).
//!
//! noodles cannot write a GT vector whose length differs between samples
//! (a missing haploid `.` next to diploid calls), so the converter writes
//! this field itself and splices it into the record noodles encoded.

/// Encode the text GT values of all samples as the typed block of one
/// FORMAT field: the type/count descriptor, then the data of every sample.
/// The FORMAT key that precedes the block is the caller's to write.
///
/// Each allele is `(allele + 1) << 1 | phased`; the phase bit comes from the
/// separator before the allele, so the first allele never has it, and a
/// missing allele is 0. Samples with fewer alleles than the largest ploidy
/// are padded with the end-of-vector value. The element width is the
/// smallest of int8, int16 and int32 that holds every value.
///
/// A bare `.` or an empty value (no GT for the sample) is one missing allele. A VCF 4.4 leading
/// separator (`|0|1`) is an error, as it is for bcftools 1.21.
pub fn encode_gt(samples: &[&[u8]], out: &mut Vec<u8>) -> Result<(), String> {
    GtEncoder::default().encode(samples.iter().copied(), out)
}

/// [`encode_gt`] with its working buffers kept between calls, for encoding
/// many records.
#[derive(Default)]
pub struct GtEncoder {
    vals: Vec<i32>,
    lens: Vec<usize>,
}

impl GtEncoder {
    /// Appends the typed GT block of `samples` to `out`; see [`encode_gt`].
    pub fn encode<'a>(
        &mut self,
        samples: impl Iterator<Item = &'a [u8]>,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        let (vals, lens) = (&mut self.vals, &mut self.lens);
        vals.clear();
        lens.clear();
        let mut max_val = 0i32;
        let mut n_samples = 0;
        for (i, s) in samples.enumerate() {
            let before = vals.len();
            parse_sample(s, vals).map_err(|e| format!("sample {}: {e}", i + 1))?;
            lens.push(vals.len() - before);
            max_val = max_val.max(vals[before..].iter().copied().max().unwrap_or(0));
            n_samples = i + 1;
        }
        write_block(vals, lens, max_val, n_samples, out);
        Ok(())
    }
}

fn write_block(vals: &[i32], lens: &[usize], max_val: i32, n_samples: usize, out: &mut Vec<u8>) {
    let ploidy = lens.iter().copied().max().unwrap_or(0);
    let (code, width) = if max_val <= i8::MAX as i32 {
        (1u8, 1usize)
    } else if max_val <= i16::MAX as i32 {
        (2, 2)
    } else {
        (3, 4)
    };
    // Descriptor: count in the high nibble, or 15 and the count as a typed int.
    if ploidy < 15 {
        out.push((ploidy as u8) << 4 | code);
    } else {
        out.push(15 << 4 | code);
        if let Ok(n) = i8::try_from(ploidy) {
            out.extend_from_slice(&[1 << 4 | 1, n as u8]);
        } else if let Ok(n) = i16::try_from(ploidy) {
            out.push(1 << 4 | 2);
            out.extend_from_slice(&n.to_le_bytes());
        } else {
            out.push(1 << 4 | 3);
            out.extend_from_slice(&(ploidy as i32).to_le_bytes());
        }
    }
    out.reserve(n_samples * ploidy * width);
    let mut at = 0;
    for &n in lens {
        for j in 0..ploidy {
            if j < n {
                put(out, code, vals[at + j]);
            } else {
                put(out, code, END[usize::from(code) - 1]);
            }
        }
        at += n;
    }
}

/// End-of-vector values (int8, int16, int32).
const END: [i32; 3] = [-127, -32767, i32::MIN + 1];

/// Marks the integer-missing value in an encoded vector; `put` writes the
/// minimum of the element width (-128, -32768 or `i32::MIN`).
const INT_MISSING: i32 = i32::MIN;

fn put(out: &mut Vec<u8>, code: u8, v: i32) {
    let v = match (v, code) {
        (INT_MISSING, 1) => i32::from(i8::MIN),
        (INT_MISSING, 2) => i32::from(i16::MIN),
        _ => v,
    };
    match code {
        1 => out.push(v as i8 as u8),
        2 => out.extend_from_slice(&(v as i16).to_le_bytes()),
        _ => out.extend_from_slice(&v.to_le_bytes()),
    }
}

/// Append the encoded alleles of one sample's text GT to `vals`.
fn parse_sample(v: &[u8], vals: &mut Vec<i32>) -> Result<(), String> {
    // An empty value (a sample cut short before its GT) is written as htslib
    // does: the first element is the integer-missing value, not the missing
    // allele 0 that `.` gives, and readers see one missing allele either way.
    if v.is_empty() {
        vals.push(INT_MISSING);
        return Ok(());
    }
    if v == b"." {
        vals.push(0);
        return Ok(());
    }
    let bad = |what: &str| {
        format!(
            "GT value {:?} {what}; use alleles separated by / or |, for example 0/1",
            String::from_utf8_lossy(v)
        )
    };
    // htslib (bcftools 1.21) rejects a VCF 4.4 leading separator, so the BCF
    // writer refuses it too rather than guess at the first allele's phase.
    if v.len() > 1 && matches!(v[0], b'|' | b'/') {
        return Err(bad(
            "starts with a phase separator, which BCF output (like bcftools 1.21) does not accept; \
remove the leading separator, or write VCF output instead",
        ));
    }
    let mut rest = v;
    let mut phased = false;
    let mut first = true;
    loop {
        let end = rest
            .iter()
            .position(|&b| b == b'/' || b == b'|')
            .unwrap_or(rest.len());
        let a = &rest[..end];
        let allele: i64 = if a == b"." {
            0
        } else if a.is_empty() {
            return Err(bad("has an empty allele"));
        } else {
            let mut n: i64 = 0;
            for &b in a {
                if !b.is_ascii_digit() {
                    return Err(bad("has an allele that is not a non-negative integer"));
                }
                n = n * 10 + i64::from(b - b'0');
                if n > i64::from(i32::MAX >> 1) - 1 {
                    return Err(bad("has an allele index too large for BCF"));
                }
            }
            n + 1
        };
        let phase = i64::from(!first && phased);
        vals.push(((allele << 1) | phase) as i32);
        first = false;
        if end == rest.len() {
            return Ok(());
        }
        phased = rest[end] == b'|';
        rest = &rest[end + 1..];
    }
}

#[cfg(test)]
mod tests {
    use super::encode_gt;

    fn enc(s: &[&str]) -> Vec<u8> {
        let v: Vec<&[u8]> = s.iter().map(|x| x.as_bytes()).collect();
        let mut out = Vec::new();
        encode_gt(&v, &mut out).unwrap();
        out
    }

    #[test]
    fn mixed_ploidy_pads_with_end_of_vector() {
        // int8, 2 elements: 0/1 -> 2,4; "." -> 0, end.
        assert_eq!(enc(&["0/1", "."]), [0x21, 2, 4, 0, 0x81]);
        assert_eq!(enc(&["1|0"]), [0x21, 4, 3]);
    }

    #[test]
    fn wide_alleles_switch_to_int16() {
        // allele 62 unphased is 126, phased 127: int8. Allele 63 unphased is 128.
        assert_eq!(enc(&["62|62"]), [0x21, 126, 127]);
        assert_eq!(enc(&["63/0", "0"]), [0x22, 128, 0, 2, 0, 2, 0, 0x01, 0x80]);
    }

    #[test]
    fn empty_gt_is_integer_missing_then_padding() {
        // As htslib writes a sample cut short before its GT: 0x80, not 0.
        assert_eq!(enc(&["0/1", ""]), [0x21, 2, 4, 0x80, 0x81]);
    }

    #[test]
    fn invalid_values_are_errors() {
        let mut o = Vec::new();
        assert!(encode_gt(&[b"0/x"], &mut o).is_err());
        assert!(encode_gt(&[b"0//1"], &mut o).is_err());
        assert!(encode_gt(&[b"0/"], &mut o).is_err());
        for lead in [&b"|1"[..], b"|0|1", b"/0|1", b"/1", b"|."] {
            let e = encode_gt(&[lead], &mut o).unwrap_err();
            assert!(e.contains("phase separator") && e.contains("remove"), "{e}");
        }
    }
}
