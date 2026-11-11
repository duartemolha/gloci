//! `gloci vcf split`: multi-allelic records to biallelic ones, as
//! `bcftools norm -m-` does without a reference (no trimming, no
//! left-alignment). Works on VCF text lines; BCF is rendered first.

use crate::extract::declarations;
use gloci_core::{VariantKind, alt_kind};
use std::collections::{HashMap, HashSet};

/// Which records are split (`bcftools norm -m-snps|indels|both|any`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplitTypes {
    /// Records with at least one SNP ALT.
    Snps,
    /// Records with at least one insertion or deletion ALT.
    Indels,
    /// Every multi-allelic record (as `any`: bcftools splits both alike).
    Both,
    /// Every multi-allelic record.
    Any,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Num {
    A,
    R,
    G,
    Other,
}

fn num_of(number: Option<&[u8]>) -> Num {
    match number {
        Some(b"A") => Num::A,
        Some(b"R") => Num::R,
        Some(b"G") => Num::G,
        _ => Num::Other,
    }
}

/// What happened to a record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// Not split: the caller writes the input record as it is.
    Unchanged,
    /// Split into this many records, appended to the output buffer.
    Split(usize),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Gt,
    Num(Num),
    Drop,
}

/// Splits records against one header's INFO and FORMAT declarations.
pub struct Splitter {
    types: SplitTypes,
    nsamples: usize,
    info: HashMap<Vec<u8>, Num>,
    format: HashMap<Vec<u8>, Num>,
    /// FORMAT keys declared Type=String or Type=Character.
    format_strings: HashSet<Vec<u8>>,
}

fn count_tokens(f: &[u8]) -> usize {
    f.iter().filter(|&&b| b == b',').count() + 1
}

fn nth_token(f: &[u8], n: usize) -> Option<&[u8]> {
    f.split(|&b| b == b',').nth(n)
}

/// Index of genotype `(j, k)`, `j <= k`, in the VCF `Number=G` order.
fn gt_index(j: usize, k: usize) -> usize {
    k * (k + 1) / 2 + j
}

fn field(f: Option<&[u8]>) -> &[u8] {
    match f {
        Some(f) if !f.is_empty() => f,
        _ => b".",
    }
}

impl Splitter {
    /// `header` is the raw VCF header text.
    pub fn new(header: &[u8], types: SplitTypes, nsamples: usize) -> Splitter {
        let map = |kind: &'static [u8]| {
            declarations(header, kind)
                .map(|d| (d.id.to_vec(), num_of(d.number)))
                .collect::<HashMap<_, _>>()
        };
        Splitter {
            types,
            nsamples,
            info: map(b"INFO"),
            format: map(b"FORMAT"),
            format_strings: declarations(header, b"FORMAT")
                .filter(|d| matches!(d.ty, Some(b"String" | b"Character")))
                .map(|d| d.id.to_vec())
                .collect(),
        }
    }

    /// Split `line` (a VCF record without its newline). Split records are
    /// appended to `out` with no newline; `ends` is cleared and gets the end
    /// offset of each. An error says what is wrong with the record.
    pub fn split(
        &self,
        line: &[u8],
        out: &mut Vec<u8>,
        ends: &mut Vec<usize>,
    ) -> Result<Outcome, String> {
        ends.clear();
        let mut cols = line.splitn(10, |&b| b == b'\t');
        let mut fixed: [&[u8]; 9] = [b""; 9];
        let mut nfixed = 0;
        for slot in fixed.iter_mut() {
            match cols.next() {
                Some(c) => {
                    *slot = c;
                    nfixed += 1;
                }
                None => break,
            }
        }
        let rest = cols.next();
        if nfixed < 8 {
            return Err(format!(
                "the record has {nfixed} tab-separated columns, but a VCF record needs at least 8"
            ));
        }
        let alt_col = fixed[4];
        if alt_col == b"." || !alt_col.contains(&b',') {
            return Ok(Outcome::Unchanged);
        }
        let alts: Vec<&[u8]> = alt_col.split(|&b| b == b',').collect();
        let nalt = alts.len();
        let split_it = match self.types {
            SplitTypes::Both | SplitTypes::Any => true,
            SplitTypes::Snps => alts
                .iter()
                .any(|a| alt_kind(0, fixed[3], a) == VariantKind::Snv),
            SplitTypes::Indels => alts.iter().any(|a| {
                matches!(
                    alt_kind(0, fixed[3], a),
                    VariantKind::Ins | VariantKind::Del
                )
            }),
        };
        if !split_it {
            return Ok(Outcome::Unchanged);
        }
        let nals = nalt + 1;
        let dip = nals * (nals + 1) / 2;

        // INFO: validate the A/R/G counts.
        let info_items: Vec<&[u8]> = if fixed[7] == b"." {
            Vec::new()
        } else {
            fixed[7].split(|&b| b == b';').collect()
        };
        let mut info_plan: Vec<Option<Num>> = Vec::with_capacity(info_items.len());
        for item in &info_items {
            let (key, val) = match item.iter().position(|&b| b == b'=') {
                Some(eq) => (&item[..eq], Some(&item[eq + 1..])),
                None => (*item, None),
            };
            let Some(&num) = self.info.get(key) else {
                info_plan.push(None);
                continue;
            };
            if let Some(v) = val {
                let (expected, ok) = match num {
                    Num::A => (nalt, count_tokens(v) == nalt),
                    Num::R => (nals, count_tokens(v) == nals),
                    Num::G => (dip, count_tokens(v) == dip),
                    Num::Other => (0, true),
                };
                if !ok {
                    return Err(format!(
                        "INFO/{} has {} values, but the record has {nalt} ALT alleles and its \
Number={} needs {expected}; correct the record, or declare the field Number=. in the header",
                        String::from_utf8_lossy(key),
                        count_tokens(v),
                        match num {
                            Num::A => "A",
                            Num::R => "R",
                            _ => "G",
                        }
                    ));
                }
            }
            info_plan.push(Some(num));
        }

        // FORMAT plan and per-key row widths.
        let mut fkeys: Vec<(&[u8], Kind)> = Vec::new();
        let mut widths: Vec<usize> = Vec::new();
        let has_format = nfixed == 9;
        if has_format {
            for k in fixed[8].split(|&b| b == b':') {
                let kind = if k == b"GT" {
                    Kind::Gt
                } else {
                    match self.format.get(k) {
                        Some(&n) => Kind::Num(n),
                        None => Kind::Drop,
                    }
                };
                fkeys.push((k, kind));
            }
            widths = vec![0; fkeys.len()];
            let mut count = 0;
            if let Some(rest) = rest {
                for s in rest.split(|&b| b == b'\t') {
                    count += 1;
                    for (i, f) in s.split(|&b| b == b':').enumerate().take(fkeys.len()) {
                        if !f.is_empty() {
                            widths[i] = widths[i].max(count_tokens(f));
                        }
                    }
                }
            }
            if count != self.nsamples {
                return Err(format!(
                    "the record has {count} sample columns, but the header lists {}",
                    self.nsamples
                ));
            }
            for (i, (k, kind)) in fkeys.iter().enumerate() {
                let w = widths[i];
                let bad = match kind {
                    Kind::Num(Num::A) => w != 1 && w != nalt,
                    Kind::Num(Num::R) => w != 1 && w != nals,
                    Kind::Num(Num::G) => w != 1 && w != nals && w != dip,
                    _ => false,
                };
                if bad {
                    return Err(format!(
                        "FORMAT/{} has up to {w} values per sample, but the record has {nalt} \
ALT alleles and its Number={} needs {}; correct the record, or declare the field Number=. in the header",
                        String::from_utf8_lossy(k),
                        match kind {
                            Kind::Num(Num::A) => "A",
                            Kind::Num(Num::R) => "R",
                            _ => "G",
                        },
                        match kind {
                            Kind::Num(Num::A) => nalt.to_string(),
                            Kind::Num(Num::R) => nals.to_string(),
                            _ => format!("{nals} or {dip}"),
                        }
                    ));
                }
            }
        }

        for (i, alt) in alts.iter().enumerate() {
            let a = i + 1;
            for c in &fixed[..4] {
                out.extend_from_slice(c);
                out.push(b'\t');
            }
            out.extend_from_slice(alt);
            out.push(b'\t');
            out.extend_from_slice(fixed[5]);
            out.push(b'\t');
            out.extend_from_slice(fixed[6]);
            out.push(b'\t');
            // INFO
            let mut wrote = false;
            for (item, plan) in info_items.iter().zip(&info_plan) {
                let Some(num) = plan else { continue };
                if wrote {
                    out.push(b';');
                }
                wrote = true;
                match (num, item.iter().position(|&b| b == b'=')) {
                    (Num::Other, _) | (_, None) => out.extend_from_slice(item),
                    (n, Some(eq)) => {
                        out.extend_from_slice(&item[..=eq]);
                        let v = &item[eq + 1..];
                        match n {
                            Num::A => out.extend_from_slice(nth_token(v, i).unwrap_or(b".")),
                            Num::R => {
                                out.extend_from_slice(nth_token(v, 0).unwrap_or(b"."));
                                out.push(b',');
                                out.extend_from_slice(nth_token(v, a).unwrap_or(b"."));
                            }
                            _ => {
                                out.extend_from_slice(nth_token(v, 0).unwrap_or(b"."));
                                for x in [gt_index(0, a), gt_index(a, a)] {
                                    out.push(b',');
                                    out.extend_from_slice(nth_token(v, x).unwrap_or(b"."));
                                }
                            }
                        }
                    }
                }
            }
            if !wrote {
                out.push(b'.');
            }
            if has_format {
                out.push(b'\t');
                let mut first = true;
                for (k, kind) in &fkeys {
                    if *kind != Kind::Drop {
                        if !first {
                            out.push(b':');
                        }
                        first = false;
                        out.extend_from_slice(k);
                    }
                }
                if first {
                    out.push(b'.');
                }
                if let Some(rest) = rest {
                    for s in rest.split(|&b| b == b'\t') {
                        out.push(b'\t');
                        let mut fields = s.split(|&b| b == b':');
                        let mut first = true;
                        for (ki, (key, kind)) in fkeys.iter().enumerate() {
                            let f = fields.next();
                            if *kind == Kind::Drop {
                                continue;
                            }
                            if !first {
                                out.push(b':');
                            }
                            first = false;
                            let f = field(f);
                            let w = widths[ki];
                            match kind {
                                Kind::Gt => recode_gt(f, a, out)?,
                                Kind::Num(_) if w <= 1 => out.extend_from_slice(f),
                                Kind::Num(Num::A) => {
                                    out.extend_from_slice(nth_token(f, i).unwrap_or(b"."))
                                }
                                Kind::Num(Num::R) if self.format_strings.contains(*key) => {
                                    // bcftools keeps a String R value only when
                                    // the sample has this ALT's value.
                                    match (nth_token(f, 0), nth_token(f, a)) {
                                        (Some(r), Some(v)) if a > 0 && count_tokens(f) > a => {
                                            out.extend_from_slice(r);
                                            out.push(b',');
                                            out.extend_from_slice(v);
                                        }
                                        _ => out.push(b'.'),
                                    }
                                }
                                Kind::Num(Num::G) if self.format_strings.contains(*key) => {
                                    let k = count_tokens(f);
                                    if w == dip && k == dip {
                                        out.extend_from_slice(nth_token(f, 0).unwrap_or(b"."));
                                        for x in [gt_index(0, a), gt_index(a, a)] {
                                            out.push(b',');
                                            out.extend_from_slice(nth_token(f, x).unwrap_or(b"."));
                                        }
                                    } else if w == dip {
                                        // Not a diploid vector: copied, except that
                                        // nals missing values are a haploid vector.
                                        if k == nals && f.split(|&b| b == b',').all(|t| t == b".") {
                                            out.extend_from_slice(b".,.");
                                        } else {
                                            out.extend_from_slice(f);
                                        }
                                    } else {
                                        // Haploid row (width nals).
                                        out.extend_from_slice(nth_token(f, 0).unwrap_or(b"."));
                                        if a < k {
                                            out.push(b',');
                                            out.extend_from_slice(nth_token(f, a).unwrap_or(b"."));
                                        }
                                    }
                                }
                                Kind::Num(Num::R) => {
                                    out.extend_from_slice(nth_token(f, 0).unwrap_or(b"."));
                                    out.push(b',');
                                    out.extend_from_slice(nth_token(f, a).unwrap_or(b"."));
                                }
                                Kind::Num(Num::G) => {
                                    let k = count_tokens(f);
                                    let idx = if k == dip {
                                        [0, gt_index(0, a), gt_index(a, a)]
                                    } else {
                                        [0, a, usize::MAX]
                                    };
                                    out.extend_from_slice(nth_token(f, 0).unwrap_or(b"."));
                                    for &x in &idx[1..] {
                                        if x < k {
                                            out.push(b',');
                                            out.extend_from_slice(nth_token(f, x).unwrap_or(b"."));
                                        }
                                    }
                                }
                                Kind::Num(Num::Other) => out.extend_from_slice(f),
                                Kind::Drop => {}
                            }
                        }
                        if first {
                            out.push(b'.');
                        }
                    }
                }
            }
            ends.push(out.len());
        }
        Ok(Outcome::Split(nalt))
    }
}

/// The genotype `gt` with allele `keep` as 1 and every other called allele
/// as 0; missing alleles, separators and so phase and ploidy are kept.
fn recode_gt(gt: &[u8], keep: usize, out: &mut Vec<u8>) -> Result<(), String> {
    let mut i = 0;
    while i < gt.len() {
        let b = gt[i];
        if b == b'/' || b == b'|' || b == b'.' {
            out.push(b);
            i += 1;
        } else if b.is_ascii_digit() {
            let mut n = 0usize;
            while i < gt.len() && gt[i].is_ascii_digit() {
                n = n.saturating_mul(10).saturating_add((gt[i] - b'0') as usize);
                i += 1;
            }
            out.push(if n == keep { b'1' } else { b'0' });
        } else {
            return Err(format!(
                "the genotype `{}` is not alleles joined by / or |; correct the record",
                String::from_utf8_lossy(gt)
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(header: &str, line: &str, t: SplitTypes) -> Vec<String> {
        let ns = header
            .lines()
            .last()
            .unwrap()
            .split('\t')
            .count()
            .saturating_sub(9);
        let s = Splitter::new(header.as_bytes(), t, ns);
        let (mut out, mut ends) = (Vec::new(), Vec::new());
        match s.split(line.as_bytes(), &mut out, &mut ends).unwrap() {
            Outcome::Unchanged => vec![line.to_string()],
            Outcome::Split(_) => {
                let mut prev = 0;
                ends.iter()
                    .map(|&e| {
                        let r = String::from_utf8(out[prev..e].to_vec()).unwrap();
                        prev = e;
                        r
                    })
                    .collect()
            }
        }
    }

    const H: &str = "##INFO=<ID=A1,Number=A,Type=Integer,Description=\"a\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n";

    #[test]
    fn splits_gt_and_pl() {
        let r = run(
            H,
            "c\t1\t.\tA\tC,G\t.\t.\tA1=1,2;U=3\tGT:PL\t1/2:1,2,3,4,5,6",
            SplitTypes::Any,
        );
        assert_eq!(r[0], "c\t1\t.\tA\tC\t.\t.\tA1=1\tGT:PL\t1/0:1,2,3");
        assert_eq!(r[1], "c\t1\t.\tA\tG\t.\t.\tA1=2\tGT:PL\t0/1:1,4,6");
    }

    #[test]
    fn types_gate() {
        let l = "c\t1\t.\tAT\tCG,GA\t.\t.\t.";
        assert_eq!(run(H, l, SplitTypes::Snps).len(), 1);
        assert_eq!(run(H, l, SplitTypes::Both).len(), 2);
    }
}
