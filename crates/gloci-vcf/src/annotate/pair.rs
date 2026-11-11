//! Pairing of the records of two files at one position, as htslib's synced
//! reader pairs them under `--pair-logic some` (bcftools annotate's default):
//! two records pair when they share a `REF>ALT` allele. The algorithm is
//! `bcf_sr_sort_set` of htslib 1.21, including the order of the rows, which
//! decides the order of records written at one position.

use std::collections::HashMap;

pub(crate) const SR_REF: u8 = 1;
pub(crate) const SR_SNP: u8 = 2;
pub(crate) const SR_INDEL: u8 = 4;
pub(crate) const SR_OTHER: u8 = 8;

// htslib's variant type bits.
const VCF_SNP: u32 = 1;
const VCF_MNP: u32 = 2;
const VCF_INDEL: u32 = 4;
const VCF_OTHER: u32 = 8;
const VCF_BND: u32 = 16;
const VCF_OVERLAP: u32 = 32;
const VCF_INS: u32 = 64;
const VCF_DEL: u32 = 128;
/// `bcf_get_variant_types` hides INS and DEL.
const ORIG_TYPES: u32 = VCF_SNP | VCF_MNP | VCF_INDEL | VCF_OTHER | VCF_BND | VCF_OVERLAP;

fn up(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

/// htslib's `bcf_set_variant_type` of one ALT against REF.
pub(crate) fn alt_type(reference: &[u8], alt: &[u8]) -> u32 {
    if alt == b"*" {
        return VCF_OVERLAP;
    }
    if reference.len() == 1 && alt.len() == 1 {
        if alt[0] == b'.' || reference[0] == alt[0] || alt[0] == b'X' {
            return 0;
        }
        return VCF_SNP;
    }
    if alt.first() == Some(&b'<') {
        if alt.starts_with(b"<X>") || alt.starts_with(b"<*>") || alt == b"<NON_REF>" {
            return 0;
        }
        return VCF_OTHER;
    }
    if matches!(alt.first(), Some(b']' | b'[')) {
        return VCF_BND;
    }
    let (mut r, mut a) = (0usize, 0usize);
    while r < reference.len() && a < alt.len() && up(reference[r]) == up(alt[a]) {
        r += 1;
        a += 1;
    }
    let (rest_r, rest_a) = (&reference[r..], &alt[a..]);
    if !rest_a.is_empty() && rest_r.is_empty() {
        if matches!(rest_a[0], b']' | b'[') {
            return VCF_BND;
        }
        return VCF_INDEL | VCF_INS;
    } else if !rest_r.is_empty() && rest_a.is_empty() {
        return VCF_INDEL | VCF_DEL;
    } else if rest_r.is_empty() && rest_a.is_empty() {
        return 0;
    }
    // Both have unmatched bases: compare from the end.
    let (mut re, mut ae) = (reference.len() - 1, alt.len() - 1);
    while re > r && ae > a && up(reference[re]) == up(alt[ae]) {
        re -= 1;
        ae -= 1;
    }
    if ae == a {
        if re == r {
            return VCF_SNP;
        }
        if up(reference[re]) == up(alt[ae]) {
            return VCF_INDEL | VCF_DEL;
        }
        return VCF_OTHER;
    } else if re == r {
        if up(reference[re]) == up(alt[ae]) {
            return VCF_INDEL | VCF_INS;
        }
        return VCF_OTHER;
    }
    if re - r == ae - a { VCF_MNP } else { VCF_OTHER }
}

/// What the pairing needs of one record.
#[derive(Clone, Debug)]
pub struct PairRec {
    /// `REF>ALT` per ALT, comma-separated (`REF>.` without an ALT); a symbolic
    /// ALT is followed by `/END` when INFO/END is present.
    pub key: Vec<u8>,
    /// The SR type bits.
    pub ty: u8,
}

/// `info_end` is INFO/END of the record, if any.
pub fn pair_rec(reference: &[u8], alts: &[&[u8]], info_end: Option<i64>) -> PairRec {
    let mut key = Vec::new();
    let mut vt = 0u32;
    for (i, a) in alts.iter().enumerate() {
        if i > 0 {
            key.push(b',');
        }
        key.extend_from_slice(reference);
        key.push(b'>');
        key.extend_from_slice(a);
        if a.first() == Some(&b'<')
            && let Some(e) = info_end
            && e != 0
        {
            key.push(b'/');
            key.extend_from_slice(e.to_string().as_bytes());
        }
        vt |= alt_type(reference, a);
    }
    if alts.is_empty() {
        key.extend_from_slice(reference);
        key.extend_from_slice(b">.");
    }
    let vt = vt & ORIG_TYPES;
    let ty = if vt == 0 {
        SR_REF
    } else {
        let mut t = 0;
        if vt & (VCF_SNP | VCF_MNP) != 0 {
            t |= SR_SNP;
        }
        if vt & VCF_INDEL != 0 {
            t |= SR_INDEL;
        }
        if vt & VCF_OTHER != 0 {
            t |= SR_OTHER;
        }
        t
    };
    PairRec { key, ty }
}

struct Var {
    ty: u8,
    key: Vec<u8>,
    /// `(reader, record index)` of every record with this variant.
    recs: Vec<(usize, usize)>,
    /// The groups (one per distinct set of records of a reader) it is in.
    mask: u64,
}

struct VSet {
    vars: Vec<usize>,
    mask: u64,
    cnt: usize,
}

fn is_subset(a: &[u8], b: &[u8]) -> bool {
    a.split(|&c| c == b',').filter(|p| !p.is_empty()).any(|pa| {
        b.split(|&c| c == b',')
            .any(|pb| pa.eq_ignore_ascii_case(pb))
    })
}

/// `pairing_score` under `some`: pair on a shared allele, else not. htslib
/// scores every other combination zero and stops at the first one it meets,
/// so only the first variant of each set is compared.
fn paired(vars: &[Var], iv: &VSet, jv: &VSet) -> bool {
    let (Some(&i), Some(&j)) = (iv.vars.first(), jv.vars.first()) else {
        return false;
    };
    let (a, b) = (&vars[i], &vars[j]);
    (a.ty == b.ty && a.key == b.key) || (a.ty & b.ty != 0 && is_subset(&a.key, &b.key))
}

/// The rows of a position that both files have records at: each row has the
/// index of a record of `file`, of `src`, or both, in bcftools' order.
pub fn pair_rows(file: &[PairRec], src: &[PairRec]) -> Vec<(Option<usize>, Option<usize>)> {
    let readers = [file, src];
    let mut vars: Vec<Var> = Vec::new();
    let mut by_key: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut groups: HashMap<Vec<Vec<u8>>, u32> = HashMap::new();
    let mut reader_vars: Vec<Vec<usize>> = Vec::new();
    let mut reader_group: Vec<u32> = Vec::new();
    for (ireader, recs) in readers.iter().enumerate() {
        let mut vs = Vec::new();
        let mut keys: Vec<Vec<u8>> = Vec::new();
        for (irec, r) in recs.iter().enumerate() {
            let mut key = r.key.clone();
            let mut idx = 0;
            let ivar = loop {
                match by_key.get(&key) {
                    None => break None,
                    Some(&v) => {
                        if vars[v].recs.last().is_some_and(|x| x.0 != ireader) {
                            break Some(v);
                        }
                        key = r.key.clone();
                        key.extend_from_slice(idx.to_string().as_bytes());
                        idx += 1;
                    }
                }
            };
            let ivar = match ivar {
                Some(v) => v,
                None => {
                    vars.push(Var {
                        ty: r.ty,
                        key: r.key.clone(),
                        recs: Vec::new(),
                        mask: 0,
                    });
                    by_key.insert(key.clone(), vars.len() - 1);
                    vars.len() - 1
                }
            };
            vars[ivar].recs.push((ireader, irec));
            vs.push(ivar);
            keys.push(key);
        }
        keys.sort();
        let next = groups.len() as u32;
        let g = *groups.entry(keys).or_insert(next);
        reader_vars.push(vs);
        reader_group.push(g);
    }
    for (ireader, vs) in reader_vars.iter().enumerate() {
        for &v in vs {
            vars[v].mask |= 1 << reader_group[ireader];
        }
    }
    let mut vsets: Vec<VSet> = vars
        .iter()
        .enumerate()
        .map(|(i, v)| VSet {
            vars: vec![i],
            mask: v.mask,
            cnt: v.recs.len(),
        })
        .collect();
    let mut rows = Vec::new();
    while !vsets.is_empty() {
        let mut imax = 0;
        for i in 1..vsets.len() {
            if vsets[imax].cnt < vsets[i].cnt {
                imax = i;
            }
        }
        let mut ipair = None;
        for i in 0..vsets.len() {
            if vsets[imax].mask & vsets[i].mask != 0 {
                continue;
            }
            if paired(&vars, &vsets[imax], &vsets[i]) {
                ipair = Some(i);
                break;
            }
        }
        if let Some(j) = ipair {
            let (lo, hi) = if imax < j { (imax, j) } else { (j, imax) };
            let removed = vsets.remove(hi);
            let keep = &mut vsets[lo];
            keep.mask |= removed.mask;
            keep.vars.extend(removed.vars);
            keep.cnt += removed.cnt;
            continue;
        }
        let vs = vsets.remove(imax);
        let mut row = (None, None);
        for v in vs.vars {
            for &(reader, rec) in &vars[v].recs {
                if reader == 0 {
                    row.0 = Some(rec);
                } else {
                    row.1 = Some(rec);
                }
            }
        }
        rows.push(row);
    }
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(r: &[u8], alts: &[&[u8]]) -> PairRec {
        pair_rec(r, alts, None)
    }

    #[test]
    fn exact_match_goes_first_then_partial() {
        // FILE: C>T,G and C>G; SRC: C>G,T and C>G: the exact pair is a row of
        // its own and comes first, then the partial pair.
        let file = [rec(b"C", &[b"T", b"G"]), rec(b"C", &[b"G"])];
        let src = [rec(b"C", &[b"G", b"T"]), rec(b"C", &[b"G"])];
        assert_eq!(
            pair_rows(&file, &src),
            vec![(Some(1), Some(1)), (Some(0), Some(0))]
        );
    }

    #[test]
    fn unrelated_records_do_not_pair() {
        let file = [rec(b"A", &[b"G"])];
        let src = [rec(b"A", &[b"T"])];
        assert_eq!(
            pair_rows(&file, &src),
            vec![(Some(0), None), (None, Some(0))]
        );
    }

    #[test]
    fn types() {
        assert_eq!(alt_type(b"A", b"G"), VCF_SNP);
        assert_eq!(alt_type(b"A", b"AT"), VCF_INDEL | VCF_INS);
        assert_eq!(alt_type(b"AT", b"A"), VCF_INDEL | VCF_DEL);
        assert_eq!(alt_type(b"AC", b"GT"), VCF_MNP);
        assert_eq!(alt_type(b"AC", b"GTT"), VCF_OTHER);
        assert_eq!(alt_type(b"A", b"<DEL>"), VCF_OTHER);
        assert_eq!(alt_type(b"A", b"*"), VCF_OVERLAP);
        assert_eq!(alt_type(b"A", b"."), 0);
        assert_eq!(alt_type(b"A", b"A]1:5]"), VCF_BND);
    }
}
