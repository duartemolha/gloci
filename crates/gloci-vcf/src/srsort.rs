//! htslib's `bcf_sr_sort` (synced_bcf_reader, 1.21): the pairing of the
//! records of several files at one position into rows, as `bcftools annotate`
//! (`--pair-logic some`) uses it for the order of its output and as
//! `bcftools merge` uses it (exact pairing) to pick the first record of
//! each file. Shared by `annotate` and `merge`.

use std::collections::HashMap;

pub(crate) const SR_REF: u8 = 1;
pub(crate) const SR_SNP: u8 = 2;
pub(crate) const SR_INDEL: u8 = 4;
pub(crate) const SR_OTHER: u8 = 8;

// htslib's variant type bits.
pub(crate) const VCF_SNP: u32 = 1;
pub(crate) const VCF_MNP: u32 = 2;
pub(crate) const VCF_INDEL: u32 = 4;
pub(crate) const VCF_OTHER: u32 = 8;
pub(crate) const VCF_BND: u32 = 16;
pub(crate) const VCF_OVERLAP: u32 = 32;
pub(crate) const VCF_INS: u32 = 64;
pub(crate) const VCF_DEL: u32 = 128;
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
    /// Number of ALT alleles.
    pub nalt: usize,
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
    PairRec {
        nalt: alts.len(),
        key,
        ty,
    }
}

/// htslib's `pair` logic as far as it is used: `Exact` (bcftools merge) and
/// `Some` (`--pair-logic some`, annotate).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pair {
    Exact,
    Some,
}

struct Var {
    ty: u8,
    nalt: usize,
    key: Vec<u8>,
    /// `(reader, record index)` of every record with this variant.
    recs: Vec<(usize, usize)>,
    /// The groups (one per distinct set of records of a reader) it is in.
    mask: Mask,
}

struct VSet {
    vars: Vec<usize>,
    mask: Mask,
    cnt: usize,
}

#[derive(Clone, Debug, Default)]
struct Mask(Vec<u64>);

impl Mask {
    fn bit(i: usize) -> Mask {
        let mut v = vec![0u64; i / 64 + 1];
        v[i / 64] = 1 << (i % 64);
        Mask(v)
    }

    fn or(&mut self, o: &Mask) {
        if self.0.len() < o.0.len() {
            self.0.resize(o.0.len(), 0);
        }
        for (a, b) in self.0.iter_mut().zip(&o.0) {
            *a |= b;
        }
    }

    fn meets(&self, o: &Mask) -> bool {
        self.0.iter().zip(&o.0).any(|(a, b)| a & b != 0)
    }
}

fn alleles(k: &[u8]) -> impl Iterator<Item = &[u8]> {
    k.split(|&c| c == b',')
}

/// `multi_is_subset`: some allele of `a` is one of `b`.
fn is_subset(a: &[u8], b: &[u8]) -> bool {
    alleles(a)
        .filter(|p| !p.is_empty())
        .any(|pa| alleles(b).any(|pb| pa.eq_ignore_ascii_case(pb)))
}

/// `multi_is_exact`: the same number of alleles, each of `a` among those of `b`.
fn is_exact(a: &Var, b: &Var) -> bool {
    a.nalt == b.nalt
        && a.key.len() == b.key.len()
        && alleles(&a.key)
            .filter(|p| !p.is_empty())
            .all(|pa| alleles(&b.key).any(|pb| pa.eq_ignore_ascii_case(pb)))
}

/// `pairing_score`: 0 for no pairing, anything else for the best kind (only the
/// `exact` and `some` logics exist here, whose scores are all or nothing).
fn pairing_score(vars: &[Var], iv: &VSet, jv: &VSet, pair: Pair) -> u32 {
    for &i in &iv.vars {
        for &j in &jv.vars {
            let (a, b) = (&vars[i], &vars[j]);
            match pair {
                Pair::Exact => {
                    if a.ty != b.ty {
                        continue;
                    }
                    if a.key == b.key || is_exact(a, b) {
                        return u32::MAX;
                    }
                }
                Pair::Some => {
                    if a.ty == b.ty && a.key == b.key {
                        return u32::MAX;
                    }
                    if a.ty & b.ty != 0 && is_subset(&a.key, &b.key) {
                        return u32::MAX;
                    }
                    // every other combination scores zero: the sets do not pair
                    return 0;
                }
            }
        }
    }
    0
}

/// The rows of one position, in the order htslib writes them: each row holds
/// the `(reader, record)` pairs paired together (at most one record per reader).
pub fn sort_rows(readers: &[&[PairRec]], pair: Pair) -> Vec<Vec<(usize, usize)>> {
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
                        nalt: r.nalt,
                        key: r.key.clone(),
                        recs: Vec::new(),
                        mask: Mask::default(),
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
            vars[v].mask.or(&Mask::bit(reader_group[ireader] as usize));
        }
    }
    let mut vsets: Vec<VSet> = vars
        .iter()
        .enumerate()
        .map(|(i, v)| VSet {
            vars: vec![i],
            mask: v.mask.clone(),
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
        let mut max_score = 0;
        for i in 0..vsets.len() {
            if vsets[imax].mask.meets(&vsets[i].mask) {
                continue;
            }
            let score = pairing_score(&vars, &vsets[imax], &vsets[i], pair);
            if max_score < score {
                max_score = score;
                ipair = Some(i);
            }
        }
        if let Some(j) = ipair
            && j != imax
        {
            let (lo, hi) = if imax < j { (imax, j) } else { (j, imax) };
            let removed = vsets.remove(hi);
            let keep = &mut vsets[lo];
            keep.mask.or(&removed.mask);
            keep.vars.extend(removed.vars);
            keep.cnt += removed.cnt;
            continue;
        }
        let vs = vsets.remove(imax);
        let mut row = Vec::new();
        for v in vs.vars {
            row.extend(vars[v].recs.iter().copied());
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
    fn exact_needs_the_same_alleles_some_a_shared_one() {
        let f0 = [rec(b"A", &[b"C", b"G"])];
        let f1 = [rec(b"A", &[b"G"])];
        let f2 = [rec(b"A", &[b"G", b"C"])];
        assert_eq!(sort_rows(&[&f0, &f1], Pair::Exact).len(), 2);
        assert_eq!(sort_rows(&[&f0, &f1], Pair::Some).len(), 1);
        // the order of the alleles does not matter for exact
        assert_eq!(sort_rows(&[&f0, &f2], Pair::Exact).len(), 1);
    }

    #[test]
    fn the_set_with_most_records_is_the_first_row() {
        let f0 = [rec(b"A", &[b"T"]), rec(b"A", &[b"C"])];
        let f1 = [rec(b"A", &[b"C"])];
        let f2 = [rec(b"A", &[b"C"])];
        let rows = sort_rows(&[&f0, &f1, &f2], Pair::Exact);
        assert_eq!(rows[0], vec![(0, 1), (1, 0), (2, 0)]);
        assert_eq!(rows[1], vec![(0, 0)]);
    }
}
