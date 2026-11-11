//! Which records of the inputs at one position are merged into one record,
//! and in what order the merged records are written: htslib's
//! `bcf_sr_sort_set`/`bcf_sr_sort_next` (the synced reader bcftools merge
//! uses) with the `-m` mode as its pairing logic. Generalises the
//! two-reader port in `annotate::pair` to any number of readers.

use crate::annotate::pair::{SR_INDEL, SR_REF, SR_SNP};
use crate::annotate::{PairRec, pair_rec};
use std::collections::HashMap;

/// `bcftools merge -m`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MergeMode {
    /// Only records with the same alleles merge.
    None,
    Snps,
    Indels,
    Both,
    All,
    /// Records with the same ID merge.
    Id,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
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

struct Var {
    ty: u8,
    key: Vec<u8>,
    nal: Vec<Vec<u8>>,
    /// The same alleles in the record's order.
    ordered: Vec<Vec<u8>>,
    /// `(reader, record index)` of every record with this variant.
    recs: Vec<(usize, usize)>,
    mask: Mask,
}

struct VSet {
    vars: Vec<usize>,
    mask: Mask,
    cnt: usize,
}

/// A record's pairing key and its normalised alleles, sorted and in order.
type PairNorm = (PairRec, (Vec<Vec<u8>>, Vec<Vec<u8>>));

/// What pairing needs of one record.
pub struct PairIn {
    pub reference: Vec<u8>,
    pub alts: Vec<Vec<u8>>,
    pub id: Vec<u8>,
    /// INFO/END, which goes into the key of a symbolic ALT.
    pub end: Option<i64>,
}

fn symbolic(a: &[u8]) -> bool {
    a.first() == Some(&b'<') || a == b"*" || a.contains(&b'[') || a.contains(&b']')
}

/// The alleles of a record with REF extended to `longest` (`REF>ALT` each),
/// so that `A>AT` and `AT>ATT` compare equal.
fn normalised(r: &PairIn, longest: &[u8]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let suffix = if longest.len() >= r.reference.len()
        && longest[..r.reference.len()].eq_ignore_ascii_case(&r.reference)
    {
        &longest[r.reference.len()..]
    } else {
        &[][..]
    };
    let mut v: Vec<Vec<u8>> = r
        .alts
        .iter()
        .map(|a| {
            // A symbolic ALT keeps its own REF: `A><DEL>` and `AT><DEL>` differ.
            let mut x = if symbolic(a) {
                r.reference.clone()
            } else {
                longest.to_vec()
            };
            x.push(b'>');
            x.extend_from_slice(a);
            if !symbolic(a) {
                x.extend_from_slice(suffix);
            }
            x
        })
        .collect();
    let ordered = v.clone();
    v.sort();
    (v, ordered)
}

/// `pairing_score` of two variants, 0 for no pairing: 5 the same record
/// text, 4 the same alleles once REF is extended, 3 an allele in the same
/// place of the lists, 2 any shared allele (all with one type), 1 what
/// `--merge` allows on top.
fn score(mode: MergeMode, a: &Var, b: &Var) -> u32 {
    if a.ty == b.ty && a.key == b.key {
        return 5;
    }
    if mode == MergeMode::Id {
        return 0;
    }
    if a.ty == b.ty {
        if a.nal == b.nal {
            return 4;
        }
        // The same allele in the same place of the list pairs before any other shared allele.
        if a.ordered.iter().zip(&b.ordered).any(|(x, y)| x == y) {
            return 3;
        }
        if a.nal.iter().any(|x| b.nal.contains(x)) {
            return 2;
        }
    }
    let (ta, tb) = (a.ty, b.ty);
    // A record without ALT pairs with any other, a SNP before the rest.
    if ta == SR_REF || tb == SR_REF {
        let other = if ta == SR_REF { tb } else { ta };
        return if other & SR_SNP != 0 { 2 } else { 1 };
    }
    let snp = ta & tb & SR_SNP != 0;
    let indel = ta & tb & SR_INDEL != 0;
    let pure_indel = ta & SR_SNP == 0 && tb & SR_SNP == 0;
    match mode {
        MergeMode::None | MergeMode::Id => 0,
        MergeMode::All => 1,
        MergeMode::Snps => u32::from(snp),
        MergeMode::Indels => u32::from(indel),
        MergeMode::Both => u32::from(snp || (pure_indel && indel)),
    }
}

/// The rows of one position: for each, the record index each reader
/// contributes (at most one per reader), in the order htslib writes them.
pub fn rows(inputs: &[Vec<PairIn>], mode: MergeMode) -> Vec<Vec<Option<usize>>> {
    let nreaders = inputs.len();
    let mut longest: &[u8] = b"";
    for r in inputs.iter().flatten() {
        if r.reference.len() > longest.len() {
            longest = &r.reference;
        }
    }
    let recs: Vec<Vec<PairNorm>> = inputs
        .iter()
        .map(|rs| {
            rs.iter()
                .map(|r| {
                    let alts: Vec<&[u8]> = r.alts.iter().map(Vec::as_slice).collect();
                    let mut p = pair_rec(&r.reference, &alts, r.end);
                    if mode == MergeMode::Id {
                        p = id_pair_rec(&r.id);
                    }
                    (p, normalised(r, longest))
                })
                .collect()
        })
        .collect();
    let mut vars: Vec<Var> = Vec::new();
    let mut by_key: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut groups: HashMap<Vec<Vec<u8>>, usize> = HashMap::new();
    let mut reader_vars: Vec<Vec<usize>> = Vec::new();
    let mut reader_group: Vec<usize> = Vec::new();
    for (ireader, rs) in recs.iter().enumerate() {
        let mut vs = Vec::new();
        let mut keys: Vec<Vec<u8>> = Vec::new();
        for (irec, (r, (nal, ordered))) in rs.iter().enumerate() {
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
                        nal: nal.clone(),
                        ordered: ordered.clone(),
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
        let next = groups.len();
        let g = *groups.entry(keys).or_insert(next);
        reader_vars.push(vs);
        reader_group.push(g);
    }
    for (ireader, vs) in reader_vars.iter().enumerate() {
        for &v in vs {
            vars[v].mask.or(&Mask::bit(reader_group[ireader]));
        }
    }
    // Records without an ALT are paired last, each with a SNP set if there is
    // one, else the first set; they are not partners while the others pair.
    let split_refs = !matches!(mode, MergeMode::All | MergeMode::Id);
    let mut vsets: Vec<VSet> = Vec::new();
    let mut refs: Vec<usize> = Vec::new();
    for (i, v) in vars.iter().enumerate() {
        if split_refs && v.ty == SR_REF {
            refs.push(i);
            continue;
        }
        vsets.push(VSet {
            vars: vec![i],
            mask: v.mask.clone(),
            cnt: v.recs.len(),
        });
    }
    let mut done: Vec<VSet> = Vec::new();
    while !vsets.is_empty() {
        let mut imax = 0;
        for i in 1..vsets.len() {
            if vsets[imax].cnt < vsets[i].cnt {
                imax = i;
            }
        }
        let mut ipair = None;
        let mut best = 0;
        for i in 0..vsets.len() {
            if i == imax || vsets[imax].mask.meets(&vsets[i].mask) {
                continue;
            }
            // Two sets pair as well as their worst pair of variants does.
            let mut s = u32::MAX;
            for &a in &vsets[imax].vars {
                for &b in &vsets[i].vars {
                    s = s.min(score(mode, &vars[a], &vars[b]));
                }
            }
            if s > best {
                best = s;
                ipair = Some(i);
            }
        }
        if let Some(j) = ipair {
            let (lo, hi) = if imax < j { (imax, j) } else { (j, imax) };
            let removed = vsets.remove(hi);
            let keep = &mut vsets[lo];
            keep.mask.or(&removed.mask);
            keep.vars.extend(removed.vars);
            keep.cnt += removed.cnt;
            continue;
        }
        done.push(vsets.remove(imax));
    }
    if split_refs {
        done.sort_by_key(|s| s.vars.iter().copied().min().unwrap_or(0));
    }
    for r in refs {
        let m = &vars[r].mask;
        let first = done.iter().position(|s| !s.mask.meets(m));
        let pick = if mode == MergeMode::Indels {
            first
        } else {
            done.iter()
                .position(|s| !s.mask.meets(m) && s.vars.iter().any(|&v| vars[v].ty & SR_SNP != 0))
                .or(first)
        };
        match pick {
            Some(k) => {
                done[k].vars.push(r);
                done[k].mask.or(m);
            }
            None => done.push(VSet {
                vars: vec![r],
                mask: m.clone(),
                cnt: 1,
            }),
        }
    }
    let mut rows = Vec::new();
    for vs in done {
        let mut row = vec![None; nreaders];
        let lo = vs.vars.iter().copied().min().unwrap_or(0);
        let snp = vs.vars.iter().any(|&v| vars[v].ty & SR_SNP != 0);
        for v in vs.vars {
            for &(reader, rec) in &vars[v].recs {
                row[reader] = Some(rec);
            }
        }
        rows.push((lo, snp, row));
    }
    // Records with a SNP go first, in the order their variants were first
    // seen; with `-m indels` the order is just that.
    if mode == MergeMode::All {
        // As the sets were formed.
    } else if mode == MergeMode::Indels {
        rows.sort_by_key(|r| r.0);
    } else {
        rows.sort_by_key(|r| (!r.1, r.0));
    }
    rows.into_iter().map(|r| r.2).collect()
}

/// What the pairing needs of a record under `-m id`: its ID alone.
pub fn id_pair_rec(id: &[u8]) -> PairRec {
    PairRec {
        key: id.to_vec(),
        ty: SR_REF,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(r: &str, alts: &[&str]) -> PairIn {
        PairIn {
            reference: r.as_bytes().to_vec(),
            alts: alts.iter().map(|a| a.as_bytes().to_vec()).collect(),
            id: b".".to_vec(),
            end: None,
        }
    }

    #[test]
    fn modes_decide_what_pairs() {
        let files = vec![vec![rec("A", &["C"])], vec![rec("A", &["G"])]];
        assert_eq!(rows(&files, MergeMode::None).len(), 2);
        assert_eq!(rows(&files, MergeMode::Indels).len(), 2);
        for m in [MergeMode::Snps, MergeMode::Both, MergeMode::All] {
            assert_eq!(rows(&files, m), vec![vec![Some(0), Some(0)]]);
        }
        let mixed = vec![vec![rec("A", &["G"])], vec![rec("AT", &["A"])]];
        assert_eq!(rows(&mixed, MergeMode::Both).len(), 2);
        assert_eq!(rows(&mixed, MergeMode::All).len(), 1);
    }

    #[test]
    fn snp_rows_come_first_and_a_shared_allele_pairs() {
        let files = vec![
            vec![rec("A", &["<DEL>"])],
            vec![rec("A", &["T", "G"])],
            vec![rec("A", &["T", "C"])],
            vec![rec("A", &["G"])],
        ];
        let r = rows(&files, MergeMode::None);
        // T,G pairs with T,C (same first ALT); G stays alone; the SNP rows go before <DEL>
        assert_eq!(
            r,
            vec![
                vec![None, Some(0), Some(0), None],
                vec![None, None, None, Some(0)],
                vec![Some(0), None, None, None],
            ]
        );
    }

    #[test]
    fn same_reference_extended_alleles_compare_equal() {
        // A>AT and AT>ATT are one variant: the exact pair is taken first.
        let files = vec![
            vec![rec("A", &["AT"]), rec("AT", &["ATT"])],
            vec![rec("AT", &["ATT"]), rec("A", &["AT"])],
        ];
        let r = rows(&files, MergeMode::None);
        assert_eq!(r, vec![vec![Some(0), Some(1)], vec![Some(1), Some(0)]]);
    }

    #[test]
    fn many_groups_do_not_overflow_the_masks() {
        // 200 files with 200 different records at one position
        let files: Vec<Vec<PairIn>> = (0..200)
            .map(|i| {
                let alt = ["C", "G", "T"][i % 3].to_string();
                let mut r = rec("A", &[&alt]);
                r.id = format!("id{i}").into_bytes();
                vec![r]
            })
            .collect();
        let r = rows(&files, MergeMode::Id);
        assert_eq!(r.len(), 200);
    }
}
