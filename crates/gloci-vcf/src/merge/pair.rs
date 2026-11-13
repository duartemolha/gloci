//! Which records of the inputs at one position are merged into one record,
//! and in what order: bcftools merge 1.21's own loop (`maux_reset`,
//! `can_merge`, `stage_line`, `types_compatible`, `merge_alleles`) after the
//! synced reader's exact pairing ([`crate::srsort`]) has chosen the first
//! record of each input.

use super::Alleles;
use crate::annotate::vcmp::Vcmp;
use crate::srsort::{Pair, VCF_INDEL, VCF_MNP, VCF_SNP, alt_type, pair_rec, sort_rows};

/// `bcftools merge -m`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MergeMode {
    /// Only records with the same alleles (or a shared allele of one type) merge.
    None,
    Snps,
    Indels,
    Both,
    All,
    /// Records with the same ID merge.
    Id,
}

/// What the merge needs of one record.
pub struct PairIn {
    pub reference: Vec<u8>,
    pub alts: Vec<Vec<u8>>,
    pub id: Vec<u8>,
    /// INFO/END, which goes into the key of a symbolic ALT.
    pub end: Option<i64>,
}

/// One merged record to write: the record each input contributes and the alleles.
pub struct StagedRow {
    pub members: Vec<Option<usize>>,
    pub alleles: Alleles,
}

const SKIP_DONE: u8 = 1;
const SKIP_DIFF: u8 = 2;
const SNP_MASK: u32 = (VCF_SNP << 1) | (VCF_MNP << 1);
const INDEL_MASK: u32 = VCF_INDEL << 1;
const REF_MASK: u32 = 1;
const COLLAPSE_SNPS: u8 = 1;
const COLLAPSE_INDELS: u8 = 2;
const COLLAPSE_ANY: u8 = 4;

struct Rec<'a> {
    pin: &'a PairIn,
    orig: usize,
    /// REF first.
    alleles: Vec<&'a [u8]>,
    /// Variant type bits per allele (0 for REF).
    types: Vec<u32>,
    /// The record's types shifted by one, or `REF_MASK`.
    ty: u32,
    skip: u8,
    /// Allele map into the merged list (valid for candidates).
    map: Vec<usize>,
}

struct Buf<'a> {
    recs: Vec<Rec<'a>>,
    /// OR of the types of its records not yet done.
    var_types: u32,
    cur: Option<usize>,
}

/// `merge_alleles` of vcfmerge.c: merge `a` into `b`, returning the map.
fn merge_alleles(a: &[Vec<u8>], b: &mut Vec<Vec<u8>>) -> Result<Vec<usize>, String> {
    let na = a.len();
    let mut map = vec![0usize; na];
    let rla = a[0].len().max(1);
    let rlb = b[0].len().max(1);
    if na == 2
        && b.len() == 2
        && rla == 1
        && rlb == 1
        && a[1].first() == b[1].first()
        && a[1].len() == 1
        && b[1].len() == 1
    {
        map[1] = 1;
        return Ok(map);
    }
    let n = rla.min(rlb);
    let mut a: Vec<Vec<u8>> = a.to_vec();
    if a[0][..n.min(a[0].len())] != b[0][..n.min(b[0].len())] {
        if !a[0][..n.min(a[0].len())].eq_ignore_ascii_case(&b[0][..n.min(b[0].len())]) {
            return Err(format!(
                "The REF prefixes differ: {} vs {} ({rla},{rlb})",
                String::from_utf8_lossy(&a[0]),
                String::from_utf8_lossy(&b[0])
            ));
        }
        for x in a.iter_mut().chain(b.iter_mut()) {
            x.make_ascii_uppercase();
        }
    }
    if rla > rlb {
        let ext = a[0][rlb.min(a[0].len())..].to_vec();
        for x in b.iter_mut() {
            if x.first() == Some(&b'<') || x.first() == Some(&b'*') {
                continue;
            }
            x.extend_from_slice(&ext);
        }
    }
    for i in 1..na {
        let mut ai = a[i].clone();
        if rlb > rla && ai.first() != Some(&b'<') && ai.first() != Some(&b'*') {
            ai.extend_from_slice(&b[0][rla.min(b[0].len())..]);
        }
        match (1..b.len()).find(|&j| ai.eq_ignore_ascii_case(&b[j])) {
            Some(j) => map[i] = j,
            None => {
                map[i] = b.len();
                b.push(ai);
            }
        }
    }
    Ok(map)
}

/// `normalize_alleles`: drop the last base while every allele ends with it
/// and keeps at least one base.
pub(crate) fn trim_common_suffix(alleles: &mut [Vec<u8>]) {
    while alleles[0].len() > 1 {
        let last = alleles[0][alleles[0].len() - 1];
        if alleles
            .iter()
            .any(|a| a.len() < 2 || a[a.len() - 1] != last)
        {
            return;
        }
        for a in alleles.iter_mut() {
            a.pop();
        }
    }
}

struct Maux {
    als: Vec<Vec<u8>>,
    cnt: Vec<i32>,
    als_types: Vec<u32>,
    var_types: u32,
}

struct Stager<'a> {
    bufs: Vec<Buf<'a>>,
    collapse: u8,
    by_id: bool,
    maux: Maux,
    vcmp: Vcmp,
}

fn allele_type_mask(t: u32) -> u32 {
    if t != 0 { t << 1 } else { REF_MASK }
}

impl Stager<'_> {
    /// `maux_update_alleles`.
    fn update_alleles(&mut self, i: usize, j: usize) -> Result<(), String> {
        let rec = &mut self.bufs[i].recs[j];
        let n = rec.alleles.len();
        if self.maux.als.is_empty() {
            self.maux.als = rec.alleles.iter().map(|a| a.to_vec()).collect();
            rec.map = (0..n).collect();
            self.maux.cnt = vec![1; n];
            self.maux.als_types = (0..n).map(|k| allele_type_mask(rec.types[k])).collect();
            return Ok(());
        }
        let a: Vec<Vec<u8>> = rec.alleles.iter().map(|x| x.to_vec()).collect();
        rec.map = merge_alleles(&a, &mut self.maux.als)?;
        self.maux.cnt.resize(self.maux.als.len(), 0);
        self.maux.als_types.resize(self.maux.als.len(), 0);
        for k in 1..n {
            let ik = rec.map[k];
            self.maux.als_types[ik] = allele_type_mask(rec.types[k]);
            self.maux.cnt[ik] += 1;
        }
        self.maux.cnt[0] += 1;
        Ok(())
    }

    /// `types_compatible`.
    fn types_compatible(&mut self, selected: u32, i: usize, j: usize) -> bool {
        let rec = &self.bufs[i].recs[j];
        let rec_types = rec.ty;
        if self.collapse & COLLAPSE_ANY != 0 {
            return true;
        }
        if selected & REF_MASK != 0 && selected & !REF_MASK == 0 {
            return true;
        }
        if rec_types & REF_MASK != 0 && rec_types & !REF_MASK == 0 {
            return true;
        }
        if self.collapse != 0 {
            if self.collapse & COLLAPSE_SNPS != 0
                && rec_types & SNP_MASK != 0
                && selected & SNP_MASK != 0
            {
                return true;
            }
            if self.collapse & COLLAPSE_INDELS != 0
                && rec_types & INDEL_MASK != 0
                && selected & INDEL_MASK != 0
            {
                return true;
            }
        }
        let (mut x, mut y) = (selected >> 1, rec_types >> 1);
        while x != 0 && y != 0 {
            x >>= 1;
            y >>= 1;
        }
        if x != 0 || y != 0 {
            return false;
        }
        if self
            .vcmp
            .set_ref(&self.maux.als[0], rec.alleles[0])
            .is_none()
        {
            return false;
        }
        for k in 1..rec.alleles.len() {
            if rec.types[k] == 0 {
                continue;
            }
            if self
                .vcmp
                .find_allele(&self.maux.als[1..], rec.alleles[k])
                .is_some()
            {
                return true;
            }
        }
        false
    }

    /// `can_merge`: choose the candidates of this round; false when nothing is left.
    fn can_merge(&mut self) -> Result<bool, String> {
        self.maux.als.clear();
        self.maux.cnt.clear();
        self.maux.als_types.clear();
        self.maux.var_types = 0;
        let mut ntodo = 0;
        let mut id: Option<Vec<u8>> = None;
        for buf in &mut self.bufs {
            buf.var_types = 0;
            for rec in &mut buf.recs {
                if rec.skip & SKIP_DONE != 0 {
                    continue;
                }
                rec.skip = SKIP_DIFF;
                ntodo += 1;
                if self.by_id && id.is_none() {
                    id = Some(rec.pin.id.clone());
                    continue;
                }
                self.maux.var_types |= rec.ty;
                buf.var_types |= rec.ty;
            }
        }
        if ntodo == 0 {
            return Ok(false);
        }
        let mut selected = 0u32;
        for i in 0..self.bufs.len() {
            for j in 0..self.bufs[i].recs.len() {
                if self.bufs[i].recs[j].skip & SKIP_DONE != 0 {
                    continue;
                }
                let line_types = self.bufs[i].recs[j].ty;
                if self.by_id {
                    if id.as_deref() != Some(self.bufs[i].recs[j].pin.id.as_slice()) {
                        continue;
                    }
                } else if selected != 0 && !self.types_compatible(selected, i, j) {
                    continue;
                } else if (self.collapse & COLLAPSE_SNPS != 0 || self.collapse == 0)
                    && self.maux.var_types & SNP_MASK != 0
                    && line_types & (SNP_MASK | REF_MASK) == 0
                {
                    // SNPs first when there are any
                    continue;
                }
                selected |= line_types;
                self.bufs[i].recs[j].skip = 0;
                self.update_alleles(i, j)?;
            }
        }
        Ok(true)
    }

    /// `stage_line`: one record per input for the next merged record.
    fn stage_line(&mut self) {
        let m = &self.maux;
        let mut icnt: Option<usize> = None;
        for i in 1..m.als.len() {
            if m.als_types[i] & REF_MASK != 0 {
                continue;
            }
            if icnt.is_none_or(|c| m.cnt[c] < m.cnt[i]) {
                icnt = Some(i);
            }
        }
        let selected_type = icnt.map_or(REF_MASK, |c| m.als_types[c]);
        let icnt_v = icnt.unwrap_or(0);
        let all_var = m.var_types;
        let collapse = self.collapse;
        for buf in &mut self.bufs {
            buf.cur = None;
            let n = buf.recs.len();
            let mut found = None;
            for j in 0..n {
                let rec = &buf.recs[j];
                if rec.skip != 0 {
                    continue;
                }
                if self.by_id
                    || (!(selected_type & buf.var_types != 0 && selected_type & rec.ty == 0)
                        && ((selected_type & buf.var_types == 0 && rec.ty & REF_MASK != 0)
                            || selected_type == REF_MASK
                            || rec.map.iter().take(rec.alleles.len()).any(|&x| x == icnt_v)))
                {
                    found = Some(j);
                    break;
                }
            }
            if found.is_none() && collapse != 0 {
                // something compatible to make a multiallelic record of
                for (j, rec) in buf.recs.iter().enumerate().take(n) {
                    if rec.skip != 0 {
                        continue;
                    }
                    let lt = rec.ty;
                    let snp = collapse & COLLAPSE_SNPS != 0;
                    let indel = collapse & COLLAPSE_INDELS != 0;
                    let ok = collapse & COLLAPSE_ANY != 0
                        || (all_var & SNP_MASK != 0 && lt & SNP_MASK != 0 && snp)
                        || (all_var & INDEL_MASK != 0 && lt & INDEL_MASK != 0 && indel)
                        || if lt & REF_MASK != 0 {
                            (all_var & SNP_MASK != 0 && snp)
                                || (all_var & INDEL_MASK != 0 && indel)
                                || all_var & REF_MASK != 0
                        } else {
                            all_var & REF_MASK != 0
                                && ((lt & SNP_MASK != 0 && snp) || (lt & INDEL_MASK != 0 && indel))
                        };
                    if ok {
                        found = Some(j);
                        break;
                    }
                }
            }
            buf.cur = found;
        }
        for buf in &mut self.bufs {
            if let Some(c) = buf.cur {
                buf.recs[c].skip = SKIP_DONE;
            }
        }
    }

    /// The alleles of the staged records (`merge_chrom2qual`): only the merged
    /// alleles a staged record uses are written, REF is the longest seen.
    fn row(&self) -> StagedRow {
        let nals = self.maux.als.len();
        let mut used = vec![false; nals];
        for buf in &self.bufs {
            if let Some(c) = buf.cur {
                for &m in &buf.recs[c].map[1..] {
                    used[m] = true;
                }
            }
        }
        let mut newidx = vec![0usize; nals];
        let mut out: Vec<Vec<u8>> = vec![self.maux.als[0].clone()];
        for i in 1..nals {
            if used[i] {
                out.push(self.maux.als[i].clone());
                newidx[i] = out.len() - 1;
            }
        }
        trim_common_suffix(&mut out);
        let reference = out.remove(0);
        let mut members = Vec::with_capacity(self.bufs.len());
        let mut maps = Vec::new();
        let mut differs = Vec::new();
        for buf in &self.bufs {
            match buf.cur {
                Some(c) => {
                    let rec = &buf.recs[c];
                    members.push(Some(rec.orig));
                    let map: Vec<usize> = rec
                        .map
                        .iter()
                        .enumerate()
                        .map(|(k, &m)| if k == 0 { 0 } else { newidx[m] })
                        .collect();
                    differs.push(
                        map.len() != out.len() + 1 || map.iter().enumerate().any(|(i, &x)| i != x),
                    );
                    maps.push(map);
                }
                None => members.push(None),
            }
        }
        StagedRow {
            members,
            alleles: Alleles {
                reference,
                alts: out,
                maps,
                differs,
            },
        }
    }
}

/// The merged records of one position, in bcftools' order. `inputs[i]` are
/// the records of input `i` at the position, in file order.
pub fn stage_position(inputs: &[Vec<PairIn>], mode: MergeMode) -> Result<Vec<StagedRow>, String> {
    let n = inputs.len();
    // The synced reader's exact pairing decides the first record of each input.
    let pairs: Vec<Vec<_>> = inputs
        .iter()
        .map(|rs| {
            rs.iter()
                .map(|r| {
                    let alts: Vec<&[u8]> = r.alts.iter().map(Vec::as_slice).collect();
                    pair_rec(&r.reference, &alts, r.end)
                })
                .collect()
        })
        .collect();
    let active = inputs.iter().filter(|r| !r.is_empty()).count();
    let mut first: Vec<Option<usize>> = vec![None; n];
    if active == 1 {
        for (i, rs) in inputs.iter().enumerate() {
            if !rs.is_empty() {
                first[i] = Some(0);
            }
        }
    } else if active > 1 {
        let refs: Vec<&[_]> = pairs.iter().map(Vec::as_slice).collect();
        if let Some(row) = sort_rows(&refs, Pair::Exact).into_iter().next() {
            for (reader, rec) in row {
                first[reader] = Some(rec);
            }
        }
    }
    let (collapse, by_id) = match mode {
        MergeMode::None => (0, false),
        MergeMode::Snps => (COLLAPSE_SNPS, false),
        MergeMode::Indels => (COLLAPSE_INDELS, false),
        MergeMode::Both => (COLLAPSE_SNPS | COLLAPSE_INDELS, false),
        MergeMode::All => (COLLAPSE_ANY, false),
        MergeMode::Id => (0, true),
    };
    let mut bufs = Vec::with_capacity(n);
    for (i, rs) in inputs.iter().enumerate() {
        let mut order: Vec<usize> = (0..rs.len()).collect();
        if let Some(f) = first[i] {
            order.retain(|&k| k != f);
            order.insert(0, f);
        }
        let recs = order
            .into_iter()
            .map(|k| {
                let pin = &rs[k];
                let mut alleles: Vec<&[u8]> = vec![pin.reference.as_slice()];
                alleles.extend(pin.alts.iter().map(Vec::as_slice));
                let types: Vec<u32> = std::iter::once(0)
                    .chain(pin.alts.iter().map(|a| alt_type(&pin.reference, a)))
                    .collect();
                let all = types.iter().fold(0, |a, b| a | b);
                Rec {
                    pin,
                    orig: k,
                    alleles,
                    types,
                    ty: allele_type_mask(all),
                    skip: 0,
                    map: Vec::new(),
                }
            })
            .collect();
        bufs.push(Buf {
            recs,
            var_types: 0,
            cur: None,
        });
    }
    let mut st = Stager {
        bufs,
        collapse,
        by_id,
        maux: Maux {
            als: Vec::new(),
            cnt: Vec::new(),
            als_types: Vec::new(),
            var_types: 0,
        },
        vcmp: Vcmp::default(),
    };
    let mut rows = Vec::new();
    while st.can_merge()? {
        st.stage_line();
        rows.push(st.row());
    }
    Ok(rows)
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

    fn rows(files: &[Vec<PairIn>], m: MergeMode) -> Vec<Vec<Option<usize>>> {
        stage_position(files, m)
            .unwrap()
            .into_iter()
            .map(|r| r.members)
            .collect()
    }

    #[test]
    fn modes_decide_what_merges() {
        let files = vec![vec![rec("A", &["C"])], vec![rec("A", &["G"])]];
        assert_eq!(rows(&files, MergeMode::None).len(), 2);
        assert_eq!(rows(&files, MergeMode::Indels).len(), 2);
        for m in [MergeMode::Snps, MergeMode::Both, MergeMode::All] {
            assert_eq!(rows(&files, m), vec![vec![Some(0), Some(0)]]);
        }
        let mixed = vec![vec![rec("A", &["G"])], vec![rec("AT", &["A"])]];
        assert_eq!(rows(&mixed, MergeMode::Both).len(), 2);
        let all = stage_position(&mixed, MergeMode::All).unwrap();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0].alleles.reference, b"AT");
        assert_eq!(all[0].alleles.alts, vec![b"GT".to_vec(), b"A".to_vec()]);
    }

    #[test]
    fn shared_allele_merges_even_with_none() {
        // A>G merges with A>G,C (they share G) but not with A>T
        let files = vec![
            vec![rec("A", &["G"])],
            vec![rec("A", &["G", "C"])],
            vec![rec("A", &["T"])],
        ];
        let r = rows(&files, MergeMode::None);
        assert_eq!(
            r,
            vec![vec![Some(0), Some(0), None], vec![None, None, Some(0)]]
        );
    }

    #[test]
    fn snps_go_before_other_records() {
        // <DEL> is staged after the SNP (a round starts over, so its REF stays A)
        let files = vec![vec![rec("A", &["<DEL>"]), rec("AT", &["GC"])]];
        let r = stage_position(&files, MergeMode::None).unwrap();
        assert_eq!(r[0].members, vec![Some(1)]);
        assert_eq!(r[1].members, vec![Some(0)]);
        // each round starts over: the <DEL> row is not extended
        assert_eq!(r[1].alleles.reference, b"A");
    }

    #[test]
    fn id_mode_pairs_equal_ids_only() {
        let mut a = rec("A", &["C"]);
        a.id = b"x".to_vec();
        let mut b = rec("A", &["G"]);
        b.id = b"x".to_vec();
        let mut c = rec("A", &["T"]);
        c.id = b"y".to_vec();
        let files = vec![vec![a], vec![b], vec![c]];
        let r = rows(&files, MergeMode::Id);
        assert_eq!(
            r,
            vec![vec![Some(0), Some(0), None], vec![None, None, Some(0)]]
        );
    }

    #[test]
    fn two_hundred_files_at_one_position() {
        let files: Vec<Vec<PairIn>> = (0..200)
            .map(|i| vec![rec("A", &[["C", "G", "T"][i % 3]])])
            .collect();
        let r = stage_position(&files, MergeMode::Both).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].alleles.alts.len(), 3);
    }
}
