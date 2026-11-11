//! The FORMAT columns of a merged record: the union of the inputs' keys,
//! each sample's values re-indexed to the merged allele list.

use super::{Alleles, Merger, P, parse_f32, put_f32};
use crate::annotate::hdr::{Number, Ty};

/// The called alleles of the merged genotypes: AN and AC per ALT.
#[derive(Clone, Debug)]
pub struct Calls {
    pub an: u64,
    pub ac: Vec<u64>,
}

#[derive(Default)]
pub struct Formatted {
    /// FORMAT and the samples, with a leading tab.
    pub text: Vec<u8>,
    pub calls: Option<Calls>,
}

/// Write one value token as htslib re-serialises it by type.
pub(crate) fn put_val(out: &mut Vec<u8>, ty: Ty, tok: &[u8]) {
    match ty {
        Ty::Int => match std::str::from_utf8(tok)
            .ok()
            .and_then(|s| s.parse::<i64>().ok())
        {
            Some(v) => out.extend_from_slice(itoa::Buffer::new().format(v).as_bytes()),
            None => out.extend_from_slice(if tok.is_empty() { b"." } else { tok }),
        },
        Ty::Float => match parse_f32(tok) {
            Some(v) if v.is_finite() => put_f32(out, v),
            _ => out.extend_from_slice(if tok.is_empty() { b"." } else { tok }),
        },
        _ => out.extend_from_slice(if tok.is_empty() { b"." } else { tok }),
    }
}

fn put_list(out: &mut Vec<u8>, ty: Ty, toks: &[&[u8]]) {
    for (i, t) in toks.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        put_val(out, ty, t);
    }
}

/// The alleles of a GT cell and the separator before each (`/` or `|`; the first has none).
fn gt_parts(cell: &[u8]) -> (Vec<&[u8]>, Vec<u8>) {
    let mut alleles = Vec::new();
    let mut seps = Vec::new();
    let mut start = 0;
    for (i, &b) in cell.iter().enumerate() {
        if b == b'/' || b == b'|' {
            alleles.push(&cell[start..i]);
            seps.push(b);
            start = i + 1;
        }
    }
    alleles.push(&cell[start..]);
    (alleles, seps)
}

#[allow(clippy::needless_range_loop)]
pub(crate) fn merge_format(m: &Merger, recs: &[P<'_>], al: &Alleles) -> Result<Formatted, String> {
    let nin = m.inputs.len();
    // Row member of each input.
    let mut member: Vec<Option<usize>> = vec![None; nin];
    for (k, r) in recs.iter().enumerate() {
        member[r.reader] = Some(k);
    }
    let mut keys: Vec<&[u8]> = Vec::new();
    for r in recs {
        for k in &r.fmt {
            // A key the merged header does not declare is dropped, as bcftools does.
            if !keys.contains(k) && m.defs.format.contains_key(*k) {
                keys.push(k);
            }
        }
    }
    if let Some(p) = keys.iter().position(|k| *k == b"GT") {
        let gt = keys.remove(p);
        keys.insert(0, gt);
    }
    let total = m.total_samples;
    let mut out = Vec::new();
    if keys.is_empty() {
        out.extend_from_slice(b"\t.");
        for _ in 0..total {
            out.extend_from_slice(b"\t.");
        }
        return Ok(Formatted {
            text: out,
            calls: None,
        });
    }
    let mut cells: Vec<Vec<Vec<u8>>> = vec![Vec::with_capacity(keys.len()); total];
    let mut calls = None;
    for key in &keys {
        // The index of the key in each input's FORMAT.
        let idx: Vec<Option<usize>> = (0..nin)
            .map(|i| member[i].and_then(|k| recs[k].fmt.iter().position(|x| x == key)))
            .collect();
        let field = |i: usize, s: usize| -> Option<&[u8]> {
            let k = member[i]?;
            let fi = idx[i]?;
            recs[k].samples[s].split(|&b| b == b':').nth(fi)
        };
        let mut col: Vec<Vec<u8>> = Vec::with_capacity(total);
        if *key == b"GT" {
            let mut nsize = 1;
            for (i, input) in m.inputs.iter().enumerate() {
                for s in 0..input.samples {
                    if let Some(f) = field(i, s) {
                        nsize = nsize.max(gt_parts(f).0.len());
                    }
                }
            }
            let mut an = 0u64;
            let mut ac = vec![0u64; al.alts.len()];
            for i in 0..nin {
                for s in 0..m.inputs[i].samples {
                    let mut cell = Vec::new();
                    match field(i, s) {
                        Some(f) => {
                            let k = member[i].unwrap_or(0);
                            let map = &al.maps[k];
                            let (alleles, seps) = gt_parts(f);
                            for (j, a) in alleles.iter().enumerate() {
                                if j > 0 {
                                    cell.push(seps[j - 1]);
                                }
                                match std::str::from_utf8(a)
                                    .ok()
                                    .and_then(|x| x.parse::<usize>().ok())
                                {
                                    Some(v) => {
                                        let nv = *map.get(v).ok_or_else(|| {
                                            format!(
                                                "{}: the genotype {} at {} uses allele {v} but the record has only {} alleles",
                                                m.inputs[i].name,
                                                String::from_utf8_lossy(f),
                                                super::pos_label(&recs[k]),
                                                map.len()
                                            )
                                        })?;
                                        cell.extend_from_slice(
                                            itoa::Buffer::new().format(nv).as_bytes(),
                                        );
                                        an += 1;
                                        if nv > 0 {
                                            ac[nv - 1] += 1;
                                        }
                                    }
                                    None => cell.push(b'.'),
                                }
                            }
                        }
                        None => {
                            let fill = if m.opts.missing_to_ref { b'0' } else { b'.' };
                            for j in 0..nsize {
                                if j > 0 {
                                    cell.push(b'/');
                                }
                                cell.push(fill);
                            }
                            if m.opts.missing_to_ref {
                                an += nsize as u64;
                            }
                        }
                    }
                    col.push(cell);
                }
            }
            calls = Some(Calls { an, ac });
        } else {
            let def = m.out_def(false, key);
            let n_alt = al.alts.len();
            for i in 0..nin {
                if let Some(k) = member[i]
                    && idx[i].is_some()
                {
                    m.check_type(false, i, key, &super::pos_label(&recs[k]))?;
                }
                for s in 0..m.inputs[i].samples {
                    let mut cell = Vec::new();
                    match field(i, s) {
                        None => cell.push(b'.'),
                        Some(f) => {
                            let k = member[i].unwrap_or(0);
                            let toks: Vec<&[u8]> = f.split(|&b| b == b',').collect();
                            let remap = al.differs[k]
                                && matches!(def.number, Number::A | Number::R | Number::G);
                            if !remap {
                                put_list(&mut cell, def.ty, &toks);
                            } else {
                                let map = &al.maps[k];
                                let n_ori = map.len();
                                let (len, slots): (usize, Vec<Option<usize>>) = match def.number {
                                    Number::A => (
                                        n_alt,
                                        toks.iter()
                                            .enumerate()
                                            .map(|(j, _)| (j < n_ori - 1).then(|| map[j + 1] - 1))
                                            .collect(),
                                    ),
                                    Number::R => (
                                        n_alt + 1,
                                        toks.iter()
                                            .enumerate()
                                            .map(|(j, _)| (j < n_ori).then(|| map[j]))
                                            .collect(),
                                    ),
                                    _ => {
                                        let nd = n_ori * (n_ori + 1) / 2;
                                        if toks.len() == nd {
                                            let mut v = Vec::new();
                                            for b in 0..n_ori {
                                                for a in 0..=b {
                                                    let (x, y) = (map[a], map[b]);
                                                    let (lo, hi) =
                                                        if x <= y { (x, y) } else { (y, x) };
                                                    v.push(Some(hi * (hi + 1) / 2 + lo));
                                                }
                                            }
                                            ((n_alt + 1) * (n_alt + 2) / 2, v)
                                        } else if toks.len() == n_ori {
                                            (n_alt + 1, (0..n_ori).map(|j| Some(map[j])).collect())
                                        } else {
                                            return Err(format!(
                                                "{}: FORMAT/{} at {} has {} values; its Number=G needs {} (diploid) or {} (haploid) for {} alleles",
                                                m.inputs[i].name,
                                                String::from_utf8_lossy(key),
                                                super::pos_label(&recs[k]),
                                                toks.len(),
                                                nd,
                                                n_ori,
                                                n_ori
                                            ));
                                        }
                                    }
                                };
                                let mut v: Vec<Option<&[u8]>> = vec![None; len];
                                for (t, sl) in toks.iter().zip(slots) {
                                    if let Some(sl) = sl {
                                        v[sl] = Some(t);
                                    }
                                }
                                for (j, x) in v.iter().enumerate() {
                                    if j > 0 {
                                        cell.push(b',');
                                    }
                                    match x {
                                        Some(t) => put_val(&mut cell, def.ty, t),
                                        None => cell.push(b'.'),
                                    }
                                }
                            }
                        }
                    }
                    col.push(cell);
                }
            }
        }
        for (s, c) in col.into_iter().enumerate() {
            cells[s].push(c);
        }
    }
    out.push(b'\t');
    for (i, k) in keys.iter().enumerate() {
        if i > 0 {
            out.push(b':');
        }
        out.extend_from_slice(k);
    }
    for sample in &cells {
        out.push(b'\t');
        for (i, c) in sample.iter().enumerate() {
            if i > 0 {
                out.push(b':');
            }
            out.extend_from_slice(c);
        }
    }
    Ok(Formatted { text: out, calls })
}
