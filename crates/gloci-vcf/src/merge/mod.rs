//! `gloci vcf merge`: the record and header rules of `bcftools merge`
//! (1.21): which records merge is in [`pair`]; this module merges the
//! records of one row into one. The streaming read of the inputs is in
//! gloci-ops.

mod format;
pub mod header;
mod info;
pub mod pair;

use crate::annotate::hdr::{InfoDef, Number, Ty};
use header::{Defs, Merged};
pub use pair::{MergeMode, PairIn, rows};

/// Options of one merge.
#[derive(Clone, Debug)]
pub struct MergeOptions {
    pub mode: MergeMode,
    /// A sample absent from a record gets `0/0` instead of `./.`.
    pub missing_to_ref: bool,
    /// `--info-rules`: `None` is bcftools' default, `-` none at all.
    pub info_rules: Option<String>,
    /// `-F x`: a PASS in any input makes the merged FILTER PASS.
    pub filter_x: bool,
}

/// The columns of one input record.
pub(crate) struct P<'a> {
    pub reader: usize,
    pub chrom: &'a [u8],
    pub pos: &'a [u8],
    pub id: &'a [u8],
    pub reference: &'a [u8],
    pub alt: &'a [u8],
    pub qual: &'a [u8],
    pub filter: &'a [u8],
    pub info: &'a [u8],
    pub fmt: Vec<&'a [u8]>,
    pub samples: Vec<&'a [u8]>,
}

impl<'a> P<'a> {
    fn parse(reader: usize, line: &'a [u8]) -> Result<P<'a>, String> {
        let cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        if cols.len() < 8 {
            return Err(format!(
                "the record has {} tab-separated columns; a VCF record needs at least 8",
                cols.len()
            ));
        }
        let fmt = if cols.len() > 8 && cols[8] != b"." {
            cols[8].split(|&b| b == b':').collect()
        } else {
            Vec::new()
        };
        Ok(P {
            reader,
            chrom: cols[0],
            pos: cols[1],
            id: cols[2],
            reference: cols[3],
            alt: cols[4],
            qual: cols[5],
            filter: cols[6],
            info: cols[7],
            fmt,
            samples: cols.get(9..).unwrap_or(&[]).to_vec(),
        })
    }

    pub(crate) fn alts(&self) -> Vec<&'a [u8]> {
        if self.alt == b"." || self.alt.is_empty() {
            Vec::new()
        } else {
            self.alt.split(|&b| b == b',').collect()
        }
    }
}

/// The merged allele list and each input's map into it.
pub(crate) struct Alleles {
    pub reference: Vec<u8>,
    pub alts: Vec<Vec<u8>>,
    /// Per row member (same order as the records): input allele index to merged index.
    pub maps: Vec<Vec<usize>>,
    /// Per row member: its alleles are not the merged list as it is.
    pub differs: Vec<bool>,
}

fn symbolic(a: &[u8]) -> bool {
    a.first() == Some(&b'<') || a == b"*" || a.contains(&b'[') || a.contains(&b']')
}

/// `normalize_alleles`: drop the last base while every allele ends with it
/// and keeps at least one base (a symbolic ALT, `*` or a one-base REF stops it).
fn trim_common_suffix(alleles: &mut [Vec<u8>]) {
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

/// The record's alleles, REF first.
fn record_alleles(r: &P<'_>) -> Vec<Vec<u8>> {
    let mut v = vec![r.reference.to_vec()];
    v.extend(r.alts().iter().map(|a| a.to_vec()));
    v
}

fn merge_alleles(recs: &[P<'_>]) -> Result<Alleles, String> {
    let all: Vec<Vec<Vec<u8>>> = recs.iter().map(record_alleles).collect();
    let mut reference: &[u8] = &all[0][0];
    for al in &all {
        if al[0].len() > reference.len() {
            reference = &al[0];
        }
    }
    let reference = reference.to_vec();
    let mut alts: Vec<Vec<u8>> = Vec::new();
    let mut maps = Vec::new();
    let mut ori_counts = Vec::new();
    for (r, al) in recs.iter().zip(&all) {
        let rl = al[0].len();
        if !reference[..rl].eq_ignore_ascii_case(&al[0]) {
            return Err(format!(
                "the REF alleles at {}:{} differ ({} and {}); records with different REF cannot be merged",
                String::from_utf8_lossy(r.chrom),
                String::from_utf8_lossy(r.pos),
                String::from_utf8_lossy(&reference),
                String::from_utf8_lossy(&al[0])
            ));
        }
        let suffix = &reference[rl..];
        let mut map = vec![0usize];
        for a in &al[1..] {
            let mut norm = a.clone();
            if !suffix.is_empty() && !symbolic(a) {
                norm.extend_from_slice(suffix);
            }
            let idx = match alts.iter().position(|x| *x == norm) {
                Some(i) => i,
                None => {
                    alts.push(norm);
                    alts.len() - 1
                }
            };
            map.push(idx + 1);
        }
        ori_counts.push(al.len());
        maps.push(map);
    }
    let mut reference = reference;
    {
        let mut v = vec![reference.clone()];
        v.extend(alts.iter().cloned());
        trim_common_suffix(&mut v);
        reference = v.remove(0);
        alts = v;
    }
    let differs = maps
        .iter()
        .zip(&ori_counts)
        .map(|(m, &n)| n != alts.len() + 1 || m.iter().enumerate().any(|(i, &x)| i != x))
        .collect();
    Ok(Alleles {
        reference,
        alts,
        maps,
        differs,
    })
}

/// What a merge needs of one input.
struct Input {
    name: String,
    samples: usize,
}

/// A configured merge of the inputs whose headers were merged.
pub struct Merger {
    opts: MergeOptions,
    defs: Defs,
    own: Vec<Defs>,
    inputs: Vec<Input>,
    total_samples: usize,
    rules: Vec<info::Rule>,
    /// AC and AN are removed and recomputed from the merged genotypes.
    calc_acan: bool,
    scratch: Vec<u8>,
}

fn pos_label(r: &P<'_>) -> String {
    format!(
        "{}:{}",
        String::from_utf8_lossy(r.chrom),
        String::from_utf8_lossy(r.pos)
    )
}

impl Merger {
    pub fn new(opts: MergeOptions, merged: &Merged, names: &[String]) -> Result<Merger, String> {
        let (rules, calc_acan) = info::build_rules(&opts, merged)?;
        Ok(Merger {
            opts,
            defs: merged.defs.clone(),
            own: merged.own.clone(),
            inputs: names
                .iter()
                .zip(&merged.counts)
                .map(|(n, &c)| Input {
                    name: n.clone(),
                    samples: c,
                })
                .collect(),
            total_samples: merged.samples.len(),
            rules,
            calc_acan,
            scratch: Vec::new(),
        })
    }

    /// What pairing needs of `line`.
    pub fn pair_in(&self, line: &[u8]) -> Result<PairIn, String> {
        let p = P::parse(0, line)?;
        Ok(PairIn {
            reference: p.reference.to_vec(),
            alts: p.alts().iter().map(|a| a.to_vec()).collect(),
            id: p.id.to_vec(),
            end: info::end_value(p.info),
        })
    }

    pub fn options(&self) -> &MergeOptions {
        &self.opts
    }

    /// Merge the records of one row (`row[i]` is input i's record line, no
    /// terminator) and append the merged line to `out`.
    pub fn merge_row(&mut self, row: &[Option<&[u8]>], out: &mut Vec<u8>) -> Result<(), String> {
        let mut recs = Vec::new();
        for (i, l) in row.iter().enumerate() {
            if let Some(l) = l {
                let p = P::parse(i, l).map_err(|e| format!("{}: {e}", self.inputs[i].name))?;
                if p.samples.len() != self.inputs[i].samples {
                    return Err(format!(
                        "{}: the record at {} has {} sample columns but the header has {}",
                        self.inputs[i].name,
                        pos_label(&p),
                        p.samples.len(),
                        self.inputs[i].samples
                    ));
                }
                recs.push(p);
            }
        }
        if recs.is_empty() {
            return Ok(());
        }
        let al = merge_alleles(&recs)?;
        let first = &recs[0];
        out.extend_from_slice(first.chrom);
        out.push(b'\t');
        out.extend_from_slice(first.pos);
        out.push(b'\t');
        merge_ids(&recs, out);
        out.push(b'\t');
        out.extend_from_slice(&al.reference);
        out.push(b'\t');
        if al.alts.is_empty() {
            out.push(b'.');
        } else {
            for (i, a) in al.alts.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(a);
            }
        }
        out.push(b'\t');
        merge_qual(&recs, out);
        out.push(b'\t');
        merge_filter(&recs, out);
        out.push(b'\t');
        let mut fmt = format::Formatted::default();
        if self.total_samples > 0 {
            fmt = format::merge_format(self, &recs, &al)?;
        }
        let mut info_out = std::mem::take(&mut self.scratch);
        info_out.clear();
        info::merge_info(self, &recs, &al, &fmt.calls, &mut info_out)?;
        if info_out.is_empty() {
            out.push(b'.');
        } else {
            out.extend_from_slice(&info_out);
        }
        self.scratch = info_out;
        if self.total_samples > 0 {
            out.extend_from_slice(&fmt.text);
        }
        Ok(())
    }

    pub(crate) fn out_def(&self, kind_info: bool, key: &[u8]) -> InfoDef {
        let m = if kind_info {
            &self.defs.info
        } else {
            &self.defs.format
        };
        m.get(key).map_or(
            InfoDef {
                ty: Ty::Str,
                number: Number::Var,
            },
            |d| d.1,
        )
    }

    /// An error when input `reader` declares `key` with another type than the merged header.
    pub(crate) fn check_type(
        &self,
        kind_info: bool,
        reader: usize,
        key: &[u8],
        at: &str,
    ) -> Result<(), String> {
        let own = if kind_info {
            self.own[reader].info.get(key)
        } else {
            self.own[reader].format.get(key)
        };
        let Some((_, d)) = own else { return Ok(()) };
        let out = self.out_def(kind_info, key);
        if d.ty != out.ty {
            return Err(format!(
                "{}: {}/{} at {at} is declared with a different type than in the first input that declares it; make the ##{} lines agree (for example with bcftools annotate -h) before merging",
                self.inputs[reader].name,
                if kind_info { "INFO" } else { "FORMAT" },
                String::from_utf8_lossy(key),
                if kind_info { "INFO" } else { "FORMAT" },
            ));
        }
        Ok(())
    }
}

fn merge_ids(recs: &[P<'_>], out: &mut Vec<u8>) {
    let mut seen: Vec<&[u8]> = Vec::new();
    for r in recs {
        if r.id == b"." || r.id.is_empty() || seen.contains(&r.id) {
            continue;
        }
        seen.push(r.id);
    }
    if seen.is_empty() {
        out.push(b'.');
    }
    for (i, s) in seen.iter().enumerate() {
        if i > 0 {
            out.push(b';');
        }
        out.extend_from_slice(s);
    }
}

pub(crate) fn parse_f32(tok: &[u8]) -> Option<f32> {
    std::str::from_utf8(tok)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map(|v| v as f32)
}

pub(crate) fn put_f32(out: &mut Vec<u8>, v: f32) {
    gloci_core::cfmt::write_g6(out, f64::from(v));
}

fn merge_qual(recs: &[P<'_>], out: &mut Vec<u8>) {
    let mut best: Option<f32> = None;
    for r in recs {
        if let Some(v) = parse_f32(r.qual)
            && best.is_none_or(|b| v > b)
        {
            best = Some(v);
        }
    }
    match best {
        Some(v) => put_f32(out, v),
        None => out.push(b'.'),
    }
}

fn merge_filter(recs: &[P<'_>], out: &mut Vec<u8>) {
    let mut names: Vec<&[u8]> = Vec::new();
    let mut pass = false;
    for r in recs {
        if r.filter == b"." || r.filter.is_empty() {
            continue;
        }
        for f in r.filter.split(|&b| b == b';') {
            if f == b"PASS" {
                pass = true;
            } else if !names.contains(&f) {
                names.push(f);
            }
        }
    }
    if names.is_empty() {
        out.extend_from_slice(if pass { b"PASS" } else { b"." });
        return;
    }
    for (i, f) in names.iter().enumerate() {
        if i > 0 {
            out.push(b';');
        }
        out.extend_from_slice(f);
    }
}
