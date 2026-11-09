//! Evaluation of a parsed expression on one record.

use super::lexer::Op;
use super::parser::{HAP_ALT, HAP_REF, HET, HOM_ALT, HOM_REF, MIS};
use super::{EvalScratch, Field, InfoTy, Node, Rhs, TypeKw};
use crate::filter::SiteView;
use crate::gt::{Genotype, GenotypeSource, SampleFields};
use crate::stats::SiteStats;
use gloci_core::{VariantKind, alt_kind};

/// Genotypes with the samples the mask rules masked read as missing (D50):
/// no GT, DP or GQ. `masked` is the vector `site_stats` filled.
pub struct MaskedGenotypes<'a> {
    pub inner: &'a dyn GenotypeSource,
    pub masked: &'a [bool],
}

impl GenotypeSource for MaskedGenotypes<'_> {
    fn n_samples(&self) -> usize {
        self.inner.n_samples()
    }

    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        self.inner.sample(i, out)?;
        if self.masked.get(i).copied().unwrap_or(false) {
            out.clear_gt();
            out.dp = None;
            out.gq = None;
        }
        Ok(())
    }
}

/// One value of a field.
#[derive(Clone, Copy)]
enum Val<'a> {
    Num(f64),
    Str(&'a [u8]),
    Miss,
}

struct Cx<'a, 'b> {
    site: &'a SiteView<'a>,
    gt: Option<&'a dyn GenotypeSource>,
    stats: Option<&'a SiteStats>,
    scratch: &'b mut EvalScratch,
    /// Inside `N_PASS`: the sample whose fields are loaded in `scratch`.
    sample: Option<usize>,
}

pub(super) fn eval(
    node: &Node,
    site: &SiteView,
    gt: Option<&dyn GenotypeSource>,
    stats: Option<&SiteStats>,
    scratch: &mut EvalScratch,
) -> Result<bool, String> {
    let mut cx = Cx {
        site,
        gt,
        stats,
        scratch,
        sample: None,
    };
    node_holds(node, &mut cx)
}

fn node_holds(node: &Node, cx: &mut Cx) -> Result<bool, String> {
    match node {
        Node::Or(v) => {
            for n in v {
                if node_holds(n, cx)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        Node::And(v) => {
            for n in v {
                if !node_holds(n, cx)? {
                    return Ok(false);
                }
            }
            Ok(true)
        }
        Node::Not(n) => Ok(!node_holds(n, cx)?),
        Node::Truthy(f) => {
            let mut any = false;
            visit(f, cx, &mut |v| {
                any = matches!(v, Val::Num(n) if n != 0.0);
                true
            })?;
            Ok(any)
        }
        Node::Cmp { field, op, rhs } => compare(field, *op, rhs, cx),
    }
}

fn compare(field: &Field, op: Op, rhs: &Rhs, cx: &mut Cx) -> Result<bool, String> {
    match (field, rhs) {
        (Field::Type, Rhs::TypeKw(k)) => {
            let m = type_matches(cx.site, *k);
            Ok(m == (op == Op::Eq))
        }
        (Field::Gt, Rhs::GtKw(mask)) => {
            let want_in = op == Op::Eq;
            let mut hit = false;
            each_sample(cx, &mut |_, f| {
                let class = f.gt.as_ref().map_or(MIS, gt_class);
                hit = (class & mask != 0) == want_in;
                hit
            })?;
            Ok(hit)
        }
        (Field::Gt, _) => {
            let mut hit = false;
            each_sample(cx, &mut |buf, f| {
                render_gt(f.gt.as_ref(), buf);
                hit = test(Val::Str(buf), op, rhs);
                hit
            })?;
            Ok(hit)
        }
        _ => {
            let mut hit = false;
            visit(field, cx, &mut |v| {
                hit = test(v, op, rhs);
                hit
            })?;
            Ok(hit)
        }
    }
}

/// One value against the literal. A missing value fails every comparison
/// except `!=` and `!~` (bcftools).
fn test(v: Val, op: Op, rhs: &Rhs) -> bool {
    match (v, rhs) {
        (Val::Num(x), Rhs::Num(y)) => match op {
            Op::Eq => x == *y,
            Op::Ne => x != *y,
            Op::Lt => x < *y,
            Op::Le => x <= *y,
            Op::Gt => x > *y,
            Op::Ge => x >= *y,
            _ => false,
        },
        (Val::Miss, Rhs::Num(_)) => op == Op::Ne,
        (Val::Miss, Rhs::Miss) => op == Op::Eq,
        (Val::Num(_), Rhs::Miss) => op == Op::Ne,
        (Val::Str(s), Rhs::Str(t)) => match op {
            Op::Eq => s == &**t,
            _ => s != &**t,
        },
        (Val::Miss, Rhs::Str(_)) => op == Op::Ne,
        (Val::Str(s), Rhs::Re(re)) => re.is_match(s) == (op == Op::Match),
        (Val::Miss, Rhs::Re(_)) => op == Op::NoMatch,
        _ => false,
    }
}

fn type_matches(site: &SiteView, k: TypeKw) -> bool {
    // True when every ALT is of kind `k` (bcftools compares the record's
    // whole type set); `*` has no TYPE.
    let mut n = 0;
    for alt in site.alts.split(|&b| b == b',') {
        let kind = if alt == b"*" {
            return false;
        } else {
            alt_kind(site.pos.saturating_sub(1), site.ref_, alt)
        };
        let alt_k = match kind {
            VariantKind::Snv => TypeKw::Snp,
            VariantKind::Mnv => TypeKw::Mnp,
            VariantKind::Ins | VariantKind::Del => TypeKw::Indel,
            VariantKind::Sv | VariantKind::Other => TypeKw::Other,
            VariantKind::Ref => TypeKw::Ref,
        };
        if alt_k != k {
            return false;
        }
        n += 1;
    }
    n > 0
}

pub(crate) fn gt_class(g: &Genotype) -> u8 {
    let a = &g.alleles;
    if a.is_empty() || a.iter().any(Option::is_none) {
        return MIS;
    }
    let first = a[0];
    match (a.len(), first) {
        (1, Some(0)) => HAP_REF,
        (1, _) => HAP_ALT,
        _ if a.iter().all(|x| *x == first) => {
            if first == Some(0) {
                HOM_REF
            } else {
                HOM_ALT
            }
        }
        _ => HET,
    }
}

/// `0/1`, `1|0`, `0`, `./.`; a mixed-separator genotype renders unphased.
fn render_gt(g: Option<&Genotype>, out: &mut Vec<u8>) {
    out.clear();
    let Some(g) = g else {
        out.push(b'.');
        return;
    };
    let sep = if g.phased { b'|' } else { b'/' };
    for (i, a) in g.alleles.iter().enumerate() {
        if i > 0 {
            out.push(sep);
        }
        match a {
            Some(n) => out.extend_from_slice(itoa::Buffer::new().format(*n).as_bytes()),
            None => out.push(b'.'),
        }
    }
}

fn need_gt<'a>(cx: &Cx<'a, '_>) -> Result<&'a dyn GenotypeSource, String> {
    cx.gt.ok_or_else(|| {
        "the expression reads genotypes (FMT/, GT or N_PASS) but none were supplied; \
         this is a bug, please report it"
            .to_string()
    })
}

fn need_stats<'a>(cx: &Cx<'a, '_>, what: &str) -> Result<&'a SiteStats, String> {
    cx.stats.ok_or_else(|| {
        format!(
            "{what} needs genotype statistics but none were supplied; \
             this is a bug, please report it"
        )
    })
}

/// Call `f` with each considered sample's fields (just the current one inside
/// `N_PASS`), stopping when it returns true. The first argument is a scratch
/// text buffer.
fn each_sample(
    cx: &mut Cx,
    f: &mut dyn FnMut(&mut Vec<u8>, &SampleFields) -> bool,
) -> Result<bool, String> {
    if cx.sample.is_some() {
        let s = &mut *cx.scratch;
        return Ok(f(&mut s.gt_text, &s.fields));
    }
    let src = need_gt(cx)?;
    for i in 0..src.n_samples() {
        if cx
            .scratch
            .keep
            .as_ref()
            .is_some_and(|k| !k.get(i).copied().unwrap_or(true))
        {
            continue;
        }
        src.sample(i, &mut cx.scratch.fields)?;
        let s = &mut *cx.scratch;
        if f(&mut s.gt_text, &s.fields) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Call `cb` with each value of `field`, stopping when it returns true.
fn visit(field: &Field, cx: &mut Cx, cb: &mut dyn FnMut(Val) -> bool) -> Result<bool, String> {
    let site = cx.site;
    Ok(match field {
        Field::Qual => cb(site.qual.map_or(Val::Miss, Val::Num)),
        Field::Pos => cb(Val::Num(site.pos as f64)),
        Field::Chrom => cb(Val::Str(site.chrom)),
        Field::Ref => cb(Val::Str(site.ref_)),
        Field::Id => site.id.split(|&b| b == b';').any(|s| cb(Val::Str(s))),
        Field::Alt => site.alts.split(|&b| b == b',').any(|s| cb(Val::Str(s))),
        Field::Filter => cb(Val::Str(if site.filter.is_empty() {
            b"."
        } else {
            site.filter
        })),
        Field::NAlt => {
            let n = if site.alts.is_empty() || site.alts == b"." {
                0
            } else {
                site.alts.iter().filter(|&&b| b == b',').count() + 1
            };
            cb(Val::Num(n as f64))
        }
        Field::NSamples(n) => {
            let n = match (cx.stats, cx.gt) {
                (Some(s), _) => s.samples as usize,
                (None, Some(g)) => g.n_samples(),
                _ => *n,
            };
            cb(Val::Num(n as f64))
        }
        Field::Info { key, ty, idx } => visit_info(site.info, key, *ty, *idx, cb)?,
        Field::FmtDp(idx) | Field::FmtGq(idx) => {
            let dp = matches!(field, Field::FmtDp(_));
            let mut stop = false;
            each_sample(cx, &mut |_, f| {
                let v = if dp { f.dp } else { f.gq };
                let val = match (v, idx) {
                    (Some(n), None | Some(0)) => Val::Num(f64::from(n)),
                    _ => Val::Miss,
                };
                stop = cb(val);
                stop
            })?;
            stop
        }
        Field::Gt => unreachable!("GT is compared through compare()"),
        Field::Type => unreachable!("TYPE is compared through compare()"),
        Field::Maf => cb(need_stats(cx, "MAF")?
            .expr_maf()
            .map_or(Val::Miss, Val::Num)),
        Field::An => cb(Val::Num(need_stats(cx, "AN")?.an() as f64)),
        Field::FMissing => cb(Val::Num(need_stats(cx, "F_MISSING")?.missing_fraction())),
        Field::NMissing => cb(Val::Num(need_stats(cx, "N_MISSING")?.missing as f64)),
        Field::Ac(idx) => {
            // Unindexed, the computed AC is the first ALT's count only (as in
            // bcftools), missing when no allele is called; AC[i] picks
            // another. INFO/AC matches any value.
            let st = need_stats(cx, "AC")?;
            let ac = if st.an() == 0 { &[][..] } else { st.ac() };
            cb(ac
                .get(idx.unwrap_or(0))
                .map_or(Val::Miss, |&c| Val::Num(c as f64)))
        }
        Field::NPass(inner) => {
            let src = need_gt(cx)?;
            let mut count = 0usize;
            for i in 0..src.n_samples() {
                if cx
                    .scratch
                    .keep
                    .as_ref()
                    .is_some_and(|k| !k.get(i).copied().unwrap_or(true))
                {
                    continue;
                }
                src.sample(i, &mut cx.scratch.fields)?;
                cx.sample = Some(i);
                let r = node_holds(inner, cx);
                cx.sample = None;
                if r? {
                    count += 1;
                }
            }
            cb(Val::Num(count as f64))
        }
    })
}

fn visit_info(
    info: &[u8],
    key: &[u8],
    ty: InfoTy,
    idx: Option<usize>,
    cb: &mut dyn FnMut(Val) -> bool,
) -> Result<bool, String> {
    let item = info
        .split(|&b| b == b';')
        .find(|item| item.split(|&b| b == b'=').next() == Some(key));
    if ty == InfoTy::Flag {
        return Ok(cb(Val::Num(f64::from(u8::from(item.is_some())))));
    }
    let value = item.and_then(|i| i.iter().position(|&b| b == b'=').map(|p| &i[p + 1..]));
    let Some(value) = value else {
        return Ok(cb(Val::Miss));
    };
    let one = |v: &[u8], cb: &mut dyn FnMut(Val) -> bool| -> Result<bool, String> {
        Ok(match ty {
            InfoTy::Str => cb(Val::Str(v)),
            _ if v == b"." => cb(Val::Miss),
            _ => {
                let n = std::str::from_utf8(v)
                    .ok()
                    .and_then(|s| s.parse::<f64>().ok())
                    .ok_or_else(|| {
                        format!(
                            "INFO/{} value {:?} is not a number",
                            String::from_utf8_lossy(key),
                            String::from_utf8_lossy(v)
                        )
                    })?;
                cb(Val::Num(n))
            }
        })
    };
    match idx {
        Some(i) => match value.split(|&b| b == b',').nth(i) {
            Some(v) => one(v, cb),
            None => Ok(cb(Val::Miss)),
        },
        None => {
            for v in value.split(|&b| b == b',') {
                if one(v, cb)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    }
}
