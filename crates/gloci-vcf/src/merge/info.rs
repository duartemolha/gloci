//! The INFO column of a merged record: bcftools merge's `--info-rules`, the
//! re-indexing of Number=A and Number=R tags, and AC/AN.

use super::format::Calls;
use super::header::Merged;
use super::{Alleles, MergeOptions, Merger, P, parse_f32, put_f32};
use crate::annotate::hdr::{InfoDef, Number, Ty};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Method {
    Sum,
    Avg,
    Min,
    Max,
    Join,
}

#[derive(Clone, Debug)]
pub(crate) struct Rule {
    pub tag: Vec<u8>,
    pub method: Method,
    pub def: InfoDef,
}

/// INFO/END of a record, if any.
pub(crate) fn end_value(info: &[u8]) -> Option<i64> {
    items(info)
        .into_iter()
        .find(|(k, _)| *k == b"END")
        .and_then(|(_, v)| std::str::from_utf8(v?).ok()?.parse().ok())
}

/// One `KEY` or `KEY=VALUE` of an INFO column.
type Item<'a> = (&'a [u8], Option<&'a [u8]>);

fn items(info: &[u8]) -> Vec<Item<'_>> {
    if info == b"." || info.is_empty() {
        return Vec::new();
    }
    info.split(|&b| b == b';')
        .filter(|i| !i.is_empty())
        .map(|i| match i.iter().position(|&b| b == b'=') {
            Some(p) => (&i[..p], Some(&i[p + 1..])),
            None => (i, None),
        })
        .collect()
}

fn method_name(m: Method) -> &'static str {
    match m {
        Method::Sum => "sum",
        Method::Avg => "avg",
        Method::Min => "min",
        Method::Max => "max",
        Method::Join => "join",
    }
}

/// The rules in effect (in header order) and whether AC and AN are recomputed.
pub(crate) fn build_rules(
    opts: &MergeOptions,
    merged: &Merged,
) -> Result<(Vec<Rule>, bool), String> {
    let defs = &merged.defs;
    let sites_only = merged.samples.is_empty();
    let mut rules: Vec<Rule> = Vec::new();
    let mut user_named = false;
    match opts.info_rules.as_deref() {
        Some("-") => {}
        Some(list) => {
            for item in list.split(',').filter(|s| !s.is_empty()) {
                let (tag, method) = item.split_once(':').ok_or_else(|| {
                    format!(
                        "the --info-rules item \"{item}\" is not TAG:METHOD; use for example DP:sum (methods: sum, avg, min, max, join)"
                    )
                })?;
                let method = match method {
                    "sum" => Method::Sum,
                    "avg" => Method::Avg,
                    "min" => Method::Min,
                    "max" => Method::Max,
                    "join" => Method::Join,
                    _ => {
                        return Err(format!(
                            "the --info-rules method \"{method}\" is not recognised; use sum, avg, min, max or join"
                        ));
                    }
                };
                let Some((_, def)) = defs.info.get(tag.as_bytes()) else {
                    return Err(format!(
                        "the INFO tag \"{tag}\" of --info-rules is not defined in the header; add its ##INFO line to an input"
                    ));
                };
                if def.ty == Ty::Flag {
                    return Err(format!(
                        "the --info-rules item \"{item}\" is not supported: INFO/{tag} is a Flag"
                    ));
                }
                if method != Method::Join
                    && !matches!(def.number, Number::Fixed(_) | Number::A | Number::R)
                {
                    return Err(format!(
                        "--info-rules {}:{tag} needs a tag of fixed length (Number=1, a number, A or R); INFO/{tag} has a variable length",
                        method_name(method)
                    ));
                }
                if method != Method::Join && def.ty == Ty::Str {
                    return Err(format!(
                        "--info-rules {}:{tag} needs a numeric tag; INFO/{tag} is a String (only join works on strings)",
                        method_name(method)
                    ));
                }
                if tag == "AC" || tag == "AN" {
                    user_named = true;
                }
                rules.retain(|r| r.tag != tag.as_bytes());
                rules.push(Rule {
                    tag: tag.as_bytes().to_vec(),
                    method,
                    def: *def,
                });
            }
        }
        None => {
            let mut tags: Vec<&[u8]> = Vec::new();
            if sites_only {
                tags.extend([&b"AC"[..], b"AN"]);
            }
            tags.extend([&b"DP"[..], b"DP4"]);
            for t in tags {
                if let Some((_, def)) = defs.info.get(t)
                    && matches!(def.number, Number::Fixed(_) | Number::A | Number::R)
                    && matches!(def.ty, Ty::Int | Ty::Float)
                {
                    rules.push(Rule {
                        tag: t.to_vec(),
                        method: Method::Sum,
                        def: *def,
                    });
                }
            }
        }
    }
    rules.sort_by(|a, b| a.tag.cmp(&b.tag));
    let calc = !(user_named || opts.info_rules.is_none() && sites_only);
    Ok((rules, calc))
}

fn put_value(out: &mut Vec<u8>, ty: Ty, val: &[u8]) {
    for (i, t) in val.split(|&b| b == b',').enumerate() {
        if i > 0 {
            out.push(b',');
        }
        super::format::put_val(out, ty, t);
    }
}

fn push_item(out: &mut Vec<u8>, key: &[u8], val: Option<&[u8]>, ty: Ty) {
    if !out.is_empty() {
        out.push(b';');
    }
    out.extend_from_slice(key);
    if let Some(v) = val
        && ty != Ty::Flag
    {
        out.push(b'=');
        put_value(out, ty, v);
    }
}

pub(crate) fn merge_info(
    m: &Merger,
    recs: &[P<'_>],
    al: &Alleles,
    calls: &Option<Calls>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let calc = m.calc_acan;
    let declared = |k: &[u8]| m.defs.info.contains_key(k);
    let is_rule = |k: &[u8]| m.rules.iter().any(|r| r.tag == k);
    let is_ar =
        |k: &[u8]| declared(k) && matches!(m.out_def(true, k).number, Number::A | Number::R);
    let removed = |k: &[u8]| is_rule(k) || is_ar(k) || (calc && (k == b"AC" || k == b"AN"));
    let at = super::pos_label(&recs[0]);
    let all: Vec<Vec<Item<'_>>> = recs.iter().map(|r| items(r.info)).collect();

    // Tags of the first record keep their place.
    let mut written: Vec<&[u8]> = Vec::new();
    for (k, v) in &all[0] {
        written.push(k);
        if removed(k) {
            continue;
        }
        m.check_type(true, recs[0].reader, k, &at)?;
        push_item(out, k, *v, m.out_def(true, k).ty);
    }
    // Tags only later records have, as they come.
    for (ri, its) in all.iter().enumerate().skip(1) {
        for (k, v) in its {
            if written.contains(k) || removed(k) {
                continue;
            }
            written.push(k);
            m.check_type(true, recs[ri].reader, k, &at)?;
            push_item(out, k, *v, m.out_def(true, k).ty);
        }
    }
    // Tags with a rule, in header order.
    for rule in &m.rules {
        let mut vals: Vec<(usize, Vec<&[u8]>)> = Vec::new();
        for (ri, its) in all.iter().enumerate() {
            if let Some((_, Some(v))) = its.iter().find(|(k, _)| *k == rule.tag.as_slice()) {
                m.check_type(true, recs[ri].reader, &rule.tag, &at)?;
                vals.push((ri, v.split(|&b| b == b',').collect()));
            }
        }
        if vals.is_empty() {
            continue;
        }
        if let Some(text) = apply_rule(rule, al, &vals, &at)? {
            if !out.is_empty() {
                out.push(b';');
            }
            out.extend_from_slice(&rule.tag);
            out.push(b'=');
            out.extend_from_slice(&text);
        }
    }
    // Number=A and Number=R tags without a rule, in header order.
    let mut ar_keys: Vec<&[u8]> = Vec::new();
    for its in &all {
        for (k, _) in its {
            if is_ar(k)
                && !is_rule(k)
                && !(calc && (*k == b"AC" || *k == b"AN"))
                && !ar_keys.contains(k)
            {
                ar_keys.push(k);
            }
        }
    }
    for key in ar_keys {
        let def = m.out_def(true, key);
        let len = if def.number == Number::A {
            al.alts.len()
        } else {
            al.alts.len() + 1
        };
        let mut v: Vec<Option<&[u8]>> = vec![None; len];
        let mut any = false;
        for (ri, its) in all.iter().enumerate() {
            let Some((_, Some(val))) = its.iter().find(|(k, _)| *k == key) else {
                continue;
            };
            any = true;
            m.check_type(true, recs[ri].reader, key, &at)?;
            let map = &al.maps[ri];
            for (j, t) in val.split(|&b| b == b',').enumerate() {
                let slot = if def.number == Number::A {
                    map.get(j + 1).map(|s| s - 1)
                } else {
                    map.get(j).copied()
                };
                match slot {
                    // A missing value does not replace one already there.
                    Some(s) => {
                        if t != b"." || v[s].is_none() {
                            v[s] = Some(t);
                        }
                    }
                    None => {
                        return Err(format!(
                            "{}: INFO/{} at {at} has more values than the record has alleles",
                            m.inputs[recs[ri].reader].name,
                            String::from_utf8_lossy(key)
                        ));
                    }
                }
            }
        }
        if !any {
            continue;
        }
        if !out.is_empty() {
            out.push(b';');
        }
        out.extend_from_slice(key);
        out.push(b'=');
        for (j, x) in v.iter().enumerate() {
            if j > 0 {
                out.push(b',');
            }
            match x {
                Some(t) => super::format::put_val(out, def.ty, t),
                None => out.push(b'.'),
            }
        }
    }
    if calc && let Some(c) = calls {
        if declared(b"AN") {
            if !out.is_empty() {
                out.push(b';');
            }
            out.extend_from_slice(b"AN=");
            out.extend_from_slice(itoa::Buffer::new().format(c.an).as_bytes());
        }
        if declared(b"AC") && !c.ac.is_empty() {
            if !out.is_empty() {
                out.push(b';');
            }
            out.extend_from_slice(b"AC=");
            for (i, a) in c.ac.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(itoa::Buffer::new().format(*a).as_bytes());
            }
        }
    }
    Ok(())
}

fn apply_rule(
    rule: &Rule,
    al: &Alleles,
    vals: &[(usize, Vec<&[u8]>)],
    at: &str,
) -> Result<Option<Vec<u8>>, String> {
    let ty = rule.def.ty;
    let tag = String::from_utf8_lossy(&rule.tag).into_owned();
    if rule.method == Method::Join {
        if ty != Ty::Str && matches!(rule.def.number, Number::Var | Number::G) {
            let n = vals[0].1.len();
            if vals.iter().any(|(_, v)| v.len() != n) {
                return Err(format!(
                    "Mismatch in number of values for INFO/{tag} at {at}"
                ));
            }
        }
        let mut out = Vec::new();
        for (_, toks) in vals {
            for t in toks {
                if *t == b"." || t.is_empty() {
                    continue;
                }
                if !out.is_empty() {
                    out.push(b',');
                }
                super::format::put_val(&mut out, ty, t);
            }
        }
        return Ok((!out.is_empty()).then_some(out));
    }
    let nslots = match rule.def.number {
        Number::A => al.alts.len(),
        Number::R => al.alts.len() + 1,
        Number::Fixed(n) => n,
        _ => 1,
    };
    let mut slots: Vec<Vec<Option<&[u8]>>> = vec![Vec::new(); nslots];
    for (ri, toks) in vals {
        let map = &al.maps[*ri];
        for (j, t) in toks.iter().enumerate() {
            let s = match rule.def.number {
                Number::A => map.get(j + 1).map(|s| s - 1),
                Number::R => map.get(j).copied(),
                _ => (j < nslots).then_some(j),
            };
            let Some(s) = s else {
                return Err(format!(
                    "Mismatch in number of values for INFO/{tag} at {at}"
                ));
            };
            slots[s].push((*t != b"." && !t.is_empty()).then_some(*t));
        }
    }
    let mut out = Vec::new();
    let mut any = false;
    for (si, s) in slots.iter().enumerate() {
        if si > 0 {
            out.push(b',');
        }
        let present: Vec<&[u8]> = s.iter().flatten().copied().collect();
        if present.is_empty() {
            // Nothing to add up is 0; nothing to compare is missing.
            if matches!(rule.method, Method::Sum | Method::Avg) {
                out.push(b'0');
                any = true;
            } else {
                out.push(b'.');
            }
            continue;
        }
        any = true;
        if ty == Ty::Int {
            let v: Vec<i64> = present
                .iter()
                .map(|t| {
                    std::str::from_utf8(t)
                        .ok()
                        .and_then(|x| x.parse::<i64>().ok())
                        .ok_or_else(|| {
                            format!(
                                "INFO/{tag} at {at}: \"{}\" is not an integer",
                                String::from_utf8_lossy(t)
                            )
                        })
                })
                .collect::<Result<_, _>>()?;
            let r = match rule.method {
                Method::Sum => v.iter().fold(0i64, |a, b| a.saturating_add(*b)),
                Method::Avg => v.iter().fold(0i64, |a, b| a.saturating_add(*b)) / vals.len() as i64,
                Method::Min => *v.iter().min().unwrap_or(&0),
                _ => *v.iter().max().unwrap_or(&0),
            };
            out.extend_from_slice(itoa::Buffer::new().format(r).as_bytes());
        } else {
            let v: Vec<f32> = present
                .iter()
                .map(|t| {
                    parse_f32(t).ok_or_else(|| {
                        format!(
                            "INFO/{tag} at {at}: \"{}\" is not a number",
                            String::from_utf8_lossy(t)
                        )
                    })
                })
                .collect::<Result<_, _>>()?;
            let r = match rule.method {
                Method::Sum => v.iter().fold(0f32, |a, b| a + b),
                Method::Avg => v.iter().fold(0f32, |a, b| a + b) / vals.len() as f32,
                Method::Min => v.iter().copied().fold(f32::INFINITY, f32::min),
                _ => v.iter().copied().fold(f32::NEG_INFINITY, f32::max),
            };
            put_f32(&mut out, r);
        }
    }
    let _ = any;
    Ok(Some(out))
}
