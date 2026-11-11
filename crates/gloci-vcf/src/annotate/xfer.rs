//! Column resolution and the per-record setters of `vcf annotate`, after
//! bcftools' `init_columns` and its `vcf_setter_*` and `setter_*` functions.

use super::hdr::{HeaderEdit, InfoDef, Number, Ty};
use super::rec::{Rec, info_get, info_remove, info_set};
use super::spec::{ColKind, ColumnSpec};
use super::tab::{TabLayout, TabRow};
use super::vcmp::Vcmp;

#[derive(Clone, Debug)]
pub struct InfoCol {
    dst: Vec<u8>,
    src: Vec<u8>,
    /// The destination's declaration (type and Number decide the transfer).
    def: InfoDef,
    /// The source VCF's declaration, if it has one.
    src_def: Option<InfoDef>,
}

#[derive(Clone, Debug)]
enum CKind {
    Placeholder,
    Id,
    Qual,
    Filter,
    Info(Box<InfoCol>),
}

/// A column to transfer.
#[derive(Clone, Debug)]
pub struct Col {
    kind: CKind,
    plus: bool,
    carry: bool,
    icol: usize,
}

impl Col {
    pub fn placeholder() -> Col {
        Col {
            kind: CKind::Placeholder,
            plus: false,
            carry: false,
            icol: 0,
        }
    }
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

fn at(rec: &Rec) -> String {
    format!("{}:{}", lossy(&rec.chrom), lossy(&rec.pos))
}

pub fn resolve_vcf(
    spec: &ColumnSpec,
    hdr_out: &mut HeaderEdit,
    src: &HeaderEdit,
) -> Result<Vec<Col>, String> {
    let mut cols = Vec::new();
    let mk = |kind, it: &super::spec::ColItem| Col {
        kind,
        plus: it.plus,
        carry: it.carry,
        icol: it.icol,
    };
    for it in &spec.items {
        match &it.kind {
            ColKind::Skip | ColKind::Chrom | ColKind::Pos | ColKind::To => {}
            ColKind::Ref | ColKind::Alt => unreachable!("rejected when parsing"),
            ColKind::Id => cols.push(mk(CKind::Id, it)),
            ColKind::Qual => cols.push(mk(CKind::Qual, it)),
            ColKind::Filter => {
                for id in src.ids("FILTER") {
                    if let Some(line) = src.line("FILTER", &id) {
                        hdr_out.append(line)?;
                    }
                }
                cols.push(mk(CKind::Filter, it));
            }
            ColKind::InfoAll => {
                for id in src.ids("INFO") {
                    if spec.skip_info.contains(&id) {
                        continue;
                    }
                    if !hdr_out.has("INFO", &id)
                        && let Some(line) = src.line("INFO", &id)
                    {
                        hdr_out.append(line)?;
                    }
                    cols.push(info_col(it, &id, &id, hdr_out, src)?);
                }
            }
            ColKind::Info { dst, src: s } => {
                if !hdr_out.has("INFO", dst) {
                    let Some(line) = src.line("INFO", s) else {
                        return Err(format!(
                            "--columns: the INFO tag {} is not defined in the annotation file; declare it with --header-lines or give a tag it declares",
                            lossy(s)
                        ));
                    };
                    let line = line.to_vec();
                    hdr_out.append_renamed("INFO", &line, dst);
                }
                cols.push(info_col(it, dst, s, hdr_out, src)?);
            }
        }
    }
    Ok(cols)
}

fn info_col(
    it: &super::spec::ColItem,
    dst: &[u8],
    src_tag: &[u8],
    hdr_out: &HeaderEdit,
    src: &HeaderEdit,
) -> Result<Col, String> {
    let def = hdr_out
        .info_def(dst)
        .ok_or_else(|| format!("the INFO tag {} could not be declared", lossy(dst)))?;
    Ok(Col {
        kind: CKind::Info(Box::new(InfoCol {
            dst: dst.to_vec(),
            src: src_tag.to_vec(),
            def,
            src_def: src.info_def(src_tag),
        })),
        plus: it.plus,
        carry: it.carry,
        icol: it.icol,
    })
}

pub fn resolve_tab(
    spec: &ColumnSpec,
    hdr_out: &mut HeaderEdit,
    bed: bool,
) -> Result<(Vec<Col>, TabLayout), String> {
    let (mut chr, mut beg, mut end, mut ref_, mut alt) = (None, None, None, None, None);
    let mut cols = Vec::new();
    let mut info_ar = false;
    for it in &spec.items {
        let mk = |kind| Col {
            kind,
            plus: it.plus,
            carry: it.carry,
            icol: it.icol,
        };
        match &it.kind {
            ColKind::Skip => {}
            ColKind::Chrom => chr = Some(it.icol),
            ColKind::Pos => beg = Some(it.icol),
            ColKind::To => end = Some(it.icol),
            ColKind::Ref => ref_ = Some(it.icol),
            ColKind::Alt => alt = Some(it.icol),
            ColKind::Id => cols.push(mk(CKind::Id)),
            ColKind::Qual => cols.push(mk(CKind::Qual)),
            ColKind::Filter => cols.push(mk(CKind::Filter)),
            ColKind::InfoAll => {
                return Err(
                    "--columns INFO (every tag) needs a VCF/BCF annotation file; name the tags for a tab or BED file"
                        .to_string(),
                );
            }
            ColKind::Info { dst, src } => {
                if dst != src {
                    // NEW:=OLD names a column of a VCF; in a tab file the
                    // position in the list is the column, so only NEW matters.
                    return Err(format!(
                        "--columns {}:={}: renaming with := needs a VCF/BCF annotation file; for a tab or BED file name the destination tag, the column is its position in the list",
                        lossy(dst),
                        lossy(src)
                    ));
                }
                let Some(def) = hdr_out.info_def(dst) else {
                    return Err(format!(
                        "--columns: the INFO tag {} is not defined in the input header; declare it with --header-lines, e.g. a file with ##INFO=<ID={},Number=1,Type=String,Description=\"...\">",
                        lossy(dst),
                        lossy(dst)
                    ));
                };
                if matches!(def.number, Number::A | Number::R) && def.ty != Ty::Flag {
                    info_ar = true;
                }
                cols.push(Col {
                    kind: CKind::Info(Box::new(InfoCol {
                        dst: dst.clone(),
                        src: src.clone(),
                        def,
                        src_def: None,
                    })),
                    plus: it.plus,
                    carry: it.carry,
                    icol: it.icol,
                });
            }
        }
    }
    let Some(chr) = chr else {
        return Err("--columns needs CHROM to say which column holds the chromosome".to_string());
    };
    let Some(beg) = beg else {
        return Err(
            "--columns needs POS (or FROM and TO) to say which column holds the position"
                .to_string(),
        );
    };
    if end.is_some() && ref_.is_some() {
        return Err(
            "--columns: REF is ignored when FROM,TO columns are present; replace TO with - or drop REF and ALT"
                .to_string(),
        );
    }
    if ref_.is_some() != alt.is_some() {
        return Err("--columns: REF and ALT must be given together".to_string());
    }
    if end.is_some() && (chr, beg, end) != (0, 1, Some(2)) {
        return Err(
            "--columns: with FROM,TO the first three columns must be CHROM,FROM,TO (as bcftools requires)"
                .to_string(),
        );
    }
    if info_ar && ref_.is_none() {
        return Err(
            "--columns: a Number=A or Number=R tag needs REF and ALT columns in the annotation file, to know which allele each value belongs to"
                .to_string(),
        );
    }
    Ok((
        cols,
        TabLayout {
            chr,
            beg,
            end,
            ref_,
            alt,
            bed,
        },
    ))
}

fn set_flag(rec: &mut Rec, key: &[u8], on: bool) {
    let present = info_get(&rec.info, key).is_some();
    if on && !present {
        info_set(&mut rec.info, key, None);
        rec.dirty = true;
    } else if !on && info_remove(&mut rec.info, key) {
        rec.dirty = true;
    }
}

/// Whether the record already has a value for the tag: its first value is
/// not missing (a string: the whole value is not `.`).
fn has_value(rec: &Rec, key: &[u8], ty: Ty) -> bool {
    match info_get(&rec.info, key) {
        Some(Some(v)) if !v.is_empty() => {
            if ty == Ty::Str {
                v != b"."
            } else {
                v.split(|&b| b == b',').next() != Some(b".")
            }
        }
        _ => false,
    }
}

fn put(rec: &mut Rec, key: &[u8], val: &[u8]) {
    if info_get(&rec.info, key) != Some(Some(val)) {
        info_set(&mut rec.info, key, Some(val));
        rec.dirty = true;
    }
}

fn set_text(dst: &mut Vec<u8>, v: &[u8], dirty: &mut bool) {
    if dst != v {
        *dst = v.to_vec();
        *dirty = true;
    }
}

fn alleles(rec: &Rec) -> Vec<Vec<u8>> {
    let mut v = vec![rec.reference.clone()];
    v.extend(rec.alts().into_iter().map(<[u8]>::to_vec));
    v
}

/// The Number=A/R transfer of bcftools' `setter_ARinfo_*`: the values of the
/// source record go to the destination's alleles by allele, not by position.
fn ar_transfer(
    ic: &InfoCol,
    col: &Col,
    rec: &mut Rec,
    src_vals: &[Vec<u8>],
    src_als: &[Vec<u8>],
    vcmp: &mut Vcmp,
) -> Result<(), String> {
    let is_a = ic.def.number == Number::A;
    let nals = src_als.len();
    let want = if is_a { nals - 1 } else { nals };
    let lone_missing = src_vals.len() == 1 && src_vals[0] == b".";
    // bcftools lets a lone `.` through for strings only; for numbers every
    // other count (a lone `.` included) is an error.
    if src_vals.len() != want && !(ic.def.ty == Ty::Str && lone_missing) {
        return Err(format!(
            "incorrect number of values ({}) for the {} tag at {}: its Number is {}, which needs {want} for a record with {} ALT allele(s)",
            src_vals.len(),
            lossy(&ic.src),
            at(rec),
            if is_a { "A" } else { "R" },
            nals - 1,
        ));
    }
    let dst_als = alleles(rec);
    let ndst = if is_a {
        dst_als.len() - 1
    } else {
        dst_als.len()
    };
    let s: Vec<&[u8]> = src_als.iter().map(Vec::as_slice).collect();
    let d: Vec<&[u8]> = dst_als.iter().map(Vec::as_slice).collect();
    let Some(map) = vcmp.map_ar_values(ndst, &s, &d) else {
        return Err(format!("REF alleles not compatible at {}", at(rec)));
    };
    let split =
        |v: &[u8]| -> Vec<Vec<u8>> { v.split(|&b| b == b',').map(<[u8]>::to_vec).collect() };
    let mut out: Vec<Vec<u8>>;
    if ic.def.ty == Ty::Str {
        let existing = match info_get(&rec.info, &ic.dst) {
            Some(Some(v)) if v != b"." => split(v),
            _ => vec![b".".to_vec(); ndst],
        };
        out = existing.clone();
        for (i, m) in map.iter().enumerate().take(ndst) {
            let Some(m) = m else { continue };
            // bcftools fills only missing values of a String A/R tag, with or
            // without +.
            if out.get(i).is_none_or(|t| t != b".") {
                continue;
            }
            let v = src_vals.get(*m).cloned().unwrap_or_else(|| b".".to_vec());
            while out.len() <= i {
                out.push(b".".to_vec());
            }
            out[i] = v;
        }
    } else {
        // bcftools reads an Integer tag's existing values as floats, which
        // fails, so for Integer the existing values never count.
        let existing: Vec<Vec<u8>> = if ic.def.ty == Ty::Float {
            match info_get(&rec.info, &ic.dst) {
                Some(Some(v)) if !v.is_empty() => split(v),
                _ => Vec::new(),
            }
        } else {
            Vec::new()
        };
        let n2 = existing.len();
        out = existing.clone();
        while out.len() < ndst {
            out.push(b".".to_vec());
        }
        for i in 0..ndst {
            match map[i] {
                None => {
                    if n2 < ndst {
                        out[i] = b".".to_vec();
                    }
                }
                Some(m) => {
                    if n2 == ndst && col.plus && existing[i] != b"." {
                        continue;
                    }
                    out[i] = src_vals[m].clone();
                }
            }
        }
    }
    if ndst == 0 && ic.def.ty != Ty::Str {
        if info_remove(&mut rec.info, &ic.dst) {
            rec.dirty = true;
        }
        return Ok(());
    }
    if ic.def.ty != Ty::Str {
        out.truncate(ndst);
    }
    put(rec, &ic.dst, &out.join(&b','));
    Ok(())
}

/// Transfer one column from the paired source VCF record.
pub fn vcf_set(col: &Col, rec: &mut Rec, src: &Rec, vcmp: &mut Vcmp) -> Result<(), String> {
    match &col.kind {
        CKind::Placeholder => {}
        CKind::Id => {
            if src.id == b"." || src.id.is_empty() {
                return Ok(());
            }
            if col.plus && rec.id != b"." && !rec.id.is_empty() {
                return Ok(());
            }
            set_text(&mut rec.id, &src.id, &mut rec.dirty);
        }
        CKind::Qual => {
            if src.qual == b"." {
                return Ok(());
            }
            if col.plus && rec.qual != b"." {
                return Ok(());
            }
            set_text(&mut rec.qual, &src.qual, &mut rec.dirty);
        }
        CKind::Filter => {
            if src.filter == b"." {
                return Ok(());
            }
            if col.plus && rec.filter != b"." {
                return Ok(());
            }
            set_text(&mut rec.filter, &src.filter, &mut rec.dirty);
        }
        CKind::Info(ic) => {
            let Some(sd) = ic.src_def else {
                return Ok(());
            };
            if ic.def.ty == Ty::Flag {
                let on = info_get(&src.info, &ic.src).is_some();
                set_flag(rec, &ic.dst, on);
                return Ok(());
            }
            if sd.ty != ic.def.ty {
                // bcftools reads the tag as the destination's type, which
                // fails, and adds nothing.
                return Ok(());
            }
            let Some(Some(val)) = info_get(&src.info, &ic.src) else {
                return Ok(());
            };
            if matches!(ic.def.number, Number::A | Number::R) {
                let vals: Vec<Vec<u8>> = val.split(|&b| b == b',').map(<[u8]>::to_vec).collect();
                return ar_transfer(ic, col, rec, &vals, &alleles(src), vcmp);
            }
            if col.plus && has_value(rec, &ic.dst, ic.def.ty) {
                return Ok(());
            }
            put(rec, &ic.dst, val);
        }
    }
    Ok(())
}

fn escape(cell: &[u8]) -> Vec<u8> {
    let mut o = Vec::with_capacity(cell.len());
    for &b in cell {
        match b {
            b';' => o.extend_from_slice(b"%3B"),
            b'=' => o.extend_from_slice(b"%3D"),
            _ => o.push(b),
        }
    }
    o
}

fn parse_tokens(ic: &InfoCol, col: &Col, cell: &[u8], rec: &Rec) -> Result<Vec<Vec<u8>>, String> {
    let mut out = Vec::new();
    let bad = |tok: &[u8]| {
        format!(
            "could not parse {} at {}: `{}` in `{}` is not {}",
            lossy(&ic.src),
            at(rec),
            lossy(tok),
            lossy(cell),
            if ic.def.ty == Ty::Int {
                "an integer"
            } else {
                "a number"
            }
        )
    };
    for tok in cell.split(|&b| b == b',') {
        if tok == b"." {
            if col.carry {
                out.push(b".".to_vec());
            }
            continue;
        }
        let s = std::str::from_utf8(tok).map_err(|_| bad(tok))?;
        if ic.def.ty == Ty::Int {
            let v: i64 = s.parse().map_err(|_| bad(tok))?;
            out.push(v.to_string().into_bytes());
        } else {
            s.parse::<f64>().map_err(|_| bad(tok))?;
            out.push(tok.to_vec());
        }
    }
    Ok(out)
}

/// Transfer one column from a tab/BED row. `Ok(true)` when the column is
/// done for this record (a later row will not set it), `Ok(false)` when the
/// row gave it nothing and a later row may.
pub fn tab_set(
    col: &Col,
    rec: &mut Rec,
    row: &TabRow,
    layout: &TabLayout,
    vcmp: &mut Vcmp,
    filter_ok: &dyn Fn(&[u8]) -> bool,
) -> Result<bool, String> {
    let cell = match col.kind {
        CKind::Placeholder => return Ok(true),
        _ => row.col(col.icol).ok_or_else(|| {
            format!(
                "the line has no column {}, which --columns names",
                col.icol + 1
            )
        })?,
    };
    match &col.kind {
        CKind::Placeholder => Ok(true),
        CKind::Id => {
            if cell == b"." {
                return Ok(true);
            }
            if col.plus && rec.id != b"." && !rec.id.is_empty() {
                return Ok(true);
            }
            set_text(&mut rec.id, cell, &mut rec.dirty);
            Ok(true)
        }
        CKind::Qual => {
            if cell == b"." {
                return Ok(true);
            }
            if col.plus && rec.qual != b"." {
                return Ok(true);
            }
            if std::str::from_utf8(cell)
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
                .is_none()
            {
                return Err(format!(
                    "could not parse QUAL at {}: `{}` is not a number",
                    at(rec),
                    lossy(cell)
                ));
            }
            set_text(&mut rec.qual, cell, &mut rec.dirty);
            Ok(true)
        }
        CKind::Filter => {
            if cell == b"." {
                return Ok(true);
            }
            if !filter_ok(cell) {
                return Err(format!(
                    "the FILTER `{}` is not defined in the header; declare it with --header-lines, e.g. ##FILTER=<ID={},Description=\"...\">",
                    lossy(cell),
                    lossy(cell)
                ));
            }
            if col.plus && rec.filter != b"." {
                return Ok(true);
            }
            set_text(&mut rec.filter, cell, &mut rec.dirty);
            Ok(true)
        }
        CKind::Info(ic) => {
            let is_ar = matches!(ic.def.number, Number::A | Number::R);
            match ic.def.ty {
                Ty::Flag => {
                    match cell {
                        b"." => {}
                        b"1" => set_flag(rec, &ic.dst, true),
                        b"0" => set_flag(rec, &ic.dst, false),
                        _ => {
                            return Err(format!(
                                "could not parse {} at {}: `{}` is not 0 or 1 (it is a flag)",
                                lossy(&ic.src),
                                at(rec),
                                lossy(cell)
                            ));
                        }
                    }
                    Ok(true)
                }
                Ty::Str => {
                    if col.plus && !is_ar && has_value(rec, &ic.dst, Ty::Str) {
                        return Ok(true);
                    }
                    if cell.is_empty() {
                        return Ok(true);
                    }
                    if cell == b"." && !col.carry {
                        return Ok(false);
                    }
                    let esc = escape(cell);
                    if is_ar {
                        let (als, _) = row
                            .alleles(layout)
                            .ok_or_else(|| "the line has no REF or ALT column".to_string())?;
                        let vals: Vec<Vec<u8>> =
                            esc.split(|&b| b == b',').map(<[u8]>::to_vec).collect();
                        ar_transfer(ic, col, rec, &vals, &als, vcmp)?;
                    } else {
                        put(rec, &ic.dst, &esc);
                    }
                    Ok(true)
                }
                Ty::Int | Ty::Float => {
                    if cell == b"." && !col.carry {
                        return Ok(false);
                    }
                    let vals = parse_tokens(ic, col, cell, rec)?;
                    if is_ar {
                        let (als, _) = row
                            .alleles(layout)
                            .ok_or_else(|| "the line has no REF or ALT column".to_string())?;
                        ar_transfer(ic, col, rec, &vals, &als, vcmp)?;
                        return Ok(true);
                    }
                    if col.plus && has_value(rec, &ic.dst, ic.def.ty) {
                        return Ok(true);
                    }
                    if vals.is_empty() {
                        if info_remove(&mut rec.info, &ic.dst) {
                            rec.dirty = true;
                        }
                    } else {
                        put(rec, &ic.dst, &vals.join(&b','));
                    }
                    Ok(true)
                }
            }
        }
    }
}
