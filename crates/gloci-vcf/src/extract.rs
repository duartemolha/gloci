//! Cells of `vcf extract`: INFO values, FORMAT values of any key (as text)
//! and the 012 genotype codes, read from a VCF text line (decision D73).

/// How the header declares an INFO key.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoKind {
    /// `Type=Flag`: written `1` when present, `0` when absent.
    Flag,
    /// Any other declared type: the value, or `.` when absent.
    Value,
    /// `Number=0` with a type other than Flag: `1` when present, `.` when absent
    /// (no `0`: the header does not call it a flag).
    Presence,
}

/// An INFO key to tabulate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InfoKey {
    pub name: Vec<u8>,
    pub kind: InfoKind,
}

/// The `KEY=value` pairs of a header line's `<...>` body, with quotes
/// respected (a `Description` may hold commas and `=`).
pub(crate) fn structured_pairs(line: &[u8]) -> Vec<(&[u8], &[u8])> {
    let Some(open) = line.iter().position(|&b| b == b'<') else {
        return Vec::new();
    };
    let body = line[open + 1..]
        .strip_suffix(b">")
        .unwrap_or(&line[open + 1..]);
    let mut pairs = Vec::new();
    let (mut start, mut quoted) = (0, false);
    for i in 0..=body.len() {
        let at_end = i == body.len();
        if !at_end && body[i] == b'"' {
            quoted = !quoted;
        }
        if at_end || (body[i] == b',' && !quoted) {
            let item = &body[start..i];
            if let Some(eq) = item.iter().position(|&b| b == b'=') {
                pairs.push((&item[..eq], &item[eq + 1..]));
            }
            start = i + 1;
        }
    }
    pairs
}

/// A header declaration: `ID`, `Type` and `Number`.
pub(crate) struct Decl<'a> {
    pub(crate) id: &'a [u8],
    pub(crate) ty: Option<&'a [u8]>,
    pub(crate) number: Option<&'a [u8]>,
}

/// The header lines `##<kind>=<...>` with their `ID`, `Type` and `Number`.
pub(crate) fn declarations<'a>(
    header: &'a [u8],
    kind: &'static [u8],
) -> impl Iterator<Item = Decl<'a>> {
    header.split(|&b| b == b'\n').filter_map(move |line| {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let rest = line.strip_prefix(b"##")?.strip_prefix(kind)?;
        rest.strip_prefix(b"=<")?;
        let pairs = structured_pairs(rest);
        let get = |k: &[u8]| pairs.iter().find(|(n, _)| *n == k).map(|p| p.1);
        Some(Decl {
            id: get(b"ID")?,
            ty: get(b"Type"),
            number: get(b"Number"),
        })
    })
}

/// How the header `header` (raw VCF header text) declares INFO key `name`;
/// `None` when it does not.
pub fn info_kind(header: &[u8], name: &[u8]) -> Option<InfoKind> {
    let d = declarations(header, b"INFO").find(|d| d.id == name)?;
    Some(match (d.ty, d.number) {
        (Some(b"Flag"), _) => InfoKind::Flag,
        (_, Some(b"0")) => InfoKind::Presence,
        _ => InfoKind::Value,
    })
}

/// Whether the header declares FORMAT key `name`.
pub fn format_declared(header: &[u8], name: &[u8]) -> bool {
    declarations(header, b"FORMAT").any(|d| d.id == name)
}

/// The IDs the header declares for `kind` (`INFO` or `FORMAT`), in order.
pub fn declared_keys(header: &[u8], kind: &'static [u8]) -> Vec<String> {
    declarations(header, kind)
        .map(|d| String::from_utf8_lossy(d.id).into_owned())
        .collect()
}

/// The tab-separated columns of `line` into `out` (cleared first). A VCF
/// record has at least eight; an error says what is missing.
pub fn split_columns<'a>(line: &'a [u8], out: &mut Vec<&'a [u8]>) -> Result<(), String> {
    out.clear();
    out.extend(line.split(|&b| b == b'\t'));
    if out.len() < 8 {
        return Err(format!(
            "the record has {} tab-separated columns, but a VCF record needs at least 8",
            out.len()
        ));
    }
    Ok(())
}

/// The value of `key` in an INFO column, as the table cell shows it.
pub fn info_cell<'a>(info: &'a [u8], key: &InfoKey) -> &'a [u8] {
    if info != b"." {
        for item in info.split(|&b| b == b';') {
            let (name, value) = match item.iter().position(|&b| b == b'=') {
                Some(eq) => (&item[..eq], Some(&item[eq + 1..])),
                None => (item, None),
            };
            if name != key.name.as_slice() {
                continue;
            }
            return match (value, key.kind) {
                (Some(v), _) if !v.is_empty() => v,
                (Some(_), _) => b".",
                (None, InfoKind::Value) => b".",
                (None, _) => b"1",
            };
        }
    }
    if key.kind == InfoKind::Flag {
        b"0"
    } else {
        b"."
    }
}

/// The position of `key` among the colon-separated keys of a FORMAT column.
pub fn format_index(format: &[u8], key: &[u8]) -> Option<usize> {
    if format == b"." {
        return None;
    }
    format.split(|&b| b == b':').position(|k| k == key)
}

/// Field `idx` of a sample column as the table cell shows it: `.` when the
/// sample has fewer fields (trailing fields may be dropped) or an empty one.
pub fn format_cell(sample: &[u8], idx: Option<usize>) -> &[u8] {
    let Some(idx) = idx else {
        return b".";
    };
    match sample.split(|&b| b == b':').nth(idx) {
        Some(v) if !v.is_empty() => v,
        _ => b".",
    }
}

/// Why a genotype has no 012 code.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GtError {
    pub msg: String,
    /// The genotype names an allele the ALT column lacks (the one error
    /// `--skip-invalid` turns into `-1`).
    pub allele_out_of_range: bool,
}

/// The 012 code of a genotype on a site with `n_alt` ALT alleles: the number
/// of non-reference alleles among the called ones, `-1` when no allele is
/// called. Ploidy above 2 and an allele the ALT column lacks are errors.
pub fn genotype_012(gt: &[u8], n_alt: usize) -> Result<i8, GtError> {
    let mut called = 0u8;
    let mut non_ref = 0i8;
    let mut ploidy = 0usize;
    for allele in gt.split(|&b| b == b'/' || b == b'|') {
        ploidy += 1;
        if ploidy > 2 {
            return Err(GtError {
                msg: format!(
                    "genotype {} has ploidy above 2, which the 012 matrix does not support; \
drop the site (--exclude-positions or --chrom) or the sample",
                    String::from_utf8_lossy(gt)
                ),
                allele_out_of_range: false,
            });
        }
        if allele == b"." || allele.is_empty() {
            continue;
        }
        let n: usize = std::str::from_utf8(allele)
            .ok()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| GtError {
                msg: format!("genotype {} is not a valid GT", String::from_utf8_lossy(gt)),
                allele_out_of_range: false,
            })?;
        if n > n_alt {
            return Err(GtError {
                msg: format!(
                    "genotype {} names allele {n}, but the site has {n_alt} ALT allele(s); \
fix the record or use --skip-invalid to write -1",
                    String::from_utf8_lossy(gt)
                ),
                allele_out_of_range: true,
            });
        }
        called += 1;
        non_ref += i8::from(n != 0);
    }
    Ok(if called == 0 { -1 } else { non_ref })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(name: &str, kind: InfoKind) -> InfoKey {
        InfoKey {
            name: name.into(),
            kind,
        }
    }

    #[test]
    fn info_cells() {
        let info = b"DP=40;AF=0.5,0.2;FLAG;EMPTY=";
        assert_eq!(info_cell(info, &key("DP", InfoKind::Value)), b"40");
        assert_eq!(info_cell(info, &key("AF", InfoKind::Value)), b"0.5,0.2");
        assert_eq!(info_cell(info, &key("FLAG", InfoKind::Flag)), b"1");
        assert_eq!(info_cell(info, &key("NOPE", InfoKind::Flag)), b"0");
        assert_eq!(info_cell(info, &key("NOPE", InfoKind::Value)), b".");
        assert_eq!(info_cell(info, &key("NOPE", InfoKind::Presence)), b".");
        assert_eq!(info_cell(info, &key("FLAG", InfoKind::Presence)), b"1");
        assert_eq!(info_cell(info, &key("FLAG", InfoKind::Value)), b".");
        assert_eq!(info_cell(info, &key("EMPTY", InfoKind::Value)), b".");
        assert_eq!(info_cell(b".", &key("DP", InfoKind::Value)), b".");
        assert_eq!(info_cell(b".", &key("FLAG", InfoKind::Flag)), b"0");
        // A key that merely starts with the name is another key.
        assert_eq!(info_cell(b"DPX=3", &key("DP", InfoKind::Value)), b".");
    }

    #[test]
    fn format_cells() {
        let idx = format_index(b"GT:DP:GQ", b"DP");
        assert_eq!(idx, Some(1));
        assert_eq!(format_cell(b"0/1:12:30", idx), b"12");
        assert_eq!(format_cell(b"0/1", idx), b".");
        assert_eq!(format_cell(b"0/1::5", idx), b".");
        assert_eq!(format_cell(b"0/1:.:5", idx), b".");
        assert_eq!(format_cell(b"0/1:3", format_index(b"GT:DP", b"AD")), b".");
        assert_eq!(format_index(b".", b"GT"), None);
    }

    #[test]
    fn header_declarations() {
        let h = b"##fileformat=VCFv4.2\n##INFO=<ID=A,Number=0,Type=Flag,Description=\"x,Type=Integer\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n";
        assert_eq!(info_kind(h, b"A"), Some(InfoKind::Flag));
        assert_eq!(info_kind(h, b"DP"), Some(InfoKind::Value));
        assert_eq!(info_kind(h, b"Z"), None);
        assert_eq!(declared_keys(h, b"INFO"), ["A", "DP"]);
        assert!(format_declared(h, b"GT"));
        assert!(!format_declared(h, b"DP"));
    }

    #[test]
    fn genotype_codes() {
        let g = |s: &str, n| genotype_012(s.as_bytes(), n);
        assert_eq!(g("0/0", 1), Ok(0));
        assert_eq!(g("0/1", 1), Ok(1));
        assert_eq!(g("1|1", 1), Ok(2));
        assert_eq!(g("./.", 1), Ok(-1));
        assert_eq!(g(".", 1), Ok(-1));
        assert_eq!(g("./1", 1), Ok(1));
        assert_eq!(g("0/.", 1), Ok(0));
        assert_eq!(g("0", 1), Ok(0));
        assert_eq!(g("1", 1), Ok(1));
        assert_eq!(g("0/0", 0), Ok(0));
        assert!(g("0/1/1", 1).unwrap_err().msg.contains("ploidy"));
        let e = g("2/2", 1).unwrap_err();
        assert!(e.msg.contains("allele 2") && e.allele_out_of_range);
        let e = g("a/b", 1).unwrap_err();
        assert!(e.msg.contains("not a valid GT") && !e.allele_out_of_range);
    }

    #[test]
    fn columns() {
        let mut v = Vec::new();
        assert!(split_columns(b"a\tb", &mut v).is_err());
        split_columns(b"1\t2\t3\t4\t5\t6\t7\t8\t9", &mut v).unwrap();
        assert_eq!(v.len(), 9);
    }
}
