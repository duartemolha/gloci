//! Cells of `vcf extract`: INFO values, FORMAT values of any key (as text)
//! and the 012 genotype codes, read from a VCF text line (decision D73).

/// How the header declares an INFO key.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoKind {
    /// `Type=Flag`: written `1` when present, `0` when absent.
    Flag,
    /// Any other declared type: the value, or `.` when absent.
    Value,
    /// Not declared in the header: the value, `1` for a bare key, `.` when absent.
    Undeclared,
}

/// An INFO key to tabulate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InfoKey {
    pub name: Vec<u8>,
    pub kind: InfoKind,
}

/// The `KEY=value` pairs of a header line's `<...>` body, with quotes
/// respected (a `Description` may hold commas and `=`).
fn structured_pairs(line: &[u8]) -> Vec<(&[u8], &[u8])> {
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

/// The header lines `##<kind>=<...>` with their `ID` and `Type` values.
fn declarations<'a>(
    header: &'a [u8],
    kind: &'static [u8],
) -> impl Iterator<Item = (&'a [u8], Option<&'a [u8]>)> {
    header.split(|&b| b == b'\n').filter_map(move |line| {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let rest = line.strip_prefix(b"##")?.strip_prefix(kind)?;
        rest.strip_prefix(b"=<")?;
        let pairs = structured_pairs(rest);
        let id = pairs.iter().find(|(k, _)| *k == b"ID")?.1;
        let ty = pairs.iter().find(|(k, _)| *k == b"Type").map(|p| p.1);
        Some((id, ty))
    })
}

/// How the header `header` (raw VCF header text) declares INFO key `name`.
pub fn info_kind(header: &[u8], name: &[u8]) -> InfoKind {
    match declarations(header, b"INFO").find(|(id, _)| *id == name) {
        None => InfoKind::Undeclared,
        Some((_, Some(b"Flag"))) => InfoKind::Flag,
        Some(_) => InfoKind::Value,
    }
}

/// Whether the header declares FORMAT key `name`.
pub fn format_declared(header: &[u8], name: &[u8]) -> bool {
    declarations(header, b"FORMAT").any(|(id, _)| id == name)
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

/// The 012 code of a genotype on a site with `n_alt` ALT alleles: the number
/// of non-reference alleles among the called ones, `-1` when no allele is
/// called. Ploidy above 2 and an allele the ALT column lacks are errors.
pub fn genotype_012(gt: &[u8], n_alt: usize) -> Result<i8, String> {
    let mut called = 0u8;
    let mut non_ref = 0i8;
    let mut ploidy = 0usize;
    for allele in gt.split(|&b| b == b'/' || b == b'|') {
        ploidy += 1;
        if ploidy > 2 {
            return Err(format!(
                "genotype {} has ploidy above 2, which the 012 matrix does not support; \
drop the site (--exclude-positions or --chrom) or the sample",
                String::from_utf8_lossy(gt)
            ));
        }
        if allele == b"." || allele.is_empty() {
            continue;
        }
        let n: usize = std::str::from_utf8(allele)
            .ok()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| format!("genotype {} is not a valid GT", String::from_utf8_lossy(gt)))?;
        if n > n_alt {
            return Err(format!(
                "genotype {} names allele {n}, but the site has {n_alt} ALT allele(s); \
fix the record",
                String::from_utf8_lossy(gt)
            ));
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
        assert_eq!(info_cell(info, &key("NOPE", InfoKind::Undeclared)), b".");
        assert_eq!(info_cell(info, &key("FLAG", InfoKind::Undeclared)), b"1");
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
        assert_eq!(info_kind(h, b"A"), InfoKind::Flag);
        assert_eq!(info_kind(h, b"DP"), InfoKind::Value);
        assert_eq!(info_kind(h, b"Z"), InfoKind::Undeclared);
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
        assert!(g("0/1/1", 1).unwrap_err().contains("ploidy"));
        assert!(g("2/2", 1).unwrap_err().contains("allele 2"));
        assert!(g("a/b", 1).unwrap_err().contains("not a valid GT"));
    }

    #[test]
    fn columns() {
        let mut v = Vec::new();
        assert!(split_columns(b"a\tb", &mut v).is_err());
        split_columns(b"1\t2\t3\t4\t5\t6\t7\t8\t9", &mut v).unwrap();
        assert_eq!(v.len(), 9);
    }
}
