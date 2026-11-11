//! Parsing of `vcf annotate`'s `--columns` list and `--set-id` format.

use super::rec::{Rec, info_get};

/// What one `--columns` item names.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ColKind {
    /// `-`, or an empty item: the column is not used.
    Skip,
    Chrom,
    /// `POS` or `FROM`/`BEG`: the start, 1-based (0-based in a BED file).
    Pos,
    /// `TO`/`END`: the end of an interval.
    To,
    Ref,
    Alt,
    Id,
    Qual,
    Filter,
    /// `INFO`: every INFO tag the source VCF declares.
    InfoAll,
    /// An INFO tag, renamed from `src` to `dst` when they differ.
    Info {
        dst: Vec<u8>,
        src: Vec<u8>,
    },
}

/// One item of the list.
#[derive(Clone, Debug)]
pub struct ColItem {
    pub kind: ColKind,
    /// `+TAG`: only fill missing values.
    pub plus: bool,
    /// `.TAG`: carry the source's missing values over (tab and BED sources).
    pub carry: bool,
    /// The item's position in the list, which is its column in a tab file.
    pub icol: usize,
    /// Whether `FROM` (rather than `POS`) named the start.
    pub from_form: bool,
}

/// A parsed `--columns` list.
#[derive(Clone, Debug, Default)]
pub struct ColumnSpec {
    pub items: Vec<ColItem>,
    /// INFO tags a `^TAG` item leaves out of `INFO`.
    pub skip_info: Vec<Vec<u8>>,
}

fn strip_ci<'a>(s: &'a str, prefix: &str) -> Option<&'a str> {
    (s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix))
        .then(|| &s[prefix.len()..])
}

/// Parse `list`. `vcf_src` is true for a VCF/BCF source, false for a tab or
/// BED file. Errors are usage messages.
pub fn parse_columns(list: &str, vcf_src: bool) -> Result<ColumnSpec, String> {
    let mut spec = ColumnSpec::default();
    let mut have_skip_list = false;
    let mut icol = 0usize;
    for raw in list.split(',') {
        // `^TAG` lists (VCF source): INFO without these tags.
        if let Some(rest) = raw.strip_prefix('^') {
            if !vcf_src {
                return Err(format!(
                    "--columns {raw}: ^ lists are for a VCF/BCF source; for a tab or BED file name the columns to transfer"
                ));
            }
            if strip_ci(rest, "FORMAT/").is_some() || strip_ci(rest, "FMT/").is_some() {
                return Err(format!(
                    "--columns {raw}: transferring FORMAT fields is not available in gloci vcf annotate"
                ));
            }
            let tag = strip_ci(rest, "INFO/")
                .or_else(|| strip_ci(rest, "INF/"))
                .unwrap_or(rest);
            if !have_skip_list {
                have_skip_list = true;
                spec.items.push(ColItem {
                    kind: ColKind::InfoAll,
                    plus: false,
                    carry: false,
                    icol,
                    from_form: false,
                });
                icol += 1;
            }
            spec.skip_info.push(tag.as_bytes().to_vec());
            continue;
        }
        let mut s = raw;
        let (mut plus, mut carry, mut minus, mut eq) = (false, false, false, false);
        while let Some(c) = s.chars().next() {
            match c {
                '+' => plus = true,
                '-' => minus = true,
                '=' => eq = true,
                '.' => carry = true,
                _ => break,
            }
            s = &s[1..];
        }
        let item = |kind: ColKind, from_form: bool| ColItem {
            kind,
            plus,
            carry,
            icol,
            from_form,
        };
        let kind;
        let mut from_form = false;
        if s.is_empty() {
            kind = ColKind::Skip;
        } else if minus {
            return Err(format!(
                "--columns {raw}: the -TAG form (overwrite only existing values) is not implemented by bcftools annotate either; use +TAG to fill missing values or TAG to overwrite"
            ));
        } else if eq {
            return Err(format!(
                "--columns {raw}: the =TAG form (append to existing values) is not available in gloci vcf annotate"
            ));
        } else if s.starts_with('~') {
            return Err(format!(
                "--columns {raw}: matching on a column (~ID, ~INFO/END) is not available in gloci vcf annotate"
            ));
        } else if s.eq_ignore_ascii_case("CHROM") {
            kind = ColKind::Chrom;
        } else if s.eq_ignore_ascii_case("POS") {
            kind = ColKind::Pos;
        } else if s.eq_ignore_ascii_case("FROM") || s.eq_ignore_ascii_case("BEG") {
            kind = ColKind::Pos;
            from_form = true;
        } else if s.eq_ignore_ascii_case("TO") || s.eq_ignore_ascii_case("END") {
            kind = ColKind::To;
        } else if s.eq_ignore_ascii_case("REF") || s.eq_ignore_ascii_case("ALT") {
            if vcf_src {
                return Err(format!(
                    "--columns {raw}: replacing REF or ALT from a VCF/BCF source is not available in gloci vcf annotate"
                ));
            }
            kind = if s.eq_ignore_ascii_case("REF") {
                ColKind::Ref
            } else {
                ColKind::Alt
            };
        } else if s.eq_ignore_ascii_case("ID") {
            kind = ColKind::Id;
        } else if s.eq_ignore_ascii_case("QUAL") {
            kind = ColKind::Qual;
        } else if s.eq_ignore_ascii_case("FILTER") {
            kind = ColKind::Filter;
        } else if vcf_src && s.eq_ignore_ascii_case("INFO") {
            kind = ColKind::InfoAll;
        } else if s.eq_ignore_ascii_case("FORMAT")
            || s.eq_ignore_ascii_case("FMT")
            || strip_ci(s, "FORMAT/").is_some()
            || strip_ci(s, "FMT/").is_some()
        {
            return Err(format!(
                "--columns {raw}: transferring FORMAT fields is not available in gloci vcf annotate"
            ));
        } else if s.len() >= 3 && s[..3].eq_ignore_ascii_case("ID:") {
            return Err(format!(
                "--columns {raw}: transferring an INFO tag into ID is not available in gloci vcf annotate"
            ));
        } else {
            let (dst, src) = match s.find(":=") {
                Some(p) => (&s[..p], &s[p + 2..]),
                None => (s, s),
            };
            let dst = strip_ci(dst, "INFO/").unwrap_or(dst);
            let src = strip_ci(src, "INFO/").unwrap_or(src);
            if strip_ci(src, "FMT/").is_some() || strip_ci(src, "FORMAT/").is_some() {
                return Err(format!(
                    "--columns {raw}: transferring FORMAT fields is not available in gloci vcf annotate"
                ));
            }
            if src.eq_ignore_ascii_case("FILTER") || src.eq_ignore_ascii_case("ID") {
                return Err(format!(
                    "--columns {raw}: transferring the ID or FILTER column into an INFO tag is not available in gloci vcf annotate"
                ));
            }
            if dst.is_empty() || src.is_empty() {
                return Err(format!("--columns {raw}: the tag name is empty"));
            }
            kind = ColKind::Info {
                dst: dst.as_bytes().to_vec(),
                src: src.as_bytes().to_vec(),
            };
        }
        let it = item(kind, from_form);
        spec.items.push(it);
        icol += 1;
    }
    Ok(spec)
}

/// One piece of a `--set-id` format.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IdTok {
    Lit(Vec<u8>),
    Chrom,
    Pos,
    Ref,
    FirstAlt,
    Id,
    Info(Vec<u8>),
}

/// A `--set-id` format, `+` prefix stripped.
#[derive(Clone, Debug)]
pub struct IdFormat {
    pub toks: Vec<IdTok>,
    /// A `+` prefix: set only missing IDs.
    pub plus: bool,
}

/// Parse a `bcftools query`-style format, the subset ID columns use. A tag
/// name is letters, digits, `_` and `.` (so `%POS_x` names the tag `POS_x`,
/// as in bcftools); put another character after the tag.
pub fn parse_id_format(fmt: &str) -> Result<IdFormat, String> {
    let (plus, fmt) = match fmt.strip_prefix('+') {
        Some(r) => (true, r),
        None => (false, fmt),
    };
    let b = fmt.as_bytes();
    let mut toks = Vec::new();
    let mut lit = Vec::new();
    let mut i = 0;
    let is_tag = |c: u8| c.is_ascii_alphanumeric() || c == b'_' || c == b'.';
    while i < b.len() {
        match b[i] {
            b'%' => {
                if !lit.is_empty() {
                    toks.push(IdTok::Lit(std::mem::take(&mut lit)));
                }
                let mut j = i + 1;
                while j < b.len() && is_tag(b[j]) {
                    j += 1;
                }
                let name = &fmt[i + 1..j];
                if name.is_empty() {
                    return Err(format!(
                        "--set-id {fmt}: a % must be followed by a tag name"
                    ));
                }
                let tok = match name {
                    "CHROM" => IdTok::Chrom,
                    "POS" => IdTok::Pos,
                    "REF" => IdTok::Ref,
                    "FIRST_ALT" => IdTok::FirstAlt,
                    "ID" => IdTok::Id,
                    "INFO" => {
                        if b.get(j) != Some(&b'/') {
                            return Err(format!(
                                "--set-id {fmt}: %INFO alone (the whole column) is not available; name a tag, %INFO/TAG"
                            ));
                        }
                        let mut k = j + 1;
                        while k < b.len() && is_tag(b[k]) {
                            k += 1;
                        }
                        let tag = &fmt[j + 1..k];
                        if tag.is_empty() {
                            return Err(format!("--set-id {fmt}: %INFO/ needs a tag name"));
                        }
                        j = k;
                        IdTok::Info(tag.as_bytes().to_vec())
                    }
                    other => {
                        return Err(format!(
                            "--set-id {fmt}: `%{other}` is not available; gloci vcf annotate supports %CHROM, %POS, %REF, %FIRST_ALT, %ID and %INFO/TAG (a tag name is letters, digits, _ and . only, so write %POS:%REF rather than %POS_%REF)"
                        ));
                    }
                };
                if matches!(b.get(j), Some(b'{')) {
                    return Err(format!(
                        "--set-id {fmt}: subscripts such as %INFO/AF{{0}} are not available in gloci vcf annotate"
                    ));
                }
                toks.push(tok);
                i = j;
            }
            b'[' | b']' => {
                return Err(format!(
                    "--set-id {fmt}: per-sample blocks [ ] cannot make an ID"
                ));
            }
            b'\\' => {
                i += 1;
                match b.get(i) {
                    Some(b'n') => lit.push(b'\n'),
                    Some(b't') => lit.push(b'\t'),
                    Some(&c) => lit.push(c),
                    None => {}
                }
                i += 1;
            }
            c => {
                lit.push(c);
                i += 1;
            }
        }
    }
    if !lit.is_empty() {
        toks.push(IdTok::Lit(lit));
    }
    if toks.is_empty() {
        return Err("--set-id needs a format, e.g. '%CHROM:%POS:%REF:%FIRST_ALT'".to_string());
    }
    Ok(IdFormat { toks, plus })
}

impl IdFormat {
    /// The INFO tags the format reads.
    pub fn info_tags(&self) -> impl Iterator<Item = &[u8]> {
        self.toks.iter().filter_map(|t| match t {
            IdTok::Info(k) => Some(k.as_slice()),
            _ => None,
        })
    }

    pub fn render(&self, rec: &Rec, out: &mut Vec<u8>) {
        for t in &self.toks {
            match t {
                IdTok::Lit(l) => out.extend_from_slice(l),
                IdTok::Chrom => out.extend_from_slice(&rec.chrom),
                IdTok::Pos => out.extend_from_slice(&rec.pos),
                IdTok::Ref => out.extend_from_slice(&rec.reference),
                IdTok::FirstAlt => {
                    let first = rec.alt.split(|&b| b == b',').next().unwrap_or(b".");
                    out.extend_from_slice(if first.is_empty() { b"." } else { first });
                }
                IdTok::Id => {
                    out.extend_from_slice(if rec.id.is_empty() { b"." } else { &rec.id });
                }
                IdTok::Info(k) => match info_get(&rec.info, k) {
                    None => out.push(b'.'),
                    Some(None) => out.push(b'1'),
                    Some(Some([])) => out.push(b'.'),
                    Some(Some(v)) => out.extend_from_slice(v),
                },
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn columns() {
        let s = parse_columns("CHROM,POS,-,ID,+INFO/AF,NEW:=INFO/OLD,.DB", false).unwrap();
        let k: Vec<_> = s
            .items
            .iter()
            .map(|i| (&i.kind, i.plus, i.carry, i.icol))
            .collect();
        assert_eq!(k[2], (&ColKind::Skip, false, false, 2));
        assert_eq!(k[3].0, &ColKind::Id);
        assert_eq!(
            k[4],
            (
                &ColKind::Info {
                    dst: b"AF".to_vec(),
                    src: b"AF".to_vec()
                },
                true,
                false,
                4
            )
        );
        assert_eq!(
            k[5].0,
            &ColKind::Info {
                dst: b"NEW".to_vec(),
                src: b"OLD".to_vec()
            }
        );
        assert!(k[6].2);
        assert!(parse_columns("-INFO/AF", true).is_err());
        assert!(parse_columns("=AF", true).is_err());
        assert!(parse_columns("FORMAT/DP", true).is_err());
        assert!(parse_columns("ALT", true).is_err());
        assert!(parse_columns("^AF", false).is_err());
        let s = parse_columns("ID,^INFO/DP,^AC", true).unwrap();
        assert_eq!(s.items.len(), 2);
        assert_eq!(s.skip_info, vec![b"DP".to_vec(), b"AC".to_vec()]);
    }

    #[test]
    fn id_format() {
        let f = parse_id_format("+%CHROM:%POS:%REF:%FIRST_ALT").unwrap();
        assert!(f.plus);
        let mut r = Rec::default();
        r.load(b"1\t200\t.\tC\tT,G\t.\t.\tDP=3;DB").unwrap();
        let mut out = Vec::new();
        f.render(&r, &mut out);
        assert_eq!(out, b"1:200:C:T");
        let f = parse_id_format("%INFO/DP-%INFO/DB-%INFO/X-%ID").unwrap();
        out.clear();
        f.render(&r, &mut out);
        assert_eq!(out, b"3-1-.-.");
        assert!(parse_id_format("%POS_%REF").is_err());
        assert!(parse_id_format("%INFO").is_err());
        assert!(parse_id_format("[%GT]").is_err());
    }
}
