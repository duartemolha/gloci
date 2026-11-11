//! The VCF header as editable lines, for `vcf annotate`: declarations are
//! looked up, removed and appended the way htslib's `bcf_hdr_*` do.

use crate::extract::structured_pairs;

/// How an INFO tag is typed in a header.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ty {
    Flag,
    Int,
    Float,
    /// `String` and `Character`.
    Str,
}

/// An INFO tag's declared `Number`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Number {
    Fixed(usize),
    A,
    R,
    G,
    Var,
}

/// A declared INFO tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InfoDef {
    pub ty: Ty,
    pub number: Number,
}

/// The header lines before `#CHROM`, and the `#CHROM` line.
#[derive(Clone, Debug)]
pub struct HeaderEdit {
    /// Every line but `#CHROM`, each ending in `\n`.
    lines: Vec<Vec<u8>>,
    chrom_line: Vec<u8>,
}

/// `(kind, id)` of a structured line `##KIND=<ID=...>`.
fn structured(line: &[u8]) -> Option<(&[u8], Vec<u8>)> {
    let rest = line.strip_prefix(b"##")?;
    let eq = rest.iter().position(|&b| b == b'=')?;
    if rest.get(eq + 1) != Some(&b'<') {
        return None;
    }
    let kind = &rest[..eq];
    let id = structured_pairs(line)
        .into_iter()
        .find(|(k, _)| *k == b"ID")
        .map(|(_, v)| v.strip_prefix(b"\"").unwrap_or(v))?;
    let id = id.strip_suffix(b"\"").unwrap_or(id);
    Some((kind, id.to_vec()))
}

fn trim_eol(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

impl HeaderEdit {
    pub fn parse(raw: &[u8]) -> HeaderEdit {
        let mut lines = Vec::new();
        let mut chrom_line = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
        for line in raw.split_inclusive(|&b| b == b'\n') {
            if line.is_empty() || trim_eol(line).is_empty() {
                continue;
            }
            if line.starts_with(b"#CHROM") {
                chrom_line = line.to_vec();
                if !chrom_line.ends_with(b"\n") {
                    chrom_line.push(b'\n');
                }
            } else {
                let mut l = line.to_vec();
                if !l.ends_with(b"\n") {
                    l.push(b'\n');
                }
                lines.push(l);
            }
        }
        let mut h = HeaderEdit { lines, chrom_line };
        // htslib always declares PASS: as its second line when the file does not.
        if !h.has("FILTER", b"PASS") {
            let at = usize::from(
                h.lines
                    .first()
                    .is_some_and(|l| l.starts_with(b"##fileformat")),
            );
            h.lines.insert(
                at,
                b"##FILTER=<ID=PASS,Description=\"All filters passed\">\n".to_vec(),
            );
        }
        h
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for l in &self.lines {
            out.extend_from_slice(l);
        }
        out.extend_from_slice(&self.chrom_line);
        out
    }

    /// The line declaring `kind` (`INFO`, `FORMAT`, `FILTER`) `id`.
    pub fn line(&self, kind: &str, id: &[u8]) -> Option<&[u8]> {
        self.lines
            .iter()
            .map(|l| trim_eol(l))
            .find(|l| matches!(structured(l), Some((k, i)) if k == kind.as_bytes() && i == id))
    }

    pub fn has(&self, kind: &str, id: &[u8]) -> bool {
        self.line(kind, id).is_some()
    }

    /// The declared IDs of `kind`, in header order.
    pub fn ids(&self, kind: &str) -> Vec<Vec<u8>> {
        self.lines
            .iter()
            .filter_map(|l| structured(trim_eol(l)))
            .filter(|(k, _)| *k == kind.as_bytes())
            .map(|(_, i)| i)
            .collect()
    }

    pub fn info_def(&self, id: &[u8]) -> Option<InfoDef> {
        let line = self.line("INFO", id)?;
        Some(def_of(line))
    }

    /// Remove the `kind` declaration of `id`; whether there was one.
    pub fn remove(&mut self, kind: &str, id: &[u8]) -> bool {
        let before = self.lines.len();
        self.lines.retain(
            |l| !matches!(structured(trim_eol(l)), Some((k, i)) if k == kind.as_bytes() && i == id),
        );
        self.lines.len() != before
    }

    /// Remove every `kind` declaration (but FORMAT/GT with `keep_gt`).
    pub fn remove_all(&mut self, kind: &str, keep_gt: bool) {
        self.lines.retain(|l| match structured(trim_eol(l)) {
            Some((k, i)) if k == kind.as_bytes() => keep_gt && kind == "FORMAT" && i == b"GT",
            _ => true,
        });
    }

    /// Remove the generic lines `##key=value` (bcftools `-x '##key'`).
    pub fn remove_generic(&mut self, key: &[u8]) {
        self.lines.retain(|l| {
            let Some(rest) = l.strip_prefix(b"##") else {
                return true;
            };
            let Some(eq) = rest.iter().position(|&b| b == b'=') else {
                return true;
            };
            !(&rest[..eq] == key && rest.get(eq + 1) != Some(&b'<'))
        });
    }

    /// Remove the structured lines `##key=<...>` that are not INFO, FILTER,
    /// FORMAT, ALT or contig declarations (bcftools `-x key`).
    pub fn remove_structured_key(&mut self, key: &[u8]) {
        self.lines.retain(|l| {
            let Some(rest) = l.strip_prefix(b"##") else {
                return true;
            };
            let Some(eq) = rest.iter().position(|&b| b == b'=') else {
                return true;
            };
            !(&rest[..eq] == key && rest.get(eq + 1) == Some(&b'<'))
        });
    }

    /// Append a header line unless the same declaration (same kind and ID, or
    /// the very same generic line) is already there, as `bcf_hdr_append`.
    pub fn append(&mut self, line: &[u8]) -> Result<(), String> {
        let line = trim_eol(line);
        if !line.starts_with(b"##") {
            return Err(format!(
                "`{}` is not a header line; header lines start with ##",
                String::from_utf8_lossy(line)
            ));
        }
        let dup = match structured(line) {
            Some((k, i)) => self
                .lines
                .iter()
                .any(|l| matches!(structured(trim_eol(l)), Some((k2, i2)) if k2 == k && i2 == i)),
            None => self.lines.iter().any(|l| trim_eol(l) == line),
        };
        if !dup {
            let mut l = line.to_vec();
            l.push(b'\n');
            self.lines.push(l);
        }
        Ok(())
    }

    /// Add a copy of `src_line` (a declaration from another header) renamed
    /// to `id`, without its `IDX=`, unless `id` is already declared.
    pub fn append_renamed(&mut self, kind: &str, src_line: &[u8], id: &[u8]) {
        if self.has(kind, id) {
            return;
        }
        let line = trim_eol(src_line);
        let Some(open) = line.iter().position(|&b| b == b'<') else {
            return;
        };
        let mut out = line[..=open].to_vec();
        let body = line[open + 1..]
            .strip_suffix(b">")
            .unwrap_or(&line[open + 1..]);
        let mut first = true;
        let (mut start, mut quoted) = (0usize, false);
        for i in 0..=body.len() {
            let at_end = i == body.len();
            if !at_end && body[i] == b'"' {
                quoted = !quoted;
            }
            if at_end || (body[i] == b',' && !quoted) {
                let item = &body[start..i];
                start = i + 1;
                if item.starts_with(b"IDX=") {
                    continue;
                }
                if !first {
                    out.push(b',');
                }
                first = false;
                if item.starts_with(b"ID=") {
                    out.extend_from_slice(b"ID=");
                    out.extend_from_slice(id);
                } else {
                    out.extend_from_slice(item);
                }
            }
        }
        out.push(b'>');
        out.push(b'\n');
        self.lines.push(out);
    }
}

fn def_of(line: &[u8]) -> InfoDef {
    let pairs = structured_pairs(line);
    let get = |k: &[u8]| pairs.iter().find(|(n, _)| *n == k).map(|p| p.1);
    let ty = match get(b"Type") {
        Some(b"Flag") => Ty::Flag,
        Some(b"Integer") => Ty::Int,
        Some(b"Float") => Ty::Float,
        _ => Ty::Str,
    };
    let number = match get(b"Number") {
        Some(b"A") => Number::A,
        Some(b"R") => Number::R,
        Some(b"G") => Number::G,
        Some(b".") => Number::Var,
        Some(n) => std::str::from_utf8(n)
            .ok()
            .and_then(|s| s.parse().ok())
            .map_or(Number::Var, Number::Fixed),
        None => Number::Var,
    };
    InfoDef { ty, number }
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: &[u8] = b"##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"a,b\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #[test]
    fn lookup_remove_append() {
        let mut h = HeaderEdit::parse(H);
        assert!(h.has("FILTER", b"PASS"));
        assert_eq!(HeaderEdit::parse(&h.to_bytes()).to_bytes(), h.to_bytes());
        let d = h.info_def(b"DP").unwrap();
        assert_eq!((d.ty, d.number), (Ty::Int, Number::Fixed(1)));
        h.remove_all("FORMAT", true);
        assert!(h.has("FORMAT", b"GT") && !h.has("FORMAT", b"DP"));
        h.append(b"##INFO=<ID=DP,Number=2,Type=Float,Description=\"x\">")
            .unwrap();
        assert_eq!(h.info_def(b"DP").unwrap().ty, Ty::Int);
        h.append(b"##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">")
            .unwrap();
        assert_eq!(h.info_def(b"AF").unwrap().number, Number::A);
        assert!(h.append(b"INFO=x").is_err());
        h.append_renamed(
            "INFO",
            b"##INFO=<ID=OLD,Number=R,Type=Integer,Description=\"o\",IDX=9>",
            b"NEW",
        );
        assert!(
            h.line("INFO", b"NEW")
                .unwrap()
                .ends_with(b"Description=\"o\">")
        );
    }
}
