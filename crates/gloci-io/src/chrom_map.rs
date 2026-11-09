//! `--chrom-map FILE` (D62): a global rename table applied on input.
//!
//! The file has two tab-separated columns, `old new`; blank lines and `#`
//! comment lines are ignored. The readers rename before interning, so
//! everything downstream (order checks, `-g` lookups, indexes, output) sees
//! the new names. A name the table does not list passes through unchanged.

use crate::error::{IoError, UsageError};
use crate::vcf::header_fields;
use std::collections::HashMap;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ChromMap {
    forward: HashMap<Vec<u8>, Vec<u8>>,
    backward: HashMap<Vec<u8>, Vec<u8>>,
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

impl ChromMap {
    /// Read and parse the table in `path`.
    pub fn load(path: &str) -> Result<ChromMap, IoError> {
        let text = std::fs::read(path).map_err(|e| {
            UsageError(format!(
                "cannot read the --chrom-map file {path}: {e}; check the --chrom-map path"
            ))
        })?;
        Ok(ChromMap::parse(path, &text)?)
    }

    /// Parse table text; `path` names the file in errors (`path:LINE: ...`).
    /// Malformed lines (not two non-empty tab-separated fields) and a name
    /// that appears twice in either column are usage errors.
    pub fn parse(path: &str, text: &[u8]) -> Result<ChromMap, UsageError> {
        let mut forward: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut backward: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        // The line each name was first given on, for the duplicate messages.
        let mut old_line: HashMap<&[u8], usize> = HashMap::new();
        let mut new_line: HashMap<&[u8], usize> = HashMap::new();
        for (i, line) in text.split(|&b| b == b'\n').enumerate() {
            let n = i + 1;
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
            if fields.len() != 2 {
                return Err(UsageError(format!(
                    "{path}:{n}: expected 2 tab-separated columns (old name, new name), found {}; \
fix or remove the line",
                    fields.len()
                )));
            }
            let (old, new) = (fields[0], fields[1]);
            if old.is_empty() || new.is_empty() {
                return Err(UsageError(format!(
                    "{path}:{n}: empty chromosome name; each line needs an old name and a new name"
                )));
            }
            if let Some(first) = old_line.insert(old, n) {
                return Err(UsageError(format!(
                    "{path}:{n}: \"{}\" is already renamed on line {first}; list each old name once",
                    lossy(old)
                )));
            }
            if let Some(first) = new_line.insert(new, n) {
                return Err(UsageError(format!(
                    "{path}:{n}: \"{}\" is already the new name on line {first}; two chromosomes \
cannot be renamed to the same name",
                    lossy(new)
                )));
            }
            forward.insert(old.to_vec(), new.to_vec());
            backward.insert(new.to_vec(), old.to_vec());
        }
        Ok(ChromMap { forward, backward })
    }

    /// The mapped name, or `name` itself when the table does not list it.
    pub fn rename<'a>(&'a self, name: &'a [u8]) -> &'a [u8] {
        self.forward.get(name).map_or(name, Vec::as_slice)
    }

    /// The old name a new name was renamed from, if the table lists it.
    pub fn reverse(&self, new: &[u8]) -> Option<&[u8]> {
        self.backward.get(new).map(Vec::as_slice)
    }

    /// The name the file itself uses for the chromosome that is called
    /// `requested` after renaming, or `None` when no file chromosome gets
    /// that name: `requested` is neither a mapped new name nor a name the
    /// table leaves alone.
    pub fn file_name<'a>(&'a self, requested: &'a [u8]) -> Option<&'a [u8]> {
        match self.reverse(requested) {
            Some(old) => Some(old),
            None if self.forward.contains_key(requested) => None,
            None => Some(requested),
        }
    }

    /// Refuse a file whose own names collide under the map: it has both an
    /// old name and its new name, and the new name is not itself renamed
    /// away, so renaming would merge two chromosomes. `file` names the file
    /// in the error.
    pub fn check_file_names<'a>(
        &self,
        file: &str,
        names: impl IntoIterator<Item = &'a [u8]> + Clone,
    ) -> Result<(), UsageError> {
        let present: std::collections::HashSet<&[u8]> = names.clone().into_iter().collect();
        for old in names {
            let Some(new) = self.forward.get(old) else {
                continue;
            };
            if new != old && !self.forward.contains_key(new) && present.contains(&new[..]) {
                return Err(UsageError(format!(
                    "{file} has both \"{o}\" and \"{n}\"; --chrom-map renames \"{o}\" to \"{n}\", \
which would merge them; drop that row of the map or rename \"{n}\" too",
                    o = lossy(old),
                    n = lossy(new)
                )));
            }
        }
        Ok(())
    }

    /// `raw` (VCF header text) with the ID of every `##contig` line renamed.
    /// Two contigs that end up with one name are an error.
    pub fn rename_header(&self, raw: &[u8]) -> Result<Vec<u8>, String> {
        let mut out = Vec::with_capacity(raw.len());
        let mut seen: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        for line in raw.split_inclusive(|&b| b == b'\n') {
            // header_fields wants the line without its terminator
            let trimmed = line.strip_suffix(b"\n").unwrap_or(line);
            let trimmed = trimmed.strip_suffix(b"\r").unwrap_or(trimmed);
            let Some(body) = trimmed.strip_prefix(b"##contig=<") else {
                out.extend_from_slice(line);
                continue;
            };
            let id = header_fields(body)
                .into_iter()
                .find(|(k, _)| *k == b"ID")
                .map(|(_, v)| v);
            let Some(id) = id.filter(|v| !v.is_empty()) else {
                out.extend_from_slice(line);
                continue;
            };
            let new = self.rename(id);
            if let Some(prev) = seen.insert(new.to_vec(), id.to_vec())
                && prev != id
            {
                return Err(format!(
                    "--chrom-map gives the contigs \"{}\" and \"{}\" the same name \"{}\"; \
fix the map or the header",
                    lossy(&prev),
                    lossy(id),
                    lossy(new)
                ));
            }
            // `id` points into `line`; replace exactly those bytes.
            let at = id.as_ptr() as usize - line.as_ptr() as usize;
            out.extend_from_slice(&line[..at]);
            out.extend_from_slice(new);
            out.extend_from_slice(&line[at + id.len()..]);
        }
        Ok(out)
    }
}

/// Renames the chromosome column (the first field) of the text lines a
/// reader appends to a batch. The last name seen is cached, so a run of
/// records on one chromosome costs one comparison each, not a lookup.
#[derive(Debug)]
pub(crate) struct Renamer {
    map: std::sync::Arc<ChromMap>,
    old: Vec<u8>,
    new: Vec<u8>,
    /// `old`/`new` hold a lookup, and it changes the name.
    renames: bool,
    cached: bool,
}

impl Renamer {
    pub(crate) fn new(map: std::sync::Arc<ChromMap>) -> Renamer {
        Renamer {
            map,
            old: Vec::new(),
            new: Vec::new(),
            renames: false,
            cached: false,
        }
    }

    pub(crate) fn map(&self) -> &ChromMap {
        &self.map
    }

    /// Rewrite the first field of the line at `buf[start..]`.
    pub(crate) fn apply(&mut self, buf: &mut Vec<u8>, start: usize) {
        let end = memchr::memchr(b'\t', &buf[start..]).map_or(buf.len(), |i| start + i);
        let name = &buf[start..end];
        if !(self.cached && self.old == name) {
            let new = self.map.rename(name);
            self.renames = new != name;
            self.new.clear();
            self.new.extend_from_slice(new);
            self.old.clear();
            self.old.extend_from_slice(name);
            self.cached = true;
        }
        if self.renames {
            buf.splice(start..end, self.new.iter().copied());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(t: &str) -> Result<ChromMap, UsageError> {
        ChromMap::parse("m.tsv", t.as_bytes())
    }

    #[test]
    fn parses_comments_blanks_and_crlf() {
        let m = parse("# c\r\n\r\n1\tchr1\r\nMT\tchrM\n").unwrap();
        assert_eq!(m.rename(b"1"), b"chr1");
        assert_eq!(m.rename(b"MT"), b"chrM");
        assert_eq!(m.rename(b"X"), b"X");
        assert_eq!(m.reverse(b"chrM"), Some(&b"MT"[..]));
        assert_eq!(m.reverse(b"X"), None);
    }

    #[test]
    fn file_name_translates_back() {
        let m = parse("1\tchr1\nchr2\tz\n").unwrap();
        assert_eq!(m.file_name(b"chr1"), Some(&b"1"[..]));
        assert_eq!(m.file_name(b"chrX"), Some(&b"chrX"[..]));
        // "chr2" is renamed to "z", so no chromosome of the file is "chr2"
        assert_eq!(m.file_name(b"chr2"), None);
    }

    #[test]
    fn errors_name_the_line() {
        let e = parse("1\tchr1\n2\n").unwrap_err().0;
        assert!(
            e.starts_with("m.tsv:2: expected 2 tab-separated columns"),
            "{e}"
        );
        let e = parse("1\tchr1\t3\n").unwrap_err().0;
        assert!(e.contains("m.tsv:1:") && e.contains("found 3"), "{e}");
        let e = parse("1\t\n").unwrap_err().0;
        assert!(e.starts_with("m.tsv:1: empty chromosome name"), "{e}");
        let e = parse("1\ta\n1\tb\n").unwrap_err().0;
        assert!(e.starts_with("m.tsv:2:") && e.contains("line 1"), "{e}");
        let e = parse("1\ta\n2\ta\n").unwrap_err().0;
        assert!(
            e.starts_with("m.tsv:2:") && e.contains("new name on line 1"),
            "{e}"
        );
    }

    #[test]
    fn header_contigs_are_renamed() {
        let m = parse("1\tchr1\n").unwrap();
        let out = m
            .rename_header(b"##x=1\n##contig=<ID=1,length=5>\n##contig=<ID=2>\n#CHROM\n")
            .unwrap();
        assert_eq!(
            out,
            b"##x=1\n##contig=<ID=chr1,length=5>\n##contig=<ID=2>\n#CHROM\n"
        );
        let e = m
            .rename_header(b"##contig=<ID=1>\n##contig=<ID=chr1>\n")
            .unwrap_err();
        assert!(e.contains("same name"), "{e}");
    }

    #[test]
    fn renamer_splices_the_first_field() {
        let m = std::sync::Arc::new(parse("1\tchr1\n").unwrap());
        let mut r = Renamer::new(m);
        let mut buf = b"keep\n".to_vec();
        let s = buf.len();
        buf.extend_from_slice(b"1\t5\t6\n1\t7");
        r.apply(&mut buf, s);
        assert_eq!(buf, b"keep\nchr1\t5\t6\n1\t7");
        let mut b2 = b"2\t1\t2".to_vec();
        r.apply(&mut b2, 0);
        assert_eq!(b2, b"2\t1\t2");
        let mut b3 = b"1".to_vec();
        r.apply(&mut b3, 0);
        assert_eq!(b3, b"chr1");
    }
}
