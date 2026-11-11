//! A VCF text record as editable columns, and INFO edits on its raw text.

/// The columns of one VCF line. Buffers are reused from record to record.
#[derive(Default, Debug, Clone)]
pub struct Rec {
    pub chrom: Vec<u8>,
    pub pos: Vec<u8>,
    pub id: Vec<u8>,
    pub reference: Vec<u8>,
    pub alt: Vec<u8>,
    pub qual: Vec<u8>,
    pub filter: Vec<u8>,
    /// The INFO column; `.` when empty.
    pub info: Vec<u8>,
    /// FORMAT and the samples, with a leading tab; empty without samples.
    pub tail: Vec<u8>,
    /// Whether any column was edited (otherwise the input line is written).
    pub dirty: bool,
}

fn set(buf: &mut Vec<u8>, v: &[u8]) {
    buf.clear();
    buf.extend_from_slice(v);
}

impl Rec {
    /// Split `line` (no terminator) into columns.
    pub fn load(&mut self, line: &[u8]) -> Result<(), String> {
        let mut it = line.splitn(9, |&b| b == b'\t');
        let mut next = |what: &str| -> Result<&[u8], String> {
            it.next().ok_or_else(|| {
                format!("the record has too few tab-separated columns (no {what}); a VCF record needs at least 8")
            })
        };
        set(&mut self.chrom, next("CHROM")?);
        set(&mut self.pos, next("POS")?);
        set(&mut self.id, next("ID")?);
        set(&mut self.reference, next("REF")?);
        set(&mut self.alt, next("ALT")?);
        set(&mut self.qual, next("QUAL")?);
        set(&mut self.filter, next("FILTER")?);
        set(&mut self.info, next("INFO")?);
        self.tail.clear();
        if let Some(rest) = it.next() {
            self.tail.push(b'\t');
            self.tail.extend_from_slice(rest);
        }
        self.dirty = false;
        Ok(())
    }

    /// Append the record as a line (no terminator).
    pub fn write(&self, out: &mut Vec<u8>) {
        for (i, c) in [
            &self.chrom,
            &self.pos,
            &self.id,
            &self.reference,
            &self.alt,
            &self.qual,
            &self.filter,
            &self.info,
        ]
        .iter()
        .enumerate()
        {
            if i > 0 {
                out.push(b'\t');
            }
            if c.is_empty() {
                out.push(b'.');
            } else {
                out.extend_from_slice(c);
            }
        }
        out.extend_from_slice(&self.tail);
    }

    pub fn position(&self) -> Result<u64, String> {
        std::str::from_utf8(&self.pos)
            .ok()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| {
                format!(
                    "POS `{}` is not a number",
                    String::from_utf8_lossy(&self.pos)
                )
            })
    }

    /// The ALT alleles (none for `.`).
    pub fn alts(&self) -> Vec<&[u8]> {
        if self.alt == b"." || self.alt.is_empty() {
            Vec::new()
        } else {
            self.alt.split(|&b| b == b',').collect()
        }
    }
}

/// `(item start, item end, value start)` of `key` in an INFO column. The value
/// start is `None` for a bare key (a flag).
fn find(info: &[u8], key: &[u8]) -> Option<(usize, usize, Option<usize>)> {
    if info == b"." || info.is_empty() {
        return None;
    }
    let mut start = 0;
    loop {
        let end = info[start..]
            .iter()
            .position(|&b| b == b';')
            .map_or(info.len(), |p| start + p);
        let item = &info[start..end];
        let (name, vstart) = match item.iter().position(|&b| b == b'=') {
            Some(eq) => (&item[..eq], Some(start + eq + 1)),
            None => (item, None),
        };
        if name == key {
            return Some((start, end, vstart));
        }
        if end >= info.len() {
            return None;
        }
        start = end + 1;
    }
}

/// The value of `key`: `Some(None)` for a bare key, `None` when absent.
pub fn info_get<'a>(info: &'a [u8], key: &[u8]) -> Option<Option<&'a [u8]>> {
    let (_, end, v) = find(info, key)?;
    Some(v.map(|v| &info[v..end]))
}

/// Set `key` to `val` (a bare key for `None`): in place when present, else
/// appended.
pub fn info_set(info: &mut Vec<u8>, key: &[u8], val: Option<&[u8]>) {
    let mut item = key.to_vec();
    if let Some(v) = val {
        item.push(b'=');
        item.extend_from_slice(v);
    }
    match find(info, key) {
        Some((s, e, _)) => {
            info.splice(s..e, item);
        }
        None => {
            if info == b"." || info.is_empty() {
                *info = item;
            } else {
                info.push(b';');
                info.extend_from_slice(&item);
            }
        }
    }
}

/// Remove `key`; whether it was there.
pub fn info_remove(info: &mut Vec<u8>, key: &[u8]) -> bool {
    let Some((s, e, _)) = find(info, key) else {
        return false;
    };
    if s == 0 && e >= info.len() {
        *info = b".".to_vec();
    } else if e >= info.len() {
        info.truncate(s - 1);
    } else {
        info.drain(s..=e);
    }
    true
}

/// Rewrite FORMAT and the samples in `tail` (with its leading tab), keeping
/// the keys for which `keep` is true. With no key left FORMAT and every sample
/// become `.`, as bcftools writes them. Returns whether anything changed.
pub fn format_retain(tail: &mut Vec<u8>, keep: impl Fn(&[u8]) -> bool) -> bool {
    if tail.is_empty() {
        return false;
    }
    let body = &tail[1..];
    let mut cols = body.split(|&b| b == b'\t');
    let fmt = cols.next().unwrap_or(b".");
    let keys: Vec<&[u8]> = if fmt == b"." {
        Vec::new()
    } else {
        fmt.split(|&b| b == b':').collect()
    };
    let kept: Vec<usize> = (0..keys.len()).filter(|&i| keep(keys[i])).collect();
    if kept.len() == keys.len() {
        return false;
    }
    let mut out = vec![b'\t'];
    if kept.is_empty() {
        out.push(b'.');
        for _ in cols {
            out.extend_from_slice(b"\t.");
        }
    } else {
        for (n, &i) in kept.iter().enumerate() {
            if n > 0 {
                out.push(b':');
            }
            out.extend_from_slice(keys[i]);
        }
        for s in cols {
            out.push(b'\t');
            let fields: Vec<&[u8]> = s.split(|&b| b == b':').collect();
            let mut first = true;
            for &i in &kept {
                if let Some(f) = fields.get(i) {
                    if !first {
                        out.push(b':');
                    }
                    first = false;
                    out.extend_from_slice(f);
                }
            }
        }
    }
    *tail = out;
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_edits() {
        let mut i = b"DP=5;OLD=x;DB".to_vec();
        assert_eq!(info_get(&i, b"OLD"), Some(Some(&b"x"[..])));
        assert_eq!(info_get(&i, b"DB"), Some(None));
        assert_eq!(info_get(&i, b"O"), None);
        info_set(&mut i, b"OLD", Some(b"yy"));
        info_set(&mut i, b"AF", Some(b"0.5"));
        assert_eq!(i, b"DP=5;OLD=yy;DB;AF=0.5");
        assert!(info_remove(&mut i, b"DP"));
        assert!(info_remove(&mut i, b"AF"));
        assert_eq!(i, b"OLD=yy;DB");
        assert!(info_remove(&mut i, b"DB"));
        assert!(info_remove(&mut i, b"OLD"));
        assert_eq!(i, b".");
        assert!(!info_remove(&mut i, b"OLD"));
        info_set(&mut i, b"DB", None);
        assert_eq!(i, b"DB");
    }

    #[test]
    fn format_removal() {
        let mut t = b"\tGT:DP:GQ\t0/1:3:9\t1/1:4".to_vec();
        assert!(format_retain(&mut t, |k| k != b"DP"));
        assert_eq!(t, b"\tGT:GQ\t0/1:9\t1/1");
        assert!(format_retain(&mut t, |_| false));
        assert_eq!(t, b"\t.\t.\t.");
        assert!(!format_retain(&mut t, |_| true));
    }

    #[test]
    fn load_and_write() {
        let mut r = Rec::default();
        let line = b"1\t5\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/1";
        r.load(line).unwrap();
        let mut out = Vec::new();
        r.write(&mut out);
        assert_eq!(out, line);
        assert!(r.load(b"1\t5\t.").is_err());
    }
}
