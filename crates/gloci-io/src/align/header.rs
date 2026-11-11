//! The alignment header shared by BAM, SAM and CRAM: the SAM header text, and the
//! reference list that BAM records name their contig by (D103).

use crate::chrom_map::ChromMap;
use std::collections::{HashMap, HashSet};
use std::io::{self, Read, Write};
use std::sync::OnceLock;

/// What `@HD SO:` says about the order of the records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Queryname,
    Unsorted,
    /// No `@HD` line, no `SO:`, or a value this build does not know.
    Unknown,
}

impl SortOrder {
    fn from_text(text: &[u8]) -> SortOrder {
        for line in text.split(|&b| b == b'\n') {
            if !line.starts_with(b"@HD") {
                continue;
            }
            for f in line.split(|&b| b == b'\t') {
                if let Some(v) = f.strip_prefix(b"SO:") {
                    return match v.strip_suffix(b"\r").unwrap_or(v) {
                        b"coordinate" => SortOrder::Coordinate,
                        b"queryname" => SortOrder::Queryname,
                        b"unsorted" => SortOrder::Unsorted,
                        _ => SortOrder::Unknown,
                    };
                }
            }
            break;
        }
        SortOrder::Unknown
    }
}

/// The order `gloci sort` puts alignments in (D121).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlnSortOrder {
    /// `@HD SO:coordinate`.
    Coordinate,
    /// `-n`: `SO:queryname SS:queryname:natural`.
    Name,
    /// `-N`: `SO:queryname SS:queryname:lexicographical`.
    NameBytes,
}

/// An alignment header.
///
/// `refs` must not be changed after the first [`AlnHeader::ref_id`] call: the name map behind
/// it is built then, and not rebuilt. Equality compares `text`, `refs` and `sort_order` and
/// ignores that map, so a header renamed by `--chrom-map` (whose map is keyed by the names in
/// the file) can compare equal to a freshly parsed one whose map is keyed by the new names.
#[derive(Debug, Clone)]
pub struct AlnHeader {
    /// The SAM header text (`@HD`, `@SQ`, `@RG`, `@PG`, `@CO` lines), each line ended by `\n`.
    pub text: Vec<u8>,
    /// The references in `@SQ` order, as `(name, length)`. Record `refID`s index this list.
    pub refs: Vec<(Vec<u8>, u64)>,
    pub sort_order: SortOrder,
    /// The value of `@HD SS:`, as written (`queryname:natural`), when the line has one.
    pub sub_sort: Option<Vec<u8>>,
    /// Reference id by the name that records carry, built on first use. For a header renamed by
    /// `--chrom-map` these are the original names.
    names: OnceLock<HashMap<Vec<u8>, usize>>,
}

impl PartialEq for AlnHeader {
    fn eq(&self, other: &AlnHeader) -> bool {
        self.text == other.text && self.refs == other.refs && self.sort_order == other.sort_order
    }
}

impl Eq for AlnHeader {}

/// Longest reference name accepted in a binary header.
const MAX_NAME: usize = 1 << 20;

fn read_full(r: &mut dyn Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut got = 0;
    while got < buf.len() {
        match r.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(k) => got += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(got)
}

const TRUNCATED: &str = "BAM header is truncated; the file was probably cut off";

fn i32_le(r: &mut dyn Read) -> Result<i32, String> {
    let mut b = [0u8; 4];
    let got = read_full(r, &mut b).map_err(|e| format!("cannot read the BAM header: {e}"))?;
    if got < 4 {
        return Err(TRUNCATED.into());
    }
    Ok(i32::from_le_bytes(b))
}

/// The value of the `@HD` line's field `key` (such as `SS:`), when there is one.
fn hd_field(text: &[u8], key: &[u8]) -> Option<Vec<u8>> {
    let line = text.split(|&b| b == b'\n').next()?;
    if !line.starts_with(b"@HD") {
        return None;
    }
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    field(line, key).map(<[u8]>::to_vec)
}

/// The value of a `KEY:` field of a tab-separated header line.
fn field<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    line.split(|&b| b == b'\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key))
}

impl AlnHeader {
    /// A header from its parts: the SAM text, the references in `@SQ` order (the ids that
    /// records use), and the sort order.
    pub fn new(text: Vec<u8>, refs: Vec<(Vec<u8>, u64)>, sort_order: SortOrder) -> AlnHeader {
        let sub_sort = hd_field(&text, b"SS:");
        AlnHeader {
            text,
            refs,
            sort_order,
            sub_sort,
            names: OnceLock::new(),
        }
    }

    /// The header as `samtools sort` leaves it for `o` (D121, `sam_hdr_update_hd`): `@HD`'s `SO`
    /// is replaced where it stands, or appended; `-n` and `-N` set `SS` the same way and
    /// coordinate order removes it; `GO` is removed. With no `@HD` line, `@HD VN:1.6 SO:...
    /// [SS:...]` comes first.
    pub fn with_sort_order(&self, o: AlnSortOrder) -> AlnHeader {
        let (so, ss): (&[u8], Option<&[u8]>) = match o {
            AlnSortOrder::Coordinate => (b"coordinate", None),
            AlnSortOrder::Name => (b"queryname", Some(b"queryname:natural")),
            AlnSortOrder::NameBytes => (b"queryname", Some(b"queryname:lexicographical")),
        };
        let order = if o == AlnSortOrder::Coordinate {
            SortOrder::Coordinate
        } else {
            SortOrder::Queryname
        };
        self.with_hd_order(so, ss, order)
    }

    /// The header with `@HD SO:unsorted`, set as [`AlnHeader::with_sort_order`] sets `SO`
    /// (`SS` and `GO` removed): for records in an order no `SO` value describes.
    pub fn with_unsorted_order(&self) -> AlnHeader {
        self.with_hd_order(b"unsorted", None, SortOrder::Unsorted)
    }

    fn with_hd_order(&self, so: &[u8], ss: Option<&[u8]>, order: SortOrder) -> AlnHeader {
        let mut text = Vec::with_capacity(self.text.len() + 40);
        let mut done = false;
        for line in self.text.split_inclusive(|&b| b == b'\n') {
            if done
                || !line.starts_with(b"@HD")
                || !matches!(line.get(3), Some(b'\t' | b'\n') | None)
            {
                text.extend_from_slice(line);
                continue;
            }
            done = true;
            let (body, end) = match line.strip_suffix(b"\n") {
                Some(b) => (b, &b"\n"[..]),
                None => (line, &b""[..]),
            };
            let mut fields: Vec<Vec<u8>> =
                body.split(|&b| b == b'\t').map(<[u8]>::to_vec).collect();
            let mut set = |key: &[u8], val: Option<&[u8]>| {
                let at = fields
                    .iter()
                    .skip(1)
                    .position(|f| f.starts_with(key))
                    .map(|i| i + 1);
                match (at, val) {
                    (Some(i), Some(v)) => fields[i] = [key, v].concat(),
                    (Some(i), None) => {
                        fields.remove(i);
                    }
                    (None, Some(v)) => fields.push([key, v].concat()),
                    (None, None) => {}
                }
            };
            set(b"SO:", Some(so));
            set(b"SS:", ss);
            set(b"GO:", None);
            text.extend_from_slice(&fields.join(&b'\t'));
            text.extend_from_slice(end);
            if end.is_empty() {
                text.push(b'\n');
            }
        }
        if !done {
            let mut hd = b"@HD\tVN:1.6\tSO:".to_vec();
            hd.extend_from_slice(so);
            if let Some(ss) = ss {
                hd.extend_from_slice(b"\tSS:");
                hd.extend_from_slice(ss);
            }
            hd.push(b'\n');
            text.splice(0..0, hd);
        }
        AlnHeader {
            sort_order: order,
            sub_sort: ss.map(<[u8]>::to_vec),
            refs: self.refs.clone(),
            names: self.names.clone(),
            text,
        }
    }

    fn names(&self) -> &HashMap<Vec<u8>, usize> {
        self.names.get_or_init(|| {
            let mut m = HashMap::with_capacity(self.refs.len());
            for (i, (n, _)) in self.refs.iter().enumerate() {
                m.entry(n.clone()).or_insert(i);
            }
            m
        })
    }

    /// The id of the reference that records name `name`, in O(1). After `--chrom-map` renamed
    /// the header, `name` is the name in the file, not the new one.
    pub fn ref_id(&self, name: &[u8]) -> Option<usize> {
        self.names().get(name).copied()
    }

    /// Read the BAM magic, header text and reference list. Errors do not name the
    /// file; the caller adds it.
    pub fn read_bam(r: &mut dyn Read) -> Result<AlnHeader, String> {
        let mut magic = [0u8; 4];
        let got =
            read_full(r, &mut magic).map_err(|e| format!("cannot read the BAM header: {e}"))?;
        if got < 4 || &magic != b"BAM\x01" {
            return Err(
                "not a BAM file (it does not start with the BAM magic bytes); check --format"
                    .into(),
            );
        }
        let l_text = u64::try_from(i32_le(r)?)
            .map_err(|_| "BAM header has a negative text length".to_string())?;
        let mut text = Vec::new();
        r.take(l_text)
            .read_to_end(&mut text)
            .map_err(|e| format!("cannot read the BAM header: {e}"))?;
        if (text.len() as u64) < l_text {
            return Err(TRUNCATED.into());
        }
        if let Some(nul) = memchr::memchr(0, &text) {
            text.truncate(nul);
        }
        let n_ref = usize::try_from(i32_le(r)?)
            .map_err(|_| "BAM header has a negative reference count".to_string())?;
        // Not reserved up front: a corrupt count must not allocate gigabytes.
        let mut refs: Vec<(Vec<u8>, u64)> = Vec::new();
        for i in 0..n_ref {
            let l_name = usize::try_from(i32_le(r)?)
                .ok()
                .filter(|&n| (1..=MAX_NAME).contains(&n))
                .ok_or_else(|| format!("BAM reference {} has an invalid name length", i + 1))?;
            let mut name = Vec::new();
            r.take(l_name as u64)
                .read_to_end(&mut name)
                .map_err(|e| format!("cannot read the BAM header: {e}"))?;
            if name.len() < l_name {
                return Err(TRUNCATED.into());
            }
            if let Some(nul) = memchr::memchr(0, &name) {
                name.truncate(nul);
            }
            if name.is_empty() {
                return Err(format!("BAM reference {} has an empty name", i + 1));
            }
            let l_ref = u64::try_from(i32_le(r)?).map_err(|_| {
                format!(
                    "BAM reference {} has a negative length",
                    String::from_utf8_lossy(&name)
                )
            })?;
            refs.push((name, l_ref));
        }
        let sort_order = SortOrder::from_text(&text);
        Ok(AlnHeader::new(text, refs, sort_order))
    }

    /// The header as SAM text shows it: a BAM header whose text has no `@SQ` lines gets one
    /// per reference of its binary list, at the end of the text, as htslib shows it. BAM output
    /// that gets no `@PG` line keeps the text as stored (as samtools does), so it does not ask
    /// for this; SAM output and BAM output with a `@PG` line do.
    pub fn with_sq_lines(&self) -> AlnHeader {
        let mut h = self.clone();
        if h.text
            .split(|&b| b == b'\n')
            .any(|l| l.starts_with(b"@SQ\t"))
        {
            return h;
        }
        if !h.text.is_empty() && h.text.last() != Some(&b'\n') {
            h.text.push(b'\n');
        }
        for (name, len) in &h.refs {
            h.text.extend_from_slice(b"@SQ\tSN:");
            h.text.extend_from_slice(name);
            h.text.extend_from_slice(format!("\tLN:{len}\n").as_bytes());
        }
        h
    }

    /// Parse SAM header text: the `@SQ` lines give the references.
    pub fn parse_sam(text: Vec<u8>) -> Result<AlnHeader, String> {
        let mut refs = Vec::new();
        let mut names = HashSet::new();
        for (i, line) in text.split(|&b| b == b'\n').enumerate() {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if !line.starts_with(b"@SQ\t") {
                continue;
            }
            let at = i + 1;
            let name = field(line, b"SN:")
                .filter(|n| !n.is_empty())
                .ok_or_else(|| format!("header line {at}: the @SQ line has no SN: field"))?;
            let len = field(line, b"LN:")
                .and_then(|v| std::str::from_utf8(v).ok())
                .and_then(|v| v.parse::<u32>().ok())
                .filter(|&v| v <= i32::MAX as u32)
                .ok_or_else(|| {
                    format!(
                        "header line {at}: the @SQ line for {} has no valid LN: field \
(a whole number up to 2147483647)",
                        String::from_utf8_lossy(name)
                    )
                })?;
            if !names.insert(name) {
                return Err(format!(
                    "header line {at}: reference {} is declared twice; remove one @SQ line",
                    String::from_utf8_lossy(name)
                ));
            }
            refs.push((name.to_vec(), u64::from(len)));
        }
        let sort_order = SortOrder::from_text(&text);
        Ok(AlnHeader::new(text, refs, sort_order))
    }

    /// Write the BAM magic, the header text and the reference list.
    pub fn write_bam(&self, w: &mut dyn Write) -> io::Result<()> {
        let too_big = |what: &str| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("the BAM header's {what} is larger than 2 GiB"),
            )
        };
        w.write_all(b"BAM\x01")?;
        let l_text = i32::try_from(self.text.len()).map_err(|_| too_big("text"))?;
        w.write_all(&l_text.to_le_bytes())?;
        w.write_all(&self.text)?;
        let n_ref = i32::try_from(self.refs.len()).map_err(|_| too_big("reference list"))?;
        w.write_all(&n_ref.to_le_bytes())?;
        for (name, len) in &self.refs {
            let l_name = i32::try_from(name.len() + 1).map_err(|_| too_big("reference name"))?;
            w.write_all(&l_name.to_le_bytes())?;
            w.write_all(name)?;
            w.write_all(&[0])?;
            let l_ref = i32::try_from(*len).map_err(|_| too_big("reference length"))?;
            w.write_all(&l_ref.to_le_bytes())?;
        }
        Ok(())
    }

    /// The header with gloci's `@PG` lines added, as htslib's `sam_hdr_add_pg` adds samtools'
    /// (D104). One line goes in for each end of the `@PG` chains, that is each line that names an
    /// existing line as its `PP` and that no other line names, with that line as its `PP`.
    /// With none, one line follows the last line that no other line names (an isolated program,
    /// or a line whose `PP` is missing), or comes without `PP` when there is none. The ID is
    /// `gloci`, or `gloci.N` for the first N that no line (added ones included) has. Lines go
    /// after the last `@PG` line, or at the end of the header.
    pub fn with_program(&self, command_line: &str, version: &str) -> AlnHeader {
        let mut text = self.text.clone();
        if !text.is_empty() && text.last() != Some(&b'\n') {
            text.push(b'\n');
        }
        // `(id, pp)` of each @PG line that has an ID, and where the last @PG line ends.
        let mut pgs: Vec<(&[u8], Option<&[u8]>)> = Vec::new();
        let mut after = text.len();
        let mut at = 0;
        for line in self.text.split_inclusive(|&b| b == b'\n') {
            at += line.len() + usize::from(!line.ends_with(b"\n"));
            if line.starts_with(b"@PG\t") {
                let line = line.strip_suffix(b"\n").unwrap_or(line);
                let line = line.strip_suffix(b"\r").unwrap_or(line);
                if let Some(id) = field(line, b"ID:") {
                    pgs.push((id, field(line, b"PP:")));
                }
                after = at;
            }
        }
        let find = |id: &[u8]| pgs.iter().position(|(i, _)| *i == id);
        // A line is named when another line (not itself) has it as `PP`.
        let mut named = vec![false; pgs.len()];
        for (i, (_, pp)) in pgs.iter().enumerate() {
            if let Some(k) = pp.and_then(find).filter(|&k| k != i) {
                named[k] = true;
            }
        }
        let mut ends: Vec<Option<&[u8]>> = (0..pgs.len())
            .filter(|&i| {
                // A line that names itself is a root, not a chain end.
                !named[i] && pgs[i].1.is_some_and(|pp| find(pp).is_some_and(|k| k != i))
            })
            .map(|i| Some(pgs[i].0))
            .collect();
        if ends.is_empty() {
            ends.push((0..pgs.len()).rev().find(|&i| !named[i]).map(|i| pgs[i].0));
        }
        let mut used: Vec<Vec<u8>> = pgs.iter().map(|(i, _)| i.to_vec()).collect();
        let mut added = Vec::new();
        for pp in ends {
            let mut id = "gloci".to_string();
            let mut n = 0u64;
            while used.iter().any(|u| u == id.as_bytes()) {
                n += 1;
                id = format!("gloci.{n}");
            }
            added.extend_from_slice(format!("@PG\tID:{id}\tPN:gloci").as_bytes());
            if let Some(pp) = pp {
                added.extend_from_slice(b"\tPP:");
                added.extend_from_slice(pp);
            }
            added.extend_from_slice(format!("\tVN:{version}\tCL:{command_line}\n").as_bytes());
            used.push(id.into_bytes());
        }
        text.splice(after..after, added);
        AlnHeader {
            text,
            refs: self.refs.clone(),
            sort_order: self.sort_order,
            sub_sort: self.sub_sort.clone(),
            names: self.names.clone(),
        }
    }

    /// The header with every reference renamed by `map` (`--chrom-map`, D103): the
    /// reference list and the `SN:` of each `@SQ` line.
    pub fn renamed(&self, map: &ChromMap) -> AlnHeader {
        let refs = self
            .refs
            .iter()
            .map(|(n, l)| (map.rename(n).to_vec(), *l))
            .collect();
        let mut text = Vec::with_capacity(self.text.len());
        for line in self.text.split_inclusive(|&b| b == b'\n') {
            match field(line, b"SN:").filter(|_| line.starts_with(b"@SQ\t")) {
                Some(sn) => {
                    // `sn` runs to the next TAB or the line end.
                    let end = sn
                        .iter()
                        .position(|&b| b == b'\t' || b == b'\n' || b == b'\r');
                    let sn = &sn[..end.unwrap_or(sn.len())];
                    let at = sn.as_ptr() as usize - line.as_ptr() as usize;
                    text.extend_from_slice(&line[..at]);
                    text.extend_from_slice(map.rename(sn));
                    text.extend_from_slice(&line[at + sn.len()..]);
                }
                None => text.extend_from_slice(line),
            }
        }
        // Records still carry the names of the file.
        let names = OnceLock::new();
        let _ = names.set(
            self.names
                .get_or_init(|| {
                    let mut m = HashMap::with_capacity(self.refs.len());
                    for (i, (n, _)) in self.refs.iter().enumerate() {
                        m.entry(n.clone()).or_insert(i);
                    }
                    m
                })
                .clone(),
        );
        AlnHeader {
            text,
            refs,
            sort_order: self.sort_order,
            sub_sort: self.sub_sort.clone(),
            names,
        }
    }
}
