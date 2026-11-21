//! `gloci dict`: a SAM sequence dictionary from a FASTA (replaces
//! `samtools dict`, D132).
//!
//! The FASTA is streamed; per sequence only a running MD5 and a length are
//! held, so memory does not depend on sequence length.

use gloci_io::{FastaLine, FastaStream, IoError, Loc};
use md5::{Digest, Md5};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

/// What `dict` writes beyond the mandatory `SN`, `LN` and `M5`.
#[derive(Debug, Clone, Default)]
pub struct DictOptions {
    /// Write the `@HD VN:1.0 SO:unsorted` line.
    pub hd_line: bool,
    /// `AS:` value.
    pub assembly: Option<String>,
    /// `SP:` value.
    pub species: Option<String>,
    /// `UR:` value; when absent, `file://` plus the canonical path of the
    /// FASTA, and no `UR` for standard input.
    pub uri: Option<String>,
    /// Add `AN:` by adding or removing `chr` (samtools `-A`).
    pub alias: bool,
    /// Names that get `AH:*` (samtools `-l`).
    pub alt: Option<HashSet<Vec<u8>>>,
}

/// Feed `line` to `state` without spaces, tabs, `\r` or any other byte outside
/// `!`..=`~`, with `a`-`z` uppercased (what samtools hashes); returns how many
/// bytes were kept, which is the sequence length they add.
pub fn normalized_md5(state: &mut Md5, line: &[u8]) -> u64 {
    let mut buf = [0u8; 4096];
    let mut kept = 0u64;
    for chunk in line.chunks(buf.len()) {
        let mut n = 0;
        for &b in chunk {
            if (b'!'..=b'~').contains(&b) {
                buf[n] = b.to_ascii_uppercase();
                n += 1;
            }
        }
        state.update(&buf[..n]);
        kept += n as u64;
    }
    kept
}

/// The names of an alternate-locus list (`samtools dict -l`): the text of each
/// line up to its first tab, after dropping a final `\r`; empty lines name
/// nothing.
pub fn parse_alt_list(text: &[u8]) -> HashSet<Vec<u8>> {
    text.split(|&b| b == b'\n')
        .filter_map(|l| {
            let l = l.strip_suffix(b"\r").unwrap_or(l);
            let key = l.split(|&b| b == b'\t').next().unwrap_or(&[]);
            (!key.is_empty()).then(|| key.to_vec())
        })
        .collect()
}

/// The `AN:` value samtools gives `name`.
fn alias_of(name: &[u8]) -> Vec<u8> {
    let core = name.strip_prefix(b"chr").unwrap_or(name);
    let mut v = if core.len() < name.len() {
        core.to_vec()
    } else {
        [b"chr".as_slice(), name].concat()
    };
    match core {
        b"M" => v.extend_from_slice(b",chrMT,MT"),
        b"MT" => v.extend_from_slice(b",chrM,M"),
        _ => {}
    }
    v
}

struct Cur {
    name: Vec<u8>,
    md5: Md5,
    len: u64,
}

fn write_sq(
    c: Cur,
    o: &DictOptions,
    uri: Option<&str>,
    out: &mut dyn Write,
) -> std::io::Result<()> {
    out.write_all(b"@SQ\tSN:")?;
    out.write_all(&c.name)?;
    write!(out, "\tLN:{}\tM5:", c.len)?;
    for b in c.md5.finalize() {
        write!(out, "{b:02x}")?;
    }
    if o.alt.as_ref().is_some_and(|a| a.contains(&c.name)) {
        out.write_all(b"\tAH:*")?;
    }
    if o.alias {
        out.write_all(b"\tAN:")?;
        out.write_all(&alias_of(&c.name))?;
    }
    if let Some(u) = uri {
        write!(out, "\tUR:{u}")?;
    }
    if let Some(a) = &o.assembly {
        write!(out, "\tAS:{a}")?;
    }
    if let Some(s) = &o.species {
        write!(out, "\tSP:{s}")?;
    }
    out.write_all(b"\n")
}

/// Write the dictionary of `fasta`. `path_for_uri` is the FASTA's path (`None`
/// for standard input), used for the default `UR`.
pub fn dict(
    mut fasta: FastaStream,
    path_for_uri: Option<&Path>,
    o: &DictOptions,
    out: &mut dyn Write,
) -> Result<(), IoError> {
    let werr = |source| IoError::Write {
        path: "output".into(),
        source,
    };
    let uri: Option<String> = match (&o.uri, path_for_uri) {
        (Some(u), _) => Some(u.clone()),
        (None, Some(p)) => {
            let abs = std::fs::canonicalize(p).map_err(|source| IoError::Open {
                path: p.display().to_string(),
                source,
            })?;
            Some(format!("file://{}", abs.display()))
        }
        (None, None) => None,
    };
    if o.hd_line {
        out.write_all(b"@HD\tVN:1.0\tSO:unsorted\n").map_err(werr)?;
    }
    let mut seen: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut cur: Option<Cur> = None;
    loop {
        let name = match fasta.next()? {
            None => break,
            Some(FastaLine::Header(h)) => {
                let end = h
                    .iter()
                    .position(u8::is_ascii_whitespace)
                    .unwrap_or(h.len());
                h[..end].to_vec()
            }
            Some(FastaLine::Seq(s)) => {
                // Text before the first header belongs to no sequence.
                if let Some(c) = cur.as_mut() {
                    c.len += normalized_md5(&mut c.md5, s);
                }
                continue;
            }
        };
        if let Some(c) = cur.take() {
            write_sq(c, o, uri.as_deref(), out).map_err(werr)?;
        }
        let line = fasta.line();
        if let Some(first) = seen.insert(name.clone(), line) {
            return Err(IoError::Parse {
                path: fasta.name().to_string(),
                loc: Loc::Line(line),
                msg: format!(
                    "sequence name '{}' was already used on line {first}; sequence names must be unique \
(samtools dict writes both); rename one of them",
                    String::from_utf8_lossy(&name)
                ),
            });
        }
        cur = Some(Cur {
            name,
            md5: Md5::new(),
            len: 0,
        });
    }
    if let Some(c) = cur.take() {
        write_sq(c, o, uri.as_deref(), out).map_err(werr)?;
    }
    out.flush().map_err(werr)
}
