//! `gloci samples` (D135): the samples of the read groups of alignment files, from their
//! headers alone. Nothing but the header is read.

use gloci_io::{
    Accept, AlignOptions, AlnHeader, InvalidPolicy, IoError, ReaderOptions, find_fai, find_index,
    open_aln_source, parse_fai,
};
use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

/// What `samples` prints besides the sample and the path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamplesOptions {
    /// The `@RG` tag whose values are listed (`SM`).
    pub tag: [u8; 2],
    /// Add a column: `Y` when the file has an index, else `N`.
    pub show_index: bool,
    /// FASTA files to match against the `@SQ` lines; one column names the matching one.
    pub match_reference: Vec<PathBuf>,
    /// Start with a `#`-prefixed line naming the columns.
    pub column_names: bool,
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// The distinct values of `tag` over the `@RG` lines of the header text, in first-seen order.
fn rg_values(text: &[u8], tag: [u8; 2]) -> Vec<Vec<u8>> {
    let key = [tag[0], tag[1], b':'];
    let mut seen: HashSet<&[u8]> = HashSet::new();
    let mut out = Vec::new();
    for line in text.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if !line.starts_with(b"@RG\t") {
            continue;
        }
        let value = line
            .split(|&b| b == b'\t')
            .skip(1)
            .find_map(|f| f.strip_prefix(&key[..]));
        if let Some(v) = value
            && seen.insert(v)
        {
            out.push(v.to_vec());
        }
    }
    out
}

/// A FASTA's `(name, length)` list from its `.fai`, which must exist (D108).
fn fai_refs(fasta: &Path) -> Result<Vec<(Vec<u8>, u64)>, IoError> {
    let (fai, stale) = find_fai(fasta)?;
    if stale {
        log::warn!(
            "the index {} is older than {}; rebuild it with: gloci index {}",
            fai.display(),
            fasta.display(),
            fasta.display()
        );
    }
    let name = fai.display().to_string();
    let text = std::fs::read(&fai).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    Ok(parse_fai(&text, &name)?
        .into_iter()
        .map(|e| (e.name, e.len))
        .collect())
}

/// Write one line per sample per file: `SAMPLE TAB PATH`, then `Y`/`N` with
/// `show_index`, then the matching FASTA or `.` with `match_reference`. A file with no value
/// of the tag gets one line whose sample is `.`. Samples come in the order the header names
/// them; samtools prints them in hash order.
pub fn samples(paths: &[String], o: &SamplesOptions, out: &mut dyn Write) -> Result<(), IoError> {
    // Every FASTA's references are read once, whatever the number of files, and the check
    // for a missing index comes before the first file is opened.
    let fastas: Vec<(&PathBuf, Vec<(Vec<u8>, u64)>)> = o
        .match_reference
        .iter()
        .map(|p| fai_refs(p).map(|r| (p, r)))
        .collect::<Result<_, _>>()?;
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions::default()),
            convert_only: false,
        },
        chrom_map: None,
        threads: 1.try_into().expect("1 is not zero"),
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let mut text = Vec::new();
    if o.column_names {
        text.extend_from_slice(b"#");
        text.extend_from_slice(&o.tag);
        text.extend_from_slice(b"\tPATH");
        if o.show_index {
            text.extend_from_slice(b"\tINDEX");
        }
        if !fastas.is_empty() {
            text.extend_from_slice(b"\tREFERENCE");
        }
        text.push(b'\n');
    }
    for path in paths {
        let src = open_aln_source(path, &ro)?;
        let header: &AlnHeader = src.header();
        let mut values = rg_values(&header.text, o.tag);
        if values.is_empty() {
            values.push(b".".to_vec());
        }
        let indexed = o.show_index && {
            let p = Path::new(path);
            match find_index(p).filter(|_| path != "-" && p.is_file()) {
                Some(info) => {
                    info.warn_if_stale(path);
                    true
                }
                None => false,
            }
        };
        let reference: Option<&Path> = if fastas.is_empty() {
            None
        } else {
            // samtools tries the FASTA given last first.
            fastas
                .iter()
                .enumerate()
                .rev()
                .find(|(_, (_, refs))| *refs == header.refs)
                .map(|(_, (p, _))| p.as_path())
        };
        for v in &values {
            text.extend_from_slice(v);
            text.push(b'\t');
            text.extend_from_slice(path.as_bytes());
            if o.show_index {
                text.extend_from_slice(if indexed { b"\tY" } else { b"\tN" });
            }
            if !fastas.is_empty() {
                text.push(b'\t');
                match reference {
                    Some(p) => text.extend_from_slice(p.display().to_string().as_bytes()),
                    None => text.push(b'.'),
                }
            }
            text.push(b'\n');
        }
    }
    out.write_all(&text).map_err(write_err)?;
    out.flush().map_err(write_err)
}
