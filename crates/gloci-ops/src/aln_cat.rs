//! `gloci cat`: join BAM files by copying their compressed blocks, as `samtools cat` does
//! (D123). It means what `vcf concat --naive` means, and shares its copy loop
//! (`gloci_io::bgzf_copy`).
//!
//! The header is the first file's, plus the `@RG` lines of the later files whose ID the
//! header does not have yet (appended at the end, as samtools does), plus gloci's `@PG` line
//! (D104). Unlike samtools, the `@SQ` lists of all files must be equal: records name their
//! reference by number, so a file with another list would silently point at the wrong
//! contigs (`#alignment-cat-sq`). Of each file, the data after its header is copied: the
//! block that holds the header's end is compressed again from there and every later block is
//! copied as it is. The files' EOF markers are dropped and one is written at the end. Records
//! are never parsed and their order is not checked.
//!
//! Each header is read once before anything is written (the output header needs all of them),
//! and again as part of the copy. A regular file is opened twice; stdin or a pipe keeps the
//! blocks its header came from, so it is read once.

use gloci_io::bgzf_copy::{HeadReader, TrackedWrite, copy_from, deflate_blocks, malformed_reason};
use gloci_io::{
    AlnHeader, BGZF_EOF, Compression, Format, IoError, Output, UsageError, WriterOptions,
    detect_compression, detect_format, is_regular_file,
};
use std::fs::File;
use std::io::{self, Cursor, Read, Write};

fn field<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    line.split(|&b| b == b'\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key))
}

fn show(name: &[u8]) -> String {
    String::from_utf8_lossy(name).into_owned()
}

/// The first way the reference list `other` (of `path`) differs from `first` (of `first_path`).
pub(crate) fn sq_difference(
    first_path: &str,
    first: &AlnHeader,
    path: &str,
    other: &AlnHeader,
) -> Option<String> {
    let same = |a: &(Vec<u8>, u64), b: &(Vec<u8>, u64)| a == b;
    for (i, (a, b)) in first.refs.iter().zip(&other.refs).enumerate() {
        if !same(a, b) {
            let what = if a.0 != b.0 {
                format!(
                    "reference {} is {} in {first_path} but {} in {path}",
                    i + 1,
                    show(&a.0),
                    show(&b.0)
                )
            } else {
                format!(
                    "reference {} ({}) is {} long in {first_path} but {} in {path}",
                    i + 1,
                    show(&a.0),
                    a.1,
                    b.1
                )
            };
            return Some(what);
        }
    }
    if first.refs.len() != other.refs.len() {
        let (n, m) = (first.refs.len(), other.refs.len());
        return Some(format!(
            "{first_path} has {n} reference{} and {path} has {m}",
            if n == 1 { "" } else { "s" }
        ));
    }
    None
}

/// The header of the joined file, from the headers of the inputs, as `(name, header)` in
/// order: the first one, plus the `@RG` lines of the others with a new ID. The `@SQ` lists
/// must all equal the first's (same names and lengths, in the same order), or the result is
/// an error naming the first difference (`#alignment-cat-sq`).
pub fn cat_header(inputs: &[(&str, &AlnHeader)]) -> Result<AlnHeader, IoError> {
    let Some(&(first_name, first)) = inputs.first() else {
        return Err(UsageError("cat needs at least one BAM file".into()).into());
    };
    let mut out = first.clone();
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut read_groups = |text: &[u8], add: &mut dyn FnMut(&[u8])| {
        for line in text.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if !line.starts_with(b"@RG\t") {
                continue;
            }
            if let Some(id) = field(line, b"ID:")
                && !ids.iter().any(|i| i == id)
            {
                ids.push(id.to_vec());
                add(line);
            }
        }
    };
    read_groups(&first.text, &mut |_| {});
    for &(name, h) in &inputs[1..] {
        if let Some(diff) = sq_difference(first_name, first, name, h) {
            return Err(IoError::Detect {
                path: name.to_string(),
                msg: format!(
                    "{name}: the @SQ list differs from that of {first_name}: {diff}. Records name their \
reference by number, so joining these files would attach reads to the wrong contigs \
(samtools cat does not check this). Join files aligned to one reference, in one @SQ order \
(see docs/differences.md#alignment-cat-sq)"
                ),
            });
        }
        let mut extra: Vec<u8> = Vec::new();
        read_groups(&h.text, &mut |line| {
            extra.extend_from_slice(line);
            extra.push(b'\n');
        });
        if !extra.is_empty() {
            if !out.text.is_empty() && out.text.last() != Some(&b'\n') {
                out.text.push(b'\n');
            }
            out.text.extend_from_slice(&extra);
        }
    }
    Ok(out)
}

/// One input after its header was read.
struct Entry {
    path: String,
    header: AlnHeader,
    /// Bytes of data the header takes (the magic, text and reference list).
    header_len: u64,
    /// Where the input is read again from: `None` for a regular file, which is reopened.
    held: Option<(Vec<u8>, Box<dyn Read>)>,
}

fn detect(path: &str, msg: String) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg,
    }
}

/// Open `path` as BAM and read its header.
fn read_header(path: &str) -> Result<Entry, IoError> {
    let open_err = |source| IoError::Open {
        path: path.to_string(),
        source,
    };
    let read_err = |source| IoError::Read {
        path: path.to_string(),
        source,
    };
    let mut src: Box<dyn Read> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(open_err)?)
    };
    let name = if path == "-" { "stdin" } else { path };
    let mut head = vec![0u8; 1 << 16];
    let mut got = 0;
    while got < head.len() {
        match src.read(&mut head[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(read_err(e)),
        }
    }
    head.truncate(got);
    if head.is_empty() {
        return Err(IoError::empty_input(name));
    }
    if head.starts_with(b"CRAM") {
        return Err(UsageError(format!(
            "{name} is CRAM; gloci cat joins BAM files only. Use samtools cat for CRAM, or convert \
it first: gloci view {name} -o x.bam"
        ))
        .into());
    }
    match detect_compression(&head) {
        Compression::Bgzf => {}
        Compression::Gzip => {
            return Err(detect(
                path,
                format!(
                    "{name} is gzip but not BGZF, so it is not a BAM file; gloci cat joins BAM files"
                ),
            ));
        }
        Compression::None => {
            return Err(match detect_format(&head) {
                Ok(Format::Sam) => UsageError(format!(
                    "{name} is SAM text; gloci cat joins BAM files byte by byte. \
Convert it first: gloci view {name} -o x.bam"
                ))
                .into(),
                Ok(f) => detect(
                    path,
                    format!("{name} is {}, not BAM; gloci cat joins BAM files", f.name()),
                ),
                Err(_) => detect(
                    path,
                    format!("{name} is not a BAM file; gloci cat joins BAM files"),
                ),
            });
        }
    }
    let mut blocks = HeadReader::new(Cursor::new(head).chain(src));
    let mut magic = [0u8; 4];
    let bad = |e: io::Error| match malformed_reason(&e) {
        Some(m) => detect(
            path,
            format!("{name}: {m}; check that the file is a complete BAM file"),
        ),
        None => read_err(e),
    };
    let mut got = 0;
    while got < magic.len() {
        match blocks.read(&mut magic[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(bad(e)),
        }
    }
    match got {
        // BGZF holding no bytes (only the end-of-file block).
        0 => return Err(IoError::empty_input(name)),
        1..=3 => {
            return Err(detect(
                path,
                format!("{name} is not a BAM file; gloci cat joins BAM files"),
            ));
        }
        _ => {}
    }
    if &magic != b"BAM\x01" {
        return Err(detect(
            path,
            format!(
                "{name} is BGZF-compressed but not BAM (it does not start with the BAM magic bytes); \
gloci cat joins BAM files"
            ),
        ));
    }
    let header = AlnHeader::read_bam(&mut Cursor::new(magic).chain(&mut blocks)).map_err(|e| {
        detect(
            path,
            format!("{name}: {e}; check that the file is a complete BAM file"),
        )
    })?;
    let (stored, rest, header_len) = blocks.finish();
    let held = if is_regular_file(path) {
        None
    } else {
        Some((stored, Box::new(rest) as Box<dyn Read>))
    };
    Ok(Entry {
        path: path.to_string(),
        header,
        header_len,
        held,
    })
}

/// Join the BAM files `paths` (`-` is stdin) into `out`, which must be a raw output
/// (`gloci_io::open_raw_output`): the blocks are written as they are. Finishes `out`. The
/// output's header gets gloci's `@PG` line when `o.command_line` is set.
pub fn cat_bams(paths: &[String], mut out: Output, o: &WriterOptions) -> Result<(), IoError> {
    let mut entries = paths
        .iter()
        .map(|p| read_header(p))
        .collect::<Result<Vec<_>, _>>()?;
    let merged = {
        let named: Vec<(&str, &AlnHeader)> = entries
            .iter()
            .map(|e| (e.path.as_str(), &e.header))
            .collect();
        cat_header(&named)?
    };
    let header = match &o.command_line {
        Some(cmd) => merged.with_program(cmd, o.version),
        None => merged,
    };
    let werr = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    let mut bytes = Vec::new();
    header.write_bam(&mut bytes).map_err(werr)?;
    out.write_all(&deflate_blocks(&bytes).map_err(werr)?)
        .map_err(werr)?;
    drop(bytes);
    for e in &mut entries {
        let mut src: Box<dyn Read> = match e.held.take() {
            Some((stored, rest)) => Box::new(Cursor::new(stored).chain(rest)),
            None => Box::new(File::open(&e.path).map_err(|source| IoError::Open {
                path: e.path.clone(),
                source,
            })?),
        };
        let mut sink = TrackedWrite::new(&mut out);
        let stats = copy_from(&mut src, &mut sink, e.header_len).map_err(|err| {
            if sink.failed {
                werr(err)
            } else if let Some(m) = malformed_reason(&err) {
                detect(
                    &e.path,
                    format!(
                        "{}: {m}; check that the file is a complete BAM file",
                        e.path
                    ),
                )
            } else {
                IoError::Read {
                    path: e.path.clone(),
                    source: err,
                }
            }
        })?;
        if !stats.eof_marker {
            return Err(IoError::Truncated {
                path: e.path.clone(),
            });
        }
    }
    out.write_all(&BGZF_EOF).map_err(werr)?;
    out.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(text: &str) -> AlnHeader {
        AlnHeader::parse_sam(text.as_bytes().to_vec()).unwrap()
    }

    const SQ: &str = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n";

    #[test]
    fn cat_header_rules() {
        let a = h(&format!(
            "@HD\tVN:1.6\tSO:coordinate\n{SQ}@RG\tID:g1\tSM:s1\n@CO\tfirst\n"
        ));
        let b = h(&format!(
            "@HD\tVN:1.6\tSO:queryname\n{SQ}@RG\tID:g1\tSM:other\n@RG\tID:g2\tSM:s2\tLB:x\n@PG\tID:p\n@CO\tsecond\n"
        ));
        let c = h(&format!("{SQ}@RG\tID:g3\n@RG\tID:g2\tSM:again\n"));
        let out = cat_header(&[("a", &a), ("b", &b), ("c", &c)]).unwrap();
        // the first header, then each new @RG ID once, verbatim, at the end
        assert_eq!(
            String::from_utf8(out.text).unwrap(),
            format!(
                "@HD\tVN:1.6\tSO:coordinate\n{SQ}@RG\tID:g1\tSM:s1\n@CO\tfirst\n\
@RG\tID:g2\tSM:s2\tLB:x\n@RG\tID:g3\n"
            )
        );
        assert_eq!(out.refs, a.refs);
    }

    #[test]
    fn a_different_sq_list_names_the_first_difference() {
        let a = h(SQ);
        let err = |b: &AlnHeader| {
            cat_header(&[("a.bam", &a), ("b.bam", b)])
                .unwrap_err()
                .to_string()
        };
        let swapped = h("@SQ\tSN:chr2\tLN:50\n@SQ\tSN:chr1\tLN:100\n");
        assert!(err(&swapped).contains("reference 1 is chr1 in a.bam but chr2 in b.bam"));
        let longer = h("@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:51\n");
        assert!(err(&longer).contains("reference 2 (chr2) is 50 long in a.bam but 51 in b.bam"));
        let fewer = h("@SQ\tSN:chr1\tLN:100\n");
        assert!(err(&fewer).contains("a.bam has 2 references and b.bam has 1"));
        assert!(err(&fewer).starts_with("b.bam:"));
    }
}
