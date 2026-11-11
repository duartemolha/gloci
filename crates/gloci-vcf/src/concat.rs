//! `vcf concat --naive`: join BGZF-compressed VCF files by copying their
//! compressed blocks (D74).
//!
//! A BGZF file is a series of independent gzip members ("blocks"), so the
//! records of several files can be joined without recompressing them: the
//! first file is copied whole, and of every later file only the blocks after
//! its header are copied. The header may end in the middle of a block. That
//! block then holds the end of the header and the first records, and only
//! its record part belongs in the output: it is compressed again as a new
//! block (the only recompression), and every block after it is copied byte
//! for byte. The EOF marker (an empty block) of each file is dropped; the
//! caller writes one at the very end ([`BGZF_EOF`]).
//!
//! The copy never parses records. To let the caller check the order across
//! files it keeps the first record of each file, and the last one, which it
//! finds by inflating blocks backwards from the end of the file until the
//! record's first byte is reached, so the memory used is about the record's
//! size. The last record is only looked for when the caller asks (every file
//! but the final one).

use gloci_io::bgzf_copy::{
    Block, BlockInfo, TrackedWrite, copy_from_observed, inflate, malformed_reason, read_block,
};
use std::fmt;
use std::io::{self, Read, Seek, SeekFrom, Write};

pub use gloci_io::BGZF_EOF;

fn read_err(e: io::Error) -> NaiveError {
    match malformed_reason(&e) {
        Some(m) => NaiveError::Malformed(m.to_string()),
        None => NaiveError::Read(e),
    }
}

/// Why a block copy failed.
#[derive(Debug)]
pub enum NaiveError {
    /// Reading the input failed.
    Read(io::Error),
    /// Writing the output failed.
    Write(io::Error),
    /// The input is not a well-formed BGZF stream.
    Malformed(String),
    /// The input has no `#CHROM` line, so its header has no end.
    NoHeader,
    /// The header differs from the first file's; holds this file's header.
    HeaderDiffers(Vec<u8>),
}

impl fmt::Display for NaiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NaiveError::Read(e) | NaiveError::Write(e) => write!(f, "{e}"),
            NaiveError::Malformed(m) => f.write_str(m),
            NaiveError::NoHeader => f.write_str("the header has no #CHROM line"),
            NaiveError::HeaderDiffers(_) => f.write_str("the header differs from the first file's"),
        }
    }
}

impl std::error::Error for NaiveError {}

/// What the block copy learned about one file.
#[derive(Debug, Default)]
pub struct NaiveFile {
    /// The decompressed header, up to and including the `#CHROM` line.
    pub header: Vec<u8>,
    /// The first record line, without its newline; `None` without records.
    pub first_record: Option<Vec<u8>>,
    /// The last record line, without its newline; `None` without records.
    pub last_record: Option<Vec<u8>>,
    /// The file ended with the BGZF EOF marker.
    pub eof_marker: bool,
    /// The data ends with a newline (true for a file without records).
    pub ends_with_newline: bool,
}

/// The length of the header in `text`: up to and including the line that
/// starts with `#CHROM`; `None` while that line is incomplete.
fn header_len(text: &[u8]) -> Option<usize> {
    let mut start = 0;
    while start < text.len() {
        let nl = text[start..].iter().position(|&b| b == b'\n')?;
        if text[start..].starts_with(b"#CHROM") {
            return Some(start + nl + 1);
        }
        start += nl + 1;
    }
    None
}

fn strip_eol(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// Copy the blocks of the BGZF VCF `src` to `out`, as the module comment
/// describes (the copy itself is `gloci_io::bgzf_copy`). `want_last` asks for
/// [`NaiveFile::last_record`] (and `ends_with_newline`), which needs reading
/// blocks again. With `expected_header` (every file but the first) the header
/// must equal it and none of the header's blocks is written; without it (the
/// first file) every block is copied.
pub fn copy_blocks(
    src: &mut (impl Read + Seek),
    out: &mut dyn Write,
    expected_header: Option<&[u8]>,
    want_last: bool,
) -> Result<NaiveFile, NaiveError> {
    let mut info = NaiveFile {
        ends_with_newline: true,
        ..NaiveFile::default()
    };
    // First pass: the header, and the first record. Nothing is written yet.
    // `acc` holds the data until the header ends, then the record bytes seen
    // before the first record's newline.
    let mut acc: Vec<u8> = Vec::new();
    let mut header_len_found = None;
    let mut first_done = false;
    while let Some(block) = read_block(src).map_err(read_err)? {
        if block.isize == 0 {
            continue;
        }
        acc.extend_from_slice(&inflate(&block).map_err(read_err)?);
        if header_len_found.is_none()
            && let Some(n) = header_len(&acc)
        {
            header_len_found = Some(n);
            info.header = acc[..n].to_vec();
            if let Some(want) = expected_header
                && want != info.header.as_slice()
            {
                return Err(NaiveError::HeaderDiffers(info.header));
            }
            acc.drain(..n);
        }
        if header_len_found.is_some()
            && let Some(nl) = acc.iter().position(|&b| b == b'\n')
        {
            info.first_record = Some(strip_eol(&acc[..=nl]).to_vec());
            first_done = true;
            break;
        }
    }
    let Some(header_len) = header_len_found else {
        return Err(NaiveError::NoHeader);
    };
    if !first_done && !acc.is_empty() {
        // A last record without a newline.
        info.first_record = Some(strip_eol(&acc).to_vec());
    }
    drop(acc);

    // Second pass: the copy. The first file is copied whole; of the others,
    // everything after the header.
    src.seek(SeekFrom::Start(0)).map_err(NaiveError::Read)?;
    let mut blocks: Vec<(u64, usize, usize)> = Vec::new();
    let mut sink = TrackedWrite::new(out);
    let start = if expected_header.is_some() {
        header_len as u64
    } else {
        0
    };
    let stats = copy_from_observed(src, &mut sink, start, &mut |b: BlockInfo| {
        if want_last && b.isize > 0 {
            blocks.push((b.offset, b.len, b.isize));
        }
    })
    .map_err(|e| {
        if sink.failed {
            NaiveError::Write(e)
        } else {
            read_err(e)
        }
    })?;
    let total = stats.uncompressed as usize;
    info.eof_marker = stats.eof_marker;
    if info.first_record.is_none() || !want_last {
        return Ok(info);
    }

    // The last record: inflate blocks from the end until a newline before
    // the record (or the end of the header) is in the text. Each block is
    // tested once and the text is joined once, so a record of many blocks
    // costs linear time; memory is bounded by the record.
    let mut parts: Vec<Vec<u8>> = Vec::new();
    let mut tail_start = total;
    for &(at, len, isize) in blocks.iter().rev() {
        src.seek(SeekFrom::Start(at)).map_err(NaiveError::Read)?;
        let mut raw = vec![0u8; len];
        src.read_exact(&mut raw).map_err(NaiveError::Read)?;
        let data = inflate(&Block { raw, isize }).map_err(read_err)?;
        tail_start -= data.len();
        // The last block's own line ending does not end an earlier line.
        let has_newline = if parts.is_empty() {
            strip_eol(&data).contains(&b'\n')
        } else {
            data.contains(&b'\n')
        };
        parts.push(data);
        if tail_start <= info.header.len() || has_newline {
            break;
        }
    }
    let mut text: Vec<u8> = Vec::with_capacity(parts.iter().map(Vec::len).sum());
    for part in parts.iter().rev() {
        text.extend_from_slice(part);
    }
    drop(parts);
    info.ends_with_newline = text.ends_with(b"\n");
    let body = strip_eol(&text);
    let header_in_tail = info.header.len().saturating_sub(tail_start);
    let start = match body.iter().rposition(|&b| b == b'\n') {
        Some(p) => (p + 1).max(header_in_tail),
        None => header_in_tail,
    };
    info.last_record = Some(strip_eol(&body[start.min(body.len())..]).to_vec());
    Ok(info)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_io::bgzf_copy::deflate_blocks as deflate;

    fn bgzf(parts: &[&[u8]]) -> Vec<u8> {
        let mut v = Vec::new();
        for p in parts {
            v.extend(deflate(p).unwrap());
        }
        v.extend_from_slice(&BGZF_EOF);
        v
    }

    fn inflate_all(mut data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        while let Some(b) = read_block(&mut data).unwrap() {
            out.extend(inflate(&b).unwrap());
        }
        out
    }

    const H: &[u8] = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #[test]
    fn header_ending_mid_block_recompresses_only_that_block() {
        let a = [H, b"1\t5\t.\tA\tC\t.\t.\t.\n1\t6\t.\tA\tC\t.\t.\t.\n"].concat();
        let b = [H, b"1\t7\t.\tA\tC\t.\t.\t.\n"].concat();
        let (fa, fb) = (bgzf(&[&a]), bgzf(&[&b, b"1\t8\t.\tA\tC\t.\t.\t.\n"]));
        let mut out = Vec::new();
        let i = copy_blocks(
            &mut std::io::Cursor::new(fa.as_slice()),
            &mut out,
            None,
            true,
        )
        .unwrap();
        let j = copy_blocks(
            &mut std::io::Cursor::new(fb.as_slice()),
            &mut out,
            Some(&i.header),
            true,
        )
        .unwrap();
        assert_eq!(j.first_record.unwrap(), b"1\t7\t.\tA\tC\t.\t.\t.");
        assert_eq!(j.last_record.unwrap(), b"1\t8\t.\tA\tC\t.\t.\t.");
        assert!(i.eof_marker && j.eof_marker);
        out.extend_from_slice(&BGZF_EOF);
        let want = [
            a.as_slice(),
            b"1\t7\t.\tA\tC\t.\t.\t.\n1\t8\t.\tA\tC\t.\t.\t.\n",
        ]
        .concat();
        assert_eq!(inflate_all(&out), want);
        // the second block of file two is copied byte for byte
        let second = &fb[fb.len() - BGZF_EOF.len() - 40..fb.len() - BGZF_EOF.len()];
        assert!(out.windows(second.len()).any(|w| w == second));
    }

    #[test]
    fn header_ending_on_a_block_edge_copies_every_later_block_verbatim() {
        let rec = |p: u32| format!("1\t{p}\t.\tA\tC\t.\t.\t.\n").into_bytes();
        let fa = bgzf(&[H, &rec(1)]);
        let (r2, r3) = (rec(2), rec(3));
        let fb = bgzf(&[H, &r2, &r3]);
        let mut out = Vec::new();
        let i = copy_blocks(
            &mut std::io::Cursor::new(fa.as_slice()),
            &mut out,
            None,
            true,
        )
        .unwrap();
        let before = out.len();
        let j = copy_blocks(
            &mut std::io::Cursor::new(fb.as_slice()),
            &mut out,
            Some(&i.header),
            true,
        )
        .unwrap();
        assert_eq!(j.first_record.unwrap(), b"1\t2\t.\tA\tC\t.\t.\t.");
        // nothing was recompressed: the output is file b's blocks after the header block
        let hb = deflate(H).unwrap().len();
        assert_eq!(&out[before..], &fb[hb..fb.len() - BGZF_EOF.len()]);
    }

    #[test]
    fn a_differing_header_is_refused() {
        let a = bgzf(&[H, b"1\t5\t.\tA\tC\t.\t.\t.\n"]);
        let other = [H, b"##x=1\n"].concat();
        let mut out = Vec::new();
        let r = copy_blocks(
            &mut std::io::Cursor::new(a.as_slice()),
            &mut out,
            Some(&other),
            true,
        );
        assert!(matches!(r, Err(NaiveError::HeaderDiffers(_))));
    }

    #[test]
    fn a_missing_eof_marker_is_reported() {
        let mut a = bgzf(&[H, b"1\t5\t.\tA\tC\t.\t.\t.\n"]);
        a.truncate(a.len() - BGZF_EOF.len());
        let mut out = Vec::new();
        let i = copy_blocks(
            &mut std::io::Cursor::new(a.as_slice()),
            &mut out,
            None,
            true,
        )
        .unwrap();
        assert!(!i.eof_marker);
        assert_eq!(i.last_record.unwrap(), b"1\t5\t.\tA\tC\t.\t.\t.");
    }

    #[test]
    fn a_huge_last_record_is_found_and_the_final_file_skips_the_scan() {
        let long = format!("1\t9\t.\tA\tC\t.\t.\t{}\n", "x".repeat(2_500_000));
        let f = bgzf(&[H, b"1\t5\t.\tA\tC\t.\t.\t.\n", long.as_bytes()]);
        let mut out = Vec::new();
        let i = copy_blocks(
            &mut std::io::Cursor::new(f.as_slice()),
            &mut out,
            None,
            true,
        )
        .unwrap();
        assert_eq!(i.last_record.unwrap(), long.trim_end().as_bytes());
        assert!(i.ends_with_newline);
        let mut out = Vec::new();
        let j = copy_blocks(
            &mut std::io::Cursor::new(f.as_slice()),
            &mut out,
            None,
            false,
        )
        .unwrap();
        assert!(j.last_record.is_none());
        assert_eq!(j.first_record.unwrap(), b"1\t5\t.\tA\tC\t.\t.\t.");
    }

    #[test]
    fn a_record_spanning_blocks_is_found() {
        let long = format!("1\t9\t.\tA\tC\t.\t.\t{}\n", "x".repeat(300));
        let f = bgzf(&[
            H,
            b"1\t5\t.\tA\tC\t.\t.\t.\n",
            &long.as_bytes()[..100],
            &long.as_bytes()[100..],
        ]);
        let mut out = Vec::new();
        let i = copy_blocks(
            &mut std::io::Cursor::new(f.as_slice()),
            &mut out,
            None,
            true,
        )
        .unwrap();
        assert_eq!(i.last_record.unwrap(), long.trim_end().as_bytes());
    }
}
