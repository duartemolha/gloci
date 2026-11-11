//! The parts of `gloci reheader` (D130) that read and write bytes: the exact length of a BAM
//! header, the `@SQ` compatibility rule, and the header container of a CRAM file.

use super::{AlnHeader, read_cram_prelude};
use crate::detect::{Compression, detect_compression};
use crate::error::{IoError, UsageError};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/// A reader that counts the bytes handed out. It does not buffer, so the count is exactly
/// what the header parser consumed.
struct Counting<'a> {
    inner: &'a mut dyn Read,
    n: u64,
}

impl Read for Counting<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let k = self.inner.read(buf)?;
        self.n += k as u64;
        Ok(k)
    }
}

/// The header of the BAM stream `r` (decompressed) and the length of the header in that
/// stream: the magic, `l_text` bytes of text (NUL padding included) and the reference list.
/// `r` must not read ahead of what is asked of it.
pub fn read_bam_header_with_extent(r: &mut dyn Read) -> Result<(AlnHeader, u64), String> {
    let mut counting = Counting { inner: r, n: 0 };
    let header = AlnHeader::read_bam(&mut counting)?;
    Ok((header, counting.n))
}

/// The `M5:` of each `@SQ` line of `text`, in order.
fn sq_m5s(text: &[u8]) -> Vec<Option<&[u8]>> {
    text.split(|&b| b == b'\n')
        .filter(|l| l.starts_with(b"@SQ\t"))
        .map(|l| {
            let l = l.strip_suffix(b"\r").unwrap_or(l);
            l.split(|&b| b == b'\t')
                .skip(1)
                .find_map(|f| f.strip_prefix(b"M5:"))
        })
        .collect()
}

fn m5_of<'a>(v: &[Option<&'a [u8]>], i: usize) -> Option<&'a [u8]> {
    v.get(i).copied().flatten()
}

fn show(h: &AlnHeader, i: usize, m5: Option<&[u8]>, with_m5: bool) -> String {
    let (name, len) = &h.refs[i];
    let mut s = format!("SN:{} LN:{len}", String::from_utf8_lossy(name));
    if let (true, Some(m5)) = (with_m5, m5) {
        s.push_str(&format!(" M5:{}", String::from_utf8_lossy(m5)));
    }
    s
}

const SQ_ADVICE: &str = "reheader keeps each read's reference by its position in the @SQ list, \
so change only SN: to rename, or use gloci view --chrom-map";

/// [`check_sq_compatible`] with the names of the two headers' files for the message.
pub fn check_sq_compatible_named(
    old: &AlnHeader,
    old_name: &str,
    new: &AlnHeader,
    new_name: &str,
    cram: bool,
) -> Result<(), String> {
    let (old_m5, new_m5) = if cram {
        (sq_m5s(&old.text), sq_m5s(&new.text))
    } else {
        (Vec::new(), Vec::new())
    };
    for i in 0..old.refs.len() {
        let n = i + 1;
        if i >= new.refs.len() {
            return Err(format!(
                "@SQ {n} of {old_name} ({}) has no @SQ {n} in {new_name}, which has only {} @SQ \
line{}; {SQ_ADVICE}. Keep every @SQ line, and add new ones only after the last",
                show(old, i, m5_of(&old_m5, i), false),
                new.refs.len(),
                if new.refs.len() == 1 { "" } else { "s" },
            ));
        }
        let (a, b) = (m5_of(&old_m5, i), m5_of(&new_m5, i));
        let m5_differs = matches!((a, b), (Some(x), Some(y)) if x != y);
        if old.refs[i].1 != new.refs[i].1 || m5_differs {
            return Err(format!(
                "@SQ {n} of {new_name} ({}) does not match @SQ {n} of {old_name} ({}); {SQ_ADVICE}",
                show(new, i, b, m5_differs),
                show(old, i, a, m5_differs),
            ));
        }
    }
    Ok(())
}

/// Can `new` replace `old` as a file's header without changing what its records mean? Records
/// name their reference by its position in the `@SQ` list, so each of `old`'s references must
/// keep its position and `LN`; names may change, and `@SQ` lines may follow the last old one.
/// For CRAM (`cram`), an `M5` present on both sides must be equal.
pub fn check_sq_compatible(old: &AlnHeader, new: &AlnHeader, cram: bool) -> Result<(), String> {
    check_sq_compatible_named(old, "the input", new, "the new header", cram)
}

/// The BAM header of `h` (magic, text, reference list) as BGZF blocks at compression `level`,
/// without an end-of-file marker. The last block ends where the header does, so the data that
/// follows starts on a block edge. Only the compressed bytes are held, not the serialised header.
pub fn bam_header_blocks(h: &AlnHeader, level: u8) -> io::Result<Vec<u8>> {
    use noodles_bgzf::io::writer::{Builder, CompressionLevel};
    let level = CompressionLevel::try_from(level)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    let mut w = Builder::default()
        .set_compression_level(level)
        .build_from_writer(Vec::new());
    h.write_bam(&mut w)?;
    let mut out = w.finish()?;
    debug_assert!(out.ends_with(&crate::input::BGZF_EOF));
    out.truncate(out.len() - crate::input::BGZF_EOF.len());
    Ok(out)
}

/// The header of the BAM file `path` and its length in the decompressed stream, reading only
/// the blocks the header spans.
pub fn read_bam_file_header(path: &Path) -> Result<(AlnHeader, u64), IoError> {
    let name = path.display().to_string();
    let file = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut r = noodles_bgzf::io::Reader::new(BufReader::new(file));
    read_bam_header_with_extent(&mut r).map_err(|e| {
        let msg = if e.starts_with("not a BAM file") {
            format!(
                "{name} is BGZF-compressed but not BAM (it does not start with the BAM magic bytes)"
            )
        } else {
            format!("{name}: {e}; check that the file is a complete BAM file")
        };
        IoError::Detect {
            path: name.clone(),
            msg,
        }
    })
}

/// The header that `gloci reheader` is given as HEADER: SAM header text (the leading lines
/// that start with `@` of a SAM file, or of plain, gzip or BGZF text), the header of a BAM file,
/// or the header of a CRAM file. `-` is standard input (not for CRAM, which is read twice).
pub fn read_new_header(src: &str) -> Result<AlnHeader, IoError> {
    let name = if src == "-" { "stdin" } else { src };
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let detect = |msg: String| IoError::Detect {
        path: name.to_string(),
        msg,
    };
    let raw: Box<dyn Read> = if src == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(src).map_err(|source| IoError::Open {
            path: src.to_string(),
            source,
        })?)
    };
    let mut br = BufReader::with_capacity(1 << 16, raw);
    let peek = br.fill_buf().map_err(read_err)?;
    if peek.starts_with(b"CRAM") {
        if src == "-" {
            return Err(UsageError(
                "HEADER on standard input cannot be CRAM; give the CRAM file's name, or its \
header text: gloci view --header-only x.cram -O sam | gloci reheader - y.cram"
                    .into(),
            )
            .into());
        }
        drop(br);
        let prelude = read_cram_prelude(Path::new(src))?;
        return AlnHeader::parse_sam(prelude.header_text)
            .map_err(|e| detect(format!("{name}: {e}")));
    }
    let mut d: Box<dyn BufRead> = match detect_compression(peek) {
        Compression::None => Box::new(br),
        Compression::Gzip => Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(br))),
        Compression::Bgzf => Box::new(BufReader::new(noodles_bgzf::io::Reader::new(br))),
    };
    let head = d.fill_buf().map_err(read_err)?;
    if head.starts_with(b"BAM\x01") {
        let h = AlnHeader::read_bam(&mut d).map_err(|e| {
            detect(format!(
                "{name}: {e}; check that the file is a complete BAM file"
            ))
        })?;
        return Ok(h.with_sq_lines());
    }
    let mut text = Vec::new();
    loop {
        let start = text.len();
        if d.read_until(b'\n', &mut text).map_err(read_err)? == 0 {
            break;
        }
        if text[start] != b'@' {
            text.truncate(start);
            break;
        }
        if text.last() != Some(&b'\n') {
            text.push(b'\n');
        }
    }
    if text.is_empty() {
        return Err(detect(format!(
            "{name} has no SAM header lines (lines that start with @, such as @SQ SN:chr1 LN:1000); \
HEADER is SAM header text, or a SAM, BAM or CRAM file whose header is taken"
        )));
    }
    AlnHeader::parse_sam(text).map_err(|e| detect(format!("{name}: {e}")))
}

fn itf8(v: u32, out: &mut Vec<u8>) {
    match v {
        0..=0x7f => out.push(v as u8),
        0x80..=0x3fff => out.extend_from_slice(&[0x80 | (v >> 8) as u8, v as u8]),
        0x4000..=0x1f_ffff => {
            out.extend_from_slice(&[0xc0 | (v >> 16) as u8, (v >> 8) as u8, v as u8])
        }
        0x20_0000..=0x0fff_ffff => out.extend_from_slice(&[
            0xe0 | (v >> 24) as u8,
            (v >> 16) as u8,
            (v >> 8) as u8,
            v as u8,
        ]),
        _ => out.extend_from_slice(&[
            0xf0 | (v >> 28) as u8,
            (v >> 20) as u8,
            (v >> 12) as u8,
            (v >> 4) as u8,
            (v & 0x0f) as u8,
        ]),
    }
}

/// Write the CRAM 3.x header container: one raw `FILE_HEADER`
/// block holding the 32-bit length of `text` and `text`, with the CRC32 of the container header
/// and of the block (CRAM 3.x; gloci does not read or write 2.x). The file definition is not written.
pub fn write_cram_header_container(text: &[u8], out: &mut dyn Write) -> io::Result<()> {
    let too_big = || {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "the CRAM header text is larger than 2 GiB",
        )
    };
    let data_len = u32::try_from(text.len() + 4)
        .ok()
        .filter(|&n| n <= i32::MAX as u32)
        .ok_or_else(too_big)?;
    // The block: method 0 (raw), content type 0 (FILE_HEADER), content id 0, sizes, data.
    let mut block = vec![0u8, 0, 0];
    itf8(data_len, &mut block);
    itf8(data_len, &mut block);
    block.extend_from_slice(&(text.len() as u32).to_le_bytes());
    block.extend_from_slice(text);
    let mut crc = flate2::Crc::new();
    crc.update(&block);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    let block_len = u32::try_from(block.len()).map_err(|_| too_big())?;
    // The container header: length, reference id 0, start 0, span 0, no records, record
    // counter 0, bases 0 (a single byte in ITF8 and LTF8 alike), one block, no landmarks.
    let mut head = block_len.to_le_bytes().to_vec();
    head.extend_from_slice(&[0, 0, 0, 0, 0, 0, 1, 0]);
    out.write_all(&head)?;
    let mut crc = flate2::Crc::new();
    crc.update(&head);
    out.write_all(&crc.sum().to_le_bytes())?;
    out.write_all(&block)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(text: &str) -> AlnHeader {
        AlnHeader::parse_sam(text.as_bytes().to_vec()).unwrap()
    }

    #[test]
    fn itf8_widths() {
        let enc = |v| {
            let mut o = Vec::new();
            itf8(v, &mut o);
            o
        };
        assert_eq!(enc(0), [0]);
        assert_eq!(enc(0x7f), [0x7f]);
        assert_eq!(enc(0x80), [0x80, 0x80]);
        assert_eq!(enc(0x3fff), [0xbf, 0xff]);
        assert_eq!(enc(0x4000), [0xc0, 0x40, 0x00]);
        assert_eq!(enc(0x20_0000), [0xe0, 0x20, 0, 0]);
        assert_eq!(enc(0x1000_0000), [0xf1, 0, 0, 0, 0]);
        assert_eq!(enc(0x7fff_ffff), [0xf7, 0xff, 0xff, 0xff, 0x0f]);
    }

    #[test]
    fn sq_compatibility_table() {
        let old = h("@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n@SQ\tSN:chr3\tLN:70\n");
        let ok = |new: &str| check_sq_compatible(&old, &h(new), false).is_ok();
        // renamed
        assert!(ok(
            "@SQ\tSN:1\tLN:100\n@SQ\tSN:2\tLN:50\n@SQ\tSN:3\tLN:70\n"
        ));
        // appended
        assert!(ok(
            "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n@SQ\tSN:chr3\tLN:70\n@SQ\tSN:x\tLN:9\n"
        ));
        // a changed LN
        let e = check_sq_compatible(
            &old,
            &h("@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n@SQ\tSN:chr3\tLN:71\n"),
            false,
        )
        .unwrap_err();
        assert!(e.contains("@SQ 3 of the new header (SN:chr3 LN:71)"), "{e}");
        assert!(e.contains("@SQ 3 of the input (SN:chr3 LN:70)"), "{e}");
        assert!(e.contains("--chrom-map"), "{e}");
        // a swapped order
        let e = check_sq_compatible(
            &old,
            &h("@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr3\tLN:70\n@SQ\tSN:chr2\tLN:50\n"),
            false,
        )
        .unwrap_err();
        assert!(e.contains("@SQ 2 of the new header"), "{e}");
        assert!(e.contains("@SQ 2 of the input"), "{e}");
        // a dropped @SQ
        let e = check_sq_compatible(
            &old,
            &h("@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n"),
            false,
        )
        .unwrap_err();
        assert!(e.contains("@SQ 3 of the input (SN:chr3 LN:70)"), "{e}");
        assert!(e.contains("only 2 @SQ lines"), "{e}");
        // names in the message
        let e =
            check_sq_compatible_named(&old, "x.bam", &h("@SQ\tSN:chr1\tLN:99\n"), "new.sam", false)
                .unwrap_err();
        assert!(e.starts_with("@SQ 1 of new.sam (SN:chr1 LN:99) does not match @SQ 1 of x.bam"));
    }

    #[test]
    fn sq_compatibility_m5_for_cram() {
        let old = h("@SQ\tSN:a\tLN:10\tM5:aaaa\n@SQ\tSN:b\tLN:10\n");
        let new = |m: &str, n: &str| h(&format!("@SQ\tSN:a\tLN:10\t{m}\n@SQ\tSN:b\tLN:10\t{n}\n"));
        assert!(check_sq_compatible(&old, &new("M5:aaaa", "M5:bbbb"), true).is_ok());
        // missing on one side: fine
        assert!(check_sq_compatible(&old, &new("UR:x", "UR:y"), true).is_ok());
        let e = check_sq_compatible(&old, &new("M5:cccc", "UR:y"), true).unwrap_err();
        assert!(e.contains("M5:cccc") && e.contains("M5:aaaa"), "{e}");
        // BAM does not look at M5
        assert!(check_sq_compatible(&old, &new("M5:cccc", "UR:y"), false).is_ok());
    }

    #[test]
    fn extent_includes_nul_padding() {
        let mut h = h("@SQ\tSN:chr1\tLN:100\n");
        let mut padded = h.text.clone();
        padded.extend_from_slice(&[0, 0, 0]);
        let mut bytes = b"BAM\x01".to_vec();
        bytes.extend_from_slice(&(padded.len() as i32).to_le_bytes());
        bytes.extend_from_slice(&padded);
        bytes.extend_from_slice(&1i32.to_le_bytes());
        bytes.extend_from_slice(&5i32.to_le_bytes());
        bytes.extend_from_slice(b"chr1\0");
        bytes.extend_from_slice(&100i32.to_le_bytes());
        let header_len = bytes.len() as u64;
        bytes.extend_from_slice(b"records follow");
        let (got, extent) = read_bam_header_with_extent(&mut bytes.as_slice()).unwrap();
        assert_eq!(extent, header_len);
        h.sort_order = got.sort_order;
        assert_eq!(got.text, h.text);
    }
}
