//! `gloci head` (D134): the header and the first records of any input, as text.
//!
//! Nothing is added: no `@PG`, no `##glociCommand`. The file is shown as stored. Memory is the
//! header plus one record (one batch of at most `records` records for the text formats).

use gloci_core::ChromTable;
use gloci_io::{
    Aln, BcfText, Header, IoError, Opened, ReaderOptions, RecordBatch, open_reader_or_alignments,
};
use std::io::Write;

/// What `head` prints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HeadOptions {
    /// Header lines to print; `None` is all of them.
    pub headers: Option<u64>,
    /// Records to print, after the header lines.
    pub records: u64,
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// The first `n` lines of `text` (all when `None`), each ended by a newline.
fn write_lines(text: &[u8], n: Option<u64>, out: &mut dyn Write) -> Result<(), IoError> {
    let mut left = n.unwrap_or(u64::MAX);
    for line in text.split_inclusive(|&b| b == b'\n') {
        if left == 0 {
            break;
        }
        left -= 1;
        out.write_all(line).map_err(write_err)?;
        if line.last() != Some(&b'\n') {
            out.write_all(b"\n").map_err(write_err)?;
        }
    }
    Ok(())
}

/// Write the header lines and records of `path` as text. Whatever was read before a failure
/// has been written to `out` when the error is returned. An alignment file prints SAM, a BCF
/// prints VCF; BED, GFF and VCF print their lines as they are. `ro.accept.alignments` must
/// be set (it carries `--reference`).
pub fn head(
    path: &str,
    o: &HeadOptions,
    ro: &ReaderOptions,
    out: &mut dyn Write,
) -> Result<(), IoError> {
    match open_reader_or_alignments(path, ro)? {
        Opened::Alignments(mut src) => {
            // htslib shows the binary references of a header whose text lacks them.
            let header = src.header().with_sq_lines();
            write_lines(&header.text, o.headers, out)?;
            let name = src.name().to_string();
            let sam = src.layout() == gloci_io::Format::Sam;
            let mut buf = Vec::new();
            let mut line = Vec::new();
            for _ in 0..o.records {
                buf.clear();
                if src.next_record(&mut buf)?.is_none() {
                    break;
                }
                if sam {
                    out.write_all(&buf).map_err(write_err)?;
                } else {
                    line.clear();
                    let bad = |msg: String| IoError::Detect {
                        path: name.clone(),
                        msg: format!("{name}: cannot write a BAM record as SAM: {msg}"),
                    };
                    Aln::parse(&buf, &header)
                        .and_then(|a| a.write_sam(&header, &mut line))
                        .map_err(bad)?;
                    out.write_all(&line).map_err(write_err)?;
                }
                out.write_all(b"\n").map_err(write_err)?;
            }
        }
        Opened::Records(mut reader) => {
            match reader.header() {
                Header::Lines(l) => write_lines(l, o.headers, out)?,
                Header::Vcf(h) => write_lines(&h.raw, o.headers, out)?,
                Header::None | Header::Alignment(_) => {}
            }
            let mut bcf = BcfText::for_reader(reader.as_ref())?;
            let mut chroms = ChromTable::first_appearance();
            let mut batch = RecordBatch::default();
            let mut left = o.records;
            while left > 0 {
                batch.clear();
                let want = usize::try_from(left).unwrap_or(usize::MAX).min(1 << 14);
                let more = reader.read_batch(&mut chroms, &mut batch, want)?;
                for (i, m) in batch.recs.iter().enumerate() {
                    if left == 0 {
                        break;
                    }
                    left -= 1;
                    let line = match &mut bcf {
                        Some(b) => b.render(batch.raw(i), m.loc)?,
                        None => batch.raw(i),
                    };
                    out.write_all(line).map_err(write_err)?;
                    out.write_all(b"\n").map_err(write_err)?;
                }
                if !more {
                    break;
                }
            }
        }
    }
    out.flush().map_err(write_err)
}
