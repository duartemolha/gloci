//! `gloci reheader`: replace the header of a BAM or CRAM file without decoding its records
//! (D130).
//!
//! BAM: the new header is written as BGZF blocks that end where the header ends, then the old
//! file's blocks follow from the end of its header: the block that held that point is
//! compressed again from there, every later block is copied byte for byte (`gloci_io::bgzf_copy`).
//! CRAM: the file definition, a new header container, and every byte of the old file from its
//! first data container on. Memory is the old and new headers plus one BGZF block.
//!
//! Records name their reference by its position in the `@SQ` list, so the new list must keep the
//! old references in order with equal `LN`; names may change and `@SQ` lines may be added at
//! the end (`#reheader-sq-check`). samtools checks nothing.

use gloci_io::bgzf_copy::{TrackedWrite, copy_from, malformed_reason};
use gloci_io::{
    AlnHeader, BGZF_EOF, Compression, Format, IoError, Output, UsageError, WriterOptions,
    bam_header_blocks, check_sq_compatible_named, cram_eof_present, detect_compression,
    detect_format, read_bam_file_header, read_cram_prelude, read_new_header,
    write_cram_header_container,
};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

/// What kind of file the input is.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Bam,
    Cram,
}

fn detect_err(path: &str, msg: String) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg,
    }
}

/// Is the regular file `path` a BAM or CRAM file? Stdin, pipes, SAM text and any other format
/// detected are usage errors (exit 2); a gzip that is not BGZF, or content no format is
/// detected in, is a data error.
fn input_kind(path: &Path) -> Result<Kind, IoError> {
    let name = path.display().to_string();
    if name == "-" || name.starts_with("/dev/fd/") || name == "/dev/stdin" {
        return Err(UsageError(
            "reheader reads IN twice, so IN must be a regular BAM or CRAM file, not standard \
input; write the stream to a file first: gloci view - -o x.bam"
                .into(),
        )
        .into());
    }
    let mut f = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let meta = f.metadata().map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    if !meta.is_file() {
        return Err(UsageError(format!(
            "{name} is not a regular file; reheader reads IN twice, so give a BAM or CRAM file \
(gloci view - -o x.bam writes a stream to one)"
        ))
        .into());
    }
    let mut head = vec![0u8; 4096];
    let mut got = 0;
    while got < head.len() {
        match f.read(&mut head[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(source) => {
                return Err(IoError::Read { path: name, source });
            }
        }
    }
    head.truncate(got);
    if head.is_empty() {
        return Err(IoError::empty_input(&name));
    }
    if head.starts_with(b"CRAM") {
        return Ok(Kind::Cram);
    }
    match detect_compression(&head) {
        Compression::Bgzf => Ok(Kind::Bam),
        // A gzip header that is not BGZF's may be a damaged BAM: a data error.
        Compression::Gzip => Err(detect_err(
            &name,
            format!(
                "{name} is gzip but not BGZF, so it is not a BAM file; reheader takes BAM or CRAM"
            ),
        )),
        Compression::None => Err(match detect_format(&head) {
            Ok(Format::Sam) => UsageError(format!(
                "{name} is SAM text; reheader takes BAM or CRAM. Edit the header of a SAM file \
as text, or make BAM first: gloci view {name} -o x.bam"
            ))
            .into(),
            Ok(f) => IoError::wrong_format(format!(
                "{name} is {}, not BAM or CRAM; reheader takes BAM or CRAM",
                f.name()
            )),
            Err(_) => detect_err(
                &name,
                format!("{name} is not a BAM or CRAM file; reheader takes BAM or CRAM"),
            ),
        }),
    }
}

/// Is the BAM file `path` ended by the BGZF end-of-file marker?
fn bgzf_eof_present(path: &Path) -> std::io::Result<bool> {
    let mut f = File::open(path)?;
    let len = f.metadata()?.len();
    if len < BGZF_EOF.len() as u64 {
        return Ok(false);
    }
    f.seek(SeekFrom::End(-(BGZF_EOF.len() as i64)))?;
    let mut tail = [0u8; 28];
    f.read_exact(&mut tail)?;
    Ok(tail == BGZF_EOF)
}

/// Replace the header of the BAM or CRAM file `input` with `header_src` (SAM header text, or a
/// SAM, BAM or CRAM file whose header is taken; `-` is stdin), writing the result to `out`: a
/// raw output for BAM (`open_raw_output`), a plain one for CRAM. Finishes `out`. The new header
/// gets gloci's `@PG` line when `o.command_line` is set; `level` is the BGZF level of a BAM's
/// new header (the rest of the file is not compressed again).
pub fn reheader(
    header_src: &str,
    input: &Path,
    mut out: Output,
    o: &WriterOptions,
    level: u8,
) -> Result<(), IoError> {
    let name = input.display().to_string();
    let werr = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    let kind = input_kind(input)?;
    if kind == Kind::Cram {
        gloci_io::check_cram_version(input)?;
    }
    let eof = match kind {
        Kind::Bam => bgzf_eof_present(input),
        Kind::Cram => cram_eof_present(input),
    }
    .map_err(|source| IoError::Read {
        path: name.clone(),
        source,
    })?;
    if !eof {
        return Err(match kind {
            Kind::Bam => IoError::Truncated { path: name },
            Kind::Cram => detect_err(
                &name,
                format!(
                    "{name}: the CRAM end-of-file container is missing; the file is probably truncated"
                ),
            ),
        });
    }
    let new_name = if header_src == "-" {
        "stdin"
    } else {
        header_src
    };
    // The old header, and where its data ends (BAM) or where the first data container starts
    // (CRAM). A CRAM's file definition is kept from this one read of its prelude.
    let (mut new, old_header, resume, cram_def) = match kind {
        Kind::Bam => {
            let (h, extent) = read_bam_file_header(input)?;
            (read_new_header(header_src)?, h, extent, None)
        }
        Kind::Cram => {
            let p = read_cram_prelude(input)?;
            let h = AlnHeader::parse_sam(p.header_text)
                .map_err(|e| detect_err(&name, format!("{name}: {e}")))?;
            let def = (p.major, p.minor, p.file_id);
            (read_new_header(header_src)?, h, p.data_start, Some(def))
        }
    };
    check_sq_compatible_named(&old_header, &name, &new, new_name, kind == Kind::Cram)
        .map_err(|msg| detect_err(new_name, msg))?;
    drop(old_header);
    if let Some(cmd) = &o.command_line {
        // `with_program` copies the whole header; take the reference list out of the copy
        // (`@PG` lines are text), so that only the text exists twice.
        let refs = std::mem::take(&mut new.refs);
        let mut with = new.with_program(cmd, o.version);
        drop(new);
        with.refs = refs;
        new = with;
    }
    match cram_def {
        None => {
            let blocks = bam_header_blocks(&new, level).map_err(werr)?;
            drop(new);
            out.write_all(&blocks).map_err(werr)?;
            drop(blocks);
            let mut src = File::open(input).map_err(|source| IoError::Open {
                path: name.clone(),
                source,
            })?;
            let mut sink = TrackedWrite::new(&mut out);
            let stats = copy_from(&mut src, &mut sink, resume).map_err(|err| {
                if sink.failed {
                    werr(err)
                } else if let Some(m) = malformed_reason(&err) {
                    detect_err(
                        &name,
                        format!("{name}: {m}; check that the file is a complete BAM file"),
                    )
                } else {
                    IoError::Read {
                        path: name.clone(),
                        source: err,
                    }
                }
            })?;
            if !stats.eof_marker {
                return Err(IoError::Truncated { path: name });
            }
            out.write_all(&BGZF_EOF).map_err(werr)?;
        }
        Some((major, minor, file_id)) => {
            let mut src = File::open(input).map_err(|source| IoError::Open {
                path: name.clone(),
                source,
            })?;
            let read_err = |source| IoError::Read {
                path: name.clone(),
                source,
            };
            out.write_all(b"CRAM").map_err(werr)?;
            out.write_all(&[major, minor]).map_err(werr)?;
            out.write_all(&file_id).map_err(werr)?;
            write_cram_header_container(&new.text, &mut out).map_err(werr)?;
            drop(new);
            src.seek(SeekFrom::Start(resume)).map_err(read_err)?;
            let mut buf = vec![0u8; 1 << 16];
            loop {
                let n = match src.read(&mut buf) {
                    Ok(0) => break,
                    Ok(n) => n,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    Err(e) => return Err(read_err(e)),
                };
                out.write_all(&buf[..n]).map_err(werr)?;
            }
        }
    }
    out.finish()
}
