//! The start and end of a CRAM file, read without decoding anything: the file definition and the
//! SAM header text (`head`, `samples`, `reheader`, `check --quick`), and the end-of-file
//! container. Only CRAM 3.0 and 3.1 are read, as by every other CRAM reader: a 2.x file is
//! refused by its file definition, with the record readers' message.

use super::cram::{CRAM_EOF, unread_version};
use crate::error::IoError;
use noodles_cram as cram;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// The file definition and SAM header of a CRAM file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CramPrelude {
    pub major: u8,
    pub minor: u8,
    pub file_id: [u8; 20],
    /// The SAM header text, up to its first NUL (the container pads it with NULs).
    pub header_text: Vec<u8>,
    /// The offset of the first data container: the end of the header container.
    pub data_start: u64,
}

/// Refuse a CRAM file of a version other than 3.0 and 3.1 (read from its file definition), with
/// the message every CRAM reader gives. noodles reads every container as CRAM 3, so a 2.x header
/// would only fail as a generic decode error.
pub fn check_cram_version(path: &Path) -> Result<(), IoError> {
    let name = path.display().to_string();
    let mut start = [0u8; 6];
    let n = File::open(path)
        .and_then(|mut f| f.read(&mut start))
        .map_err(|source| IoError::Open {
            path: name.clone(),
            source,
        })?;
    if n == 6
        && &start[..4] == b"CRAM"
        && let Some(e) = unread_version(&name, start[4], start[5])
    {
        return Err(e);
    }
    Ok(())
}

/// Read the file definition and the header container of the CRAM file `path`.
pub fn read_cram_prelude(path: &Path) -> Result<CramPrelude, IoError> {
    check_cram_version(path)?;
    let name = path.display().to_string();
    let file = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut src = BufReader::new(file);
    let bad = |e: io::Error| match e.kind() {
        io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof => IoError::Detect {
            path: name.clone(),
            msg: format!(
                "{name} is not a readable CRAM file ({e}); the start of the file is cut short or \
corrupt"
            ),
        },
        _ => IoError::Read {
            path: name.clone(),
            source: e,
        },
    };
    let mut head = cram::io::Reader::new(&mut src);
    let mut hr = head.header_reader();
    hr.read_magic_number().map_err(bad)?;
    let version = hr.read_format_version().map_err(bad)?;
    let file_id = hr.read_file_id().map_err(bad)?;
    let mut container = hr.container_reader().map_err(bad)?;
    let mut text = Vec::new();
    {
        let mut raw = container.raw_sam_header_reader().map_err(bad)?;
        raw.read_to_end(&mut text).map_err(bad)?;
        raw.discard_to_end().map_err(bad)?;
    }
    container.discard_to_end().map_err(bad)?;
    if let Some(nul) = memchr::memchr(0, &text) {
        text.truncate(nul);
    }
    let data_start = src.stream_position().map_err(|source| IoError::Read {
        path: name.clone(),
        source,
    })?;
    // The container reader stops quietly where the file does, so a cut inside the header
    // container shows as a position at the end of the file: no data or end-of-file container
    // follows.
    let len = src
        .get_ref()
        .metadata()
        .map_err(|source| IoError::Read {
            path: name.clone(),
            source,
        })?
        .len();
    if data_start >= len {
        return Err(IoError::Detect {
            path: name.clone(),
            msg: format!(
                "{name} ends inside its CRAM header container (or has no end-of-file container \
after it); the file is probably truncated"
            ),
        });
    }
    Ok(CramPrelude {
        major: version.major(),
        minor: version.minor(),
        file_id,
        header_text: text,
        data_start,
    })
}

/// Does the CRAM file `path` end with the 38-byte CRAM 3 end-of-file container? Only the tail
/// is read (a 2.x file never gets this far: [`check_cram_version`]).
pub fn cram_eof_present(path: &Path) -> io::Result<bool> {
    let mut f = File::open(path)?;
    let len = f.metadata()?.len();
    if len < CRAM_EOF.len() as u64 {
        return Ok(false);
    }
    f.seek(SeekFrom::Start(len - CRAM_EOF.len() as u64))?;
    let mut buf = [0u8; CRAM_EOF.len()];
    f.read_exact(&mut buf)?;
    Ok(buf == CRAM_EOF)
}

/// Read the container at `offset` (the first after the header: [`CramPrelude::data_start`]) as
/// far as its header and blocks, decoding no slice, so that no reference is needed. A container
/// that is cut short or corrupt is an error ([`IoError::Detect`]); the end-of-file container
/// reads like any other.
pub fn cram_first_container(path: &Path, offset: u64) -> Result<(), IoError> {
    let name = path.display().to_string();
    let file = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut src = BufReader::new(file);
    let read_err = |source| IoError::Read {
        path: name.clone(),
        source,
    };
    src.seek(SeekFrom::Start(offset)).map_err(read_err)?;
    let mut container = cram::io::reader::Container::default();
    cram::io::Reader::new(src)
        .read_container(&mut container)
        .map(|_| ())
        .map_err(|e| match e.kind() {
            io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof => IoError::Detect {
                path: name.clone(),
                msg: format!(
                    "the first data container of {name} is cut short or corrupt ({e}); the file \
is probably truncated"
                ),
            },
            _ => IoError::Read {
                path: name.clone(),
                source: e,
            },
        })
}
