//! The start and end of a CRAM file, read without decoding anything: the file definition and the
//! SAM header text (`head`, `samples`, `reheader`, `check --quick`), and the end-of-file
//! container. Any version noodles reads (2.x and 3.x) is accepted, because these commands only
//! move the header and look at the tail.

use crate::error::IoError;
use noodles_cram as cram;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// The 38 bytes of a CRAM 3.x end-of-file container.
const EOF_V3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// The 30 bytes of a CRAM 2.1 end-of-file container.
const EOF_V2: [u8; 30] = [
    0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
];

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

/// Refuse a CRAM file older than 3.0 (read from its file definition), with a message that says
/// what to do. noodles reads every container as CRAM 3, so a 2.x header would only fail as a
/// generic decode error.
pub fn check_cram_version(path: &Path) -> Result<(), IoError> {
    let name = path.display().to_string();
    let mut start = [0u8; 6];
    let n = File::open(path)
        .and_then(|mut f| f.read(&mut start))
        .map_err(|source| IoError::Open {
            path: name.clone(),
            source,
        })?;
    if n == 6 && &start[..4] == b"CRAM" && start[4] < 3 {
        return Err(IoError::Detect {
            path: name.clone(),
            msg: format!(
                "{name} is CRAM {}.{}, which gloci does not read (it reads CRAM 3.0 and 3.1); \
convert it with samtools view -C -o x3.cram {name}",
                start[4], start[5]
            ),
        });
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

/// Does the CRAM file `path` end with its end-of-file container: 38 bytes for CRAM 3.x, 30 for
/// 2.1? Only the version bytes and the tail are read. CRAM 2.0 is checked like 2.1: any major
/// version 2 uses the 30-byte container (htslib reads 2.0 and 2.1 alike).
pub fn cram_eof_present(path: &Path) -> io::Result<bool> {
    let mut f = File::open(path)?;
    let len = f.metadata()?.len();
    let mut start = [0u8; 6];
    let n = f.read(&mut start)?;
    // Version 2.0 is treated like 2.1 here: both end with the 30-byte container.
    let tail: &[u8] = if n == 6 && start[4] == 2 {
        &EOF_V2
    } else {
        &EOF_V3
    };
    if len < tail.len() as u64 {
        return Ok(false);
    }
    f.seek(SeekFrom::Start(len - tail.len() as u64))?;
    let mut buf = vec![0u8; tail.len()];
    f.read_exact(&mut buf)?;
    Ok(buf == tail)
}
