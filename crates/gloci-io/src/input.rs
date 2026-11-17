use crate::detect::{Compression, Format, bbi_format, detect_compression, detect_format};
use crate::error::{IoError, UsageError};
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;

const PEEK: usize = 64 * 1024;

/// The 28-byte BGZF end-of-file block.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

#[derive(Debug, Clone)]
pub struct InputOptions {
    pub threads: NonZeroUsize,
    pub format: Option<Format>,
}

pub struct Input {
    pub name: String,
    pub format: Format,
    pub compression: Compression,
    pub reader: Box<dyn BufRead + Send>,
}

/// Read until `buf` is full or EOF; returns bytes read.
pub(crate) fn read_up_to(r: &mut dyn Read, n: usize) -> std::io::Result<Vec<u8>> {
    let mut v = vec![0u8; n];
    let mut got = 0;
    while got < n {
        match r.read(&mut v[got..]) {
            Ok(0) => break,
            Ok(k) => got += k,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    v.truncate(got);
    Ok(v)
}

/// Standard input, for an input named `-`. A run reads it for at most one input: the reader
/// that opens it second is refused with a usage error ([`crate::cancel::claim_stdin`]). Every
/// reader of `-` opens standard input here.
pub fn open_stdin() -> Result<std::io::Stdin, IoError> {
    crate::cancel::claim_stdin()?;
    Ok(std::io::stdin())
}

pub fn open_input(path: &str, opts: &InputOptions) -> Result<Input, IoError> {
    if path == "-" {
        return open_input_from("stdin", Box::new(open_stdin()?), opts);
    }
    if crate::memio::is_mem_path(path) {
        let (display, src) = crate::memio::take_input(path)?;
        return open_input_from(&display, src, opts);
    }
    let open_err = |source| IoError::Open {
        path: path.to_string(),
        source,
    };
    let mut f = File::open(path).map_err(open_err)?;
    let meta = f.metadata().map_err(open_err)?;
    if meta.is_file() {
        let read_err = |source| IoError::Read {
            path: path.to_string(),
            source,
        };
        let head = read_up_to(&mut f, 18).map_err(read_err)?;
        // bigWig and bigBed (D171) are read through their own index, by the reader, which opens
        // the file itself; a descriptor path is not a file it can reopen.
        let bbi = match opts.format {
            Some(f) if f.is_bbi() => Some(f),
            Some(_) => None,
            None => bbi_format(&head),
        };
        if let Some(format) = bbi
            && is_regular_file(path)
        {
            return Ok(Input {
                name: path.to_string(),
                format,
                compression: Compression::None,
                reader: Box::new(std::io::empty()),
            });
        }
        if detect_compression(&head) == Compression::Bgzf {
            let len = meta.len();
            let ok = len >= 28 && {
                f.seek(SeekFrom::Start(len - 28)).map_err(read_err)?;
                read_up_to(&mut f, 28).map_err(read_err)? == BGZF_EOF
            };
            if !ok {
                return Err(IoError::Truncated {
                    path: path.to_string(),
                });
            }
        }
        f.seek(SeekFrom::Start(0)).map_err(read_err)?;
    }
    open_input_from(path, Box::new(f), opts)
}

/// The bytes of `path` as they are, not decompressed or detected: `-` is stdin, an in-memory
/// name (D84) is taken (once), anything else is a file. For commands that copy or inspect
/// bytes themselves (`cat`, `check --quick`).
pub fn open_raw(path: &str) -> Result<Box<dyn Read + Send>, IoError> {
    if path == "-" {
        return Ok(Box::new(crate::cancel::CancelRead::new(open_stdin()?)));
    }
    if crate::memio::is_mem_path(path) {
        let (_, src) = crate::memio::take_input(path)?;
        return Ok(Box::new(crate::cancel::CancelRead::new(src)));
    }
    let f = File::open(path).map_err(|source| IoError::Open {
        path: path.to_string(),
        source,
    })?;
    Ok(Box::new(crate::cancel::CancelRead::new(f)))
}

pub fn open_input_from(
    name: &str,
    src: Box<dyn Read + Send>,
    opts: &InputOptions,
) -> Result<Input, IoError> {
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    // Every source fails once the host cancels the run (D91).
    let mut src: Box<dyn Read + Send> = Box::new(crate::cancel::CancelRead::new(src));
    let raw = read_up_to(&mut src, PEEK).map_err(read_err)?;
    let bbi = match opts.format {
        Some(f) if f.is_bbi() => Some(f),
        Some(_) => None,
        None => bbi_format(&raw),
    };
    if let Some(format) = bbi {
        return Err(IoError::Usage(UsageError(format!(
            "{name} is {}, which gloci reads through its own index and so only from a file, \
not from standard input, a pipe or memory; save it to a file first",
            format.name()
        ))));
    }
    let compression = detect_compression(&raw);
    let chained = Cursor::new(raw).chain(src);
    let mut decoded: Box<dyn Read + Send> = match compression {
        Compression::None => Box::new(chained),
        Compression::Gzip => Box::new(flate2::read::MultiGzDecoder::new(chained)),
        Compression::Bgzf if opts.threads.get() == 1 => {
            Box::new(noodles_bgzf::io::Reader::new(chained))
        }
        Compression::Bgzf => Box::new(
            crate::bgzf::Reader::new(chained, opts.threads).map_err(crate::bgzf::start_error)?,
        ),
    };
    let peek = read_up_to(&mut decoded, PEEK).map_err(read_err)?;
    let format = match opts.format {
        Some(f) => f,
        None => detect_format(&peek).map_err(|msg| IoError::Detect {
            path: name.to_string(),
            msg: msg.replace("{path}", name),
        })?,
    };
    let reader = BufReader::with_capacity(PEEK, Cursor::new(peek).chain(decoded));
    Ok(Input {
        name: name.to_string(),
        format,
        compression,
        reader: Box::new(reader),
    })
}

impl Input {
    /// Whether the input has no bytes after decompression ([`IoError::empty_input`]).
    pub fn is_empty(&mut self) -> Result<bool, IoError> {
        loop {
            match self.reader.fill_buf() {
                Ok(b) => return Ok(b.is_empty()),
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(source) => {
                    return Err(IoError::Read {
                        path: self.name.clone(),
                        source,
                    });
                }
            }
        }
    }
}

impl std::fmt::Debug for Input {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Input")
            .field("name", &self.name)
            .field("format", &self.format)
            .field("compression", &self.compression)
            .finish_non_exhaustive()
    }
}

/// Whether `path` names a regular file, which can be opened again and read
/// from the start. `-` (stdin), pipes, FIFOs and `/dev/fd/N` cannot (even
/// when `/dev/fd/N` is a regular file, on macOS a reopen shares the file
/// offset): reading one twice silently loses data, so they must be read
/// exactly once.
pub fn is_regular_file(path: &str) -> bool {
    path != "-" && !is_fd_path(path) && std::fs::metadata(path).is_ok_and(|m| m.is_file())
}

/// `/dev/fd/N`, `/dev/stdin`, `/proc/self/fd/N` and `/proc/PID/fd/N`.
fn is_fd_path(path: &str) -> bool {
    if path == "/dev/stdin" || path.starts_with("/dev/fd/") {
        return true;
    }
    path.strip_prefix("/proc/")
        .and_then(|r| r.split_once("/fd/"))
        .is_some_and(|(pid, _)| pid == "self" || pid.bytes().all(|b| b.is_ascii_digit()))
}

/// Copies of inputs that cannot be reopened (stdin, pipes, FIFOs), kept in a
/// registered temp directory so a command that must read an input in a later
/// pass reads its source only once. The directory is removed on drop and, via
/// the cleanup registry, on SIGINT/SIGTERM/SIGHUP or a panic.
pub struct InputSpool {
    // Dropped in this order: the directory is removed, then unregistered.
    dir: tempfile::TempDir,
    _guard: crate::atomic::TempGuard,
    next: usize,
}

impl InputSpool {
    /// Create the spool directory in `parent` (the `--tmp-dir`).
    pub fn create(parent: &std::path::Path) -> Result<InputSpool, IoError> {
        let (dir, guard) =
            crate::atomic::create_temp_dir_in(parent, "gloci-spool-").map_err(|source| {
                IoError::Write {
                    path: parent.display().to_string(),
                    source,
                }
            })?;
        Ok(InputSpool {
            dir,
            _guard: guard,
            next: 0,
        })
    }

    /// Copy all of `path` (`-` is stdin) byte for byte into a new spool file
    /// and return the file's path.
    pub fn spool(&mut self, path: &str) -> Result<std::path::PathBuf, IoError> {
        let shown = crate::memio::display_name(path);
        let name = if path == "-" { "stdin" } else { shown.as_str() };
        let dest = self.dir.path().join(format!("input-{:06}", self.next));
        self.next += 1;
        let mut out =
            crate::atomic::create_file_in_temp_dir(&dest).map_err(|source| IoError::Write {
                path: dest.display().to_string(),
                source,
            })?;
        let mut src: Box<dyn Read> = if path == "-" {
            Box::new(crate::cancel::CancelRead::new(open_stdin()?.lock()))
        } else if crate::memio::is_mem_path(path) {
            Box::new(crate::cancel::CancelRead::new(
                crate::memio::take_input(path)?.1,
            ))
        } else {
            Box::new(File::open(path).map_err(|source| IoError::Open {
                path: path.to_string(),
                source,
            })?)
        };
        let mut buf = vec![0u8; 256 << 10];
        loop {
            let n = match src.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(source) => {
                    return Err(IoError::Read {
                        path: name.to_string(),
                        source,
                    });
                }
            };
            std::io::Write::write_all(&mut out, &buf[..n]).map_err(|source| IoError::Write {
                path: dest.display().to_string(),
                source,
            })?;
        }
        Ok(dest)
    }
}

#[cfg(test)]
mod tests {
    use super::is_regular_file;

    #[test]
    fn fd_paths_and_stdin_are_not_reopenable() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("f");
        std::fs::write(&p, b"x").unwrap();
        assert!(is_regular_file(p.to_str().unwrap()));
        for s in [
            "-",
            "/dev/fd/0",
            "/dev/fd/63",
            "/dev/stdin",
            "/proc/self/fd/3",
            "/proc/1234/fd/5",
        ] {
            assert!(!is_regular_file(s), "{s}");
        }
    }
}
