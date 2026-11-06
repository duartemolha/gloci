use crate::detect::{Compression, Format, detect_compression, detect_format};
use crate::error::IoError;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;

const PEEK: usize = 64 * 1024;

/// The 28-byte BGZF end-of-file block.
const BGZF_EOF: [u8; 28] = [
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
fn read_up_to(r: &mut dyn Read, n: usize) -> std::io::Result<Vec<u8>> {
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

pub fn open_input(path: &str, opts: &InputOptions) -> Result<Input, IoError> {
    if path == "-" {
        return open_input_from("stdin", Box::new(std::io::stdin()), opts);
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

pub fn open_input_from(
    name: &str,
    mut src: Box<dyn Read + Send>,
    opts: &InputOptions,
) -> Result<Input, IoError> {
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let raw = read_up_to(&mut src, PEEK).map_err(read_err)?;
    let compression = detect_compression(&raw);
    let chained = Cursor::new(raw).chain(src);
    let mut decoded: Box<dyn Read + Send> = match compression {
        Compression::None => Box::new(chained),
        Compression::Gzip => Box::new(flate2::read::MultiGzDecoder::new(chained)),
        Compression::Bgzf if opts.threads.get() == 1 => {
            Box::new(noodles_bgzf::io::Reader::new(chained))
        }
        Compression::Bgzf => Box::new(noodles_bgzf::io::MultithreadedReader::with_worker_count(
            opts.threads,
            chained,
        )),
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

impl std::fmt::Debug for Input {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Input")
            .field("name", &self.name)
            .field("format", &self.format)
            .field("compression", &self.compression)
            .finish_non_exhaustive()
    }
}
