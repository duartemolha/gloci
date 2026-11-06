//! Byte-level fixture builders shared by the test suites.

use std::io::Write;
use std::path::{Path, PathBuf};

/// BGZF-compress `data`; the output ends with the 28-byte EOF block.
pub fn bgzip(data: &[u8]) -> Vec<u8> {
    let mut w = noodles_bgzf::io::Writer::new(Vec::new());
    w.write_all(data).expect("write to Vec");
    w.finish().expect("finish bgzf")
}

/// Plain (non-BGZF) gzip.
pub fn gzip(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).expect("write to Vec");
    e.finish().expect("finish gzip")
}

/// Write `bytes` to `dir/name` and return the path.
pub fn write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, bytes).expect("write fixture");
    p
}
