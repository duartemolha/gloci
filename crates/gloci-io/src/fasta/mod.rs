//! FASTA: `.fai`/`.gzi` indexes written byte-identical to samtools 1.21 and
//! htslib 1.21 (D108), random access with a bounded chunk size, and a
//! forward stream for commands that read a whole FASTA.
//!
//! Random access never builds an index implicitly: a missing `.fai` (and,
//! for BGZF, `.gzi`) is a usage error that names `gloci index FILE`.

mod fetch;
mod index;
mod stream;

pub use fetch::{FETCH_CHUNK, IndexedFasta};
pub use index::{FaiEntry, build_fai, fai_bytes, fasta_index_bytes, parse_fai, write_fasta_index};
pub use stream::{FastaLine, FastaStream};

use crate::error::{IoError, UsageError};
use std::ffi::OsString;
use std::path::{Path, PathBuf};

/// `PATH.fai`.
pub fn fai_path(fasta: &Path) -> PathBuf {
    with_suffix(fasta, ".fai")
}

/// `PATH.gzi`.
pub fn gzi_path(fasta: &Path) -> PathBuf {
    with_suffix(fasta, ".gzi")
}

fn with_suffix(p: &Path, suffix: &str) -> PathBuf {
    let mut s: OsString = p.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// The `.fai` next to `fasta`, and whether it is older than the FASTA.
/// A missing index is a usage error naming the fix (D108); it is never built
/// implicitly.
pub fn find_fai(fasta: &Path) -> Result<(PathBuf, bool), IoError> {
    let fai = fai_path(fasta);
    let name = fasta.display();
    let fai_meta = match std::fs::metadata(&fai) {
        Ok(m) => m,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Err(UsageError(format!(
                "{name} has no index {}; create it with: gloci index {name} \
(gloci never builds a FASTA index implicitly)",
                fai.display()
            ))
            .into());
        }
        Err(source) => {
            return Err(IoError::Open {
                path: fai.display().to_string(),
                source,
            });
        }
    };
    let stale = matches!(
        (fai_meta.modified(), std::fs::metadata(fasta).and_then(|m| m.modified())),
        (Ok(i), Ok(d)) if i < d
    );
    Ok((fai, stale))
}

/// The usage error for a plain-gzip FASTA, which cannot be read at random.
pub(crate) fn gzip_error(name: &str) -> IoError {
    let bgz = name.strip_suffix(".gz").unwrap_or(name);
    UsageError(format!(
        "{name} is gzip-compressed but not BGZF, so it cannot be read at random; recompress it \
with bgzip: gzip -dc {name} | bgzip -c > {bgz}.bgz && mv {bgz}.bgz {name}, then run: gloci index {name}"
    ))
    .into()
}
