pub mod index;
pub mod sort;

use gloci_io::IoError;

/// Libraries report output write errors under the placeholder path
/// `output`; name the real destination, keeping the `io::ErrorKind` so a
/// broken pipe is still recognised.
pub(crate) fn name_output(e: IoError, dest: &str) -> IoError {
    match e {
        IoError::Write { path, source } if path == "output" => IoError::Write {
            path: dest.to_string(),
            source,
        },
        e => e,
    }
}
