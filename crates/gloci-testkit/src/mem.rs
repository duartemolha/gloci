//! In-memory record streams: operations take readers, never only paths
//! (spec §8.1), so tests feed them text directly.

use gloci_io::{ReaderOptions, RecordReader, open_reader_from};
use std::io::Cursor;

/// A reader over `text`, named `name` in messages (and in order errors).
/// Panics when the text cannot be opened (format detection or a bad
/// header), since tests pass fixed text.
pub fn mem_reader(name: &str, text: &str, o: &ReaderOptions) -> Box<dyn RecordReader> {
    open_reader_from(name, Box::new(Cursor::new(text.as_bytes().to_vec())), o)
        .unwrap_or_else(|e| panic!("cannot open in-memory input {name}: {e}"))
}
