//! Bodies of the cargo-fuzz targets in `fuzz/`, kept here so that stable tests
//! can drive the same code with arbitrary bytes. None of these may panic.

use std::io::{BufRead, Cursor};
use std::num::NonZeroUsize;

use gloci_core::{ChromTable, Region, VcfSpan};

use crate::{
    Format, InputOptions, InvalidPolicy, ReaderOptions, RecordBatch, open_input_from,
    open_reader_from,
};

/// Open `data` as `format` (or autodetect) and read every batch.
/// Errors are fine; panics are not.
pub fn read_all(data: &[u8], format: Option<Format>, invalid: InvalidPolicy) {
    let o = ReaderOptions {
        threads: NonZeroUsize::MIN,
        format,
        vcf_span: VcfSpan::Allele,
        invalid,
    };
    let Ok(mut reader) = open_reader_from("fuzz", Box::new(Cursor::new(data.to_vec())), &o) else {
        return;
    };
    let mut chroms = ChromTable::karyotype();
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        match reader.read_batch(&mut chroms, &mut batch, 1024) {
            Ok(true) => {}
            Ok(false) | Err(_) => return,
        }
    }
}

/// Run `read_all` under both invalid-record policies.
pub fn read_all_both(data: &[u8], format: Option<Format>) {
    read_all(data, format, InvalidPolicy::Skip);
    read_all(data, format, InvalidPolicy::Fail);
}

/// Open `data` as an input (detecting compression and format) and drain it.
pub fn drain_input(data: &[u8]) {
    let o = InputOptions {
        threads: NonZeroUsize::MIN,
        format: None,
    };
    let Ok(mut input) = open_input_from("fuzz", Box::new(Cursor::new(data.to_vec())), &o) else {
        return;
    };
    loop {
        let n = match input.reader.fill_buf() {
            Ok(b) => b.len(),
            Err(_) => return,
        };
        if n == 0 {
            return;
        }
        input.reader.consume(n);
    }
}

/// Parse a region string; a success must satisfy `start <= end`.
pub fn region_parse(data: &[u8]) {
    let Ok(s) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(r) = Region::parse(s, |name| name.starts_with(b"chr")) {
        if let Some(end) = r.end {
            assert!(r.start <= end, "start > end for {s:?}");
        }
    }
}
