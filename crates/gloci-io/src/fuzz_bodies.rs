//! Bodies of the cargo-fuzz targets in `fuzz/`, kept here so that stable tests
//! can drive the same code with arbitrary bytes. None of these may panic.
//!
//! The first input byte is a config byte: bits 1-2 pick the reader thread
//! count (1..=4), bits 3-4 the batch size; bit 0 is left to the caller (the
//! targets use it to choose BGZF wrapping). The rest is the input.

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
    let (cfg, data) = split_config(data);
    let o = ReaderOptions {
        threads: threads(cfg),
        format,
        vcf_span: VcfSpan::Allele,
        invalid,
    };
    let Ok(mut reader) = open_reader_from("fuzz", Box::new(Cursor::new(data.to_vec())), &o) else {
        return;
    };
    let max = [1, 2, 7, 1024][usize::from(cfg >> 3 & 3)];
    let mut chroms = ChromTable::karyotype();
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        match reader.read_batch(&mut chroms, &mut batch, max) {
            Ok(true) => {}
            Ok(false) | Err(_) => return,
        }
    }
}

fn split_config(data: &[u8]) -> (u8, &[u8]) {
    match data.split_first() {
        Some((&c, rest)) => (c, rest),
        None => (0, data),
    }
}

fn threads(cfg: u8) -> NonZeroUsize {
    NonZeroUsize::new(1 + usize::from(cfg >> 1 & 3)).expect("1..=4")
}

/// Run `read_all` under both invalid-record policies.
pub fn read_all_both(data: &[u8], format: Option<Format>) {
    read_all(data, format, InvalidPolicy::Skip);
    read_all(data, format, InvalidPolicy::Fail);
}

/// Open `data` as an input (detecting compression and format) and drain it.
pub fn drain_input(data: &[u8]) {
    let (cfg, data) = split_config(data);
    let o = InputOptions {
        threads: threads(cfg),
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
