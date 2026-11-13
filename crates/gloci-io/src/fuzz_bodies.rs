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
        accept: crate::align::Accept {
            alignments: Some(crate::align::AlignOptions::default()),
            ..crate::align::Accept::default()
        },
        chrom_map: None,
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

/// Read `data` as BEDPE under both invalid-record policies, with and without a
/// `--chrom-map` (bit 0 of the config byte): the column rules, the absent-end
/// form and the chromosome rewrite are in this path.
pub fn read_bedpe(data: &[u8]) {
    let (cfg, body) = split_config(data);
    let map = (cfg & 1 == 1)
        .then(|| crate::ChromMap::parse("fuzz-map", b"1\tchr1\nchr2\tc2\n").ok())
        .flatten()
        .map(std::sync::Arc::new);
    for invalid in [InvalidPolicy::Skip, InvalidPolicy::Fail] {
        let o = ReaderOptions {
            accept: crate::align::Accept::default(),
            chrom_map: map.clone(),
            threads: threads(cfg),
            format: None,
            vcf_span: VcfSpan::Allele,
            invalid,
        };
        let Ok(mut r) =
            crate::BedpeReader::open_from("fuzz", Box::new(Cursor::new(body.to_vec())), &o)
        else {
            continue;
        };
        let _ = r.header().len();
        let max = [1, 2, 7, 1024][usize::from(cfg >> 3 & 3)];
        let mut chroms = ChromTable::karyotype();
        let mut batch = crate::BedpeBatch::default();
        loop {
            batch.clear();
            match r.read_batch(&mut chroms, &mut batch, max) {
                Ok(true) => {}
                Ok(false) | Err(_) => break,
            }
            for i in 0..batch.len() {
                let _ = crate::bedpe_name(batch.raw(i));
            }
        }
        let _ = r.diagnostics();
    }
}

/// Read `data` as a BAM under both invalid-record policies (the BAM header, framing
/// and every record check are in this path).
pub fn read_bam(data: &[u8]) {
    read_all_both(data, Some(Format::Bam));
}

/// Read `data` as BAM or SAM through the every-record source and convert each record to the
/// other layout, under both invalid-record policies: the SAM line checks, the SAM-to-BAM
/// encoder and the BAM-to-SAM writer are in this path.
pub fn read_aln_source(data: &[u8], format: Format) {
    for invalid in [InvalidPolicy::Skip, InvalidPolicy::Fail] {
        let (cfg, body) = split_config(data);
        let o = ReaderOptions {
            accept: crate::align::Accept {
                alignments: Some(crate::align::AlignOptions::default()),
                ..crate::align::Accept::default()
            },
            chrom_map: None,
            threads: threads(cfg),
            format: Some(format),
            vcf_span: VcfSpan::Allele,
            invalid,
        };
        let Ok(mut src) =
            crate::open_aln_source_from("fuzz", Box::new(Cursor::new(body.to_vec())), &o)
        else {
            continue;
        };
        let header = src.header().clone();
        let encoder = crate::SamToBam::new(&header);
        let (mut raw, mut out) = (Vec::new(), Vec::new());
        loop {
            raw.clear();
            out.clear();
            match src.next_record(&mut raw) {
                Ok(Some(_)) => {}
                Ok(None) | Err(_) => break,
            }
            if src.layout() == Format::Sam {
                let _ = encoder.encode(&raw, &mut out);
                let _ = crate::sam_line_to_bam(&raw, &header, &mut out);
            } else if let Ok(a) = crate::Aln::parse(&raw, &header) {
                let _ = a.write_sam(&header, &mut out);
            }
        }
    }
}

/// Read `data` as a CRAM with no `--reference` (a file with an embedded reference, or written
/// without one, decodes; any other is an error), through the interval reader and the
/// every-record source, under both invalid-record policies: gloci's framing and slice checks,
/// the decode pool, noodles' codecs and the transcoding to BAM bytes are in this path.
pub fn read_cram(data: &[u8]) {
    read_all_both(data, Some(Format::Cram));
    read_aln_source(data, Format::Cram);
}

/// Index `data` as a FASTA and parse `data` as a `.fai`. An index that was built must
/// serialise to text that parses back to the same rows.
pub fn fasta_index(data: &[u8]) {
    let (_, body) = split_config(data);
    if let Ok(entries) = crate::build_fai(&mut Cursor::new(body), "fuzz") {
        let text = crate::fai_bytes(&entries);
        let back = crate::parse_fai(&text, "fuzz").expect("a built .fai parses");
        assert_eq!(crate::fai_bytes(&back), text, "a built .fai is stable");
    }
    if let Ok(entries) = crate::parse_fai(body, "fuzz") {
        let _ = crate::fai_bytes(&entries);
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
    if let Ok(r) = Region::parse(s, |name| name.starts_with(b"chr"))
        && let Some(end) = r.end
    {
        assert!(r.start <= end, "start > end for {s:?}");
    }
}
