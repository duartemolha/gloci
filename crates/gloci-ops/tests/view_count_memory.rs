//! `view --count` over a scan holds one record and the reader's pool, however many records the
//! file has (D133). A counting allocator measures the peak of the count over a small and a large
//! generated file, SAM and BAM, with and without a regions file.
//!
//! The measurement is in-process on purpose. A child's `ru_maxrss` from `wait4` is not usable
//! here: the child shares the parent's address space until exec and the kernel keeps the parent's
//! high-water mark in the child's figure, so whatever the test process held (a generated input,
//! say) is a floor under every reading, and a buffering regression would hide beneath it.
//! The input is written to disk record by record, never held in memory.
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, OutputOptions, ReaderOptions, WriterOptions,
    plain_text_output_spec,
};
use gloci_ops::{ViewOptions, run_view};
use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

struct Counting;
static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(l) };
        if !p.is_null() {
            let now = LIVE.fetch_add(l.size(), Relaxed) + l.size();
            PEAK.fetch_max(now, Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) };
        LIVE.fetch_sub(l.size(), Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let q = unsafe { System.realloc(p, l, new) };
        if !q.is_null() {
            let now = if new >= l.size() {
                LIVE.fetch_add(new - l.size(), Relaxed) + (new - l.size())
            } else {
                LIVE.fetch_sub(l.size() - new, Relaxed) - (l.size() - new)
            };
            PEAK.fetch_max(now, Relaxed);
        }
        q
    }
}

#[global_allocator]
static A: Counting = Counting;

fn peak_of<T>(f: impl FnOnce() -> T) -> (T, usize) {
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let v = f();
    (v, PEAK.load(Relaxed).saturating_sub(base))
}

/// A SAM of `n` 160-base reads on chr1, streamed to `path`.
fn write_sam(path: &Path, n: usize) {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    writeln!(w, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000").unwrap();
    let seq = "ACGT".repeat(40);
    let qual = "I".repeat(seq.len());
    for i in 0..n {
        writeln!(
            w,
            "read{i}\t0\tchr1\t{}\t60\t160M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0",
            i * 100 + 1
        )
        .unwrap();
    }
    w.flush().unwrap();
}

fn ropts() -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
        chrom_map: None,
        accept: Accept {
            alignments: Some(AlignOptions {
                full_records: true,
                ..AlignOptions::default()
            }),
            ..Accept::default()
        },
    }
}

/// `view --count` of `input`, with `regions_file` if given: the count and the allocator peak.
fn count(input: &Path, regions_file: Option<&Path>, out: &Path) -> (u64, usize) {
    let vo = ViewOptions {
        count: true,
        regions_file: regions_file.map(|p| p.to_string_lossy().into_owned()),
        ..ViewOptions::default()
    };
    let spec = plain_text_output_spec(Some(out));
    let oo = OutputOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        level: 6,
        write_index: false,
    };
    let wo = WriterOptions {
        write_header: true,
        command_line: None,
        version: "0",
    };
    let ((), peak) = peak_of(|| {
        run_view(input.to_str().unwrap(), &vo, &ropts(), &spec, &oo, &wo).unwrap();
    });
    (
        std::fs::read_to_string(out)
            .unwrap()
            .trim()
            .parse()
            .unwrap(),
        peak,
    )
}

/// The same file as BAM, made by `view` (outside the measurement).
fn to_bam(sam: &Path, bam: &Path) {
    let vo = ViewOptions::default();
    let spec = gloci_io::parse_output_spec(Some(bam), None).unwrap();
    let oo = OutputOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        level: 1,
        write_index: false,
    };
    let wo = WriterOptions {
        write_header: true,
        command_line: None,
        version: "0",
    };
    run_view(sam.to_str().unwrap(), &vo, &ropts(), &spec, &oo, &wo).unwrap();
}

#[test]
fn view_count_scan_holds_one_record_and_the_pool() {
    const SMALL: usize = 1_000;
    const BIG: usize = 200_000;
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let bed = d.join("t.bed");
    std::fs::write(&bed, "chr1\t0\t1000\nchr1\t19000000\t19000100\n").unwrap();
    let out = d.join("n.txt");
    let mib = 1usize << 20;
    for (ext, n_small, n_big) in [("sam", SMALL, BIG), ("bam", SMALL, BIG)] {
        let mut files = Vec::new();
        for (tag, n) in [("small", n_small), ("big", n_big)] {
            let sam = d.join(format!("{tag}.sam"));
            if !sam.exists() {
                write_sam(&sam, n);
            }
            let f = if ext == "sam" {
                sam
            } else {
                let bam = d.join(format!("{tag}.bam"));
                to_bam(&sam, &bam);
                bam
            };
            files.push(f);
        }
        // The records' size, not the (very compressible) BAM's.
        let input_bytes = std::fs::metadata(d.join("big.sam")).unwrap().len() as usize;
        let (n, small) = count(&files[0], None, &out);
        assert_eq!(n, SMALL as u64);
        let (n, big) = count(&files[1], None, &out);
        assert_eq!(n, BIG as u64);
        // The big file is far larger than the slack: buffering it would show.
        let slack = 4 * mib;
        assert!(
            input_bytes > 8 * slack,
            "{ext}: input {input_bytes} too small"
        );
        assert!(
            big <= small + slack,
            "{ext}: counting {BIG} reads peaked at {} KiB, {SMALL} reads at {} KiB",
            big / 1024,
            small / 1024
        );
        // A regions file, scanned (no index): records are filtered, not held.
        let (n, rb) = count(&files[1], Some(&bed), &out);
        assert!(n > 0 && n < BIG as u64);
        assert!(
            rb <= small + slack,
            "{ext}: regions-file count peaked at {} KiB",
            rb / 1024
        );
    }
}
