//! `faidx` holds at most one fetch chunk plus a small block, however long the
//! region is. A counting allocator measures the peak while a 300 Mb sequence
//! (a sparse file) is written, forwards and reverse-complemented.
use gloci_core::Interval;
use gloci_io::{FETCH_CHUNK, IndexedFasta};
use gloci_ops::faidx::{FaidxOptions, FaidxRegion, FaidxWriter, RegionsFile, faidx};
use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Write;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

/// The tests measure one global counter, so they run one at a time.
static ONE_AT_A_TIME: Mutex<()> = Mutex::new(());

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

/// Counts what is written and keeps nothing.
struct Sink {
    bytes: u64,
    lines: u64,
}

impl Write for Sink {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes += b.len() as u64;
        self.lines += b.iter().filter(|&&c| c == b'\n').count() as u64;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn faidx_memory() {
    let _g = ONE_AT_A_TIME.lock().unwrap_or_else(|e| e.into_inner());
    const LEN: u64 = 300_000_000;
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("big.fa");
    // ">big\n" + LEN bases on one line + "\n"; the bases are NUL bytes in a
    // sparse file, which faidx passes through like any other byte.
    let f = std::fs::File::create(&p).unwrap();
    (&f).write_all(b">big\n").unwrap();
    f.set_len(5 + LEN + 1).unwrap();
    std::fs::write(
        d.path().join("big.fa.fai"),
        format!("big\t{LEN}\t5\t{LEN}\t{}\n", LEN + 1),
    )
    .unwrap();
    let mut fasta = IndexedFasta::open(&p).unwrap();
    let r = [FaidxRegion {
        label: "big".into(),
        chrom: b"big".to_vec(),
        iv: Interval {
            start: 0,
            end: u64::MAX,
        },
    }];
    let one_line = 1024usize;
    for reverse_complement in [false, true] {
        let o = FaidxOptions {
            line_width: Some(one_line as u64),
            reverse_complement,
            ..Default::default()
        };
        let mut sink = Sink { bytes: 0, lines: 0 };
        let base = LIVE.load(Relaxed);
        PEAK.store(base, Relaxed);
        faidx(&mut fasta, &r, &o, &mut sink).unwrap();
        let held = PEAK.load(Relaxed) - base;
        let header = if reverse_complement { 8 } else { 5 };
        assert_eq!(sink.bytes, header + LEN + LEN.div_ceil(one_line as u64));
        assert!(
            held <= FETCH_CHUNK as usize + one_line + 64 * 1024,
            "rc={reverse_complement}: held {held} bytes for a {LEN}-base sequence"
        );
    }
}

#[test]
fn faidx_regions_file_streams() {
    let _g = ONE_AT_A_TIME.lock().unwrap_or_else(|e| e.into_inner());
    const LINES: u64 = 1_000_000;
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("s.fa");
    std::fs::write(&p, ">s\nACGTACGTAC\nGTACGTACGT\n").unwrap();
    std::fs::write(d.path().join("s.fa.fai"), "s\t20\t3\t10\t11\n").unwrap();
    let bed = d.path().join("r.bed");
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&bed).unwrap());
        for i in 0..LINES {
            // a zero-length line every 1000, a missing sequence every 1001
            match i {
                _ if i % 1000 == 0 => writeln!(w, "s\t5\t5").unwrap(),
                _ if i % 1001 == 0 => writeln!(w, "nope\t1\t2").unwrap(),
                _ => writeln!(w, "s\t{}\t{}", i % 10, 11 + i % 9).unwrap(),
            }
        }
    }
    let ro = gloci_io::ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: std::num::NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::default(),
        invalid: gloci_io::InvalidPolicy::Fail,
    };
    let mut fasta = IndexedFasta::open(&p).unwrap();
    let o = FaidxOptions::default();
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let rf = RegionsFile::scan(bed.to_str().unwrap(), &ro, &fasta).unwrap();
    let mut sink = Sink { bytes: 0, lines: 0 };
    let mut w = FaidxWriter::new(&o);
    rf.write(&mut fasta, &ro, true, &mut w, &mut sink).unwrap();
    let held = PEAK.load(Relaxed) - base;
    assert_eq!(rf.zero_length, 1000);
    assert_eq!(rf.missing, 999);
    assert!(
        sink.lines > LINES,
        "regions were written: {} lines",
        sink.lines
    );
    // 1M regions held as (name, interval) and a label each would be >100 MB
    assert!(
        held < 8 * 1024 * 1024,
        "held {held} bytes for {LINES} BED lines"
    );
}
