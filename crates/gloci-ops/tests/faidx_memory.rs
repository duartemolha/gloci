//! `faidx` holds at most one fetch chunk plus a small block, however long the
//! region is. A counting allocator measures the peak while a 300 Mb sequence
//! (a sparse file) is written, forwards and reverse-complemented.
use gloci_core::Interval;
use gloci_io::{FETCH_CHUNK, IndexedFasta};
use gloci_ops::faidx::{FaidxOptions, FaidxRegion, faidx};
use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Write;
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
