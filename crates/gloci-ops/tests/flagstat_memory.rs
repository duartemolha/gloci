//! A `flagstat` scan holds one record and the readers' buffers, however many records the file
//! has (Global Constraints, scan memory bound). A counting allocator measures the peak of a
//! small and of a file 100 times larger, for BAM and SAM.
use gloci_core::VcfSpan;
use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions, open_aln_source};
use gloci_ops::flagstat;
use std::alloc::{GlobalAlloc, Layout, System};
use std::num::NonZeroUsize;
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

fn sam(n: usize) -> String {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:10000000\n");
    for i in 0..n {
        let flag = [0, 16, 99, 147, 4, 256, 2048, 1024][i % 8];
        s.push_str(&format!(
            "read{i}\t{flag}\tc1\t{}\t60\t10M\t=\t{}\t0\tACGTACGTAC\tIIIIIIIIII\n",
            i % 1_000_000 + 1,
            i % 1_000_000 + 1
        ));
    }
    s
}

/// Peak bytes above the baseline while counting `path`.
fn measure(path: &std::path::Path, threads: usize) -> (u64, usize) {
    let o = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions::default()),
            convert_only: false,
        },
        chrom_map: None,
        threads: NonZeroUsize::new(threads).unwrap(),
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let stats = flagstat(open_aln_source(path.to_str().unwrap(), &o).unwrap()).unwrap();
    (stats.total[0], PEAK.load(Relaxed) - base)
}

#[test]
fn flagstat_scan_memory_is_bounded() {
    let d = tempfile::tempdir().unwrap();
    let (small_n, big_n) = (2_000, 200_000);
    for (ext, bam) in [("sam", false), ("bam", true)] {
        let mut paths = Vec::new();
        for n in [small_n, big_n] {
            let text = sam(n);
            let p = d.path().join(format!("r{n}.{ext}"));
            let bytes = if bam {
                gloci_testkit::fixtures::sam_to_bam(&text)
            } else {
                text.into_bytes()
            };
            std::fs::write(&p, bytes).unwrap();
            paths.push(p);
        }
        for threads in [1, 2] {
            let (n1, small) = measure(&paths[0], threads);
            let (n2, big) = measure(&paths[1], threads);
            assert_eq!((n1, n2), (small_n as u64, big_n as u64));
            // 100 times the records, no more memory (a few KiB of slack for allocator noise
            // and the thread pool's queues).
            assert!(
                big <= small + 256 * 1024,
                "{ext} -t {threads}: {small} B for {small_n} records, {big} B for {big_n}"
            );
            assert!(big < 8 * 1024 * 1024, "{ext} -t {threads}: {big} B");
        }
    }
}
