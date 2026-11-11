//! `head` holds the header and a record, whatever the size of the file (Global Constraints,
//! `head` bound). A counting allocator measures the peak of a small file and of one 100 times
//! larger, for BAM, SAM, VCF and BED.
use gloci_core::VcfSpan;
use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions};
use gloci_ops::{HeadOptions, head};
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
        s.push_str(&format!(
            "read{i}\t0\tc1\t{}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
            i % 1_000_000 + 1
        ));
    }
    s
}

fn vcf(n: usize) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=c1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..n {
        s.push_str(&format!("c1\t{}\t.\tA\tT\t.\t.\t.\n", i + 1));
    }
    s
}

fn bed(n: usize) -> String {
    let mut s = String::from("#header\n");
    for i in 0..n {
        s.push_str(&format!("c1\t{}\t{}\n", i * 10, i * 10 + 5));
    }
    s
}

/// Peak bytes above the baseline while printing the header and two records of `path`.
fn measure(path: &std::path::Path, threads: usize) -> usize {
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
    let h = HeadOptions {
        headers: None,
        records: 2,
    };
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    head(path.to_str().unwrap(), &h, &o, &mut std::io::sink()).unwrap();
    PEAK.load(Relaxed) - base
}

#[test]
fn head_memory_is_bounded() {
    let d = tempfile::tempdir().unwrap();
    let (small_n, big_n) = (2_000, 200_000);
    type Maker = fn(usize) -> Vec<u8>;
    let makers: [(&str, Maker); 4] = [
        ("sam", |n| sam(n).into_bytes()),
        ("bam", |n| gloci_testkit::fixtures::sam_to_bam(&sam(n))),
        ("vcf", |n| vcf(n).into_bytes()),
        ("bed", |n| bed(n).into_bytes()),
    ];
    for (ext, make) in makers {
        let mut paths = Vec::new();
        for n in [small_n, big_n] {
            let p = d.path().join(format!("r{n}.{ext}"));
            std::fs::write(&p, make(n)).unwrap();
            paths.push(p);
        }
        for threads in [1, 2] {
            let small = measure(&paths[0], threads);
            let big = measure(&paths[1], threads);
            assert!(
                big <= small + 256 * 1024,
                "{ext} -t {threads}: {small} B for {small_n} records, {big} B for {big_n}"
            );
            assert!(big < 8 * 1024 * 1024, "{ext} -t {threads}: {big} B");
        }
    }
}
