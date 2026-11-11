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

/// The counters are global: one test at a time, building its files included.
static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

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

/// Peak bytes above the baseline while printing the header and `records` records of `path`.
fn measure(path: &std::path::Path, threads: usize, records: u64) -> usize {
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
        records,
    };
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    head(path.to_str().unwrap(), &h, &o, &mut std::io::sink()).unwrap();
    PEAK.load(Relaxed) - base
}

#[test]
fn head_memory_is_bounded() {
    let _one = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
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
            let small = measure(&paths[0], threads, 2);
            let big = measure(&paths[1], threads, 2);
            assert!(
                big <= small + 256 * 1024,
                "{ext} -t {threads}: {small} B for {small_n} records, {big} B for {big_n}"
            );
            assert!(big < 8 * 1024 * 1024, "{ext} -t {threads}: {big} B");
        }
    }
}

/// The text formats read records in batches of at most `next_batch_size`'s cap (16384 records, a
/// few MiB of bytes), so what `head` holds is the header plus one such batch, a constant: the
/// batch, not one record, is the bound (D134; the memory table of the docs says so). With
/// `--records` above the batch size the peak must still not grow with the input.
#[test]
fn head_large_records_is_bounded_by_the_batch() {
    let _one = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let d = tempfile::tempdir().unwrap();
    let mut peaks = Vec::new();
    for n in [200_000usize, 2_000_000] {
        let p = d.path().join(format!("b{n}.bed"));
        std::fs::write(&p, bed(n)).unwrap();
        peaks.push(measure(&p, 1, 1_000_000));
    }
    assert!(
        peaks[1] <= peaks[0] + 256 * 1024,
        "{} B for 200000 records, {} B for 2000000",
        peaks[0],
        peaks[1]
    );
    assert!(peaks[1] < 16 * 1024 * 1024, "{} B", peaks[1]);
}
