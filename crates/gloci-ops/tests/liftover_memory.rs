//! `gloci liftover` holds the chain index and one batch of records, whatever the size of the
//! input (D174, R7 and R12). A counting allocator measures the peak of lifting 20,000, 100,000 and
//! 400,000 records of the Ensembl chromosome 21 and 22 slice. The smallest run already fills one
//! batch (16,384 records, about 1.6 MB), which is the bound; a run of 1,000 records peaks lower.
use gloci_core::{ChromOrder, VcfSpan};
use gloci_io::{
    Accept, ChainIndex, InvalidPolicy, IoError, ReaderOptions, RecordWriter, chrom_table_for,
    open_reader,
};
use gloci_ops::liftover::{LiftoverOptions, liftover};
use std::alloc::{GlobalAlloc, Layout, System};
use std::io;
use std::num::NonZeroUsize;
use std::path::PathBuf;
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

/// Drops every record, keeping a count.
struct Drop(u64);

impl RecordWriter for Drop {
    fn write_record(&mut self, _raw: &[u8]) -> io::Result<()> {
        self.0 += 1;
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/liftover")
        .join(name)
}

/// Peak bytes above the baseline while lifting `path`, and the records written.
fn measure(chain: &ChainIndex, path: &std::path::Path) -> (usize, u64) {
    let ro = ReaderOptions {
        accept: Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let reader = open_reader(path.to_str().unwrap(), &ro).unwrap();
    let mut chroms = chrom_table_for(None, &[reader.as_ref()], ChromOrder::Karyotype).unwrap();
    let mut out = Drop(0);
    liftover(
        reader,
        &mut chroms,
        chain,
        &LiftoverOptions::default(),
        &mut out,
        None,
    )
    .unwrap();
    (PEAK.load(Relaxed) - base, out.0)
}

#[test]
fn liftover_memory_is_flat() {
    let chain_text = std::fs::read(fixture("ensembl-grch37-to-grch38.21-22.chain")).unwrap();
    let chain = ChainIndex::parse("ensembl.chain", &chain_text).unwrap();
    drop(chain_text);
    let bed = std::fs::read_to_string(fixture("grch37-21-22.bed")).unwrap();
    let lines: Vec<&str> = bed.lines().collect();
    let d = tempfile::tempdir().unwrap();
    let mut peaks = Vec::new();
    for n in [20_000usize, 100_000, 400_000] {
        let mut text = String::new();
        for i in 0..n {
            text.push_str(lines[i % lines.len()]);
            text.push('\n');
        }
        let p = d.path().join(format!("r{n}.bed"));
        std::fs::write(&p, text).unwrap();
        let (peak, written) = measure(&chain, &p);
        assert!(written >= n as u64 / 2, "{n} records wrote {written}");
        peaks.push(peak);
    }
    assert!(
        peaks[2] < peaks[0] + 1024 * 1024 && peaks[1] < peaks[0] + 1024 * 1024,
        "{} B for 20,000 records, {} B for 100,000, {} B for 400,000",
        peaks[0],
        peaks[1],
        peaks[2]
    );
}
