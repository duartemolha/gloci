//! A region as long as a chromosome costs what its longest read costs, not what its length does
//! (D140): a 250 Mb region over 1,000,000 reads keeps the ring at 2^12 counters and the heap
//! under a megabyte, measured by an in-process counting allocator.

use gloci_core::Interval;
use gloci_io::{AlnHeader, AlnSource, Diagnostics, Format, IoError, SortOrder};
use gloci_ops::bedcov::region_totals_in;
use gloci_ops::depth::{CountRules, DepthRing};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::Arc;
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

const LN: u64 = 250_000_000;
const READS: u64 = 1_000_000;
const STEP: u64 = 250;
const LEN: u32 = 100;

/// A million reads of 100 M, one every 250 bases, made as they are asked for.
struct Generated {
    header: Arc<AlnHeader>,
    next: u64,
    diag: Diagnostics,
}

impl AlnSource for Generated {
    fn name(&self) -> &str {
        "generated.bam"
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn header(&self) -> &Arc<AlnHeader> {
        &self.header
    }
    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        if self.next == READS {
            return Ok(None);
        }
        let pos = (self.next * STEP) as i32;
        buf.extend_from_slice(&0i32.to_le_bytes());
        buf.extend_from_slice(&pos.to_le_bytes());
        buf.push(2); // l_read_name
        buf.push(60);
        buf.extend_from_slice(&0u16.to_le_bytes()); // bin
        buf.extend_from_slice(&1u16.to_le_bytes()); // n_cigar_op
        buf.extend_from_slice(&0u16.to_le_bytes()); // flag
        buf.extend_from_slice(&LEN.to_le_bytes());
        buf.extend_from_slice(&(-1i32).to_le_bytes());
        buf.extend_from_slice(&(-1i32).to_le_bytes());
        buf.extend_from_slice(&0i32.to_le_bytes());
        buf.extend_from_slice(b"r\0");
        buf.extend_from_slice(&(LEN << 4).to_le_bytes());
        buf.extend(std::iter::repeat_n(0x11u8, LEN.div_ceil(2) as usize));
        buf.extend(std::iter::repeat_n(30u8, LEN as usize));
        self.next += 1;
        Ok(Some(self.next))
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }
}

#[test]
fn bedcov_whole_chromosome_region() {
    let header = Arc::new(AlnHeader::new(
        format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:big\tLN:{LN}\n").into_bytes(),
        vec![(b"big".to_vec(), LN)],
        SortOrder::Coordinate,
    ));
    let mut src = Generated {
        header,
        next: 0,
        diag: Diagnostics::default(),
    };
    let mut ring = DepthRing::new();
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let t = region_totals_in(
        &mut ring,
        &mut src,
        LN,
        Interval { start: 0, end: LN },
        &CountRules {
            deletions: true,
            ref_skips: true,
            min_baseq: 0,
        },
        Some(1),
        &mut None,
    )
    .unwrap();
    let held = PEAK.load(Relaxed) - base;
    assert_eq!(t.reads, READS);
    assert_eq!(t.sum, READS * u64::from(LEN));
    assert_eq!(t.at_depth, READS * u64::from(LEN));
    assert!(ring.capacity() <= 1 << 12, "{}", ring.capacity());
    assert!(held < 1 << 20, "the heap peaked at {held} bytes");
}
