//! `dict` holds about one 1 MiB buffer, however long the sequence is. A
//! counting allocator measures the peak while a 1 GB single-line sequence (a
//! sparse file of NUL bytes, which count for nothing) and a 3 MiB line that
//! straddles the stream's pieces are read.
use gloci_io::FastaStream;
use gloci_ops::dict::{DictOptions, dict};
use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Write;
use std::num::NonZeroUsize;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

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

fn measure(path: &std::path::Path) -> (String, usize) {
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let s = FastaStream::open(path.to_str().unwrap(), NonZeroUsize::new(1).unwrap()).unwrap();
    let o = DictOptions {
        uri: Some("x".into()),
        ..Default::default()
    };
    let mut out = Vec::new();
    dict(s, None, &o, &mut out).unwrap();
    let held = PEAK.load(Relaxed) - base;
    (String::from_utf8(out).unwrap(), held)
}

#[test]
fn dict_streams() {
    let _g = ONE_AT_A_TIME.lock().unwrap_or_else(|e| e.into_inner());
    const LEN: u64 = 1_000_000_000;
    let d = tempfile::tempdir().unwrap();
    let big = d.path().join("big.fa");
    let f = std::fs::File::create(&big).unwrap();
    (&f).write_all(b">big\n").unwrap();
    f.set_len(5 + LEN + 1).unwrap();
    let (got, held) = measure(&big);
    assert_eq!(
        got,
        "@SQ\tSN:big\tLN:0\tM5:d41d8cd98f00b204e9800998ecf8427e\tUR:x\n"
    );
    assert!(held < 3 << 20, "held {held} bytes for a 1 GB line");

    // Real bases over several pieces: 3 MiB + 17 lower-case A on one line.
    let n = 3 * (1 << 20) + 17;
    let mut text = b">a\n".to_vec();
    text.extend(std::iter::repeat_n(b'a', n));
    text.push(b'\n');
    let p = d.path().join("a.fa");
    std::fs::write(&p, &text).unwrap();
    let (got, held) = measure(&p);
    use md5::{Digest, Md5};
    let want = format!("{:x}", Md5::digest(vec![b'A'; n]));
    assert_eq!(got, format!("@SQ\tSN:a\tLN:{n}\tM5:{want}\tUR:x\n"));
    assert!(held < 3 << 20, "held {held} bytes");
}
