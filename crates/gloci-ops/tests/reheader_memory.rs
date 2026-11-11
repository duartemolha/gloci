//! `reheader` holds the old and new headers and about one BGZF block, however much data follows
//! the header. A counting allocator measures the peak while a BAM with a 200,000-`@SQ` header
//! and several MiB of reads is given a new (renamed) header.
use gloci_io::bgzf_copy::deflate_blocks;
use gloci_io::{
    AlnHeader, BGZF_EOF, OutFormat, OutputOptions, WriterOptions, open_raw_output,
    parse_output_spec, read_bam_file_header, read_new_header,
};
use gloci_ops::reheader::reheader;
use std::alloc::{GlobalAlloc, Layout, System};
use std::fmt::Write as _;
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

/// The peak bytes above the live bytes at the start while `f` runs.
fn peak_of<T>(f: impl FnOnce() -> T) -> (T, usize) {
    let base = LIVE.load(Relaxed);
    PEAK.store(base, Relaxed);
    let v = f();
    (v, PEAK.load(Relaxed) - base)
}

const N: usize = 200_000;
const READS: usize = 30_000;

/// An unmapped read: no reference, 100 bases and qualities from a generator.
fn read(i: usize, state: &mut u64, raw: &mut Vec<u8>) {
    let name = format!("r{i:06}\0");
    let mut seq = [0u8; 50];
    let mut qual = [0u8; 100];
    for b in &mut seq {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *b = (*state >> 33) as u8;
    }
    for q in &mut qual {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *q = ((*state >> 40) % 40) as u8;
    }
    let size = 32 + name.len() + seq.len() + qual.len();
    raw.extend_from_slice(&(size as u32).to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&[name.len() as u8, 0]);
    raw.extend_from_slice(&4680u16.to_le_bytes());
    raw.extend_from_slice(&0u16.to_le_bytes());
    raw.extend_from_slice(&4u16.to_le_bytes());
    raw.extend_from_slice(&100u32.to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&0i32.to_le_bytes());
    raw.extend_from_slice(name.as_bytes());
    raw.extend_from_slice(&seq);
    raw.extend_from_slice(&qual);
}

#[test]
fn reheader_holds_the_headers_and_a_block() {
    let dir = tempfile::tempdir().unwrap();
    let mut old = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    let mut new = old.clone();
    for i in 0..N {
        let len = 1000 + (i % 7000);
        writeln!(old, "@SQ\tSN:contig_{i:06}\tLN:{len}").unwrap();
        writeln!(new, "@SQ\tSN:c{i}\tLN:{len}").unwrap();
    }
    let header = AlnHeader::parse_sam(old.into_bytes()).unwrap();
    let mut raw = Vec::new();
    header.write_bam(&mut raw).unwrap();
    let header_len = raw.len();
    let mut state = 7u64;
    for i in 0..READS {
        read(i, &mut state, &mut raw);
    }
    let mut bam = deflate_blocks(&raw).unwrap();
    bam.extend_from_slice(&BGZF_EOF);
    let data_len = raw.len() - header_len;
    drop((raw, header));
    let input = dir.path().join("in.bam");
    std::fs::write(&input, &bam).unwrap();
    let hdr = dir.path().join("new.sam");
    std::fs::write(&hdr, &new).unwrap();
    let out_path = dir.path().join("out.bam");

    // What reading the two headers costs, as the command reads them, both held at once (the
    // parsers' own temporaries, such as the duplicate-name set, included).
    let (_, headers) = peak_of(|| {
        let (old_h, _) = read_bam_file_header(&input).unwrap();
        let new_h = read_new_header(hdr.to_str().unwrap()).unwrap();
        assert!(old_h.text.len() + new_h.text.len() > 8_000_000);
    });

    let spec = parse_output_spec(Some(&out_path), None).unwrap();
    let wo = WriterOptions {
        write_header: true,
        command_line: Some("gloci reheader".into()),
        version: "0",
    };
    let ((), peak) = peak_of(|| {
        let oo = OutputOptions {
            threads: NonZeroUsize::MIN,
            level: 6,
            write_index: false,
        };
        let out = open_raw_output(&spec, OutFormat::Bam, &oo).unwrap();
        reheader(hdr.to_str().unwrap(), &input, out, &wo, 6).unwrap();
    });
    // The two headers, the new header compressed (about 1.3 MB), the output's 1 MiB buffer and
    // a few 64 KiB blocks. Never the data.
    let slack = (2 << 20) + (1 << 20) + 4 * 65536;
    assert!(data_len > slack, "{data_len} <= {slack}");
    assert!(
        peak <= headers + slack,
        "peak {peak} bytes, headers {headers}, data {data_len}"
    );
    // And the result is right.
    let (got, _) = read_bam_file_header(&out_path).unwrap();
    assert_eq!(got.refs.len(), N);
    assert_eq!(got.refs[17].0, b"c17");
    assert!(std::fs::metadata(&out_path).unwrap().len() as usize > data_len / 2);
}
