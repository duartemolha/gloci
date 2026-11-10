//! The BGZF block copy behind `vcf concat --naive` and `gloci cat` (D74, D123).

use gloci_io::BGZF_EOF;
use gloci_io::bgzf_copy::{
    HeadReader, TrackedWrite, copy_from, deflate_blocks, inflate, malformed_reason, read_block,
};
use std::io::{self, Cursor, Read, Write};

/// A stream with no `Seek`: the copy must work on stdin.
struct Pipe(Cursor<Vec<u8>>);

impl Read for Pipe {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.0.read(buf)
    }
}

/// Each part as its own block(s), then the EOF marker.
fn bgzf(parts: &[&[u8]]) -> Vec<u8> {
    let mut v = Vec::new();
    for p in parts {
        v.extend(deflate_blocks(p).unwrap());
    }
    v.extend_from_slice(&BGZF_EOF);
    v
}

fn inflate_all(mut data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    while let Some(b) = read_block(&mut data).unwrap() {
        out.extend(inflate(&b).unwrap());
    }
    out
}

fn noise(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) as u8
        })
        .collect()
}

#[test]
fn copy_from_mid_block_and_edge() {
    let (a, b, c) = (noise(1000, 1), noise(2000, 2), noise(3000, 3));
    let file = bgzf(&[&a, &b, &c]);
    let blocks: Vec<usize> = {
        let mut r = file.as_slice();
        let mut v = Vec::new();
        while let Some(x) = read_block(&mut r).unwrap() {
            v.push(x.raw.len());
        }
        v
    };
    assert_eq!(blocks.len(), 4, "three blocks and the EOF marker");

    // Mid-block: block 1 is recompressed from byte 10, blocks 2 and 3 are copied as they are.
    let mut out = Vec::new();
    let st = copy_from(&mut Pipe(Cursor::new(file.clone())), &mut out, 10).unwrap();
    assert_eq!((st.blocks_recompressed, st.blocks_copied), (1, 2));
    assert!(st.eof_marker);
    assert_eq!(st.uncompressed, 6000);
    let expect = [&a[10..], &b, &c].concat();
    assert_eq!(inflate_all(&out), expect);
    let verbatim = &file[blocks[0]..file.len() - BGZF_EOF.len()];
    assert!(out.ends_with(verbatim), "later blocks are byte for byte");

    // On an edge: nothing is recompressed.
    let mut out = Vec::new();
    let st = copy_from(&mut Pipe(Cursor::new(file.clone())), &mut out, 1000).unwrap();
    assert_eq!((st.blocks_recompressed, st.blocks_copied), (0, 2));
    assert_eq!(out, verbatim);

    // From the start: every block, the EOF marker dropped.
    let mut out = Vec::new();
    copy_from(&mut Pipe(Cursor::new(file.clone())), &mut out, 0).unwrap();
    assert_eq!(out, &file[..file.len() - BGZF_EOF.len()]);

    // Past the end: nothing.
    let mut out = Vec::new();
    let st = copy_from(&mut Pipe(Cursor::new(file)), &mut out, 7000).unwrap();
    assert!(out.is_empty());
    assert_eq!(st.uncompressed, 6000);
}

#[test]
fn a_missing_eof_marker_is_reported_not_added() {
    let mut file = bgzf(&[&noise(500, 4)]);
    file.truncate(file.len() - BGZF_EOF.len());
    let mut out = Vec::new();
    let st = copy_from(&mut file.as_slice(), &mut out, 0).unwrap();
    assert!(!st.eof_marker);
    assert_eq!(out, file);
}

#[test]
fn memory_is_one_block() {
    // About 6.5 MB stored, over 100 blocks of noise: the copy holds one at a time.
    let big = noise(100 * 65280, 5);
    let file = bgzf(&[&big]);
    assert!(file.len() > 6_000_000);
    let mut sizes = Vec::new();
    struct Sizes<'a>(&'a mut Vec<usize>);
    impl Write for Sizes<'_> {
        fn write(&mut self, b: &[u8]) -> io::Result<usize> {
            self.0.push(b.len());
            Ok(b.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let st = copy_from(&mut Pipe(Cursor::new(file)), &mut Sizes(&mut sizes), 5).unwrap();
    // Stored, inflated and compressed again: three buffers of one block at most (plus the
    // few bytes a stored block may exceed 64 KiB by).
    assert!(
        st.peak_block_bytes <= 3 * (65536 + 64),
        "peak {}",
        st.peak_block_bytes
    );
    assert!(sizes.iter().all(|&n| n <= 65536 + 64), "{sizes:?}");
    assert!(sizes.len() >= 99, "{}", sizes.len());
}

#[test]
fn malformed_input_is_an_error_not_a_panic() {
    let file = bgzf(&[&noise(3000, 6), &noise(3000, 7)]);
    let mut out = Vec::new();
    // cut inside a block
    let e = copy_from(&mut &file[..file.len() / 2], &mut out, 0).unwrap_err();
    assert!(malformed_reason(&e).unwrap().contains("ends inside"));
    // plain gzip
    let gz = {
        let mut w = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        w.write_all(b"hello").unwrap();
        w.finish().unwrap()
    };
    let e = copy_from(&mut gz.as_slice(), &mut out, 0).unwrap_err();
    assert!(malformed_reason(&e).unwrap().contains("BGZF"));
    // text
    let e = copy_from(&mut &b"@HD\tVN:1.6\n"[..], &mut out, 0).unwrap_err();
    assert!(malformed_reason(&e).is_some());
    // a corrupt body, and a block that claims to be huge
    let mut bad = file.clone();
    bad[40] ^= 0xff;
    bad[41] ^= 0xff;
    let _ = copy_from(&mut bad.as_slice(), &mut Vec::new(), 1); // error or not, no panic
    let mut huge = file.clone();
    let first = read_block(&mut file.as_slice()).unwrap().unwrap().raw.len();
    huge[first - 1] = 0xff;
    let e = copy_from(&mut huge.as_slice(), &mut Vec::new(), 0).unwrap_err();
    assert!(malformed_reason(&e).is_some());
    // empty input is an empty stream
    let st = copy_from(&mut io::empty(), &mut out, 0).unwrap();
    assert!(!st.eof_marker && st.uncompressed == 0);
}

#[test]
fn write_errors_are_told_from_read_errors() {
    struct Full;
    impl Write for Full {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("disk full"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let file = bgzf(&[&noise(500, 8)]);
    let mut sink = TrackedWrite::new(Full);
    assert!(copy_from(&mut file.as_slice(), &mut sink, 0).is_err());
    assert!(sink.failed);
}

#[test]
fn head_reader_hands_the_whole_stream_back() {
    let (a, b) = (noise(1000, 9), noise(5000, 10));
    let file = bgzf(&[&a, &b]);
    let mut head = HeadReader::new(Cursor::new(file.clone()));
    let mut buf = [0u8; 300];
    head.read_exact(&mut buf).unwrap();
    assert_eq!(buf[..], a[..300]);
    let (stored, rest, produced) = head.finish();
    assert_eq!(produced, 300);
    let mut again = Vec::new();
    Cursor::new(stored)
        .chain(rest)
        .read_to_end(&mut again)
        .unwrap();
    assert_eq!(again, file);
}
