//! BGZF writing on workers that sleep while they wait (design §8.3, stage 1).
//!
//! [`Writer`] cuts what is written to it into blocks on the calling thread, has each block
//! compressed on a worker of a [`Pool`], and writes the blocks to its sink in order, on the
//! calling thread again.
//!
//! It does what `noodles_bgzf::io::MultithreadedWriter` does, and leans on noodles for it:
//! the compressed stream is the same byte for byte at every worker count, with the same
//! block size and each block compressed by noodles' single-threaded writer. What differs is
//! how the threads wait: noodles' writer runs a rayon pool and crossbeam channels, which
//! wait by yielding the CPU ([`crate::handoff`] has the cost).

use crate::error::IoError;
use crate::handoff::Pool;
use crate::input::BGZF_EOF;
use noodles_bgzf::io::writer::{Builder as BlockBuilder, CompressionLevel};
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex, MutexGuard};

/// The error of a writer whose worker threads cannot be started.
pub(crate) fn start_error(e: io::Error) -> IoError {
    IoError::Internal(format!("cannot start a worker thread: {e}"))
}

/// The most data one block holds: noodles' write buffer (`MAX_BUF_SIZE`, private there),
/// which is the 64 KiB a block may inflate to less the BGZF header (18), the gzip trailer
/// (8) and the most a stored deflate stream adds (15). A test compares it with noodles.
const BLOCK_DATA: usize = (1 << 16) - 18 - 8 - 15;

/// The most a block takes when it is stored: its size less one is kept in 16 bits.
const MAX_STORED: usize = 1 << 16;

/// A block to compress, and a buffer to compress it into.
type Deflate = (Vec<u8>, Vec<u8>);
/// The compressed block as it is stored, and the buffer the data came in.
type Deflated = (Vec<u8>, Vec<u8>);

/// Where a worker's block writer puts a block: a buffer the worker swaps in before a block
/// and takes out after it, while the writer keeps its end.
#[derive(Clone, Default)]
struct Stored(Arc<Mutex<Vec<u8>>>);

impl Stored {
    fn buffer(&self) -> MutexGuard<'_, Vec<u8>> {
        self.0.lock().unwrap_or_else(|e| e.into_inner())
    }
}

impl Write for Stored {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.buffer().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// What a compression worker keeps: noodles' single-threaded writer, which makes one block
/// of what it is given up to a flush, with the same deflate call and the same framing as its
/// multithreaded one. The worker keeps one writer for every block because a writer allocates
/// two 64 KiB buffers: made and freed for each block, they left about 2 MiB more resident in
/// every worker (measured with 8 workers).
struct Deflater {
    level: CompressionLevel,
    block: noodles_bgzf::io::Writer<Stored>,
    stored: Stored,
}

impl Deflater {
    fn new(level: CompressionLevel) -> Deflater {
        let stored = Stored::default();
        Deflater {
            level,
            block: BlockBuilder::default()
                .set_compression_level(level)
                .build_from_writer(stored.clone()),
            stored,
        }
    }
}

/// One block of data as one stored BGZF block.
fn deflate(deflater: &mut Deflater, (data, mut stored): Deflate) -> io::Result<Deflated> {
    stored.clear();
    *deflater.stored.buffer() = stored;
    let written = deflater
        .block
        .write_all(&data)
        .and_then(|()| deflater.block.flush());
    let stored = std::mem::take(&mut *deflater.stored.buffer());
    if let Err(e) = written {
        // The writer may still hold the block: start again with a new one.
        *deflater = Deflater::new(deflater.level);
        return Err(e);
    }
    Ok((stored, data))
}

fn ended() -> io::Error {
    io::Error::other("the BGZF output has ended")
}

/// A BGZF writer whose blocks are compressed on worker threads.
///
/// - [`Write::write`] takes at most what fills the current block, as noodles' writers do.
/// - [`Write::flush`] ends the current block, so a block may be short; it does not wait for
///   the block to reach the sink, nor flush the sink.
/// - [`Writer::finish`] writes what is left and the end-of-file block, and returns the sink.
///   A writer dropped without it finishes too, ignoring errors.
/// - An error from the sink or from a worker comes back from the `write`, `flush` or `finish`
///   that meets it, and nothing more reaches the sink after it.
/// - At most twice the worker count of blocks are compressed or wait to be written.
pub(crate) struct Writer<W: Write> {
    /// `None` once the stream is finished or has failed.
    sink: Option<W>,
    pool: Pool<Deflate, Deflated>,
    /// The blocks that may be in flight.
    cap: u64,
    /// The block being filled.
    block: Vec<u8>,
    /// Buffers of blocks that have been written, to fill and compress into again.
    spare: Vec<(Vec<u8>, Vec<u8>)>,
}

impl<W: Write> Writer<W> {
    /// A writer to `sink` with `workers` compression threads.
    pub fn new(sink: W, level: CompressionLevel, workers: NonZeroUsize) -> io::Result<Self> {
        Ok(Writer {
            sink: Some(sink),
            pool: Pool::new(
                "bgzf-deflate",
                (0..workers.get()).map(|_| Deflater::new(level)).collect(),
                deflate,
            )?,
            cap: 2 * workers.get() as u64,
            block: Vec::with_capacity(BLOCK_DATA),
            spare: Vec::new(),
        })
    }

    /// Write the blocks that are compressed to the sink, in order, waiting for the earliest
    /// while more than `allowed` are in flight.
    fn settle(&mut self, allowed: u64) -> io::Result<()> {
        loop {
            let next = if self.pool.in_flight() > allowed {
                self.pool.take()
            } else {
                self.pool.try_take()
            };
            let Some(deflated) = next else {
                return Ok(());
            };
            let sink = self.sink.as_mut().ok_or_else(ended)?;
            let written = deflated.and_then(|(stored, data)| {
                sink.write_all(&stored)?;
                Ok((data, stored))
            });
            match written {
                Ok(buffers) => self.spare.push(buffers),
                Err(e) => {
                    // Nothing follows a block that failed: not the later blocks, not the
                    // end-of-file block.
                    self.sink = None;
                    return Err(e);
                }
            }
        }
    }

    /// Write what is left, then the end-of-file block, and return the sink.
    pub fn finish(&mut self) -> io::Result<W> {
        self.flush()?;
        self.settle(0)?;
        let mut sink = self.sink.take().ok_or_else(ended)?;
        sink.write_all(&BGZF_EOF)?;
        Ok(sink)
    }
}

impl<W: Write> Write for Writer<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = (BLOCK_DATA - self.block.len()).min(buf.len());
        self.block.extend_from_slice(&buf[..n]);
        if self.block.len() == BLOCK_DATA {
            self.flush()?;
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        if self.sink.is_none() {
            return Err(ended());
        }
        // Make room first, so the blocks in flight never pass the bound.
        self.settle(self.cap - 1)?;
        // Buffers are made here, at their full size, and come back: the workers allocate
        // nothing for a block.
        let (mut next, stored) = self.spare.pop().unwrap_or_else(|| {
            (
                Vec::with_capacity(BLOCK_DATA),
                Vec::with_capacity(MAX_STORED),
            )
        });
        next.clear();
        let data = std::mem::replace(&mut self.block, next);
        self.pool.submit((data, stored));
        Ok(())
    }
}

impl<W: Write> Drop for Writer<W> {
    fn drop(&mut self) {
        if self.sink.is_some() {
            let _ = self.finish();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_bgzf::io::MultithreadedWriter;
    use noodles_bgzf::io::multithreaded_writer::Builder as NoodlesBuilder;

    fn level(n: u8) -> CompressionLevel {
        CompressionLevel::new(n).unwrap()
    }

    fn workers(n: usize) -> NonZeroUsize {
        NonZeroUsize::new(n).unwrap()
    }

    /// `n` bytes that compress, but not to nothing.
    fn text(n: usize, seed: u32) -> Vec<u8> {
        let mut x = seed.wrapping_mul(2_654_435_761).max(1);
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 17;
                x ^= x << 5;
                b"ACGTacgt\t\n0123456789chr."[(x >> 9) as usize % 24]
            })
            .collect()
    }

    /// Write `data` in pieces of `piece` bytes, flushing after the pieces listed.
    fn write_pieces(w: &mut dyn Write, data: &[u8], piece: usize, flush_after: &[usize]) {
        for (i, c) in data.chunks(piece.max(1)).enumerate() {
            w.write_all(c).unwrap();
            if flush_after.contains(&i) {
                w.flush().unwrap();
            }
        }
    }

    fn ours(data: &[u8], lvl: u8, n: usize, piece: usize, flush_after: &[usize]) -> Vec<u8> {
        let mut w = Writer::new(Vec::new(), level(lvl), workers(n)).unwrap();
        write_pieces(&mut w, data, piece, flush_after);
        w.finish().unwrap()
    }

    fn noodles(data: &[u8], lvl: u8, n: usize, piece: usize, flush_after: &[usize]) -> Vec<u8> {
        let mut w: MultithreadedWriter<Vec<u8>> = NoodlesBuilder::default()
            .set_compression_level(level(lvl))
            .set_worker_count(workers(n))
            .build_from_writer(Vec::new());
        write_pieces(&mut w, data, piece, flush_after);
        w.finish().unwrap()
    }

    /// The size and the inflated size of each block of `bgzf`.
    fn blocks(bgzf: &[u8]) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut at = 0;
        while at < bgzf.len() {
            let size = u16::from_le_bytes([bgzf[at + 16], bgzf[at + 17]]) as usize + 1;
            let isize = u32::from_le_bytes(bgzf[at + size - 4..at + size].try_into().unwrap());
            out.push((size, isize as usize));
            at += size;
        }
        out
    }

    #[test]
    fn a_block_holds_what_a_noodles_block_holds() {
        let data = text(BLOCK_DATA + 1, 1);
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(&data).unwrap();
        let bgzf = w.finish().unwrap();
        let sizes: Vec<usize> = blocks(&bgzf).iter().map(|b| b.1).collect();
        assert_eq!(sizes, [BLOCK_DATA, 1, 0]);
    }

    /// The compressed stream is the one noodles' multithreaded writer writes, whatever the
    /// worker count, the level, the sizes of the writes and the flushes between them.
    #[test]
    fn writes_the_bytes_the_noodles_multithreaded_writer_writes() {
        let big = text(5 * BLOCK_DATA + 12_345, 7);
        let cases: [(&[u8], usize, &[usize]); 9] = [
            (&[], 1, &[]),
            (b"chr1\t1\t2\n", 3, &[]),
            (&big[..BLOCK_DATA - 1], 4096, &[]),
            (&big[..BLOCK_DATA], 4096, &[]),
            (&big[..BLOCK_DATA + 1], 4096, &[]),
            (&big, 7919, &[]),
            (&big, 1 << 20, &[]),
            // Flushes cut short blocks, also right after a full one and twice in a row.
            (&big, 7919, &[0, 3, 4, 8, 16, 40]),
            (&big[..2 * BLOCK_DATA], BLOCK_DATA, &[0, 1]),
        ];
        for (data, piece, flush_after) in cases {
            for lvl in [0, 1, 6, 9, 12] {
                let want = noodles(data, lvl, 1, piece, flush_after);
                assert!(want.ends_with(&BGZF_EOF));
                for n in [1, 2, 3, 8] {
                    let got = ours(data, lvl, n, piece, flush_after);
                    assert!(
                        got == want,
                        "{} bytes in pieces of {piece}, level {lvl}, {n} workers: the blocks \
are {:?}, noodles writes {:?}",
                        data.len(),
                        blocks(&got),
                        blocks(&want)
                    );
                }
            }
        }
        // noodles writes the same at any worker count of its own.
        assert_eq!(
            noodles(&big, 6, 4, 7919, &[3]),
            noodles(&big, 6, 1, 7919, &[3])
        );
    }

    /// A sink the test can look into while the writer owns it. It refuses the write that
    /// would take it past `limit` bytes, once: whatever is written after that is taken.
    #[derive(Clone)]
    struct Shared {
        bytes: Arc<Mutex<Vec<u8>>>,
        limit: Arc<Mutex<usize>>,
    }

    impl Shared {
        fn new(limit: usize) -> Shared {
            Shared {
                bytes: Arc::default(),
                limit: Arc::new(Mutex::new(limit)),
            }
        }
        fn bytes(&self) -> Vec<u8> {
            self.bytes.lock().unwrap().clone()
        }
    }

    impl Write for Shared {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            let mut bytes = self.bytes.lock().unwrap();
            let mut limit = self.limit.lock().unwrap();
            if bytes.len() + buf.len() > *limit {
                *limit = usize::MAX;
                return Err(io::Error::new(
                    io::ErrorKind::StorageFull,
                    "the disk is full",
                ));
            }
            bytes.extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A writer that is dropped writes what it holds and the end-of-file block, as noodles'
    /// writers do; one that is finished adds nothing when it is dropped.
    #[test]
    fn a_dropped_writer_finishes_its_stream() {
        let data = text(3 * BLOCK_DATA + 99, 3);
        let want = noodles(&data, 6, 1, 4096, &[]);
        for n in [1, 4] {
            let sink = Shared::new(usize::MAX);
            let mut w = Writer::new(sink.clone(), level(6), workers(n)).unwrap();
            w.write_all(&data).unwrap();
            drop(w);
            assert_eq!(sink.bytes(), want);

            let sink = Shared::new(usize::MAX);
            let mut w = Writer::new(sink.clone(), level(6), workers(n)).unwrap();
            w.write_all(&data).unwrap();
            drop(w.finish().unwrap());
            drop(w);
            assert_eq!(sink.bytes(), want);
        }
    }

    /// A write the sink refuses comes back from the writer with the sink's error, and
    /// nothing reaches the sink after it, though the sink would take it: no later block and
    /// no end-of-file block, from a write, from `finish` or from the drop.
    #[test]
    fn a_failed_write_surfaces_and_nothing_follows_it() {
        let data = text(40 * BLOCK_DATA, 5);
        let whole = noodles(&data, 6, 1, 4096, &[]);
        for n in [1, 2, 8] {
            let sink = Shared::new(whole.len() / 2);
            let mut w = Writer::new(sink.clone(), level(6), workers(n)).unwrap();
            let failed = w.write_all(&data).and_then(|()| w.finish().map(drop));
            let e = failed.expect_err("the sink fills up half way");
            assert_eq!(e.kind(), io::ErrorKind::StorageFull, "{e}");
            assert_eq!(e.to_string(), "the disk is full");
            let held = sink.bytes();
            assert!(w.finish().is_err());
            assert!(w.write_all(&data).is_err());
            drop(w);
            assert_eq!(
                sink.bytes(),
                held,
                "bytes reached the sink after the failure"
            );
            // What reached the sink is whole blocks from the start of the stream.
            assert!(whole.starts_with(&held));
            assert!(!held.ends_with(&BGZF_EOF));
            let size: usize = blocks(&whole)
                .iter()
                .map(|b| b.0)
                .scan(0, |sum, b| {
                    *sum += b;
                    Some(*sum)
                })
                .take_while(|&end| end <= whole.len() / 2)
                .last()
                .unwrap();
            assert_eq!(held.len(), size);
        }
    }

    /// The end-of-file block that the sink refuses is the error of `finish`.
    #[test]
    fn a_failure_at_the_end_is_the_error_of_finish() {
        let data = text(2 * BLOCK_DATA, 9);
        let whole = noodles(&data, 6, 1, 4096, &[]);
        let sink = Shared::new(whole.len() - 1);
        let mut w = Writer::new(sink.clone(), level(6), workers(2)).unwrap();
        w.write_all(&data).unwrap();
        let e = w.finish().map(drop).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::StorageFull);
        drop(w);
        assert_eq!(sink.bytes(), whole[..whole.len() - BGZF_EOF.len()]);
    }

    /// The writer holds a bounded number of blocks however fast it is written to: what is
    /// not yet in the sink stays within the blocks that may be in flight.
    #[test]
    fn the_blocks_in_flight_are_bounded() {
        let data = text(BLOCK_DATA, 11);
        let stored = noodles(&data, 1, 1, BLOCK_DATA, &[]).len() - BGZF_EOF.len();
        let sink = Shared::new(usize::MAX);
        let mut w = Writer::new(sink.clone(), level(1), workers(2)).unwrap();
        for written in 1..=60 {
            w.write_all(&data).unwrap();
            // Every block but the 2 x 2 that may be in flight is in the sink.
            let in_sink = sink.bytes().len() / stored;
            assert!(
                in_sink + 4 >= written,
                "after {written} blocks the sink holds {in_sink}"
            );
        }
        w.finish().unwrap();
        assert_eq!(sink.bytes().len(), 60 * stored + BGZF_EOF.len());
    }
}
