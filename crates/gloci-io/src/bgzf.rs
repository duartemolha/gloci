//! BGZF reading and writing on workers that sleep while they wait (design §8.3, stage 1).
//!
//! [`Writer`] cuts what is written to it into blocks on the calling thread, has each block
//! compressed on a worker of a [`Pool`], and writes the blocks to its sink in order, on the
//! calling thread again. [`Reader`] reads the stored blocks of its source on the calling
//! thread, has each inflated on a worker, and hands their data out in order.
//!
//! Both do what `noodles_bgzf::io::MultithreadedWriter` and `MultithreadedReader` do, and
//! lean on noodles for it. The compressed stream is the same byte for byte at every worker
//! count: the same block size, and each block compressed by noodles' single-block writer.
//! A block is checked and inflated by noodles' single-threaded reader, so a bad block fails
//! with noodles' error. What differs is how the threads wait: noodles' pair run a rayon pool
//! and crossbeam channels each, which wait by yielding the CPU ([`crate::handoff`] has the
//! cost).

use crate::error::IoError;
use crate::handoff::Pool;
use crate::input::BGZF_EOF;
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::writer::{Builder as BlockBuilder, CompressionLevel};
use std::io::{self, BufRead, BufReader, Cursor, Read, Write};
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex, MutexGuard};

/// The error of a reader or writer whose worker threads cannot be started.
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

/// The most a block inflates to.
const MAX_ISIZE: usize = 1 << 16;
/// A block's BGZF header, which ends with the block's size, and the gzip trailer after its
/// data.
const HEADER: usize = 18;
const TRAILER: usize = 8;
/// How much of the source is read at a time.
const SOURCE_BUF: usize = 256 << 10;

/// A stored block to inflate, and a buffer to inflate it into.
type Inflate = (Vec<u8>, Vec<u8>);
/// The data of a block, and the block as it was stored.
type Inflated = (Vec<u8>, Vec<u8>);
/// noodles' reader over one stored block at a time.
type BlockReader = noodles_bgzf::io::Reader<Cursor<Vec<u8>>>;

/// The next block of `src` as it is stored, in `stored`; false at the end of the stream.
///
/// This is noodles' frame reader (private there), error for error: a stream that ends
/// inside the 18-byte header ends there without an error, a block too small to hold a
/// header and a trailer is an invalid frame size, and one that ends early fails to fill its
/// buffer. The header's other fields are checked where the block is inflated.
fn read_stored(src: &mut impl Read, stored: &mut Vec<u8>) -> io::Result<bool> {
    stored.resize(HEADER, 0);
    match src.read_exact(stored) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e),
    }
    let size = usize::from(u16::from_le_bytes([stored[16], stored[17]])) + 1;
    if size < HEADER + TRAILER {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid frame size",
        ));
    }
    stored.resize(size, 0);
    src.read_exact(&mut stored[HEADER..])?;
    Ok(true)
}

/// The data of one stored block: its header, sizes and checksum checked and its data
/// inflated by noodles' single-threaded reader, as its multithreaded one does.
fn inflate(reader: &mut BlockReader, (stored, mut data): Inflate) -> io::Result<Inflated> {
    *reader.get_mut() = Cursor::new(stored);
    // noodles inflates straight into a buffer that can hold any block.
    data.resize(MAX_ISIZE, 0);
    let n = reader.read(&mut data)?;
    data.truncate(n);
    Ok((data, std::mem::take(reader.get_mut()).into_inner()))
}

/// A BGZF reader whose blocks are inflated on worker threads.
///
/// - Blocks without data are skipped, and the stream ends where its source ends.
/// - An error is reported in its place in the stream, after the data of the blocks before
///   it. After a block that fails the reader goes on with the next one; after a source that
///   fails, or a block whose size cannot be read, the stream has ended.
/// - [`Reader::virtual_position`] is where the next byte comes from: the block's offset in
///   the file and the offset within its data, or the start of the next block once a block
///   is read to its end.
/// - At most twice the worker count of blocks are read ahead of the one being handed out.
pub(crate) struct Reader<R> {
    src: BufReader<R>,
    pool: Pool<Inflate, Inflated>,
    /// The blocks that may be in flight.
    cap: u64,
    /// No more blocks are read from the source: it has ended, or failed.
    ended: bool,
    /// What ended the source early, reported once the blocks before it are read.
    failed: Option<io::Error>,
    /// The current block's data, and how much of it is read.
    data: Vec<u8>,
    at: usize,
    /// Where the current block starts in the file, and its stored size: the block after it
    /// starts where it ends.
    block_at: u64,
    block_size: u64,
    /// Buffers of blocks that have been read, to read and inflate into again.
    spare: Vec<(Vec<u8>, Vec<u8>)>,
}

impl<R: Read> Reader<R> {
    /// A reader of `src` with `workers` inflating threads.
    pub fn new(src: R, workers: NonZeroUsize) -> io::Result<Self> {
        let readers = (0..workers.get())
            .map(|_| BlockReader::new(Cursor::new(Vec::new())))
            .collect();
        Ok(Reader {
            src: BufReader::with_capacity(SOURCE_BUF, src),
            pool: Pool::new("bgzf-inflate", readers, inflate)?,
            cap: 2 * workers.get() as u64,
            ended: false,
            failed: None,
            data: Vec::with_capacity(MAX_ISIZE),
            at: 0,
            block_at: 0,
            block_size: 0,
            spare: Vec::new(),
        })
    }

    /// Where the next byte read comes from, as a BGZF virtual position.
    pub fn virtual_position(&self) -> VirtualPosition {
        let (block, within) = if self.at < self.data.len() {
            (self.block_at, self.at as u16)
        } else {
            (self.block_at + self.block_size, 0)
        };
        VirtualPosition::try_from((block, within)).expect("a BGZF file is smaller than 256 TiB")
    }

    /// Read stored blocks from the source and queue them, up to the bound.
    fn read_ahead(&mut self) {
        while !self.ended && self.pool.in_flight() < self.cap {
            // Buffers are made here, at their full size, and come back: the workers
            // allocate nothing for a block.
            let (mut stored, data) = self.spare.pop().unwrap_or_else(|| {
                (
                    Vec::with_capacity(MAX_STORED),
                    Vec::with_capacity(MAX_ISIZE),
                )
            });
            match read_stored(&mut self.src, &mut stored) {
                Ok(true) => self.pool.submit((stored, data)),
                Ok(false) => self.ended = true,
                Err(e) => {
                    self.ended = true;
                    self.failed = Some(e);
                }
            }
        }
    }

    /// Make the next block that holds data the current one. At the end of the stream there
    /// is no data left, and the position stays after the last block.
    fn next_block(&mut self) -> io::Result<()> {
        loop {
            self.read_ahead();
            let Some(inflated) = self.pool.take() else {
                if let Some(e) = self.failed.take() {
                    return Err(e);
                }
                self.data.clear();
                self.at = 0;
                return Ok(());
            };
            let (data, stored) = inflated?;
            self.block_at += self.block_size;
            self.block_size = stored.len() as u64;
            let read = std::mem::replace(&mut self.data, data);
            self.at = 0;
            self.spare.push((stored, read));
            if !self.data.is_empty() {
                return Ok(());
            }
        }
    }
}

impl<R: Read> Read for Reader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let data = self.fill_buf()?;
        let n = data.len().min(buf.len());
        buf[..n].copy_from_slice(&data[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl<R: Read> BufRead for Reader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.at >= self.data.len() {
            self.next_block()?;
        }
        Ok(&self.data[self.at..])
    }

    fn consume(&mut self, amt: usize) {
        self.at = (self.at + amt).min(self.data.len());
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

    /// What a reader under test and noodles' multithreaded reader have in common.
    trait Positioned: BufRead {
        fn at(&self) -> u64;
    }

    impl<R: Read> Positioned for Reader<R> {
        fn at(&self) -> u64 {
            u64::from(self.virtual_position())
        }
    }

    impl<R: Read + Send + 'static> Positioned for noodles_bgzf::io::MultithreadedReader<R> {
        fn at(&self) -> u64 {
            u64::from(self.virtual_position())
        }
    }

    /// What one read gave: the bytes or the error, and the virtual position after it.
    #[derive(Debug, PartialEq)]
    struct Step {
        read: Result<Vec<u8>, (io::ErrorKind, String)>,
        at: u64,
    }

    /// Read `r` to its end with a fixed mix of calls (lines, small and large reads, a
    /// buffer consumed in part), going on after an error for as long as `after_error`
    /// allows, and return every step.
    fn trace(r: &mut dyn Positioned, after_error: usize) -> Vec<Step> {
        let mut steps = vec![Step {
            read: Ok(Vec::new()),
            at: r.at(),
        }];
        let (mut ends, mut errors) = (0, 0);
        for i in 0usize.. {
            let mut got = Vec::new();
            let read = match i % 5 {
                0 => r.read_until(b'\n', &mut got).map(drop),
                1 => {
                    got.resize(100_000, 0);
                    r.read(&mut got).map(|n| got.truncate(n))
                }
                2 => {
                    got.resize(7, 0);
                    r.read(&mut got).map(|n| got.truncate(n))
                }
                3 => r.fill_buf().map(|b| b.len() / 3).map(|n| {
                    got.extend_from_slice(&r.fill_buf().unwrap()[..n]);
                    r.consume(n);
                }),
                _ => {
                    got.resize(70_000, 0);
                    let mut n = 0;
                    let read = loop {
                        match r.read(&mut got[n..]) {
                            Ok(0) => break Ok(()),
                            Ok(k) => n += k,
                            Err(e) => break Err(e),
                        }
                        if n == got.len() {
                            break Ok(());
                        }
                    };
                    got.truncate(n);
                    read
                }
            };
            let failed = read.is_err();
            let empty = got.is_empty();
            steps.push(Step {
                read: read.map(|()| got).map_err(|e| (e.kind(), e.to_string())),
                at: r.at(),
            });
            errors += usize::from(failed);
            // A read that gives nothing is the end unless it is the third kind, which takes
            // a third of what is buffered: see the end twice in a row to be sure.
            ends = if empty && !failed { ends + 1 } else { 0 };
            if ends == 6 || errors > after_error {
                break;
            }
        }
        steps
    }

    /// The traces of noodles' multithreaded reader and of ours at several worker counts.
    fn assert_reads_like_noodles(name: &str, bgzf: &[u8]) -> Vec<Step> {
        let mut theirs = noodles_bgzf::io::MultithreadedReader::with_worker_count(
            workers(2),
            Cursor::new(bgzf.to_vec()),
        );
        let want = trace(&mut theirs, 3);
        for n in [1, 2, 3, 8] {
            let mut r = Reader::new(bgzf, workers(n)).unwrap();
            let got = trace(&mut r, 3);
            if let Some(i) = (0..got.len().max(want.len())).find(|&i| got.get(i) != want.get(i)) {
                let brief = |s: Option<&Step>| {
                    s.map(|s| (s.read.as_ref().map(Vec::len).map_err(Clone::clone), s.at))
                };
                panic!(
                    "{name}, {n} workers: step {i} is {:?}, noodles gives {:?}",
                    brief(got.get(i)),
                    brief(want.get(i))
                );
            }
        }
        want
    }

    fn noodles_block(data: &[u8]) -> Vec<u8> {
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(data).unwrap();
        w.flush().unwrap();
        w.into_inner()
    }

    /// Lines of text as BGZF blocks of several sizes, with blocks that hold no data among
    /// them (the end-of-file block of a file that another was appended to).
    fn stream() -> (Vec<u8>, Vec<u8>) {
        let mut text = Vec::new();
        for i in 0..40_000u32 {
            text.extend_from_slice(
                format!("chr{}\t{}\t{}\n", i % 23, i * 7, i * 7 + 50).as_bytes(),
            );
        }
        let text = [text, self::text(150_000, 21)].concat();
        let mut bgzf = Vec::new();
        let mut rest = &text[..];
        for (i, size) in [1usize, 65_495, 0, 30_000, 65_280, 17, 0, 0, 65_495, 4_000]
            .iter()
            .cycle()
            .enumerate()
        {
            if rest.is_empty() {
                break;
            }
            if *size == 0 {
                bgzf.extend_from_slice(&BGZF_EOF);
                continue;
            }
            let n = (size + i % 3).min(rest.len()).min(BLOCK_DATA);
            bgzf.extend_from_slice(&noodles_block(&rest[..n]));
            rest = &rest[n..];
        }
        (text, bgzf)
    }

    /// The data, and the virtual position after every read, are those of noodles'
    /// multithreaded reader, whatever the worker count and however the stream is read.
    #[test]
    fn reads_what_the_noodles_multithreaded_reader_reads() {
        let (text, mut bgzf) = stream();
        for end in [&b""[..], &BGZF_EOF[..]] {
            bgzf.extend_from_slice(end);
            let steps = assert_reads_like_noodles("a stream of blocks", &bgzf);
            let read: Vec<u8> = steps.into_iter().flat_map(|s| s.read.unwrap()).collect();
            assert!(read == text, "the data read is not the data written");
            let mut all = Vec::new();
            let mut r = Reader::new(&bgzf[..], workers(4)).unwrap();
            r.read_to_end(&mut all).unwrap();
            assert!(all == text);
            assert_eq!(
                r.virtual_position(),
                VirtualPosition::try_from((bgzf.len() as u64, 0)).unwrap()
            );
        }
        for (name, bgzf) in [
            ("an empty stream", Vec::new()),
            ("only the end-of-file block", BGZF_EOF.to_vec()),
            ("one short block", noodles_block(b"chr1\t1\t2\n")),
        ] {
            assert_reads_like_noodles(name, &bgzf);
        }
    }

    /// Where the blocks of `bgzf` start.
    fn starts(bgzf: &[u8]) -> Vec<usize> {
        blocks(bgzf)
            .iter()
            .scan(0, |at, b| {
                *at += b.0;
                Some(*at - b.0)
            })
            .collect()
    }

    /// Streams that are not sound BGZF fail as noodles' multithreaded reader fails on them:
    /// the same data before the error, the same error, and the same reads after it.
    #[test]
    fn bad_streams_fail_as_the_noodles_multithreaded_reader_fails() {
        let (_, mut good) = stream();
        good.extend_from_slice(&BGZF_EOF);
        let starts = starts(&good);
        // A block in the middle that holds data, and what follows it.
        let k = 9;
        let (b, next) = (starts[k], starts[k + 1]);
        assert!(next - b > 1_000);
        let patched = |at: usize, bytes: &[u8]| {
            let mut bad = good.clone();
            bad[at..at + bytes.len()].copy_from_slice(bytes);
            bad
        };
        let flipped = |at: usize| patched(at, &[!good[at]]);
        let gz = {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            e.write_all(&text(200_000, 2)).unwrap();
            e.finish().unwrap()
        };
        // The kind and the text of the first error; a text of "" stands for libdeflate's own.
        type FirstError = Option<(io::ErrorKind, &'static str)>;
        let cases: Vec<(&str, Vec<u8>, FirstError)> = vec![
            (
                "cut inside a block's data",
                good[..b + 500].to_vec(),
                Some((io::ErrorKind::UnexpectedEof, "failed to fill whole buffer")),
            ),
            ("cut inside a block's header", good[..b + 5].to_vec(), None),
            ("cut between two blocks", good[..b].to_vec(), None),
            (
                "cut inside the last block's trailer",
                good[..good.len() - BGZF_EOF.len() - 3].to_vec(),
                Some((io::ErrorKind::UnexpectedEof, "failed to fill whole buffer")),
            ),
            (
                "a header that does not start as gzip",
                patched(b, b"XX"),
                Some((io::ErrorKind::InvalidData, "invalid BGZF header")),
            ),
            (
                "a header without the BGZF field",
                patched(b + 12, b"XY"),
                Some((io::ErrorKind::InvalidData, "invalid BGZF header")),
            ),
            (
                "a wrong checksum",
                flipped(next - 8),
                Some((io::ErrorKind::InvalidData, "block data checksum mismatch")),
            ),
            (
                "an inflated size past 64 KiB",
                patched(next - 4, &70_000u32.to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "invalid BGZF ISIZE")),
            ),
            (
                "an inflated size that is too small",
                patched(next - 4, &100u32.to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "")),
            ),
            (
                "an inflated size that is too large",
                patched(next - 4, &65_536u32.to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "block data checksum mismatch")),
            ),
            (
                "a block size too small for a block",
                patched(b + 16, &10u16.to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "invalid frame size")),
            ),
            (
                "a block size one short of a header and a trailer",
                patched(b + 16, &24u16.to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "invalid frame size")),
            ),
            (
                "a block size of a header and a trailer",
                patched(b + 16, &25u16.to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "")),
            ),
            (
                "a block size that ends inside the block",
                patched(b + 16, &((next - b - 300) as u16).to_le_bytes()),
                Some((io::ErrorKind::InvalidData, "invalid BGZF ISIZE")),
            ),
            (
                "data that is not a deflate stream",
                patched(b + 18, &[0xff; 64]),
                Some((io::ErrorKind::InvalidData, "")),
            ),
            (
                "plain gzip after BGZF blocks",
                [&good[..b], &gz[..]].concat(),
                Some((io::ErrorKind::InvalidData, "invalid BGZF header")),
            ),
            (
                "text after the end-of-file block",
                [&good[..], &text(200_000, 4)[..]].concat(),
                Some((io::ErrorKind::InvalidData, "invalid BGZF header")),
            ),
        ];
        for (name, bad, first_error) in cases {
            let steps = assert_reads_like_noodles(name, &bad);
            let failed = steps.iter().find_map(|s| s.read.as_ref().err());
            let failed = failed.map(|(kind, msg)| match first_error {
                Some((_, "")) => (*kind, ""),
                _ => (*kind, msg.as_str()),
            });
            assert_eq!(failed, first_error, "{name}");
        }
    }

    /// A source that fails is reported after the blocks read before it failed, and then the
    /// stream has ended.
    #[test]
    fn a_source_that_fails_is_reported_in_its_place() {
        struct Failing<'a>(&'a [u8]);
        impl Read for Failing<'_> {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.0.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::BrokenPipe,
                        "the source failed",
                    ));
                }
                self.0.read(buf)
            }
        }
        let (text, bgzf) = stream();
        // The source fails 40 bytes into a block that is larger than that.
        let sizes = blocks(&bgzf);
        let k = (20..).find(|&k| sizes[k].0 > 1_000).unwrap();
        let cut = starts(&bgzf)[k] + 40;
        let whole: usize = sizes[..k].iter().map(|b| b.1).sum();
        for n in [1, 4] {
            let mut r = Reader::new(Failing(&bgzf[..cut]), workers(n)).unwrap();
            let mut read = Vec::new();
            let e = r.read_to_end(&mut read).unwrap_err();
            assert_eq!(
                (e.kind(), e.to_string().as_str()),
                (io::ErrorKind::BrokenPipe, "the source failed")
            );
            assert!(read == text[..whole]);
            assert_eq!(r.read(&mut [0; 10]).unwrap(), 0);
        }
    }

    /// The reader takes from its source what its workers and its buffer hold, not the whole
    /// stream, however little of the data is read.
    #[test]
    fn the_source_is_read_a_bounded_way_ahead() {
        struct Counting<'a>(&'a [u8], Arc<Mutex<usize>>);
        impl Read for Counting<'_> {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                let n = self.0.read(buf)?;
                *self.1.lock().unwrap() += n;
                Ok(n)
            }
        }
        let bgzf = noodles(&text(200 * BLOCK_DATA, 13), 1, 1, 1 << 20, &[]);
        let stored = blocks(&bgzf)[0].0;
        let taken = Arc::new(Mutex::new(0));
        let mut r = Reader::new(Counting(&bgzf, taken.clone()), workers(2)).unwrap();
        for block in 1..=20 {
            r.read_exact(&mut vec![0; BLOCK_DATA]).unwrap();
            // The blocks read, the 2 x 2 in flight, and one buffer of the source.
            let most = (block + 4) * (stored + 64) + SOURCE_BUF;
            let taken = *taken.lock().unwrap();
            assert!(
                taken <= most,
                "after {block} blocks, {taken} bytes of the source are read"
            );
        }
    }
}
