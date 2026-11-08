//! Sorted runs on disk: LZ4 frames of `(chrom u32, key u64, end u64,
//! seq u64, len u32, bytes)` records in a registered temp directory (D24).

use gloci_io::{IoError, TempGuard, create_registered_guarded, create_unregistered};
use lz4_flex::frame::{BlockMode, BlockSize, FrameDecoder, FrameEncoder, FrameInfo};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// A record's sort fields as stored in a run (the chromosome as its
/// `ChromId`, not its rank, since ranks can change as chromosomes appear).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct RecHead {
    pub chrom: u32,
    pub key: u64,
    pub end: u64,
    pub seq: u64,
}

const HEAD_BYTES: usize = 4 + 8 + 8 + 8 + 4;

/// LZ4 block size: the smallest, so each open run costs little memory.
const BLOCK: BlockSize = BlockSize::Max64KB;

/// The spill directory, removed when dropped (and, through the cleanup
/// registry, on SIGINT/SIGTERM/SIGHUP or a panic).
pub(crate) struct SpillDir {
    // Dropped in this order: the directory is removed, then unregistered.
    dir: TempDir,
    _guard: TempGuard,
    next: u64,
}

impl SpillDir {
    pub fn create(parent: &Path) -> Result<SpillDir, IoError> {
        // Created and registered under one lock, so a signal can't land between.
        let (dir, guard) = create_registered_guarded(|| {
            let dir = tempfile::Builder::new()
                .prefix("gloci-sort-")
                .tempdir_in(parent)?;
            let path = dir.path().to_path_buf();
            Ok((dir, path))
        })
        .map_err(|source| IoError::Write {
            path: parent.display().to_string(),
            source,
        })?;
        Ok(SpillDir {
            dir,
            _guard: guard,
            next: 0,
        })
    }

    /// Start writing a new run file.
    pub fn create_run(&mut self) -> Result<RunWriter, IoError> {
        let path = self.dir.path().join(format!("run-{:06}.lz4", self.next));
        self.next += 1;
        // Under the registry lock: the termination path removes the directory
        // while holding it, so no run file can appear mid-removal.
        let file = create_unregistered(|| File::create(&path))
            .map_err(|source| write_err(&path, source))?;
        let info = FrameInfo::new()
            .block_size(BLOCK)
            .block_mode(BlockMode::Independent);
        Ok(RunWriter {
            enc: FrameEncoder::with_frame_info(info, BufWriter::with_capacity(128 << 10, file)),
            path,
        })
    }

    /// Remove the directory now.
    pub fn close(self) {
        let SpillDir { dir, _guard, .. } = self;
        // Best effort, like TempDir's own drop.
        let _ = dir.close();
    }
}

fn write_err(path: &Path, source: io::Error) -> IoError {
    IoError::Write {
        path: path.display().to_string(),
        source,
    }
}

fn read_err(path: &Path, source: io::Error) -> IoError {
    IoError::Read {
        path: path.display().to_string(),
        source,
    }
}

/// A finished run file.
#[derive(Debug)]
pub(crate) struct Run {
    path: PathBuf,
}

impl Run {
    /// Delete the file (best effort; the directory goes at the end anyway).
    pub fn remove(self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

pub(crate) struct RunWriter {
    enc: FrameEncoder<BufWriter<File>>,
    path: PathBuf,
}

impl RunWriter {
    pub fn write(&mut self, h: &RecHead, raw: &[u8]) -> Result<(), IoError> {
        let mut buf = [0u8; HEAD_BYTES];
        buf[0..4].copy_from_slice(&h.chrom.to_le_bytes());
        buf[4..12].copy_from_slice(&h.key.to_le_bytes());
        buf[12..20].copy_from_slice(&h.end.to_le_bytes());
        buf[20..28].copy_from_slice(&h.seq.to_le_bytes());
        // Record lengths come from `RecordMeta::len`, a u32.
        buf[28..32].copy_from_slice(&(raw.len() as u32).to_le_bytes());
        self.enc
            .write_all(&buf)
            .and_then(|()| self.enc.write_all(raw))
            .map_err(|e| write_err(&self.path, e))
    }

    pub fn finish(self) -> Result<Run, IoError> {
        let RunWriter { enc, path } = self;
        let flushed = enc
            .finish()
            .map_err(io::Error::from)
            .and_then(|mut w| w.flush());
        flushed.map_err(|e| write_err(&path, e))?;
        Ok(Run { path })
    }
}

pub(crate) struct RunReader {
    dec: FrameDecoder<BufReader<File>>,
    path: PathBuf,
}

impl RunReader {
    pub fn open(run: &Run) -> Result<RunReader, IoError> {
        let file = File::open(&run.path).map_err(|source| IoError::Open {
            path: run.path.display().to_string(),
            source,
        })?;
        Ok(RunReader {
            dec: FrameDecoder::new(BufReader::with_capacity(8 << 10, file)),
            path: run.path.clone(),
        })
    }

    /// Read the next record into `head` and `bytes`; false at the end of the
    /// run.
    pub fn next(&mut self, head: &mut RecHead, bytes: &mut Vec<u8>) -> Result<bool, IoError> {
        let mut buf = [0u8; HEAD_BYTES];
        let mut got = 0;
        while got < HEAD_BYTES {
            match self.dec.read(&mut buf[got..]) {
                Ok(0) => break,
                Ok(n) => got += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(read_err(&self.path, e)),
            }
        }
        if got == 0 {
            return Ok(false);
        }
        if got < HEAD_BYTES {
            let e = io::Error::new(io::ErrorKind::UnexpectedEof, "spill file is truncated");
            return Err(read_err(&self.path, e));
        }
        let u32_at = |i: usize| u32::from_le_bytes(buf[i..i + 4].try_into().unwrap());
        let u64_at = |i: usize| u64::from_le_bytes(buf[i..i + 8].try_into().unwrap());
        *head = RecHead {
            chrom: u32_at(0),
            key: u64_at(4),
            end: u64_at(12),
            seq: u64_at(20),
        };
        bytes.clear();
        bytes.resize(u32_at(28) as usize, 0);
        self.dec
            .read_exact(bytes)
            .map_err(|e| read_err(&self.path, e))?;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_round_trip() {
        let tmp = tempfile::tempdir().unwrap();
        let mut dir = SpillDir::create(tmp.path()).unwrap();
        let recs: Vec<(RecHead, Vec<u8>)> = (0..50_000u64)
            .map(|i| {
                let h = RecHead {
                    chrom: (i % 7) as u32,
                    key: i * 3,
                    end: u64::MAX - i,
                    seq: i,
                };
                (h, format!("line {i}").into_bytes())
            })
            .chain([(RecHead::default(), vec![b'z'; 300_000])])
            .collect();
        let mut w = dir.create_run().unwrap();
        for (h, raw) in &recs {
            w.write(h, raw).unwrap();
        }
        let run = w.finish().unwrap();
        let mut r = RunReader::open(&run).unwrap();
        let (mut h, mut raw) = (RecHead::default(), Vec::new());
        for (want_h, want_raw) in &recs {
            assert!(r.next(&mut h, &mut raw).unwrap());
            assert_eq!((&h, &raw), (want_h, want_raw));
        }
        assert!(!r.next(&mut h, &mut raw).unwrap());
        let path = tmp
            .path()
            .read_dir()
            .unwrap()
            .next()
            .unwrap()
            .unwrap()
            .path();
        dir.close();
        assert!(!path.exists());
    }
}
