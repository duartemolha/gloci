//! ChunkSink: in-memory record buffer that spills to LZ4 in tmp_dir (D41).

use gloci_io::{IoError, RecordWriter, TempGuard, create_file_in_temp_dir, create_temp_dir_in};
use lz4_flex::frame::{BlockMode, BlockSize, FrameDecoder, FrameEncoder, FrameInfo};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const DEFAULT_MEM_LIMIT: usize = 64 * 1024 * 1024; // 64 MiB

struct SpillState {
    _dir: TempDir,
    _guard: TempGuard,
    enc: Option<FrameEncoder<BufWriter<File>>>,
    path: PathBuf,
}

pub struct ChunkSink {
    mem_limit: usize,
    tmp_dir: PathBuf,
    mem_buf: Vec<u8>,
    spill: Option<SpillState>,
    is_bcf: bool,
}

impl ChunkSink {
    /// Creates a sink buffering up to 64 MiB (or env `GLOCI_CHUNK_MEM_LIMIT`)
    /// in memory before spilling to LZ4 in `tmp_dir`.
    pub fn new(tmp_dir: &Path) -> Self {
        let limit = std::env::var("GLOCI_CHUNK_MEM_LIMIT")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(DEFAULT_MEM_LIMIT);
        Self::with_mem_limit(tmp_dir, limit)
    }

    /// Creates a sink with a custom in-memory limit, for testing.
    pub fn with_mem_limit(tmp_dir: &Path, mem_limit: usize) -> Self {
        Self {
            mem_limit,
            tmp_dir: tmp_dir.to_path_buf(),
            mem_buf: Vec::new(),
            spill: None,
            is_bcf: false,
        }
    }

    /// Drain all buffered records to `out` in order, then remove any spill files.
    pub fn drain_to(&mut self, out: &mut dyn RecordWriter) -> Result<(), IoError> {
        let mut offset = 0;
        while offset + 4 <= self.mem_buf.len() {
            let len =
                u32::from_le_bytes(self.mem_buf[offset..offset + 4].try_into().unwrap()) as usize;
            offset += 4;
            let rec = &self.mem_buf[offset..offset + len];
            offset += len;
            out.write_record(rec).map_err(|source| IoError::Write {
                path: "output".to_string(),
                source,
            })?;
        }
        self.mem_buf.clear();

        if let Some(mut spill) = self.spill.take() {
            if let Some(enc) = spill.enc.take() {
                let mut buf_writer = enc.finish().map_err(|source| IoError::Write {
                    path: spill.path.display().to_string(),
                    source: source.into(),
                })?;
                buf_writer.flush().map_err(|source| IoError::Write {
                    path: spill.path.display().to_string(),
                    source,
                })?;
            }

            let file = File::open(&spill.path).map_err(|source| IoError::Open {
                path: spill.path.display().to_string(),
                source,
            })?;
            let mut dec = FrameDecoder::new(BufReader::with_capacity(128 << 10, file));
            let mut len_buf = [0u8; 4];
            let mut rec_buf = Vec::new();

            loop {
                let mut n = 0;
                while n < 4 {
                    match dec.read(&mut len_buf[n..]) {
                        Ok(0) => break,
                        Ok(k) => n += k,
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                        Err(e) => {
                            return Err(IoError::Read {
                                path: spill.path.display().to_string(),
                                source: e,
                            });
                        }
                    }
                }
                if n == 0 {
                    break;
                }
                if n < 4 {
                    return Err(IoError::Read {
                        path: spill.path.display().to_string(),
                        source: io::Error::new(
                            io::ErrorKind::UnexpectedEof,
                            "truncated record length in spill file",
                        ),
                    });
                }
                let len = u32::from_le_bytes(len_buf) as usize;
                rec_buf.resize(len, 0);
                dec.read_exact(&mut rec_buf)
                    .map_err(|source| IoError::Read {
                        path: spill.path.display().to_string(),
                        source,
                    })?;
                out.write_record(&rec_buf)
                    .map_err(|source| IoError::Write {
                        path: "output".to_string(),
                        source,
                    })?;
            }
        }

        Ok(())
    }
}

impl RecordWriter for ChunkSink {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        let Ok(rec_len) = u32::try_from(raw.len()) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "record exceeds 4 GiB",
            ));
        };

        if self.spill.is_none() {
            if self.mem_buf.len() + 4 + raw.len() <= self.mem_limit {
                self.mem_buf.extend_from_slice(&rec_len.to_le_bytes());
                self.mem_buf.extend_from_slice(raw);
                return Ok(());
            }

            let (dir, guard) =
                create_temp_dir_in(&self.tmp_dir, "gloci-chunk-").map_err(|source| {
                    io::Error::new(
                        source.kind(),
                        format!(
                            "failed to create temp dir in {}: {source}",
                            self.tmp_dir.display()
                        ),
                    )
                })?;
            let path = dir.path().join("chunk-spill.lz4");
            let file = create_file_in_temp_dir(&path).map_err(|source| {
                io::Error::new(
                    source.kind(),
                    format!("failed to create spill file {}: {source}", path.display()),
                )
            })?;
            let info = FrameInfo::new()
                .block_size(BlockSize::Max64KB)
                .block_mode(BlockMode::Independent);
            let mut enc =
                FrameEncoder::with_frame_info(info, BufWriter::with_capacity(128 << 10, file));
            enc.write_all(&rec_len.to_le_bytes())?;
            enc.write_all(raw)?;
            self.spill = Some(SpillState {
                _dir: dir,
                _guard: guard,
                enc: Some(enc),
                path,
            });
            return Ok(());
        }

        let spill = self.spill.as_mut().unwrap();
        let enc = spill.enc.as_mut().unwrap();
        enc.write_all(&rec_len.to_le_bytes())?;
        enc.write_all(raw)?;
        Ok(())
    }

    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }

    fn is_bcf(&self) -> bool {
        self.is_bcf
    }
}
