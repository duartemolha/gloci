//! The BAM reader for interval commands (D101, D102, D115).
//!
//! Records are decoded by hand from the BGZF-decoded input, so the thread count of
//! `Input` applies. Each record kept is a footprint plus its raw bytes (slim unless
//! the command writes alignments); unmapped reads and reads the filters drop are
//! counted and not returned.

use super::place::{IntervalRules, Verdict};
use super::{AlignOptions, Aln, AlnHeader, AlnRef};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::Input;
use crate::record::{Header, RecordBatch, RecordMeta, RecordReader, warn_undeclared};
use gloci_core::ChromTable;
use std::io::{self, BufRead};
use std::sync::Arc;

/// Read until `buf` is full or input ends; returns the bytes read.
fn read_full(r: &mut dyn BufRead, out: &mut Vec<u8>, n: usize) -> io::Result<usize> {
    crate::bcf::append(r, out, n)
}

/// The framing of a BAM file: the header, then records cut at their `block_size`. Shared by
/// the interval reader and [`super::source::BamSource`], so both count and name records alike.
pub(crate) struct BamFrames {
    pub(crate) name: String,
    pub(crate) compression: Compression,
    pub(crate) aln: Arc<AlnHeader>,
    src: Box<dyn BufRead + Send>,
    /// Records started so far in the file (`RecordMeta::loc`).
    pub(crate) rec: u64,
    /// The `refID` of the current run of records and how many of the run were read:
    /// the "record N of CHROM" errors name, so every thread count reports one place.
    run_ref: i32,
    run_seen: u64,
    eof: bool,
}

impl BamFrames {
    pub(crate) fn open(input: Input, map: Option<&ChromMap>) -> Result<BamFrames, IoError> {
        let Input {
            name,
            compression,
            mut reader,
            ..
        } = input;
        let detect = |msg: String| IoError::Detect {
            msg: format!("{name}: {msg}"),
            path: name.clone(),
        };
        let mut h = AlnHeader::read_bam(&mut reader).map_err(detect)?;
        if let Some(m) = map {
            m.check_file_names(&name, h.refs.iter().map(|r| &r.0[..]))?;
            h = h.renamed(m);
        }
        Ok(BamFrames {
            aln: Arc::new(h),
            name,
            compression,
            src: reader,
            rec: 0,
            run_ref: i32::MIN,
            run_seen: 0,
            eof: false,
        })
    }

    /// "record N of CHROM (record M of the file)" for the record that is `seen` in its run,
    /// on `ref_id`. `self.rec` is the file-wide number, the same at every thread count.
    fn place_of(&self, ref_id: i32, seen: u64) -> String {
        let chrom = match usize::try_from(ref_id) {
            Ok(i) => match self.aln.refs.get(i) {
                Some((n, _)) => String::from_utf8_lossy(n).into_owned(),
                None => format!("refID {ref_id}"),
            },
            Err(_) if ref_id == -1 => "unplaced reads".to_string(),
            Err(_) => format!("refID {ref_id}"),
        };
        format!(
            "record {} of {chrom} (record {} of the file)",
            gloci_core::thousands(seen),
            gloci_core::thousands(self.rec)
        )
    }

    /// The place of the record last returned by [`BamFrames::next`] (or being read).
    pub(crate) fn place(&self) -> String {
        self.place_of(self.run_ref, self.run_seen)
    }

    pub(crate) fn data_err(&self, msg: &str) -> IoError {
        IoError::Detect {
            path: self.name.clone(),
            msg: format!("{} {}: {msg}", self.name, self.place()),
        }
    }

    pub(crate) fn read_err(&self, source: io::Error) -> IoError {
        IoError::Read {
            path: self.name.clone(),
            source,
        }
    }

    /// An invalid record: an error, or counted in `diag` under `--skip-invalid`.
    pub(crate) fn invalid(
        &self,
        policy: InvalidPolicy,
        diag: &mut Diagnostics,
        msg: &str,
    ) -> Result<(), IoError> {
        match policy {
            InvalidPolicy::Fail => Err(self.data_err(msg)),
            InvalidPolicy::Skip => {
                diag.record_place(&self.place(), msg);
                Ok(())
            }
        }
    }

    /// Append the next record (without its `block_size`) to `buf`. `Ok(false)` at the end of
    /// the file. The record is not checked beyond its frame.
    pub(crate) fn next(&mut self, buf: &mut Vec<u8>) -> Result<bool, IoError> {
        if self.eof {
            return Ok(false);
        }
        let start = buf.len();
        let got = read_full(&mut self.src, buf, 4).map_err(|e| self.read_err(e))?;
        if got == 0 {
            self.eof = true;
            return Ok(false);
        }
        self.rec += 1;
        // Until the record's reference is known, errors name it as the next of the last run.
        let truncated = |this: &mut Self, buf: &mut Vec<u8>| {
            buf.truncate(start);
            this.run_seen += 1;
            this.data_err(crate::bcf::RECORD_TRUNCATED)
        };
        if got < 4 {
            return Err(truncated(self, buf));
        }
        let block =
            i32::from_le_bytes([buf[start], buf[start + 1], buf[start + 2], buf[start + 3]]);
        buf.truncate(start);
        let Ok(size) = usize::try_from(block) else {
            self.run_seen += 1;
            return Err(self.data_err(&format!(
                "block_size {block} is negative, so the next record cannot be found; \
the file is corrupt"
            )));
        };
        let body = read_full(&mut self.src, buf, size).map_err(|e| self.read_err(e))?;
        if body < size {
            return Err(truncated(self, buf));
        }
        // The record's place: its number among the records of its reference.
        let ref_id = if size >= 4 {
            i32::from_le_bytes([buf[start], buf[start + 1], buf[start + 2], buf[start + 3]])
        } else {
            self.run_ref
        };
        if ref_id != self.run_ref {
            self.run_ref = ref_id;
            self.run_seen = 0;
        }
        self.run_seen += 1;
        Ok(true)
    }
}

/// Reads BAM records for interval commands.
pub struct BamReader {
    frames: BamFrames,
    header: Header,
    diag: Diagnostics,
    policy: InvalidPolicy,
    opts: AlignOptions,
    rules: IntervalRules,
    scratch: Vec<u8>,
}

impl BamReader {
    pub(crate) fn open_bam(
        input: Input,
        policy: InvalidPolicy,
        map: Option<&ChromMap>,
        opts: AlignOptions,
    ) -> Result<BamReader, IoError> {
        let frames = BamFrames::open(input, map)?;
        Ok(BamReader {
            header: Header::Alignment(frames.aln.clone()),
            rules: IntervalRules::new(frames.aln.refs.len()),
            frames,
            diag: Diagnostics::default(),
            policy,
            opts,
            scratch: Vec::new(),
        })
    }

    /// An invalid record: an error, or counted under `--skip-invalid`.
    fn bad(&mut self, msg: &str) -> Result<(), IoError> {
        self.frames.invalid(self.policy, &mut self.diag, msg)
    }
}

impl RecordReader for BamReader {
    fn name(&self) -> &str {
        &self.frames.name
    }
    fn format(&self) -> Format {
        Format::Bam
    }
    fn compression(&self) -> Compression {
        self.frames.compression
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }
    fn whole_records(&self) -> bool {
        self.opts.full_records
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        debug_assert!(max > 0, "read_batch max must be > 0");
        let hdr = self.frames.aln.clone();
        let mut n = 0;
        while n < max {
            let start = batch.buf.len();
            if !self.frames.next(&mut batch.buf)? {
                return Ok(false);
            }
            if batch.buf.len() - start < 32 {
                let msg = format!(
                    "block_size {} is below the 32 bytes of fixed fields of a BAM record; \
the file is corrupt",
                    batch.buf.len() - start
                );
                batch.buf.truncate(start);
                self.bad(&msg)?;
                continue;
            }
            let a = match Aln::parse(&batch.buf[start..], &hdr) {
                Ok(a) => a,
                Err(msg) => {
                    batch.buf.truncate(start);
                    self.bad(&msg)?;
                    continue;
                }
            };
            let verdict = self.rules.judge(
                &self.frames.name,
                &AlnRef::Bam(a),
                &hdr,
                &self.opts,
                chroms,
                &mut self.diag,
            );
            let p = match verdict {
                Verdict::Keep(p) => p,
                Verdict::Skip => {
                    batch.buf.truncate(start);
                    continue;
                }
                Verdict::Bad(msg) => {
                    batch.buf.truncate(start);
                    self.bad(&msg)?;
                    continue;
                }
            };
            if chroms.has_undeclared() {
                warn_undeclared(
                    chroms,
                    &self.frames.name,
                    Some(Loc::Record(self.frames.rec)),
                );
            }
            if !self.opts.full_records {
                self.scratch.clear();
                a.write_slim(self.opts.score_tag, &mut self.scratch);
                batch.buf.truncate(start);
                batch.buf.extend_from_slice(&self.scratch);
            }
            let len = batch.buf.len() - start;
            batch.recs.push(RecordMeta {
                chrom: p.chrom,
                key: p.pos,
                span: p.span,
                strand: p.strand,
                loc: self.frames.rec,
                off: start,
                len: len as u32,
            });
            n += 1;
        }
        Ok(true)
    }
}
