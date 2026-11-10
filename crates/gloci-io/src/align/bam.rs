//! The BAM reader for interval commands (D101, D102, D115).
//!
//! Records are decoded by hand from the BGZF-decoded input, so the thread count of
//! `Input` applies. Each record kept is a footprint plus its raw bytes (slim unless
//! the command writes alignments); unmapped reads and reads the filters drop are
//! counted and not returned.

use super::{AlignOptions, Aln, AlnHeader};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::Input;
use crate::record::{Header, RecordBatch, RecordMeta, RecordReader, intern_chrom, warn_undeclared};
use gloci_core::{ChromId, ChromTable};
use std::io::{self, BufRead};
use std::sync::Arc;

const FLAG_UNMAPPED: u16 = 0x4;

/// Where the records of one reference are counted from: "record N of CHROM".
pub struct BamReader {
    name: String,
    compression: Compression,
    header: Header,
    aln: Arc<AlnHeader>,
    diag: Diagnostics,
    policy: InvalidPolicy,
    opts: AlignOptions,
    src: Box<dyn BufRead + Send>,
    /// Records started so far in the file (`RecordMeta::loc`).
    rec: u64,
    /// The `refID` of the current run of records and how many of the run were read:
    /// the number errors name, so every thread count reports one place.
    run_ref: i32,
    run_seen: u64,
    eof: bool,
    /// Reference id to chromosome id, filled as references are first met.
    ids: Vec<Option<ChromId>>,
    scratch: Vec<u8>,
    clip_warned: bool,
}

/// Read until `buf` is full or input ends; returns the bytes read.
fn read_full(r: &mut dyn BufRead, out: &mut Vec<u8>, n: usize) -> io::Result<usize> {
    crate::bcf::append(r, out, n)
}

impl BamReader {
    pub(crate) fn open_bam(
        input: Input,
        policy: InvalidPolicy,
        map: Option<&ChromMap>,
        opts: AlignOptions,
    ) -> Result<BamReader, IoError> {
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
        let aln = Arc::new(h);
        Ok(BamReader {
            header: Header::Alignment(aln.clone()),
            ids: vec![None; aln.refs.len()],
            aln,
            name,
            compression,
            diag: Diagnostics::default(),
            policy,
            opts,
            src: reader,
            rec: 0,
            run_ref: i32::MIN,
            run_seen: 0,
            eof: false,
            scratch: Vec::new(),
            clip_warned: false,
        })
    }

    /// "record N of CHROM" for the record that is `seen` in its run, on `ref_id`.
    fn place(&self, ref_id: i32, seen: u64) -> String {
        let chrom = match usize::try_from(ref_id) {
            Ok(i) => match self.aln.refs.get(i) {
                Some((n, _)) => String::from_utf8_lossy(n).into_owned(),
                None => format!("refID {ref_id}"),
            },
            Err(_) => "unplaced reads".to_string(),
        };
        format!("record {} of {chrom}", gloci_core::thousands(seen))
    }

    fn data_err(&self, place: &str, msg: &str) -> IoError {
        IoError::Detect {
            path: self.name.clone(),
            msg: format!("{} {place}: {msg}", self.name),
        }
    }

    fn read_err(&self, source: io::Error) -> IoError {
        IoError::Read {
            path: self.name.clone(),
            source,
        }
    }

    /// An invalid record: an error, or counted under `--skip-invalid`.
    fn bad(&mut self, place: &str, msg: &str) -> Result<(), IoError> {
        match self.policy {
            InvalidPolicy::Fail => Err(self.data_err(place, msg)),
            InvalidPolicy::Skip => {
                self.diag.record_place(place, msg);
                Ok(())
            }
        }
    }

    fn dropped_by_filters(&self, flag: u16, mapq: u8) -> bool {
        mapq < self.opts.min_mapq
            || flag & self.opts.require_flags != self.opts.require_flags
            || flag & self.opts.exclude_flags != 0
    }
}

impl RecordReader for BamReader {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        Format::Bam
    }
    fn compression(&self) -> Compression {
        self.compression
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        debug_assert!(max > 0, "read_batch max must be > 0");
        let hdr = self.aln.clone();
        let mut n = 0;
        while n < max {
            if self.eof {
                return Ok(false);
            }
            let start = batch.buf.len();
            let got = read_full(&mut self.src, &mut batch.buf, 4).map_err(|e| self.read_err(e))?;
            if got == 0 {
                self.eof = true;
                return Ok(false);
            }
            self.rec += 1;
            let truncated = |this: &Self, buf: &mut Vec<u8>| {
                buf.truncate(start);
                this.data_err(
                    &this.place(this.run_ref, this.run_seen + 1),
                    crate::bcf::RECORD_TRUNCATED,
                )
            };
            if got < 4 {
                return Err(truncated(self, &mut batch.buf));
            }
            let block = i32::from_le_bytes([
                batch.buf[start],
                batch.buf[start + 1],
                batch.buf[start + 2],
                batch.buf[start + 3],
            ]);
            batch.buf.truncate(start);
            let size = match usize::try_from(block) {
                Ok(s) => s,
                Err(_) => {
                    let place = self.place(self.run_ref, self.run_seen + 1);
                    return Err(self.data_err(
                        &place,
                        &format!(
                            "block_size {block} is negative, so the next record cannot be found; \
the file is corrupt"
                        ),
                    ));
                }
            };
            let body =
                read_full(&mut self.src, &mut batch.buf, size).map_err(|e| self.read_err(e))?;
            if body < size {
                return Err(truncated(self, &mut batch.buf));
            }
            // The record's place: its number among the records of its reference.
            let ref_id = if size >= 4 {
                i32::from_le_bytes([
                    batch.buf[start],
                    batch.buf[start + 1],
                    batch.buf[start + 2],
                    batch.buf[start + 3],
                ])
            } else {
                self.run_ref
            };
            if ref_id != self.run_ref {
                self.run_ref = ref_id;
                self.run_seen = 0;
            }
            self.run_seen += 1;
            let place = self.place(ref_id, self.run_seen);
            if size < 32 {
                batch.buf.truncate(start);
                let msg = format!(
                    "block_size {size} is below the 32 bytes of fixed fields of a BAM record; \
the file is corrupt"
                );
                self.bad(&place, &msg)?;
                continue;
            }
            let a = match Aln::parse(&batch.buf[start..], &hdr) {
                Ok(a) => a,
                Err(msg) => {
                    batch.buf.truncate(start);
                    self.bad(&place, &msg)?;
                    continue;
                }
            };
            let flag = a.flag();
            if flag & FLAG_UNMAPPED != 0 {
                self.diag.unmapped += 1;
                batch.buf.truncate(start);
                continue;
            }
            if self.dropped_by_filters(flag, a.mapq()) {
                self.diag.filtered += 1;
                batch.buf.truncate(start);
                continue;
            }
            // `parse` guarantees a mapped read has a reference and a position.
            let (Some(rid), Some(pos)) = (a.ref_id(), a.pos()) else {
                batch.buf.truncate(start);
                continue;
            };
            let (span, clipped) = match a.span(hdr.refs[rid].1) {
                Some(s) => s,
                None => {
                    batch.buf.truncate(start);
                    continue;
                }
            };
            if clipped && !self.clip_warned {
                self.clip_warned = true;
                log::warn!(
                    "{}: read {} on {} extends past the end of the reference ({} bases); \
reads like it are cut at the reference length",
                    self.name,
                    String::from_utf8_lossy(a.name()),
                    String::from_utf8_lossy(&hdr.refs[rid].0),
                    hdr.refs[rid].1
                );
            }
            let strand = a.strand();
            let chrom = match self.ids[rid] {
                Some(c) => c,
                None => match intern_chrom(chroms, &hdr.refs[rid].0) {
                    Ok(c) => {
                        self.ids[rid] = Some(c);
                        c
                    }
                    Err(msg) => {
                        batch.buf.truncate(start);
                        self.bad(&place, &msg)?;
                        continue;
                    }
                },
            };
            if chroms.has_undeclared() {
                warn_undeclared(chroms, &self.name, Some(Loc::Record(self.rec)));
            }
            if !self.opts.full_records {
                self.scratch.clear();
                a.write_slim(self.opts.score_tag, &mut self.scratch);
                batch.buf.truncate(start);
                batch.buf.extend_from_slice(&self.scratch);
            }
            let len = batch.buf.len() - start;
            batch.recs.push(RecordMeta {
                chrom,
                key: pos,
                span,
                strand,
                loc: self.rec,
                off: start,
                len: len as u32,
            });
            n += 1;
        }
        Ok(true)
    }
}
