//! CRAM input (D109, D110): the file's framing, the reference, the decode pool, and the records
//! transcoded to BAM bytes.
//!
//! - **Framing** is gloci's own ([`frame`]): the file definition, the header container and each data
//!   container are read in order and CRC-checked, and every slice header is checked (landmarks,
//!   block counts, reference ids, spans) before noodles-cram decodes anything, so its `expect`s on
//!   a bad reference id or an empty repository stay unreachable. A file that ends without the
//!   end-of-file container is an error naming it (`#cram-eof-container`).
//! - **The reference** is `--reference FASTA` only (bedtools `--cram-ref`): no `REF_PATH`, no
//!   cache directory, no network. Each `@SQ` the `.fai` names must have its `LN`, checked before
//!   any record. [`RefAdapter`](reference::RefAdapter) reads whole sequences through
//!   [`IndexedFasta`](crate::fasta::IndexedFasta), uppercased as
//!   htslib reads them, and never answers `None`: a missing name is `Some(Err)` with the D109
//!   message.
//! - **Each decoder** owns a noodles repository and clears it when the slice's reference id
//!   changes, so it holds one sequence (a slice holding several references holds those).
//! - **The pool** (D110): containers are read in order, decoded on up to `threads` workers with
//!   at most `2 × threads` in flight, and released in order, so the output does not depend on
//!   the thread count.
//! - **Transcoding**: each record is encoded as BAM bytes once, with htslib's `bin` and MAPQ 0
//!   for unmapped reads (CRAM does not store their MAPQ), and the bytes go through the BAM
//!   reader's framing ([`BamFrames`]), so interval commands and the every-record source treat a
//!   CRAM exactly as the BAM it encodes. MD and NM are not regenerated (`#cram-md-nm`).

mod decode;
mod frame;
mod reference;

pub(crate) use decode::ExtentDecoder;
pub(crate) use frame::{CramFrames, RawContainer, slice_heads};

use super::AlignOptions;
use super::bam::{BamFrames, BamReader};
use super::header::AlnHeader;
use super::source::{AlnSource, BamSource};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::Diagnostics;
use crate::error::IoError;
use crate::input::{Input, is_regular_file};
use crate::pool::Ordered;
use crate::record::{Header, ReaderOptions, RecordBatch, RecordReader};
use decode::Decoder;
use frame::{CRAM_EOF, CramHead, detect_err, read_head, take_into, truncated};
use gloci_core::ChromTable;
use reference::{CramInfo, RefState};
use std::io::{self, BufRead, Read};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::Ordering;

fn decode_job(d: &mut Decoder, c: RawContainer) -> Result<Vec<u8>, IoError> {
    let mut out = Vec::new();
    d.decode(&c, &mut out)?;
    Ok(out)
}

enum Mode {
    Inline(Box<Decoder>),
    Pool(Ordered<RawContainer, Vec<u8>>, u64),
}

/// The decoded records of a CRAM as the byte stream of a BAM body: each record with its
/// `block_size`, in file order. A decoding failure is an `io::Error` holding the gloci error
/// ([`unwrap_failure`] takes it out).
struct CramStream {
    frames: CramFrames,
    mode: Mode,
    cur: Vec<u8>,
    at: usize,
    reading_done: bool,
    failed: bool,
    /// Each decoder's reference counters (the bounded-cache test reads them).
    states: Vec<Arc<RefState>>,
}

impl CramStream {
    fn new(
        frames: CramFrames,
        info: Arc<CramInfo>,
        threads: NonZeroUsize,
    ) -> Result<Self, IoError> {
        let n = threads.get();
        let (mode, states) = if n == 1 {
            let d = Decoder::new(info)?;
            let states = vec![d.state.clone()];
            (Mode::Inline(Box::new(d)), states)
        } else {
            let decoders = (0..n)
                .map(|_| Decoder::new(info.clone()))
                .collect::<Result<Vec<_>, _>>()?;
            let states = decoders.iter().map(|d| d.state.clone()).collect();
            (
                Mode::Pool(Ordered::new(decoders, decode_job)?, 2 * n as u64),
                states,
            )
        };
        Ok(CramStream {
            frames,
            mode,
            cur: Vec::new(),
            at: 0,
            reading_done: false,
            failed: false,
            states,
        })
    }

    /// The next decoded container, in file order.
    fn next_container(&mut self) -> Result<Option<Vec<u8>>, IoError> {
        match &mut self.mode {
            Mode::Inline(d) => match self.frames.next()? {
                None => Ok(None),
                Some(c) => {
                    let mut out = std::mem::take(&mut self.cur);
                    out.clear();
                    d.decode(&c, &mut out)?;
                    Ok(Some(out))
                }
            },
            Mode::Pool(pool, cap) => {
                while !self.reading_done && pool.in_flight() < *cap {
                    match self.frames.next() {
                        Ok(Some(c)) => pool.submit(c),
                        Ok(None) => self.reading_done = true,
                        Err(e) => {
                            pool.submit_result(Err(e));
                            self.reading_done = true;
                        }
                    }
                }
                pool.next().transpose()
            }
        }
    }
}

/// A CRAM decoding failure carried through `io::Error`.
#[derive(Debug)]
struct CramFailure(Box<IoError>);

impl std::fmt::Display for CramFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for CramFailure {}

impl Read for CramStream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let avail = self.fill_buf()?;
        let n = avail.len().min(buf.len());
        buf[..n].copy_from_slice(&avail[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for CramStream {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.at >= self.cur.len() {
            if self.failed {
                return Err(io::Error::other("the CRAM input failed earlier"));
            }
            match self.next_container() {
                Ok(Some(b)) => {
                    self.cur = b;
                    self.at = 0;
                }
                Ok(None) => return Ok(&[]),
                Err(e) => {
                    self.failed = true;
                    return Err(io::Error::other(CramFailure(Box::new(e))));
                }
            }
        }
        Ok(&self.cur[self.at..])
    }

    fn consume(&mut self, n: usize) {
        self.at = (self.at + n).min(self.cur.len());
    }
}

/// The gloci error inside a read error from a CRAM stream, or `e` itself.
pub(crate) fn unwrap_failure(e: IoError) -> IoError {
    match e {
        IoError::Read { path, source }
            if source
                .get_ref()
                .is_some_and(|inner| inner.is::<CramFailure>()) =>
        {
            match source.into_inner().map(|i| i.downcast::<CramFailure>()) {
                Some(Ok(f)) => *f.0,
                _ => IoError::Internal(format!("{path}: a CRAM failure was lost")),
            }
        }
        e => e,
    }
}

// ---- opening ----

/// Check that a seekable CRAM ends with the end-of-file container, before any record is read
/// (as gloci checks a BGZF file's end-of-file block).
fn check_eof_container(name: &str) -> Result<(), IoError> {
    use std::io::{Seek, SeekFrom};
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let mut f = std::fs::File::open(name).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let len = f.metadata().map_err(read_err)?.len();
    let mut tail = Vec::new();
    if len >= CRAM_EOF.len() as u64 {
        f.seek(SeekFrom::Start(len - CRAM_EOF.len() as u64))
            .map_err(read_err)?;
        take_into(&mut f, &mut tail, CRAM_EOF.len()).map_err(read_err)?;
    }
    if tail != CRAM_EOF {
        return Err(truncated(name, "the file does not end with one"));
    }
    Ok(())
}

/// A CRAM input's header, its decoded records as a BAM body, and the counters of its decoders.
struct Opened {
    frames: BamFrames,
    states: Vec<Arc<RefState>>,
}

/// Open the CRAM `input` for the BAM reader's framing: the header (renamed by `map`), and the
/// records decoded on `threads` workers with the reference `reference`.
fn open(
    input: Input,
    map: Option<&ChromMap>,
    reference: Option<&Path>,
    threads: NonZeroUsize,
) -> Result<Opened, IoError> {
    let Input {
        name, mut reader, ..
    } = input;
    let head = read_head(&mut reader, &name)?;
    if is_regular_file(&name) {
        check_eof_container(&name)?;
    }
    let h = parse_header(&name, &head)?;
    let info = Arc::new(CramInfo::new(&name, h.clone(), reference)?);
    let frames = CramFrames::new(&name, reader, &head);
    let stream = CramStream::new(frames, info, threads)?;
    let states = stream.states.clone();
    let mut aln = h;
    if let Some(m) = map {
        m.check_file_names(&name, aln.refs.iter().map(|r| &r.0[..]))?;
        aln = aln.renamed(m);
    }
    Ok(Opened {
        frames: BamFrames::from_parts(name, Compression::None, aln, Box::new(stream)),
        states,
    })
}

fn reference_of(o: &ReaderOptions) -> Option<&Path> {
    o.accept
        .alignments
        .as_ref()
        .and_then(|a| a.reference.as_deref())
}

/// Reads CRAM records for interval commands: the BAM reader over the transcoded records, so
/// spans, skips, filters and warnings are the BAM's (D101, D102).
pub struct CramReader {
    inner: BamReader,
    states: Vec<Arc<RefState>>,
}

impl CramReader {
    pub(crate) fn open_cram(
        input: Input,
        o: &ReaderOptions,
        opts: AlignOptions,
    ) -> Result<CramReader, IoError> {
        let Opened { frames, states } =
            open(input, o.chrom_map.as_deref(), reference_of(o), o.threads)?;
        Ok(CramReader {
            inner: BamReader::from_frames(frames, o.invalid, opts),
            states,
        })
    }

    /// Sequences each decoder read from the FASTA, and the most it held at once.
    pub fn reference_counts(&self) -> Vec<(usize, usize)> {
        counts(&self.states)
    }
}

fn counts(states: &[Arc<RefState>]) -> Vec<(usize, usize)> {
    states
        .iter()
        .map(|s| {
            (
                s.fetches.load(Ordering::Relaxed),
                s.max_live.load(Ordering::Relaxed),
            )
        })
        .collect()
}

impl RecordReader for CramReader {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        Format::Cram
    }
    fn compression(&self) -> Compression {
        Compression::None
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        self.inner
            .read_batch(chroms, batch, max)
            .map_err(unwrap_failure)
    }
}

/// Every record of a CRAM (D120): the BAM source over the transcoded records.
pub(crate) struct CramSource {
    inner: BamSource,
}

impl CramSource {
    pub(crate) fn open(input: Input, o: &ReaderOptions) -> Result<CramSource, IoError> {
        let Opened { frames, .. } =
            open(input, o.chrom_map.as_deref(), reference_of(o), o.threads)?;
        let opts = o.accept.alignments.clone().unwrap_or_default();
        Ok(CramSource {
            inner: BamSource::new(frames, o.invalid, opts),
        })
    }
}

impl AlnSource for CramSource {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn format(&self) -> Format {
        Format::Cram
    }
    fn header(&self) -> &Arc<AlnHeader> {
        self.inner.header()
    }
    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        self.inner.next_record(buf).map_err(unwrap_failure)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
}

fn parse_header(name: &str, head: &CramHead) -> Result<AlnHeader, IoError> {
    AlnHeader::parse_sam(head.text.clone()).map_err(|m| {
        detect_err(
            name,
            format!("{name}: the SAM header of the CRAM file: {m}"),
        )
    })
}

/// The alignment header of the CRAM file `name`, read from `r` (for `-g x.cram`); no reference
/// is needed.
pub fn read_cram_header(r: &mut dyn Read, name: &str) -> Result<AlnHeader, IoError> {
    parse_header(name, &read_head(r, name)?)
}

/// Warn that `--reference` is not used when no input is CRAM (D109).
pub fn warn_unused_reference(reference: Option<&Path>, any_cram: bool) {
    if let Some(r) = reference
        && !any_cram
    {
        log::warn!(
            "--reference {} is not used: no input is CRAM (only CRAM input needs a reference)",
            r.display()
        );
    }
}

/// The header of the CRAM `name` as decoders share it, and its frames.
pub(crate) fn open_for_index(
    name: &str,
    src: Box<dyn BufRead + Send>,
    reference: Option<&Path>,
) -> Result<(Arc<CramInfo>, CramFrames), IoError> {
    let mut src = src;
    let head = read_head(&mut src, name)?;
    if is_regular_file(name) {
        check_eof_container(name)?;
    }
    let h = parse_header(name, &head)?;
    let info = Arc::new(CramInfo::new(name, h, reference)?);
    Ok((info, CramFrames::new(name, src, &head)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diag::InvalidPolicy;
    use crate::input::{InputOptions, open_input};
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../xtask/diff/inputs")
            .join(name)
    }

    fn options(reference: Option<PathBuf>, threads: usize) -> ReaderOptions {
        ReaderOptions {
            threads: NonZeroUsize::new(threads).unwrap(),
            format: None,
            vcf_span: gloci_core::VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
            chrom_map: None,
            accept: crate::align::Accept {
                alignments: Some(AlignOptions {
                    full_records: true,
                    reference,
                    ..AlignOptions::default()
                }),
                ..crate::align::Accept::default()
            },
        }
    }

    fn drain(path: &Path, threads: usize) -> (Vec<Vec<u8>>, Vec<(usize, usize)>) {
        let o = options(Some(fixture("fa/ref.fa")), threads);
        let input = open_input(
            &path.to_string_lossy(),
            &InputOptions {
                threads: o.threads,
                format: None,
            },
        )
        .unwrap();
        let opts = o.accept.alignments.clone().unwrap();
        let mut r = CramReader::open_cram(input, &o, opts).unwrap();
        let Header::Alignment(h) = r.header().clone() else {
            panic!("no alignment header");
        };
        let mut chroms = ChromTable::first_appearance();
        let mut recs = Vec::new();
        let mut batch = RecordBatch::default();
        loop {
            batch.clear();
            let more = r.read_batch(&mut chroms, &mut batch, 3).unwrap();
            for i in 0..batch.len() {
                assert!(crate::align::Aln::parse(batch.raw(i), &h).is_ok());
                recs.push(batch.raw(i).to_vec());
            }
            if !more {
                break;
            }
        }
        (recs, r.reference_counts())
    }

    /// A counting adapter: each decoder holds at most one sequence of a single-reference run, and
    /// reads it once per run (D109).
    #[test]
    fn cram_reference_cache_bounded() {
        // flags.multi.cram: slices of several references (as multi_seq_per_slice=1 writes every
        // slice), five on chr1, then chr1 (an unmapped read) + chr2, then chr2 + unplaced reads:
        // chr1 and chr2 are each read once. tags.cram: chrM, chr1, then chr2 + chr10 (an
        // unmapped read).
        for (file, fetches) in [("aln/flags.multi.cram", 2), ("aln/tags.cram", 3)] {
            let (recs, counts) = drain(&fixture(file), 1);
            assert!(!recs.is_empty());
            assert_eq!(counts.len(), 1);
            let (got, live) = counts[0];
            assert_eq!(got, fetches, "{file}: sequences read");
            assert!(live <= 1, "{file}: {live} sequences held at once");
            // On every worker of a pool: one sequence each, at most.
            for t in [2, 8] {
                let (again, counts) = drain(&fixture(file), t);
                assert_eq!(again, recs, "{file} at -t {t}");
                assert_eq!(counts.len(), t);
                for (_, live) in counts {
                    assert!(live <= 1, "{file} at -t {t}: {live} sequences held at once");
                }
            }
        }
    }
}
