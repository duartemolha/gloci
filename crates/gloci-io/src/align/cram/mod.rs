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
//! - **The reference memory** is shared by every decoder ([`reference::RefCache`]): two slots,
//!   each one reference, for slices whose reads need one (sorted input reads each reference
//!   once and holds at most two whatever the thread count), and shared generations of
//!   reference bases and of placeholders for slices that need several (unsorted or
//!   name-sorted input), reads past a reference's end, or a CRAM stored without the reference.
//!   A slice leases every reference it needs at once, each loaded once; there are no
//!   per-thread copies. A file not sorted by coordinate uses the generations for every slice.
//! - **The pool** (D110, [`crate::pool`]): containers are read in order, decoded on up to
//!   `threads` workers (with 16 MiB stacks for noodles' codecs) with at most `2 × threads` in
//!   flight, and released in order, so the output does not depend on the thread count;
//!   dropping the reader cancels the containers not started.
//! - **Before noodles decodes a slice**, the counts it allocates by are checked against the
//!   slice's bytes, and reads may run at most 1 Mb (or four read lengths) past a reference's
//!   end (`#cram-limits`).
//! - **Transcoding**: each record is encoded as BAM bytes once, with htslib's `bin` and MAPQ 0
//!   for unmapped reads (CRAM does not store their MAPQ), and the bytes go through the BAM
//!   reader's framing ([`BamFrames`]), so interval commands and the every-record source treat a
//!   CRAM exactly as the BAM it encodes. MD and NM are not regenerated (`#cram-md-nm`).
//! - **Flagstat's fields** ([`fields`]): the every-record source gives the flag, reference ids
//!   and MAPQ of each record from the data series they come from alone, and decodes a container
//!   in full where that would not be exact.

mod decode;
mod fields;
mod file;
mod frame;
mod reference;

pub(crate) use decode::ExtentDecoder;
pub(crate) use file::CramFile;
pub(crate) use frame::{CRAM_EOF, CramFrames, RawContainer, slice_heads, unread_version};
pub(crate) use reference::CramInfo;

use super::bam::{BamFrames, BamReader, bam_record_place};
use super::header::AlnHeader;
use super::source::{AlnSource, BamSource, FlagFields, FlagFieldsVisitor, flag_fields_of_records};
use super::{AlignOptions, Aln};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::IoError;
use crate::input::{Input, is_regular_file};
use crate::pool::Ordered;
use crate::record::{Header, ReaderOptions, RecordBatch, RecordReader};
use decode::Decoder;
use fields::{Counted, FieldsCounts, FieldsDecoder};
use frame::{CramHead, detect_err, read_head, take_into, truncated};
use gloci_core::ChromTable;
pub use reference::ReferenceCounts;
use reference::{Counters, RefCache};
use std::io::{self, BufRead, Read};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Arc;

fn decode_job(d: &mut Decoder, c: RawContainer) -> Result<Vec<u8>, IoError> {
    let mut out = Vec::new();
    d.decode(&c, &mut out)?;
    Ok(out)
}

/// Where containers are decoded: on the reading thread (one thread), or on the pool, in order.
enum Mode<S, R: Send + 'static> {
    Inline(Box<S>),
    Pool(Ordered<RawContainer, R>, u64),
}

impl<S: Send + 'static, R: Send + 'static> Mode<S, R> {
    /// One state per thread: inline for one, else a pool running `work`, with at most twice as
    /// many containers in flight as workers.
    fn new(
        mut states: Vec<S>,
        work: fn(&mut S, RawContainer) -> Result<R, IoError>,
    ) -> Result<Self, IoError> {
        if states.len() == 1
            && let Some(one) = states.pop()
        {
            return Ok(Mode::Inline(Box::new(one)));
        }
        let cap = 2 * states.len() as u64;
        Ok(Mode::Pool(Ordered::new(states, work)?, cap))
    }

    /// The next container of `frames` decoded, in file order; `inline` decodes it on this
    /// thread. `reading_done` is set once the pool has been given the last container.
    fn next(
        &mut self,
        frames: &mut CramFrames,
        reading_done: &mut bool,
        inline: impl FnOnce(&mut S, RawContainer) -> Result<R, IoError>,
    ) -> Result<Option<R>, IoError> {
        match self {
            Mode::Inline(d) => match frames.next()? {
                None => Ok(None),
                Some(c) => inline(d, c).map(Some),
            },
            Mode::Pool(pool, cap) => {
                while !*reading_done && pool.in_flight() < *cap {
                    match frames.next() {
                        Ok(Some(c)) => pool.submit(c),
                        Ok(None) => *reading_done = true,
                        Err(e) => {
                            pool.submit_result(Err(e));
                            *reading_done = true;
                        }
                    }
                }
                pool.next().transpose()
            }
        }
    }
}

/// The decoded records of a CRAM as the byte stream of a BAM body: each record with its
/// `block_size`, in file order. A decoding failure is an `io::Error` holding the gloci error
/// ([`unwrap_failure`] takes it out).
struct CramStream {
    frames: CramFrames,
    mode: Mode<Decoder, Vec<u8>>,
    cur: Vec<u8>,
    at: usize,
    reading_done: bool,
    failed: bool,
    /// What the decoders read from the FASTA and held (the bounded-cache tests read it).
    counters: Arc<Counters>,
}

impl CramStream {
    fn new(
        frames: CramFrames,
        info: Arc<CramInfo>,
        threads: NonZeroUsize,
    ) -> Result<Self, IoError> {
        let counters = Arc::new(Counters::default());
        let cache = Arc::new(RefCache::new(info.clone(), counters.clone()));
        let decoders = (0..threads.get())
            .map(|_| Decoder::new(info.clone(), cache.clone()))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(CramStream {
            frames,
            mode: Mode::new(decoders, decode_job)?,
            cur: Vec::new(),
            at: 0,
            reading_done: false,
            failed: false,
            counters,
        })
    }

    /// The next decoded container, in file order.
    fn next_container(&mut self) -> Result<Option<Vec<u8>>, IoError> {
        let cur = &mut self.cur;
        self.mode
            .next(&mut self.frames, &mut self.reading_done, |d, c| {
                let mut out = std::mem::take(cur);
                out.clear();
                d.decode(&c, &mut out)?;
                Ok(out)
            })
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

/// A CRAM input before any record is read: its name, its header (renamed by `--chrom-map`), the
/// header as its decoders share it, and its data containers.
struct Parts {
    name: String,
    aln: AlnHeader,
    info: Arc<CramInfo>,
    frames: CramFrames,
}

/// Open the CRAM `input`: its header (renamed by `map`), checked against the reference
/// `reference`, and its frames.
fn open_parts(
    input: Input,
    map: Option<&ChromMap>,
    reference: Option<&Path>,
) -> Result<Parts, IoError> {
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
    let mut aln = h;
    if let Some(m) = map {
        m.check_file_names(&name, aln.refs.iter().map(|r| &r.0[..]))?;
        aln = aln.renamed(m);
    }
    Ok(Parts {
        name,
        aln,
        info,
        frames,
    })
}

/// A CRAM input's header, its decoded records as a BAM body, and the counters of its decoders.
struct Opened {
    frames: BamFrames,
    counters: Arc<Counters>,
}

/// The BAM reader's framing over the records of `frames`, decoded on `threads` workers.
fn bam_frames(
    name: String,
    aln: impl Into<Arc<AlnHeader>>,
    info: Arc<CramInfo>,
    frames: CramFrames,
    threads: NonZeroUsize,
) -> Result<Opened, IoError> {
    let stream = CramStream::new(frames, info, threads)?;
    let counters = stream.counters.clone();
    Ok(Opened {
        frames: BamFrames::from_parts(name, Compression::None, aln, Box::new(stream))
            .without_file_numbers(),
        counters,
    })
}

/// Open the CRAM `input` for the BAM reader's framing: the header (renamed by `map`), and the
/// records decoded on `threads` workers with the reference `reference`.
fn open(
    input: Input,
    map: Option<&ChromMap>,
    reference: Option<&Path>,
    threads: NonZeroUsize,
) -> Result<Opened, IoError> {
    let Parts {
        name,
        aln,
        info,
        frames,
    } = open_parts(input, map, reference)?;
    bam_frames(name, aln, info, frames, threads)
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
    counters: Arc<Counters>,
}

impl CramReader {
    pub(crate) fn open_cram(
        input: Input,
        o: &ReaderOptions,
        opts: AlignOptions,
    ) -> Result<CramReader, IoError> {
        let Opened { frames, counters } =
            open(input, o.chrom_map.as_deref(), reference_of(o), o.threads)?;
        Ok(CramReader {
            inner: BamReader::from_frames(frames, o.invalid, opts),
            counters,
        })
    }

    /// What the decoders read from the FASTA so far, and the most they held at once.
    pub fn reference_counts(&self) -> ReferenceCounts {
        self.counters.snapshot()
    }
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
    fn whole_records(&self) -> bool {
        self.inner.whole_records()
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

/// Every record of a CRAM (D120): the BAM source over the transcoded records, started at the
/// first record read. Flagstat's fields come from the data series they need alone
/// ([`fields`]), with the record path's filters, errors and places.
pub(crate) struct CramSource {
    name: String,
    header: Arc<AlnHeader>,
    policy: InvalidPolicy,
    opts: AlignOptions,
    threads: NonZeroUsize,
    state: SourceState,
    /// What the fields path dropped or skipped.
    diag: Diagnostics,
    /// The reference id of the fields path's current run of records, and how many of the run
    /// it has read: its errors name "record N of CHROM" as the record path does.
    run: (i32, u64),
    /// Containers the fields path read alone, and those it decoded in full.
    counts: Arc<FieldsCounts>,
}

enum SourceState {
    /// Nothing read yet.
    Unread {
        info: Arc<CramInfo>,
        frames: CramFrames,
    },
    /// Records are read through the BAM source over the transcoded records.
    Records(Box<BamSource>),
    /// The fields path has read the file.
    Done,
}

impl CramSource {
    pub(crate) fn open(input: Input, o: &ReaderOptions) -> Result<CramSource, IoError> {
        let Parts {
            name,
            aln,
            info,
            frames,
        } = open_parts(input, o.chrom_map.as_deref(), reference_of(o))?;
        Ok(CramSource {
            name,
            header: Arc::new(aln),
            policy: o.invalid,
            opts: o.accept.alignments.clone().unwrap_or_default(),
            threads: o.threads,
            state: SourceState::Unread { info, frames },
            diag: Diagnostics::default(),
            run: (i32::MIN, 0),
            counts: Arc::default(),
        })
    }

    /// The containers the fields path read alone and those it decoded in full, so far.
    #[cfg(test)]
    pub(crate) fn fields_counts(&self) -> (u64, u64) {
        self.counts.get()
    }

    /// The record path's source, started on first use; `None` once the fields path has read
    /// the file.
    fn records(&mut self) -> Result<Option<&mut BamSource>, IoError> {
        if matches!(self.state, SourceState::Unread { .. }) {
            let SourceState::Unread { info, frames } =
                std::mem::replace(&mut self.state, SourceState::Done)
            else {
                unreachable!("the state was just matched");
            };
            let Opened { frames, .. } = bam_frames(
                self.name.clone(),
                self.header.clone(),
                info,
                frames,
                self.threads,
            )?;
            let records = BamSource::new(frames, self.policy, self.opts.clone());
            self.state = SourceState::Records(Box::new(records));
        }
        Ok(match &mut self.state {
            SourceState::Records(s) => Some(&mut **s),
            _ => None,
        })
    }

    /// Count a record of reference `ref_id` in the run of records the fields path reads.
    fn see(&mut self, ref_id: i32) {
        if ref_id != self.run.0 {
            self.run = (ref_id, 0);
        }
        self.run.1 += 1;
    }

    /// Give `f` the record `fields`, unless the read filters drop it.
    fn deliver(
        &mut self,
        fields: FlagFields,
        f: &mut FlagFieldsVisitor<'_>,
    ) -> Result<(), IoError> {
        if self.opts.drops(fields.flag, fields.mapq) {
            self.diag.filtered += 1;
            return Ok(());
        }
        f(fields)
    }

    /// Give `f` the fields of the BAM records `bytes` (each with its `block_size`) of a container
    /// decoded in full, checked, filtered and reported as the record path's BAM source does.
    fn walk_records(
        &mut self,
        mut bytes: &[u8],
        f: &mut FlagFieldsVisitor<'_>,
    ) -> Result<(), IoError> {
        while let Some((head, rest)) = bytes.split_first_chunk::<4>() {
            let size = u32::from_le_bytes(*head) as usize;
            let Some(body) = rest.get(..size) else {
                return Err(IoError::Internal(format!(
                    "{}: a decoded CRAM record runs past its container",
                    self.name
                )));
            };
            bytes = &rest[size..];
            if let Some(r) = body.first_chunk::<4>() {
                self.see(i32::from_le_bytes(*r));
            }
            // As BamSource::next_record: the same check, filters, message and place.
            let problem = match Aln::parse_framed(body, &self.header) {
                Ok(_) => {
                    if let Some(fields) = FlagFields::of_bam(body) {
                        self.deliver(fields, f)?;
                    }
                    continue;
                }
                Err(msg) => msg,
            };
            let place = bam_record_place(&self.header.refs, self.run.0, self.run.1, None);
            match self.policy {
                InvalidPolicy::Fail => {
                    return Err(IoError::Detect {
                        path: self.name.clone(),
                        msg: format!("{} {place}: {problem}", self.name),
                    });
                }
                InvalidPolicy::Skip => self.diag.record_place(&place, &problem),
            }
        }
        Ok(())
    }
}

fn fields_job(d: &mut FieldsDecoder, c: RawContainer) -> Result<Counted, IoError> {
    d.decode(c)
}

impl AlnSource for CramSource {
    fn name(&self) -> &str {
        &self.name
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn format(&self) -> Format {
        Format::Cram
    }
    fn header(&self) -> &Arc<AlnHeader> {
        &self.header
    }
    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        match self.records()? {
            Some(s) => s.next_record(buf).map_err(unwrap_failure),
            None => Ok(None),
        }
    }
    fn for_each_flag_fields(&mut self, f: &mut FlagFieldsVisitor<'_>) -> Result<(), IoError> {
        let (info, mut frames) = match std::mem::replace(&mut self.state, SourceState::Done) {
            SourceState::Unread { info, frames } => (info, frames),
            other => {
                // Records were read already: the rest come the same way.
                self.state = other;
                return flag_fields_of_records(self, f);
            }
        };
        // The run that reads the file owns its decoders (D91); the pool's workers take the
        // run's token, and this loop polls it at each container.
        let watch = crate::cancel::Watch::current();
        let cache = Arc::new(RefCache::new(info.clone(), Arc::new(Counters::default())));
        let decoders = (0..self.threads.get())
            .map(|_| {
                let full = Decoder::new(info.clone(), cache.clone())?;
                Ok(FieldsDecoder::new(full, self.counts.clone()))
            })
            .collect::<Result<Vec<_>, IoError>>()?;
        let mut mode = Mode::new(decoders, fields_job)?;
        let mut reading_done = false;
        while let Some(counted) = mode.next(&mut frames, &mut reading_done, |d, c| d.decode(c))? {
            watch.check()?;
            match counted {
                Counted::Fields(records) => {
                    for r in records {
                        self.see(r.ref_id);
                        self.deliver(r.into(), f)?;
                    }
                }
                Counted::Records(bytes) => self.walk_records(&bytes, f)?,
            }
        }
        Ok(())
    }
    fn diagnostics(&self) -> &Diagnostics {
        match &self.state {
            SourceState::Records(s) => s.diagnostics(),
            _ => &self.diag,
        }
    }
    fn place(&self) -> String {
        match &self.state {
            SourceState::Records(s) => s.place(),
            _ => bam_record_place(&self.header.refs, self.run.0, self.run.1, None),
        }
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

    fn drain(path: &Path, threads: usize) -> (Vec<Vec<u8>>, ReferenceCounts) {
        let noref = path.to_string_lossy().contains("noref");
        let o = options((!noref).then(|| fixture("fa/ref.fa")), threads);
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

    /// Without the reference, a slice of several references gets a placeholder per reference,
    /// as long as the reference, shared by every decoder: the counts are totals across the
    /// decoders (D109). contigs.noref.cram (chr1 of 10 Mb with reads near its end, then 200
    /// contigs of 50 kb) holds at most one placeholder of each, 20 MB in all;
    /// pairs.namesorted.noref.cram (two references of 10 Mb, every slice needing both) holds
    /// exactly one of each whatever the thread count.
    #[test]
    fn cram_placeholders_bounded() {
        // mixed.namesorted.noref.cram: chr1 of 10 Mb and 300 contigs of 5 kb, every slice needing
        // chr1 and other contigs. The small ones stay (they are few bytes); chr1's placeholder is
        // made once: one generation, whatever the thread count.
        let mixed = drain(&fixture("aln/mixed.namesorted.noref.cram"), 1);
        assert_eq!(mixed.0.len(), 600);
        for t in [1, 2, 8] {
            let (recs, counts) = drain(&fixture("aln/mixed.namesorted.noref.cram"), t);
            assert_eq!(recs, mixed.0, "-t {t}");
            assert_eq!(counts.generations, 1, "mixed -t {t}: {counts:?}");
            assert!(
                counts.max_placeholder_bytes <= 10_000_000 + 300 * 5_000,
                "mixed -t {t}: {counts:?}"
            );
        }
        let one = drain(&fixture("aln/contigs.noref.cram"), 1);
        let pairs = drain(&fixture("aln/pairs.namesorted.noref.cram"), 1);
        assert_eq!(pairs.0.len(), 300);
        for t in [1, 2, 8] {
            let (recs, counts) = drain(&fixture("aln/contigs.noref.cram"), t);
            assert_eq!(recs, one.0, "-t {t}");
            assert_eq!(recs.len(), 430);
            assert_eq!(counts.fetches, 0);
            assert!(
                (10_000_000..=20_000_000).contains(&counts.max_placeholder_bytes),
                "contigs -t {t}: {counts:?}"
            );
            let (recs, counts) = drain(&fixture("aln/pairs.namesorted.noref.cram"), t);
            assert_eq!(recs, pairs.0, "-t {t}");
            assert_eq!(
                counts.max_placeholder_bytes, 20_000_000,
                "pairs -t {t}: {counts:?}"
            );
            assert_eq!(counts.max_live, 2, "pairs -t {t}: {counts:?}");
        }
    }

    /// A name-tokeniser block whose clear-text header claims a decoded length other than its
    /// block's, or more names than the slice has records, is refused before noodles reserves by
    /// them (the reviewer's 4 GB and 10.7 M-name reproducers, rebuilt from cat-big.v31.cram with
    /// the block checksum recomputed).
    #[test]
    fn tokeniser_header_is_checked() {
        let path = fixture("aln/cat-big.v31.cram");
        let bytes = std::fs::read(&path).unwrap();
        let mut f = std::io::BufReader::new(std::fs::File::open(&path).unwrap());
        let head = read_head(&mut f, "x").unwrap();
        let c = CramFrames::new("x", Box::new(f), &head)
            .next()
            .unwrap()
            .unwrap();
        // (start of the block, start of its data, end of its data) of the tokeniser block.
        let base = c.offset as usize + c.header.header_len;
        let mut rest = c.data();
        let mut found = None;
        while !rest.is_empty() {
            let start = base + (c.data().len() - rest.len());
            let b = frame::read_block(&mut rest).unwrap();
            let data = b.data.as_ptr() as usize - c.bytes.as_ptr() as usize + c.offset as usize;
            if b.method == 8 {
                found = Some((start, data, data + b.data.len()));
            }
        }
        let (start, data, end) = found.expect("a tokeniser block");
        for (at, value, want) in [
            (
                data,
                4_009_763_914u32,
                "says it decodes to 4,009,763,914 bytes",
            ),
            (
                data + 4,
                10_683_668,
                "claims 10,683,668 names in 9,290 bytes",
            ),
        ] {
            let mut b = bytes.clone();
            b[at..at + 4].copy_from_slice(&value.to_le_bytes());
            let mut crc = flate2::Crc::new();
            crc.update(&b[start..end]);
            b[end..end + 4].copy_from_slice(&crc.sum().to_le_bytes());
            let o = options(None, 2);
            let e =
                crate::align::open_aln_source_from("x.cram", Box::new(std::io::Cursor::new(b)), &o)
                    .and_then(|mut s| {
                        let mut buf = Vec::new();
                        while s.next_record(&mut buf)?.is_some() {
                            buf.clear();
                        }
                        Ok(())
                    })
                    .unwrap_err()
                    .to_string();
            assert!(e.contains(want) && e.contains("x.cram"), "{e}");
        }
    }

    /// A counting adapter: the decoders share the references, read each from the FASTA once and
    /// hold at most two at a time in all, whatever the thread count (D109).
    #[test]
    fn cram_reference_cache_bounded() {
        // flags.multi.cram: slices of several references (as multi_seq_per_slice=1 writes every
        // slice), five on chr1, then chr1 (an unmapped read) + chr2, then chr2 + unplaced reads;
        // their reads' references are known from their RI and BF series, so each slice that
        // needs one reference takes it from the shared cache. tags.cram: chrM, chr1, then
        // chr2 + chr10 (an unmapped read).
        // flags.namesorted.cram: name-sorted, so every slice needs chr1 and chr2; the decoders
        // share one copy of each, read once.
        for (file, fetches) in [
            ("aln/flags.multi.cram", 2),
            ("aln/tags.cram", 3),
            ("aln/flags.namesorted.cram", 2),
        ] {
            let (recs, counts) = drain(&fixture(file), 1);
            assert!(!recs.is_empty());
            assert_eq!(counts.fetches, fetches, "{file}: sequences read");
            assert!(counts.max_live <= 2, "{file}: {counts:?}");
            for t in [2, 8] {
                let (again, counts) = drain(&fixture(file), t);
                assert_eq!(again, recs, "{file} at -t {t}");
                assert_eq!(counts.fetches, fetches, "{file} at -t {t}: sequences read");
                assert!(counts.max_live <= 2, "{file} at -t {t}: {counts:?}");
            }
        }
    }
}
