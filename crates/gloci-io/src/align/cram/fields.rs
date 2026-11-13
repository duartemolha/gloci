//! The fields `flagstat` counts, decoded from a CRAM without the rest of each record (Extra X1):
//! a read's BAM flag, reference id, mate reference id and MAPQ come from the data series BF, CF,
//! RI (in a slice of several references), MQ (mapped reads), MF and NS (detached mates) and NF
//! (mates further on in the slice), as `samtools flagstat` reads them through htslib.
//!
//! - **Exactly the full decode.** Each value is read as noodles-cram reads it, and mates within a
//!   slice are resolved as its `resolve_mates` does (a read gets its mate's reference id, and
//!   0x20 and 0x8 from its mate's 0x10 and 0x4), so the fields equal those of the records the full
//!   decode transcodes: MAPQ 0 for an unmapped read, 255 kept for a mapped one.
//! - **Only those blocks.** A series must be in an external block no other encoding reads, or
//!   coded as a constant (a one-symbol Huffman table, which reads no bits): anything else would
//!   interleave it with series flagstat does not read. noodles-cram 0.102 keeps its block codecs
//!   to itself, so the container goes to it reframed: its compression header and slice headers
//!   as stored, the blocks of those series as stored, and every other block empty. noodles
//!   parses the headers as it does for the full decode and decompresses those blocks alone.
//! - **Fallback.** A container the fields cannot be read from exactly (an encoding as above, a
//!   value out of range, a count that disagrees with the slice's, a block that does not frame,
//!   check or decompress) is decoded in full instead, with a debug line: what the full decode
//!   gives, records or an error, is the answer, so errors keep their wording and place.
//! - **Checked as before:** every block's CRC32 and stated sizes, the slice headers, `--reference`
//!   for a file that needs one (D109), and the slice counts on the series read. **Not checked:**
//!   the series flagstat does not read (positions, names, bases, qualities, tags), so a read past
//!   gloci's overhang limit, a reference whose bases differ from the file's checksums, or damage
//!   confined to those series in blocks whose checksums still match (a crafted file) goes
//!   unreported, as in `samtools flagstat` (#cram-limits).

use super::decode::{Bases, Decoder, SLACK_RECORDS, decoder_note};
use super::frame::{
    Block, CompressionInfo, RawContainer, SliceHead, check_block, check_compression_block,
    compression, crc32, itf8, read_block, slice_heads,
};
use crate::align::FlagFields;
use crate::align::reheader::itf8 as put_itf8;
use crate::error::IoError;
use noodles_cram as cram;
use std::borrow::Cow;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

/// Block content types (CRAM 3 §8).
const CORE: u8 = 5;
const EXTERNAL: u8 = 4;
/// The highest block compression method noodles reads (the name tokeniser).
const LAST_METHOD: u8 = 8;
/// CRAM flags (CF): the mate is stored with the read, or is a later read of the slice.
const DETACHED: u8 = 0x2;
const DOWNSTREAM: u8 = 0x4;
/// Mate flags (MF): the mate is on the reverse strand, or unmapped.
const MATE_REVERSE: u8 = 0x1;
const MATE_UNMAPPED: u8 = 0x2;
/// BAM flags.
const REVERSE: u16 = 0x10;
const UNMAPPED: u16 = 0x4;
const MREVERSE: u16 = 0x20;
const MUNMAP: u16 = 0x8;

/// The fields of one record `flagstat` counts, as the full decode transcodes them: the BAM flag,
/// MAPQ, and reference ids (-1 for none).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Fields {
    pub(crate) flag: u16,
    pub(crate) mapq: u8,
    pub(crate) ref_id: i32,
    pub(crate) mate_ref_id: i32,
}

impl From<Fields> for FlagFields {
    fn from(r: Fields) -> FlagFields {
        FlagFields {
            flag: r.flag,
            ref_id: usize::try_from(r.ref_id).ok(),
            mate_ref_id: usize::try_from(r.mate_ref_id).ok(),
            mapq: r.mapq,
        }
    }
}

/// A container decoded for `flagstat`.
pub(crate) enum Counted {
    /// The fields of its records, read alone.
    Fields(Vec<Fields>),
    /// Its records decoded in full: BAM records, each with its `block_size`.
    Records(Vec<u8>),
}

/// How many containers were read for their fields alone, and how many were decoded in full.
#[derive(Debug, Default)]
pub(crate) struct FieldsCounts {
    fast: AtomicU64,
    full: AtomicU64,
}

impl FieldsCounts {
    #[cfg(test)]
    pub(crate) fn get(&self) -> (u64, u64) {
        (
            self.fast.load(Ordering::Relaxed),
            self.full.load(Ordering::Relaxed),
        )
    }
}

/// Decodes containers for `flagstat`: the fields alone where that is exact, else in full. Each
/// worker of the pool has one.
pub(crate) struct FieldsDecoder {
    full: Decoder,
    counts: Arc<FieldsCounts>,
}

impl FieldsDecoder {
    pub(crate) fn new(full: Decoder, counts: Arc<FieldsCounts>) -> FieldsDecoder {
        FieldsDecoder { full, counts }
    }

    pub(crate) fn decode(&mut self, c: RawContainer) -> Result<Counted, IoError> {
        match self.fields(&c) {
            Ok(records) => {
                self.counts.fast.fetch_add(1, Ordering::Relaxed);
                Ok(Counted::Fields(records))
            }
            Err(why) => {
                log::debug!(
                    "{}: {} is decoded in full for flagstat: {why}",
                    self.full.info.name,
                    c.place()
                );
                self.counts.full.fetch_add(1, Ordering::Relaxed);
                let mut out = Vec::new();
                self.full.decode(&c, &mut out)?;
                Ok(Counted::Records(out))
            }
        }
    }

    /// The fields of the records of `c`, in order, or why they cannot be read alone exactly.
    fn fields(&self, c: &RawContainer) -> Result<Vec<Fields>, String> {
        let info = &self.full.info;
        // noodles' codecs trust what is inside a compressed stream: a panic there names the
        // container, as in the full decode (#cram-limits).
        let _note = decoder_note(&info.name, c);
        let n_refs = info.aln.refs.len();
        let heads = slice_heads(c, &info.name, n_refs).map_err(|e| e.to_string())?;
        let comp = compression(c);
        if let Some(what) = &comp.unsupported {
            return Err(format!("it uses the encoding {what}"));
        }
        check_compression_block(c)?;
        // Where each slice's reads would find their bases, checked as the full decode checks it
        // (a file that needs the reference needs --reference), though no base is read.
        let bases = heads
            .iter()
            .map(|head| self.full.bases(c, head, comp.required))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| e.to_string())?;
        let plan = Plan::new(&comp);
        let framed = reframe(c, &heads, &plan, comp.external.get(b"RN").copied())?;
        for ((head, base), ids) in heads.iter().zip(&bases).zip(&framed.ids) {
            if *base == Bases::Embedded && !ids.contains(&head.embedded) {
                return Err("a slice names an embedded reference block it does not have".into());
            }
        }
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&framed.bytes[..])
            .read_container(&mut container)
            .map_err(|e| e.to_string())?;
        container.compression_header().map_err(|e| e.to_string())?;
        let mut out = Vec::new();
        for (head, slice) in heads.iter().zip(container.slices()) {
            let slice = slice.map_err(|e| e.to_string())?;
            let (_, ext) = slice.decode_blocks().map_err(|e| e.to_string())?;
            read_slice(head, n_refs, &plan, &ext, &mut out)?;
        }
        Ok(out)
    }
}

/// Where a data series `flagstat` reads comes from in a container.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Source {
    /// Every record has this value (a one-symbol Huffman table).
    Constant(i32),
    /// The external block with this content id, which no other encoding reads.
    Block(i32),
    /// Anywhere else: reading it would mean decoding other series.
    Elsewhere,
}

/// The series `flagstat` reads, and where each comes from.
struct Plan {
    bf: Source,
    cf: Source,
    ri: Source,
    mq: Source,
    mf: Source,
    ns: Source,
    nf: Source,
}

impl Plan {
    fn new(comp: &CompressionInfo) -> Plan {
        let of = |key: &[u8; 2]| match (comp.constant.get(key), comp.own_block(key)) {
            (Some(&v), None) => Source::Constant(v),
            (None, Some(id)) => Source::Block(id),
            _ => Source::Elsewhere,
        };
        Plan {
            bf: of(b"BF"),
            cf: of(b"CF"),
            ri: of(b"RI"),
            mq: of(b"MQ"),
            mf: of(b"MF"),
            ns: of(b"NS"),
            nf: of(b"NF"),
        }
    }

    /// Whether slice `head` reads the external block `id`.
    fn reads(&self, head: &SliceHead, id: i32) -> bool {
        let block = Source::Block(id);
        [self.bf, self.cf, self.mq, self.mf, self.ns, self.nf].contains(&block)
            || head.ref_id == -2 && self.ri == block
    }
}

/// A container framed for noodles to decode the blocks `flagstat` reads ([`reframe`]), and the
/// content ids of each slice's external blocks, in order.
struct Framed {
    bytes: Vec<u8>,
    ids: Vec<Vec<i32>>,
}

/// The container `c` reframed: its compression header and slice headers as stored, in each slice
/// the external blocks `plan` reads as stored, and every other block (the core block included)
/// empty, with its content type and id. Every block is CRC-checked and size-checked as the full
/// decode checks it (`rn`: the read names' block), and each slice must hold the blocks its header
/// lists, and no more, and no more records than the full decode accepts.
fn reframe(
    c: &RawContainer,
    heads: &[SliceHead],
    plan: &Plan,
    rn: Option<i32>,
) -> Result<Framed, String> {
    let data = c.data();
    let first = c.header.landmarks.first().copied().unwrap_or(data.len());
    let mut body = data[..first].to_vec();
    let mut landmarks = Vec::with_capacity(heads.len());
    let mut ids = Vec::with_capacity(heads.len());
    for head in heads {
        landmarks.push(body.len());
        let mut src = &data[head.landmark..head.landmark + head.size];
        let all = src;
        read_block(&mut src)?;
        body.extend_from_slice(&all[..all.len() - src.len()]);
        let check = |b: &Block<'_>, want: u8| -> Result<usize, String> {
            if b.content_type != want || b.method > LAST_METHOD {
                return Err(format!(
                    "a slice block has content type {} and compression method {}",
                    b.content_type, b.method
                ));
            }
            check_block(b, (rn == Some(b.content_id)).then_some(head.records))?;
            // The bytes noodles decodes it to, for the records-against-bytes limit: a raw block,
            // or one of raw size 0 whatever its method (noodles reads that as raw), is its stored
            // bytes; any other is its raw size, which check_block holds the codec's own header to
            // (gzip, bzip2 and lzma decode into a buffer of that size). The full decode counts
            // what noodles returned; the two agree but where a codec returns other than its
            // stated size, which is an error in both decodes.
            Ok(if b.method == 0 || b.raw_size == 0 {
                b.data.len()
            } else {
                b.raw_size
            })
        };
        let core = read_block(&mut src)?;
        let mut bytes = check(&core, CORE)?;
        empty_block(&mut body, CORE, core.content_id);
        // No capacity from the header's block count: a corrupt one must not ask for gigabytes.
        let mut slice_ids = Vec::new();
        for _ in 1..head.blocks {
            let at = src;
            let b = read_block(&mut src)?;
            bytes += check(&b, EXTERNAL)?;
            if plan.reads(head, b.content_id) {
                body.extend_from_slice(&at[..at.len() - src.len()]);
            } else {
                empty_block(&mut body, EXTERNAL, b.content_id);
            }
            slice_ids.push(b.content_id);
        }
        if !src.is_empty() {
            return Err("a slice holds more blocks than its header lists".into());
        }
        if head.records > SLACK_RECORDS + bytes {
            return Err(format!(
                "a slice claims {} records in {bytes} bytes",
                head.records
            ));
        }
        ids.push(slice_ids);
    }
    let mut bytes = Vec::with_capacity(body.len() + 64);
    container_header(&mut bytes, body.len(), &landmarks)?;
    bytes.extend_from_slice(&body);
    Ok(Framed { bytes, ids })
}

/// Append a raw block of no bytes with this content type and id.
fn empty_block(out: &mut Vec<u8>, content_type: u8, id: i32) {
    let start = out.len();
    out.extend_from_slice(&[0, content_type]);
    put_itf8(id as u32, out);
    out.extend_from_slice(&[0, 0]);
    let crc = crc32(&out[start..]);
    out.extend_from_slice(&crc.to_le_bytes());
}

/// Append the header of a container of `len` bytes whose slices start at `landmarks`, with no
/// reference context and no counts (noodles reads neither to decode blocks).
fn container_header(out: &mut Vec<u8>, len: usize, landmarks: &[usize]) -> Result<(), String> {
    let too_big = || "a container is too big to reframe".to_string();
    let start = out.len();
    let len = i32::try_from(len).map_err(|_| too_big())?;
    out.extend_from_slice(&len.to_le_bytes());
    put_itf8(u32::MAX, out); // reference id -1
    out.extend_from_slice(&[0, 0, 0, 0, 0, 0]); // start, span, records, counter, bases, blocks
    put_itf8(landmarks.len() as u32, out);
    for &l in landmarks {
        put_itf8(u32::try_from(l).map_err(|_| too_big())?, out);
    }
    let crc = crc32(&out[start..]);
    out.extend_from_slice(&crc.to_le_bytes());
    Ok(())
}

/// The values of one data series in a slice, read in record order.
enum Series<'a> {
    Constant(i32),
    Values(&'a [u8]),
    Elsewhere(&'static str),
}

impl<'a> Series<'a> {
    /// Series `key`, from `source`, among the decoded external blocks `ext` of a slice. A block
    /// the slice lacks holds no values: reading one fails, as noodles' read does.
    fn of(
        source: Source,
        ext: &'a [(i32, Cow<'_, [u8]>)],
        key: &'static str,
    ) -> Result<Series<'a>, String> {
        Ok(match source {
            Source::Constant(v) => Series::Constant(v),
            Source::Elsewhere => Series::Elsewhere(key),
            Source::Block(id) => {
                let mut found = ext.iter().filter(|(i, _)| *i == id);
                match (found.next(), found.next()) {
                    (Some((_, b)), None) => Series::Values(b),
                    (None, _) => Series::Values(&[]),
                    (Some(_), Some(_)) => {
                        return Err(format!("a slice has two blocks with content id {id}"));
                    }
                }
            }
        })
    }

    fn next(&mut self) -> Result<i32, String> {
        match self {
            Series::Constant(v) => Ok(*v),
            Series::Values(s) => {
                itf8(s).ok_or_else(|| "a data series holds fewer values than its reads".into())
            }
            Series::Elsewhere(key) => Err(format!(
                "its {key} data series is not in an external block of its own"
            )),
        }
    }

    /// Whether every value was read (a constant has as many as are asked for).
    fn exhausted(&self) -> bool {
        !matches!(self, Series::Values(s) if !s.is_empty())
    }
}

/// Append the fields of the records of slice `head` to `out`, from its decoded external blocks
/// `ext`, with the mates within the slice resolved.
fn read_slice(
    head: &SliceHead,
    n_refs: usize,
    plan: &Plan,
    ext: &[(i32, Cow<'_, [u8]>)],
    out: &mut Vec<Fields>,
) -> Result<(), String> {
    let multi = head.ref_id == -2;
    let mut bf = Series::of(plan.bf, ext, "BF")?;
    let mut cf = Series::of(plan.cf, ext, "CF")?;
    let mut ri = Series::of(
        if multi { plan.ri } else { Source::Constant(-1) },
        ext,
        "RI",
    )?;
    let mut mq = Series::of(plan.mq, ext, "MQ")?;
    let mut mf = Series::of(plan.mf, ext, "MF")?;
    let mut ns = Series::of(plan.ns, ext, "NS")?;
    let mut nf = Series::of(plan.nf, ext, "NF")?;
    let reference = |v: i32| {
        if v == -1 || usize::try_from(v).is_ok_and(|r| r < n_refs) {
            Ok(v)
        } else {
            Err(format!("a read names reference id {v}"))
        }
    };
    let n = head.records;
    let first = out.len();
    out.reserve(n);
    // The later read of the slice each read names as its mate (NF).
    let mut mates: Vec<Option<usize>> = Vec::new();
    for i in 0..n {
        let flag = bf.next()?;
        let flag = u16::try_from(flag).map_err(|_| format!("a read has BAM flag {flag}"))?;
        let cram_flags = cf.next()?;
        let cram_flags =
            u8::try_from(cram_flags).map_err(|_| format!("a read has CRAM flags {cram_flags}"))?;
        let ref_id = if multi {
            reference(ri.next()?)?
        } else {
            head.ref_id
        };
        let mut r = Fields {
            flag,
            mapq: 0,
            ref_id,
            mate_ref_id: -1,
        };
        if cram_flags & DETACHED != 0 {
            let m = mf.next()?;
            let m = u8::try_from(m).map_err(|_| format!("a read has mate flags {m}"))?;
            if m & MATE_REVERSE != 0 {
                r.flag |= MREVERSE;
            }
            if m & MATE_UNMAPPED != 0 {
                r.flag |= MUNMAP;
            }
            r.mate_ref_id = reference(ns.next()?)?;
        } else if cram_flags & DOWNSTREAM != 0 {
            let d = nf.next()?;
            let mate = usize::try_from(d)
                .ok()
                .map(|d| i + d + 1)
                .filter(|&j| j < n)
                .ok_or_else(|| format!("read {} names its mate {d} reads further on", i + 1))?;
            if mates.is_empty() {
                mates.resize(n, None);
            }
            mates[i] = Some(mate);
        }
        if r.flag & UNMAPPED == 0 {
            if ref_id < 0 {
                return Err("a mapped read has no reference".into());
            }
            let q = mq.next()?;
            r.mapq = u8::try_from(q).map_err(|_| format!("a read has MAPQ {q}"))?;
        }
        out.push(r);
    }
    // BF, CF and (in a slice of several references) RI hold one value per read.
    if !bf.exhausted() || !cf.exhausted() || multi && !ri.exhausted() {
        return Err("a data series holds more values than the slice has reads".into());
    }
    if !mates.is_empty() {
        resolve_mates(&mut out[first..], &mut mates);
    }
    Ok(())
}

/// Give each read whose mate is a later read of the slice its mate's reference id, and 0x20 and
/// 0x8 from its mate's 0x10 and 0x4, as noodles-cram 0.102's `resolve_mates` does: along each
/// chain of reads, each takes the next, and the last takes the first.
fn resolve_mates(recs: &mut [Fields], mates: &mut [Option<usize>]) {
    let set_mate = |recs: &mut [Fields], to: usize, from: usize| {
        let (flag, ref_id) = (recs[from].flag, recs[from].ref_id);
        if flag & REVERSE != 0 {
            recs[to].flag |= MREVERSE;
        }
        if flag & UNMAPPED != 0 {
            recs[to].flag |= MUNMAP;
        }
        recs[to].mate_ref_id = ref_id;
    };
    for i in 0..recs.len() {
        if mates[i].is_none() {
            continue;
        }
        let mut j = i;
        while let Some(m) = mates[j] {
            set_mate(recs, j, m);
            j = m;
        }
        set_mate(recs, j, i);
        let mut j = i;
        while let Some(m) = mates[j].take() {
            j = m;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::CramSource;
    use super::super::frame::{CramFrames, ltf8, read_head};
    use super::*;
    use crate::align::{AlignOptions, AlnSource};
    use crate::cancel::{RunGuard, RunToken};
    use crate::diag::InvalidPolicy;
    use crate::input::{InputOptions, open_input_from};
    use crate::record::ReaderOptions;
    use std::collections::HashSet;
    use std::io::{Cursor, Read};
    use std::num::NonZeroUsize;
    use std::path::{Path, PathBuf};

    fn root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    fn aln(name: &str) -> PathBuf {
        root().join("xtask/diff/inputs/aln").join(name)
    }

    fn fasta(name: &str) -> PathBuf {
        root().join("xtask/diff/inputs/fa").join(name)
    }

    /// flagstat's reader options: no read filters.
    fn options(reference: Option<&Path>, threads: usize) -> ReaderOptions {
        ReaderOptions {
            threads: NonZeroUsize::new(threads).unwrap(),
            format: None,
            vcf_span: gloci_core::VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
            chrom_map: None,
            accept: crate::align::Accept {
                alignments: Some(AlignOptions {
                    reference: reference.map(Path::to_path_buf),
                    ..AlignOptions::default()
                }),
                ..crate::align::Accept::default()
            },
        }
    }

    /// The CRAM `bytes`, named `name`, as an every-record source.
    fn source(
        name: &str,
        src: Box<dyn Read + Send>,
        reference: Option<&Path>,
        threads: usize,
    ) -> Result<CramSource, IoError> {
        let o = options(reference, threads);
        let opts = InputOptions {
            threads: o.threads,
            format: None,
        };
        CramSource::open(open_input_from(name, src, &opts)?, &o)
    }

    /// What a read of a CRAM gave: the fields of the records it got to, and the error that
    /// ended it, if any.
    type Run = (Vec<FlagFields>, Option<String>);

    /// The fields of every record of the record path (the full decode), read in place as
    /// flagstat read them before it had a fields path.
    fn full_run(name: &str, bytes: &[u8], reference: Option<&Path>) -> Run {
        let mut out = Vec::new();
        let end =
            source(name, Box::new(Cursor::new(bytes.to_vec())), reference, 1).and_then(|mut s| {
                let mut buf = Vec::new();
                while s.next_record(&mut buf)?.is_some() {
                    out.push(FlagFields::of_bam(&buf).unwrap());
                    buf.clear();
                }
                Ok(())
            });
        (out, end.err().map(|e| e.to_string()))
    }

    /// The fields path's run at `threads`, and the containers it read alone and in full.
    fn fields_run(
        name: &str,
        bytes: &[u8],
        reference: Option<&Path>,
        threads: usize,
    ) -> (Run, (u64, u64)) {
        let mut out = Vec::new();
        let mut counts = (0, 0);
        let end = source(
            name,
            Box::new(Cursor::new(bytes.to_vec())),
            reference,
            threads,
        )
        .and_then(|mut s| {
            let end = s.for_each_flag_fields(&mut |r| {
                out.push(r);
                Ok(())
            });
            counts = s.fields_counts();
            end
        });
        ((out, end.err().map(|e| e.to_string())), counts)
    }

    /// Every CRAM file among the fixtures and the goldens, once each (by content).
    fn fixtures() -> Vec<PathBuf> {
        fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
            for e in std::fs::read_dir(dir).unwrap() {
                let p = e.unwrap().path();
                if p.is_dir() {
                    walk(&p, out);
                } else if p.extension().is_some_and(|x| x == "cram") {
                    out.push(p);
                }
            }
        }
        let mut all = Vec::new();
        walk(&root().join("xtask/diff/inputs"), &mut all);
        walk(&root().join("crates/gloci/tests/golden"), &mut all);
        all.sort();
        let mut seen = HashSet::new();
        all.retain(|p| seen.insert(std::fs::read(p).unwrap()));
        all
    }

    /// The FASTA a fixture was written against (its @SQ lines check against it), if any.
    fn reference_for(name: &str, bytes: &[u8]) -> Option<PathBuf> {
        [fasta("ref.fa"), fasta("manyref.fa")]
            .into_iter()
            .find(|r| source(name, Box::new(Cursor::new(bytes.to_vec())), Some(r), 1).is_ok())
    }

    /// What the fixtures hold that the fields path has to get right.
    #[derive(Debug, Default)]
    struct Census {
        downstream: u64,
        detached: u64,
        multi_ref_slices: u64,
        unmapped: u64,
        unplaced: u64,
        secondary: u64,
        supplementary: u64,
        v30: bool,
        v31: bool,
    }

    /// Count, in the CRAM `bytes`, what [`Census`] counts, from the series as stored (BF, CF), in
    /// the containers whose fields can be read alone.
    fn census(bytes: &[u8], c: &mut Census) {
        c.v30 |= bytes[4..6] == [3, 0];
        c.v31 |= bytes[4..6] == [3, 1];
        let mut r = Cursor::new(bytes.to_vec());
        let Ok(head) = read_head(&mut r, "x") else {
            return;
        };
        let mut frames = CramFrames::new("x", Box::new(r), &head);
        while let Ok(Some(cont)) = frames.next() {
            let Ok(heads) = slice_heads(&cont, "x", 1 << 20) else {
                continue;
            };
            let comp = compression(&cont);
            let plan = Plan::new(&comp);
            let Ok(framed) = reframe(&cont, &heads, &plan, None) else {
                continue;
            };
            let mut container = cram::io::reader::Container::default();
            if cram::io::Reader::new(&framed.bytes[..])
                .read_container(&mut container)
                .is_err()
            {
                continue;
            }
            for (head, slice) in heads.iter().zip(container.slices()) {
                let Ok((_, ext)) = slice.and_then(|s| s.decode_blocks()) else {
                    continue;
                };
                let (Ok(mut bf), Ok(mut cf)) = (
                    Series::of(plan.bf, &ext, "BF"),
                    Series::of(plan.cf, &ext, "CF"),
                ) else {
                    continue;
                };
                c.multi_ref_slices += u64::from(head.ref_id == -2);
                for _ in 0..head.records {
                    let (Ok(flag), Ok(cram_flags)) = (bf.next(), cf.next()) else {
                        break;
                    };
                    let cram_flags = cram_flags as u8;
                    c.downstream += u64::from(cram_flags & (DETACHED | DOWNSTREAM) == DOWNSTREAM);
                    c.detached += u64::from(cram_flags & DETACHED != 0);
                    c.unmapped += u64::from(flag & 0x4 != 0);
                    c.unplaced += u64::from(head.ref_id == -1);
                    c.secondary += u64::from(flag & 0x100 != 0);
                    c.supplementary += u64::from(flag & 0x800 != 0);
                }
            }
        }
    }

    /// Every CRAM fixture of the repo: the fields path gives the full decode's fields, record by
    /// record, and the same error after the same records, at -t 1, 2 and 8; it reads every
    /// container samtools wrote alone. The fixtures hold mates further on in a slice and detached
    /// ones, unmapped and unplaced reads, slices of several references, secondary and
    /// supplementary reads, in CRAM 3.0 and 3.1.
    #[test]
    fn cram_flag_fields_match_full_decode() {
        let overhang = std::fs::read(aln("overhang.cram")).unwrap();
        let mut c = Census::default();
        let mut checked = 0;
        for path in fixtures() {
            let name = path.to_string_lossy().into_owned();
            let bytes = std::fs::read(&path).unwrap();
            census(&bytes, &mut c);
            let reference = reference_for(&name, &bytes);
            let want = full_run(&name, &bytes, reference.as_deref());
            for t in [1, 2, 8] {
                let (got, (_, full)) = fields_run(&name, &bytes, reference.as_deref(), t);
                if bytes == overhang {
                    // A read 1.5 Mb past the end of chr1 is past gloci's limit for decoding its
                    // bases (#cram-limits); flagstat reads none, so it counts both reads, as
                    // samtools does.
                    assert!(
                        want.1
                            .as_deref()
                            .is_some_and(|e| e.contains("gloci's limit")),
                        "{name}: {want:?}"
                    );
                    assert_eq!((got.0.len(), got.1.as_deref()), (2, None), "{name} -t {t}");
                    continue;
                }
                assert_eq!(got, want, "{name} -t {t}");
                if want.1.is_none() {
                    assert_eq!(full, 0, "{name} -t {t}: containers decoded in full");
                }
            }
            checked += 1;
        }
        assert!(checked >= 40, "{checked} CRAM fixtures");
        for (what, n) in [
            ("mates further on in a slice", c.downstream),
            ("detached mates", c.detached),
            ("slices of several references", c.multi_ref_slices),
            ("unmapped reads", c.unmapped),
            ("unplaced reads", c.unplaced),
            ("secondary reads", c.secondary),
            ("supplementary reads", c.supplementary),
        ] {
            assert!(n > 0, "no fixture has {what}: {c:?}");
        }
        assert!(c.v30 && c.v31, "{c:?}");
    }

    // ---- crafted encodings ----

    /// A CRAM whose blocks are all stored raw (`level=0`), split up to be rewritten: the file
    /// definition and header container, each data container, and the end-of-file container.
    struct RawCram {
        head: Vec<u8>,
        containers: Vec<RawContainerParts>,
        eof: Vec<u8>,
    }

    struct RawContainerParts {
        /// The container header's fields between its length and its landmarks.
        mid: Vec<u8>,
        /// The compression header's preservation map, data series encodings (key, encoding)
        /// and tag encoding map.
        preservation: Vec<u8>,
        series: Vec<([u8; 2], Vec<u8>)>,
        tags: Vec<u8>,
        slices: Vec<RawSlice>,
    }

    struct RawSlice {
        /// The slice header block, as stored.
        header: Vec<u8>,
        records: usize,
        /// (content type, content id, contents) of each block after the header, core first.
        blocks: Vec<(u8, i32, Vec<u8>)>,
    }

    fn put_raw_block(out: &mut Vec<u8>, content_type: u8, id: i32, data: &[u8]) {
        let start = out.len();
        out.extend_from_slice(&[0, content_type]);
        put_itf8(id as u32, out);
        put_itf8(data.len() as u32, out);
        put_itf8(data.len() as u32, out);
        out.extend_from_slice(data);
        let crc = crc32(&out[start..]);
        out.extend_from_slice(&crc.to_le_bytes());
    }

    impl RawCram {
        fn parse(bytes: &[u8]) -> RawCram {
            let head_len = read_head(&mut &bytes[..], "x").unwrap().len as usize;
            let mut rest = &bytes[head_len..];
            let mut containers = Vec::new();
            loop {
                let all = rest;
                let len = u32_le_at(rest) as usize;
                let mut s = &rest[4..];
                let mid_start = s;
                for _ in 0..4 {
                    itf8(&mut s).unwrap();
                }
                ltf8(&mut s).unwrap();
                ltf8(&mut s).unwrap();
                itf8(&mut s).unwrap();
                let mid = mid_start[..mid_start.len() - s.len()].to_vec();
                let n = itf8(&mut s).unwrap() as usize;
                let landmarks: Vec<usize> =
                    (0..n).map(|_| itf8(&mut s).unwrap() as usize).collect();
                let s = &s[4..];
                if n == 0 && len == 15 {
                    return RawCram {
                        head: bytes[..head_len].to_vec(),
                        containers,
                        eof: all.to_vec(),
                    };
                }
                let data = &s[..len];
                rest = &s[len..];
                let mut d = data;
                let comp = read_block(&mut d).unwrap();
                assert_eq!(comp.method, 0, "the fixture is stored raw");
                let mut m = comp.data;
                let mut map = || {
                    let n = itf8(&mut m).unwrap() as usize;
                    let (x, r) = m.split_at(n);
                    m = r;
                    x.to_vec()
                };
                let (preservation, ds, tags) = (map(), map(), map());
                let mut e = &ds[..];
                let n = itf8(&mut e).unwrap();
                let mut series = Vec::new();
                for _ in 0..n {
                    let key = [e[0], e[1]];
                    e = &e[2..];
                    let start = e;
                    itf8(&mut e).unwrap();
                    let len = itf8(&mut e).unwrap() as usize;
                    e = &e[len..];
                    series.push((key, start[..start.len() - e.len()].to_vec()));
                }
                let heads = slice_heads_of(data, &landmarks);
                let slices = heads
                    .into_iter()
                    .map(|(mut s, records)| {
                        let all = s;
                        read_block(&mut s).unwrap();
                        let header = all[..all.len() - s.len()].to_vec();
                        let mut blocks = Vec::new();
                        while !s.is_empty() {
                            let b = read_block(&mut s).unwrap();
                            assert_eq!(b.method, 0, "the fixture is stored raw");
                            blocks.push((b.content_type, b.content_id, b.data.to_vec()));
                        }
                        RawSlice {
                            header,
                            records,
                            blocks,
                        }
                    })
                    .collect();
                containers.push(RawContainerParts {
                    mid,
                    preservation,
                    series,
                    tags,
                    slices,
                });
            }
        }

        fn write(&self) -> Vec<u8> {
            let mut out = self.head.clone();
            for c in &self.containers {
                let mut ds = Vec::new();
                put_itf8(c.series.len() as u32, &mut ds);
                for (key, enc) in &c.series {
                    ds.extend_from_slice(key);
                    ds.extend_from_slice(enc);
                }
                let mut comp = Vec::new();
                for m in [&c.preservation, &ds, &c.tags] {
                    put_itf8(m.len() as u32, &mut comp);
                    comp.extend_from_slice(m);
                }
                let mut data = Vec::new();
                put_raw_block(&mut data, 1, 0, &comp);
                let mut landmarks = Vec::new();
                for s in &c.slices {
                    landmarks.push(data.len());
                    data.extend_from_slice(&s.header);
                    for (t, id, d) in &s.blocks {
                        put_raw_block(&mut data, *t, *id, d);
                    }
                }
                let start = out.len();
                out.extend_from_slice(&(data.len() as u32).to_le_bytes());
                out.extend_from_slice(&c.mid);
                put_itf8(landmarks.len() as u32, &mut out);
                for l in landmarks {
                    put_itf8(l as u32, &mut out);
                }
                let crc = crc32(&out[start..]);
                out.extend_from_slice(&crc.to_le_bytes());
                out.extend_from_slice(&data);
            }
            out.extend_from_slice(&self.eof);
            out
        }
    }

    fn u32_le_at(b: &[u8]) -> u32 {
        u32::from_le_bytes([b[0], b[1], b[2], b[3]])
    }

    /// Each slice of a container's `data` (slices at `landmarks`), and its record count.
    fn slice_heads_of<'a>(data: &'a [u8], landmarks: &[usize]) -> Vec<(&'a [u8], usize)> {
        landmarks
            .iter()
            .enumerate()
            .map(|(i, &at)| {
                let end = landmarks.get(i + 1).copied().unwrap_or(data.len());
                let mut s = &data[at..end];
                let h = read_block(&mut s).unwrap();
                let mut c = h.data;
                for _ in 0..3 {
                    itf8(&mut c).unwrap();
                }
                (&data[at..end], itf8(&mut c).unwrap() as usize)
            })
            .collect()
    }

    /// How a data series is coded in a container.
    #[derive(Clone, Copy)]
    enum Coding {
        External(i32),
        Constant(i32),
    }

    impl RawContainerParts {
        /// How data series `key` is coded: EXTERNAL, or a one-symbol Huffman table.
        fn coding(&self, key: &[u8; 2]) -> Coding {
            let (_, enc) = self.series.iter().find(|(k, _)| k == key).unwrap();
            let mut e = &enc[..];
            let codec = itf8(&mut e).unwrap();
            itf8(&mut e).unwrap();
            match codec {
                1 => Coding::External(itf8(&mut e).unwrap()),
                3 => {
                    assert_eq!(itf8(&mut e), Some(1), "a one-symbol table");
                    Coding::Constant(itf8(&mut e).unwrap())
                }
                _ => panic!("{} uses codec {codec}", String::from_utf8_lossy(key)),
            }
        }

        fn set_encoding(&mut self, key: &[u8; 2], codec: u32, params: &[u8]) {
            let mut enc = Vec::new();
            put_itf8(codec, &mut enc);
            put_itf8(params.len() as u32, &mut enc);
            enc.extend_from_slice(params);
            self.series.iter_mut().find(|(k, _)| k == key).unwrap().1 = enc;
        }
    }

    impl RawSlice {
        fn block(&self, content_type: u8, id: i32) -> &[u8] {
            let found = self
                .blocks
                .iter()
                .find(|(t, i, _)| *t == content_type && *i == id);
            &found.unwrap().2
        }

        fn set_block(&mut self, content_type: u8, id: i32, data: Vec<u8>) {
            self.blocks
                .iter_mut()
                .find(|(t, i, _)| *t == content_type && *i == id)
                .unwrap()
                .2 = data;
        }

        /// The values of a series coded as `coding`, for the `n` reads that read it.
        fn read(&self, coding: Coding, n: usize) -> Vec<i32> {
            let v = match coding {
                Coding::External(id) => itf8_values(self.block(EXTERNAL, id)),
                Coding::Constant(v) => vec![v; n],
            };
            assert_eq!(v.len(), n);
            v
        }

        /// Each read's BF, and whether it reads MQ (it is mapped).
        fn flags(&self, c: &RawContainerParts) -> Vec<(i32, bool)> {
            let bf = self.read(c.coding(b"BF"), self.records);
            bf.into_iter().map(|f| (f, f & 0x4 == 0)).collect()
        }

        /// The values of `key` (BF or MQ), for the reads that read it.
        fn values_of(&self, c: &RawContainerParts, key: &[u8; 2]) -> Vec<i32> {
            let flags = self.flags(c);
            match key {
                b"BF" => flags.into_iter().map(|(f, _)| f).collect(),
                _ => self.read(c.coding(key), flags.iter().filter(|(_, m)| *m).count()),
            }
        }
    }

    fn itf8_values(mut d: &[u8]) -> Vec<i32> {
        let mut v = Vec::new();
        while !d.is_empty() {
            v.push(itf8(&mut d).unwrap());
        }
        v
    }

    fn itf8_bytes(values: impl IntoIterator<Item = i32>) -> Vec<u8> {
        let mut out = Vec::new();
        for v in values {
            put_itf8(v as u32, &mut out);
        }
        out
    }

    /// Whether a container reads series `key`: every one reads BF, one with a mapped read MQ.
    fn reads(c: &RawContainerParts, key: &[u8; 2]) -> bool {
        key == b"BF"
            || c.slices
                .iter()
                .any(|s| s.flags(c).iter().any(|&(_, mapped)| mapped))
    }

    /// Series `key` (BF or MQ) coded in the core block: a Huffman table of every value it has in
    /// the container, codes of one length (two symbols at least, so each read takes bits).
    /// Returns whether the container reads it there.
    fn in_core(c: &mut RawContainerParts, key: &[u8; 2]) -> bool {
        let reading = reads(c, key);
        let per_slice: Vec<Vec<i32>> = c.slices.iter().map(|s| s.values_of(c, key)).collect();
        let mut symbols: Vec<i32> = per_slice.iter().flatten().copied().collect();
        symbols.sort_unstable();
        symbols.dedup();
        while symbols.len() < 2 {
            symbols.push(symbols.last().map_or(0, |v| v + 1));
        }
        let bits = usize::BITS - (symbols.len() - 1).leading_zeros();
        let mut params = itf8_bytes([symbols.len() as i32]);
        params.extend(itf8_bytes(symbols.iter().copied()));
        params.extend(itf8_bytes([symbols.len() as i32]));
        params.extend(itf8_bytes(symbols.iter().map(|_| bits as i32)));
        c.set_encoding(key, 3, &params);
        for (s, values) in c.slices.iter_mut().zip(per_slice) {
            let (content_type, id, core) = &s.blocks[0];
            assert!(
                *content_type == CORE && core.is_empty(),
                "the core block is empty"
            );
            let id = *id;
            let (mut bytes, mut acc, mut n) = (Vec::new(), 0u32, 0u32);
            for v in values {
                let code = symbols.binary_search(&v).unwrap() as u32;
                for k in (0..bits).rev() {
                    acc = acc << 1 | (code >> k & 1);
                    n += 1;
                    if n == 8 {
                        bytes.push(acc as u8);
                        (acc, n) = (0, 0);
                    }
                }
            }
            if n > 0 {
                bytes.push((acc << (8 - n)) as u8);
            }
            s.set_block(CORE, id, bytes);
        }
        reading
    }

    /// Series `key` (BF or MQ) read from AP's block, which flagstat skips: each read's BF then
    /// its AP, or its AP then (a mapped read) its MQ, as noodles reads them. A container whose
    /// AP is a constant is left as it is. Returns whether the container reads the series there.
    fn in_ap_block(c: &mut RawContainerParts, key: &[u8; 2]) -> bool {
        let Coding::External(ap) = c.coding(b"AP") else {
            return false;
        };
        let reading = reads(c, key);
        let per_slice: Vec<Vec<i32>> = c.slices.iter().map(|s| s.values_of(c, key)).collect();
        let flags: Vec<Vec<(i32, bool)>> = c.slices.iter().map(|s| s.flags(c)).collect();
        c.set_encoding(key, 1, &itf8_bytes([ap]));
        for ((s, values), flags) in c.slices.iter_mut().zip(per_slice).zip(flags) {
            let aps = itf8_values(s.block(EXTERNAL, ap));
            let mut values = values.into_iter();
            let mut data = Vec::new();
            for (a, (_, mapped)) in aps.into_iter().zip(flags) {
                if key == b"BF" {
                    put_itf8(values.next().unwrap() as u32, &mut data);
                }
                put_itf8(a as u32, &mut data);
                if key == b"MQ" && mapped {
                    put_itf8(values.next().unwrap() as u32, &mut data);
                }
            }
            assert!(values.next().is_none());
            s.set_block(EXTERNAL, ap, data);
        }
        reading
    }

    /// BF declared twice: first EXTERNAL in its own block, whose values are now wrong (every read
    /// unmapped), then in the core block. noodles keeps a key's last encoding, so the full decode
    /// reads the core; the fields path must not read the stale block. Returns whether the
    /// container reads BF where the fields path cannot read it alone.
    fn bf_declared_twice(c: &mut RawContainerParts, key: &[u8; 2]) -> bool {
        let Coding::External(stale) = c.coding(key) else {
            return in_core(c, key);
        };
        let first = c.series.iter().find(|(k, _)| k == key).unwrap().clone();
        let reading = in_core(c, key);
        for s in &mut c.slices {
            let wrong = itf8_bytes(std::iter::repeat_n(0x4, s.records));
            s.set_block(EXTERNAL, stale, wrong);
        }
        c.series.insert(0, first);
        reading
    }

    /// RN, the read names (BYTE_ARRAY_STOP), read from BF's block: each read's BF then its name
    /// and the stop byte, as noodles reads them. A container whose BF is a constant is left as
    /// it is. Returns whether the container reads BF from a block another series reads.
    fn rn_in_bf_block(c: &mut RawContainerParts, key: &[u8; 2]) -> bool {
        let Coding::External(bf) = c.coding(key) else {
            return false;
        };
        let (_, enc) = c.series.iter().find(|(k, _)| k == b"RN").unwrap();
        let mut e = &enc[..];
        assert_eq!(itf8(&mut e), Some(5), "RN is BYTE_ARRAY_STOP");
        itf8(&mut e).unwrap();
        let stop = e[0];
        e = &e[1..];
        let rn = itf8(&mut e).unwrap();
        let mut params = vec![stop];
        put_itf8(bf as u32, &mut params);
        c.set_encoding(b"RN", 5, &params);
        for s in &mut c.slices {
            let flags = itf8_values(s.block(EXTERNAL, bf));
            let names = s.block(EXTERNAL, rn).to_vec();
            let names: Vec<&[u8]> = names.split(|&b| b == stop).take(s.records).collect();
            assert_eq!(names.len(), flags.len());
            let mut data = Vec::new();
            for (flag, name) in flags.into_iter().zip(names) {
                put_itf8(flag as u32, &mut data);
                data.extend_from_slice(name);
                data.push(stop);
            }
            s.set_block(EXTERNAL, bf, data);
        }
        true
    }

    /// A CRAM whose BF or MQ is in the core block, or shares a block with a series flagstat skips
    /// (AP, the read names), or whose BF is declared twice (noodles keeps the last): every
    /// container that reads the series there is decoded in full, and the fields are those of the
    /// file as samtools wrote it, at -t 1 and 2.
    #[test]
    fn cram_flag_fields_fall_back() {
        let bytes = std::fs::read(aln("mates.raw.cram")).unwrap();
        let reference = fasta("ref.fa");
        let want = full_run("mates.raw.cram", &bytes, Some(&reference));
        assert!(want.1.is_none() && want.0.len() > 100, "{want:?}");
        let (got, (fast, full)) = fields_run("mates.raw.cram", &bytes, Some(&reference), 1);
        assert_eq!(got, want);
        assert!(fast > 1 && full == 0, "{fast} {full}");
        type Edit = fn(&mut RawContainerParts, &[u8; 2]) -> bool;
        let edits: [(&str, Edit, &[u8; 2]); 6] = [
            ("BF in the core block", in_core, b"BF"),
            ("MQ in the core block", in_core, b"MQ"),
            ("BF in AP's block", in_ap_block, b"BF"),
            ("MQ in AP's block", in_ap_block, b"MQ"),
            ("RN in BF's block", rn_in_bf_block, b"BF"),
            (
                "BF declared twice, last in the core block",
                bf_declared_twice,
                b"BF",
            ),
        ];
        for (what, edit, key) in edits {
            let mut raw = RawCram::parse(&bytes);
            assert_eq!(raw.write(), bytes, "the raw CRAM is rewritten as it was");
            // The containers that read the series where it now is: a container that never reads
            // it (MQ without a mapped read) still has its fields read alone.
            let reading = raw
                .containers
                .iter_mut()
                .map(|c| edit(c, key))
                .filter(|&r| r)
                .count() as u64;
            assert!(reading > 1, "{what}");
            let crafted = raw.write();
            assert_eq!(
                full_run("crafted.cram", &crafted, Some(&reference)),
                want,
                "{what}"
            );
            for t in [1, 2] {
                let (got, counts) = fields_run("crafted.cram", &crafted, Some(&reference), t);
                assert_eq!(got, want, "{what} -t {t}");
                assert_eq!(
                    counts,
                    (fast - reading, reading),
                    "{what} -t {t}: the containers that read {} are decoded in full",
                    String::from_utf8_lossy(key)
                );
            }
        }
    }

    /// A run cancelled while flagstat reads its CRAM (by the visitor, at the first record, so
    /// every byte is in memory by then) stops with the cancellation: at -t 1 the loop polls the
    /// run's token at each container, at -t 2 the pool's workers poll it too.
    #[test]
    fn cram_flag_fields_are_cancelled() {
        let bytes = std::fs::read(aln("mates.noref.cram")).unwrap();
        for t in [1, 2] {
            let token = RunToken::new();
            let _run = RunGuard::enter(&token);
            let mut s =
                source("mates.cram", Box::new(Cursor::new(bytes.clone())), None, t).unwrap();
            let mut seen = 0u64;
            let e = s
                .for_each_flag_fields(&mut |_| {
                    token.cancel();
                    seen += 1;
                    Ok(())
                })
                .unwrap_err();
            assert!(e.is_cancelled(), "-t {t}: {e} after {seen} records");
            assert!(seen < 100, "-t {t}: {seen} records after the cancel");
        }
    }
}
