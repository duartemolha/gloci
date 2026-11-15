//! `convert --bam`: BED, GFF and VCF features as BAM records (bedtools `bedtobam`), and BEDPE
//! pairs as two records each (bedtools `bedpetobam`).
//!
//! Each record is the SAM line bedtools 2.31.1's record reads as, encoded by the BAM writer's own
//! encoder: no bases or qualities, the MAPQ given, and bedtools' place, flags, CIGAR and mate
//! fields. The records are written in input order under an `SO:unsorted` header whose `@SQ`
//! lines are the genome's chromosomes in its order.

use crate::browser::for_each_feature;
use gloci_core::{ChromTable, Genome, OrderSource};
use gloci_io::{
    AlnHeader, BedpeBatch, BedpeReader, Format, IoError, Loc, RecordReader, RecordWriter,
    SortOrder, bed12_blocks, bedpe_name, column_count, field, sam_line_to_bam, unknown_chrom_err,
};
use std::io::Write;

/// What `convert --bam` writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BamOptions {
    /// The MAPQ of every record (bedtools `-mapq`; 255 by default).
    pub mapq: u8,
    /// BED12 input: each block is an `M` operation and each gap between blocks an `N`
    /// (bedtools `-bed12`).
    pub bed12: bool,
}

/// The header bedtools writes for `-g GENOME`: `@HD VN:1.0 SO:unsorted`, then one `@SQ` line
/// per chromosome of `genome`, in its order, with `AS:` naming the genome as `-g` gave it
/// (`assembly`). A name with a control character in it would break the line and is left out.
pub fn bam_header(genome: &Genome, assembly: &str) -> Result<AlnHeader, String> {
    let mut text = b"@HD\tVN:1.0\tSO:unsorted\n".to_vec();
    let assembly = (!assembly.is_empty() && !assembly.bytes().any(|b| b.is_ascii_control()))
        .then_some(assembly);
    let mut refs = Vec::with_capacity(genome.len());
    for i in 0..genome.len() {
        let name = genome.name(i);
        let shown = String::from_utf8_lossy(name);
        let len = genome.length(i).ok_or_else(|| {
            format!(
                "the genome gives no length for {shown}, and a BAM header needs one; pass -g a \
.fai or .genome file, or a built-in genome"
            )
        })?;
        if len > i32::MAX as u64 {
            return Err(format!(
                "{shown} is {len} bases long, and BAM holds references of up to 2147483647 bases"
            ));
        }
        text.extend_from_slice(b"@SQ\tSN:");
        text.extend_from_slice(name);
        if let Some(a) = assembly {
            text.extend_from_slice(b"\tAS:");
            text.extend_from_slice(a.as_bytes());
        }
        text.extend_from_slice(format!("\tLN:{len}\n").as_bytes());
        refs.push((name.to_vec(), len));
    }
    Ok(AlnHeader::new(text, refs, SortOrder::Unsorted))
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// Why a chromosome cannot be written, and the fix.
fn missing(genome: &OrderSource, chrom: &[u8]) -> String {
    let c = String::from_utf8_lossy(chrom);
    format!(
        "{}; a BAM record needs a reference of the header, which -g gives: add {c} to the genome, \
pass -g a genome that lists it, or rename the chromosomes with --chrom-map",
        unknown_chrom_err(genome, &c)
    )
}

/// The SAM line of one record and its BAM bytes, reused from record to record.
#[derive(Default)]
struct Line {
    sam: Vec<u8>,
    bam: Vec<u8>,
}

/// Where a record is placed: its reference name and 0-based position, or neither.
type Place<'a> = Option<(&'a [u8], i64)>;

impl Line {
    /// Encode `self.sam` and write the record. `fail` turns an encoding error into the input's.
    fn write(
        &mut self,
        header: &AlnHeader,
        w: &mut dyn RecordWriter,
        fail: impl Fn(String) -> IoError,
    ) -> Result<(), IoError> {
        self.bam.clear();
        sam_line_to_bam(&self.sam, header, &mut self.bam)
            .map_err(|e| fail(format!("cannot be written as BAM: {e}")))?;
        w.write_alignment(&self.sam, Some(&self.bam))
            .map_err(write_err)
    }

    /// Start the line: QNAME to CIGAR.
    fn start(&mut self, name: &[u8], flag: u16, place: Place, mapq: u8, cigar: &[u8]) {
        self.sam.clear();
        self.sam.extend_from_slice(name);
        let (rname, pos) = place.unwrap_or((&b"*"[..], -1));
        let _ = write!(self.sam, "\t{flag}\t");
        self.sam.extend_from_slice(rname);
        let _ = write!(self.sam, "\t{}\t{mapq}\t", pos + 1);
        self.sam.extend_from_slice(cigar);
    }

    /// End the line: RNEXT to QUAL, with no bases and no qualities.
    fn end(&mut self, own: Place, mate: Place, tlen: i64) {
        self.sam.push(b'\t');
        match mate {
            None => self.sam.extend_from_slice(b"*\t0"),
            Some((name, pos)) => {
                if own.is_some_and(|(n, _)| n == name) {
                    self.sam.push(b'=');
                } else {
                    self.sam.extend_from_slice(name);
                }
                let _ = write!(self.sam, "\t{}", pos + 1);
            }
        }
        let _ = write!(self.sam, "\t{tlen}\t*\t*");
    }
}

/// Write each feature of `input` (BED, GFF or VCF) as one BAM record to `w`, a BAM writer for
/// SAM-layout records under `header`. Returns the number of records.
///
/// bedtools' reading of a feature (shared with `igv`): a BED record needs a name (BED4 or
/// more), a zero-length BED record is one base wider on each side (only the base after it at
/// the start of a chromosome), a GFF record is its feature type at `start-1`, a VCF record is
/// `REF/ALT[_ID]` over its REF bases. The record is `NAME FLAG CHROM START+1 MAPQ (END-START)M
/// * 0 0 * *`, with FLAG 16 when the strand is `-`; with `o.bed12` the CIGAR is the blocks.
pub fn features_to_bam(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    header: &AlnHeader,
    genome: &OrderSource,
    o: &BamOptions,
    w: &mut dyn RecordWriter,
) -> Result<u64, IoError> {
    let path = input.name().to_string();
    let format = input.format();
    let mut line = Line::default();
    let mut cigar = Vec::new();
    let mut n = 0u64;
    for_each_feature("convert --bam", input, chroms, |f, loc| {
        let fail = |msg: String| IoError::Parse {
            path: path.clone(),
            loc: Loc::Line(loc),
            msg,
        };
        let name = f.name.as_deref().ok_or_else(|| {
            fail(
                "the record has 3 columns; a BAM record needs a read name, so the input must \
have a name column (BED4 or more)"
                    .into(),
            )
        })?;
        if header.ref_id(f.chrom).is_none() {
            return Err(fail(missing(genome, f.chrom)));
        }
        // A zero-length record at the start of a chromosome has no base before it.
        let start = f.start.max(0);
        cigar.clear();
        if o.bed12 {
            let bed = matches!(format, Format::Bed | Format::BedGraph);
            if !bed || f.columns != 12 {
                let what = if bed {
                    format!("this one has {}", f.columns)
                } else {
                    format!("this input is {}", format.name())
                };
                return Err(fail(format!(
                    "--bed12 needs BED12 records with 12 columns; {what}"
                )));
            }
            let blocks = bed12_blocks(f.raw)
                .map_err(|e| fail(format!("invalid BED12: {e}")))?
                .unwrap_or_default();
            let mut at = blocks.first().map_or(0, |b| b.start);
            for b in &blocks {
                if b.start > at {
                    let _ = write!(cigar, "{}N", b.start - at);
                }
                let _ = write!(cigar, "{}M", b.end - b.start);
                at = b.end;
            }
        } else {
            let _ = write!(cigar, "{}M", f.end - start);
        }
        let flag = if f.minus { 0x10 } else { 0 };
        let place = Some((f.chrom, start as i64));
        line.start(name, flag, place, o.mapq, &cigar);
        line.end(place, None, 0);
        line.write(header, w, fail)?;
        n += 1;
        Ok(())
    })?;
    Ok(n)
}

/// Write each pair of `input` as two BAM records to `w`, a BAM writer for SAM-layout records
/// under `header`: the first end as READ1, the second as READ2. Returns the number of records.
///
/// As bedtools: both records are paired (0x1) and take the MAPQ given; each has 0x10 when its
/// own strand is `-` and 0x20 when its mate's is; both have 0x2 when the ends are on one
/// chromosome with the first on `+` and the second on `-`; the mate fields point at the other
/// end, and TLEN on both is `start2 - start1` on one chromosome, else 0. The CIGAR of an end is
/// `(end-start)M`, `0M` when it is zero-length. An absent end (chromosome `.`) is an unmapped
/// record (0x4, and 0x8 on its mate) with no CIGAR and MAPQ 0, placed at its mate as the SAM
/// specification recommends, or unplaced when both ends are absent; TLEN is then 0. Every line
/// needs the first line's number of columns, as bedtools requires.
pub fn bedpe_to_bam(
    mut input: BedpeReader,
    chroms: &mut ChromTable,
    header: &AlnHeader,
    genome: &OrderSource,
    mapq: u8,
    w: &mut dyn RecordWriter,
) -> Result<u64, IoError> {
    let path = input.name().to_string();
    let mut batch = BedpeBatch::default();
    let mut line = Line::default();
    let mut cigar: [Vec<u8>; 2] = Default::default();
    let mut first_columns: Option<usize> = None;
    let mut n = 0u64;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, 4096)?;
        for i in 0..batch.len() {
            let m = batch.recs[i];
            let raw = batch.raw(i);
            let fail = |msg: String| IoError::Parse {
                path: path.clone(),
                loc: Loc::Line(m.loc),
                msg,
            };
            let columns = column_count(raw);
            let first = *first_columns.get_or_insert(columns);
            if columns != first {
                return Err(fail(format!(
                    "{columns} columns, but the first line has {first}; every line must have the \
same number of columns"
                )));
            }
            let name = bedpe_name(raw).ok_or_else(|| {
                fail(
                    "the line has 6 columns; a BAM record needs a read name, so the line must \
have a name column (column 7)"
                        .into(),
                )
            })?;
            // The strand columns as written (an absent end keeps its own).
            let strand = |k: usize| field(raw, 8 + k);
            let minus = [strand(0) == Some(b"-"), strand(1) == Some(b"-")];
            let mut ends: [Option<(&[u8], i64, i64)>; 2] = [None, None];
            for (k, e) in m.ends.iter().enumerate() {
                if let Some(e) = e {
                    let c = chroms.name(e.chrom);
                    if header.ref_id(c).is_none() {
                        return Err(fail(missing(genome, c)));
                    }
                    ends[k] = Some((c, e.span.start as i64, e.span.end as i64));
                }
            }
            let same = matches!((ends[0], ends[1]), (Some(a), Some(b)) if a.0 == b.0);
            let proper = same && strand(0) == Some(b"+") && strand(1) == Some(b"-");
            let tlen = match (ends[0], ends[1]) {
                (Some(a), Some(b)) if same => b.1 - a.1,
                _ => 0,
            };
            // An absent end sits at its mate; with both absent neither is placed.
            let place = |k: usize| ends[k].or(ends[1 - k]).map(|(c, s, _)| (c, s));
            for k in 0..2 {
                let o = 1 - k;
                let mut flag: u16 = 0x1 | if k == 0 { 0x40 } else { 0x80 };
                if minus[k] {
                    flag |= 0x10;
                }
                if minus[o] {
                    flag |= 0x20;
                }
                if proper {
                    flag |= 0x2;
                }
                if ends[k].is_none() {
                    flag |= 0x4;
                }
                if ends[o].is_none() {
                    flag |= 0x8;
                }
                cigar[k].clear();
                let q = match ends[k] {
                    Some((_, s, e)) => {
                        let _ = write!(cigar[k], "{}M", e - s);
                        mapq
                    }
                    None => {
                        cigar[k].push(b'*');
                        0
                    }
                };
                line.start(name, flag, place(k), q, &cigar[k]);
                line.end(place(k), place(o), tlen);
                line.write(header, w, fail)?;
                n += 1;
            }
        }
        if !more {
            return Ok(n);
        }
    }
}
