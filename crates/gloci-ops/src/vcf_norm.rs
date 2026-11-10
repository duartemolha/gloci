//! `gloci vcf norm`: left-align and trim VCF/BCF records against a reference
//! FASTA, as `bcftools norm -f` does.

use crate::vcf_select::vcf_header;
use gloci_core::{ChromTable, Interval};
use gloci_io::{
    BcfText, Format, IndexedFasta, IoError, Loc, RecordBatch, RecordReader, RecordWriter,
    TextToBcf, next_batch_size,
};
use gloci_vcf::{
    CheckRef, Dedup, Normalizer, Reference, RmDup, SplitOutcome, SplitTypes, Splitter, Verdict,
    clean_reference_bases,
};
use std::collections::VecDeque;

pub use crate::vcf_split::vcf_rewrite_output as vcf_norm_output;

const COMMAND: &str = "vcf norm";

/// Bases read from the FASTA at a time.
const WINDOW_BASES: u64 = 1 << 16;

/// Options for `vcf norm`.
#[derive(Clone, Copy, Debug)]
pub struct NormOptions {
    /// Split multi-allelic records first (`bcftools norm -m-`).
    pub split: Option<SplitTypes>,
    pub check_ref: CheckRef,
    pub rm_dup: Option<RmDup>,
    /// Records that move left are re-sorted within this many bases.
    pub window: u64,
    /// Skip a record that cannot be processed, with a warning, instead of failing.
    pub skip_invalid: bool,
}

/// Records read, written and changed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NormSummary {
    pub records_in: u64,
    pub records_split: u64,
    pub realigned: u64,
    pub excluded: u64,
    pub duplicates: u64,
    pub records_out: u64,
    pub records_skipped: u64,
    /// The most reference bases held at once (bounded, whatever the contig size).
    pub peak_reference_bases: u64,
}

/// A FASTA read through one bounded window of bases.
struct FastaRef {
    fa: IndexedFasta,
    chrom: Vec<u8>,
    lo: u64,
    buf: Vec<u8>,
    peak: u64,
    fatal: Option<IoError>,
}

impl FastaRef {
    fn load(&mut self, chrom: &[u8], a: u64, b: u64) -> Result<(), IoError> {
        self.buf.clear();
        self.chrom.clear();
        self.lo = a;
        let buf = &mut self.buf;
        self.fa
            .fetch(chrom, Interval { start: a, end: b }, &mut |c| {
                clean_reference_bases(c, buf);
                Ok(())
            })?;
        self.chrom.extend_from_slice(chrom);
        self.peak = self.peak.max(self.buf.len() as u64);
        Ok(())
    }
}

impl Reference for FastaRef {
    fn get(
        &mut self,
        chrom: &[u8],
        start0: u64,
        len: u64,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        let Some(clen) = self.fa.len_of(chrom) else {
            let msg = format!(
                "the contig {} is not in the FASTA index of {}; use the reference the VCF was \
called against (the contigs of the FASTA are listed in its .fai)",
                String::from_utf8_lossy(chrom),
                self.fa.name()
            );
            self.fatal = Some(IoError::Detect {
                path: self.fa.name().to_string(),
                msg: msg.clone(),
            });
            return Err(msg);
        };
        if start0 >= clen {
            return Ok(());
        }
        let end = (start0 + len).min(clen);
        let held =
            self.chrom == chrom && start0 >= self.lo && end <= self.lo + self.buf.len() as u64;
        if !held {
            let (a, b) = if end - start0 >= WINDOW_BASES {
                (start0, end)
            } else if self.chrom == chrom && start0 < self.lo {
                (end.saturating_sub(WINDOW_BASES), end)
            } else {
                (start0, (start0 + WINDOW_BASES).min(clen))
            };
            if let Err(e) = self.load(chrom, a, b) {
                let msg = e.to_string();
                self.fatal = Some(e);
                return Err(msg);
            }
        }
        out.extend_from_slice(&self.buf[(start0 - self.lo) as usize..(end - self.lo) as usize]);
        Ok(())
    }
}

struct Pending {
    chrom: Vec<u8>,
    pos: u64,
    line: Vec<u8>,
    loc: u64,
}

struct Sink<'a> {
    w: &'a mut dyn RecordWriter,
    text_to_bcf: Option<TextToBcf>,
    enc: Vec<u8>,
    dedup: Option<Dedup>,
    sum: NormSummary,
}

impl Sink<'_> {
    fn emit(&mut self, p: &Pending) -> Result<(), IoError> {
        if let Some(d) = self.dedup.as_mut()
            && !d.keep(&p.line)
        {
            self.sum.duplicates += 1;
            return Ok(());
        }
        let write_err = |source| IoError::Write {
            path: "output".to_string(),
            source,
        };
        if let Some(t) = self.text_to_bcf.as_mut() {
            self.enc.clear();
            t.encode(&p.line, p.loc, &mut self.enc)?;
            self.w.write_record(&self.enc).map_err(write_err)?;
        } else {
            self.w.write_record(&p.line).map_err(write_err)?;
        }
        self.sum.records_out += 1;
        Ok(())
    }
}

/// Normalise every record of `input` against `fasta` and write the records
/// to `w` (opened with [`vcf_norm_output`]'s format and header), sorted by
/// position within `--window` bases as `bcftools norm` does.
pub fn vcf_norm(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    fasta: IndexedFasta,
    o: &NormOptions,
    w: &mut dyn RecordWriter,
) -> Result<NormSummary, IoError> {
    let h = vcf_header(COMMAND, input.as_ref())?;
    let name = input.name().to_string();
    let splitter = o.split.map(|t| Splitter::new(&h.raw, t, h.samples.len()));
    let norm = Normalizer::new(o.check_ref);
    let mut refs = FastaRef {
        fa: fasta,
        chrom: Vec::new(),
        lo: 0,
        buf: Vec::new(),
        peak: 0,
        fatal: None,
    };
    let from_bcf = input.format() == Format::Bcf;
    let mut bcf_text = if from_bcf {
        Some(BcfText::for_header(&name, &h)?)
    } else {
        None
    };
    let text_to_bcf = if w.is_bcf() {
        Some(TextToBcf::new(&name, &h)?)
    } else {
        None
    };
    let mut sink = Sink {
        w,
        text_to_bcf,
        enc: Vec::new(),
        dedup: o.rm_dup.map(Dedup::new),
        sum: NormSummary::default(),
    };
    let mut buf: VecDeque<Pending> = VecDeque::new();
    let (mut sout, mut ends, mut nout) = (Vec::new(), Vec::new(), Vec::new());
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for i in 0..batch.len() {
            let raw = batch.raw(i);
            let loc = batch.recs[i].loc;
            sink.sum.records_in += 1;
            let text = match bcf_text.as_mut() {
                Some(t) => t.render(raw, loc)?,
                None => raw,
            };
            let parse_err = |msg: &str| {
                let l = if from_bcf {
                    Loc::Record(loc)
                } else {
                    Loc::Line(loc)
                };
                IoError::Parse {
                    path: name.clone(),
                    loc: l,
                    msg: format!("{msg} (pass --skip-invalid to skip such records)"),
                }
            };
            // Split first, then normalise each piece, as bcftools norm -m- -f.
            sout.clear();
            let outcome = match &splitter {
                Some(s) => match s.split(text, &mut sout, &mut ends) {
                    Ok(x) => x,
                    Err(msg) => {
                        if o.skip_invalid {
                            sink.sum.records_skipped += 1;
                            log::warn!("{COMMAND}: {name}: record {loc} skipped: {msg}");
                            continue;
                        }
                        return Err(parse_err(&msg));
                    }
                },
                None => SplitOutcome::Unchanged,
            };
            let mut pieces: Vec<&[u8]> = Vec::new();
            match outcome {
                SplitOutcome::Unchanged => pieces.push(text),
                SplitOutcome::Split(_) => {
                    sink.sum.records_split += 1;
                    let mut prev = 0;
                    for &e in &ends {
                        pieces.push(&sout[prev..e]);
                        prev = e;
                    }
                }
            }
            for piece in pieces {
                let verdict = match norm.process(piece, &mut refs, &mut nout) {
                    Ok(v) => v,
                    Err(msg) => {
                        if let Some(e) = refs.fatal.take() {
                            return Err(e);
                        }
                        if o.skip_invalid {
                            sink.sum.records_skipped += 1;
                            log::warn!("{COMMAND}: {name}: record {loc} skipped: {msg}");
                            continue;
                        }
                        return Err(parse_err(&msg));
                    }
                };
                let Verdict::Keep { pos, changed } = verdict else {
                    sink.sum.excluded += 1;
                    continue;
                };
                if changed {
                    sink.sum.realigned += 1;
                }
                let chrom = nout
                    .split(|&b| b == b'\t')
                    .next()
                    .unwrap_or_default()
                    .to_vec();
                // Write what is far enough behind, then insert after any
                // record at the same position.
                while let Some(head) = buf.front() {
                    if head.chrom != chrom || head.pos + o.window <= pos {
                        let head = buf.pop_front().expect("front exists");
                        sink.emit(&head)?;
                    } else {
                        break;
                    }
                }
                let mut at = buf.len();
                while at > 0 && buf[at - 1].pos > pos {
                    at -= 1;
                }
                buf.insert(
                    at,
                    Pending {
                        chrom,
                        pos,
                        line: nout.clone(),
                        loc,
                    },
                );
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    while let Some(head) = buf.pop_front() {
        sink.emit(&head)?;
    }
    if let Some(s) = input.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    sink.sum.peak_reference_bases = refs.peak;
    Ok(sink.sum)
}
