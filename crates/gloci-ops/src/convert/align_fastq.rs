//! BAM, CRAM and SAM to FASTQ (D166, bedtools `bamtofastq`, `-fq2` and `-tags`).
//!
//! Without a second file every record is written, in file order, under its own name. With one,
//! the records of a read name that are next to each other are the pair: READ1 goes to the first
//! file as `NAME/1` and READ2 to the second as `NAME/2`. With `tags` each pair is written from
//! its first record alone, whose `R2` and `Q2` tags hold the other end. A read on the reverse
//! strand is written reverse-complemented with its qualities reversed.

use crate::seq::revcomp_in_place;
use gloci_io::{AlnRef, AlnSource, IoError, SortOrder, UsageError, aln_ref};
use std::io;

const PAIRED: u16 = 0x1;
const REVERSE: u16 = 0x10;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;

#[derive(Debug, Clone, Copy, Default)]
pub struct FastqOptions {
    /// Pair the records by name and write the ends to two outputs (bedtools `-fq2`).
    pub paired: bool,
    /// With `paired`: take the second end from the `R2` and `Q2` tags of the first record
    /// (bedtools `-tags`).
    pub tags: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastqSummary {
    /// Reads written, or pairs with `paired`.
    pub written: u64,
    /// Read names with one record only, which were not written, and the first of them.
    pub no_mate: u64,
    pub first_no_mate: Option<String>,
    /// Read names with two records that are not a READ1 and a READ2 (with `tags`: not both
    /// marked paired), which were not written, and the first of them.
    pub not_a_pair: u64,
    pub first_not_a_pair: Option<String>,
}

/// Where the FASTQ records go: the second end of a pair may share the first one's output.
pub trait FastqSink {
    fn first(&mut self, record: &[u8]) -> io::Result<()>;
    fn second(&mut self, record: &[u8]) -> io::Result<()>;
}

/// The `R2` and `Q2` tags of a first record, or why they are missing.
type MateTags = Result<(Vec<u8>, Vec<u8>), String>;

/// A read held until its name group is complete.
#[derive(Default)]
struct Read {
    flag: u16,
    seq: Vec<u8>,
    /// Raw Phred values, `0xFF` for each base when QUAL is missing.
    qual: Vec<u8>,
    /// `R2` and `Q2` (with `tags`, for the first record of a name), or why they are missing.
    mate: Option<MateTags>,
}

impl Read {
    fn is_first(&self) -> bool {
        self.flag & READ1 != 0
    }
    fn is_second(&self) -> bool {
        self.flag & READ2 != 0
    }
    fn is_reverse(&self) -> bool {
        self.flag & REVERSE != 0
    }
}

/// `@name[/end]`, the bases and the qualities of a read, reverse-complemented when `reverse`.
fn push_record(
    out: &mut Vec<u8>,
    name: &[u8],
    end: Option<u8>,
    seq: &[u8],
    qual: &[u8],
    reverse: bool,
) {
    out.push(b'@');
    out.extend_from_slice(name);
    if let Some(e) = end {
        out.extend_from_slice(&[b'/', e]);
    }
    out.push(b'\n');
    let at = out.len();
    out.extend_from_slice(seq);
    if reverse {
        revcomp_in_place(&mut out[at..], false);
    }
    out.extend_from_slice(b"\n+\n");
    let at = out.len();
    if qual.first() == Some(&0xFF) {
        out.resize(at + qual.len(), b'!');
    } else {
        out.extend(qual.iter().map(|q| q.wrapping_add(33)));
        if reverse {
            out[at..].reverse();
        }
    }
    out.push(b'\n');
}

/// Qualities given as text (a `Q2` tag) are written as they are.
fn push_text_record(out: &mut Vec<u8>, name: &[u8], end: u8, seq: &[u8], qual: &[u8]) {
    out.push(b'@');
    out.extend_from_slice(name);
    out.extend_from_slice(&[b'/', end, b'\n']);
    out.extend_from_slice(seq);
    out.extend_from_slice(b"\n+\n");
    out.extend_from_slice(qual);
    out.push(b'\n');
}

fn plural(n: u64, one: &str, many: &str) -> String {
    format!(
        "{} {}",
        gloci_core::thousands(n),
        if n == 1 { one } else { many }
    )
}

struct Pairing<'a> {
    o: &'a FastqOptions,
    sink: &'a mut dyn FastqSink,
    name: Vec<u8>,
    group: Vec<Read>,
    /// How many of `group` are in use (the others keep their buffers for the next name).
    used: usize,
    sum: FastqSummary,
    a: Vec<u8>,
    b: Vec<u8>,
}

impl Pairing<'_> {
    fn slot(&mut self) -> &mut Read {
        if self.used == self.group.len() {
            self.group.push(Read::default());
        }
        let r = &mut self.group[self.used];
        self.used += 1;
        r.seq.clear();
        r.qual.clear();
        r.mate = None;
        r
    }

    fn name_string(&self) -> String {
        String::from_utf8_lossy(&self.name).into_owned()
    }

    fn skip_no_mate(&mut self) {
        self.sum.no_mate += 1;
        if self.sum.first_no_mate.is_none() {
            self.sum.first_no_mate = Some(self.name_string());
        }
    }

    fn skip_not_a_pair(&mut self) {
        self.sum.not_a_pair += 1;
        if self.sum.first_not_a_pair.is_none() {
            self.sum.first_not_a_pair = Some(self.name_string());
        }
    }

    /// Write the pair a name group makes, or count it as skipped. The group has at most two
    /// reads here.
    fn flush(&mut self, file: &str) -> Result<(), IoError> {
        let n = std::mem::take(&mut self.used);
        if n == 0 {
            return Ok(());
        }
        if n == 1 {
            self.skip_no_mate();
            return Ok(());
        }
        let g = &self.group;
        self.a.clear();
        self.b.clear();
        if self.o.tags {
            // Both records are marked paired; the first one carries both ends.
            if g[0].flag & g[1].flag & PAIRED == 0 {
                self.skip_not_a_pair();
                return Ok(());
            }
            let first = &g[0];
            let (mseq, mqual) = match first.mate.as_ref().expect("tags are read with the record") {
                Ok(m) => m,
                Err(msg) => {
                    return Err(IoError::Detect {
                        path: file.to_string(),
                        msg: msg.clone(),
                    });
                }
            };
            let (mine, other, me_end, other_end) = if first.is_first() {
                (&mut self.a, &mut self.b, b'1', b'2')
            } else {
                (&mut self.b, &mut self.a, b'2', b'1')
            };
            push_record(
                mine,
                &self.name,
                Some(me_end),
                &first.seq,
                &first.qual,
                first.is_reverse(),
            );
            push_text_record(other, &self.name, other_end, mseq, mqual);
        } else {
            let (x, y) = if g[0].is_first() && g[1].is_second() {
                (&g[0], &g[1])
            } else if g[0].is_second() && g[1].is_first() {
                (&g[1], &g[0])
            } else {
                self.skip_not_a_pair();
                return Ok(());
            };
            push_record(
                &mut self.a,
                &self.name,
                Some(b'1'),
                &x.seq,
                &x.qual,
                x.is_reverse(),
            );
            push_record(
                &mut self.b,
                &self.name,
                Some(b'2'),
                &y.seq,
                &y.qual,
                y.is_reverse(),
            );
        }
        self.sum.written += 1;
        let w = |source| IoError::Write {
            path: "output".into(),
            source,
        };
        // The order the ends are written in matters when they share an output: bedtools writes
        // the end that comes from the tags second.
        if self.o.tags && !self.group[0].is_first() {
            self.sink.second(&self.b).map_err(w)?;
            self.sink.first(&self.a).map_err(w)
        } else {
            self.sink.first(&self.a).map_err(w)?;
            self.sink.second(&self.b).map_err(w)
        }
    }
}

/// Write the reads of `src` as FASTQ on `sink`: each record under its own name, or with
/// `o.paired` the pairs of name-grouped records, READ1 to `first` and READ2 to `second` (D166).
/// Returns what was skipped; the caller says so.
pub fn alignments_to_fastq(
    mut src: Box<dyn AlnSource>,
    o: &FastqOptions,
    sink: &mut dyn FastqSink,
) -> Result<FastqSummary, IoError> {
    let name = src.name().to_string();
    let hdr = src.header().clone();
    if o.paired && hdr.sort_order == SortOrder::Coordinate {
        return Err(UsageError(format!(
            "{name} is sorted by coordinate (@HD SO:coordinate), and --fastq2 pairs the reads that \
are next to each other by name; sort it by name first: gloci sort -n {name} -o x.byname.bam"
        ))
        .into());
    }
    let layout = src.layout();
    let mut p = Pairing {
        o,
        sink,
        name: Vec::new(),
        group: Vec::new(),
        used: 0,
        sum: FastqSummary::default(),
        a: Vec::new(),
        b: Vec::new(),
    };
    let mut buf = Vec::new();
    let mut out = Vec::new();
    loop {
        buf.clear();
        if src.next_record(&mut buf)?.is_none() {
            break;
        }
        let fail = |src: &dyn AlnSource, msg: &str| IoError::Detect {
            path: name.clone(),
            msg: format!("{name} {}: {msg}", src.place()),
        };
        let a: AlnRef<'_> = aln_ref(layout, &buf, &hdr).map_err(|m| fail(&*src, &m))?;
        if !o.paired {
            out.clear();
            let (mut seq, qual) = (Vec::new(), a.quals());
            a.seq_into(&mut seq);
            push_record(
                &mut out,
                a.name(),
                None,
                &seq,
                &qual,
                a.flag() & REVERSE != 0,
            );
            p.sink.first(&out).map_err(|source| IoError::Write {
                path: "output".into(),
                source,
            })?;
            p.sum.written += 1;
            continue;
        }
        let flag = a.flag();
        // A name group is its first record, whatever its flags, and the later records of the
        // same name; without --tags an unpaired record is not a mate and is passed over.
        if p.used > 0 && a.name() == &p.name[..] {
            if !o.tags && flag & PAIRED == 0 {
                continue;
            }
        } else {
            p.flush(&name)?;
            p.name.clear();
            p.name.extend_from_slice(a.name());
        }
        if p.used == 2 {
            return Err(fail(
                &*src,
                &format!(
                    "the read name {} has more than two records, and --fastq2 needs one READ1 \
and one READ2 per name; drop the secondary and supplementary alignments with \
--exclude-flags 0x900 (and duplicates with 0x400 if they are marked)",
                    String::from_utf8_lossy(a.name())
                ),
            ));
        }
        // The first record of a name carries the mate's bases in its tags; the others are
        // never read for them.
        let tags = (o.tags && p.used == 0).then(|| {
            let text = |t: [u8; 2]| {
                a.str_tag(t).map(<[u8]>::to_vec).ok_or_else(|| {
                    let t = String::from_utf8_lossy(&t).into_owned();
                    format!(
                        "{name} {}: read {} has no {t} string tag, which --tags takes the mate \
from (the R2 and Q2 tags); drop --tags to pair the records by name instead",
                        src.place(),
                        String::from_utf8_lossy(a.name())
                    )
                })
            };
            Ok::<_, String>((text(*b"R2")?, text(*b"Q2")?))
        });
        let r = p.slot();
        r.flag = flag;
        a.seq_into(&mut r.seq);
        r.qual.extend_from_slice(&a.quals());
        r.mate = tags;
    }
    if o.paired {
        p.flush(&name)?;
    }
    let sum = p.sum;
    src.diagnostics().log_alignment_counts(&name);
    if sum.no_mate > 0 {
        log::warn!(
            "{name}: {} one record and no pair was written for them (the first is {}); --fastq2 \
needs the records of a name next to each other (gloci sort -n) and every mate present",
            plural(sum.no_mate, "read name has", "read names have"),
            sum.first_no_mate.as_deref().unwrap_or("?")
        );
    }
    if sum.not_a_pair > 0 {
        log::warn!(
            "{name}: {} two records that {} not {}, and no pair was written for them (the \
first is {})",
            plural(sum.not_a_pair, "read name has", "read names have"),
            "are",
            if o.tags {
                "both marked paired (flag 0x1)"
            } else {
                "a READ1 and a READ2 (flags 0x40 and 0x80)"
            },
            sum.first_not_a_pair.as_deref().unwrap_or("?")
        );
    }
    if let Some(s) = src.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    Ok(sum)
}
