//! BAM, CRAM and SAM to BEDPE (D156, bedtools `bamtobed -bedpe` and `-mate1`): reads that are
//! adjacent by name make a pair, the memory held is the primary mates of one read name.

use gloci_core::Strand;
use gloci_io::{AlnHeader, AlnRef, AlnSource, IoError, SortOrder, UsageError, aln_ref};
use std::io::Write;

#[derive(Debug, Clone, Default)]
pub struct BedpeOptions {
    /// Keep READ1 as the first end (bedtools `-mate1`); otherwise the end that sorts first by
    /// `(chromosome as bytes, start)` comes first, an absent end before any.
    pub mate1_first: bool,
    /// The integer tag whose sum over the mapped ends is the score (`NM` is bedtools `-ed`);
    /// otherwise the smaller MAPQ when both ends are mapped, else 0.
    pub score_tag: Option<[u8; 2]>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BedpeSummary {
    /// Pairs written.
    pub pairs: u64,
    /// Read names with one primary mate only (not written).
    pub missing_mate: u64,
    /// Reads that are neither READ1 nor READ2 (not written).
    pub unpaired: u64,
}

/// One end of a pair: the reference and span of a mapped read, or an absent end.
#[derive(Clone, Copy)]
struct End {
    /// `Some((reference id, start, end))` for a mapped read.
    at: Option<(usize, u64, u64)>,
    minus: bool,
    mapq: u8,
    /// The `--score-tag` value (0 when unmapped or no tag asked for).
    tag: i64,
}

struct Group {
    name: Vec<u8>,
    mate1: Option<End>,
    mate2: Option<End>,
}

struct Pairing<'a> {
    o: &'a BedpeOptions,
    hdr: &'a AlnHeader,
    out: &'a mut dyn Write,
    row: Vec<u8>,
    sum: BedpeSummary,
}

impl Pairing<'_> {
    /// The sort key of an end: absent first, then chromosome bytes, then start.
    fn key(&self, e: &End) -> (bool, &[u8], u64) {
        match e.at {
            None => (false, &[][..], 0),
            Some((rid, s, _)) => (true, &self.hdr.refs[rid].0[..], s),
        }
    }

    fn end_columns(&mut self, e: &End) {
        match e.at {
            None => self.row.extend_from_slice(b".\t-1\t-1"),
            Some((rid, s, t)) => {
                self.row.extend_from_slice(&self.hdr.refs[rid].0);
                self.row.push(b'\t');
                self.row.extend_from_slice(s.to_string().as_bytes());
                self.row.push(b'\t');
                self.row.extend_from_slice(t.to_string().as_bytes());
            }
        }
    }

    fn strand(e: &End) -> u8 {
        match (e.at, e.minus) {
            (None, _) => b'.',
            (Some(_), true) => b'-',
            (Some(_), false) => b'+',
        }
    }

    fn flush(&mut self, g: &mut Group) -> Result<(), IoError> {
        match (g.mate1.take(), g.mate2.take()) {
            (Some(a), Some(b)) => {
                let (first, second) = if !self.o.mate1_first && self.key(&a) > self.key(&b) {
                    (b, a)
                } else {
                    (a, b)
                };
                let score = if self.o.score_tag.is_some() {
                    first.tag + second.tag
                } else if first.at.is_some() && second.at.is_some() {
                    i64::from(first.mapq.min(second.mapq))
                } else {
                    0
                };
                self.row.clear();
                self.end_columns(&first);
                self.row.push(b'\t');
                self.end_columns(&second);
                self.row.push(b'\t');
                self.row.extend_from_slice(&g.name);
                self.row.push(b'\t');
                self.row.extend_from_slice(score.to_string().as_bytes());
                self.row.push(b'\t');
                self.row.push(Self::strand(&first));
                self.row.push(b'\t');
                self.row.push(Self::strand(&second));
                self.row.push(b'\n');
                self.sum.pairs += 1;
                self.out
                    .write_all(&self.row)
                    .map_err(|source| IoError::Write {
                        path: "output".into(),
                        source,
                    })
            }
            (Some(_), None) | (None, Some(_)) => {
                self.sum.missing_mate += 1;
                Ok(())
            }
            (None, None) => Ok(()),
        }
    }
}

fn score_of(a: &AlnRef<'_>, tag: [u8; 2]) -> Result<i64, String> {
    let t = String::from_utf8_lossy(&tag).into_owned();
    let read = String::from_utf8_lossy(a.name()).into_owned();
    if let Some(v) = a.int_tag(tag) {
        return Ok(v);
    }
    Err(match a.tag_type(tag) {
        None => format!(
            "read {read} has no {t} tag, which --score-tag {t} sums into the score; \
use a tag every mapped read has (NM, AS, ...) or drop --score-tag"
        ),
        Some(ty) => format!(
            "the {t} tag of read {read} has type {}, not an integer; --score-tag needs \
an integer (i) tag such as NM or AS",
            ty as char
        ),
    })
}

fn plural(n: u64, one: &str, many: &str) -> String {
    format!(
        "{} {}",
        gloci_core::thousands(n),
        if n == 1 { one } else { many }
    )
}

/// Write the pairs of the name-grouped alignments `src` as BEDPE rows on `out`, in file
/// order, one per read name that has a primary READ1 and READ2 (D156).
pub fn alignments_to_bedpe(
    mut src: Box<dyn AlnSource>,
    o: &BedpeOptions,
    out: &mut dyn Write,
) -> Result<BedpeSummary, IoError> {
    let name = src.name().to_string();
    let hdr = src.header().clone();
    if hdr.sort_order == SortOrder::Coordinate {
        return Err(UsageError(format!(
            "{name} is sorted by coordinate (@HD SO:coordinate), and --bedpe pairs the reads \
that are next to each other by name; sort it by name first: gloci sort -n {name} -o \
x.byname.bam"
        ))
        .into());
    }
    let layout = src.layout();
    let mut p = Pairing {
        o,
        hdr: &hdr,
        out,
        row: Vec::new(),
        sum: BedpeSummary::default(),
    };
    let mut g = Group {
        name: Vec::new(),
        mate1: None,
        mate2: None,
    };
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if src.next_record(&mut buf)?.is_none() {
            break;
        }
        let fail = |src: &dyn AlnSource, msg: &str| IoError::Detect {
            path: name.clone(),
            msg: format!("{name} {}: {msg}", src.place()),
        };
        let a = aln_ref(layout, &buf, &hdr).map_err(|m| fail(&*src, &m))?;
        let flag = a.flag();
        // Secondary and supplementary alignments are not mates (`#bamtobed-bedpe`).
        if flag & 0x900 != 0 {
            continue;
        }
        let first = flag & 0x40 != 0;
        if first == (flag & 0x80 != 0) {
            p.sum.unpaired += 1;
            continue;
        }
        if a.name() != &g.name[..] {
            p.flush(&mut g)?;
            g.name.clear();
            g.name.extend_from_slice(a.name());
        }
        let at = match (a.ref_id(), flag & 0x4 == 0) {
            (Some(rid), true) => a
                .span(hdr.refs[rid].1)
                .map(|(iv, _)| (rid, iv.start, iv.end)),
            _ => None,
        };
        let tag = match (o.score_tag, at) {
            (Some(t), Some(_)) => score_of(&a, t).map_err(|m| fail(&*src, &m))?,
            _ => 0,
        };
        let end = End {
            at,
            minus: a.strand() == Strand::Minus,
            mapq: a.mapq(),
            tag,
        };
        let slot = if first { &mut g.mate1 } else { &mut g.mate2 };
        slot.get_or_insert(end);
    }
    p.flush(&mut g)?;
    let sum = p.sum;
    src.diagnostics().log_alignment_counts(&name);
    if sum.missing_mate > 0 {
        log::warn!(
            "{name}: {} one mate and no pair was written for them; --bedpe needs the reads of \
a name next to each other (gloci sort -n) and every mate present",
            plural(sum.missing_mate, "read name has", "read names have")
        );
    }
    if sum.unpaired > 0 {
        log::warn!(
            "{name}: {} neither READ1 nor READ2 (flag 0x40 or 0x80) and {} left out; --bedpe \
writes paired-end reads",
            plural(sum.unpaired, "read is", "reads are"),
            if sum.unpaired == 1 { "was" } else { "were" }
        );
    }
    if let Some(s) = src.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    Ok(sum)
}
