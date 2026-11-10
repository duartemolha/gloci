//! `gloci vcf norm`: left-align and trim records against a reference, as
//! `bcftools norm -f` does (without `-a`), plus `--check-ref` and `--rm-dup`.
//! Works on VCF text lines; the reference is read through [`Reference`].

use std::collections::VecDeque;

/// What to do with a record whose REF (or other alleles) disagree with the
/// reference (`bcftools norm -c`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CheckRef {
    /// Fail on the first such record.
    Exit,
    /// Warn and write the record as it is.
    Warn,
    /// Drop the record.
    Exclude,
    /// Set REF from the reference (swapping with a matching ALT).
    Set,
}

/// Which records `--rm-dup` removes (`bcftools norm -d`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RmDup {
    /// A SNP or MNP record at the position of an earlier kept one.
    Snps,
    /// An indel record at the position of an earlier kept indel record.
    Indels,
    /// Both of the above, each within its own type.
    Both,
    /// Any record at the position of an earlier kept one.
    All,
    /// A record with the same REF and ALT alleles (any order, any case).
    Exact,
}

/// Reference bases, read in bounded windows by the caller.
pub trait Reference {
    /// Append the bases of `chrom:[start0, start0+len)` to `out`, clipped at
    /// the end of the contig. Letters other than ACGTN (either case) are
    /// returned as `N`. The error says what is wrong (for example, an
    /// unknown contig).
    fn get(&mut self, chrom: &[u8], start0: u64, len: u64, out: &mut Vec<u8>)
    -> Result<(), String>;
}

/// What happened to a record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// Dropped (`--check-ref x`).
    Drop,
    /// Written at this 1-based position; `changed` is true when the alleles
    /// or position were rewritten.
    Keep { pos: u64, changed: bool },
}

/// Per-record rules of `vcf norm`.
#[derive(Clone, Copy, Debug)]
pub struct Normalizer {
    check: CheckRef,
}

/// Append `chunk` to `out` the way [`Reference::get`] returns bases: ACGTN in
/// either case as they are, every other byte as `N`.
pub fn clean_reference_bases(chunk: &[u8], out: &mut Vec<u8>) {
    out.extend(chunk.iter().map(|&b| match b {
        b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' => b,
        _ => b'N',
    }));
}

fn acgtn(a: &[u8]) -> bool {
    a.iter().all(|&b| {
        matches!(
            b,
            b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
        )
    })
}

fn is_symbolic(a: &[u8]) -> bool {
    a.first() == Some(&b'<')
}

fn is_del(a: &[u8]) -> bool {
    a.starts_with(b"<DEL")
}

fn is_breakend(a: &[u8]) -> bool {
    a.iter().any(|&b| b == b'[' || b == b']')
}

/// An ALT written as bases (letters only).
fn is_bases(a: &[u8]) -> bool {
    !a.is_empty() && a.iter().all(|b| b.is_ascii_alphabetic())
}

struct Cols<'a> {
    chrom: &'a [u8],
    pos: u64,
    id: &'a [u8],
    reference: &'a [u8],
    alt: &'a [u8],
    qual: &'a [u8],
    filter: &'a [u8],
    info: &'a [u8],
    /// FORMAT and the samples, with the tab before them removed.
    rest: Option<&'a [u8]>,
}

fn parse(line: &[u8]) -> Result<Cols<'_>, String> {
    let mut it = line.splitn(9, |&b| b == b'\t');
    let mut c: [&[u8]; 8] = [b""; 8];
    let mut n = 0;
    for slot in c.iter_mut() {
        match it.next() {
            Some(x) => {
                *slot = x;
                n += 1;
            }
            None => break,
        }
    }
    if n < 8 {
        return Err(format!(
            "the record has {n} tab-separated columns, but a VCF record needs at least 8"
        ));
    }
    let pos = std::str::from_utf8(c[1])
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .filter(|&p| p >= 1)
        .ok_or_else(|| {
            format!(
                "POS `{}` is not a positive integer",
                String::from_utf8_lossy(c[1])
            )
        })?;
    Ok(Cols {
        chrom: c[0],
        pos,
        id: c[2],
        reference: c[3],
        alt: c[4],
        qual: c[5],
        filter: c[6],
        info: c[7],
        rest: it.next(),
    })
}

fn show(a: &[u8]) -> String {
    String::from_utf8_lossy(a).into_owned()
}

/// Swap alleles 0 and `i` in every GT of `rest` (FORMAT and samples).
fn swap_gt(rest: &[u8], i: usize, out: &mut Vec<u8>) {
    let mut lines = rest.split(|&b| b == b'\t');
    let Some(format) = lines.next() else { return };
    out.extend_from_slice(format);
    let gt = format.split(|&b| b == b':').position(|k| k == b"GT");
    for s in lines {
        out.push(b'\t');
        let Some(gt) = gt else {
            out.extend_from_slice(s);
            continue;
        };
        for (k, f) in s.split(|&b| b == b':').enumerate() {
            if k > 0 {
                out.push(b':');
            }
            if k != gt {
                out.extend_from_slice(f);
                continue;
            }
            let mut num = Vec::new();
            let flush = |num: &mut Vec<u8>, out: &mut Vec<u8>| {
                if let Some(v) = std::str::from_utf8(num)
                    .ok()
                    .and_then(|s| s.parse::<usize>().ok())
                {
                    let v = if v == 0 {
                        i
                    } else if v == i {
                        0
                    } else {
                        v
                    };
                    out.extend_from_slice(v.to_string().as_bytes());
                } else {
                    out.extend_from_slice(num);
                }
                num.clear();
            };
            for &b in f {
                if b == b'/' || b == b'|' {
                    flush(&mut num, out);
                    out.push(b);
                } else {
                    num.push(b);
                }
            }
            flush(&mut num, out);
        }
    }
}

/// `info` with the value of END replaced by `end` (unchanged without END).
fn set_end(info: &[u8], end: u64, out: &mut Vec<u8>) {
    let mut first = true;
    for kv in info.split(|&b| b == b';') {
        if !first {
            out.push(b';');
        }
        first = false;
        if kv.starts_with(b"END=") {
            out.extend_from_slice(b"END=");
            out.extend_from_slice(end.to_string().as_bytes());
        } else {
            out.extend_from_slice(kv);
        }
    }
}

impl Normalizer {
    pub fn new(check: CheckRef) -> Normalizer {
        Normalizer { check }
    }

    /// Check `line` against the reference and left-align and trim it. The
    /// record to write is left in `out` (the input, or the rewritten line).
    /// An error says what is wrong with the record and how to fix it.
    pub fn process(
        &self,
        line: &[u8],
        refs: &mut dyn Reference,
        out: &mut Vec<u8>,
    ) -> Result<Verdict, String> {
        out.clear();
        let c = parse(line)?;
        let unchanged = |out: &mut Vec<u8>| {
            out.extend_from_slice(line);
            Ok(Verdict::Keep {
                pos: c.pos,
                changed: false,
            })
        };
        let where_ = format!("{}:{}", show(c.chrom), c.pos);
        let mut alleles: Vec<Vec<u8>> = Vec::new();
        alleles.push(c.reference.to_vec());
        if c.alt != b"." || c.alt.is_empty() {
            alleles.extend(c.alt.split(|&b| b == b',').map(|a| a.to_vec()));
        } else {
            alleles.push(b".".to_vec());
        }
        let pos0 = c.pos - 1;
        let reflen = c.reference.len() as u64;
        let mut fasta = Vec::new();
        refs.get(c.chrom, pos0, reflen, &mut fasta)?;

        // Each failed check: error, warn and keep, drop, or (set) fix.
        let mut fixed = false;
        let mut gt_swap: Option<usize> = None;
        let bad_ref = !acgtn(&alleles[0]);
        let mismatch = !alleles[0].eq_ignore_ascii_case(&fasta);
        if bad_ref || mismatch {
            let what = if bad_ref {
                "REF contains a base other than ACGTN"
            } else {
                "REF does not match the reference"
            };
            match self.check {
                CheckRef::Exit => {
                    return Err(format!(
                        "{where_}: {what}: the VCF has `{}`, the FASTA has `{}`; \
fix the VCF, or pass --check-ref w (warn), x (drop the record) or s (set REF from the FASTA)",
                        show(&alleles[0]),
                        show(&fasta)
                    ));
                }
                CheckRef::Warn => {
                    log::warn!(
                        "vcf norm: {where_}: {what}: the VCF has `{}`, the FASTA has `{}`; record written as it is",
                        show(&alleles[0]),
                        show(&fasta)
                    );
                    return unchanged(out);
                }
                CheckRef::Exclude => return Ok(Verdict::Drop),
                CheckRef::Set => {
                    let maxlen = alleles
                        .iter()
                        .filter(|a| is_bases(a))
                        .map(|a| a.len() as u64)
                        .max()
                        .unwrap_or(reflen)
                        .max(reflen);
                    let mut long = Vec::new();
                    refs.get(c.chrom, pos0, maxlen, &mut long)?;
                    let hit = (1..alleles.len()).find(|&i| {
                        let a = &alleles[i];
                        is_bases(a)
                            && a.len() <= long.len()
                            && a.eq_ignore_ascii_case(&long[..a.len()])
                    });
                    match hit {
                        Some(i) => {
                            alleles.swap(0, i);
                            gt_swap = Some(i);
                        }
                        None => {
                            if fasta.is_empty() {
                                return Err(format!(
                                    "{where_}: the position is past the end of the contig in the FASTA, so REF cannot be set; \
check that the VCF and the FASTA use the same assembly"
                                ));
                            }
                            alleles[0] = fasta.clone();
                        }
                    }
                    fixed = true;
                }
            }
        }

        // Other alleles that are not bases.
        let has_bases = alleles[1..].iter().any(|a| is_bases(a));
        let bad_alt = alleles[1..].iter().find(|a| {
            if is_symbolic(a) || a.as_slice() == b"*" {
                return false;
            }
            if is_bases(a) {
                !acgtn(a)
            } else {
                has_bases && (a.as_slice() == b"." || is_breakend(a))
            }
        });
        if let Some(a) = bad_alt {
            let msg = format!(
                "{where_}: the ALT `{}` contains a base other than ACGTN",
                show(a)
            );
            match self.check {
                CheckRef::Exit => {
                    return Err(format!(
                        "{msg}; fix the VCF, or pass --check-ref w (warn), x (drop the record) or s (keep it as it is)"
                    ));
                }
                CheckRef::Warn => {
                    log::warn!("vcf norm: {msg}; record written as it is");
                    return self.finish_fixed(&c, &alleles, fixed, gt_swap, out);
                }
                CheckRef::Exclude => return Ok(Verdict::Drop),
                CheckRef::Set => {
                    return self.finish_fixed(&c, &alleles, fixed, gt_swap, out);
                }
            }
        }

        // An ALT equal to REF.
        if alleles[1..]
            .iter()
            .any(|a| a.eq_ignore_ascii_case(&alleles[0]))
        {
            match self.check {
                CheckRef::Exit => {
                    return Err(format!(
                        "{where_}: an ALT allele is the same as REF (`{}`); fix the VCF, or pass \
--check-ref w, x or s to keep the record as it is",
                        show(&alleles[0])
                    ));
                }
                CheckRef::Warn | CheckRef::Set => {
                    log::warn!(
                        "vcf norm: {where_}: an ALT allele is the same as REF (`{}`); record written without left-alignment",
                        show(&alleles[0])
                    );
                }
                CheckRef::Exclude => {}
            }
            return self.finish_fixed(&c, &alleles, fixed, gt_swap, out);
        }

        // Left-align and trim.
        match realign(c.chrom, pos0, &alleles, refs)? {
            Some((p0, a)) => {
                self.write(&c, p0 + 1, &a, gt_swap, true, out);
                Ok(Verdict::Keep {
                    pos: p0 + 1,
                    changed: true,
                })
            }
            None => self.finish_fixed(&c, &alleles, fixed, gt_swap, out),
        }
    }

    /// Write the record with its (possibly fixed) alleles at its own position.
    fn finish_fixed(
        &self,
        c: &Cols<'_>,
        alleles: &[Vec<u8>],
        fixed: bool,
        gt_swap: Option<usize>,
        out: &mut Vec<u8>,
    ) -> Result<Verdict, String> {
        self.write(c, c.pos, alleles, gt_swap, false, out);
        Ok(Verdict::Keep {
            pos: c.pos,
            changed: fixed,
        })
    }

    fn write(
        &self,
        c: &Cols<'_>,
        pos: u64,
        alleles: &[Vec<u8>],
        gt_swap: Option<usize>,
        update_end: bool,
        out: &mut Vec<u8>,
    ) {
        out.extend_from_slice(c.chrom);
        out.push(b'\t');
        out.extend_from_slice(pos.to_string().as_bytes());
        out.push(b'\t');
        out.extend_from_slice(c.id);
        out.push(b'\t');
        out.extend_from_slice(&alleles[0]);
        out.push(b'\t');
        for (i, a) in alleles[1..].iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            out.extend_from_slice(a);
        }
        out.push(b'\t');
        out.extend_from_slice(c.qual);
        out.push(b'\t');
        out.extend_from_slice(c.filter);
        out.push(b'\t');
        if update_end {
            set_end(c.info, pos + alleles[0].len() as u64 - 1, out);
        } else {
            out.extend_from_slice(c.info);
        }
        if let Some(rest) = c.rest {
            out.push(b'\t');
            match gt_swap {
                Some(i) => swap_gt(rest, i, out),
                None => out.extend_from_slice(rest),
            }
        }
    }
}

/// A new 0-based position and the alleles written there.
type Moved = (u64, Vec<Vec<u8>>);

/// Left-align and trim the alleles of one record (`alleles[0]` is REF), as
/// bcftools does: trim the base every allele ends with, extending all alleles
/// to the left by one reference base whenever one would become empty, then
/// trim the bases every allele starts with. `<DEL...>` ALTs are ignored; a
/// record with another symbolic ALT, or without a base ALT, is left alone.
/// Returns the new 0-based position and alleles, or `None` when nothing
/// changes. Rewritten alleles are upper case, except bases taken from the
/// reference, which keep its case.
fn realign(
    chrom: &[u8],
    pos0: u64,
    alleles: &[Vec<u8>],
    refs: &mut dyn Reference,
) -> Result<Option<Moved>, String> {
    let mut idx = vec![0usize];
    for (i, a) in alleles.iter().enumerate().skip(1) {
        if is_del(a) {
            continue;
        }
        if is_symbolic(a) {
            return Ok(None);
        }
        idx.push(i);
    }
    if idx.len() < 2 {
        return Ok(None);
    }
    let mut w: Vec<VecDeque<u8>> = idx
        .iter()
        .map(|&i| alleles[i].iter().map(|b| b.to_ascii_uppercase()).collect())
        .collect();
    let lens: Vec<usize> = w.iter().map(|a| a.len()).collect();
    let mut pos = pos0;
    let mut base = Vec::new();
    while let Some(&last) = w[0].back() {
        if !w
            .iter()
            .all(|a| a.back().is_some_and(|b| b.eq_ignore_ascii_case(&last)))
        {
            break;
        }
        if w.iter().any(|a| a.len() == 1) {
            if pos == 0 {
                break;
            }
            base.clear();
            refs.get(chrom, pos - 1, 1, &mut base)?;
            let Some(&b) = base.first() else { break };
            for a in w.iter_mut() {
                a.pop_back();
                a.push_front(b);
            }
            pos -= 1;
        } else {
            for a in w.iter_mut() {
                a.pop_back();
            }
        }
    }
    while let Some(&first) = w[0].front() {
        if w.iter().all(|a| a.len() > 1)
            && w.iter()
                .all(|a| a.front().is_some_and(|b| b.eq_ignore_ascii_case(&first)))
        {
            for a in w.iter_mut() {
                a.pop_front();
            }
            pos += 1;
        } else {
            break;
        }
    }
    if pos == pos0 && w.iter().map(|a| a.len()).eq(lens.iter().copied()) {
        return Ok(None);
    }
    let mut out: Vec<Vec<u8>> = alleles.to_vec();
    for (k, &i) in idx.iter().enumerate() {
        out[i] = w[k].iter().copied().collect();
    }
    Ok(Some((pos, out)))
}

/// Applies `--rm-dup` to records in output order.
pub struct Dedup {
    mode: RmDup,
    chrom: Vec<u8>,
    pos: u64,
    /// The kept records at the current position: type bits and the allele set.
    kept: Vec<(u8, Vec<Vec<u8>>)>,
}

const SNP: u8 = 1;
const INDEL: u8 = 2;

impl Dedup {
    pub fn new(mode: RmDup) -> Dedup {
        Dedup {
            mode,
            chrom: Vec::new(),
            pos: 0,
            kept: Vec::new(),
        }
    }

    /// True when `line` is to be written (and is remembered as kept).
    pub fn keep(&mut self, line: &[u8]) -> bool {
        let Ok(c) = parse(line) else { return true };
        if c.chrom != self.chrom.as_slice() || c.pos != self.pos {
            self.chrom.clear();
            self.chrom.extend_from_slice(c.chrom);
            self.pos = c.pos;
            self.kept.clear();
        }
        let reflen = c.reference.len();
        let mut ty = 0u8;
        let mut set: Vec<Vec<u8>> = Vec::with_capacity(4);
        set.push(c.reference.to_ascii_uppercase());
        for a in c.alt.split(|&b| b == b',') {
            set.push(a.to_ascii_uppercase());
            if is_del(a) {
                ty |= INDEL;
            } else if is_bases(a) {
                ty |= if a.len() == reflen { SNP } else { INDEL };
            }
        }
        set[1..].sort();
        let dup = match self.mode {
            RmDup::All => !self.kept.is_empty(),
            RmDup::Exact => self.kept.iter().any(|(_, s)| *s == set),
            RmDup::Snps => self.kept.iter().any(|(t, _)| t & ty & SNP != 0),
            RmDup::Indels => self.kept.iter().any(|(t, _)| t & ty & INDEL != 0),
            RmDup::Both => self.kept.iter().any(|(t, _)| t & ty != 0),
        };
        if dup {
            return false;
        }
        self.kept.push((ty, set));
        true
    }
}
