//! Chain files for `gloci liftover` (D174, R1, R2): read as CrossMap reads them, into aligned
//! blocks indexed per source chromosome.
//!
//! A chain is a header (`chain score tName tSize + tStart tEnd qName qSize qStrand qStart qEnd
//! [id]`) and block lines (`size dt dq`, then a last `size`). The source is the `t` side and
//! must be `+`; the target is the `q` side, whose coordinates count from the end of its
//! chromosome for a `-` chain. A [`Segment`] always carries forward target coordinates.

use crate::error::{IoError, Loc};
use gloci_core::Interval;
use rust_lapper::{Interval as LapInterval, Lapper};
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// What a chain's header says, kept for the chain's place in the file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainInfo {
    /// The header's line number.
    pub line: u64,
    pub source: String,
    pub target: String,
    pub target_size: u64,
    /// The target strand is `-`.
    pub reverse: bool,
}

/// One piece of a record: the part of an aligned block that overlaps it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment {
    /// The chain's place in the file, from 0 ([`ChainIndex::chain`]).
    pub chain: u32,
    pub source: Interval,
    /// In forward target coordinates, also for a `-` chain.
    pub target: Interval,
    pub reverse: bool,
}

/// A source chromosome of a [`ChainIndex`].
#[derive(Debug, Clone)]
pub struct SourceRef(usize);

#[derive(Debug, Clone, PartialEq, Eq)]
struct Block {
    /// Where the block starts in the chain's own target coordinates.
    target_start: u64,
    chain: u32,
}

/// How liftover treats the mitochondrial chromosome for a chain (D174, R13).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MtRule {
    /// Names as they are in the chain: `chrM` finds no `M` and is never sent to `MT` (R4), and
    /// a target `MT` is styled like any other name (`chrMT` for a `chr` input under `as-input`).
    #[default]
    Plain,
    /// The built-in GRCh37-to-GRCh38 chain: hg19's `chrM` is not GRCh37's `MT`, so `chrM` and
    /// `M` stay unmapped.
    Hg19,
}

/// The aligned blocks of a chain file, one interval index per source chromosome.
pub struct ChainIndex {
    pub(crate) mt: MtRule,
    names: HashMap<String, usize>,
    lappers: Vec<Lapper<u64, Block>>,
    chains: Vec<ChainInfo>,
    blocks: usize,
}

/// The chain being read.
struct Open {
    line: u64,
    t_end: u64,
    q_end: u64,
    t: u64,
    q: u64,
    chrom: usize,
    chain: u32,
}

/// A chromosome size and the line that gave it.
type Sizes = HashMap<String, (u64, u64)>;

impl ChainIndex {
    /// Read the chain file at `path` (`-` is stdin, and an in-memory source is read once, D84):
    /// plain, gzip or BGZF, told by magic bytes. Errors call it by its display name.
    pub fn load(path: &str) -> Result<ChainIndex, IoError> {
        let raw = crate::input::open_raw(path)?;
        Self::from_reader(&crate::memio::display_name(path), BufReader::new(raw))
    }

    /// Read chain text held in memory (plain, gzip or BGZF); `name` is what errors call it.
    pub fn parse(name: &str, text: &[u8]) -> Result<ChainIndex, IoError> {
        Self::from_reader(name, text)
    }

    fn from_reader<R: BufRead>(name: &str, mut raw: R) -> Result<ChainIndex, IoError> {
        let read_err = |source| IoError::Read {
            path: name.to_string(),
            source,
        };
        let gz = raw.fill_buf().map_err(read_err)?.starts_with(&[0x1f, 0x8b]);
        if gz {
            build(
                name,
                BufReader::new(flate2::bufread::MultiGzDecoder::new(raw)),
            )
        } else {
            build(name, raw)
        }
    }

    /// The source chromosome called exactly `chrom`.
    pub fn source(&self, chrom: &str) -> Option<SourceRef> {
        self.names.get(chrom).map(|&i| SourceRef(i))
    }

    /// The pieces of `iv`, in `out` (cleared first): the blocks with `block.start < iv.end &&
    /// block.end > iv.start`, clipped to `iv`, ordered by source start and then the chain's
    /// place in the file (R2).
    pub fn segments(&self, src: &SourceRef, iv: Interval, out: &mut Vec<Segment>) {
        out.clear();
        for b in self.lappers[src.0].find(iv.start, iv.end) {
            let (start, end) = (b.start.max(iv.start), b.stop.min(iv.end));
            let info = &self.chains[b.val.chain as usize];
            let q0 = b.val.target_start + (start - b.start);
            let q1 = b.val.target_start + (end - b.start);
            let target = if info.reverse {
                Interval {
                    start: info.target_size - q1,
                    end: info.target_size - q0,
                }
            } else {
                Interval { start: q0, end: q1 }
            };
            out.push(Segment {
                chain: b.val.chain,
                source: Interval { start, end },
                target,
                reverse: info.reverse,
            });
        }
        out.sort_by_key(|s| (s.source.start, s.chain));
    }

    /// The header of chain `i` (a [`Segment::chain`]).
    pub fn chain(&self, i: u32) -> &ChainInfo {
        &self.chains[i as usize]
    }

    pub fn chains(&self) -> usize {
        self.chains.len()
    }

    pub fn blocks(&self) -> usize {
        self.blocks
    }

    /// How the mitochondrial chromosome is treated; only a built-in chain changes it.
    pub fn mt_rule(&self) -> MtRule {
        self.mt
    }
}

fn build<R: BufRead>(name: &str, mut r: R) -> Result<ChainIndex, IoError> {
    let err = |line: u64, msg: String| IoError::Parse {
        path: name.to_string(),
        loc: Loc::Line(line),
        msg,
    };
    let mut names: HashMap<String, usize> = HashMap::new();
    let mut lists: Vec<Vec<LapInterval<u64, Block>>> = vec![];
    let mut chains: Vec<ChainInfo> = vec![];
    let (mut t_sizes, mut q_sizes): (Sizes, Sizes) = (HashMap::new(), HashMap::new());
    let mut blocks = 0usize;
    let mut open: Option<Open> = None;
    let (mut buf, mut n, mut last) = (Vec::new(), 0u64, 0u64);
    loop {
        buf.clear();
        let got = r
            .read_until(b'\n', &mut buf)
            .map_err(|source| IoError::Read {
                path: name.to_string(),
                source,
            })?;
        if got == 0 {
            break;
        }
        n += 1;
        let f: Vec<&[u8]> = buf
            .split(|b| b.is_ascii_whitespace())
            .filter(|t| !t.is_empty())
            .collect();
        if f.is_empty() || f[0][0] == b'#' {
            continue;
        }
        last = n;
        if f[0] == b"chain" {
            if let Some(o) = &open {
                return Err(err(
                    n,
                    format!(
                        "the chain at line {} ends only at its line of one number, but a new chain header came first",
                        o.line
                    ),
                ));
            }
            if f.len() != 12 && f.len() != 13 {
                return Err(err(
                    n,
                    format!("a chain header has 12 or 13 fields, found {}", f.len()),
                ));
            }
            let text = |i: usize| {
                std::str::from_utf8(f[i])
                    .map_err(|_| err(n, "a chain header is not UTF-8 text".into()))
            };
            let num = |i: usize, what: &str| {
                std::str::from_utf8(f[i])
                    .ok()
                    .and_then(|s| s.parse::<u64>().ok())
                    .ok_or_else(|| {
                        err(
                            n,
                            format!("{what} is not a number: {}", String::from_utf8_lossy(f[i])),
                        )
                    })
            };
            let (t_name, q_name) = (text(2)?, text(7)?);
            let (t_size, t_start, t_end) = (num(3, "tSize")?, num(5, "tStart")?, num(6, "tEnd")?);
            let (q_size, q_start, q_end) = (num(8, "qSize")?, num(10, "qStart")?, num(11, "qEnd")?);
            if f[4] != b"+" {
                return Err(err(
                    n,
                    format!(
                        "the source strand is {}; only `+` is supported",
                        String::from_utf8_lossy(f[4])
                    ),
                ));
            }
            let reverse = match f[9] {
                b"+" => false,
                b"-" => true,
                s => {
                    return Err(err(
                        n,
                        format!(
                            "the target strand is {}; it must be `+` or `-`",
                            String::from_utf8_lossy(s)
                        ),
                    ));
                }
            };
            for (side, start, end, size) in [
                ("tStart/tEnd", t_start, t_end, t_size),
                ("qStart/qEnd", q_start, q_end, q_size),
            ] {
                if start > end || end > size {
                    return Err(err(
                        n,
                        format!("{side} {start}..{end} must lie inside the size {size}"),
                    ));
                }
            }
            for (sizes, chrom, size) in [
                (&mut t_sizes, t_name, t_size),
                (&mut q_sizes, q_name, q_size),
            ] {
                let seen = *sizes.entry(chrom.to_string()).or_insert((size, n));
                if seen.0 != size {
                    return Err(err(
                        n,
                        format!(
                            "{chrom} has size {size} here but size {} on line {}",
                            seen.0, seen.1
                        ),
                    ));
                }
            }
            let chrom = *names.entry(t_name.to_string()).or_insert_with(|| {
                lists.push(vec![]);
                lists.len() - 1
            });
            chains.push(ChainInfo {
                line: n,
                source: t_name.to_string(),
                target: q_name.to_string(),
                target_size: q_size,
                reverse,
            });
            open = Some(Open {
                line: n,
                t_end,
                q_end,
                t: t_start,
                q: q_start,
                chrom,
                chain: (chains.len() - 1) as u32,
            });
            continue;
        }
        let Some(o) = open.as_mut() else {
            return Err(err(n, "a block line before any chain header".into()));
        };
        if f.len() != 1 && f.len() != 3 {
            return Err(err(
                n,
                format!("a block line has 1 or 3 numbers, found {}", f.len()),
            ));
        }
        let mut v = [0u64; 3];
        for (i, tok) in f.iter().enumerate() {
            v[i] = std::str::from_utf8(tok)
                .ok()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| err(n, format!("not a number: {}", String::from_utf8_lossy(tok))))?;
        }
        let size = v[0];
        let end_of = |a: u64, b: u64| {
            a.checked_add(b)
                .ok_or_else(|| err(n, "block sizes overflow".into()))
        };
        if size > 0 {
            lists[o.chrom].push(LapInterval {
                start: o.t,
                stop: end_of(o.t, size)?,
                val: Block {
                    target_start: o.q,
                    chain: o.chain,
                },
            });
            blocks += 1;
        }
        o.t = end_of(o.t, size)?;
        o.q = end_of(o.q, size)?;
        if f.len() == 3 {
            o.t = end_of(o.t, v[1])?;
            o.q = end_of(o.q, v[2])?;
            continue;
        }
        if (o.t, o.q) != (o.t_end, o.q_end) {
            return Err(err(
                n,
                format!(
                    "the blocks of the chain at line {} end at {} (source) and {} (target), but its header says tEnd {} and qEnd {}",
                    o.line, o.t, o.q, o.t_end, o.q_end
                ),
            ));
        }
        open = None;
    }
    if let Some(o) = open {
        return Err(err(
            last,
            format!("the file ends inside the chain at line {}", o.line),
        ));
    }
    if chains.is_empty() {
        return Err(IoError::Detect {
            path: name.to_string(),
            msg: format!("{name}: no chain in the file"),
        });
    }
    let lappers = lists.into_iter().map(Lapper::new).collect();
    Ok(ChainIndex {
        mt: MtRule::Plain,
        names,
        lappers,
        chains,
        blocks,
    })
}
