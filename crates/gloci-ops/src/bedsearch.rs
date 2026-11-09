//! Region queries on uncompressed BED without an index under
//! `--assume-sorted` (D11): find the first line of a chromosome by bisecting
//! byte offsets, as bedops `bedextract` does, then scan that chromosome, so
//! features that start before a region and reach into it are found. Only
//! the blocks read are checked for order; the file is trusted to be sorted.

use gloci_core::{ChromTable, RegionSet};
use gloci_io::{
    ChromMap, Diagnostics, Format, IoError, Loc, OrderGuard, ReaderOptions, RecordBatch,
    RecordMeta, RecordReader, is_bed_record_line, open_reader_from,
};
use std::cmp::Ordering;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// Where bisection found a chromosome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChromBlock {
    /// The byte offset of the chromosome's first line.
    At(u64),
    /// The search ended on another chromosome, at the end of the file, or on
    /// a name the table cannot place: the file lacks the chromosome, or is
    /// not sorted in the table's order.
    NotFoundByBisection,
}

/// The first BED record line at or after a line start.
struct Probe {
    start: u64,
    chrom: Vec<u8>,
}

/// The first record line starting at or after the first line start at or
/// after `x` (offset 0, or just past a `\n`). Blank, `#`, `track` and
/// `browser` lines are skipped, as the BED reader skips them.
fn probe(file: &mut (impl Read + Seek), x: u64, len: u64) -> io::Result<Option<Probe>> {
    let from = x.saturating_sub(1);
    file.seek(SeekFrom::Start(from))?;
    let mut r = BufReader::new(&mut *file);
    let mut pos = from;
    if x > 0 {
        // the rest of the line holding byte x-1 (just the `\n` if x starts a line)
        let n = r.skip_until(b'\n')?;
        if n == 0 {
            return Ok(None);
        }
        pos += n as u64;
    }
    let mut line = Vec::new();
    while pos < len {
        line.clear();
        let n = r.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        let start = pos;
        pos += n as u64;
        let body = line.strip_suffix(b"\n").unwrap_or(&line);
        let body = body.strip_suffix(b"\r").unwrap_or(body);
        if is_bed_record_line(body) {
            let end = body.iter().position(|&b| b == b'\t').unwrap_or(body.len());
            return Ok(Some(Probe {
                start,
                chrom: body[..end].to_vec(),
            }));
        }
    }
    Ok(None)
}

/// Find the first line of `chrom` in the BED text `file` of `len` bytes.
///
/// Bisects byte offsets for the first record line whose chromosome sorts at
/// or after `chrom` in `table`'s order, and accepts it only if that line is
/// on `chrom`. The line before it is then on a chromosome that sorts before
/// `chrom`, so the line starts `chrom`'s block. A file sorted in another
/// chromosome order may still be found (any block start is a block start)
/// or not (`NotFoundByBisection`); the caller then scans the whole file.
///
/// Read errors name the placeholder path `input`
/// (see [`IoError::with_input_name`]).
pub fn find_chrom_block(
    file: &mut (impl Read + Seek),
    len: u64,
    chrom: &[u8],
    table: &ChromTable,
) -> Result<ChromBlock, IoError> {
    find_chrom_block_mapped(file, len, chrom, table, None)
}

/// [`find_chrom_block`] for a file read with `--chrom-map`: `chrom` and
/// `table` use the renamed names, so each line's chromosome is renamed
/// before it is compared. A chromosome the map renames away is not found.
pub fn find_chrom_block_mapped(
    file: &mut (impl Read + Seek),
    len: u64,
    chrom: &[u8],
    table: &ChromTable,
    map: Option<&ChromMap>,
) -> Result<ChromBlock, IoError> {
    let named =
        |p: &Probe| -> Vec<u8> { map.map_or(p.chrom.clone(), |m| m.rename(&p.chrom).to_vec()) };
    let read_err = |source| IoError::Read {
        path: "input".to_string(),
        source,
    };
    // Lower bound over offsets x in [0, len] of "the first record at or
    // after x's line is on `chrom` or a later chromosome, or there is none".
    // Invariant: that is false at lo - 1 (when lo > 0) and true at hi.
    let (mut lo, mut hi) = (0u64, len);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        match probe(file, mid, len).map_err(read_err)? {
            None => hi = mid,
            Some(p) => match table.cmp_names(&named(&p), chrom) {
                None => return Ok(ChromBlock::NotFoundByBisection),
                // every offset up to p.start finds p as well
                Some(Ordering::Less) => lo = p.start + 1,
                Some(_) => hi = mid,
            },
        }
    }
    Ok(match probe(file, lo, len).map_err(read_err)? {
        Some(p) if named(&p) == chrom => ChromBlock::At(p.start),
        _ => ChromBlock::NotFoundByBisection,
    })
}

/// A chromosome's block and the regions wanted on it.
pub(crate) struct Block {
    pub(crate) offset: u64,
    pub(crate) chrom: Vec<u8>,
    pub(crate) set: RegionSet,
}

/// What locating every region chromosome found.
pub(crate) enum Located {
    /// One block per chromosome, in file order.
    Blocks(Vec<Block>),
    /// The chromosomes bisection could not find.
    Missing(Vec<String>),
}

/// Bisect for the block of every chromosome of `sets`.
pub(crate) fn locate_blocks(
    path: &Path,
    name: &str,
    sets: Vec<(Vec<u8>, RegionSet)>,
    table: &ChromTable,
    map: Option<&ChromMap>,
) -> Result<Located, IoError> {
    let mut f = File::open(path).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let len = f
        .metadata()
        .map_err(|source| IoError::Read {
            path: name.to_string(),
            source,
        })?
        .len();
    let (mut blocks, mut missing) = (Vec::new(), Vec::new());
    for (chrom, set) in sets {
        match find_chrom_block_mapped(&mut f, len, &chrom, table, map)
            .map_err(|e| e.with_input_name(name))?
        {
            ChromBlock::At(offset) => blocks.push(Block { offset, chrom, set }),
            ChromBlock::NotFoundByBisection => {
                missing.push(String::from_utf8_lossy(&chrom).into_owned());
            }
        }
    }
    if !missing.is_empty() {
        return Ok(Located::Missing(missing));
    }
    blocks.sort_by_key(|b| b.offset);
    Ok(Located::Blocks(blocks))
}

/// Number of `\n` bytes before `offset`: the lines a reader that starts at
/// `offset` does not count.
fn lines_before(path: &Path, name: &str, offset: u64) -> Result<u64, IoError> {
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let f = File::open(path).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let mut r = f.take(offset);
    let mut buf = vec![0u8; 1 << 18];
    let mut n = 0u64;
    loop {
        let got = match r.read(&mut buf) {
            Ok(0) => return Ok(n),
            Ok(k) => k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(read_err(e)),
        };
        n += buf[..got].iter().filter(|&&b| b == b'\n').count() as u64;
    }
}

/// Renumber the lines an error names from a reader that started at
/// `offset` to lines of the whole file. Only done on error, so the search
/// never pays for counting lines.
fn with_file_lines(e: IoError, path: &Path, name: &str, offset: u64) -> IoError {
    let shift = |by: u64, loc: Loc| match loc {
        Loc::Line(n) => Loc::Line(n + by),
        other => other,
    };
    match e {
        IoError::Parse {
            path: p,
            loc: Loc::Line(n),
            msg,
        } => match lines_before(path, name, offset) {
            Ok(by) => IoError::Parse {
                path: p,
                loc: Loc::Line(n + by),
                msg,
            },
            Err(count_err) => count_err,
        },
        IoError::Order(mut o) => match lines_before(path, name, offset) {
            Ok(by) => {
                o.loc = shift(by, o.loc);
                o.prev_loc = shift(by, o.prev_loc);
                IoError::Order(o)
            }
            Err(count_err) => count_err,
        },
        e => e,
    }
}

/// The block being read.
struct Current {
    block: usize,
    reader: Box<dyn RecordReader>,
    guard: OrderGuard,
}

/// Reads the blocks found by bisection, in file order. Each block is read
/// from its first line until the first record on another chromosome or past
/// the block's last region end, with that chromosome's order checked.
pub(crate) struct BlockScan {
    path: PathBuf,
    name: String,
    format: Format,
    opts: ReaderOptions,
    blocks: Vec<Block>,
    next: usize,
    cur: Option<Current>,
    /// Invalid records skipped so far, numbered by line of the whole file.
    pub(crate) diag: Diagnostics,
}

impl BlockScan {
    pub(crate) fn new(
        path: &Path,
        name: &str,
        format: Format,
        ro: &ReaderOptions,
        blocks: Vec<Block>,
    ) -> BlockScan {
        BlockScan {
            path: path.to_path_buf(),
            name: name.to_string(),
            format,
            opts: ReaderOptions {
                threads: NonZeroUsize::MIN,
                format: Some(format),
                ..ro.clone()
            },
            blocks,
            next: 0,
            cur: None,
            diag: Diagnostics::default(),
        }
    }

    fn open_next(&mut self) -> Result<bool, IoError> {
        let Some(b) = self.blocks.get(self.next) else {
            return Ok(false);
        };
        let mut f = File::open(&self.path).map_err(|source| IoError::Open {
            path: self.name.clone(),
            source,
        })?;
        f.seek(SeekFrom::Start(b.offset))
            .map_err(|source| IoError::Read {
                path: self.name.clone(),
                source,
            })?;
        let reader = open_reader_from(&self.name, Box::new(f), &self.opts)?;
        self.cur = Some(Current {
            block: self.next,
            reader,
            guard: OrderGuard::contiguous(&self.name, self.format),
        });
        self.next += 1;
        Ok(true)
    }

    /// Append to `out` the next records that overlap the current block's
    /// regions (possibly none); false once every block has been read.
    pub(crate) fn read(
        &mut self,
        chroms: &mut ChromTable,
        scratch: &mut RecordBatch,
        out: &mut RecordBatch,
        max: usize,
        inclusive: bool,
    ) -> Result<bool, IoError> {
        if self.cur.is_none() && !self.open_next()? {
            return Ok(false);
        }
        let cur = self.cur.as_mut().expect("a block is open");
        let block = &self.blocks[cur.block];
        let fix = |e| with_file_lines(e, &self.path, &self.name, block.offset);
        scratch.clear();
        let more = cur.reader.read_batch(chroms, scratch, max).map_err(fix)?;
        let id = chroms.get(&block.chrom);
        let (mut end, mut done) = (scratch.len(), !more);
        for (i, m) in scratch.recs.iter().enumerate() {
            if Some(m.chrom) != id {
                end = i;
                done = true;
                break;
            }
            if m.key > block.set.end() {
                // checked for order, but nothing past here can match
                end = i + 1;
                done = true;
                break;
            }
        }
        scratch.recs.truncate(end);
        cur.guard.check(chroms, scratch).map_err(fix)?;
        append_matching(out, scratch, |m| block.set.matches(m.span, inclusive));
        if done {
            let cur = self.cur.take().expect("a block is open");
            let mut d = cur.reader.diagnostics().clone();
            if d.skipped > 0 {
                d.shift_lines(lines_before(&self.path, &self.name, block.offset)?);
            }
            self.diag.absorb(d);
        }
        Ok(true)
    }
}

/// Append the records of `from` for which `keep` holds to `out`.
pub(crate) fn append_matching(
    out: &mut RecordBatch,
    from: &RecordBatch,
    mut keep: impl FnMut(&RecordMeta) -> bool,
) {
    for (i, m) in from.recs.iter().enumerate() {
        if keep(m) {
            let off = out.buf.len();
            out.buf.extend_from_slice(from.raw(i));
            out.recs.push(RecordMeta { off, ..*m });
        }
    }
}
