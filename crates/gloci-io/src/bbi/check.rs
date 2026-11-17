//! Reading back a bigWig or bigBed that bigtools has just written, before the rename (D171).
//!
//! bigtools spools each chromosome's sections and zoom data to temporary files, and loses
//! some of their write errors (a buffer flushed when it is dropped): when it cannot create or
//! fill one, the file comes out with sections missing while bigtools reports success. So the
//! temporary file is checked against the format (the bigWig/bigBed format paper,
//! doi:10.1093/bioinformatics/btq351) before it is renamed into place:
//! - the header, the zoom headers and the chromosome tree parse;
//! - the blocks of the full-data index, and of each zoom level's index, lie in their data
//!   region in order and fill it exactly;
//! - every block decompresses and agrees with its index entry;
//! - the full-data blocks hold as many records as were written.
//!
//! The same pass finds the index ends to correct. bigtools 0.5.8 gives a bigBed section, and
//! so its index entry and the nodes above it, the end of its last record rather than its
//! largest end, so a region read after that end skips a record that encloses the ones after
//! it (a gene around shorter genes). The pass computes each block's largest end, and each
//! node's as the largest of its children's, and returns the index entries whose end differs
//! as [`EndPatch`]es, which the writer applies before the rename.

use super::write::BbiKind;
use crate::detect::{BIGBED_MAGIC, BIGWIG_MAGIC};
use flate2::{Decompress, FlushDecompress, Status};
use gloci_core::thousands;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;
const CIR_TREE_MAGIC: u32 = 0x2468_ACE0;
/// The fixed header, before the zoom headers.
const HEADER: usize = 64;
const ZOOM_HEADER: usize = 24;
const CHROM_TREE_HEADER: usize = 32;
/// An R-tree's own header, before its root node.
const CIR_TREE_HEADER: usize = 48;
const CIR_LEAF_ITEM: usize = 32;
const CIR_NODE_ITEM: usize = 24;
const WIG_SECTION_HEADER: usize = 24;
const ZOOM_RECORD: usize = 32;

/// What the check went through, and the index ends to correct.
#[derive(Debug)]
pub(crate) struct Checked {
    pub sections: u64,
    pub zoom_levels: usize,
    pub patches: Vec<EndPatch>,
}

/// An index end to rewrite: endChromIx then endBase, as two 32-bit numbers at `at`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct EndPatch {
    pub at: u64,
    pub chrom: u32,
    pub end: u32,
}

/// A position as an index orders ends: chromosome id, then base.
type End = (u32, u32);

/// Checks the bigWig or bigBed at `path` (see the module docs), which must hold `records`
/// records. The error says what is wrong with it.
pub(crate) fn check_written(path: &Path, kind: BbiKind, records: u64) -> Result<Checked, String> {
    let file = File::open(path).map_err(|e| format!("it cannot be reopened to check it: {e}"))?;
    let len = file
        .metadata()
        .map_err(|e| format!("its size cannot be read: {e}"))?
        .len();
    let mut r = Reader {
        file,
        pos: 0,
        len,
        little: true,
        buf: Vec::new(),
    };
    let head = r.read(0, HEADER, "the header")?.to_vec();
    let magic = match kind {
        BbiKind::BigWig => BIGWIG_MAGIC,
        BbiKind::BigBed => BIGBED_MAGIC,
    };
    let first = [head[0], head[1], head[2], head[3]];
    let little = u32::from_le_bytes(first) == magic;
    if !little && u32::from_be_bytes(first) != magic {
        return Err("it does not start with its magic number".into());
    }
    r.little = little;
    let mut f = Fields::new(&head, little);
    f.skip(6);
    let zoom_levels = usize::from(f.u16());
    let chrom_tree = f.u64();
    let full_data = f.u64();
    let full_index = f.u64();
    f.skip(4);
    let autosql = f.u64();
    let total_summary = f.u64();
    let uncompress_buf = f.u32();
    let zh = r.read(HEADER as u64, ZOOM_HEADER * zoom_levels, "the zoom headers")?;
    let zooms: Vec<(u64, u64)> = zh
        .as_chunks::<ZOOM_HEADER>()
        .0
        .iter()
        .map(|z| {
            let mut f = Fields::new(z, little);
            f.skip(8);
            (f.u64(), f.u64())
        })
        .collect();
    // Where each part of the file starts: a data region ends where the next part starts.
    let mut parts = vec![chrom_tree, full_data, full_index, len.saturating_sub(4)];
    parts.extend([autosql, total_summary].into_iter().filter(|&o| o != 0));
    parts.extend(zooms.iter().flat_map(|&(d, i)| [d, i]));
    if let Some(&o) = parts.iter().find(|&&o| o < HEADER as u64 || o >= len) {
        return Err(format!(
            "its header points at byte {}, outside the file ({} bytes)",
            thousands(o),
            thousands(len)
        ));
    }
    let chroms = chrom_lengths(&mut r, chrom_tree)?;
    let data_count = {
        let b = r.read(full_data, 8, "the data count")?;
        Fields::new(b, little).u64()
    };
    let mut inflate = Inflater::new(uncompress_buf);
    let mut items = 0u64;
    let mut patches = Vec::new();
    let start = full_data + 8;
    let sections = walk_index(
        &mut r,
        full_index,
        "the full-data index",
        (start, region_end(start, &parts)),
        &mut patches,
        |r, leaf| {
            leaf.check_chrom(&chroms)?;
            let block = inflate.load(r, leaf, "a full-data block")?;
            let (n, end) = block_records(kind, block, leaf, little)?;
            items += n;
            Ok(end)
        },
    )?;
    // gloci refuses overlapping bigWig records, so a bigWig's ends ascend and every index end
    // is already its largest: a correction would mean a damaged index.
    if kind == BbiKind::BigWig && !patches.is_empty() {
        return Err(format!(
            "its full-data index ends at byte {} do not match its sections",
            thousands(patches[0].at)
        ));
    }
    if items != records {
        return Err(format!(
            "its blocks hold {} records, but {} were written",
            thousands(items),
            thousands(records)
        ));
    }
    let (count, of) = match kind {
        BbiKind::BigWig => (sections, "sections"),
        BbiKind::BigBed => (records, "records"),
    };
    if data_count != count {
        return Err(format!(
            "its data count is {}, but it has {} {of}",
            thousands(data_count),
            thousands(count)
        ));
    }
    // bigtools builds zoom records as tiles of ascending, non-overlapping spans, so their
    // index ends are already the largest; the same walk would correct them if they were not.
    for (level, &(data, index)) in zooms.iter().enumerate() {
        let block = format!("a block of zoom level {}", level + 1);
        walk_index(
            &mut r,
            index,
            &format!("the index of zoom level {}", level + 1),
            (data, region_end(data, &parts)),
            &mut patches,
            |r, leaf| {
                leaf.check_chrom(&chroms)?;
                let b = inflate.load(r, leaf, &block)?;
                check_zoom_block(b, leaf, little)
            },
        )?;
    }
    for p in &patches {
        if chroms.get(&p.chrom).is_none_or(|&len| p.end > len) {
            return Err(format!(
                "a block ends at {} on chromosome id {}, past its chromosome",
                thousands(u64::from(p.end)),
                p.chrom
            ));
        }
    }
    Ok(Checked {
        sections,
        zoom_levels: zooms.len(),
        patches,
    })
}

/// Where the data region starting at `start` ends: where the next part of the file starts.
fn region_end(start: u64, parts: &[u64]) -> u64 {
    parts
        .iter()
        .copied()
        .filter(|&p| p > start)
        .min()
        .unwrap_or(start)
}

/// The file, read at offsets, in its byte order.
struct Reader {
    file: File,
    /// Where the next read starts without a seek; reads in file order need none.
    pos: u64,
    len: u64,
    little: bool,
    buf: Vec<u8>,
}

impl Reader {
    /// The `n` bytes at `off`, which must lie inside the file.
    fn read(&mut self, off: u64, n: usize, what: &str) -> Result<&[u8], String> {
        let end = u64::try_from(n)
            .ok()
            .and_then(|n| off.checked_add(n))
            .filter(|&e| e <= self.len)
            .ok_or_else(|| {
                format!(
                    "{what} at byte {} runs past the end of the file ({} bytes)",
                    thousands(off),
                    thousands(self.len)
                )
            })?;
        if self.pos != off {
            self.file
                .seek(SeekFrom::Start(off))
                .map_err(|e| format!("{what} cannot be read: {e}"))?;
        }
        self.buf.resize(n, 0);
        self.pos = u64::MAX;
        self.file
            .read_exact(&mut self.buf)
            .map_err(|e| format!("{what} cannot be read: {e}"))?;
        self.pos = end;
        Ok(&self.buf)
    }
}

/// Numbers read in order from a slice that holds them all.
struct Fields<'a> {
    b: &'a [u8],
    at: usize,
    little: bool,
}

impl<'a> Fields<'a> {
    fn new(b: &'a [u8], little: bool) -> Fields<'a> {
        Fields { b, at: 0, little }
    }

    fn skip(&mut self, n: usize) {
        self.at += n;
    }

    fn take<const N: usize>(&mut self) -> [u8; N] {
        let mut a = [0u8; N];
        a.copy_from_slice(&self.b[self.at..self.at + N]);
        self.at += N;
        a
    }

    fn u8(&mut self) -> u8 {
        self.take::<1>()[0]
    }

    fn u16(&mut self) -> u16 {
        let a = self.take();
        if self.little {
            u16::from_le_bytes(a)
        } else {
            u16::from_be_bytes(a)
        }
    }

    fn u32(&mut self) -> u32 {
        let a = self.take();
        if self.little {
            u32::from_le_bytes(a)
        } else {
            u32::from_be_bytes(a)
        }
    }

    fn u64(&mut self) -> u64 {
        let a = self.take();
        if self.little {
            u64::from_le_bytes(a)
        } else {
            u64::from_be_bytes(a)
        }
    }
}

/// The node count past which a tree is taken to loop: no file has more nodes than 4-byte words.
fn too_many_nodes(nodes: u64, len: u64) -> bool {
    nodes > len / 4
}

/// The chromosome lengths by id, from the B+ tree at `at`.
fn chrom_lengths(r: &mut Reader, at: u64) -> Result<HashMap<u32, u32>, String> {
    let what = "the chromosome tree";
    let little = r.little;
    let h = r.read(at, CHROM_TREE_HEADER, what)?;
    let mut f = Fields::new(h, little);
    if f.u32() != CHROM_TREE_MAGIC {
        return Err(format!("{what} does not start with its magic number"));
    }
    f.skip(4);
    let key = f.u32() as usize;
    let val = f.u32();
    let count = f.u64();
    if val != 8 {
        return Err(format!(
            "{what} gives {val}-byte values; a chromosome's id and length take 8"
        ));
    }
    let item = key
        .checked_add(8)
        .ok_or_else(|| format!("{what} has a key size too large"))?;
    let mut lengths = HashMap::new();
    let mut stack = vec![at + CHROM_TREE_HEADER as u64];
    let mut nodes = 0u64;
    while let Some(node) = stack.pop() {
        nodes += 1;
        if too_many_nodes(nodes, r.len) {
            return Err(format!("{what} loops"));
        }
        let h = r.read(node, 4, what)?;
        let is_leaf = h[0] == 1;
        let n = usize::from(Fields::new(&h[2..], little).u16());
        let size = n
            .checked_mul(item)
            .ok_or_else(|| format!("{what} has a node too large"))?;
        let body = r.read(node + 4, size, what)?;
        let mut children = Vec::new();
        for e in body.chunks_exact(item) {
            let mut f = Fields::new(&e[key..], little);
            if is_leaf {
                let (id, length) = (f.u32(), f.u32());
                if lengths.insert(id, length).is_some() {
                    return Err(format!("{what} lists chromosome id {id} twice"));
                }
            } else {
                children.push(f.u64());
            }
        }
        stack.extend(children.into_iter().rev());
    }
    if lengths.len() as u64 != count {
        return Err(format!(
            "{what} lists {} chromosomes, but its header says {}",
            thousands(lengths.len() as u64),
            thousands(count)
        ));
    }
    Ok(lengths)
}

/// One block of an R-tree index.
#[derive(Clone, Copy, Debug)]
struct Leaf {
    chrom: u32,
    start: u32,
    end_chrom: u32,
    end: u32,
    offset: u64,
    size: u64,
}

impl Leaf {
    fn describe(&self) -> String {
        format!(
            "at byte {} ({} bytes, chromosome id {} from {})",
            thousands(self.offset),
            thousands(self.size),
            self.chrom,
            thousands(u64::from(self.start))
        )
    }

    /// Its chromosomes are in the chromosome tree, and it ends within the last one.
    fn check_chrom(&self, chroms: &HashMap<u32, u32>) -> Result<(), String> {
        let unknown =
            |id| format!("a block names chromosome id {id}, which the chromosome tree lacks");
        chroms.get(&self.chrom).ok_or_else(|| unknown(self.chrom))?;
        let length = *chroms
            .get(&self.end_chrom)
            .ok_or_else(|| unknown(self.end_chrom))?;
        if self.end_chrom < self.chrom || self.end > length {
            return Err(format!(
                "the block {} does not fit its chromosome",
                self.describe()
            ));
        }
        Ok(())
    }
}

/// Where an R-tree walk stands: the next block's offset, the blocks and nodes seen.
struct Walk {
    next: u64,
    end: u64,
    blocks: u64,
    nodes: u64,
}

/// Visits the blocks of the R-tree index at `at` in order, checking that they fill `region`
/// (`[start, end)`) exactly and in order; returns how many there are. `visit` returns a
/// block's largest end; every index entry, node and header whose end is not the largest of
/// what it covers gets an [`EndPatch`].
fn walk_index(
    r: &mut Reader,
    at: u64,
    what: &str,
    region: (u64, u64),
    patches: &mut Vec<EndPatch>,
    mut visit: impl FnMut(&mut Reader, &Leaf) -> Result<End, String>,
) -> Result<u64, String> {
    let little = r.little;
    let h = r.read(at, CIR_TREE_HEADER, what)?;
    let mut f = Fields::new(h, little);
    if f.u32() != CIR_TREE_MAGIC {
        return Err(format!(
            "{what} at byte {} does not start with its magic number",
            thousands(at)
        ));
    }
    f.skip(4);
    let count = f.u64();
    f.skip(8);
    let stored = (f.u32(), f.u32());
    let (start, end) = region;
    let mut w = Walk {
        next: start,
        end,
        blocks: 0,
        nodes: 0,
    };
    let root = walk_node(
        r,
        at + CIR_TREE_HEADER as u64,
        what,
        0,
        &mut w,
        patches,
        &mut visit,
    )?;
    // The header repeats the root's bounds.
    if let Some(e) = root
        && e != stored
    {
        patches.push(EndPatch {
            at: at + 24,
            chrom: e.0,
            end: e.1,
        });
    }
    if w.blocks != count {
        return Err(format!(
            "{what} lists {} blocks, but its header says {}",
            thousands(w.blocks),
            thousands(count)
        ));
    }
    if w.next != end {
        return Err(format!(
            "{what}'s blocks end at byte {}, but its data runs to byte {}",
            thousands(w.next),
            thousands(end)
        ));
    }
    Ok(w.blocks)
}

/// The deepest an R-tree goes: 256 entries a node make 32 levels far more than any file.
const MAX_DEPTH: usize = 32;

/// Walks the R-tree node at `node` and what is below it, in order; returns the largest end
/// under it (`None` for an empty node).
fn walk_node<F>(
    r: &mut Reader,
    node: u64,
    what: &str,
    depth: usize,
    w: &mut Walk,
    patches: &mut Vec<EndPatch>,
    visit: &mut F,
) -> Result<Option<End>, String>
where
    F: FnMut(&mut Reader, &Leaf) -> Result<End, String>,
{
    w.nodes += 1;
    if depth > MAX_DEPTH || too_many_nodes(w.nodes, r.len) {
        return Err(format!("{what} loops"));
    }
    let little = r.little;
    let h = r.read(node, 4, what)?;
    let is_leaf = h[0] == 1;
    let n = usize::from(Fields::new(&h[2..], little).u16());
    let items = node + 4;
    let mut largest: Option<End> = None;
    if !is_leaf {
        let body = r.read(items, n * CIR_NODE_ITEM, what)?;
        let children: Vec<(End, u64)> = body
            .as_chunks::<CIR_NODE_ITEM>()
            .0
            .iter()
            .map(|e| {
                let mut f = Fields::new(e, little);
                f.skip(8);
                ((f.u32(), f.u32()), f.u64())
            })
            .collect();
        for (i, (stored, child)) in children.into_iter().enumerate() {
            let e = walk_node(r, child, what, depth + 1, w, patches, visit)?
                .ok_or_else(|| format!("{what} has an empty node at byte {}", thousands(child)))?;
            if e != stored {
                patches.push(EndPatch {
                    at: items + (i * CIR_NODE_ITEM) as u64 + 8,
                    chrom: e.0,
                    end: e.1,
                });
            }
            largest = largest.max(Some(e));
        }
        return Ok(largest);
    }
    let body = r.read(items, n * CIR_LEAF_ITEM, what)?;
    let leaves: Vec<Leaf> = body
        .as_chunks::<CIR_LEAF_ITEM>()
        .0
        .iter()
        .map(|e| {
            let mut f = Fields::new(e, little);
            Leaf {
                chrom: f.u32(),
                start: f.u32(),
                end_chrom: f.u32(),
                end: f.u32(),
                offset: f.u64(),
                size: f.u64(),
            }
        })
        .collect();
    for (i, leaf) in leaves.iter().enumerate() {
        if leaf.offset != w.next {
            return Err(format!(
                "{what} has a block {}, where the block before it ends at byte {}: blocks are \
missing or out of place",
                leaf.describe(),
                thousands(w.next)
            ));
        }
        match leaf.offset.checked_add(leaf.size) {
            Some(e) if leaf.size > 0 && e <= w.end => w.next = e,
            _ => {
                return Err(format!(
                    "{what} has a block {} that reaches past the end of its data, at byte {}",
                    leaf.describe(),
                    thousands(w.end)
                ));
            }
        }
        let e = visit(r, leaf)?;
        if e != (leaf.end_chrom, leaf.end) {
            patches.push(EndPatch {
                at: items + (i * CIR_LEAF_ITEM) as u64 + 8,
                chrom: e.0,
                end: e.1,
            });
        }
        largest = largest.max(Some(e));
        w.blocks += 1;
    }
    Ok(largest)
}

/// Reads blocks and inflates them, in buffers kept from one block to the next.
struct Inflater {
    /// The header's largest uncompressed block; 0 when blocks are stored as they are.
    max: usize,
    z: Decompress,
    out: Vec<u8>,
}

impl Inflater {
    fn new(max: u32) -> Inflater {
        let max = max as usize;
        Inflater {
            max,
            z: Decompress::new(true),
            out: Vec::with_capacity(max),
        }
    }

    /// The contents of the block `leaf` points at.
    fn load<'a>(
        &'a mut self,
        r: &'a mut Reader,
        leaf: &Leaf,
        what: &str,
    ) -> Result<&'a [u8], String> {
        let size = usize::try_from(leaf.size)
            .map_err(|_| format!("{what} {} is too large to read", leaf.describe()))?;
        let raw = r.read(leaf.offset, size, what)?;
        if self.max == 0 {
            return Ok(raw);
        }
        self.out.clear();
        self.z.reset(true);
        let status = self
            .z
            .decompress_vec(raw, &mut self.out, FlushDecompress::Finish)
            .map_err(|e| format!("{what} {} does not decompress: {e}", leaf.describe()))?;
        if status != Status::StreamEnd || self.z.total_in() != raw.len() as u64 {
            return Err(format!(
                "{what} {} does not decompress to a whole block of at most {} bytes",
                leaf.describe(),
                thousands(self.max as u64)
            ));
        }
        Ok(&self.out)
    }
}

/// The records a full-data block holds, checked against its index entry, and their largest
/// end.
fn block_records(kind: BbiKind, b: &[u8], leaf: &Leaf, little: bool) -> Result<(u64, End), String> {
    let bad = |why: String| format!("the full-data block {} {why}", leaf.describe());
    match kind {
        BbiKind::BigWig => {
            if b.len() < WIG_SECTION_HEADER {
                return Err(bad("is shorter than a section header".into()));
            }
            let mut f = Fields::new(&b[..WIG_SECTION_HEADER], little);
            let (chrom, start, end) = (f.u32(), f.u32(), f.u32());
            f.skip(8);
            let kind = f.u8();
            f.skip(1);
            let n = usize::from(f.u16());
            if (chrom, start, end) != (leaf.chrom, leaf.start, leaf.end) {
                return Err(bad(format!(
                    "has a section header for chromosome id {chrom} {start}-{end}, where the \
index says id {} {}-{}",
                    leaf.chrom, leaf.start, leaf.end
                )));
            }
            let item = match kind {
                1 => 12,
                2 => 8,
                3 => 4,
                _ => return Err(bad(format!("has section type {kind}, not 1, 2 or 3"))),
            };
            if b.len() != WIG_SECTION_HEADER + n * item {
                return Err(bad(format!(
                    "is {} bytes, where its {n} items take {}",
                    b.len(),
                    WIG_SECTION_HEADER + n * item
                )));
            }
            // The section's records do not overlap, so the last one's end is the largest.
            Ok((n as u64, (chrom, end)))
        }
        BbiKind::BigBed => {
            let mut at = 0;
            let mut n = 0u64;
            // The first record's chromosome and start, and the last's chromosome and end, which
            // the index entry repeats.
            let mut first = None;
            let mut last = (0, 0);
            let mut largest = (0, 0);
            while at < b.len() {
                if b.len() - at < 12 {
                    return Err(bad("ends inside a record".into()));
                }
                let mut f = Fields::new(&b[at..at + 12], little);
                let (chrom, start, end) = (f.u32(), f.u32(), f.u32());
                first.get_or_insert((chrom, start));
                last = (chrom, end);
                largest = largest.max((chrom, end));
                let rest = memchr::memchr(0, &b[at + 12..])
                    .ok_or_else(|| bad("ends inside a record".into()))?;
                at += 12 + rest + 1;
                n += 1;
            }
            let Some(first) = first else {
                return Err(bad("holds no records".into()));
            };
            // bigtools ends the entry at the last record; once corrected, at the largest end.
            let entry_end = (leaf.end_chrom, leaf.end);
            if first != (leaf.chrom, leaf.start) || (last != entry_end && largest != entry_end) {
                return Err(bad(format!(
                    "holds records from chromosome id {} {} to id {} {}, where the index says \
id {} {} to id {} {}",
                    first.0,
                    first.1,
                    last.0,
                    last.1,
                    leaf.chrom,
                    leaf.start,
                    leaf.end_chrom,
                    leaf.end
                )));
            }
            Ok((n, largest))
        }
    }
}

/// A zoom block holds whole zoom records, from its index entry's start to its end; their
/// largest end.
fn check_zoom_block(b: &[u8], leaf: &Leaf, little: bool) -> Result<End, String> {
    if b.is_empty() || !b.len().is_multiple_of(ZOOM_RECORD) {
        return Err(format!(
            "the zoom block {} is {} bytes, not whole {ZOOM_RECORD}-byte records",
            leaf.describe(),
            b.len()
        ));
    }
    let mut first = Fields::new(&b[..12], little);
    let (chrom, start) = (first.u32(), first.u32());
    let mut last = Fields::new(&b[b.len() - ZOOM_RECORD..][..12], little);
    let (end_chrom, _, end) = (last.u32(), last.u32(), last.u32());
    if (chrom, start, end_chrom, end) != (leaf.chrom, leaf.start, leaf.end_chrom, leaf.end) {
        return Err(format!(
            "the zoom block {} holds records from chromosome id {chrom} {start} to id \
{end_chrom} {end}, not what its index entry says",
            leaf.describe()
        ));
    }
    let mut largest = (0, 0);
    for rec in b.as_chunks::<ZOOM_RECORD>().0 {
        let mut f = Fields::new(rec, little);
        let c = f.u32();
        f.skip(4);
        largest = largest.max((c, f.u32()));
    }
    Ok(largest)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bbi::{BbiSizes, BbiWriteOptions, BbiWriter};
    use gloci_core::Genome;
    use std::collections::HashMap;
    use std::num::NonZeroUsize;

    const RECORDS: u64 = 3010;

    /// A file gloci wrote: 3,000 records on chr1 and 10 on chr2, so four sections and some
    /// zoom levels; its bytes.
    fn written(kind: BbiKind) -> Vec<u8> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(match kind {
            BbiKind::BigWig => "w.bw",
            BbiKind::BigBed => "w.bb",
        });
        let genome = Genome::new(vec![
            (b"chr1".to_vec(), Some(1_000_000)),
            (b"chr2".to_vec(), Some(50_000)),
        ])
        .unwrap();
        let spec = crate::output::parse_track_output_spec(Some(&path), None).unwrap();
        let mut w = BbiWriter::create(
            &spec,
            BbiWriteOptions {
                kind,
                sizes: BbiSizes {
                    genome,
                    source: "-g test".into(),
                },
                threads: NonZeroUsize::new(1).unwrap(),
                autosql: None,
                autosql_name: String::new(),
                defined_field_count: None,
            },
        )
        .unwrap();
        for (chrom, n) in [("chr1", 3000u64), ("chr2", 10)] {
            for i in 0..n {
                let (s, e) = (i * 10, i * 10 + 5);
                let line = match kind {
                    BbiKind::BigWig => format!("{chrom}\t{s}\t{e}\t{}", i % 7),
                    BbiKind::BigBed => format!("{chrom}\t{s}\t{e}\tn{i}\t{}\t+", i % 1000),
                };
                w.push(line.as_bytes()).unwrap();
            }
        }
        w.finish().unwrap();
        std::fs::read(&path).unwrap()
    }

    fn check(bytes: &[u8], kind: BbiKind, records: u64) -> Result<Checked, String> {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("c");
        std::fs::write(&p, bytes).unwrap();
        check_written(&p, kind, records)
    }

    fn u64_at(b: &[u8], at: usize) -> u64 {
        u64::from_le_bytes(b[at..at + 8].try_into().unwrap())
    }

    fn put_u64(b: &mut [u8], at: usize, v: u64) {
        b[at..at + 8].copy_from_slice(&v.to_le_bytes());
    }

    /// Where the entries of the R-tree index at `index` start (its root is a leaf in these
    /// small files), and how many there are.
    fn leaf_entries(b: &[u8], index: usize) -> (usize, usize) {
        let root = index + CIR_TREE_HEADER;
        assert_eq!(b[root], 1, "the root is a leaf");
        let n = u16::from_le_bytes([b[root + 2], b[root + 3]]) as usize;
        (root + 4, n)
    }

    /// The full-data index's entries.
    fn full_entries(b: &[u8]) -> (usize, usize) {
        leaf_entries(b, u64_at(b, 24) as usize)
    }

    /// The first zoom level's index entries.
    fn zoom_entries(b: &[u8]) -> (usize, usize) {
        leaf_entries(b, u64_at(b, HEADER + 16) as usize)
    }

    fn refused(bytes: &[u8], kind: BbiKind, needle: &str) {
        match check(bytes, kind, RECORDS) {
            Err(e) => assert!(e.contains(needle), "{e}"),
            Ok(c) => panic!("a damaged file passed: {c:?}"),
        }
    }

    #[test]
    fn a_written_file_passes() {
        for kind in [BbiKind::BigWig, BbiKind::BigBed] {
            let b = written(kind);
            let c = check(&b, kind, RECORDS).unwrap();
            assert_eq!(c.sections, 4, "{kind:?}");
            assert!(c.zoom_levels >= 1, "{kind:?}");
        }
    }

    #[test]
    fn a_record_count_that_differs_is_refused() {
        for kind in [BbiKind::BigWig, BbiKind::BigBed] {
            let e = check(&written(kind), kind, RECORDS + 1).unwrap_err();
            assert!(e.contains("3,010 records, but 3,011 were written"), "{e}");
        }
    }

    #[test]
    fn a_wrong_data_count_is_refused() {
        let mut b = written(BbiKind::BigWig);
        let at = u64_at(&b, 16) as usize;
        put_u64(&mut b, at, 5);
        refused(
            &b,
            BbiKind::BigWig,
            "its data count is 5, but it has 4 sections",
        );
    }

    /// One byte changed inside a block: it no longer decompresses.
    #[test]
    fn a_corrupt_block_is_refused() {
        for kind in [BbiKind::BigWig, BbiKind::BigBed] {
            let mut b = written(kind);
            let (e, _) = full_entries(&b);
            let (off, size) = (u64_at(&b, e + 16), u64_at(&b, e + 24));
            b[(off + size / 2) as usize] ^= 0x5a;
            refused(&b, kind, "does not decompress");
        }
    }

    /// What a section lost to a temporary file leaves: the index lists blocks that reach past
    /// the data the file holds.
    #[test]
    fn blocks_past_their_data_are_refused() {
        let mut b = written(BbiKind::BigWig);
        let (e, n) = full_entries(&b);
        let last = e + (n - 1) * CIR_LEAF_ITEM;
        let size = u64_at(&b, last + 24);
        put_u64(&mut b, last + 24, size + 16);
        refused(&b, BbiKind::BigWig, "reaches past the end of its data");
    }

    #[test]
    fn a_block_out_of_place_is_refused() {
        let mut b = written(BbiKind::BigBed);
        let (e, _) = full_entries(&b);
        let second = e + CIR_LEAF_ITEM;
        let off = u64_at(&b, second + 16);
        put_u64(&mut b, second + 16, off + 1);
        refused(&b, BbiKind::BigBed, "missing or out of place");
    }

    #[test]
    fn a_corrupt_zoom_block_is_refused() {
        let mut b = written(BbiKind::BigWig);
        let (e, _) = zoom_entries(&b);
        let (off, size) = (u64_at(&b, e + 16), u64_at(&b, e + 24));
        b[(off + size / 2) as usize] ^= 0x5a;
        refused(&b, BbiKind::BigWig, "a block of zoom level 1");
    }

    #[test]
    fn zoom_blocks_past_their_data_are_refused() {
        let mut b = written(BbiKind::BigBed);
        let (e, n) = zoom_entries(&b);
        let last = e + (n - 1) * CIR_LEAF_ITEM;
        let size = u64_at(&b, last + 24);
        put_u64(&mut b, last + 24, size + 8);
        refused(&b, BbiKind::BigBed, "the index of zoom level 1");
    }

    #[test]
    fn written_files_need_no_index_correction() {
        for kind in [BbiKind::BigWig, BbiKind::BigBed] {
            let c = check(&written(kind), kind, RECORDS).unwrap();
            assert!(c.patches.is_empty(), "{kind:?}: {:?}", c.patches);
        }
    }

    /// A bigBed as bigtools alone writes it, with 4 records a section and 4 entries a node,
    /// so that the index has leaf nodes under a root: an enclosing record and 40 inside it.
    fn raw_nested_bigbed(path: &Path) {
        use bigtools::bed::bedparser::BedValueError;
        use bigtools::beddata::BedParserStreamingIterator;
        use bigtools::{BedEntry, BigBedWrite, InputSortType};
        let entry = |start, end, rest: &str| {
            Ok::<_, BedValueError>((
                "chr1".to_string(),
                BedEntry {
                    start,
                    end,
                    rest: rest.to_string(),
                },
            ))
        };
        let mut items = vec![entry(0, 100_000, "outer")];
        items.extend((0..40u32).map(|i| entry(10 + i * 20, 15 + i * 20, "n")));
        let sizes = HashMap::from([("chr1".to_string(), 200_000u32)]);
        let mut w = BigBedWrite::new(File::create(path).unwrap(), sizes);
        w.options.input_sort_type = InputSortType::START;
        w.options.items_per_slot = 4;
        w.options.block_size = 4;
        let rt = tokio::runtime::Builder::new_current_thread()
            .build()
            .unwrap();
        w.write(
            BedParserStreamingIterator::wrap_iter(items.into_iter(), true),
            rt,
        )
        .unwrap();
    }

    /// The records a region read through bigtools' index returns.
    fn region(path: &Path, start: u32, end: u32) -> Vec<String> {
        let mut r = bigtools::BigBedRead::open_file(path).unwrap();
        r.get_interval("chr1", start, end)
            .unwrap()
            .map(|e| e.unwrap().rest)
            .collect()
    }

    #[test]
    fn nested_bigbed_index_ends_are_corrected() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("n.bb");
        raw_nested_bigbed(&path);
        // bigtools' index ends the first section at its last record, 55: a read past the
        // nested records misses the enclosing one.
        assert!(region(&path, 80_000, 80_100).is_empty());
        let c = check_written(&path, BbiKind::BigBed, 41).unwrap();
        let b = std::fs::read(&path).unwrap();
        let index = u64_at(&b, 24);
        let root = index + CIR_TREE_HEADER as u64;
        assert_eq!(b[root as usize], 0, "the root is a node above leaves");
        let first_child = u64_at(&b, root as usize + 4 + 16);
        let outer = EndPatch {
            chrom: 0,
            end: 100_000,
            at: 0,
        };
        let at = |at| EndPatch { at, ..outer };
        // The first leaf entry, its node's entry in the root, and the tree header.
        assert!(
            c.patches.contains(&at(first_child + 4 + 8)),
            "{:?}",
            c.patches
        );
        assert!(c.patches.contains(&at(root + 4 + 8)), "{:?}", c.patches);
        assert!(c.patches.contains(&at(index + 24)), "{:?}", c.patches);
        let mut f = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        for p in &c.patches {
            super::super::write::write_index_end(&mut f, true, p).unwrap();
        }
        drop(f);
        assert_eq!(region(&path, 80_000, 80_100), ["outer"]);
        assert_eq!(region(&path, 10, 16), ["outer", "n"]);
        let again = check_written(&path, BbiKind::BigBed, 41).unwrap();
        assert!(again.patches.is_empty(), "{:?}", again.patches);
    }

    #[test]
    fn a_truncated_file_is_refused() {
        let b = written(BbiKind::BigWig);
        let e = check(&b[..b.len() / 2], BbiKind::BigWig, RECORDS).unwrap_err();
        assert!(e.contains("outside the file"), "{e}");
    }
}
