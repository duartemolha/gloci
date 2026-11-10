//! BAM index building (D107): `.bai`, or `.csi` for `--csi` and for references longer than
//! 2^29. The spans are htslib's: a mapped read covers `[pos, pos + reference length)` (at
//! least one base), and a placed unmapped read `[pos, pos + 1)`, counted in the metadata's
//! unmapped total. Reads with no reference go to `n_no_coor`.

use super::build::{Linear, check_bgzf, csi_bytes, csi_depth, with_htslib_loffsets};
use super::*;
use crate::align::{Aln, AlnHeader, index_span};
use crate::error::OrderError;
use noodles_csi::binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex};

/// What a pass over the records learned.
struct BamPass<I> {
    /// `None` when a record ended beyond the range of the chosen depth.
    index: Option<BinningIndexOf<I>>,
    linear: Vec<Linear>,
    max_end: u64,
}

/// A read that breaks the coordinate order an index needs.
fn order_error(name: &str, (prev_n, prev): (u64, String), (n, rec): (u64, String)) -> IoError {
    IoError::Order(Box::new(OrderError {
        path: name.to_string(),
        prev_loc: Loc::Record(prev_n),
        prev,
        loc: Loc::Record(n),
        rec,
        hint: format!(
            "an index needs the reads in coordinate order, with the unplaced reads last; \
run: gloci sort {name} -o sorted.bam"
        ),
    }))
}

/// Read the header and every record of `data`, feeding an index of `depth`. Records must be
/// in coordinate order, with the unplaced reads last.
fn bam_pass<I: BinIndexKind + Default>(
    data: &Path,
    name: &str,
    threads: NonZeroUsize,
    depth: u8,
) -> Result<BamPass<I>, IoError> {
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let file = File::open(data).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let mut r = MultithreadedReader::with_worker_count(threads, file);
    let hdr = AlnHeader::read_bam(&mut r).map_err(|msg| IoError::Detect {
        path: name.to_string(),
        msg: format!("{name}: {msg}"),
    })?;
    let mut indexer = Indexer::<I>::new(MIN_SHIFT, depth).expect("depth is within 5..=9");
    let mut linear: Vec<Linear> = Vec::new();
    let mut full = false;
    let (mut max_end, mut n) = (0u64, 0u64);
    let (mut last_rid, mut last_pos, mut last_n) = (0usize, 0u64, 0u64);
    let place = |rid: usize, pos: u64| {
        format!(
            "read at {}:{}",
            String::from_utf8_lossy(&hdr.refs[rid].0),
            pos + 1
        )
    };
    let mut seen_placed = false;
    let mut unplaced = false;
    let mut buf = Vec::new();
    loop {
        let start = r.virtual_position();
        buf.clear();
        let got = bcf::append(&mut r, &mut buf, 4).map_err(read_err)?;
        if got == 0 {
            break;
        }
        n += 1;
        let err = |msg: String| IoError::Parse {
            path: name.to_string(),
            loc: Loc::Record(n),
            msg,
        };
        if got < 4 {
            return Err(err(bcf::RECORD_TRUNCATED.into()));
        }
        let block = i32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
        buf.clear();
        let size = usize::try_from(block).map_err(|_| {
            err(format!(
                "block_size {block} is negative, so the next record cannot be found; \
the file is corrupt"
            ))
        })?;
        if bcf::append(&mut r, &mut buf, size).map_err(read_err)? < size {
            return Err(err(bcf::RECORD_TRUNCATED.into()));
        }
        let end = r.virtual_position();
        let a = Aln::parse(&buf, &hdr).map_err(&err)?;
        let chunk = Chunk::new(start, end);
        let Some((rid, iv, mapped)) = index_span(&a) else {
            // No place: unplaced reads (refID -1) must come last, as htslib requires.
            if a.ref_id().is_none() {
                unplaced = true;
                last_n = n;
            }
            if !full {
                indexer
                    .add_record(None, chunk)
                    .expect("an unplaced record is always accepted");
            }
            continue;
        };
        if unplaced {
            return Err(order_error(
                name,
                (last_n, "a read with no reference".to_string()),
                (n, place(rid, iv.start)),
            ));
        }
        if seen_placed && (rid < last_rid || (rid == last_rid && iv.start < last_pos)) {
            return Err(order_error(
                name,
                (last_n, place(last_rid, last_pos)),
                (n, place(rid, iv.start)),
            ));
        }
        (seen_placed, last_rid, last_pos, last_n) = (true, rid, iv.start, n);
        if iv.end > MAX_INDEX_END {
            return Err(err(format!(
                "the read ends at {}, beyond {}, the largest position an index can hold",
                thousands(iv.end),
                thousands(MAX_INDEX_END)
            )));
        }
        max_end = max_end.max(iv.end);
        if full {
            continue;
        }
        let ctx = Some((rid, position(iv.start + 1), position(iv.end), mapped));
        match indexer.add_record(ctx, chunk) {
            Ok(()) => {}
            Err(AddRecordError::EndExceedsMaxPosition { .. }) => {
                full = true;
                continue;
            }
            Err(e) => return Err(err(format!("cannot index the read: {e}"))),
        }
        if linear.len() <= rid {
            linear.resize_with(rid + 1, Linear::default);
        }
        linear[rid].add(iv.end, chunk.start());
    }
    let index = (!full).then(|| indexer.build(hdr.refs.len()));
    Ok(BamPass {
        index,
        linear,
        max_end,
    })
}

fn bai_bytes(index: &noodles_bam::bai::Index) -> io::Result<Vec<u8>> {
    let mut w = noodles_bam::bai::io::Writer::new(Vec::new());
    w.write_index(index)?;
    Ok(w.into_inner())
}

/// Build the index of the BAM file `data` (BGZF with its EOF block) and return its kind and
/// bytes, without writing anything. `.bai`, unless `force_csi` or a reference is longer than
/// 2^29 (a read that ends beyond 2^29 on a shorter one also gets a `.csi`). `data` may be a
/// temp file that becomes `final_data_path`; messages name `final_data_path`.
pub(crate) fn build_bam_index(
    data: &Path,
    final_data_path: &Path,
    force_csi: bool,
    threads: NonZeroUsize,
) -> Result<(IndexKind, Vec<u8>), IoError> {
    let name = final_data_path.display().to_string();
    check_bgzf(data, &name)?;
    let werr = |kind| {
        move |source| IoError::Write {
            path: index_path(final_data_path, kind).display().to_string(),
            source,
        }
    };
    let max_ln = longest_reference(data, &name, threads)?;
    let mut csi_depth_for = None;
    if !force_csi && max_ln <= 1 << 29 {
        let p = bam_pass::<LinearIndex>(data, &name, threads, TBI_DEPTH)?;
        match p.index {
            Some(index) => {
                let bytes = bai_bytes(&index).map_err(werr(IndexKind::Bai))?;
                return Ok((IndexKind::Bai, bytes));
            }
            // A read ends beyond what a BAI holds: a CSI deep enough for it.
            None => csi_depth_for = Some(p.max_end),
        }
    }
    let mut depth = csi_depth(csi_depth_for.unwrap_or(max_ln));
    let p = loop {
        let p = bam_pass::<BinnedIndex>(data, &name, threads, depth)?;
        if p.index.is_some() {
            break p;
        }
        let deeper = csi_depth(p.max_end);
        if deeper <= depth {
            return Err(IoError::Detect {
                path: name.clone(),
                msg: format!("{name} changed while it was being indexed; run the command again"),
            });
        }
        depth = deeper;
    };
    let index = with_htslib_loffsets(p.index.expect("checked above"), &p.linear);
    let bytes = csi_bytes(&index).map_err(werr(IndexKind::Csi))?;
    Ok((IndexKind::Csi, bytes))
}

/// The longest reference of the BAM header, which picks `.bai` or `.csi`.
fn longest_reference(data: &Path, name: &str, threads: NonZeroUsize) -> Result<u64, IoError> {
    let file = File::open(data).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let mut r = MultithreadedReader::with_worker_count(threads, file);
    let hdr = AlnHeader::read_bam(&mut r).map_err(|msg| IoError::Detect {
        path: name.to_string(),
        msg: format!("{name}: {msg}"),
    })?;
    Ok(hdr.refs.iter().map(|(_, l)| *l).max().unwrap_or(0))
}

/// Write the BAM index of `data` (`.bai`, or `.csi` with `csi` or for long references)
/// atomically beside `final_data_path`, remove any other index of it, and return the path.
pub fn write_alignment_index(
    data: &Path,
    final_data_path: &Path,
    csi: bool,
    threads: NonZeroUsize,
) -> Result<PathBuf, IoError> {
    let (kind, bytes) = build_bam_index(data, final_data_path, csi, threads)?;
    super::build::publish_index(final_data_path, kind, &bytes)
}
