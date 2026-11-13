//! `gloci idxstats` (D127): per-reference read counts of a BAM, CRAM or SAM.

use gloci_io::{
    AlnIndexCounts, Format, Header, IndexKind, IndexedFile, IoError, ReaderOptions, RefCounts,
    find_index, open_aln_source,
};
use std::io::Write;
use std::path::Path;

/// Where the counts came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdxstatsPath {
    /// The metadata of a BAI or CSI: no record was read.
    Index,
    /// Every record was counted.
    Scan,
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// Write `name TAB LN TAB mapped TAB unmapped` per reference, then `* 0 0 unplaced`. Counts
/// come from the index of a BAM when it has them, else from reading every record: its flag and
/// reference id, which is all a CRAM decodes ([`AlnSource::for_each_flag_fields`]). `o` carries
/// `--reference` (for a CRAM that needs one) and the threads, and no read filters.
///
/// [`AlnSource::for_each_flag_fields`]: gloci_io::AlnSource::for_each_flag_fields
pub fn idxstats(
    path: &str,
    o: &ReaderOptions,
    out: &mut dyn Write,
) -> Result<IdxstatsPath, IoError> {
    let src = open_aln_source(path, o)?;
    let format = src.format();
    let name = src.name().to_string();
    let header = src.header().clone();
    let file = Some(Path::new(path)).filter(|p| path != "-" && p.is_file());
    let index = file
        .filter(|_| format == Format::Bam)
        .and_then(|p| find_index(p).map(|i| (p, i)))
        .filter(|(_, i)| matches!(i.kind, IndexKind::Bai | IndexKind::Csi));
    let mut src = Some(src);
    if let Some((p, info)) = index {
        info.warn_if_stale(&name);
        src = None;
        let ix = IndexedFile::open(p, format, Header::Alignment(header.clone()), o)?;
        if let Some(c) = ix.alignment_counts() {
            if c.unplaced.is_none() {
                log::warn!(
                    "the index {} has no count of unplaced reads; counting every record instead",
                    info.path.display()
                );
            } else if let Some(i) = c.per_ref.iter().position(Option::is_none) {
                log::warn!(
                    "the index {} has no read counts for {}; counting every record instead",
                    info.path.display(),
                    String::from_utf8_lossy(&header.refs[i].0)
                );
            } else {
                write_counts(&header.refs, &c, out)?;
                return Ok(IdxstatsPath::Index);
            }
        }
    } else if format == Format::Bam && file.is_some() {
        log::warn!(
            "{name} has no index; counting every record (index it for speed: gloci index {name})"
        );
    } else if format == Format::Bam {
        // stdin, a FIFO, /dev/fd/N: an index beside it could not be used, so none is suggested.
        log::warn!("{name} is not a regular file, so no index can be read; counting every record");
    } else {
        log::warn!(
            "{name}: a {} holds no per-reference counts; counting every record",
            format.name()
        );
    }
    let mut src = match src {
        Some(s) => s,
        None => open_aln_source(path, o)?,
    };
    let n = header.refs.len();
    let mut per_ref = vec![RefCounts::default(); n];
    let mut unplaced = 0u64;
    // Only the flag and the reference id count, so a CRAM decodes those data series alone
    // (as flagstat does) and not the names, bases, qualities and tags of every read.
    src.for_each_flag_fields(&mut |r| {
        match r.ref_id {
            Some(rid) if rid < n => {
                if r.flag & 0x4 != 0 {
                    per_ref[rid].unmapped += 1;
                } else {
                    per_ref[rid].mapped += 1;
                }
            }
            // The sources check refID against the header, so this is a bug.
            Some(rid) => {
                return Err(IoError::Internal(format!(
                    "{name}: record with refID {rid} passed the header's {n} references"
                )));
            }
            None => unplaced += 1,
        }
        Ok(())
    })?;
    let counts = AlnIndexCounts {
        per_ref: per_ref.into_iter().map(Some).collect(),
        unplaced: Some(unplaced),
    };
    write_counts(&header.refs, &counts, out)?;
    Ok(IdxstatsPath::Scan)
}

fn write_counts(
    refs: &[(Vec<u8>, u64)],
    c: &AlnIndexCounts,
    out: &mut dyn Write,
) -> Result<(), IoError> {
    let mut line = Vec::with_capacity(64);
    for ((name, len), rc) in refs.iter().zip(&c.per_ref) {
        let rc = rc.unwrap_or_default();
        line.clear();
        line.extend_from_slice(name);
        writeln!(line, "\t{len}\t{}\t{}", rc.mapped, rc.unmapped).map_err(write_err)?;
        out.write_all(&line).map_err(write_err)?;
    }
    writeln!(out, "*\t0\t0\t{}", c.unplaced.unwrap_or(0)).map_err(write_err)?;
    out.flush().map_err(write_err)
}
