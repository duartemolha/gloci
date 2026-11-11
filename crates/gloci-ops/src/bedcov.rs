//! `gloci bedcov`: the depth of coverage summed over each BED region, as `samtools bedcov`
//! (D140). Each region is a random-access query through the alignment index, counted by
//! [`DepthRing`] (D137's rules), so memory follows the longest read and not the region.
//!
//! Batches of BED lines are read in order, counted on worker threads, and written in order,
//! so the output and what stderr says do not depend on the number of threads.

use crate::depth::{CountRules, DepthRing, check_same_references, log_diagnostics};
use gloci_core::{ChromTable, Interval};
use gloci_io::{
    AlnHeader, AlnSource, Diagnostics, Format, Header, IndexedFile, InvalidPolicy, IoError, Loc,
    OrderError, ReaderOptions, RecordBatch, RecordReader, UsageError, aln_ref, find_index,
    next_batch_size, open_aln_source,
};
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver};
use std::sync::{Arc, Mutex};

const FLAG_UNMAPPED: u16 = 0x4;
/// BED lines in a batch, at most (and at most the 4 MiB `next_batch_size` allows).
const BATCH_LINES: usize = 256;
/// Bytes of BED lines in a batch, at most.

/// What one input contributes to one region.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RegionTotals {
    /// The depth summed over the region's bases.
    pub sum: u64,
    /// The bases of the region whose depth is at least the threshold (0 without one).
    pub at_depth: u64,
    /// The passing reads whose span overlaps the region.
    pub reads: u64,
}

#[derive(Debug, Clone)]
pub struct BedcovOptions {
    /// What counts toward a base's depth (`--split` drops `ref_skips`, `--split-deletions`
    /// drops `deletions`).
    pub rules: CountRules,
    /// `--bases-at-depth N`: add a column of the bases with depth at least N.
    pub bases_at_depth: Option<u32>,
    /// `--read-count`: add a column of the reads over each region.
    pub read_count: bool,
    /// Write `#chrom start end <path>_cov...` first.
    pub column_names: bool,
}

/// The totals of the region `iv` from `src`, whose records are the reads that overlap it, in
/// file order (see `IndexedFile::aln_regions_source`). `ref_len` is the length of the region's
/// reference: a read is counted up to it. A placed unmapped read never counts. A zero-length
/// region has no base and, as in samtools 1.21, no read: all three totals are 0.
pub fn region_totals(
    src: &mut dyn AlnSource,
    ref_len: u64,
    iv: Interval,
    rules: &CountRules,
    at_depth: Option<u32>,
) -> Result<RegionTotals, IoError> {
    region_totals_in(
        &mut DepthRing::new(),
        src,
        ref_len,
        iv,
        rules,
        at_depth,
        &AtomicBool::new(false),
    )
}

/// [`region_totals`] with the caller's ring (reset here, so its memory is reused) and the flag
/// that makes the past-`LN` warning once per run.
pub fn region_totals_in(
    ring: &mut DepthRing,
    src: &mut dyn AlnSource,
    ref_len: u64,
    iv: Interval,
    rules: &CountRules,
    at_depth: Option<u32>,
    warned_past: &AtomicBool,
) -> Result<RegionTotals, IoError> {
    ring.reset();
    if iv.start >= iv.end {
        return Ok(RegionTotals::default());
    }
    let (layout, hdr, name) = (src.layout(), src.header().clone(), src.name().to_string());
    let mut t = RegionTotals::default();
    let mut buf = Vec::new();
    // The first base of the region not summed yet.
    let mut cursor = iv.start;
    let mut sweep = |ring: &mut DepthRing, upto: u64, t: &mut RegionTotals| {
        let stop = upto.min(iv.end);
        let hi = stop.min(ring.covered_end());
        for p in cursor..hi {
            let d = ring.take(p);
            t.sum += u64::from(d);
            if at_depth.is_some_and(|n| d >= n) {
                t.at_depth += 1;
            }
        }
        cursor = cursor.max(stop);
    };
    let mut last: Option<(u64, u64)> = None;
    loop {
        buf.clear();
        let Some(rec) = src.next_record(&mut buf)? else {
            break;
        };
        let a = aln_ref(layout, &buf, &hdr).map_err(|msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(rec),
            msg,
        })?;
        // The pileup skips unmapped reads, so they never count.
        if a.flag() & FLAG_UNMAPPED != 0 {
            continue;
        }
        let Some(pos) = a.pos() else { continue };
        let span = match a.span(ref_len) {
            Some((span, clipped)) => {
                if clipped && !warned_past.swap(true, Ordering::Relaxed) {
                    log::warn!(
                        "{name}: a read ends past the end of its reference (LN {ref_len}); \
bedcov counts only the bases up to LN, where samtools counts those beyond it too (this warning \
is given once)"
                    );
                }
                span
            }
            None => continue,
        };
        if span.start >= span.end {
            continue;
        }
        if let Some((prec, ppos)) = last
            && pos < ppos
        {
            let at = |p: u64| format!("{}:{}", iv_name(&hdr, &a), p + 1);
            return Err(OrderError {
                path: name.clone(),
                prev_loc: Loc::Record(prec),
                prev: at(ppos),
                loc: Loc::Record(rec),
                rec: at(pos),
                hint: format!("run: gloci sort {name} -o sorted.bam, and use that file"),
            }
            .into());
        }
        last = Some((rec, pos));
        t.reads += 1;
        sweep(ring, span.start, &mut t);
        ring.add(&a, span, rules, 0);
    }
    sweep(ring, iv.end, &mut t);
    Ok(t)
}

fn iv_name(hdr: &AlnHeader, a: &gloci_io::AlnRef<'_>) -> String {
    a.ref_id()
        .and_then(|r| hdr.refs.get(r))
        .map(|(n, _)| String::from_utf8_lossy(n).into_owned())
        .unwrap_or_default()
}

/// One batch of BED lines on its way to a worker.
struct Job {
    id: u64,
    batch: RecordBatch,
    /// The reference id of each record's chromosome.
    rids: Vec<usize>,
}

struct Done {
    out: Vec<u8>,
    lines: u64,
    diags: Vec<Diagnostics>,
}

struct Shared<'a> {
    files: &'a [IndexedFile],
    refs: &'a [(Vec<u8>, u64)],
    o: &'a BedcovOptions,
    warned_past: &'a AtomicBool,
}

fn push_u64(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(v.to_string().as_bytes());
}

/// Count every line of one batch, writing the lines with their columns.
fn run_batch(sh: &Shared<'_>, job: &Job, ring: &mut DepthRing) -> Result<Done, IoError> {
    let n_files = sh.files.len();
    let mut diags = vec![Diagnostics::default(); n_files];
    let mut out = Vec::with_capacity(job.batch.buf.len() + job.batch.len() * 24 * n_files);
    let mut totals = vec![RegionTotals::default(); n_files];
    for (i, m) in job.batch.recs.iter().enumerate() {
        let rid = job.rids[i];
        let (name, ref_len) = (&sh.refs[rid].0, sh.refs[rid].1);
        let region = [(name.clone(), m.span)];
        for (k, f) in sh.files.iter().enumerate() {
            if m.span.start >= m.span.end {
                totals[k] = RegionTotals::default();
                continue;
            }
            let mut src = f.aln_regions_source(&region, false)?;
            totals[k] = region_totals_in(
                ring,
                &mut *src,
                ref_len,
                m.span,
                &sh.o.rules,
                sh.o.bases_at_depth,
                sh.warned_past,
            )?;
            diags[k].absorb(src.diagnostics().clone());
        }
        out.extend_from_slice(job.batch.raw(i));
        for t in &totals {
            out.push(b'\t');
            push_u64(&mut out, t.sum);
        }
        if sh.o.bases_at_depth.is_some() {
            for t in &totals {
                out.push(b'\t');
                push_u64(&mut out, t.at_depth);
            }
        }
        if sh.o.read_count {
            for t in &totals {
                out.push(b'\t');
                push_u64(&mut out, t.reads);
            }
        }
        out.push(b'\n');
    }
    Ok(Done {
        out,
        lines: job.batch.len() as u64,
        diags,
    })
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// The header line `--column-names` writes.
fn header_line(files: &[String], o: &BedcovOptions) -> Vec<u8> {
    let mut line = b"#chrom\tstart\tend".to_vec();
    let mut group = |suffix: &str| {
        for p in files {
            line.push(b'\t');
            line.extend_from_slice(p.as_bytes());
            line.extend_from_slice(suffix.as_bytes());
        }
    };
    group("_cov");
    if o.bases_at_depth.is_some() {
        group("_depth");
    }
    if o.read_count {
        group("_count");
    }
    line.push(b'\n');
    line
}

/// `samtools bedcov` of the regions of `bed` over `files`, written to `out`. Returns the number of
/// lines written.
pub fn bedcov(
    mut bed: Box<dyn RecordReader>,
    files: &[String],
    o: &BedcovOptions,
    ro: &ReaderOptions,
    threads: NonZeroUsize,
    out: &mut dyn Write,
) -> Result<u64, IoError> {
    if files.is_empty() {
        return Err(UsageError("bedcov needs at least one alignment file".into()).into());
    }
    let mut headers: Vec<Arc<AlnHeader>> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut indexed: Vec<IndexedFile> = Vec::new();
    for p in files {
        let src = open_aln_source(p, ro)?;
        let name = src.name().to_string();
        let format = src.format();
        let hdr = src.header().clone();
        drop(src);
        let file = Some(Path::new(p)).filter(|f| p != "-" && f.is_file());
        let info = file
            .filter(|_| matches!(format, Format::Bam | Format::Cram))
            .and_then(|f| find_index(f).map(|i| (f, i)));
        let Some((f, info)) = info else {
            let how = match format {
                Format::Bam | Format::Cram if file.is_some() => {
                    format!("create an index with: gloci index {name}")
                }
                Format::Bam | Format::Cram => {
                    "write the input to a file and index it with: gloci index FILE".to_string()
                }
                _ => format!(
                    "make a sorted BAM with gloci sort {name} -o sorted.bam, then index it with: \
gloci index sorted.bam"
                ),
            };
            return Err(UsageError(format!(
                "bedcov reads each region through an index, and {name} has none; {how}"
            ))
            .into());
        };
        info.warn_if_stale(&name);
        indexed.push(IndexedFile::open(
            f,
            format,
            Header::Alignment(hdr.clone()),
            ro,
        )?);
        headers.push(hdr);
        names.push(name);
    }
    check_same_references(&names, &headers)?;
    let refs = &headers[0].refs;
    let rid_of: HashMap<&[u8], usize> = refs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (&n[..], i))
        .collect();
    if o.column_names {
        out.write_all(&header_line(files, o)).map_err(write_err)?;
    }

    let warned_past = AtomicBool::new(false);
    let sh = Shared {
        files: &indexed,
        refs,
        o,
        warned_past: &warned_past,
    };
    let workers = threads.get();
    let bed_name = bed.name().to_string();
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut max = BATCH_LINES;
    let mut next_id = 0u64;
    let mut more = true;
    let mut read_err: Option<IoError> = None;
    let (mut skipped, mut first_skipped) = (0u64, None::<String>);
    let mut merged: Vec<Diagnostics> = vec![Diagnostics::default(); files.len()];
    let mut lines = 0u64;

    // Read the next batch of BED lines whose chromosomes the references list.
    let mut next_job = |bed: &mut Box<dyn RecordReader>,
                        read_err: &mut Option<IoError>,
                        more: &mut bool|
     -> Result<Option<Job>, IoError> {
        loop {
            if !*more || read_err.is_some() {
                return Ok(None);
            }
            batch.clear();
            *more = bed.read_batch(&mut chroms, &mut batch, max)?;
            max = next_batch_size(&batch, max).min(BATCH_LINES);
            let mut rids = Vec::with_capacity(batch.len());
            let mut keep = RecordBatch::default();
            for (i, m) in batch.recs.iter().enumerate() {
                let chrom = chroms.name(m.chrom);
                match rid_of.get(chrom) {
                    Some(&rid) => {
                        rids.push(rid);
                        let off = keep.buf.len();
                        keep.buf.extend_from_slice(batch.raw(i));
                        let mut m = *m;
                        m.off = off;
                        keep.recs.push(m);
                    }
                    None => {
                        let shown = String::from_utf8_lossy(chrom);
                        if ro.invalid == InvalidPolicy::Skip {
                            skipped += 1;
                            first_skipped.get_or_insert_with(|| {
                                format!("{bed_name} line {}: {shown}", m.loc)
                            });
                            continue;
                        }
                        *read_err = Some(IoError::Parse {
                            path: bed_name.clone(),
                            loc: Loc::Line(m.loc),
                            msg: format!(
                                "{shown} is not a reference of {}; fix the BED line, rename \
the chromosome with --chrom-map, or pass --skip-invalid to skip lines on unknown references",
                                names[0]
                            ),
                        });
                        break;
                    }
                }
            }
            if keep.is_empty() && read_err.is_none() {
                continue;
            }
            let id = next_id;
            next_id += 1;
            return Ok(Some(Job {
                id,
                batch: keep,
                rids,
            }));
        }
    };

    let absorb = |d: Done, merged: &mut Vec<Diagnostics>, lines: &mut u64| {
        *lines += d.lines;
        for (m, x) in merged.iter_mut().zip(d.diags) {
            m.absorb(x);
        }
    };

    let result: Result<(), IoError> = if workers == 1 {
        let mut ring = DepthRing::new();
        loop {
            let Some(job) = next_job(&mut bed, &mut read_err, &mut more)? else {
                break;
            };
            let done = run_batch(&sh, &job, &mut ring)?;
            out.write_all(&done.out).map_err(write_err)?;
            absorb(done, &mut merged, &mut lines);
        }
        Ok(())
    } else {
        let in_flight = 2 * workers;
        std::thread::scope(|s| {
            let (job_tx, job_rx) = mpsc::channel::<Job>();
            let job_rx: Arc<Mutex<Receiver<Job>>> = Arc::new(Mutex::new(job_rx));
            let (res_tx, res_rx) = mpsc::channel::<(u64, Result<Done, IoError>)>();
            for _ in 0..workers {
                let (rx, tx, sh) = (job_rx.clone(), res_tx.clone(), &sh);
                s.spawn(move || {
                    let mut ring = DepthRing::new();
                    loop {
                        let job = match rx.lock().unwrap().recv() {
                            Ok(j) => j,
                            Err(_) => break,
                        };
                        let r = run_batch(sh, &job, &mut ring);
                        if tx.send((job.id, r)).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(res_tx);
            let mut pending: BTreeMap<u64, Result<Done, IoError>> = BTreeMap::new();
            let (mut outstanding, mut next_out) = (0usize, 0u64);
            let mut failure: Option<IoError> = None;
            // Write what is ready, in order; the first error stops the writing.
            let mut drain = |pending: &mut BTreeMap<u64, Result<Done, IoError>>,
                             failure: &mut Option<IoError>,
                             next_out: &mut u64|
             -> Result<(), IoError> {
                while failure.is_none()
                    && let Some(r) = pending.remove(next_out)
                {
                    *next_out += 1;
                    match r {
                        Ok(d) => {
                            out.write_all(&d.out).map_err(write_err)?;
                            absorb(d, &mut merged, &mut lines);
                        }
                        Err(e) => *failure = Some(e),
                    }
                }
                Ok(())
            };
            let mut fatal: Option<IoError> = None;
            loop {
                if failure.is_some() || fatal.is_some() {
                    break;
                }
                while outstanding >= in_flight {
                    let Ok((id, r)) = res_rx.recv() else { break };
                    outstanding -= 1;
                    pending.insert(id, r);
                    if let Err(e) = drain(&mut pending, &mut failure, &mut next_out) {
                        fatal = Some(e);
                        break;
                    }
                }
                if failure.is_some() || fatal.is_some() {
                    break;
                }
                match next_job(&mut bed, &mut read_err, &mut more) {
                    Ok(Some(job)) => {
                        outstanding += 1;
                        let _ = job_tx.send(job);
                    }
                    Ok(None) => break,
                    Err(e) => {
                        fatal = Some(e);
                        break;
                    }
                }
            }
            drop(job_tx);
            while outstanding > 0 {
                let Ok((id, r)) = res_rx.recv() else { break };
                outstanding -= 1;
                pending.insert(id, r);
                if fatal.is_none()
                    && let Err(e) = drain(&mut pending, &mut failure, &mut next_out)
                {
                    fatal = Some(e);
                }
            }
            if let Some(e) = failure {
                return Err(e);
            }
            match fatal {
                Some(e) => Err(e),
                None => Ok(()),
            }
        })
    };
    result?;
    if let Some(e) = read_err {
        return Err(e);
    }
    out.flush().map_err(write_err)?;
    if let Some(s) = bed.diagnostics().summary(&bed_name) {
        log::warn!("{s}");
    }
    if skipped > 0 {
        log::warn!(
            "skipped {skipped} BED line{} on references {} does not list (first: {})",
            if skipped == 1 { "" } else { "s" },
            names[0],
            first_skipped.unwrap_or_default()
        );
    }
    for (n, d) in names.iter().zip(&merged) {
        log_diagnostics(n, d);
    }
    Ok(lines)
}
