//! Per-chromosome reading of plain, seekable files (D76): one pass over an
//! uncompressed BED, bedGraph or VCF finds where each chromosome's records
//! lie, and readers then take byte ranges of the file, so chromosomes can be
//! processed in parallel as they are for indexed files. A range reader gives
//! the records, line numbers and errors the whole-file reader gives for the
//! same lines.

use crate::bed::{Bed, is_bed_record_line};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::IoError;
use crate::record::{Header, ReaderOptions, RecordBatch, RecordReader, TextKind, TextReader};
use crate::vcf::Vcf;
use gloci_core::ChromTable;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Bytes read per step of the scan.
const SCAN_BUF: usize = 4 << 20;
/// Lines after which a file with more chromosomes than half its lines so far
/// is read whole.
const PATHOLOGICAL_AFTER: u64 = 10_000;
/// Read buffer of a range reader.
const RANGE_BUF: usize = 256 << 10;

/// One chromosome's records in a file: the bytes from its first record to
/// the first record of the next chromosome (or the end of the file), which
/// may include blank and comment lines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChromRange {
    /// The name as the file writes it (before `--chrom-map`).
    pub name: Vec<u8>,
    pub bytes: Range<u64>,
    /// 1-based number of the line the range starts on.
    pub first_line: u64,
}

/// Whether `line` (no terminator) holds a record of `format`.
fn is_record(format: Format, line: &[u8]) -> bool {
    match format {
        Format::Vcf => line.first().is_some_and(|&b| b != b'#'),
        _ => is_bed_record_line(line),
    }
}

/// Whether the plain-file scan supports `format`: BED, bedGraph and VCF.
/// GFF may carry a `##FASTA` section after its records and WIG names its
/// chromosome on declaration lines, so those are read whole.
pub fn range_format(format: Format) -> bool {
    matches!(format, Format::Bed | Format::BedGraph | Format::Vcf)
}

/// One pass over `path` (BED, bedGraph or VCF text) to find where each
/// chromosome's records lie, in file order, with the line each range starts
/// on. `None` when a chromosome's records are not contiguous: the caller
/// reads the file whole, which reports that as the sortedness error with
/// its line numbers. A record line's chromosome is the text before its
/// first tab. `None` also for a record line with no tab (which no reader
/// accepts) and for a file with more chromosomes than half its lines (after
/// the first 10,000 lines): neither gains from being split.
pub fn chrom_ranges(path: &Path, format: Format) -> Result<Option<Vec<ChromRange>>, IoError> {
    scan(path, format, SCAN_BUF)
}

fn scan(path: &Path, format: Format, buf_size: usize) -> Result<Option<Vec<ChromRange>>, IoError> {
    debug_assert!(range_format(format));
    let name = path.display().to_string();
    let read_err = |source| IoError::Read {
        path: name.clone(),
        source,
    };
    let mut file = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let len = file.metadata().map_err(read_err)?.len();
    let mut ranges: Vec<ChromRange> = Vec::new();
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut buf = vec![0u8; buf_size];
    // `buf[pos..filled]` is unread; `base` is the file offset of `buf[0]`.
    let (mut filled, mut pos, mut base, mut line_no) = (0usize, 0usize, 0u64, 0u64);
    let mut eof = false;
    let mut line = |at: u64, no: u64, text: &[u8]| -> bool {
        let text = text.strip_suffix(b"\r").unwrap_or(text);
        if !is_record(format, text) {
            return true;
        }
        // A record line with no tab is malformed: no reader accepts it, and
        // as a chromosome of its own it could make the ranges as many as the
        // lines. Read the file whole and let the reader report it.
        let Some(tab) = memchr::memchr(b'\t', text) else {
            return false;
        };
        let chrom = &text[..tab];
        if ranges.last().is_some_and(|r| r.name == chrom) {
            return true;
        }
        if !seen.insert(chrom.to_vec()) {
            return false;
        }
        // Past the first lines, more chromosomes than half the lines means a
        // file that is not a few long runs: nothing to gain from splitting it.
        if no > PATHOLOGICAL_AFTER && ranges.len() as u64 * 2 > no {
            return false;
        }
        if let Some(prev) = ranges.last_mut() {
            prev.bytes.end = at;
        }
        ranges.push(ChromRange {
            name: chrom.to_vec(),
            bytes: at..len,
            first_line: no,
        });
        true
    };
    loop {
        while let Some(i) = memchr::memchr(b'\n', &buf[pos..filled]) {
            line_no += 1;
            if !line(base + pos as u64, line_no, &buf[pos..pos + i]) {
                return Ok(None);
            }
            pos += i + 1;
        }
        if eof {
            if pos < filled && !line(base + pos as u64, line_no + 1, &buf[pos..filled]) {
                return Ok(None);
            }
            return Ok(Some(ranges));
        }
        buf.copy_within(pos..filled, 0);
        base += pos as u64;
        filled -= pos;
        pos = 0;
        if filled == buf.len() {
            // One line longer than the buffer.
            buf.resize(buf.len() * 2, 0);
        }
        loop {
            match (&mut file)
                .take(len - base - filled as u64)
                .read(&mut buf[filled..])
            {
                Ok(0) => eof = true,
                Ok(n) => filled += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(read_err(e)),
            }
            break;
        }
    }
}

/// The records of the byte ranges of `path`, one after another, as the
/// whole-file reader gives them: `ranges` must be whole-record ranges of a
/// file in `format` that has `header`, as [`chrom_ranges`] finds them.
/// `--chrom-map` and `--skip-invalid` act as for the whole-file reader.
pub fn open_range_reader(
    path: &Path,
    ranges: &[ChromRange],
    format: Format,
    header: &Header,
    o: &ReaderOptions,
) -> Result<Box<dyn RecordReader>, IoError> {
    let ranges = ranges.to_vec();
    let name = path.display().to_string();
    let (header, policy, map) = (header.clone(), o.invalid, o.chrom_map.clone());
    Ok(match format {
        Format::Vcf => Box::new(RangeReader::new(
            path,
            name,
            format,
            header,
            ranges,
            Vcf::new(o.vcf_span),
            policy,
            map,
        )),
        _ => Box::new(RangeReader::new(
            path, name, format, header, ranges, Bed, policy, map,
        )),
    })
}

struct RangeReader<K: TextKind + Clone> {
    path: PathBuf,
    name: String,
    format: Format,
    header: Header,
    policy: InvalidPolicy,
    map: Option<Arc<ChromMap>>,
    kind: K,
    ranges: Vec<ChromRange>,
    next: usize,
    cur: Option<TextReader<K>>,
    /// The diagnostics of the ranges read to the end, and with the current
    /// one's added.
    done: Diagnostics,
    diag: Diagnostics,
    /// An error met after records of earlier ranges were put in the batch:
    /// they go out first, and the next call raises it.
    deferred: Option<IoError>,
}

impl<K: TextKind + Clone> RangeReader<K> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        path: &Path,
        name: String,
        format: Format,
        header: Header,
        ranges: Vec<ChromRange>,
        kind: K,
        policy: InvalidPolicy,
        map: Option<Arc<ChromMap>>,
    ) -> Self {
        RangeReader {
            path: path.to_path_buf(),
            name,
            format,
            header,
            policy,
            map,
            kind,
            ranges,
            next: 0,
            cur: None,
            done: Diagnostics::default(),
            diag: Diagnostics::default(),
            deferred: None,
        }
    }

    /// Start the next range; false when there is none.
    fn open_next(&mut self) -> Result<bool, IoError> {
        let Some(r) = self.ranges.get(self.next) else {
            return Ok(false);
        };
        self.next += 1;
        let io = |source| IoError::Read {
            path: self.name.clone(),
            source,
        };
        let mut f = File::open(&self.path).map_err(|source| IoError::Open {
            path: self.name.clone(),
            source,
        })?;
        f.seek(SeekFrom::Start(r.bytes.start)).map_err(io)?;
        let src: Box<dyn BufRead + Send> = Box::new(BufReader::with_capacity(
            RANGE_BUF,
            f.take(r.bytes.end - r.bytes.start),
        ));
        self.cur = Some(TextReader::over_range(
            &self.name,
            self.format,
            self.header.clone(),
            self.policy,
            src,
            r.first_line,
            self.kind.clone(),
            self.map.clone(),
        ));
        Ok(true)
    }

    fn sync_diag(&mut self) {
        self.diag = self.done.clone();
        if let Some(c) = &self.cur {
            self.diag.absorb(c.diagnostics().clone());
        }
    }
}

impl<K: TextKind + Clone> RecordReader for RangeReader<K> {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        self.format
    }
    fn compression(&self) -> Compression {
        Compression::None
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }
    fn consume_to_end(&self) -> bool {
        true
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        if let Some(e) = self.deferred.take() {
            return Err(e);
        }
        let before = batch.recs.len();
        let result = loop {
            if self.cur.is_none() && !self.open_next()? {
                break Ok(false);
            }
            let want = max - (batch.recs.len() - before);
            let cur = self.cur.as_mut().expect("opened above");
            match cur.read_batch(chroms, batch, want) {
                Err(e) if batch.recs.len() > before => {
                    self.deferred = Some(e);
                    break Ok(true);
                }
                Err(e) => break Err(e),
                Ok(true) => break Ok(true),
                Ok(false) => {
                    self.done.absorb(cur.diagnostics().clone());
                    self.cur = None;
                    if batch.recs.len() - before >= max {
                        break Ok(self.next < self.ranges.len());
                    }
                }
            }
        };
        self.sync_diag();
        result
    }
}

/// A plain file whose chromosomes' byte ranges are known, from which
/// readers of whole chromosomes are made cheaply (D76); the counterpart of
/// [`crate::IndexedFile`].
pub struct SeekableFile {
    path: PathBuf,
    name: String,
    format: Format,
    header: Header,
    opts: ReaderOptions,
    /// In file order; `names[i]` is `ranges[i]`'s name after `--chrom-map`.
    ranges: Vec<ChromRange>,
    names: Vec<Vec<u8>>,
    by_name: HashMap<Vec<u8>, usize>,
}

impl SeekableFile {
    /// Scan `path` (`format` and `header` already read from it, for example
    /// by the reader the command opened). `None` when it cannot be read per
    /// chromosome: its chromosomes are not contiguous, or `--chrom-map` gives
    /// two of them one name. The whole-file reader then reports whatever is
    /// wrong, or reads it as it is.
    pub fn open(
        path: &Path,
        format: Format,
        header: Header,
        o: &ReaderOptions,
    ) -> Result<Option<SeekableFile>, IoError> {
        let Some(ranges) = chrom_ranges(path, format)? else {
            return Ok(None);
        };
        let map = o.chrom_map.as_deref();
        let names: Vec<Vec<u8>> = ranges
            .iter()
            .map(|r| map.map_or(&r.name[..], |m| m.rename(&r.name)).to_vec())
            .collect();
        let by_name: HashMap<Vec<u8>, usize> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        if by_name.len() != names.len() {
            return Ok(None);
        }
        Ok(Some(SeekableFile {
            path: path.to_path_buf(),
            name: path.display().to_string(),
            format,
            header,
            opts: o.clone(),
            ranges,
            names,
            by_name,
        }))
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// The chromosomes with records, as the records name them (after
    /// `--chrom-map`), in file order.
    pub fn reference_names(&self) -> Vec<Vec<u8>> {
        self.names.clone()
    }

    /// How much of the file `chrom`'s records take, in bytes (0 when it has
    /// none).
    pub fn weight(&self, chrom: &[u8]) -> u64 {
        self.by_name.get(chrom).map_or(0, |&i| {
            self.ranges[i].bytes.end - self.ranges[i].bytes.start
        })
    }

    /// Every record on `chroms` (named as records are after `--chrom-map`),
    /// chromosome by chromosome in the order given. Chromosomes the file
    /// lacks yield nothing.
    pub fn chromosomes_reader(&self, chroms: &[Vec<u8>]) -> Result<Box<dyn RecordReader>, IoError> {
        let ranges: Vec<ChromRange> = chroms
            .iter()
            .filter_map(|c| self.by_name.get(c))
            .map(|&i| self.ranges[i].clone())
            .collect();
        open_range_reader(&self.path, &ranges, self.format, &self.header, &self.opts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::{next_batch_size, open_reader};
    use gloci_core::{ChromTable, VcfSpan};
    use std::num::NonZeroUsize;

    fn opts() -> ReaderOptions {
        ReaderOptions {
            threads: NonZeroUsize::new(1).unwrap(),
            format: None,
            vcf_span: VcfSpan::Allele,
            invalid: InvalidPolicy::Skip,
            chrom_map: None,
        }
    }

    /// Every record as `(chrom name, key, loc, raw)`.
    fn read_all(mut r: Box<dyn RecordReader>) -> Vec<(Vec<u8>, u64, u64, Vec<u8>)> {
        let mut chroms = ChromTable::first_appearance();
        let (mut batch, mut max, mut out) = (RecordBatch::default(), 3, Vec::new());
        loop {
            batch.clear();
            let more = r.read_batch(&mut chroms, &mut batch, max).unwrap();
            for (i, m) in batch.recs.iter().enumerate() {
                out.push((
                    chroms.name(m.chrom).to_vec(),
                    m.key,
                    m.loc,
                    batch.raw(i).to_vec(),
                ));
            }
            max = next_batch_size(&batch, max).min(5);
            if !more {
                return out;
            }
        }
    }

    const BED: &str = "track name=t\n# c\n\nchr1\t1\t5\nchr1\t2\t6\r\na\tx\ty\n\n#x\nchr2\t0\t3\nchr2\t4\t9\nchr3\t7\t8";
    const VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n\nchr1\t6\t.\tA\tC\t.\t.\t.\n#mid\nchr2\t7\t.\tA\tC\t.\t.\t.\r\nchr2\t8\t.\tA\tC\t.\t.\t.\n";

    /// Ranges found with any buffer size, including one smaller than a
    /// line, read back as the whole-file reader reads the file: the same
    /// records, line numbers (blank, comment and unterminated lines
    /// included) and skipped-record diagnostics.
    #[test]
    fn ranges_read_as_the_whole_file() {
        for (text, format, name) in [(BED, Format::Bed, "t.bed"), (VCF, Format::Vcf, "t.vcf")] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join(name);
            std::fs::write(&path, text).unwrap();
            let whole = open_reader(path.to_str().unwrap(), &opts()).unwrap();
            let header = whole.header().clone();
            let want = read_all(whole);
            let want_diag = {
                let mut r = open_reader(path.to_str().unwrap(), &opts()).unwrap();
                let mut c = ChromTable::first_appearance();
                let mut b = RecordBatch::default();
                while r.read_batch(&mut c, &mut b, 100).unwrap() {}
                r.diagnostics().clone()
            };
            let reference = scan(&path, format, SCAN_BUF).unwrap().unwrap();
            for buf in [1, 2, 7, 16, 64, SCAN_BUF] {
                assert_eq!(
                    scan(&path, format, buf).unwrap().unwrap(),
                    reference,
                    "{buf}"
                );
            }
            // (BED has a line `a x y`, which no reader accepts: it is a range
            // of its own, and `--skip-invalid` skips it there.)
            assert_eq!(reference.len(), if format == Format::Bed { 4 } else { 2 });
            // The ranges tile the records: each starts where the last ends.
            for w in reference.windows(2) {
                assert_eq!(w[0].bytes.end, w[1].bytes.start);
            }
            assert_eq!(reference.last().unwrap().bytes.end, text.len() as u64);
            let got =
                read_all(open_range_reader(&path, &reference, format, &header, &opts()).unwrap());
            assert_eq!(got, want, "{name}");
            // Any one range gives its chromosome's records alone.
            for r in &reference {
                let one = read_all(
                    open_range_reader(&path, std::slice::from_ref(r), format, &header, &opts())
                        .unwrap(),
                );
                let expect: Vec<_> = want.iter().filter(|x| x.0 == r.name).cloned().collect();
                assert_eq!(one, expect, "{name} {:?}", r.name);
            }
            // The diagnostics, summed over ranges.
            let mut r = open_range_reader(&path, &reference, format, &header, &opts()).unwrap();
            let mut c = ChromTable::first_appearance();
            let mut b = RecordBatch::default();
            while r.read_batch(&mut c, &mut b, 2).unwrap() {}
            assert_eq!(r.diagnostics().skipped, want_diag.skipped, "{name}");
            assert_eq!(r.diagnostics().examples, want_diag.examples, "{name}");
        }
    }

    #[test]
    fn non_contiguous_chromosomes_have_no_ranges() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bed");
        std::fs::write(&path, "chr1\t1\t2\nchr2\t1\t2\nchr1\t3\t4\n").unwrap();
        for buf in [1, 5, SCAN_BUF] {
            assert_eq!(scan(&path, Format::Bed, buf).unwrap(), None);
        }
        // Comments and blanks between the runs do not split a chromosome.
        std::fs::write(&path, "chr1\t1\t2\n\n# c\nchr1\t3\t4\n").unwrap();
        assert_eq!(scan(&path, Format::Bed, 4).unwrap().unwrap().len(), 1);
        // A record line with no tab: not splittable.
        std::fs::write(&path, "chr1\t1\t2\nchr1 3 4\n").unwrap();
        assert_eq!(scan(&path, Format::Bed, 4).unwrap(), None);
        // Many more chromosomes than lines/2 after 10,000 lines: not worth it.
        let many: String = (0..20_000).map(|i| format!("c{i:06}\t1\t2\n")).collect();
        std::fs::write(&path, &many).unwrap();
        assert_eq!(scan(&path, Format::Bed, SCAN_BUF).unwrap(), None);
        // Two lines a chromosome is fine.
        let two: String = (0..20_000)
            .map(|i| format!("c{:06}\t{}\t2\n", i / 2, i % 2))
            .collect();
        std::fs::write(&path, &two).unwrap();
        assert_eq!(
            scan(&path, Format::Bed, SCAN_BUF).unwrap().unwrap().len(),
            10_000
        );
        // No records, no ranges.
        std::fs::write(&path, "# only a comment\n\n").unwrap();
        assert_eq!(scan(&path, Format::Bed, 4).unwrap().unwrap(), vec![]);
    }
}
