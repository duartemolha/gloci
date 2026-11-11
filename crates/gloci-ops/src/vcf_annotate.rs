//! `gloci vcf annotate`: a streaming merge of a VCF/BCF with an annotation
//! file (a VCF/BCF, a tab-delimited file or a BED file), both sorted in one
//! chromosome order and never loaded whole. The record edits are in
//! `gloci_vcf::annotate`; this module reads the files and keeps the rows
//! that can still match.

use crate::vcf_split::vcf_rewrite_output;
use gloci_core::{ChromId, ChromTable, OrderChecker, OrderViolation};
use gloci_io::{
    BcfText, Format, Header, InputOptions, IoError, Loc, OrderError, OrderGuard, OutFormat,
    RecordBatch, RecordReader, RecordWriter, TextToBcf, VcfHeader, next_batch_size, open_input,
    sort_hint,
};
use gloci_vcf::annotate::{
    Annotator, PairRec, Rec, TabLayout, TabRow, info_get, pair_rec, pair_rows,
};
use std::cmp::Ordering;
use std::io::BufRead;

const COMMAND: &str = "vcf annotate";

/// Counters of one run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AnnotateSummary {
    pub records_in: u64,
    pub records_out: u64,
    /// Records that had a matching annotation.
    pub matched: u64,
    /// The most records held at once, of FILE and of the annotation file together.
    pub peak_buffered: u64,
}

/// The annotation file.
pub enum AnnotationSource {
    None,
    Vcf(Box<dyn RecordReader>),
    /// A tab or BED file, opened by the command.
    Tab(Box<dyn BufRead + Send>, String),
}

/// The output format (`requested`, else FILE's) and the header with the
/// annotator's changes.
pub fn vcf_annotate_output(
    input: &dyn RecordReader,
    requested: Option<OutFormat>,
    header_text: &[u8],
) -> Result<(OutFormat, Header), IoError> {
    let (to, _) = vcf_rewrite_output(COMMAND, input, requested)?;
    let header = VcfHeader::parse(header_text.to_vec()).map_err(|e| IoError::Detect {
        path: input.name().to_string(),
        msg: e,
    })?;
    Ok((to, Header::Vcf(header)))
}

/// Whether `path` is read as BED (0-based starts) by its name, as bcftools does.
pub fn is_bed_name(path: &str) -> bool {
    let p = path.to_ascii_lowercase();
    p.ends_with(".bed") || p.ends_with(".bed.gz") || p.ends_with(".bed.bgz")
}

/// Open a tab or BED file for line reading (plain, gzip or BGZF).
pub fn open_tab(
    path: &str,
    threads: std::num::NonZeroUsize,
) -> Result<Box<dyn BufRead + Send>, IoError> {
    let input = open_input(
        path,
        &InputOptions {
            threads,
            format: Some(Format::Bed),
        },
    )?;
    Ok(input.reader)
}

/// FILE's or the source's records, a batch at a time.
struct Stream {
    reader: Box<dyn RecordReader>,
    name: String,
    batch: RecordBatch,
    i: usize,
    more: bool,
    max: usize,
    guard: OrderGuard,
    bcf: Option<BcfText>,
}

impl Stream {
    fn new(reader: Box<dyn RecordReader>) -> Result<Stream, IoError> {
        let name = reader.name().to_string();
        let guard = OrderGuard::new(&name, reader.format());
        let bcf = BcfText::for_reader(reader.as_ref())?;
        Ok(Stream {
            reader,
            name,
            batch: RecordBatch::default(),
            i: 0,
            more: true,
            max: 1,
            guard,
            bcf,
        })
    }

    /// The `(chrom, key)` of the next record, reading ahead when needed.
    fn peek(&mut self, chroms: &mut ChromTable) -> Result<Option<(ChromId, u64)>, IoError> {
        while self.i >= self.batch.len() {
            if !self.more {
                return Ok(None);
            }
            self.batch.clear();
            self.i = 0;
            self.more = self.reader.read_batch(chroms, &mut self.batch, self.max)?;
            self.guard.check(chroms, &self.batch)?;
            self.max = next_batch_size(&self.batch, self.max);
        }
        let m = &self.batch.recs[self.i];
        Ok(Some((m.chrom, m.key)))
    }

    /// Take the next record as a VCF line (appended to `out`); its number.
    fn take(&mut self, out: &mut Vec<u8>) -> Result<u64, IoError> {
        let loc = self.batch.recs[self.i].loc;
        let raw = self.batch.raw(self.i);
        match self.bcf.as_mut() {
            Some(t) => out.extend_from_slice(t.render(raw, loc)?),
            None => out.extend_from_slice(raw),
        }
        self.i += 1;
        Ok(loc)
    }

    fn loc(&self, n: u64) -> Loc {
        if self.reader.format() == Format::Bcf {
            Loc::Record(n)
        } else {
            Loc::Line(n)
        }
    }
}

fn cmp_pos(chroms: &ChromTable, a: (ChromId, u64), b: (ChromId, u64)) -> Ordering {
    chroms.cmp(a.0, b.0).then(a.1.cmp(&b.1))
}

/// One tab or BED row, with its chromosome in the shared table.
struct Row {
    row: TabRow,
    chrom: ChromId,
    loc: u64,
}

/// Reads annotation rows in order, checking the order.
struct TabStream {
    reader: Box<dyn BufRead + Send>,
    name: String,
    layout: TabLayout,
    line: Vec<u8>,
    lineno: u64,
    checker: OrderChecker,
    prev: Vec<u8>,
    prev_loc: u64,
    pending: Option<Row>,
    done: bool,
}

fn shown(b: &[u8]) -> String {
    let s = String::from_utf8_lossy(b);
    match s.char_indices().nth(120) {
        Some((cut, _)) => format!("{}...", &s[..cut]),
        None => s.into_owned(),
    }
}

impl TabStream {
    fn read_row(&mut self, chroms: &mut ChromTable) -> Result<Option<Row>, IoError> {
        loop {
            self.line.clear();
            let n = self
                .reader
                .read_until(b'\n', &mut self.line)
                .map_err(|source| IoError::Read {
                    path: self.name.clone(),
                    source,
                })?;
            if n == 0 {
                self.done = true;
                return Ok(None);
            }
            self.lineno += 1;
            while matches!(self.line.last(), Some(b'\n' | b'\r')) {
                self.line.pop();
            }
            let row = TabRow::parse(&self.line, &self.layout).map_err(|msg| IoError::Parse {
                path: self.name.clone(),
                loc: Loc::Line(self.lineno),
                msg,
            })?;
            let Some(row) = row else { continue };
            let chrom = chroms.intern(row.chrom()).map_err(|e| IoError::Parse {
                path: self.name.clone(),
                loc: Loc::Line(self.lineno),
                msg: e.to_string(),
            })?;
            let key = row.start.max(0) as u64;
            if let Err(v) = self.checker.check(chroms, chrom, key) {
                let hint = sort_hint(&self.name, chroms.order_source());
                let _: OrderViolation = v;
                return Err(OrderError {
                    path: self.name.clone(),
                    prev_loc: Loc::Line(self.prev_loc),
                    prev: shown(&self.prev),
                    loc: Loc::Line(self.lineno),
                    rec: shown(&self.line),
                    hint,
                }
                .into());
            }
            self.prev.clear();
            self.prev.extend_from_slice(&self.line);
            self.prev_loc = self.lineno;
            return Ok(Some(Row {
                row,
                chrom,
                loc: self.lineno,
            }));
        }
    }

    fn peek(&mut self, chroms: &mut ChromTable) -> Result<Option<&Row>, IoError> {
        if self.pending.is_none() && !self.done {
            self.pending = self.read_row(chroms)?;
        }
        Ok(self.pending.as_ref())
    }

    fn pop(&mut self) -> Option<Row> {
        self.pending.take()
    }
}

/// bcftools' region reader over a tab file with one position per row
/// (`bcf_sr_regions_overlap`/`_next`), made to stream.
struct PosReader {
    stream: TabStream,
    /// The current region: the last row read (`tgts` in bcftools).
    cur: Option<Row>,
    cur_chrom: Option<ChromId>,
    prev_chrom: Option<ChromId>,
    /// Rows that may still match: bcftools' `alines`.
    buf: Vec<Row>,
}

impl PosReader {
    /// `bcf_sr_regions_next` within chromosome `c`.
    fn next(&mut self, chroms: &mut ChromTable, c: ChromId) -> Result<bool, IoError> {
        match self.stream.peek(chroms)? {
            Some(r) if r.chrom == c => {
                self.cur = self.stream.pop();
                self.cur_chrom = Some(c);
                Ok(true)
            }
            _ => {
                self.cur = None;
                self.cur_chrom = None;
                Ok(false)
            }
        }
    }

    /// `bcf_sr_regions_overlap`: 0 when the current row overlaps `[start, end]`.
    fn overlap(
        &mut self,
        chroms: &mut ChromTable,
        c: ChromId,
        start: i64,
        end: i64,
    ) -> Result<bool, IoError> {
        if self.prev_chrom != Some(c) {
            // Seek to the chromosome: skip the rows before it.
            while let Some(r) = self.stream.peek(chroms)? {
                if chroms.cmp(r.chrom, c) != Ordering::Less {
                    break;
                }
                self.stream.pop();
            }
            match self.stream.peek(chroms)? {
                Some(r) if r.chrom == c => {}
                _ => return Ok(false), // the chromosome is not in the file
            }
            self.cur = None;
            self.cur_chrom = Some(c);
            self.prev_chrom = Some(c);
        }
        if self.cur_chrom != Some(c) {
            return Ok(false);
        }
        while self.cur.as_ref().is_none_or(|r| r.row.end < start) {
            if !self.next(chroms, c)? {
                return Ok(false);
            }
        }
        Ok(self.cur.as_ref().is_some_and(|r| r.row.start <= end))
    }

    /// `buffer_annot_lines`.
    fn fill(
        &mut self,
        chroms: &mut ChromTable,
        c: ChromId,
        beg: i64,
        end: i64,
        with_alleles: bool,
    ) -> Result<(), IoError> {
        if self.buf.first().is_some_and(|r| r.chrom != c) {
            self.buf.clear();
        }
        self.buf.retain(|r| beg <= r.row.end);
        if !with_alleles && !self.buf.is_empty() {
            return Ok(());
        }
        while self.overlap(chroms, c, beg, end)? {
            if self.buf.len() + 1 >= 0xffff {
                break;
            }
            let r = self.cur.take().expect("an overlapping current row");
            self.buf.push(Row {
                row: r.row.clone(),
                chrom: r.chrom,
                loc: r.loc,
            });
            self.cur = Some(r);
            if !with_alleles {
                break;
            }
            if !self.next(chroms, c)? {
                break;
            }
        }
        Ok(())
    }
}

/// Interval rows (a BED file or `FROM,TO` columns): the rows that can still overlap.
struct RegionReader {
    stream: TabStream,
    /// Sorted by start, then longer first, then arrival.
    active: Vec<Row>,
}

impl RegionReader {
    fn fill(
        &mut self,
        chroms: &mut ChromTable,
        c: ChromId,
        beg: i64,
        end: i64,
    ) -> Result<(), IoError> {
        self.active.retain(|r| r.chrom == c && r.row.end >= beg);
        while let Some(r) = self.stream.peek(chroms)? {
            let before = chroms.cmp(r.chrom, c) == Ordering::Less;
            let within = r.chrom == c && r.row.start <= end;
            if !before && !within {
                break;
            }
            let r = self.stream.pop().expect("a peeked row");
            if before || r.row.end < beg {
                continue;
            }
            let at = self
                .active
                .iter()
                .rposition(|a| {
                    a.row.start < r.row.start
                        || (a.row.start == r.row.start && a.row.end >= r.row.end)
                })
                .map_or(0, |p| p + 1);
            self.active.insert(at, r);
        }
        Ok(())
    }
}

fn parse_err(name: &str, loc: Loc, msg: String) -> IoError {
    IoError::Parse {
        path: name.to_string(),
        loc,
        msg,
    }
}

/// The records of FILE at one position.
struct Group {
    chrom: ChromId,
    key: u64,
    /// `(line, record number)`.
    recs: Vec<(Vec<u8>, u64)>,
}

/// INFO/END of a record, for the pairing key of symbolic alleles.
fn info_end(r: &Rec) -> Option<i64> {
    match info_get(&r.info, b"END") {
        Some(Some(v)) => std::str::from_utf8(v).ok()?.parse().ok(),
        _ => None,
    }
}

fn to_pair(r: &Rec) -> PairRec {
    pair_rec(&r.reference, &r.alts(), info_end(r))
}

struct Out<'a> {
    w: &'a mut dyn RecordWriter,
    enc: Option<TextToBcf>,
    buf: Vec<u8>,
    bin: Vec<u8>,
}

impl Out<'_> {
    fn put(&mut self, line: &[u8], n: u64) -> Result<(), IoError> {
        let werr = |source| IoError::Write {
            path: "output".to_string(),
            source,
        };
        match self.enc.as_mut() {
            Some(e) => {
                self.bin.clear();
                e.encode(line, n, &mut self.bin)?;
                self.w.write_record(&self.bin).map_err(werr)
            }
            None => self.w.write_record(line).map_err(werr),
        }
    }
}

/// Annotate `input` and write the records to `w` (opened with the output of
/// [`vcf_annotate_output`]).
pub fn vcf_annotate(
    input: Box<dyn RecordReader>,
    source: AnnotationSource,
    chroms: &mut ChromTable,
    annotator: &mut Annotator,
    out_header: &VcfHeader,
    w: &mut dyn RecordWriter,
) -> Result<AnnotateSummary, IoError> {
    let name = input.name().to_string();
    let enc = if w.is_bcf() {
        Some(TextToBcf::new(&name, out_header)?)
    } else {
        None
    };
    let mut out = Out {
        w,
        enc,
        buf: Vec::new(),
        bin: Vec::new(),
    };
    let mut file = Stream::new(input)?;
    let mut sum = AnnotateSummary::default();
    let editing = annotator.edits_records();
    let mut rec = Rec::default();

    enum Src {
        None,
        Vcf(Box<Stream>),
        Pos(Box<PosReader>),
        Regions(Box<RegionReader>),
    }
    let mut src = match source {
        AnnotationSource::None => Src::None,
        AnnotationSource::Vcf(r) => Src::Vcf(Box::new(Stream::new(r)?)),
        AnnotationSource::Tab(reader, sname) => {
            let layout = annotator
                .layout()
                .cloned()
                .ok_or_else(|| IoError::Internal("no tab layout".into()))?;
            let is_regions = layout.end.is_some();
            let stream = TabStream {
                reader,
                name: sname,
                layout,
                line: Vec::new(),
                lineno: 0,
                checker: OrderChecker::new(),
                prev: Vec::new(),
                prev_loc: 0,
                pending: None,
                done: false,
            };
            if is_regions {
                Src::Regions(Box::new(RegionReader {
                    stream,
                    active: Vec::new(),
                }))
            } else {
                Src::Pos(Box::new(PosReader {
                    stream,
                    cur: None,
                    cur_chrom: None,
                    prev_chrom: None,
                    buf: Vec::new(),
                }))
            }
        }
    };
    let mut src_line = Vec::new();

    // Records at one position are handled together: the source decides the
    // order they are written in (bcftools' synced reader).
    loop {
        let Some((chrom, key)) = file.peek(chroms)? else {
            break;
        };
        let mut group = Group {
            chrom,
            key,
            recs: Vec::new(),
        };
        while let Some(p) = file.peek(chroms)? {
            if p != (group.chrom, group.key) {
                break;
            }
            let mut l = Vec::new();
            let n = file.take(&mut l)?;
            group.recs.push((l, n));
        }
        let mut peak = group.recs.len() as u64;
        // (record index, source record) rows in write order.
        let mut vcf_rows: Vec<(usize, Option<Rec>)> = Vec::new();
        match &mut src {
            Src::Vcf(s) => {
                loop {
                    match s.peek(chroms)? {
                        Some(p)
                            if cmp_pos(chroms, p, (group.chrom, group.key)) == Ordering::Less =>
                        {
                            src_line.clear();
                            s.take(&mut src_line)?;
                        }
                        _ => break,
                    }
                }
                let mut srcs: Vec<Rec> = Vec::new();
                let mut locs: Vec<u64> = Vec::new();
                while let Some(p) = s.peek(chroms)? {
                    if p != (group.chrom, group.key) {
                        break;
                    }
                    src_line.clear();
                    let n = s.take(&mut src_line)?;
                    let mut r = Rec::default();
                    r.load(&src_line)
                        .map_err(|m| parse_err(&s.name, s.loc(n), m))?;
                    srcs.push(r);
                    locs.push(n);
                }
                peak += srcs.len() as u64;
                if srcs.is_empty() {
                    vcf_rows.extend((0..group.recs.len()).map(|i| (i, None)));
                } else {
                    let mut frecs = Vec::with_capacity(group.recs.len());
                    let mut fpair = Vec::with_capacity(group.recs.len());
                    for (l, n) in &group.recs {
                        let mut r = Rec::default();
                        r.load(l)
                            .map_err(|m| parse_err(&file.name, file.loc(*n), m))?;
                        fpair.push(to_pair(&r));
                        frecs.push(r);
                    }
                    let spair: Vec<PairRec> = srcs.iter().map(to_pair).collect();
                    let mut srcs: Vec<Option<Rec>> = srcs.into_iter().map(Some).collect();
                    for (fi, si) in pair_rows(&fpair, &spair) {
                        if let Some(fi) = fi {
                            vcf_rows.push((fi, si.and_then(|s| srcs[s].take())));
                        }
                    }
                }
            }
            Src::Pos(_) | Src::Regions(_) | Src::None => {}
        }
        sum.peak_buffered = sum.peak_buffered.max(peak);

        let order: Vec<usize> = if matches!(src, Src::Vcf(_)) {
            vcf_rows.iter().map(|r| r.0).collect()
        } else {
            (0..group.recs.len()).collect()
        };
        let mut vcf_src: Vec<Option<Rec>> = vcf_rows.into_iter().map(|r| r.1).collect();
        for (k, &ri) in order.iter().enumerate() {
            let (text, n) = &group.recs[ri];
            sum.records_in += 1;
            let mut has_overlap = false;
            if editing {
                rec.load(text)
                    .map_err(|m| parse_err(&file.name, file.loc(*n), m))?;
                annotator.apply_remove(&mut rec);
                let fail = |m: String| parse_err(&file.name, file.loc(*n), m);
                match &mut src {
                    Src::None => {}
                    Src::Vcf(_) => {
                        if let Some(s) = vcf_src[k].take() {
                            has_overlap = true;
                            annotator.apply_vcf(&mut rec, &s).map_err(fail)?;
                        }
                    }
                    Src::Pos(p) => {
                        let pos0 = rec.position().map_err(fail)? as i64 - 1;
                        let end0 = pos0 + rec.reference.len().max(1) as i64 - 1;
                        p.fill(
                            chroms,
                            group.chrom,
                            pos0,
                            end0,
                            annotator.layout().is_some_and(|l| l.ref_.is_some()),
                        )?;
                        sum.peak_buffered = sum
                            .peak_buffered
                            .max(group.recs.len() as u64 + p.buf.len() as u64);
                        let mut hit: Vec<usize> = Vec::new();
                        for (i, r) in p.buf.iter().enumerate() {
                            if pos0 > r.row.end || end0 < r.row.start {
                                continue;
                            }
                            if annotator.tab_matches(&rec, pos0, &r.row) {
                                hit.push(i);
                            }
                        }
                        if !hit.is_empty() {
                            has_overlap = true;
                            let rows: Vec<&TabRow> = hit.iter().map(|&i| &p.buf[i].row).collect();
                            let loc = hit.first().map(|&i| p.buf[i].loc).unwrap_or(0);
                            let sname = p.stream.name.clone();
                            annotator
                                .apply_rows(&mut rec, &rows)
                                .map_err(|m| parse_err(&sname, Loc::Line(loc), m))?;
                        }
                    }
                    Src::Regions(rg) => {
                        let pos0 = rec.position().map_err(fail)? as i64 - 1;
                        let len = rec.reference.len().max(1) as i64;
                        let end0 = pos0 + len - 1;
                        rg.fill(chroms, group.chrom, pos0, end0)?;
                        sum.peak_buffered = sum
                            .peak_buffered
                            .max(group.recs.len() as u64 + rg.active.len() as u64);
                        let mo = annotator.min_overlap();
                        let mut rows: Vec<&TabRow> = Vec::new();
                        for r in &rg.active {
                            if r.row.start > end0 || r.row.end < pos0 {
                                continue;
                            }
                            let isec = (r.row.end.min(end0) - r.row.start.max(pos0) + 1) as f32;
                            let len_ann = (r.row.end - r.row.start + 1) as f32;
                            if mo.ann > 0.0 && mo.ann > isec / len_ann {
                                continue;
                            }
                            if mo.vcf > 0.0 && mo.vcf > isec / len as f32 {
                                continue;
                            }
                            rows.push(&r.row);
                        }
                        if !rows.is_empty() {
                            has_overlap = true;
                            let sname = rg.stream.name.clone();
                            annotator
                                .apply_rows(&mut rec, &rows)
                                .map_err(|m| parse_err(&sname, Loc::Line(0), m))?;
                        }
                    }
                }
                annotator.finish(&mut rec, has_overlap);
                if has_overlap {
                    sum.matched += 1;
                }
                if rec.dirty {
                    out.buf.clear();
                    rec.write(&mut out.buf);
                    let b = std::mem::take(&mut out.buf);
                    out.put(&b, *n)?;
                    out.buf = b;
                } else {
                    out.put(text, *n)?;
                }
            } else {
                out.put(text, *n)?;
            }
            sum.records_out += 1;
        }
    }
    // Both inputs are read to the end, so that every record is validated.
    match &mut src {
        Src::Vcf(s) => {
            while s.peek(chroms)?.is_some() {
                src_line.clear();
                s.take(&mut src_line)?;
            }
        }
        Src::Pos(p) => {
            while p.stream.peek(chroms)?.is_some() {
                p.stream.pop();
            }
        }
        Src::Regions(rg) => {
            while rg.stream.peek(chroms)?.is_some() {
                rg.stream.pop();
            }
        }
        Src::None => {}
    }
    if let Some(s) = file.reader.diagnostics().summary(&file.name) {
        log::warn!("{s}");
    }
    Ok(sum)
}
