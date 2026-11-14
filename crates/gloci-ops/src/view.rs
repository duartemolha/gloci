//! `gloci view`: stream a file unchanged or re-compressed, convert VCF and
//! BCF to each other, write only the header or only the records, and
//! restrict the output to regions (§7.4).
//!
//! A region query returns the records whose footprint (per `--vcf-span`)
//! overlaps a region under §4.3, each once, in file order (D9, D10). It
//! reads through an index when the file has one, and otherwise scans the
//! whole file with a warning, checking the order of every record, so an
//! unsorted file is an error, never a wrong answer. With `assume_sorted`,
//! an uncompressed BED file without an index is binary-searched instead
//! (D11): fast, but only the chromosome block it reads is checked, so the
//! answer is right only if the file really is sorted.

use crate::bedsearch::{BlockScan, Located, append_matching, locate_blocks};
use gloci_core::{
    ChromId, ChromOrder, ChromTable, Genome, Interval, OrderSource, Region, RegionSet,
};
use gloci_io::{
    AlnSource, Compression, Diagnostics, Format, Header, HeaderMode, IndexInfo, IndexedFile,
    InvalidPolicy, IoError, Opened, OrderGuard, OutFormat, OutputOptions, OutputSpec,
    ReaderOptions, RecordBatch, RecordMeta, RecordReader, UsageError, WriterOptions,
    alignment_query_span, aln_ref, check_family, chrom_table_lenient, convert_variants, find_index,
    index_reference_names, next_batch_size, open_headerless_writer, open_output, open_reader,
    open_reader_or_alignments, open_record_writer, open_region_reader, refuse_binary_to_tty,
    resolve_genome,
};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::IsTerminal;
use std::path::Path;

#[derive(Debug, Clone, Default)]
pub struct ViewOptions {
    /// Write the header and no records.
    pub header_only: bool,
    /// Write the records and no header.
    pub no_header: bool,
    /// Region strings (D8: 1-based, closed).
    pub regions: Vec<String>,
    /// A BED file of regions (0-based, half-open).
    pub regions_file: Option<String>,
    /// A point on a region's edge overlaps it (§4.3).
    pub inclusive_points: bool,
    /// `-g`: a built-in genome id or a file giving the chromosome names and
    /// the order binary search assumes.
    pub genome: Option<String>,
    /// The order binary search assumes without `-g` (D32).
    pub chrom_order: ChromOrder,
    /// Binary-search an uncompressed BED file without an index instead of
    /// scanning it: trusts that the file is sorted (D11).
    pub assume_sorted: bool,
    /// The input need not be in order (a relation with `--unsorted` or
    /// `--sort`): never read it through an index and do not check its order,
    /// so the records come out in input order.
    pub unordered: bool,
    /// Write the number of records `view` would write, as `N\n`, instead of the records (D133).
    pub count: bool,
}

impl ViewOptions {
    pub fn has_regions(&self) -> bool {
        !self.regions.is_empty() || self.regions_file.is_some()
    }
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Open `path` for viewing, with the chromosome table to read it with.
/// Without regions every record is read, in any order. With regions only
/// the records that overlap one are read, each once, in file order, and
/// the order of the input is checked.
pub fn open_view(
    path: &str,
    opts: &ViewOptions,
    ro: &ReaderOptions,
) -> Result<(Box<dyn RecordReader>, ChromTable), IoError> {
    open_view_with(open_reader(path, ro)?, path, opts, ro)
}

/// [`open_view`] for a file already opened as `reader`.
fn open_view_with(
    reader: Box<dyn RecordReader>,
    path: &str,
    opts: &ViewOptions,
    ro: &ReaderOptions,
) -> Result<(Box<dyn RecordReader>, ChromTable), IoError> {
    let genome = match &opts.genome {
        Some(arg) => Some(resolve_genome(arg)?),
        None => None,
    };
    let chroms = match &genome {
        Some((g, source)) => ChromTable::closed(g, source.clone()),
        None => chrom_table_lenient(None, &[reader.as_ref()], opts.chrom_order)?,
    };
    if !opts.has_regions() {
        return Ok((reader, chroms));
    }
    let genome = genome.as_ref().map(|(g, source)| (source, g));
    let selection = select(path, reader, opts, ro, genome, &chroms)?;
    Ok((Box::new(selection), chroms))
}

/// The records of the stream `reader` that overlap `opts`' regions, each
/// once, in input order, as [`open_view`] selects them from a file without
/// an index; `reader` itself when there are no regions.
pub(crate) fn select_stream(
    reader: Box<dyn RecordReader>,
    opts: &ViewOptions,
    ro: &ReaderOptions,
    chroms: &ChromTable,
) -> Result<Box<dyn RecordReader>, IoError> {
    if !opts.has_regions() {
        return Ok(reader);
    }
    Ok(Box::new(select("-", reader, opts, ro, None, chroms)?))
}

/// Regions as `(chromosome, interval)` pairs.
pub type WantedRegions = Vec<(Vec<u8>, Interval)>;

/// Restrict a relation's A to the records that overlap `opts`' regions: as
/// [`open_view`] selects them from the file at `path` (read through its
/// index when it has one), `reader` being that file opened. Also returns
/// the regions as `(chromosome, interval)` pairs in new-name space, minus
/// those `view` warns are on a chromosome the file lacks.
pub fn select_a(
    path: &str,
    reader: Box<dyn RecordReader>,
    opts: &ViewOptions,
    ro: &ReaderOptions,
) -> Result<(Box<dyn RecordReader>, WantedRegions), IoError> {
    if !opts.has_regions() {
        return Ok((reader, Vec::new()));
    }
    // Only binary search reads the table, and a relation never asks for it.
    let chroms = ChromTable::first_appearance();
    let sel = select(path, reader, opts, ro, None, &chroms)?;
    let wanted = sel.wanted.clone();
    Ok((Box::new(sel), wanted))
}

/// `chromosome "1" not found in x.vcf.gz (did you mean "chr1"?)`
fn not_found(chrom: &[u8], place: &str, has: impl Fn(&[u8]) -> bool) -> String {
    let alt = match chrom.strip_prefix(b"chr") {
        Some(rest) if !rest.is_empty() => rest.to_vec(),
        _ => [b"chr", chrom].concat(),
    };
    let hint = if has(&alt) {
        format!(" (did you mean \"{}\"?)", String::from_utf8_lossy(&alt))
    } else {
        String::new()
    };
    format!(
        "chromosome \"{}\" not found in {place}{hint}",
        String::from_utf8_lossy(chrom)
    )
}

/// The regions of a BED file: 0-based, half-open, as written.
pub fn read_regions_file(
    path: &str,
    ro: &ReaderOptions,
) -> Result<Vec<(Vec<u8>, Interval)>, IoError> {
    let mut out = Vec::new();
    for_each_region(path, ro, &mut |c, iv| {
        out.push((c.to_vec(), iv));
        Ok(())
    })?;
    Ok(out)
}

/// What [`for_each_region`] calls for each region: name and span.
pub type RegionSink<'a> = dyn FnMut(&[u8], Interval) -> Result<(), IoError> + 'a;

/// [`read_regions_file`] without keeping the regions: each is passed to `f`
/// as it is read.
pub fn for_each_region(path: &str, ro: &ReaderOptions, f: &mut RegionSink) -> Result<(), IoError> {
    // Regions are written with the chromosome names the records have after
    // `--chrom-map`, so the table is not applied to the regions file.
    let o = ReaderOptions {
        format: Some(Format::Bed),
        invalid: InvalidPolicy::Fail,
        chrom_map: None,
        ..ro.clone()
    };
    let mut r = open_reader(path, &o)?;
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, max)?;
        for m in &batch.recs {
            f(chroms.name(m.chrom), m.span)?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            return Ok(());
        }
    }
}

/// How to fix a file that must be scanned for want of an index.
fn index_hint(name: &str, format: Format, compression: Compression) -> String {
    if compression == Compression::Bgzf {
        return format!("create one with: gloci index {name}");
    }
    let out = if name.ends_with(".gz") {
        name.to_string()
    } else {
        format!("{name}.gz")
    };
    let mut hint = format!("sort and index it with: gloci sort {name} -o {out} --write-index");
    if compression == Compression::None && matches!(format, Format::Bed | Format::BedGraph) {
        hint.push_str("; or, if it is sorted, pass --assume-sorted to binary-search it");
    }
    hint
}

/// The requested regions, resolved against what is known of the file.
struct Wanted {
    /// Regions on chromosomes the file can have, in new-name space.
    wanted: Vec<(Vec<u8>, Interval)>,
    /// `wanted` by chromosome.
    sets: Vec<(Vec<u8>, RegionSet)>,
    /// The file's chromosome list is complete (an index or genome gave it).
    complete: bool,
}

/// Parse `opts`' regions and drop (with a warning) those on chromosomes the
/// file does not have, when its chromosomes are known up front.
fn resolve_regions(
    name: &str,
    header: &Header,
    opts: &ViewOptions,
    ro: &ReaderOptions,
    genome: Option<(&OrderSource, &Genome)>,
    file: Option<&Path>,
    index: Option<&IndexInfo>,
) -> Result<Wanted, IoError> {
    // The chromosomes the file can have, when that is known before reading
    // it, and what to name in warnings about the others.
    let mut known: HashSet<Vec<u8>> = HashSet::new();
    let mut known_in: Option<String> = None;
    // The list is complete when an index or a genome gives it; a VCF header
    // may omit contigs that records use, which htslib accepts.
    let mut complete = genome.is_some();
    if let (Some(p), Some(info)) = (file, index) {
        // The index lists the file's own names; regions use the renamed ones.
        let map = ro.chrom_map.as_deref();
        known.extend(
            index_reference_names(p, info)?
                .into_iter()
                .map(|n| map.map_or(n.clone(), |m| m.rename(&n).to_vec())),
        );
        known_in = Some(name.to_string());
        complete = true;
    }
    if let Header::Alignment(h) = &header {
        // The `@SQ` lines list every reference a record can name.
        known.extend(h.refs.iter().map(|(n, _)| n.clone()));
        known_in = Some(name.to_string());
        complete = true;
    }
    if let Header::Vcf(h) = &header {
        let before = known.len();
        match &h.bcf {
            Some(d) => known.extend(d.contigs.iter().filter(|c| !c.is_empty()).cloned()),
            None => known.extend(h.contigs.iter().map(|(c, _)| c.clone())),
        }
        if known.len() > before {
            known_in = Some(name.to_string());
        }
    }
    if let Some((source, g)) = genome {
        known.extend((0..g.len()).map(|i| g.name(i).to_vec()));
        known_in.get_or_insert_with(|| match source {
            OrderSource::BuiltinGenome(id) => format!("the built-in genome {id}"),
            OrderSource::GenomeFile(p) => p.display().to_string(),
            _ => "the genome".to_string(),
        });
    }

    let mut wanted: Vec<(Vec<u8>, Interval)> = Vec::new();
    for s in &opts.regions {
        let r = Region::parse(s, |n| known.contains(n)).map_err(|e| usage(e.to_string()))?;
        // An open end means the whole rest of the chromosome, whatever
        // length a header or genome declares.
        let iv = r.interval(None);
        wanted.push((r.chrom, iv));
    }
    if let Some(f) = &opts.regions_file {
        wanted.extend(read_regions_file(f, ro)?);
    }
    if let Some(place) = known_in.as_ref().filter(|_| complete) {
        let mut warned = HashSet::new();
        wanted.retain(|(c, _)| {
            if known.contains(c) {
                return true;
            }
            if warned.insert(c.clone()) {
                log::warn!("{}", not_found(c, place, |n| known.contains(n)));
            }
            false
        });
    }
    let mut by_chrom: BTreeMap<Vec<u8>, Vec<Interval>> = BTreeMap::new();
    for (c, iv) in &wanted {
        by_chrom.entry(c.clone()).or_default().push(*iv);
    }
    let sets: Vec<(Vec<u8>, RegionSet)> = by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, RegionSet::new(ivs)))
        .collect();
    Ok(Wanted {
        wanted,
        sets,
        complete,
    })
}

/// Choose how to read the records of `head` (an open reader of `path`)
/// that overlap the requested regions.
fn select(
    path: &str,
    head: Box<dyn RecordReader>,
    opts: &ViewOptions,
    ro: &ReaderOptions,
    genome: Option<(&OrderSource, &Genome)>,
    chroms: &ChromTable,
) -> Result<Selection, IoError> {
    let name = head.name().to_string();
    let (format, compression) = (head.format(), head.compression());
    let header = head.header().clone();
    let file = Some(Path::new(path)).filter(|p| path != "-" && p.is_file());
    let index = match file.and_then(find_index).filter(|_| !opts.unordered) {
        Some(info) if compression == Compression::Bgzf => {
            info.warn_if_stale(&name);
            Some(info)
        }
        Some(info) => {
            log::warn!(
                "ignoring the index {}: {name} is not BGZF-compressed",
                info.path.display()
            );
            None
        }
        None => None,
    };

    let Wanted {
        wanted,
        sets,
        complete,
    } = resolve_regions(&name, &header, opts, ro, genome, file, index.as_ref())?;

    let mut warn_unseen = false;
    let src = if wanted.is_empty() {
        Src::Empty
    } else if let (Some(p), Some(_)) = (file, &index) {
        Src::Indexed(open_region_reader(p, &wanted, ro, opts.inclusive_points)?)
    } else if let Some(p) = file.filter(|_| {
        opts.assume_sorted
            && compression == Compression::None
            && matches!(format, Format::Bed | Format::BedGraph)
    }) {
        match locate_blocks(p, &name, sets.clone(), chroms, ro.chrom_map.as_deref())? {
            Located::Blocks(blocks) => {
                Src::Blocks(Box::new(BlockScan::new(p, &name, format, ro, blocks)))
            }
            Located::Missing(missing) => {
                log::warn!(
                    "could not locate {} by binary search; scanning the whole file (sort it with gloci sort, or index it)",
                    missing.join(", ")
                );
                warn_unseen = !complete;
                Src::Scan(head)
            }
        }
    } else {
        if gloci_io::memio::is_mem_path(path) {
            // Nothing to index: it is read once, as it is.
            log::debug!("{name} is in memory; scanning it whole for the regions");
        } else if !gloci_io::is_regular_file(path) {
            // Stdin, a pipe or a descriptor path (`<(...)`, /dev/stdin): no index can
            // exist for it, so there is nothing to advise.
            log::debug!("{name} is not a regular file; scanning it whole for the regions");
        } else if !opts.unordered {
            log::warn!(
                "no index for {name}; scanning the whole file ({})",
                index_hint(&name, format, compression)
            );
        }
        warn_unseen = !complete;
        Src::Scan(head)
    };

    Ok(Selection {
        guard: OrderGuard::contiguous(&name, format),
        check_order: !opts.unordered,
        wanted: wanted.clone(),
        name,
        format,
        compression,
        header,
        inclusive: opts.inclusive_points,
        src,
        matcher: Matcher::new(sets),
        scratch: RecordBatch::default(),
        read_max: 1,
        pending: RecordBatch::default(),
        taken: 0,
        warn_unseen,
        done: false,
        no_diag: Diagnostics::default(),
    })
}

/// Where a selection's records come from.
enum Src {
    /// No region can match.
    Empty,
    /// The index reader, which already selects records.
    Indexed(Box<dyn RecordReader>),
    /// Every record of the file, filtered here.
    Scan(Box<dyn RecordReader>),
    /// The chromosome blocks binary search found.
    Blocks(Box<BlockScan>),
}

/// The regions by chromosome name, looked up once per run of records on
/// one chromosome.
struct Matcher {
    sets: Vec<(Vec<u8>, RegionSet)>,
    by_name: HashMap<Vec<u8>, usize>,
    last: Option<(ChromId, Option<usize>)>,
}

impl Matcher {
    fn new(sets: Vec<(Vec<u8>, RegionSet)>) -> Matcher {
        let by_name = sets
            .iter()
            .enumerate()
            .map(|(i, (c, _))| (c.clone(), i))
            .collect();
        Matcher {
            sets,
            by_name,
            last: None,
        }
    }

    fn matches(&mut self, chroms: &ChromTable, m: &RecordMeta, inclusive: bool) -> bool {
        let i = match self.last {
            Some((id, i)) if id == m.chrom => i,
            _ => {
                let i = self.by_name.get(chroms.name(m.chrom)).copied();
                self.last = Some((m.chrom, i));
                i
            }
        };
        i.is_some_and(|i| self.sets[i].1.matches(m.span, inclusive))
    }
}

/// The records of one input that overlap the requested regions.
struct Selection {
    name: String,
    format: Format,
    compression: Compression,
    header: Header,
    inclusive: bool,
    /// Index and scan paths: chromosomes contiguous, each sorted by key.
    guard: OrderGuard,
    /// Whether to run `guard` (not for an input that may be unordered).
    check_order: bool,
    /// The regions, for [`select_a`].
    wanted: Vec<(Vec<u8>, Interval)>,
    src: Src,
    matcher: Matcher,
    /// The source's records, read `read_max` at a time.
    scratch: RecordBatch,
    read_max: usize,
    /// Selected records, of which `taken` have been handed out.
    pending: RecordBatch,
    taken: usize,
    /// After a scan, warn about region chromosomes the file never showed
    /// (when its chromosomes were not known up front).
    warn_unseen: bool,
    /// The source is exhausted.
    done: bool,
    no_diag: Diagnostics,
}

impl Selection {
    /// Read the next batch of the source (sized by `next_batch_size`, not by
    /// the caller, so a run of non-matching records is read in large
    /// batches) and queue the records that overlap a region.
    fn fill(&mut self, chroms: &mut ChromTable) -> Result<(), IoError> {
        let (scratch, max) = (&mut self.scratch, self.read_max);
        scratch.clear();
        let more = match &mut self.src {
            Src::Empty => false,
            Src::Indexed(r) => {
                let more = r.read_batch(chroms, scratch, max)?;
                self.guard.check(chroms, scratch)?;
                append_matching(&mut self.pending, scratch, |_| true);
                more
            }
            Src::Scan(r) => {
                let more = r.read_batch(chroms, scratch, max)?;
                if self.check_order {
                    self.guard.check(chroms, scratch)?;
                }
                let (matcher, inclusive) = (&mut self.matcher, self.inclusive);
                append_matching(&mut self.pending, scratch, |m| {
                    matcher.matches(chroms, m, inclusive)
                });
                more
            }
            Src::Blocks(b) => b.read(chroms, scratch, &mut self.pending, max, self.inclusive)?,
        };
        self.read_max = next_batch_size(scratch, max);
        if !more {
            self.done = true;
            if self.warn_unseen {
                for (c, _) in &self.matcher.sets {
                    if chroms.get(c).is_none() {
                        log::warn!("{}", not_found(c, &self.name, |n| chroms.get(n).is_some()));
                    }
                }
            }
        }
        Ok(())
    }
}

impl RecordReader for Selection {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        self.format
    }
    fn compression(&self) -> Compression {
        self.compression
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        match &self.src {
            Src::Indexed(r) | Src::Scan(r) => r.diagnostics(),
            Src::Blocks(b) => &b.diag,
            Src::Empty => &self.no_diag,
        }
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        loop {
            let queued = self.pending.len() - self.taken;
            if queued > 0 {
                let n = queued.min(max);
                for i in self.taken..self.taken + n {
                    let off = batch.buf.len();
                    batch.buf.extend_from_slice(self.pending.raw(i));
                    batch.recs.push(RecordMeta {
                        off,
                        ..self.pending.recs[i]
                    });
                }
                self.taken += n;
                return Ok(self.taken < self.pending.len() || !self.done);
            }
            if self.done {
                return Ok(false);
            }
            self.pending.clear();
            self.taken = 0;
            self.fill(chroms)?;
        }
    }
}

/// View `path` into the output `out_spec`. With no regions the file is
/// streamed in file order; sortedness is not required.
pub fn run_view(
    path: &str,
    opts: &ViewOptions,
    ro: &ReaderOptions,
    out_spec: &OutputSpec,
    oo: &OutputOptions,
    wo: &WriterOptions,
) -> Result<(), IoError> {
    if opts.header_only && opts.no_header {
        return Err(usage("--header-only and --no-header cannot be combined"));
    }
    if opts.count {
        let clash = if opts.header_only {
            Some("--header-only")
        } else if oo.write_index {
            Some("--write-index")
        } else {
            None
        };
        if let Some(flag) = clash {
            return Err(usage(format!(
                "--count cannot be combined with {flag}; it writes a number, not records"
            )));
        }
    }
    let reference = ro
        .accept
        .alignments
        .as_ref()
        .and_then(|a| a.reference.as_deref());
    let reader = match open_reader_or_alignments(path, ro)? {
        Opened::Records(r) => {
            gloci_io::warn_unused_reference(reference, false);
            gloci_io::warn_unused_read_filters(ro.accept.alignments.as_ref(), false);
            r
        }
        Opened::Alignments(src) => {
            gloci_io::warn_unused_reference(reference, src.format() == Format::Cram);
            return view_alignments(src, path, opts, ro, out_spec, oo, wo);
        }
    };
    let (mut reader, mut chroms) = if opts.header_only {
        // No records are read, so regions do not matter.
        (reader, ChromTable::first_appearance())
    } else {
        open_view_with(reader, path, opts, ro)?
    };
    if opts.count {
        let mut batch = RecordBatch::default();
        let (mut n, mut max) = (0u64, 1);
        loop {
            batch.clear();
            let more = reader.read_batch(&mut chroms, &mut batch, max)?;
            n += batch.len() as u64;
            max = next_batch_size(&batch, max);
            if !more {
                break;
            }
        }
        write_count(n, out_spec, oo)?;
        if let Some(s) = reader.diagnostics().summary(reader.name()) {
            log::warn!("{s}");
        }
        return Ok(());
    }
    let format = reader.format();
    let family = OutFormat::family(format);
    let out_format = out_spec.format.unwrap_or(family);
    let is_variant = |f: OutFormat| matches!(f, OutFormat::Vcf | OutFormat::Bcf);
    let converting = is_variant(family) && is_variant(out_format) && out_format != family;
    if !converting {
        check_family(format, out_format)?;
    }
    if opts.no_header && out_format == OutFormat::Bcf {
        return Err(usage(
            "--no-header cannot be used with BCF output; BCF records need their header",
        ));
    }
    refuse_binary_to_tty(
        out_spec.path.is_none() && std::io::stdout().is_terminal(),
        out_format,
        out_spec.bgzf,
    )?;
    let dest = out_spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    // Errors of the output side name the destination.
    let named = |e: IoError| e.with_output_name(&dest);

    let mut wo = wo.clone();
    // View keeps whatever header the file has unless told otherwise.
    wo.write_header = !opts.no_header;
    let mut out = open_output(out_spec, out_format, oo)?;
    // An order error found while indexing names the input, which is the file that needs sorting.
    out.set_index_input(reader.name());

    if opts.header_only {
        // No records are read, so nothing is checked either.
        let mode = HeaderMode::Only;
        if converting {
            convert_variants(
                reader.as_mut(),
                &mut chroms,
                out,
                out_format,
                &wo,
                mode,
                None,
            )
            .map_err(named)?;
        } else {
            let header = reader.header().clone();
            open_record_writer(out, format, out_format, &header, &wo)
                .map_err(named)?
                .finish()
                .map_err(named)?;
        }
        return Ok(());
    }

    if converting {
        let mode = if opts.no_header {
            HeaderMode::None
        } else {
            HeaderMode::Full
        };
        // A region view checks order itself.
        convert_variants(
            reader.as_mut(),
            &mut chroms,
            out,
            out_format,
            &wo,
            mode,
            None,
        )
        .map_err(named)?;
    } else {
        let header = reader.header().clone();
        let mut w = if opts.no_header {
            open_headerless_writer(out, format, out_format)
        } else {
            open_record_writer(out, format, out_format, &header, &wo)
        }
        .map_err(named)?;
        let mut batch = RecordBatch::default();
        let mut max = 1;
        loop {
            batch.clear();
            let more = reader.read_batch(&mut chroms, &mut batch, max)?;
            for i in 0..batch.len() {
                w.write_record(batch.raw(i))
                    .map_err(|source| IoError::Write {
                        path: dest.clone(),
                        source,
                    })?;
            }
            max = next_batch_size(&batch, max);
            if !more {
                break;
            }
        }
        w.finish().map_err(named)?;
    }
    if let Some(s) = reader.diagnostics().summary(reader.name()) {
        log::warn!("{s}");
    }
    Ok(())
}

/// Write `n` and a newline to the output (`view --count`).
fn write_count(n: u64, out_spec: &OutputSpec, oo: &OutputOptions) -> Result<(), IoError> {
    use std::io::Write;
    let dest = out_spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let mut out = open_output(out_spec, OutFormat::Tsv, oo)?;
    writeln!(out, "{n}").map_err(|source| IoError::Write {
        path: dest.clone(),
        source,
    })?;
    out.finish().map_err(|e| e.with_output_name(&dest))
}

/// `view` of BAM, CRAM or SAM: every record is written, unmapped included (D116, D120), through
/// the read filters; with regions, those that overlap one, found through the BAM's index
/// when it has one and by scanning the file otherwise (a CRAM is scanned). Both select with samtools' span
/// (`alignment_query_span`) and give the same records.
fn view_alignments(
    mut src: Box<dyn AlnSource>,
    path: &str,
    opts: &ViewOptions,
    ro: &ReaderOptions,
    out_spec: &OutputSpec,
    oo: &OutputOptions,
    wo: &WriterOptions,
) -> Result<(), IoError> {
    let count = opts.count;
    let layout = src.layout();
    // A CRAM's records are BAM bytes, and its regions are found through its `.crai`.
    let input_format = src.format();
    let out_format = out_spec.format.unwrap_or_else(|| OutFormat::family(layout));
    if !count {
        check_family(layout, out_format)?;
    }
    if !count && opts.no_header && out_format == OutFormat::Bam {
        return Err(usage(
            "--no-header cannot be used with BAM output; BAM records need their header; \
write SAM instead (-O sam)",
        ));
    }
    if !count {
        refuse_binary_to_tty(
            out_spec.path.is_none() && std::io::stdout().is_terminal(),
            out_format,
            out_spec.bgzf,
        )?;
    }
    let dest = out_spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let header = src.header().clone();
    let name = src.name().to_string();
    let named = |e: IoError| e.with_output_name(&dest);

    // Regions, as the chromosomes the header lists (a name it lacks is warned about and
    // dropped); `by_ref[rid]` is the region set of reference `rid`.
    let hdr = Header::Alignment(header.clone());
    let mut by_ref: Vec<Option<RegionSet>> = Vec::new();
    let regions = opts.has_regions() && !opts.header_only;
    // The index reader already selects the records.
    let mut indexed = false;
    if regions {
        let Wanted { wanted, sets, .. } = resolve_regions(&name, &hdr, opts, ro, None, None, None)?;
        let file = Some(Path::new(path)).filter(|p| path != "-" && p.is_file());
        let index = file
            .filter(|_| matches!(input_format, Format::Bam | Format::Cram) && !wanted.is_empty())
            .and_then(|p| find_index(p).map(|i| (p, i)));
        if let Some((p, info)) = index {
            info.warn_if_stale(&name);
            let ix = IndexedFile::open(p, input_format, hdr.clone(), ro)?;
            src = ix.aln_regions_source(&wanted, opts.inclusive_points)?;
            indexed = true;
        } else {
            let mut by_name: HashMap<Vec<u8>, RegionSet> = sets.into_iter().collect();
            by_ref = header.refs.iter().map(|(n, _)| by_name.remove(n)).collect();
            if path != "-" && by_ref.iter().any(Option::is_some) {
                // Only a regular file can be indexed (not a FIFO such as `<(cat x.bam)`).
                let hint = match input_format {
                    Format::Bam | Format::Cram if !wanted.is_empty() && file.is_some() => {
                        format!("; create an index with: gloci index {name}")
                    }
                    _ => String::new(),
                };
                log::warn!("scanning the whole file {name} to find the regions{hint}");
            }
        }
    }

    let mut wo = wo.clone();
    wo.write_header = !opts.no_header;
    let mut w = if count {
        None
    } else {
        let mut out = open_output(out_spec, out_format, oo)?;
        // An order error found while indexing names the input, which is the file to sort.
        out.set_index_input(&name);
        Some(open_record_writer(out, layout, out_format, &hdr, &wo).map_err(named)?)
    };
    let mut counted = 0u64;
    // No region on any chromosome of the file: nothing to read.
    let filter = regions && !indexed;
    let nothing = filter && by_ref.iter().all(Option::is_none);
    if count && !regions {
        // Only the read filters decide what counts, and they read the flag and the MAPQ: a
        // CRAM decodes those data series alone, as flagstat does, and no other record field.
        src.for_each_flag_fields(&mut |_| {
            counted += 1;
            Ok(())
        })?;
    } else if !opts.header_only && !nothing {
        let mut buf = Vec::new();
        while src.next_record(&mut buf)?.is_some() {
            let keep = !filter
                || aln_ref(layout, &buf, &header)
                    .ok()
                    .and_then(|a| alignment_query_span(&a))
                    .is_some_and(|(rid, iv)| {
                        by_ref[rid]
                            .as_ref()
                            .is_some_and(|set| set.matches(iv, opts.inclusive_points))
                    });
            if keep && let Some(w) = w.as_mut() {
                w.write_alignment(&buf, src.last_bam()).map_err(|source| {
                    if source.kind() == std::io::ErrorKind::InvalidData {
                        // The record cannot be written in this format; not a write failure.
                        IoError::Detect {
                            path: name.clone(),
                            msg: format!("{name}: {source}"),
                        }
                    } else {
                        IoError::Write {
                            path: dest.clone(),
                            source,
                        }
                    }
                })?;
            } else if keep {
                counted += 1;
            }
            buf.clear();
        }
    }
    match w {
        Some(w) => w.finish().map_err(named)?,
        None => write_count(counted, out_spec, oo)?,
    }
    src.diagnostics().log_alignment_counts(&name);
    if let Some(s) = src.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    Ok(())
}
