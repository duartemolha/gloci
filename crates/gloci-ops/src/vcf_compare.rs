//! `gloci vcf compare`: genotype concordance of two sorted VCF/BCF files, the
//! counterpart of vcftools `--diff-site-discordance`, `--diff-indv-discordance`
//! and `--diff-discordance-matrix` (decision D163).
//!
//! Each file goes through the shared selection stage on its own thread and
//! hands owned records (a few hundred at a time, through a bounded channel) to
//! the walk on the calling thread. The walk pairs the records of the two files
//! (`gloci_vcf::compare`), checks each file's order, and counts concordance
//! over the samples both files have. Memory: the records of one position in
//! both files, the channel's few batches, and 24 bytes per shared sample of
//! each record held.

use crate::vcf_select::{SelectStage, SelectedInput, VcfSelect, rewrite_plan, run_selection};
use crossbeam_channel::{Receiver, Sender};
use gloci_core::{ChromId, ChromOrder, ChromTable, OrderChecker};
use gloci_io::{
    AtomicFile, BcfDicts, IoError, Loc, OrderError, StagedSet, UsageError, VcfHeader,
    chrom_table_lenient, sort_hint,
};
use gloci_vcf::compare::genotypes::alt_count;
use gloci_vcf::compare::{
    AlleleMap, DiscordanceMatrix, Gt, GtCompare, MatchRule, MatrixCell, Side, Site, SiteClass,
    SiteSource, Walk,
};
use gloci_vcf::{BcfSite, GtFields, RewritePlan, SampleFields, SiteView};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::ffi::OsString;
use std::io::{self, BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;

const COMMAND: &str = "vcf compare";

/// Which files to write; none requested means all three.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompareRequest {
    /// `PREFIX.diff.sites`.
    pub site_discordance: bool,
    /// `PREFIX.diff.indv`.
    pub sample_discordance: bool,
    /// `PREFIX.diff.discordance_matrix`.
    pub discordance_matrix: bool,
}

impl CompareRequest {
    fn resolved(&self) -> CompareRequest {
        if *self == CompareRequest::default() {
            CompareRequest {
                site_discordance: true,
                sample_discordance: true,
                discordance_matrix: true,
            }
        } else {
            self.clone()
        }
    }
}

#[derive(Clone, Debug)]
pub struct CompareOptions {
    pub select: VcfSelect,
    pub matching: MatchRule,
    /// Two columns: a sample's name in the second file, then its name in the
    /// first (vcftools `--diff-indv-map`).
    pub sample_map: Option<PathBuf>,
    pub request: CompareRequest,
    pub out_prefix: PathBuf,
    /// `-g`: the chromosome order, else the first file's header, else
    /// `chrom_order`.
    pub genome: Option<String>,
    pub chrom_order: ChromOrder,
}

/// What a run did: rows of the walk by class, and the files written.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompareSummary {
    pub both: u64,
    pub first: u64,
    pub second: u64,
    pub overlap: u64,
    pub files: Vec<PathBuf>,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn write_err(path: &str, source: io::Error) -> IoError {
    IoError::Write {
        path: path.to_string(),
        source,
    }
}

/// Check the options before any input is read.
pub fn check_options(o: &CompareOptions) -> Result<(), IoError> {
    o.select.validate()?;
    if o.out_prefix.as_os_str().is_empty() {
        return Err(usage(
            "--out-prefix is empty; give a prefix such as --out-prefix results/cmp",
        ));
    }
    Ok(())
}

/// `--sample-map`: (name in the second file, name in the first).
fn read_sample_map(path: &std::path::Path) -> Result<HashMap<Vec<u8>, Vec<u8>>, IoError> {
    let name = path.display().to_string();
    let text = std::fs::read(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut map = HashMap::new();
    for (i, line) in text.split(|&b| b == b'\n').enumerate() {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() || line.starts_with(b"#") {
            continue;
        }
        let cols: Vec<&[u8]> = line
            .split(|b| b.is_ascii_whitespace())
            .filter(|c| !c.is_empty())
            .collect();
        if cols.len() != 2 {
            return Err(IoError::Parse {
                path: name,
                loc: Loc::Line(i as u64 + 1),
                msg: format!(
                    "expected two columns (a sample's name in the second file, then in the first), found {}",
                    cols.len()
                ),
            });
        }
        if map.insert(cols[0].to_vec(), cols[1].to_vec()).is_some() {
            return Err(IoError::Parse {
                path: name,
                loc: Loc::Line(i as u64 + 1),
                msg: format!(
                    "{} is named twice; give each sample of the second file one new name",
                    String::from_utf8_lossy(cols[0])
                ),
            });
        }
    }
    Ok(map)
}

/// One file's selected samples: (index in the record, name after mapping).
struct Selected {
    list: Vec<(usize, Vec<u8>)>,
}

/// The names of `list` that `own` has, by raw or by mapped name, as raw names.
fn own_names(list: &[String], raw: &[Vec<u8>], mapped: &[Vec<u8>]) -> Vec<String> {
    let want: HashSet<&[u8]> = list.iter().map(|s| s.as_bytes()).collect();
    raw.iter()
        .zip(mapped)
        .filter(|(r, m)| want.contains(r.as_slice()) || want.contains(m.as_slice()))
        .map(|(r, _)| String::from_utf8_lossy(r).into_owned())
        .collect()
}

/// The selection of one file: `--samples` and `--exclude-samples` apply by
/// name to each file that has the sample.
fn side_select(sel: &VcfSelect, raw: &[Vec<u8>], mapped: &[Vec<u8>]) -> VcfSelect {
    let mut s = sel.clone();
    s.include_samples = sel
        .include_samples
        .as_ref()
        .map(|l| own_names(l, raw, mapped));
    s.exclude_samples = own_names(&sel.exclude_samples, raw, mapped);
    s
}

/// A listed name that neither file has is a usage error.
fn check_listed(
    flag: &str,
    list: &[String],
    names: &[&[Vec<u8>]],
    files: [&str; 2],
) -> Result<(), IoError> {
    let have: HashSet<&[u8]> = names
        .iter()
        .flat_map(|n| n.iter().map(|s| s.as_slice()))
        .collect();
    for n in list {
        if !have.contains(n.as_bytes()) {
            return Err(usage(format!(
                "{flag} names {n}, which is in neither {} nor {}; check the spelling against the #CHROM line of each file (with --sample-map, either name of a renamed sample works)",
                files[0], files[1]
            )));
        }
    }
    Ok(())
}

fn selected(h: &VcfHeader, sel: &VcfSelect, mapped: &[Vec<u8>]) -> Result<Selected, IoError> {
    let plan = rewrite_plan(h, sel, &RewritePlan::default())?;
    let keep = plan.selection.as_ref().map(|s| &s.keep);
    Ok(Selected {
        list: (0..h.samples.len())
            .filter(|&i| keep.is_none_or(|k| k[i]))
            .map(|i| (i, mapped[i].clone()))
            .collect(),
    })
}

/// The optional FORMAT fields the selection reads (the comparison reads GT only).
fn decode_fields(sel: &VcfSelect) -> GtFields {
    if sel.expr.is_some() {
        return GtFields::ALL;
    }
    let (m, s) = (&sel.mask, &sel.sites);
    GtFields {
        dp: m.min_dp.is_some()
            || m.max_dp.is_some()
            || s.min_mean_dp.is_some()
            || s.max_mean_dp.is_some(),
        gq: m.min_gq.is_some(),
        ft: m.mask_filtered,
    }
}

enum Msg {
    Batch(Vec<Site>),
    End(Result<(), IoError>),
}

/// Read one file through the selection stage and send its records.
fn produce(
    input: SelectedInput,
    sel: &VcfSelect,
    shared: &[usize],
    threads: NonZeroUsize,
    tx: &Sender<Msg>,
) -> Result<(), IoError> {
    let (reader, mut chroms) = input.into_parts();
    let base = RewritePlan::default();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), sel, &base, !shared.is_empty())?;
    stage.set_threads(threads);
    if !shared.is_empty() {
        stage.share_decoded(decode_fields(sel));
    }
    let name = stage.name().to_string();
    let header = stage.header().clone();
    let is_bcf = stage.input_is_bcf();
    let dicts = if is_bcf {
        match &header.bcf {
            Some(d) => d.clone(),
            None => BcfDicts::from_header(&header.raw).map_err(|e| IoError::Detect {
                path: name.clone(),
                msg: e,
            })?,
        }
    } else {
        BcfDicts::default()
    };
    let mut bcf_site = BcfSite::default();
    let cap = (262_144 / (shared.len() + 16)).clamp(1, 1024);
    let mut batch: Vec<Site> = Vec::with_capacity(cap);
    let mut fields = SampleFields::default();
    run_selection(reader, &mut chroms, &mut stage, |hit| {
        let parse_err = |msg: String| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(hit.rec_no),
            msg,
        };
        let view: SiteView = if is_bcf {
            bcf_site.view(hit.raw, &dicts, hit.meta.span)
        } else {
            SiteView::from_text(hit.raw, hit.meta.span)
        }
        .map_err(parse_err)?;
        let n_alts = alt_count(view.alts) as u32;
        let mut gts = Vec::with_capacity(shared.len());
        if !shared.is_empty() {
            let src = hit.genotypes().ok_or_else(|| {
                IoError::Internal("the selection kept no genotypes for a record".into())
            })?;
            for &i in shared {
                src.sample(i, &mut fields).map_err(parse_err)?;
                let gt = match &fields.gt {
                    None => Gt::Missing,
                    Some(g) => Gt::from_alleles(&g.alleles),
                };
                if let Some(m) = gt.max_allele()
                    && m > n_alts
                {
                    return Err(parse_err(format!(
                        "a genotype names allele {m}, but the record has {n_alts} ALT allele(s); fix the GT or the ALT column"
                    )));
                }
                gts.push(gt);
            }
        }
        batch.push(Site {
            chrom_name: view.chrom.into(),
            chrom: ChromId(0),
            pos: view.pos,
            ref_: view.ref_.to_ascii_uppercase(),
            alts: view.alts.to_ascii_uppercase(),
            rec_no: hit.rec_no,
            gts,
        });
        if batch.len() >= cap {
            let full = std::mem::replace(&mut batch, Vec::with_capacity(cap));
            tx.send(Msg::Batch(full))
                .map_err(|_| IoError::Internal("the comparison stopped reading".into()))?;
        }
        Ok(())
    })?;
    if !batch.is_empty() {
        tx.send(Msg::Batch(batch))
            .map_err(|_| IoError::Internal("the comparison stopped reading".into()))?;
    }
    Ok(())
}

/// One file's records as the walk takes them: ids resolved and order checked.
struct Stream {
    rx: Receiver<Msg>,
    name: String,
    buf: VecDeque<Site>,
    done: bool,
    checker: OrderChecker,
    prev: Option<(Box<[u8]>, u64)>,
    last: Option<(Box<[u8]>, ChromId)>,
}

struct Sources {
    streams: [Stream; 2],
    chroms: ChromTable,
}

impl Stream {
    fn next(&mut self, chroms: &mut ChromTable) -> Result<Option<Site>, IoError> {
        loop {
            if let Some(mut s) = self.buf.pop_front() {
                self.resolve(chroms, &mut s)?;
                return Ok(Some(s));
            }
            if self.done {
                return Ok(None);
            }
            match self.rx.recv() {
                Ok(Msg::Batch(v)) => self.buf = v.into(),
                Ok(Msg::End(r)) => {
                    self.done = true;
                    r?;
                }
                Err(_) => {
                    return Err(IoError::Internal(format!(
                        "the reader of {} stopped before the end of the file",
                        self.name
                    )));
                }
            }
        }
    }

    fn resolve(&mut self, chroms: &mut ChromTable, s: &mut Site) -> Result<(), IoError> {
        let id = match &self.last {
            Some((n, id)) if *n == s.chrom_name => *id,
            _ => {
                let id = chroms.intern(&s.chrom_name).map_err(|e| IoError::Parse {
                    path: self.name.clone(),
                    loc: Loc::Record(s.rec_no),
                    msg: format!(
                        "chromosome {} is not in the order given by -g ({e}); add it to the genome file or drop -g",
                        String::from_utf8_lossy(&s.chrom_name)
                    ),
                })?;
                if chroms.has_undeclared() {
                    for u in chroms.take_undeclared() {
                        log::warn!(
                            "{}: contig {} is not declared in the header that sets the chromosome order; it is ordered after the declared ones, as first seen",
                            self.name,
                            String::from_utf8_lossy(chroms.name(u))
                        );
                    }
                }
                self.last = Some((s.chrom_name.clone(), id));
                id
            }
        };
        s.chrom = id;
        if self.checker.check(chroms, id, s.pos).is_err() {
            let (pc, pp) = self.prev.clone().unwrap_or_default();
            let show = |c: &[u8], p: u64| format!("{}:{}", String::from_utf8_lossy(c), p);
            return Err(OrderError {
                path: self.name.clone(),
                prev_loc: Loc::Record(s.rec_no.saturating_sub(1)),
                prev: show(&pc, pp),
                loc: Loc::Record(s.rec_no),
                rec: show(&s.chrom_name, s.pos),
                hint: format!(
                    "{}; both files must be sorted in one chromosome order",
                    sort_hint(&self.name, chroms.order_source())
                ),
            }
            .into());
        }
        self.prev = Some((s.chrom_name.clone(), s.pos));
        Ok(())
    }
}

impl SiteSource for Sources {
    type Err = IoError;
    fn next(&mut self, side: Side) -> Result<Option<Site>, IoError> {
        let i = usize::from(side == Side::Second);
        self.streams[i].next(&mut self.chroms)
    }
    fn cmp_chrom(&self, a: ChromId, b: ChromId) -> Ordering {
        self.chroms.cmp(a, b)
    }
}

fn dest_of(prefix: &std::path::Path, ext: &str) -> PathBuf {
    let mut name: OsString = prefix.to_path_buf().into_os_string();
    name.push(".");
    name.push(ext);
    PathBuf::from(name)
}

struct Out {
    w: BufWriter<AtomicFile>,
    dest: PathBuf,
}

impl Out {
    fn create(prefix: &std::path::Path, ext: &str) -> Result<Out, IoError> {
        let dest = dest_of(prefix, ext);
        let file = AtomicFile::create(&dest).map_err(|e| IoError::Create {
            path: dest.display().to_string(),
            source: e,
        })?;
        Ok(Out {
            w: BufWriter::new(file),
            dest,
        })
    }

    fn put(&mut self, bytes: &[u8]) -> Result<(), IoError> {
        self.w
            .write_all(bytes)
            .map_err(|e| write_err(&self.dest.display().to_string(), e))
    }

    fn stage(self, set: &mut StagedSet) -> Result<(), IoError> {
        let dest = self.dest.display().to_string();
        let file = self
            .w
            .into_inner()
            .map_err(|e| write_err(&dest, e.into_error()))?;
        set.push(file, self.dest);
        Ok(())
    }
}

/// `x` as C++ prints a double: six significant digits, `nan` for 0/0.
fn ratio(out: &mut Vec<u8>, discord: u64, common: u64) {
    if common == 0 {
        out.extend_from_slice(b"nan");
    } else {
        gloci_core::cfmt::write_g6(out, discord as f64 / common as f64);
    }
}

fn push_u64(out: &mut Vec<u8>, v: u64) {
    gloci_core::cfmt::write_u64(out, v);
}

/// Compare the genotypes of the samples `a` and `b` share, by site, by
/// sample and as a matrix, and write the files `o.request` names.
pub fn vcf_compare(
    a: SelectedInput,
    b: SelectedInput,
    o: &CompareOptions,
    threads: NonZeroUsize,
) -> Result<CompareSummary, IoError> {
    check_options(o)?;
    let request = o.request.resolved();
    let ha = crate::vcf_select::vcf_header(COMMAND, a.reader())?;
    let hb = crate::vcf_select::vcf_header(COMMAND, b.reader())?;
    let (name_a, name_b) = (a.reader().name().to_string(), b.reader().name().to_string());
    let map = match &o.sample_map {
        Some(p) => read_sample_map(p)?,
        None => HashMap::new(),
    };
    let mapped_b: Vec<Vec<u8>> = hb
        .samples
        .iter()
        .map(|s| map.get(s).cloned().unwrap_or_else(|| s.clone()))
        .collect();
    let unused = map.keys().filter(|k| !hb.samples.contains(k)).count();
    if unused > 0 {
        log::warn!(
            "--sample-map names {unused} sample(s) that {name_b} does not have; they are ignored"
        );
    }
    for flag_list in [
        (
            "--samples",
            o.select.include_samples.as_deref().unwrap_or(&[]),
        ),
        ("--exclude-samples", o.select.exclude_samples.as_slice()),
    ] {
        check_listed(
            flag_list.0,
            flag_list.1,
            &[&ha.samples, &hb.samples, &mapped_b],
            [&name_a, &name_b],
        )?;
    }
    let sel_a = side_select(&o.select, &ha.samples, &ha.samples);
    let sel_b = side_select(&o.select, &hb.samples, &mapped_b);
    let chosen_a = selected(&ha, &sel_a, &ha.samples)?;
    let chosen_b = selected(&hb, &sel_b, &mapped_b)?;
    // The union by name, in byte order; the samples in both files are compared.
    let mut union: BTreeMap<Vec<u8>, (Option<usize>, Option<usize>)> = BTreeMap::new();
    for (i, n) in &chosen_a.list {
        if union.insert(n.clone(), (Some(*i), None)).is_some() {
            return Err(usage(format!(
                "{name_a} has the sample {} twice; sample names must be distinct",
                String::from_utf8_lossy(n)
            )));
        }
    }
    for (i, n) in &chosen_b.list {
        let e = union.entry(n.clone()).or_insert((None, None));
        if e.1.is_some() {
            return Err(usage(format!(
                "{name_b} has the sample {} twice{}; sample names must be distinct",
                String::from_utf8_lossy(n),
                if o.sample_map.is_some() {
                    " after --sample-map"
                } else {
                    ""
                }
            )));
        }
        e.1 = Some(*i);
    }
    let names: Vec<&Vec<u8>> = union.keys().collect();
    let shared_pos: Vec<usize> = union
        .values()
        .enumerate()
        .filter(|(_, (x, y))| x.is_some() && y.is_some())
        .map(|(k, _)| k)
        .collect();
    let shared_a: Vec<usize> = shared_pos
        .iter()
        .map(|&k| union[names[k]].0.unwrap())
        .collect();
    let shared_b: Vec<usize> = shared_pos
        .iter()
        .map(|&k| union[names[k]].1.unwrap())
        .collect();
    let n_shared = shared_pos.len();

    let mut chroms = chrom_table_lenient(
        o.genome.as_deref(),
        &[a.reader(), b.reader()],
        o.chrom_order,
    )?;
    chroms.take_undeclared();
    let per_side =
        NonZeroUsize::new(threads.get().div_ceil(2)).expect("a thread count is at least 1");
    let mut staged = StagedSet::default();
    let mut sites_out = if request.site_discordance {
        let mut w = Out::create(&o.out_prefix, "diff.sites")?;
        w.put(b"CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n")?;
        Some(w)
    } else {
        None
    };

    let mut sum = CompareSummary::default();
    let mut common = vec![0u64; n_shared];
    let mut discord = vec![0u64; n_shared];
    let mut matrix = DiscordanceMatrix::default();
    let rule = o.matching;

    let walked = std::thread::scope(|scope| -> Result<(), IoError> {
        let (tx_a, rx_a) = crossbeam_channel::bounded::<Msg>(4);
        let (tx_b, rx_b) = crossbeam_channel::bounded::<Msg>(4);
        let (sel_a, sel_b) = (&sel_a, &sel_b);
        let (sa, sb) = (&shared_a, &shared_b);
        scope.spawn(move || {
            let r = produce(a, sel_a, sa, per_side, &tx_a);
            let _ = tx_a.send(Msg::End(r));
        });
        scope.spawn(move || {
            let r = produce(b, sel_b, sb, per_side, &tx_b);
            let _ = tx_b.send(Msg::End(r));
        });
        let stream = |rx: Receiver<Msg>, name: &str| Stream {
            rx,
            name: name.to_string(),
            buf: VecDeque::new(),
            done: false,
            checker: OrderChecker::new(),
            prev: None,
            last: None,
        };
        let mut walk = Walk::new(
            Sources {
                streams: [stream(rx_a, &name_a), stream(rx_b, &name_b)],
                chroms,
            },
            rule,
        );
        let mut map = AlleleMap::default();
        let mut line = Vec::new();
        while let Some(row) = walk.next_row()? {
            let (first, second) = (row.first.as_ref(), row.second.as_ref());
            match row.class {
                SiteClass::First => sum.first += 1,
                SiteClass::Second => sum.second += 1,
                SiteClass::Overlap => sum.overlap += 1,
                SiteClass::Both => sum.both += 1,
            }
            let one = |s: &Site, letter: &[u8], w: &mut Option<Out>, line: &mut Vec<u8>| {
                if let Some(w) = w {
                    line.clear();
                    line.extend_from_slice(&s.chrom_name);
                    line.push(b'\t');
                    push_u64(line, s.pos);
                    line.push(b'\t');
                    line.extend_from_slice(letter);
                    line.extend_from_slice(b"\t0\t0\t0\tnan\n");
                    return w.put(line);
                }
                Ok(())
            };
            match (row.class, first, second) {
                (SiteClass::First, Some(s), _) => one(s, b"1", &mut sites_out, &mut line)?,
                (SiteClass::Second, _, Some(s)) => one(s, b"2", &mut sites_out, &mut line)?,
                (SiteClass::Both, Some(x), Some(y)) => {
                    map.build(x, y, rule);
                    let (mut n_common, mut n_discord) = (0u64, 0u64);
                    for k in 0..n_shared {
                        let (gx, gy) = (&x.gts[k], &y.gts[k]);
                        match map.compare(gx, gy) {
                            GtCompare::Concordant => {
                                n_common += 1;
                                common[k] += 1;
                            }
                            GtCompare::Discordant => {
                                n_common += 1;
                                n_discord += 1;
                                common[k] += 1;
                                discord[k] += 1;
                            }
                            _ => {}
                        }
                        if map.biallelic_same
                            && let (MatrixCell::Cell(c1), MatrixCell::Cell(c2)) =
                                (AlleleMap::cell(gx), AlleleMap::cell(gy))
                        {
                            matrix.add(c1, c2);
                        }
                    }
                    if let Some(w) = sites_out.as_mut() {
                        line.clear();
                        line.extend_from_slice(&x.chrom_name);
                        line.push(b'\t');
                        push_u64(&mut line, x.pos);
                        line.extend_from_slice(b"\tB\t");
                        line.push(if map.matching { b'1' } else { b'0' });
                        line.push(b'\t');
                        push_u64(&mut line, n_common);
                        line.push(b'\t');
                        push_u64(&mut line, n_discord);
                        line.push(b'\t');
                        ratio(&mut line, n_discord, n_common);
                        line.push(b'\n');
                        w.put(&line)?;
                    }
                }
                _ => {}
            }
        }
        Ok(())
    });
    walked?;
    if let Some(w) = sites_out {
        w.stage(&mut staged)?;
    }
    if request.sample_discordance {
        let mut w = Out::create(&o.out_prefix, "diff.indv")?;
        let mut text = b"INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n".to_vec();
        let mut shared_k = 0usize;
        for (k, name) in names.iter().enumerate() {
            let (c, d) = if shared_pos.get(shared_k) == Some(&k) {
                shared_k += 1;
                (common[shared_k - 1], discord[shared_k - 1])
            } else {
                (0, 0)
            };
            text.extend_from_slice(name);
            text.push(b'\t');
            push_u64(&mut text, c);
            text.push(b'\t');
            push_u64(&mut text, d);
            text.push(b'\t');
            ratio(&mut text, d, c);
            text.push(b'\n');
        }
        w.put(&text)?;
        w.stage(&mut staged)?;
    }
    if request.discordance_matrix {
        let mut w = Out::create(&o.out_prefix, "diff.discordance_matrix")?;
        const G: [&str; 4] = ["0/0", "0/1", "1/1", "./."];
        let mut text = String::from("-");
        for g in G {
            text.push_str(&format!("\tN_{g}_file1"));
        }
        text.push('\n');
        for (r, g) in G.iter().enumerate() {
            text.push_str(&format!("N_{g}_file2"));
            for c in 0..4 {
                text.push_str(&format!("\t{}", matrix.cells[r][c]));
            }
            text.push('\n');
        }
        w.put(text.as_bytes())?;
        w.stage(&mut staged)?;
    }
    sum.files = staged.publish()?;
    Ok(sum)
}
