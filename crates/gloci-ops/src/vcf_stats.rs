//! `gloci vcf stats`: every requested statistic in one pass over the records
//! that pass the shared selection stage (decision D72).
//!
//! Each statistic writes its own temporary file next to its destination.
//! Only when the whole input has been read and every statistic has finished
//! are the files renamed into place, so an error publishes nothing; the
//! outputs of a previous run that are replaced are put back if publishing
//! fails part-way.

use crate::vcf_select::{SelectStage, SelectedInput, VcfSelect, run_selection};
use gloci_io::{AtomicFile, BcfDicts, IoError, Loc, StagedSet, UsageError, VcfHeader};
use gloci_vcf::{
    BcfSite, FilterSummary, Freq, GenoDepth, GtFields, Hardy, Heterozygosity, IndelHist,
    MissingSample, MissingSite, RewritePlan, SampleDepth, Singletons, SiteDepth, SiteMeanDepth,
    SiteQuality, SiteView, SnpDensity, Stat, StatSink, StatWrite, StatsEngine, Summary, TsTv,
};
use std::cell::RefCell;
use std::ffi::OsString;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::rc::Rc;

const COMMAND: &str = "vcf stats";

/// Which statistics to compute; each names the file it writes, as vcftools
/// does (`PREFIX.frq`, `PREFIX.idepth`, ...).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StatsRequest {
    pub freq: bool,
    pub counts: bool,
    pub depth: bool,
    pub site_depth: bool,
    pub site_mean_depth: bool,
    pub geno_depth: bool,
    pub site_quality: bool,
    pub missing_site: bool,
    pub missing_indv: bool,
    pub het: bool,
    pub hardy: bool,
    pub singletons: bool,
    pub filter_summary: bool,
    pub tstv_summary: bool,
    pub indel_hist: bool,
    /// Bin size of `--snp-density`.
    pub snp_density: Option<u64>,
    /// `PREFIX.summary.txt`.
    pub summary: bool,
    /// `PREFIX.summary.json`.
    pub json: bool,
}

impl StatsRequest {
    /// True when no statistic is requested.
    pub fn is_empty(&self) -> bool {
        *self == StatsRequest::default()
    }

    /// True when a requested statistic reads genotypes. Only quality, the
    /// FILTER summary, Ts/Tv, the indel histogram and SNP density do not.
    pub fn needs_genotypes(&self) -> bool {
        let site_only = StatsRequest {
            site_quality: self.site_quality,
            filter_summary: self.filter_summary,
            tstv_summary: self.tstv_summary,
            indel_hist: self.indel_hist,
            snp_density: self.snp_density,
            ..StatsRequest::default()
        };
        *self != site_only
    }
}

/// Options for `vcf stats`.
#[derive(Clone, Debug)]
pub struct StatsOptions {
    pub select: VcfSelect,
    /// Threads the selection runs on (`None` = 1); see `SelectStage::set_threads`.
    pub threads: Option<NonZeroUsize>,
    pub request: StatsRequest,
    /// Output files are `PREFIX.<extension>`.
    pub out_prefix: PathBuf,
}

/// What a `vcf stats` run did.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct StatsSummary {
    pub records_in: u64,
    pub records_selected: u64,
    /// The files written, in publication order.
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

/// Check the request before any input is read.
pub fn check_request(o: &StatsOptions) -> Result<(), IoError> {
    o.select.validate()?;
    if o.request.is_empty() {
        return Err(usage(
            "no statistic requested; add at least one of --freq, --counts, --depth, \
             --site-depth, --site-mean-depth, --geno-depth, --site-quality, --missing-site, \
             --missing-indv, --het, --hardy, --singletons, --filter-summary, --tstv-summary, \
             --indel-hist, --snp-density BIN, --summary or --json",
        ));
    }
    if o.request.snp_density == Some(0) {
        return Err(usage(
            "--snp-density needs a bin size above 0; give the bin width in bases, such as \
             --snp-density 10000",
        ));
    }
    if o.out_prefix.as_os_str().is_empty() {
        return Err(usage(
            "--out-prefix is empty; give a prefix such as --out-prefix results/cohort",
        ));
    }
    Ok(())
}

/// The statistics a request asks for, over the selected samples
/// (`samples`: index in the record and name).
fn build_stats(r: &StatsRequest, samples: &[(usize, String)]) -> Vec<Box<dyn Stat>> {
    let mut v: Vec<Box<dyn Stat>> = Vec::new();
    if r.freq {
        v.push(Box::new(Freq::frequencies()));
    }
    if r.counts {
        v.push(Box::new(Freq::counts()));
    }
    if r.depth {
        v.push(Box::new(SampleDepth::new(samples.to_vec())));
    }
    if r.site_depth {
        v.push(Box::new(SiteDepth::new()));
    }
    if r.site_mean_depth {
        v.push(Box::new(SiteMeanDepth::new()));
    }
    if r.geno_depth {
        v.push(Box::new(GenoDepth::new(samples.to_vec())));
    }
    if r.site_quality {
        v.push(Box::new(SiteQuality::new()));
    }
    if r.missing_site {
        v.push(Box::new(MissingSite::new()));
    }
    if r.missing_indv {
        v.push(Box::new(MissingSample::new(samples.to_vec())));
    }
    if r.het {
        v.push(Box::new(Heterozygosity::new(samples.to_vec())));
    }
    if r.hardy {
        v.push(Box::new(Hardy::new()));
    }
    if r.singletons {
        v.push(Box::new(Singletons::new(samples.to_vec())));
    }
    if r.filter_summary {
        v.push(Box::new(FilterSummary::new()));
    }
    if r.tstv_summary {
        v.push(Box::new(TsTv::new()));
    }
    if r.indel_hist {
        v.push(Box::new(IndelHist::new()));
    }
    if let Some(bin) = r.snp_density {
        v.push(Box::new(SnpDensity::new(bin)));
    }
    if r.summary || r.json {
        v.push(Box::new(Summary::new(
            samples.len() as u64,
            r.summary,
            r.json,
        )));
    }
    v
}

struct PrefixSink {
    prefix: PathBuf,
    staged: Rc<RefCell<StagedSet>>,
}

struct StagedFile {
    file: AtomicFile,
    dest: PathBuf,
    staged: Rc<RefCell<StagedSet>>,
}

impl Write for StagedFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.file.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

impl StatWrite for StagedFile {
    fn commit(mut self: Box<Self>) -> io::Result<()> {
        self.file.flush()?;
        let s = *self;
        s.staged.borrow_mut().push(s.file, s.dest);
        Ok(())
    }
}

impl StatSink for PrefixSink {
    fn create(&mut self, ext: &str) -> io::Result<Box<dyn StatWrite>> {
        let mut name: OsString = self.prefix.clone().into_os_string();
        name.push(".");
        name.push(ext);
        let dest = PathBuf::from(name);
        let file = AtomicFile::create(&dest).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!(
                    "cannot create {}: {e}; check that the directory exists and is writable",
                    dest.display()
                ),
            )
        })?;
        Ok(Box::new(StagedFile {
            file,
            dest,
            staged: Rc::clone(&self.staged),
        }))
    }
}

/// The optional FORMAT fields `vcf stats` decodes: the selection's, plus DP
/// for a depth statistic or a summary.
fn decode_fields(o: &StatsOptions) -> GtFields {
    let r = &o.request;
    let mut f = selection_fields(&o.select);
    f.dp |= r.depth || r.site_depth || r.site_mean_depth || r.geno_depth || r.summary || r.json;
    f
}

/// The selected samples as (index in the record, name).
fn selected_samples(stage: &SelectStage<'_>) -> Vec<(usize, String)> {
    let all = &stage.header().samples;
    let keep = stage.plan().selection.as_ref().map(|s| &s.keep);
    all.iter()
        .enumerate()
        .filter(|(i, _)| keep.is_none_or(|k| k[*i]))
        .map(|(i, s)| (i, String::from_utf8_lossy(s).into_owned()))
        .collect()
}

/// The optional FORMAT fields the selection reads: DP for a DP mask or a
/// mean-depth site filter, GQ for a GQ mask, FT for an FT mask. An expression
/// may read any, so it takes all.
pub(crate) fn selection_fields(sel: &VcfSelect) -> GtFields {
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

/// Compute the requested statistics of the records of `input` and write
/// `PREFIX.<extension>` for each.
pub fn vcf_stats(input: SelectedInput, o: &StatsOptions) -> Result<StatsSummary, IoError> {
    check_request(o)?;
    run_stats(
        input,
        &StatsRun {
            select: &o.select,
            prefix: &o.out_prefix,
            fields: decode_fields(o),
            need_genotypes: o.request.needs_genotypes(),
            threads: o.threads.unwrap_or(NonZeroUsize::MIN),
            command: COMMAND,
        },
        &|samples, _| Ok(build_stats(&o.request, samples)),
    )
}

/// What builds a run's statistics from the selected samples (index in the
/// record and name) and the header.
pub(crate) type BuildStats<'a> =
    dyn Fn(&[(usize, String)], &VcfHeader) -> Result<Vec<Box<dyn Stat>>, IoError> + 'a;

/// How [`run_stats`] runs: the selection, where the files go, what to decode.
pub(crate) struct StatsRun<'a> {
    pub select: &'a VcfSelect,
    /// Each statistic writes `prefix.<extension>`.
    pub prefix: &'a Path,
    /// The optional FORMAT fields to decode.
    pub fields: GtFields,
    /// True when a statistic needs genotypes or [`gloci_vcf::SiteStats`].
    pub need_genotypes: bool,
    /// Threads the selection runs on.
    pub threads: NonZeroUsize,
    /// The command, for messages.
    pub command: &'a str,
}

/// Run the selection stage on `input` and feed every selected record to the
/// statistics `build` makes from the selected samples (index in the record
/// and name) and the header, in one pass, as `run` says. Each statistic
/// writes `PREFIX.<extension>`; the files are published together once every
/// statistic has finished.
pub(crate) fn run_stats(
    input: SelectedInput,
    run: &StatsRun<'_>,
    build: &BuildStats<'_>,
) -> Result<StatsSummary, IoError> {
    let &StatsRun {
        select,
        prefix,
        fields,
        need_genotypes,
        threads,
        command,
    } = run;
    let (reader, mut chroms) = input.into_parts();
    let base = RewritePlan::default();
    // The selection computes the site statistics (under the masking rules)
    // whenever a statistic needs them.
    let mut stage = SelectStage::new(command, reader.as_ref(), select, &base, need_genotypes)?;
    stage.set_threads(threads);
    stage.decode_fields(fields);
    let name = stage.name().to_string();
    let samples = selected_samples(&stage);
    let header = stage.header().clone();
    let mut engine = StatsEngine::new(build(&samples, &header)?);
    debug_assert_eq!(engine.needs_genotypes(), need_genotypes);
    // On several threads the workers keep the per-sample counters of --depth,
    // --missing-indv and --het; the records then carry no samples unless a
    // statistic still reads them.
    let indices: Vec<usize> = samples.iter().map(|(i, _)| *i).collect();
    let tallied = stage.enable_tally(engine.tally_wants(), &indices);
    if tallied {
        engine.use_tally();
    }
    if engine.reads_genotypes() {
        stage.share_decoded(fields);
    }

    let staged = Rc::new(RefCell::new(StagedSet::default()));
    let mut sink = PrefixSink {
        prefix: prefix.to_path_buf(),
        staged: Rc::clone(&staged),
    };
    let dest = prefix.display().to_string();
    engine.start(&mut sink).map_err(|e| write_err(&dest, e))?;

    let dicts = if stage.input_is_bcf() {
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
    let is_bcf = stage.input_is_bcf();
    let mut bcf_site = BcfSite::default();
    let mut selected = 0u64;
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        selected += 1;
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
        if !need_genotypes {
            return engine.site(&view, None, None).map_err(parse_err);
        }
        // Every selected sample was decoded by the selection stage.
        engine
            .site(&view, hit.genotypes(), hit.stats)
            .map_err(parse_err)?;
        if !hit.tally_called().is_empty() {
            engine.called(hit.tally_called());
        }
        Ok(())
    })?;
    if tallied {
        engine.absorb(stage.tally());
    }
    engine.finish().map_err(|e| write_err(&dest, e))?;
    drop(engine);
    drop(sink);
    let staged = Rc::try_unwrap(staged)
        .map_err(|_| IoError::Internal("the staged outputs are still shared".into()))?
        .into_inner();
    let files = staged.publish()?;
    Ok(StatsSummary {
        records_in,
        records_selected: selected,
        files,
    })
}
