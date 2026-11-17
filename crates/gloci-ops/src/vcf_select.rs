//! The selection stage shared by the `vcf` commands: sample selection,
//! genotype masking, site filters, the include/exclude expression, regions
//! and thinning (spec §7.8, decisions D50 and D71).
//!
//! `vcf filter` writes the records that pass; `vcf extract` tabulates them.
//! A command opens its input with [`open_selected`], builds a
//! [`SelectStage`] and drives it with [`run_selection`].

use crate::view::{ViewOptions, open_view, select_stream};
use gloci_core::{ChromOrder, ChromTable, VcfSpan};
use gloci_io::handoff::{Driving, Present, Sequenced};
use gloci_io::{
    BcfDicts, BcfText, Format, Header, IoError, Loc, OrderGuard, ReaderOptions, RecordBatch,
    RecordMeta, RecordReader, UsageError, VcfHeader, next_batch_size,
};
use gloci_vcf::{
    BcfSite, EvalScratch, Expr, GenotypeSource, GtFields, MaskRules, MaskedGenotypes,
    PackedGenotypes, QueryPass, RecordSlot, RewritePlan, SampleCache, SampleFields, SampleTally,
    SiteFilters, SiteStats, SiteView, TallyWants, Thinner, bcf_genotypes_with, parse_expr,
    rewrite_text, select_samples, site_stats, text_genotypes_with,
};
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// A `--include` or `--exclude` expression (decision D71).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FilterExpr {
    pub src: String,
    /// `--exclude`: keep the records the expression does not hold on.
    pub exclude: bool,
}

/// The selection options of a `vcf` command: which samples, which genotypes
/// are masked, which sites and regions pass.
#[derive(Clone, Debug, Default)]
pub struct VcfSelect {
    pub include_samples: Option<Vec<String>>,
    pub exclude_samples: Vec<String>,
    pub mask: MaskRules,
    pub sites: SiteFilters,
    pub expr: Option<FilterExpr>,
    pub regions: Vec<String>,
    pub regions_file: Option<String>,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Check the numeric thresholds, which needs no input: callers run it before
/// reading any list file, so a usage error is never reported as a file error.
pub fn validate_thresholds(
    min_maf: Option<f64>,
    max_maf: Option<f64>,
    max_missing: Option<f64>,
    min_alleles: Option<usize>,
    max_alleles: Option<usize>,
) -> Result<(), IoError> {
    let unit = [
        ("--min-maf", min_maf),
        ("--max-maf", max_maf),
        ("--max-missing", max_missing),
    ];
    for (flag, v) in unit {
        if let Some(v) = v
            && !(0.0..=1.0).contains(&v)
        {
            return Err(usage(format!(
                "{flag} {v} is not in [0, 1]; give a fraction between 0 and 1"
            )));
        }
    }
    if let (Some(lo), Some(hi)) = (min_alleles, max_alleles)
        && lo > hi
    {
        return Err(usage(format!(
            "--min-alleles {lo} is more than --max-alleles {hi}; no site could pass"
        )));
    }
    Ok(())
}

impl VcfSelect {
    /// The checks that need no input: thresholds within range.
    pub fn validate(&self) -> Result<(), IoError> {
        validate_thresholds(
            self.sites.min_maf,
            self.sites.max_maf,
            self.sites.max_missing,
            self.sites.min_alleles,
            self.sites.max_alleles,
        )
    }

    /// The region selection, as `view` does it.
    fn view_options(&self) -> ViewOptions {
        ViewOptions {
            header_only: false,
            no_header: false,
            regions: self.regions.clone(),
            regions_file: self.regions_file.clone(),
            inclusive_points: self.sites.inclusive_points,
            genome: None,
            chrom_order: ChromOrder::default(),
            assume_sorted: false,
            unordered: false,
            count: false,
        }
    }
}

/// `input`'s VCF header, or an error when it is not VCF/BCF. `command` names
/// the command for the message.
pub fn vcf_header(command: &str, input: &dyn RecordReader) -> Result<VcfHeader, IoError> {
    let format = input.format();
    if format != Format::Vcf && format != Format::Bcf {
        return Err(usage(format!(
            "{command} works on VCF and BCF input, but {} is {}; give it a VCF or BCF file",
            input.name(),
            format.name()
        )));
    }
    match input.header() {
        Header::Vcf(h) => Ok(h.clone()),
        _ => Err(IoError::Detect {
            path: input.name().to_string(),
            msg: "a VCF/BCF input has no VCF header".into(),
        }),
    }
}

/// The rewrite of each record: the sample selection (from the include and
/// exclude lists, else `base`'s) and the masking rules.
pub fn rewrite_plan(
    h: &VcfHeader,
    sel: &VcfSelect,
    base: &RewritePlan,
) -> Result<RewritePlan, IoError> {
    let mut plan = base.clone();
    if sel.include_samples.is_some() || !sel.exclude_samples.is_empty() {
        plan.selection = Some(select_samples(
            &h.samples,
            sel.include_samples.as_deref(),
            &sel.exclude_samples,
        )?);
    }
    plan.mask = sel.mask;
    Ok(plan)
}

/// Parse the expression of `sel` against `h`, if there is one.
pub fn parse_select_expr(h: &VcfHeader, sel: &VcfSelect) -> Result<Option<Expr>, IoError> {
    let Some(e) = &sel.expr else {
        return Ok(None);
    };
    parse_expr(&e.src, h).map(Some).map_err(IoError::from)
}

/// An input opened by [`open_selected`], its regions already selected.
pub struct SelectedInput {
    pub(crate) reader: Box<dyn RecordReader>,
    pub(crate) chroms: ChromTable,
}

impl SelectedInput {
    pub fn reader(&self) -> &dyn RecordReader {
        self.reader.as_ref()
    }

    /// The reader and the chromosome table, to hand to [`run_selection`].
    pub fn into_parts(self) -> (Box<dyn RecordReader>, ChromTable) {
        (self.reader, self.chroms)
    }
}

/// Open `path` for a `vcf` command. With `--region`/`--regions-file`, only the
/// records that overlap a region are read: through the index when the file
/// has one, else by scanning it (as `view` does).
pub fn open_selected(
    path: &str,
    ro: &ReaderOptions,
    sel: &VcfSelect,
) -> Result<SelectedInput, IoError> {
    let (reader, chroms) = open_view(path, &sel.view_options(), ro)?;
    Ok(SelectedInput { reader, chroms })
}

/// Apply the regions of `sel` to an already open stream.
pub fn select_regions(
    command: &str,
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    sel: &VcfSelect,
) -> Result<Box<dyn RecordReader>, IoError> {
    sel.validate()?;
    vcf_header(command, input.as_ref())?;
    let region_ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Ref,
        invalid: gloci_io::InvalidPolicy::Fail,
        threads: std::num::NonZero::new(1).unwrap(),
    };
    select_stream(input, &sel.view_options(), &region_ro, chroms)
}

/// A record that passed the selection.
pub struct Hit<'a> {
    /// The raw record: a VCF line, or a BCF record.
    pub raw: &'a [u8],
    /// 1-based record number within the stream read.
    pub rec_no: u64,
    pub meta: &'a RecordMeta,
    /// Site statistics, when the selection computed them.
    pub stats: Option<&'a SiteStats>,
    /// Per-genotype mask flags (empty when nothing is masked).
    pub masked: &'a [bool],
    genotypes: Option<&'a dyn GenotypeSource>,
    called: &'a [u64],
    /// With `SelectStage::filter_samples`: the header samples that pass a per-sample filter
    /// (D172); `None`: every selected sample.
    pub sample_pass: Option<&'a [bool]>,
}

impl Hit<'_> {
    /// With [`SelectStage::enable_tally`]: the bitset of selected samples
    /// that are fully called diploid at this site, when the site contributes
    /// to `--het` (else empty).
    pub fn tally_called(&self) -> &[u64] {
        self.called
    }

    /// The record's decoded samples, the masked ones reading as missing. `None`
    /// unless [`SelectStage::share_decoded`] was called and the selection
    /// computed site statistics for the record. Only the samples the
    /// selection keeps are decoded: reading another sample is an error (a
    /// bug in the caller), the same on one thread and on many.
    pub fn genotypes(&self) -> Option<&dyn GenotypeSource> {
        self.genotypes
    }
}

/// The records of a stage as VCF text with the selection's sample choice and
/// masking applied: a BCF record is rendered by the converter `view -O vcf`
/// uses, so VCF and BCF input give identical lines.
pub(crate) struct Lines {
    bcf: Option<BcfText>,
    plan: RewritePlan,
    rewritten: Vec<u8>,
}

impl Lines {
    pub(crate) fn for_stage(stage: &SelectStage<'_>) -> Result<Lines, IoError> {
        Ok(Lines {
            bcf: if stage.input_is_bcf() {
                Some(BcfText::for_header(stage.name(), stage.header())?)
            } else {
                None
            },
            plan: stage.plan().clone(),
            rewritten: Vec::new(),
        })
    }

    pub(crate) fn line<'b>(&'b mut self, name: &str, hit: &Hit<'b>) -> Result<&'b [u8], IoError> {
        let text: &[u8] = match self.bcf.as_mut() {
            Some(b) => b.render(hit.raw, hit.rec_no)?,
            None => hit.raw,
        };
        if self.plan.is_passthrough() {
            return Ok(text);
        }
        self.rewritten.clear();
        rewrite_text(text, &self.plan, hit.masked, hit.stats, &mut self.rewritten).map_err(
            |msg| IoError::Parse {
                path: name.to_string(),
                loc: Loc::Record(hit.rec_no),
                msg,
            },
        )?;
        Ok(&self.rewritten)
    }
}

/// Header indices of the samples a stage keeps, in header order.
pub(crate) fn kept_indices(stage: &SelectStage<'_>) -> Vec<usize> {
    match &stage.plan().selection {
        Some(sel) => (0..sel.keep.len()).filter(|&i| sel.keep[i]).collect(),
        None => (0..stage.header().samples.len()).collect(),
    }
}

/// Records one batch holds at most, so a batch decodes about this many genotypes.
const BATCH_GENOTYPES: usize = 524_288;

/// What the selection's static settings are: read-only, shared by the
/// workers of a parallel run.
struct Spec<'a> {
    sel: &'a VcfSelect,
    name: String,
    header: VcfHeader,
    in_format: Format,
    plan: RewritePlan,
    exclude: bool,
    dicts: BcfDicts,
    need_stats: bool,
    /// Statistics are computed only to decode the samples: no filter, mask
    /// rule or expression reads them.
    stats_only_for_decode: bool,
    /// Decode the samples for the consumers without computing statistics
    /// (see [`SelectStage::decode_only`]).
    decode_only: bool,
    share_decoded: bool,
    fields: GtFields,
    tally: Option<TallySpec>,
    /// Evaluate the expression with bcftools query's rule (D172).
    filter_samples: bool,
}

/// Per-sample counters the workers keep (see [`SelectStage::enable_tally`]).
struct TallySpec {
    wants: TallyWants,
    /// The selected samples' indices in the record.
    samples: Vec<usize>,
}

/// The per-thread working state of the selection.
struct Core {
    expr: Option<Expr>,
    scratch: EvalScratch,
    bcf_site: BcfSite,
    masked_buf: Vec<bool>,
    decoded: SampleCache,
    fields_buf: SampleFields,
    /// The per-sample passes of the last record (see [`Passed::per_sample`]).
    sample_buf: Vec<bool>,
}

/// A record the selection kept, before thinning.
struct Passed {
    stats: Option<SiteStats>,
    pos: u64,
    /// The record's per-sample passes are in `Core::sample_buf`.
    per_sample: bool,
}

impl Core {
    fn new(spec: &Spec<'_>) -> Result<Core, IoError> {
        let expr = parse_select_expr(&spec.header, spec.sel)?;
        let mut scratch = EvalScratch::default();
        scratch.set_selection(spec.plan.selection.as_ref().map(|s| s.keep.as_slice()));
        let bcf_site = if expr.as_ref().is_some_and(Expr::needs_info_values) {
            BcfSite::with_info_values()
        } else {
            BcfSite::default()
        };
        Ok(Core {
            expr,
            scratch,
            bcf_site,
            masked_buf: Vec::new(),
            decoded: SampleCache::default(),
            fields_buf: SampleFields::default(),
            sample_buf: Vec::new(),
        })
    }

    /// Run the record through the selection (all but thinning).
    fn pass(
        &mut self,
        spec: &Spec<'_>,
        raw: &[u8],
        m: &RecordMeta,
        rec_no: u64,
    ) -> Result<Option<Passed>, IoError> {
        let parse_err = |msg| IoError::Parse {
            path: spec.name.clone(),
            loc: Loc::Record(rec_no),
            msg,
        };
        let input_is_bcf = spec.in_format == Format::Bcf;
        let view = if input_is_bcf {
            self.bcf_site.view(raw, &spec.dicts, m.span)
        } else {
            SiteView::from_text(raw, m.span)
        }
        .map_err(parse_err)?;
        let sel = spec.sel;

        // 1. Static site filters
        if !sel.sites.static_pass(&view) {
            return Ok(None);
        }

        let mut per_sample = false;
        // 2. Stats filters, masking and the expression
        let stats = if spec.need_stats {
            let keep = spec.plan.selection.as_ref().map(|s| s.keep.as_slice());
            let n_alts = if view.alts.is_empty() || view.alts == b"." {
                0
            } else {
                view.alts.split(|&b| b == b',').count()
            };
            let n_alleles = if n_alts == 0 { 1 } else { 1 + n_alts };
            let fields = spec.fields;
            let bcf_gt;
            let text_gt;
            // A BCF record's statistics are counted from its bytes, unless
            // its samples are decoded anyway for the stages after this one.
            let mut counted = None;
            let gt: &dyn GenotypeSource = if input_is_bcf {
                bcf_gt = bcf_genotypes_with(raw, &spec.dicts, fields).map_err(parse_err)?;
                counted = Some(&bcf_gt);
                &bcf_gt
            } else {
                text_gt = text_genotypes_with(raw, fields).map_err(parse_err)?;
                &text_gt
            };
            let fast = spec.decode_only && spec.share_decoded && spec.tally.is_none();
            let cached;
            let gt: &dyn GenotypeSource = if !fast && (spec.share_decoded || spec.tally.is_some()) {
                self.decoded.start_record(gt.n_samples());
                cached = self.decoded.wrap(gt);
                counted = None;
                &cached
            } else {
                gt
            };
            let st = if fast {
                // Decode only; the checks are those of `site_stats_into`.
                let n = gt.n_samples();
                if let Some(k) = keep
                    && k.len() != n
                {
                    return Err(parse_err(format!(
                        "the sample selection has {} entries; the record has {n} samples",
                        k.len()
                    )));
                }
                self.masked_buf.clear();
                self.masked_buf.resize(n, false);
                self.decoded
                    .decode_record(gt, keep, |_, f| {
                        for a in f.gt.iter().flat_map(|g| g.alleles.iter().flatten()) {
                            if *a as usize >= n_alleles {
                                return Err(format!(
                                    "genotype allele {a} but the site lists {n_alleles} alleles"
                                ));
                            }
                        }
                        Ok(())
                    })
                    .map_err(parse_err)?;
                SiteStats::default()
            } else if let Some(bcf) = counted {
                bcf.site_stats(keep, &spec.plan.mask, n_alleles, &mut self.masked_buf)
                    .map_err(parse_err)?
            } else {
                site_stats(gt, keep, &spec.plan.mask, n_alleles, &mut self.masked_buf)
                    .map_err(parse_err)?
            };

            if !sel.sites.stats_pass(&st) {
                return Ok(None);
            }

            // 3. The include/exclude expression sees the masked genotypes of the
            // selected samples.
            if let Some(e) = self
                .expr
                .as_ref()
                .filter(|e| e.needs_genotypes() || e.needs_stats())
            {
                let mg = MaskedGenotypes {
                    inner: gt,
                    masked: &self.masked_buf,
                };
                let gt_arg = e.needs_genotypes().then_some(&mg as &dyn GenotypeSource);
                let st_arg = e.needs_stats().then_some(&st);
                if spec.filter_samples {
                    match e
                        .eval_query(
                            spec.exclude,
                            &view,
                            gt_arg,
                            st_arg,
                            &mut self.scratch,
                            &mut self.sample_buf,
                        )
                        .map_err(parse_err)?
                    {
                        QueryPass::Drop => return Ok(None),
                        QueryPass::All => {}
                        QueryPass::Samples => per_sample = true,
                    }
                } else {
                    let hold = e
                        .eval(&view, gt_arg, st_arg, &mut self.scratch)
                        .map_err(parse_err)?;
                    if hold == spec.exclude {
                        return Ok(None);
                    }
                }
            }
            Some(st)
        } else {
            None
        };

        // 3b. An expression over the site fields alone needs no genotypes.
        if let Some(e) = self
            .expr
            .as_ref()
            .filter(|e| !e.needs_genotypes() && !e.needs_stats())
        {
            let hold = e
                .eval(&view, None, None, &mut self.scratch)
                .map_err(parse_err)?;
            if hold == spec.exclude {
                return Ok(None);
            }
        }
        Ok(Some(Passed {
            stats,
            pos: view.pos,
            per_sample,
        }))
    }
}

/// The state of the selection over one input.
pub struct SelectStage<'a> {
    spec: Spec<'a>,
    core: Core,
    guard: Option<OrderGuard>,
    thinner: Option<Thinner>,
    threads: NonZeroUsize,
    tally: SampleTally,
}

impl<'a> SelectStage<'a> {
    /// The stage for `input`. `command` names the command in messages, `base`
    /// carries the rewrite settings the selection does not own (INFO pruning,
    /// AC recalculation), and `also_stats` asks for site statistics even when
    /// no filter needs them. The stage runs on one thread until
    /// [`Self::set_threads`].
    pub fn new(
        command: &str,
        input: &dyn RecordReader,
        sel: &'a VcfSelect,
        base: &RewritePlan,
        also_stats: bool,
    ) -> Result<Self, IoError> {
        sel.validate()?;
        let header = vcf_header(command, input)?;
        let in_format = input.format();
        let plan = rewrite_plan(&header, sel, base)?;
        let expr = parse_select_expr(&header, sel)?;
        let exclude = sel.expr.as_ref().is_some_and(|e| e.exclude);
        let dicts = if in_format == Format::Bcf {
            match &header.bcf {
                Some(d) => d.clone(),
                None => BcfDicts::from_header(&header.raw).map_err(|e| IoError::Detect {
                    path: input.name().to_string(),
                    msg: e,
                })?,
            }
        } else {
            BcfDicts::default()
        };
        // --thin needs each chromosome in one block with ascending positions,
        // in any chromosome order (what tabix needs, and what the windowed
        // statistics check): it keeps a last position per chromosome.
        let guard = sel
            .sites
            .thin
            .is_some()
            .then(|| OrderGuard::contiguous(input.name(), in_format));
        let selection_needs_stats = sel.sites.needs_genotypes()
            || !plan.mask.is_empty()
            || expr
                .as_ref()
                .is_some_and(|e| e.needs_genotypes() || e.needs_stats());
        let need_stats = selection_needs_stats || also_stats;
        let spec = Spec {
            sel,
            name: input.name().to_string(),
            header,
            in_format,
            plan,
            exclude,
            dicts,
            need_stats,
            stats_only_for_decode: !selection_needs_stats,
            decode_only: false,
            share_decoded: false,
            fields: GtFields::ALL,
            tally: None,
            filter_samples: false,
        };
        let core = Core::new(&spec)?;
        Ok(SelectStage {
            spec,
            core,
            guard,
            thinner: sel.sites.thinner(),
            threads: NonZeroUsize::MIN,
            tally: SampleTally::default(),
        })
    }

    /// Evaluate --include/--exclude with bcftools query's rule (D172) instead of D71's.
    pub fn filter_samples(&mut self) {
        self.spec.filter_samples = true;
    }

    /// Whether --include/--exclude reads genotypes (`Expr::needs_genotypes`).
    pub fn expr_reads_samples(&self) -> bool {
        self.core.expr.as_ref().is_some_and(Expr::needs_genotypes)
    }

    /// The first sample subscript of --include/--exclude (`GT[0]`), if any.
    pub fn expr_sample_subscript(&self) -> Option<String> {
        self.core.expr.as_ref().and_then(Expr::sample_subscript)
    }

    /// Keep each decoded sample for the record's consumers (see
    /// [`Hit::genotypes`]), so a later stage reads the sample from the
    /// selection instead of parsing it again. `fields` names the optional
    /// FORMAT fields (besides GT) the selection and its consumers read; the
    /// others are not decoded. The selection's own reads (masking rules, site
    /// filters, an expression) must be in `fields`.
    pub fn share_decoded(&mut self, fields: GtFields) {
        self.spec.share_decoded = true;
        self.spec.fields = fields;
    }

    /// For a consumer that reads the samples and nothing computed from them
    /// ([`Hit::stats`] is then an empty placeholder): when no filter, mask
    /// rule or expression needs the site statistics, skip computing them and
    /// only decode the samples (after [`Self::share_decoded`]).
    pub fn decode_only(&mut self) {
        self.spec.decode_only = self.spec.stats_only_for_decode;
    }

    /// Decode only the optional FORMAT fields `fields` names, without
    /// keeping the samples for the consumers (see [`Self::share_decoded`],
    /// which also sets this).
    pub fn decode_fields(&mut self, fields: GtFields) {
        self.spec.fields = fields;
    }

    /// Let the workers keep the per-sample counters `wants` names for the
    /// `samples` (their indices in the record) and sum them per batch, so a
    /// consumer needs no samples for them. Returns false, and does nothing,
    /// when that cannot be: on one thread, with `--thin` (a thinned record
    /// must not count, and thinning runs after the workers), or when the
    /// selection computes no site statistics. [`Self::tally`] then holds the
    /// totals once [`run_selection`] returns; [`Hit::tally_called`] carries
    /// each record's bitset.
    pub fn enable_tally(&mut self, wants: TallyWants, samples: &[usize]) -> bool {
        if wants.is_empty()
            || self.threads.get() < 2
            || self.thinner.is_some()
            || !self.spec.need_stats
        {
            return false;
        }
        self.spec.tally = Some(TallySpec {
            wants,
            samples: samples.to_vec(),
        });
        true
    }

    /// The counters workers summed over every selected record.
    pub fn tally(&self) -> &SampleTally {
        &self.tally
    }

    /// Run the selection on `n` threads: the calling thread reads, orders and
    /// hands records to the sink, and `n - 1` workers decode and filter them.
    /// The sink sees the same records in the same order at every `n`.
    pub fn set_threads(&mut self, n: NonZeroUsize) {
        self.threads = n;
    }

    pub fn name(&self) -> &str {
        &self.spec.name
    }

    pub fn header(&self) -> &VcfHeader {
        &self.spec.header
    }

    pub fn input_is_bcf(&self) -> bool {
        self.spec.in_format == Format::Bcf
    }

    /// The rewrite of each passing record: sample selection and masking
    /// (plus whatever `base` carried).
    pub fn plan(&self) -> &RewritePlan {
        &self.spec.plan
    }

    /// Records a batch may hold: about `BATCH_GENOTYPES` decoded genotypes.
    fn batch_limit(&self) -> usize {
        (BATCH_GENOTYPES / self.spec.header.samples.len().max(1)).max(1)
    }
}

/// Test instrumentation: the most decoded genotypes the parallel stage held
/// in flight since [`reset_inflight_peak`]. (The running count is local to
/// each run, so an error cannot leave it off.)
static INFLIGHT_PEAK: AtomicUsize = AtomicUsize::new(0);

/// Reset the in-flight peak counter (see [`inflight_peak`]).
#[doc(hidden)]
pub fn reset_inflight_peak() {
    INFLIGHT_PEAK.store(0, Ordering::SeqCst);
}

/// The most decoded genotypes held in flight by the parallel selection
/// since [`reset_inflight_peak`]. For tests.
#[doc(hidden)]
pub fn inflight_peak() -> usize {
    INFLIGHT_PEAK.load(Ordering::SeqCst)
}

/// Test hook: the batch number a worker panics on (`u64::MAX` = never).
#[cfg(feature = "test-hooks")]
static PANIC_AT_BATCH: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(u64::MAX);

/// Make the worker that takes batch `seq` (0-based) panic, to test that the
/// stage reports the panic instead of hanging. `None` turns it off. Only with
/// the `test-hooks` feature.
#[cfg(feature = "test-hooks")]
pub fn set_worker_panic_at(seq: Option<u64>) {
    PANIC_AT_BATCH.store(seq.unwrap_or(u64::MAX), Ordering::SeqCst);
}

/// What a worker keeps of one passing record.
struct HitRec {
    /// Index within the batch.
    idx: usize,
    pos: u64,
    stats: Option<SiteStats>,
    /// The record's mask flags in `Hits::masks`.
    masks: std::ops::Range<usize>,
    slot: Option<RecordSlot>,
    /// The record's `--het` bitset in `Hits::called`.
    called: std::ops::Range<usize>,
    /// The record's per-sample passes in `Hits::sample_pass`.
    sample_pass: Option<std::ops::Range<usize>>,
}

/// The passing records of one batch.
#[derive(Default)]
struct Hits {
    recs: Vec<HitRec>,
    masks: Vec<bool>,
    packed: PackedGenotypes,
    called: Vec<u64>,
    sample_pass: Vec<bool>,
    /// This batch's per-sample counters.
    tally: SampleTally,
}

impl Hits {
    fn clear(&mut self) {
        self.recs.clear();
        self.masks.clear();
        self.packed.clear();
        self.called.clear();
        self.sample_pass.clear();
        self.tally.clear();
    }
}

struct Job {
    seq: u64,
    /// Records before this batch in the stream.
    base: u64,
    batch: RecordBatch,
    hits: Hits,
}

/// A batch back from a worker: its passing records, then the error that
/// ended it (if any).
struct Done {
    job: Job,
    err: Option<IoError>,
}

impl Sequenced for Done {
    fn seq(&self) -> u64 {
        self.job.seq
    }
}

/// Where the driver and the workers of a parallel run meet: the batches to
/// work on, and the batches done. Whoever waits there sleeps until it is
/// woken ([`gloci_io::handoff`]).
type Handoff = gloci_io::handoff::Handoff<Job, Done>;

fn work_batch(spec: &Spec<'_>, core: &mut Core, mut job: Job, flight: &AtomicUsize) -> Done {
    let mut err = None;
    #[cfg(feature = "test-hooks")]
    if job.seq == PANIC_AT_BATCH.load(Ordering::SeqCst) {
        // `resume_unwind` panics without the panic hook's message.
        std::panic::resume_unwind(Box::new("test hook: worker panic"));
    }
    for i in 0..job.batch.len() {
        let rec_no = job.base + i as u64 + 1;
        let m = &job.batch.recs[i];
        let passed = match core.pass(spec, job.batch.raw(i), m, rec_no) {
            Ok(Some(p)) => p,
            Ok(None) => continue,
            Err(e) => {
                err = Some(e);
                break;
            }
        };
        let hits = &mut job.hits;
        let mask_start = hits.masks.len();
        let mut slot = None;
        let called_start = hits.called.len();
        if let (Some(t), Some(st)) = (&spec.tally, &passed.stats) {
            let reader = core.decoded.reader();
            let src = MaskedGenotypes {
                inner: &reader,
                masked: &core.masked_buf,
            };
            if let Err(msg) = hits.tally.add_record(
                t.wants,
                &t.samples,
                &src,
                &mut core.fields_buf,
                st,
                &mut hits.called,
            ) {
                err = Some(IoError::Parse {
                    path: spec.name.clone(),
                    loc: Loc::Record(rec_no),
                    msg,
                });
                break;
            }
        }
        if passed.stats.is_some() {
            hits.masks.extend_from_slice(&core.masked_buf);
            if spec.share_decoded {
                let keep = spec.plan.selection.as_ref().map(|s| s.keep.as_slice());
                match hits.packed.push_record_selected(
                    &core.decoded.reader(),
                    spec.fields,
                    &core.masked_buf,
                    keep,
                ) {
                    Ok(s) => slot = Some(s),
                    Err(msg) => {
                        err = Some(IoError::Parse {
                            path: spec.name.clone(),
                            loc: Loc::Record(rec_no),
                            msg,
                        });
                        break;
                    }
                }
            }
        }
        let sample_pass = passed.per_sample.then(|| {
            let start = hits.sample_pass.len();
            hits.sample_pass.extend_from_slice(&core.sample_buf);
            start..hits.sample_pass.len()
        });
        hits.recs.push(HitRec {
            idx: i,
            pos: passed.pos,
            stats: passed.stats,
            masks: mask_start..hits.masks.len(),
            slot,
            called: called_start..hits.called.len(),
            sample_pass,
        });
    }
    let held = job.hits.packed.n_samples_held();
    let now = flight.fetch_add(held, Ordering::SeqCst) + held;
    INFLIGHT_PEAK.fetch_max(now, Ordering::SeqCst);
    Done { job, err }
}

fn worker(
    spec: &Spec<'_>,
    handoff: &Handoff,
    abort: &AtomicBool,
    flight: &AtomicUsize,
    run: &gloci_io::cancel::Watch,
) {
    let _present = Present(handoff);
    let _run = run.enter();
    // The expression parsed for the driver's core, so this cannot fail; a
    // worker that cannot start leaves its jobs to the others, and the driver
    // reports a stage that has no workers left.
    let Ok(mut core) = Core::new(spec) else {
        return;
    };
    while let Some(job) = handoff.next_job() {
        if abort.load(Ordering::Relaxed) {
            break;
        }
        let (seq, base) = (job.seq, job.base);
        // The host cancelled the run (D91): the driver waits for this batch, so answer it.
        if run.is_cancelled() {
            let job = Job {
                seq,
                base,
                batch: RecordBatch::default(),
                hits: Hits::default(),
            };
            let err = Some(IoError::Cancelled);
            handoff.finish(Done { job, err });
            break;
        }
        let ran = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            work_batch(spec, &mut core, job, flight)
        }));
        // A panic must reach the driver, which waits for this batch.
        let (d, panicked) = match ran {
            Ok(d) => (d, false),
            Err(_) => {
                let err = IoError::Internal(format!(
                    "a selection worker panicked on the batch starting at record {}; this is a bug in gloci, please report it",
                    base + 1
                ));
                let job = Job {
                    seq,
                    base,
                    batch: RecordBatch::default(),
                    hits: Hits::default(),
                };
                (
                    Done {
                        job,
                        err: Some(err),
                    },
                    true,
                )
            }
        };
        handoff.finish(d);
        if panicked {
            break;
        }
    }
}

/// Read `input` to its end, calling `sink` on each record that passes `stage`,
/// in input order. Returns the number of records read. With more than one
/// thread the records are decoded and filtered on workers; the sink and the
/// result are the same as on one thread.
pub fn run_selection<F>(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    stage: &mut SelectStage<'_>,
    sink: F,
) -> Result<u64, IoError>
where
    F: FnMut(Hit<'_>) -> Result<(), IoError>,
{
    let records_in = if stage.threads.get() > 1 {
        run_parallel(input.as_mut(), chroms, stage, sink)?
    } else {
        run_sequential(input.as_mut(), chroms, stage, sink)?
    };
    if let Some(s) = input.diagnostics().summary(input.name()) {
        log::warn!("{s}");
    }
    Ok(records_in)
}

/// What ended a batch early: an error in reading it or in the record order.
/// The batch holds the records before it, which still run (an earlier
/// record's own error then comes first).
fn read_guarded(
    input: &mut dyn RecordReader,
    chroms: &mut ChromTable,
    guard: &mut Option<OrderGuard>,
    batch: &mut RecordBatch,
    max: usize,
) -> (bool, Option<IoError>) {
    let (more, read_err) = match input.read_batch(chroms, batch, max) {
        Ok(more) => (more, None),
        Err(e) => (false, Some(e)),
    };
    if let Some(g) = guard.as_mut() {
        for i in 0..batch.len() {
            if let Err(e) = g.check_range(chroms, batch, i..i + 1) {
                batch.recs.truncate(i);
                return (false, Some(e));
            }
        }
        g.settle(chroms, batch, batch.len());
    }
    (more, read_err)
}

fn run_sequential<F>(
    input: &mut dyn RecordReader,
    chroms: &mut ChromTable,
    stage: &mut SelectStage<'_>,
    mut sink: F,
) -> Result<u64, IoError>
where
    F: FnMut(Hit<'_>) -> Result<(), IoError>,
{
    let limit = stage.batch_limit();
    let run = gloci_io::cancel::Watch::current();
    let mut batch = RecordBatch::default();
    let mut max = 1;
    let mut records_in = 0u64;
    loop {
        batch.clear();
        let (more, ended) =
            read_guarded(input, chroms, &mut stage.guard, &mut batch, max.min(limit));
        for i in 0..batch.len() {
            run.check()?;
            records_in += 1;
            let m = &batch.recs[i];
            let raw = batch.raw(i);
            let Some(p) = stage.core.pass(&stage.spec, raw, m, records_in)? else {
                continue;
            };
            if let Some(t) = stage.thinner.as_mut()
                && !t.keep(m.chrom, p.pos)
            {
                continue;
            }
            // The selected samples were decoded by the pass; masked ones read as missing.
            let shared = stage.core.decoded.reader();
            let shared = MaskedGenotypes {
                inner: &shared,
                masked: &stage.core.masked_buf,
            };
            sink(Hit {
                raw,
                rec_no: records_in,
                meta: m,
                stats: p.stats.as_ref(),
                masked: &stage.core.masked_buf,
                genotypes: (stage.spec.share_decoded && p.stats.is_some())
                    .then_some(&shared as &dyn GenotypeSource),
                called: &[],
                sample_pass: p.per_sample.then_some(stage.core.sample_buf.as_slice()),
            })?;
        }
        if let Some(e) = ended {
            return Err(e);
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    Ok(records_in)
}

fn run_parallel<F>(
    input: &mut dyn RecordReader,
    chroms: &mut ChromTable,
    stage: &mut SelectStage<'_>,
    mut sink: F,
) -> Result<u64, IoError>
where
    F: FnMut(Hit<'_>) -> Result<(), IoError>,
{
    let workers = stage.threads.get() - 1;
    // At most this many batches are in flight (read, being worked on, or
    // finished but not yet handed to the sink).
    let cap = 2 * workers;
    let limit = stage.batch_limit();
    let SelectStage {
        spec,
        guard,
        thinner,
        tally,
        ..
    } = stage;
    let spec: &Spec<'_> = spec;
    let abort = AtomicBool::new(false);
    let run = gloci_io::cancel::Watch::current();
    let run = &run;
    let flight = AtomicUsize::new(0);
    let flight = &flight;
    let handoff = Handoff::new(workers);
    let handoff = &handoff;
    std::thread::scope(|scope| {
        // Before the first worker starts: whatever ends this closure (even a
        // worker that cannot be spawned) lets the workers go.
        let _driving = Driving {
            handoff,
            abort: &abort,
        };
        for _ in 0..workers {
            let abort = &abort;
            scope.spawn(move || worker(spec, handoff, abort, flight, run));
        }
        let no_workers =
            || IoError::Internal("the selection workers stopped before the input was read".into());

        let mut spare: Vec<(RecordBatch, Hits)> = Vec::new();
        let (mut sent, mut taken) = (0u64, 0u64);
        let mut records_in = 0u64;
        let mut max = 1;
        let mut eof = false;
        let mut workers_gone = false;
        // A read or order error ends the input: the batches before it are
        // still taken (so an earlier record's error comes first).
        let mut input_err: Option<IoError> = None;
        (|| -> Result<u64, IoError> {
            loop {
                run.check()?;
                while !eof && input_err.is_none() && sent - taken < cap as u64 {
                    let (mut batch, mut hits) = spare.pop().unwrap_or_default();
                    batch.clear();
                    hits.clear();
                    let (more, ended) =
                        read_guarded(input, chroms, guard, &mut batch, max.min(limit));
                    input_err = ended;
                    max = next_batch_size(&batch, max);
                    eof = !more;
                    if batch.is_empty() {
                        continue;
                    }
                    let n = batch.len() as u64;
                    let job = Job {
                        seq: sent,
                        base: records_in,
                        batch,
                        hits,
                    };
                    if !handoff.send(job) {
                        // Every worker has exited. One may have left an
                        // error (a panic) for a batch already sent: stop
                        // reading and take the batches in order, so that
                        // error is the one reported. With none, the wait
                        // below fails and says the workers stopped.
                        eof = true;
                        workers_gone = true;
                        break;
                    }
                    sent += 1;
                    records_in += n;
                }
                if sent == taken {
                    break;
                }
                let done = handoff.take(taken).ok_or_else(|| {
                    if run.is_cancelled() {
                        IoError::Cancelled
                    } else {
                        no_workers()
                    }
                })?;
                taken += 1;
                let Done { job, err } = done;
                let held = job.hits.packed.n_samples_held();
                for h in &job.hits.recs {
                    let m = &job.batch.recs[h.idx];
                    if let Some(t) = thinner.as_mut()
                        && !t.keep(m.chrom, h.pos)
                    {
                        continue;
                    }
                    let packed = h.slot.map(|s| job.hits.packed.record(s));
                    sink(Hit {
                        raw: job.batch.raw(h.idx),
                        rec_no: job.base + h.idx as u64 + 1,
                        meta: m,
                        stats: h.stats.as_ref(),
                        masked: &job.hits.masks[h.masks.clone()],
                        genotypes: packed.as_ref().map(|p| p as &dyn GenotypeSource),
                        called: &job.hits.called[h.called.clone()],
                        sample_pass: h.sample_pass.clone().map(|r| &job.hits.sample_pass[r]),
                    })?;
                }
                if spec.tally.is_some() {
                    tally.add(&job.hits.tally);
                }
                flight.fetch_sub(held, Ordering::SeqCst);
                if let Some(e) = err {
                    return Err(e);
                }
                spare.push((job.batch, job.hits));
            }
            if workers_gone {
                return Err(no_workers());
            }
            match input_err.take() {
                Some(e) => Err(e),
                None => Ok(records_in),
            }
        })()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn job(seq: u64) -> Job {
        Job {
            seq,
            base: 0,
            batch: RecordBatch::default(),
            hits: Hits::default(),
        }
    }

    /// The driver sees what the workers hand back in the order it asks for,
    /// and learns that the last worker has left.
    #[test]
    fn handoff_returns_batches_by_number_and_notices_the_last_worker_leaving() {
        let handoff = Handoff::new(2);
        let stop = AtomicBool::new(false);
        std::thread::scope(|s| {
            // The test is the driver. Its end of the hand-off is dropped
            // however this closure ends, a failed assertion included: the
            // scope waits for the workers, and they wait for the hand-off.
            let _driving = Driving {
                handoff: &handoff,
                abort: &stop,
            };
            for _ in 0..2 {
                s.spawn(|| {
                    let _present = Present(&handoff);
                    while let Some(job) = handoff.next_job() {
                        // The later batches come back first.
                        std::thread::sleep(std::time::Duration::from_millis(20 - 5 * job.seq));
                        handoff.finish(Done { job, err: None });
                    }
                });
            }
            for seq in 0..4 {
                assert!(handoff.send(job(seq)));
            }
            for seq in 0..4 {
                assert_eq!(handoff.take(seq).map(|d| d.job.seq), Some(seq));
            }
        });
        // Both workers have returned: nothing more comes back, nothing is taken.
        assert!(handoff.take(4).is_none());
        assert!(!handoff.send(job(4)));
    }

    /// How often the kernel has switched this thread out while it could run.
    #[cfg(target_os = "linux")]
    fn switched_out() -> u64 {
        let status = std::fs::read_to_string("/proc/thread-self/status").unwrap();
        let line = status
            .lines()
            .find_map(|l| l.strip_prefix("nonvoluntary_ctxt_switches:"))
            .expect("the kernel reports context switches");
        line.trim().parse().unwrap()
    }

    /// A worker with nothing to do sleeps until a job wakes it; it does not
    /// retry around `sched_yield`. With every CPU busy, a thread that yields
    /// is switched out at each yield while it could still run, and the kernel
    /// counts that; one that sleeps is not. (A channel's blocking receive
    /// yields several times before it sleeps: on Linux each yield on a busy
    /// CPU costs a scheduler slice, which is what made the selection slower
    /// with more threads than cores.)
    #[cfg(target_os = "linux")]
    #[test]
    fn a_waiting_worker_sleeps_instead_of_yielding() {
        let rounds = 100u64;
        let stop = AtomicBool::new(false);
        let handoff = Handoff::new(1);
        let (taken, switched) = std::thread::scope(|s| {
            // The test is the driver. Dropping its end of the hand-off stops
            // the busy threads and lets the worker go, however this closure
            // ends: the scope waits for every thread, so a failed assertion
            // that left them running would hang the test with every CPU busy.
            let driving = Driving {
                handoff: &handoff,
                abort: &stop,
            };
            let cpus = std::thread::available_parallelism().map_or(4, |n| n.get());
            for _ in 0..cpus {
                s.spawn(|| {
                    while !stop.load(Ordering::Relaxed) {
                        std::hint::spin_loop();
                    }
                });
            }
            let worker = s.spawn(|| {
                let _present = Present(&handoff);
                let before = switched_out();
                let mut taken = 0;
                while let Some(job) = handoff.next_job() {
                    taken += 1;
                    handoff.finish(Done { job, err: None });
                }
                (taken, switched_out() - before)
            });
            for seq in 0..rounds {
                // Long enough for the worker to be waiting again.
                std::thread::sleep(std::time::Duration::from_millis(2));
                assert!(handoff.send(job(seq)));
                assert!(handoff.take(seq).is_some());
            }
            drop(driving);
            worker.join().unwrap()
        });
        assert_eq!(taken, rounds);
        assert!(
            switched < rounds / 2,
            "the waiting worker was switched out {switched} times while runnable in {rounds} waits"
        );
    }
}
