//! The selection stage shared by the `vcf` commands: sample selection,
//! genotype masking, site filters, the include/exclude expression, regions
//! and thinning (spec §7.8, decisions D50 and D71).
//!
//! `vcf filter` writes the records that pass; `vcf extract` tabulates them.
//! A command opens its input with [`open_selected`], builds a
//! [`SelectStage`] and drives it with [`run_selection`].

use crate::view::{ViewOptions, open_view, select_stream};
use gloci_core::{ChromOrder, ChromTable, VcfSpan};
use gloci_io::{
    BcfDicts, Format, Header, IoError, Loc, OrderGuard, ReaderOptions, RecordBatch, RecordMeta,
    RecordReader, UsageError, VcfHeader, next_batch_size,
};
use gloci_vcf::{
    BcfSite, EvalScratch, Expr, GenotypeSource, GtFields, MaskRules, MaskedGenotypes, RewritePlan,
    SampleCache, SiteFilters, SiteStats, SiteView, Thinner, bcf_genotypes_with, parse_expr,
    select_samples, site_stats, text_genotypes_with,
};

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
            chrom_order: ChromOrder::Karyotype,
            assume_sorted: false,
            unordered: false,
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
    /// The samples the selection decoded, when [`SelectStage::share_decoded`]
    /// is on: wrap the record's genotype source with it to avoid decoding twice.
    pub decoded: &'a SampleCache,
}

/// The state of the selection over one input.
pub struct SelectStage<'a> {
    sel: &'a VcfSelect,
    name: String,
    header: VcfHeader,
    in_format: Format,
    plan: RewritePlan,
    expr: Option<Expr>,
    exclude: bool,
    scratch: EvalScratch,
    dicts: BcfDicts,
    bcf_site: BcfSite,
    guard: Option<OrderGuard>,
    thinner: Option<Thinner>,
    need_stats: bool,
    masked_buf: Vec<bool>,
    decoded: SampleCache,
    share_decoded: bool,
    fields: GtFields,
}

impl<'a> SelectStage<'a> {
    /// The stage for `input`. `command` names the command in messages, `base`
    /// carries the rewrite settings the selection does not own (INFO pruning,
    /// AC recalculation), and `also_stats` asks for site statistics even when
    /// no filter needs them.
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
        let mut scratch = EvalScratch::default();
        scratch.set_selection(plan.selection.as_ref().map(|s| s.keep.as_slice()));
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
        let bcf_site = if expr.as_ref().is_some_and(Expr::needs_info_values) {
            BcfSite::with_info_values()
        } else {
            BcfSite::default()
        };
        let guard = sel
            .sites
            .thin
            .is_some()
            .then(|| OrderGuard::new(input.name(), in_format));
        let need_stats = sel.sites.needs_genotypes()
            || !plan.mask.is_empty()
            || also_stats
            || expr
                .as_ref()
                .is_some_and(|e| e.needs_genotypes() || e.needs_stats());
        Ok(SelectStage {
            sel,
            name: input.name().to_string(),
            header,
            in_format,
            plan,
            expr,
            exclude,
            scratch,
            dicts,
            bcf_site,
            guard,
            thinner: sel.sites.thinner(),
            need_stats,
            masked_buf: Vec::new(),
            decoded: SampleCache::default(),
            share_decoded: false,
            fields: GtFields::ALL,
        })
    }

    /// Keep each decoded sample for the record's consumers (see [`Hit::decoded`]),
    /// so a later stage reads the sample from the cache instead of parsing it again.
    /// `fields` names the optional FORMAT fields (besides GT) the selection
    /// and its consumers read; the others are not decoded. The selection's own
    /// reads (masking rules, site filters, an expression) must be in `fields`.
    pub fn share_decoded(&mut self, fields: GtFields) {
        self.share_decoded = true;
        self.fields = fields;
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    pub fn input_is_bcf(&self) -> bool {
        self.in_format == Format::Bcf
    }

    /// The rewrite of each passing record: sample selection and masking
    /// (plus whatever `base` carried).
    pub fn plan(&self) -> &RewritePlan {
        &self.plan
    }

    fn check_batch(&mut self, chroms: &ChromTable, batch: &RecordBatch) -> Result<(), IoError> {
        if let Some(g) = self.guard.as_mut() {
            g.check(chroms, batch)?;
        }
        Ok(())
    }

    /// Run the record through the selection: `Some(stats)` when it passes.
    fn pass(
        &mut self,
        raw: &[u8],
        m: &RecordMeta,
        rec_no: u64,
    ) -> Result<Option<Option<SiteStats>>, IoError> {
        let parse_err = |msg| IoError::Parse {
            path: self.name.clone(),
            loc: Loc::Record(rec_no),
            msg,
        };
        let input_is_bcf = self.in_format == Format::Bcf;
        let view = if input_is_bcf {
            self.bcf_site.view(raw, &self.dicts, m.span)
        } else {
            SiteView::from_text(raw, m.span)
        }
        .map_err(parse_err)?;
        let sel = self.sel;

        // 1. Static site filters
        if !sel.sites.static_pass(&view) {
            return Ok(None);
        }

        // 2. Stats filters, masking and the expression
        let stats = if self.need_stats {
            let keep = self.plan.selection.as_ref().map(|s| s.keep.as_slice());
            let n_alts = if view.alts.is_empty() || view.alts == b"." {
                0
            } else {
                view.alts.split(|&b| b == b',').count()
            };
            let n_alleles = if n_alts == 0 { 1 } else { 1 + n_alts };
            let fields = if self.share_decoded {
                self.fields
            } else {
                GtFields::ALL
            };
            let bcf_gt;
            let text_gt;
            let gt: &dyn GenotypeSource = if input_is_bcf {
                bcf_gt = bcf_genotypes_with(raw, &self.dicts, fields).map_err(parse_err)?;
                &bcf_gt
            } else {
                text_gt = text_genotypes_with(raw, fields).map_err(parse_err)?;
                &text_gt
            };
            let cached;
            let gt: &dyn GenotypeSource = if self.share_decoded {
                self.decoded.start_record(gt.n_samples());
                cached = self.decoded.wrap(gt);
                &cached
            } else {
                gt
            };
            let st = site_stats(gt, keep, &self.plan.mask, n_alleles, &mut self.masked_buf)
                .map_err(parse_err)?;

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
                let hold = e
                    .eval(
                        &view,
                        e.needs_genotypes().then_some(&mg as &dyn GenotypeSource),
                        e.needs_stats().then_some(&st),
                        &mut self.scratch,
                    )
                    .map_err(parse_err)?;
                if hold == self.exclude {
                    return Ok(None);
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
            if hold == self.exclude {
                return Ok(None);
            }
        }

        // 4. Thinning (runs last on passing sites)
        if let Some(t) = self.thinner.as_mut()
            && !t.keep(m.chrom, view.pos)
        {
            return Ok(None);
        }
        Ok(Some(stats))
    }
}

/// Read `input` to its end, calling `sink` on each record that passes `stage`.
/// Returns the number of records read.
pub fn run_selection<F>(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    stage: &mut SelectStage<'_>,
    mut sink: F,
) -> Result<u64, IoError>
where
    F: FnMut(Hit<'_>) -> Result<(), IoError>,
{
    let mut batch = RecordBatch::default();
    let mut max = 1;
    let mut records_in = 0u64;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        stage.check_batch(chroms, &batch)?;
        for i in 0..batch.len() {
            records_in += 1;
            let m = &batch.recs[i];
            let raw = batch.raw(i);
            let Some(stats) = stage.pass(raw, m, records_in)? else {
                continue;
            };
            sink(Hit {
                raw,
                rec_no: records_in,
                meta: m,
                stats: stats.as_ref(),
                masked: &stage.masked_buf,
                decoded: &stage.decoded,
            })?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    if let Some(s) = input.diagnostics().summary(input.name()) {
        log::warn!("{s}");
    }
    Ok(records_in)
}
