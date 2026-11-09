//! `gloci vcf filter`: the selection pipeline (spec §7.8, decision D50).

use crate::view::{ViewOptions, open_view, select_stream};
use gloci_core::{ChromOrder, ChromTable, VcfSpan};
use gloci_io::{
    BcfDicts, BcfText, Format, Header, IoError, Loc, OrderGuard, OutFormat, ReaderOptions,
    RecordBatch, RecordReader, RecordWriter, TextToBcf, UsageError, VcfHeader, next_batch_size,
};
use gloci_vcf::{
    BcfSite, MaskRules, RewritePlan, SiteFilters, SiteView, Thinner, bcf_genotypes, rewrite_header,
    rewrite_text, select_samples, site_stats, text_genotypes,
};

/// Options for `vcf filter`.
#[derive(Clone, Debug, Default)]
pub struct VcfFilterOptions {
    pub include_samples: Option<Vec<String>>,
    pub exclude_samples: Vec<String>,
    pub mask: MaskRules,
    pub sites: SiteFilters,
    pub plan: RewritePlan,
    pub regions: Vec<String>,
    pub regions_file: Option<String>,
}

/// Summary of records read and written by `vcf filter`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct FilterSummary {
    pub records_in: u64,
    pub records_out: u64,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

impl VcfFilterOptions {
    /// The checks that need no input: thresholds within range.
    pub fn validate(&self) -> Result<(), IoError> {
        let unit = [
            ("--min-maf", self.sites.min_maf),
            ("--max-maf", self.sites.max_maf),
            ("--max-missing", self.sites.max_missing),
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
        if let (Some(lo), Some(hi)) = (self.sites.min_alleles, self.sites.max_alleles)
            && lo > hi
        {
            return Err(usage(format!(
                "--min-alleles {lo} is more than --max-alleles {hi}; no site could pass"
            )));
        }
        Ok(())
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
        }
    }
}

/// `input`'s VCF header, or an error when it is not VCF/BCF.
fn vcf_header(input: &dyn RecordReader) -> Result<VcfHeader, IoError> {
    let format = input.format();
    if format != Format::Vcf && format != Format::Bcf {
        return Err(usage(format!(
            "vcf filter works on VCF and BCF input, but {} is {}; give it a VCF or BCF file",
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
/// exclude lists, else `o.plan`'s) and the masking rules.
fn rewrite_plan(h: &VcfHeader, o: &VcfFilterOptions) -> Result<RewritePlan, IoError> {
    let mut plan = o.plan.clone();
    if o.include_samples.is_some() || !o.exclude_samples.is_empty() {
        plan.selection = Some(select_samples(
            &h.samples,
            o.include_samples.as_deref(),
            &o.exclude_samples,
        )?);
    }
    plan.mask = o.mask;
    Ok(plan)
}

/// What `vcf filter` writes for `input`: the output format (`requested`,
/// else the input's) and the header (samples, INFO lines and so on as the
/// rewrite leaves them).
pub fn vcf_filter_output(
    input: &dyn RecordReader,
    requested: Option<OutFormat>,
    o: &VcfFilterOptions,
) -> Result<(OutFormat, Header), IoError> {
    let h = vcf_header(input)?;
    let to = requested.unwrap_or(OutFormat::family(input.format()));
    if to != OutFormat::Vcf && to != OutFormat::Bcf {
        return Err(usage(format!(
            "vcf filter writes VCF or BCF, not {}; name the output .vcf, .vcf.gz or .bcf, \
or pass -O vcf or -O bcf",
            to.name()
        )));
    }
    let plan = rewrite_plan(&h, o)?;
    let header = VcfHeader::parse(rewrite_header(&h, &plan)).map_err(|e| IoError::Detect {
        path: input.name().to_string(),
        msg: e,
    })?;
    Ok((to, Header::Vcf(header)))
}

/// An input opened by [`open_vcf_filter`], its regions already selected.
pub struct VcfFilterInput {
    reader: Box<dyn RecordReader>,
    chroms: ChromTable,
}

impl VcfFilterInput {
    pub fn reader(&self) -> &dyn RecordReader {
        self.reader.as_ref()
    }
}

/// Open `path` for `vcf filter`. With `--region`/`--regions-file`, only the
/// records that overlap a region are read: through the index when the file
/// has one, else by scanning it (as `view` does).
pub fn open_vcf_filter(
    path: &str,
    ro: &ReaderOptions,
    o: &VcfFilterOptions,
) -> Result<VcfFilterInput, IoError> {
    let (reader, chroms) = open_view(path, &o.view_options(), ro)?;
    Ok(VcfFilterInput { reader, chroms })
}

/// Filter an input opened by [`open_vcf_filter`] and write the passing
/// records to `w` (opened with [`vcf_filter_output`]'s format and header).
pub fn vcf_filter_file(
    input: VcfFilterInput,
    o: &VcfFilterOptions,
    w: &mut dyn RecordWriter,
) -> Result<FilterSummary, IoError> {
    let VcfFilterInput { reader, mut chroms } = input;
    filter(reader, &mut chroms, o, w)
}

/// Filter the VCF/BCF stream `input` and write passing records to `w`. The
/// regions of `o` are selected from the stream as it is read.
pub fn vcf_filter(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &VcfFilterOptions,
    w: &mut dyn RecordWriter,
) -> Result<FilterSummary, IoError> {
    o.validate()?;
    vcf_header(input.as_ref())?;
    let region_ro = ReaderOptions {
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Ref,
        invalid: gloci_io::InvalidPolicy::Fail,
        threads: std::num::NonZero::new(1).unwrap(),
    };
    let input = select_stream(input, &o.view_options(), &region_ro, chroms)?;
    filter(input, chroms, o, w)
}

/// The pipeline: site filters, masking, thinning, rewriting.
fn filter(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &VcfFilterOptions,
    w: &mut dyn RecordWriter,
) -> Result<FilterSummary, IoError> {
    o.validate()?;
    let vcf_header = vcf_header(input.as_ref())?;
    let in_format = input.format();
    let plan = rewrite_plan(&vcf_header, o)?;

    let input_is_bcf = in_format == Format::Bcf;
    let output_is_bcf = w.is_bcf();

    let mut text_to_bcf = if output_is_bcf && (!input_is_bcf || !plan.is_passthrough()) {
        let rewritten_header_raw = rewrite_header(&vcf_header, &plan);
        let out_vcf_header =
            VcfHeader::parse(rewritten_header_raw).map_err(|e| IoError::Detect {
                path: input.name().to_string(),
                msg: e,
            })?;
        Some(TextToBcf::new(input.name(), &out_vcf_header)?)
    } else {
        None
    };

    let dicts = if input_is_bcf {
        match &vcf_header.bcf {
            Some(d) => d.clone(),
            None => BcfDicts::from_header(&vcf_header.raw).map_err(|e| IoError::Detect {
                path: input.name().to_string(),
                msg: e,
            })?,
        }
    } else {
        BcfDicts::default()
    };

    let mut bcf_site = BcfSite::default();
    let mut bcf_text = if input_is_bcf {
        Some(BcfText::for_header(input.name(), &vcf_header)?)
    } else {
        None
    };

    let mut guard: Option<OrderGuard> = if o.sites.thin.is_some() {
        Some(OrderGuard::new(input.name(), in_format))
    } else {
        None
    };
    let mut thinner: Option<Thinner> = o.sites.thinner();

    let need_stats = o.sites.needs_genotypes() || !plan.mask.is_empty() || plan.recalc_ac;
    let mut masked_buf = Vec::new();
    let mut rewritten_buf = Vec::new();
    let mut bcf_enc_buf = Vec::new();

    let mut batch = RecordBatch::default();
    let mut max = 1;
    let mut records_in = 0u64;
    let mut records_out = 0u64;

    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        if let Some(g) = guard.as_mut() {
            g.check(chroms, &batch)?;
        }
        for i in 0..batch.len() {
            records_in += 1;
            let m = &batch.recs[i];
            let raw = batch.raw(i);

            let view = if input_is_bcf {
                bcf_site.view(raw, &dicts, m.span)
            } else {
                SiteView::from_text(raw, m.span)
            }
            .map_err(|msg| IoError::Parse {
                path: input.name().to_string(),
                loc: Loc::Record(records_in),
                msg,
            })?;

            // 1. Static site filters
            if !o.sites.static_pass(&view) {
                continue;
            }

            // 2. Stats filters & masking
            let stats = if need_stats {
                let keep = plan.selection.as_ref().map(|s| s.keep.as_slice());
                let n_alts = if view.alts.is_empty() || view.alts == b"." {
                    0
                } else {
                    view.alts.split(|&b| b == b',').count()
                };
                let n_alleles = if n_alts == 0 { 1 } else { 1 + n_alts };
                let st = if input_is_bcf {
                    let gt = bcf_genotypes(raw, &dicts).map_err(|msg| IoError::Parse {
                        path: input.name().to_string(),
                        loc: Loc::Record(records_in),
                        msg,
                    })?;
                    site_stats(&gt, keep, &plan.mask, n_alleles, &mut masked_buf)
                } else {
                    let gt = text_genotypes(raw).map_err(|msg| IoError::Parse {
                        path: input.name().to_string(),
                        loc: Loc::Record(records_in),
                        msg,
                    })?;
                    site_stats(&gt, keep, &plan.mask, n_alleles, &mut masked_buf)
                }
                .map_err(|msg| IoError::Parse {
                    path: input.name().to_string(),
                    loc: Loc::Record(records_in),
                    msg,
                })?;

                if !o.sites.stats_pass(&st) {
                    continue;
                }
                Some(st)
            } else {
                None
            };

            // 3. Thinning (runs last on passing sites)
            if let Some(t) = thinner.as_mut()
                && !t.keep(m.chrom, view.pos)
            {
                continue;
            }

            // 4. Output record
            if input_is_bcf && output_is_bcf && plan.is_passthrough() {
                w.write_record(raw).map_err(|source| IoError::Write {
                    path: "output".to_string(),
                    source,
                })?;
            } else if !input_is_bcf && !output_is_bcf && plan.is_passthrough() {
                w.write_record(raw).map_err(|source| IoError::Write {
                    path: "output".to_string(),
                    source,
                })?;
            } else {
                let text_line = if input_is_bcf {
                    bcf_text.as_mut().unwrap().render(raw, records_in)?
                } else {
                    raw
                };

                let out_text = if !plan.is_passthrough() {
                    rewritten_buf.clear();
                    rewrite_text(
                        text_line,
                        &plan,
                        &masked_buf,
                        stats.as_ref(),
                        &mut rewritten_buf,
                    )
                    .map_err(|msg| IoError::Parse {
                        path: input.name().to_string(),
                        loc: Loc::Record(records_in),
                        msg,
                    })?;
                    &rewritten_buf[..]
                } else {
                    text_line
                };

                if output_is_bcf {
                    bcf_enc_buf.clear();
                    text_to_bcf
                        .as_mut()
                        .unwrap()
                        .encode(out_text, records_in, &mut bcf_enc_buf)?;
                    w.write_record(&bcf_enc_buf)
                        .map_err(|source| IoError::Write {
                            path: "output".to_string(),
                            source,
                        })?;
                } else {
                    w.write_record(out_text).map_err(|source| IoError::Write {
                        path: "output".to_string(),
                        source,
                    })?;
                }
            }
            records_out += 1;
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }

    if let Some(s) = input.diagnostics().summary(input.name()) {
        log::warn!("{s}");
    }

    Ok(FilterSummary {
        records_in,
        records_out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_io::open_reader;

    struct VecWriter(Vec<u8>);
    impl RecordWriter for VecWriter {
        fn write_record(&mut self, record: &[u8]) -> Result<(), std::io::Error> {
            self.0.extend_from_slice(record);
            self.0.push(b'\n');
            Ok(())
        }
        fn finish(self: Box<Self>) -> Result<(), IoError> {
            Ok(())
        }
    }

    #[test]
    fn vcf_filter_library_caller_regions() {
        let vcf_data = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\trs1\tA\tC\t50\tPASS\t.\n\
chr1\t200\trs2\tG\tT\t50\tPASS\t.\n\
chr2\t300\trs3\tT\tA\t50\tPASS\t.\n";
        let temp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), vcf_data).unwrap();
        let ro = ReaderOptions {
            format: Some(Format::Vcf),
            vcf_span: VcfSpan::Ref,
            invalid: gloci_io::InvalidPolicy::Fail,
            threads: std::num::NonZero::new(1).unwrap(),
        };
        let reader = open_reader(temp.path().to_str().unwrap(), &ro).unwrap();
        let mut chroms = ChromTable::first_appearance();
        let opts = VcfFilterOptions {
            regions: vec!["chr1:150-250".to_string()],
            ..Default::default()
        };
        let mut writer = VecWriter(Vec::new());
        let summary = vcf_filter(reader, &mut chroms, &opts, &mut writer).unwrap();
        // The region is selected before the filters, which see one record.
        assert_eq!(summary.records_in, 1);
        assert_eq!(summary.records_out, 1);
        let out_str = String::from_utf8(writer.0).unwrap();
        assert!(out_str.contains("rs2"));
        assert!(!out_str.contains("rs1"));
        assert!(!out_str.contains("rs3"));
    }
}
