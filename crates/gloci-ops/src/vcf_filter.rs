//! `gloci vcf filter`: the selection pipeline (spec §7.8, decision D50).

use gloci_core::{ChromTable, Interval, Region, RegionSet, VcfSpan};
use gloci_io::{
    BcfDicts, BcfText, Format, Header, IoError, Loc, OrderGuard, ReaderOptions, RecordBatch,
    RecordReader, RecordWriter, TextToBcf, UsageError, VcfHeader, next_batch_size,
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

/// Filter VCF/BCF `input` and write passing records to `w`.
pub fn vcf_filter(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &VcfFilterOptions,
    w: &mut dyn RecordWriter,
) -> Result<FilterSummary, IoError> {
    let in_format = input.format();
    if in_format != Format::Vcf && in_format != Format::Bcf {
        return Err(UsageError(format!(
            "vcf filter works on VCF and BCF input, but {} is {:?}",
            input.name(),
            in_format
        ))
        .into());
    }

    let Header::Vcf(vcf_header) = input.header().clone() else {
        return Err(IoError::Detect {
            path: input.name().to_string(),
            msg: "a VCF/BCF input has no VCF header".into(),
        });
    };

    let selection = if o.include_samples.is_some() || !o.exclude_samples.is_empty() {
        Some(select_samples(
            &vcf_header.samples,
            o.include_samples.as_deref(),
            &o.exclude_samples,
        )?)
    } else {
        o.plan.selection.clone()
    };

    let mut plan = o.plan.clone();
    plan.selection = selection;
    plan.mask = o.mask;

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

    let include_regions: Option<std::collections::HashMap<Vec<u8>, RegionSet>> =
        if !o.regions.is_empty() || o.regions_file.is_some() {
            let mut known = std::collections::HashSet::new();
            if let Some(d) = &vcf_header.bcf {
                known.extend(d.contigs.iter().filter(|c| !c.is_empty()).cloned());
            } else {
                known.extend(vcf_header.contigs.iter().map(|(c, _)| c.clone()));
            }
            let mut wanted: Vec<(Vec<u8>, Interval)> = Vec::new();
            for s in &o.regions {
                let r = Region::parse(s, |n| {
                    known.contains(n) || (known.is_empty() && !n.contains(&b':'))
                })
                .map_err(|e| IoError::Usage(UsageError(e.to_string())))?;
                let iv = r.interval(None);
                wanted.push((r.chrom, iv));
            }
            if let Some(f) = &o.regions_file {
                let ro = ReaderOptions {
                    format: Some(Format::Bed),
                    vcf_span: VcfSpan::Ref,
                    invalid: gloci_io::InvalidPolicy::Fail,
                    threads: std::num::NonZero::new(1).unwrap(),
                };
                wanted.extend(crate::view::read_regions_file(f, &ro)?);
            }
            let mut by_chrom: std::collections::HashMap<Vec<u8>, Vec<Interval>> =
                std::collections::HashMap::new();
            for (c, iv) in wanted {
                by_chrom.entry(c).or_default().push(iv);
            }
            Some(
                by_chrom
                    .into_iter()
                    .map(|(c, ivs)| (c, RegionSet::new(ivs)))
                    .collect(),
            )
        } else {
            None
        };

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

            if let Some(regions) = &include_regions {
                match regions.get(view.chrom) {
                    Some(r) => {
                        if !r.matches(view.span, o.sites.inclusive_points) {
                            continue;
                        }
                    }
                    None => continue,
                }
            }

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
        assert_eq!(summary.records_in, 3);
        assert_eq!(summary.records_out, 1);
        let out_str = String::from_utf8(writer.0).unwrap();
        assert!(out_str.contains("rs2"));
        assert!(!out_str.contains("rs1"));
        assert!(!out_str.contains("rs3"));
    }
}
