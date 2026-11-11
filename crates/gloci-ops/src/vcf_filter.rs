//! `gloci vcf filter`: the selection pipeline (spec §7.8, decision D50).

use crate::vcf_select::{
    SelectStage, SelectedInput, VcfSelect, open_selected, parse_select_expr, rewrite_plan,
    run_selection, select_regions, vcf_header,
};
use gloci_core::ChromTable;
use gloci_io::{
    BcfText, Header, IoError, Loc, OutFormat, ReaderOptions, RecordReader, RecordWriter, TextToBcf,
    UsageError, VcfHeader,
};
use gloci_vcf::{RewritePlan, rewrite_header, rewrite_text};
use std::num::NonZeroUsize;

const COMMAND: &str = "vcf filter";

/// Options for `vcf filter`: the selection shared with the other `vcf`
/// commands, and how the passing records are rewritten.
#[derive(Clone, Debug, Default)]
pub struct VcfFilterOptions {
    pub select: VcfSelect,
    /// Threads the selection runs on (`None` = 1); see `SelectStage::set_threads`.
    pub threads: Option<NonZeroUsize>,
    /// INFO pruning and AC recalculation; `selection` and `mask` come from
    /// `select`.
    pub plan: RewritePlan,
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
        self.select.validate()
    }
}

/// What `vcf filter` writes for `input`: the output format (`requested`,
/// else the input's) and the header (samples, INFO lines and so on as the
/// rewrite leaves them).
pub fn vcf_filter_output(
    input: &dyn RecordReader,
    requested: Option<OutFormat>,
    o: &VcfFilterOptions,
) -> Result<(OutFormat, Header), IoError> {
    let h = vcf_header(COMMAND, input)?;
    let to = requested.unwrap_or(OutFormat::family(input.format()));
    if to != OutFormat::Vcf && to != OutFormat::Bcf {
        return Err(usage(format!(
            "vcf filter writes VCF or BCF, not {}; name the output .vcf, .vcf.gz or .bcf, \
or pass -O vcf or -O bcf",
            to.name()
        )));
    }
    let plan = rewrite_plan(&h, &o.select, &o.plan)?;
    // A bad expression is a usage error before any record is read.
    parse_select_expr(&h, &o.select)?;
    let header = VcfHeader::parse(rewrite_header(&h, &plan)).map_err(|e| IoError::Detect {
        path: input.name().to_string(),
        msg: e,
    })?;
    Ok((to, Header::Vcf(header)))
}

/// An input opened by [`open_vcf_filter`], its regions already selected.
pub type VcfFilterInput = SelectedInput;

/// Open `path` for `vcf filter`. With `--region`/`--regions-file`, only the
/// records that overlap a region are read: through the index when the file
/// has one, else by scanning it (as `view` does).
pub fn open_vcf_filter(
    path: &str,
    ro: &ReaderOptions,
    o: &VcfFilterOptions,
) -> Result<VcfFilterInput, IoError> {
    open_selected(path, ro, &o.select)
}

/// Filter an input opened by [`open_vcf_filter`] and write the passing
/// records to `w` (opened with [`vcf_filter_output`]'s format and header).
pub fn vcf_filter_file(
    input: VcfFilterInput,
    o: &VcfFilterOptions,
    w: &mut dyn RecordWriter,
) -> Result<FilterSummary, IoError> {
    let (reader, mut chroms) = input.into_parts();
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
    let input = select_regions(COMMAND, input, chroms, &o.select)?;
    filter(input, chroms, o, w)
}

/// The pipeline: the shared selection, then rewriting.
fn filter(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &VcfFilterOptions,
    w: &mut dyn RecordWriter,
) -> Result<FilterSummary, IoError> {
    let mut stage = SelectStage::new(
        COMMAND,
        input.as_ref(),
        &o.select,
        &o.plan,
        o.plan.recalc_ac,
    )?;
    stage.set_threads(o.threads.unwrap_or(NonZeroUsize::MIN));
    let plan = stage.plan().clone();
    let input_is_bcf = stage.input_is_bcf();
    let output_is_bcf = w.is_bcf();
    let name = stage.name().to_string();

    let mut text_to_bcf = if output_is_bcf && (!input_is_bcf || !plan.is_passthrough()) {
        let rewritten_header_raw = rewrite_header(stage.header(), &plan);
        let out_vcf_header =
            VcfHeader::parse(rewritten_header_raw).map_err(|e| IoError::Detect {
                path: name.clone(),
                msg: e,
            })?;
        Some(TextToBcf::new(&name, &out_vcf_header)?)
    } else {
        None
    };
    let mut bcf_text = if input_is_bcf {
        Some(BcfText::for_header(&name, stage.header())?)
    } else {
        None
    };

    let mut rewritten_buf = Vec::new();
    let mut bcf_enc_buf = Vec::new();
    let mut records_out = 0u64;
    let write_err = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };

    let records_in = run_selection(input, chroms, &mut stage, |hit| {
        let raw = hit.raw;
        // Output record
        if input_is_bcf == output_is_bcf && plan.is_passthrough() {
            w.write_record(raw).map_err(write_err)?;
        } else {
            let text_line = if input_is_bcf {
                bcf_text.as_mut().unwrap().render(raw, hit.rec_no)?
            } else {
                raw
            };

            let out_text = if !plan.is_passthrough() {
                rewritten_buf.clear();
                rewrite_text(text_line, &plan, hit.masked, hit.stats, &mut rewritten_buf).map_err(
                    |msg| IoError::Parse {
                        path: name.clone(),
                        loc: Loc::Record(hit.rec_no),
                        msg,
                    },
                )?;
                &rewritten_buf[..]
            } else {
                text_line
            };

            if output_is_bcf {
                bcf_enc_buf.clear();
                text_to_bcf
                    .as_mut()
                    .unwrap()
                    .encode(out_text, hit.rec_no, &mut bcf_enc_buf)?;
                w.write_record(&bcf_enc_buf).map_err(write_err)?;
            } else {
                w.write_record(out_text).map_err(write_err)?;
            }
        }
        records_out += 1;
        Ok(())
    })?;

    Ok(FilterSummary {
        records_in,
        records_out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_core::VcfSpan;
    use gloci_io::{Format, open_reader};

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
            accept: gloci_io::Accept::default(),
            chrom_map: None,
            format: Some(Format::Vcf),
            vcf_span: VcfSpan::Ref,
            invalid: gloci_io::InvalidPolicy::Fail,
            threads: std::num::NonZero::new(1).unwrap(),
        };
        let reader = open_reader(temp.path().to_str().unwrap(), &ro).unwrap();
        let mut chroms = ChromTable::first_appearance();
        let opts = VcfFilterOptions {
            select: VcfSelect {
                regions: vec!["chr1:150-250".to_string()],
                ..Default::default()
            },
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
