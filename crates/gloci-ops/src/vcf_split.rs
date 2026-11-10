//! `gloci vcf split`: multi-allelic records to biallelic ones.

use crate::vcf_select::vcf_header;
use gloci_core::ChromTable;
use gloci_io::{
    BcfText, Format, Header, IoError, Loc, OutFormat, RecordBatch, RecordReader, RecordWriter,
    TextToBcf, UsageError, VcfHeader, next_batch_size,
};
use gloci_vcf::{SplitOutcome, SplitTypes, Splitter};

const COMMAND: &str = "vcf split";

/// Options for `vcf split`.
#[derive(Clone, Copy, Debug)]
pub struct SplitOptions {
    pub types: SplitTypes,
    /// Skip a record that cannot be split, with a warning, instead of failing.
    pub skip_invalid: bool,
}

/// Records read, split and written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SplitSummary {
    pub records_in: u64,
    pub records_split: u64,
    pub records_out: u64,
    pub records_skipped: u64,
}

/// The output format (`requested`, else the input's) and header: the input's.
pub fn vcf_split_output(
    input: &dyn RecordReader,
    requested: Option<OutFormat>,
) -> Result<(OutFormat, Header), IoError> {
    vcf_rewrite_output(COMMAND, input, requested)
}

/// [`vcf_split_output`] for any command that rewrites VCF/BCF records and
/// keeps the input's header (`command` names it in errors).
pub fn vcf_rewrite_output(
    command: &str,
    input: &dyn RecordReader,
    requested: Option<OutFormat>,
) -> Result<(OutFormat, Header), IoError> {
    let h = vcf_header(command, input)?;
    let to = requested.unwrap_or(OutFormat::family(input.format()));
    if to != OutFormat::Vcf && to != OutFormat::Bcf {
        return Err(UsageError(format!(
            "{command} writes VCF or BCF, not {}; name the output .vcf, .vcf.gz or .bcf, \
or pass -O vcf or -O bcf",
            to.name()
        ))
        .into());
    }
    let header = VcfHeader::parse(h.raw).map_err(|e| IoError::Detect {
        path: input.name().to_string(),
        msg: e,
    })?;
    Ok((to, Header::Vcf(header)))
}

/// Split every multi-allelic record of `input` and write the records to `w`
/// (opened with [`vcf_split_output`]'s format and header), in input order.
pub fn vcf_split(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &SplitOptions,
    w: &mut dyn RecordWriter,
) -> Result<SplitSummary, IoError> {
    let h = vcf_header(COMMAND, input.as_ref())?;
    let name = input.name().to_string();
    let splitter = Splitter::new(&h.raw, o.types, h.samples.len());
    let from_bcf = input.format() == Format::Bcf;
    let to_bcf = w.is_bcf();
    let mut bcf_text = if from_bcf {
        Some(BcfText::for_header(&name, &h)?)
    } else {
        None
    };
    let mut text_to_bcf = if to_bcf {
        Some(TextToBcf::new(&name, &h)?)
    } else {
        None
    };
    let write_err = |source| IoError::Write {
        path: "output".to_string(),
        source,
    };
    let mut sum = SplitSummary::default();
    let (mut out, mut ends, mut enc) = (Vec::new(), Vec::new(), Vec::new());
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        for i in 0..batch.len() {
            let raw = batch.raw(i);
            let loc = batch.recs[i].loc;
            sum.records_in += 1;
            let text = match bcf_text.as_mut() {
                Some(t) => t.render(raw, loc)?,
                None => raw,
            };
            out.clear();
            let outcome = match splitter.split(text, &mut out, &mut ends) {
                Ok(x) => x,
                Err(msg) => {
                    if o.skip_invalid {
                        sum.records_skipped += 1;
                        log::warn!("{COMMAND}: {name}: record {loc} skipped: {msg}");
                        continue;
                    }
                    let l = if from_bcf {
                        Loc::Record(loc)
                    } else {
                        Loc::Line(loc)
                    };
                    return Err(IoError::Parse {
                        path: name.clone(),
                        loc: l,
                        msg: format!("{msg} (pass --skip-invalid to skip such records)"),
                    });
                }
            };
            let mut emit = |line: &[u8]| -> Result<(), IoError> {
                if let Some(t) = text_to_bcf.as_mut() {
                    enc.clear();
                    t.encode(line, loc, &mut enc)?;
                    w.write_record(&enc).map_err(write_err)
                } else {
                    w.write_record(line).map_err(write_err)
                }
            };
            match outcome {
                SplitOutcome::Unchanged => {
                    emit(text)?;
                    sum.records_out += 1;
                }
                SplitOutcome::Split(n) => {
                    let mut prev = 0;
                    for &e in &ends {
                        emit(&out[prev..e])?;
                        prev = e;
                    }
                    sum.records_split += 1;
                    sum.records_out += n as u64;
                }
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    if let Some(s) = input.diagnostics().summary(&name) {
        log::warn!("{s}");
    }
    Ok(sum)
}
