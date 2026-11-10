//! Opening the writer for a relation report.

use gloci_io::{
    Format, Header, IoError, OutFormat, OutputOptions, OutputSpec, RecordReader, RecordWriter,
    TsvHeader, UsageError, WriterOptions, open_output, open_record_writer, open_tsv_writer,
    refuse_binary_to_tty, tsv_output_ok,
};
use std::io::IsTerminal;

/// What a report looks like.
#[derive(Debug, Clone)]
pub enum ReportOutput {
    /// A's records, in A's format (`intersect` without `--mode pairs`).
    SameAsA,
    /// A's records with these columns appended, as TSV (D42).
    Tsv { appended: Vec<String> },
}

/// Open the destination and write the header. `header` asks for BED/GFF
/// header lines (`--header`); a VCF always gets its header (D42).
pub fn open_report(
    a: &dyn RecordReader,
    kind: &ReportOutput,
    spec: &OutputSpec,
    oo: &OutputOptions,
    wo: &WriterOptions,
    header: bool,
) -> Result<Box<dyn RecordWriter>, IoError> {
    let a_format = a.format();
    let to_tty = spec.path.is_none() && std::io::stdout().is_terminal();
    match kind {
        ReportOutput::SameAsA => {
            let to = spec.format.unwrap_or(OutFormat::family(a_format));
            if a_format.is_alignment() && !matches!(to, OutFormat::Bam | OutFormat::Sam) {
                return Err(UsageError(format!(
                    "A is {} alignments, which are written as BAM or SAM, not {}; use -O bam \
or -O sam (-o x.bam, -o x.sam), or for BED write BAM and convert: \
gloci intersect -a reads.bam -b targets.bed -O bam | gloci convert - -o x.bed",
                    a_format.name(),
                    to.name()
                ))
                .into());
            }
            refuse_binary_to_tty(to_tty, to, spec.bgzf)?;
            let out = open_output(spec, to, oo)?;
            let wo = WriterOptions {
                // Alignments always carry their header (BAM cannot do without it).
                write_header: header || a_format.is_alignment(),
                ..wo.clone()
            };
            open_record_writer(out, a_format, to, a.header(), &wo)
        }
        ReportOutput::Tsv { appended } => {
            let to = spec.format.unwrap_or(OutFormat::Tsv);
            tsv_output_ok(a_format, to)?;
            let tsv_header = match (a.header(), a_format) {
                (Header::Vcf(h), Format::Vcf | Format::Bcf) => TsvHeader::VcfColumns {
                    samples: h.samples.clone(),
                },
                (Header::Lines(l), _) if header => TsvHeader::Lines(l.clone()),
                (_, Format::Vcf | Format::Bcf) => {
                    return Err(IoError::Internal(
                        "a VCF/BCF input has no VCF header".into(),
                    ));
                }
                _ => TsvHeader::None,
            };
            let out = open_output(spec, to, oo)?;
            open_tsv_writer(out, a_format, appended, tsv_header, oo.write_index)
        }
    }
}
