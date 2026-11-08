//! `gloci merge`: cluster the footprints of sorted inputs.
use crate::Ctx;
use crate::cli::{MergeArgs, VcfSpan};
use gloci_io::{
    Format, Header, InvalidPolicy, OutFormat, OutputOptions, ReaderOptions, RecordReader,
    TsvHeader, UsageError, WriterOptions, chrom_table_for, open_output, open_reader,
    open_record_writer, open_tsv_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::agg::parse_aggs;
use gloci_ops::merge::{MergeOptions, merge, validate};
use std::io::IsTerminal;

pub fn run(ctx: &Ctx, a: &MergeArgs) -> anyhow::Result<()> {
    let aggs = parse_aggs(&a.columns)?;
    let o = MergeOptions {
        max_gap: a.max_gap,
        per_strand: a.strand,
        count: a.count,
        aggs,
        empty: a.empty.clone(),
        delim: a.delim.clone(),
    };
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let ro = ReaderOptions {
        threads: ctx.threads,
        format: match a.format.format.as_deref() {
            None => None,
            Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: match a.vcf_span.vcf_span {
            VcfSpan::Allele => gloci_core::VcfSpan::Allele,
            VcfSpan::Ref => gloci_core::VcfSpan::Ref,
        },
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let readers = a
        .files
        .iter()
        .map(|f| open_reader(f, &ro))
        .collect::<Result<Vec<_>, _>>()?;
    let formats: Vec<Format> = readers.iter().map(|r| r.format()).collect();
    validate(&formats, &o)?;
    let refs: Vec<&dyn RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let mut chroms = chrom_table_for(a.genome.genome.as_deref(), &refs, a.chrom_order.chrom_order)?;

    // The rows are new BED lines whatever the input format.
    let to = spec.format.unwrap_or(OutFormat::Bed);
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        to,
        spec.bgzf,
    )?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let out = open_output(&spec, to, &oo)?;
    let mut w = if to == OutFormat::Tsv {
        open_tsv_writer(out, Format::Bed, &[], TsvHeader::None, oo.write_index)
    } else {
        let wo = WriterOptions {
            write_header: false,
            command_line: Some(ctx.command_line.clone()),
            version: ctx.version,
        };
        open_record_writer(out, Format::Bed, to, &Header::None, &wo)
    }
    .map_err(|e| e.with_output_name(&dest))?;
    merge(readers, &mut chroms, &o, w.as_mut()).map_err(|e| e.with_output_name(&dest))?;
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
