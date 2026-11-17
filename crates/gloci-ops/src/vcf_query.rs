//! `gloci vcf query`: the records of the shared selection stage written through
//! a format string (decision D172, after bcftools 1.21 `query`).

use crate::vcf_select::{
    Lines, SelectStage, SelectedInput, VcfSelect, kept_indices, run_selection,
};
use gloci_io::{IoError, Loc};
use gloci_vcf::{QueryFormat, RewritePlan, SampleCols};
use std::collections::HashMap;
use std::io::Write;
use std::num::NonZeroUsize;

const COMMAND: &str = "vcf query";

/// Options for `vcf query`.
#[derive(Clone, Debug)]
pub struct QueryOptions {
    pub select: VcfSelect,
    /// Threads the selection runs on (`None` = 1).
    pub threads: Option<NonZeroUsize>,
    pub format: QueryFormat,
    /// `-H`.
    pub print_header: bool,
}

/// Counts from a query run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QuerySummary {
    pub records_in: u64,
    pub sites_out: u64,
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".to_string(),
        source,
    }
}

/// Positions in `kept` in the order the samples print: the order `--samples` gave
/// them (a name excluded or repeated is skipped), else `0..kept.len()` (D172).
fn print_order(sel: &VcfSelect, names: &[Vec<u8>], kept: &[usize]) -> Vec<usize> {
    let Some(given) = sel.include_samples.as_ref() else {
        return (0..kept.len()).collect();
    };
    let mut at: HashMap<&[u8], usize> = HashMap::with_capacity(kept.len());
    for (pos, &h) in kept.iter().enumerate() {
        at.entry(names[h].as_slice()).or_insert(pos);
    }
    let mut seen = vec![false; kept.len()];
    let mut order = Vec::with_capacity(kept.len());
    for name in given {
        if let Some(&pos) = at.get(name.as_bytes())
            && !std::mem::replace(&mut seen[pos], true)
        {
            order.push(pos);
        }
    }
    order
}

/// Write each selected record of `input` through the format to `w`, after the -H line when asked.
/// Write errors name the placeholder path `output`.
pub fn vcf_query(
    input: SelectedInput,
    o: &QueryOptions,
    w: &mut dyn Write,
) -> Result<QuerySummary, IoError> {
    o.select.validate()?;
    let base = RewritePlan::default();
    let (reader, mut chroms) = input.into_parts();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), &o.select, &base, false)?;
    stage.set_threads(o.threads.unwrap_or(NonZeroUsize::MIN));
    stage.filter_samples();
    let name = stage.name().to_string();
    let query = o.format.resolve(stage.header(), o.print_header)?;
    let kept = kept_indices(&stage);
    let header_samples = stage.header().samples.clone();
    let order = print_order(&o.select, &header_samples, &kept);

    if o.print_header {
        let names: Vec<&[u8]> = order
            .iter()
            .map(|&p| header_samples[kept[p]].as_slice())
            .collect();
        let mut head = Vec::new();
        query.write_header(&names, stage.expr_reads_samples(), &mut head);
        w.write_all(&head).map_err(write_err)?;
    }

    let mut lines = Lines::for_stage(&stage)?;
    let mut out = Vec::new();
    let mut sites_out = 0u64;
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        let line = lines.line(&name, &hit)?;
        let cols = SampleCols {
            kept: &kept,
            order: &order,
            names: &header_samples,
            pass: hit.sample_pass,
        };
        out.clear();
        query
            .write_record(line, &cols, &mut out)
            .map_err(|msg| IoError::Parse {
                path: name.clone(),
                loc: Loc::Record(hit.rec_no),
                msg,
            })?;
        w.write_all(&out).map_err(write_err)?;
        sites_out += 1;
        Ok(())
    })?;
    Ok(QuerySummary {
        records_in,
        sites_out,
    })
}
