//! Relations over files: the inputs as a command opened them, and the one
//! place that decides between streaming and per-chromosome parallelism.

use super::peek::peek_first;
use super::{RelateInputs, RelateOptions};
use crate::view::{SelectedA, ViewOptions, select_a};
use gloci_core::{ChromTable, Interval};
use gloci_engine::{SweepSummary, all_indexed};
use gloci_io::{Compression, is_regular_file, seekable::range_format};
use gloci_io::{
    Format, Header, IoError, ReaderOptions, RecordReader, RecordWriter, chrom_table_for,
};
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// Where a relation's inputs come from, and the options they were opened
/// with.
#[derive(Debug, Clone)]
pub struct InputPaths {
    pub a: String,
    pub b: Vec<String>,
    pub ro_a: ReaderOptions,
    pub ro_b: ReaderOptions,
}

/// How a relation over files runs.
#[derive(Debug, Clone)]
pub struct RunOptions {
    pub threads: NonZeroUsize,
    /// Where per-chromosome chunks spill (`--tmp-dir`).
    pub tmp_dir: PathBuf,
}

/// A relation's inputs as a command opened them: one reader per input, each
/// opened once, plus the paths, so that indexed inputs can be read per
/// chromosome (D41). Nothing here reopens a path that is not a regular file.
pub struct RelationFiles {
    pub(super) inputs: RelateInputs,
    pub(super) paths: InputPaths,
    /// The relation's chromosome table, once something was read with it.
    pub(super) chroms: Option<ChromTable>,
    b_width: Option<usize>,
    /// `--region`/`--regions-file`: the regions A is restricted to.
    pub(super) a_regions: Option<ARegions>,
}

/// The regions A is restricted to (D75).
pub(super) struct ARegions {
    /// In new-name space (after `--chrom-map`).
    pub(super) wanted: Vec<(Vec<u8>, Interval)>,
    pub(super) inclusive_points: bool,
    /// See [`crate::SelectedA::warn_unseen`].
    pub(super) warn_unseen: bool,
}

impl RelationFiles {
    /// `inputs` must have been opened from `paths` with their options.
    pub fn new(inputs: RelateInputs, paths: InputPaths) -> RelationFiles {
        RelationFiles {
            inputs,
            paths,
            chroms: None,
            b_width: None,
            a_regions: None,
        }
    }

    /// Restrict A to the records that overlap `regions` (D75), as `view`
    /// selects them: through A's index when it has one, else by scanning.
    /// B is not restricted: an A record that overlaps a region may reach
    /// past it, and `closest` must find B anywhere. With `unordered` (A may
    /// be out of order) A is scanned and keeps its order. No regions leaves
    /// the files as they are.
    pub fn restrict_a(self, regions: &ViewOptions, unordered: bool) -> Result<Self, IoError> {
        if !regions.has_regions() {
            return Ok(self);
        }
        let opts = ViewOptions {
            unordered,
            assume_sorted: false,
            genome: None,
            ..regions.clone()
        };
        let RelationFiles {
            inputs: RelateInputs { a, b },
            paths,
            chroms,
            b_width,
            ..
        } = self;
        let SelectedA {
            reader: a,
            wanted,
            warn_unseen,
        } = select_a(&paths.a, a, &opts, &paths.ro_a)?;
        Ok(RelationFiles {
            inputs: RelateInputs { a, b },
            paths,
            chroms,
            b_width,
            a_regions: (!unordered).then_some(ARegions {
                wanted,
                inclusive_points: opts.inclusive_points,
                warn_unseen,
            }),
        })
    }

    /// Open A and every B, each once.
    pub fn open(paths: InputPaths) -> Result<RelationFiles, IoError> {
        let a = gloci_io::open_reader(&paths.a, &paths.ro_a)?;
        let b = paths
            .b
            .iter()
            .map(|p| gloci_io::open_reader(p, &paths.ro_b))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(RelationFiles::new(RelateInputs { a, b }, paths))
    }

    pub fn a(&self) -> &dyn RecordReader {
        self.inputs.a.as_ref()
    }

    pub fn b_formats(&self) -> Vec<Format> {
        self.inputs.b.iter().map(|r| r.format()).collect()
    }

    pub fn b_count(&self) -> usize {
        self.inputs.b.len()
    }

    /// How many columns B file 1 has, for its missing-B placeholder (D45)
    /// and column names: fixed for VCF/BCF and GFF/GTF, else those of B's
    /// first record (3 when B has none). The first record is read ahead,
    /// once, from the reader already open, with the relation's chromosome
    /// table (D51).
    pub fn b_width(&mut self, relate: &RelateOptions) -> Result<usize, IoError> {
        if let Some(w) = self.b_width {
            return Ok(w);
        }
        let chroms = match &mut self.chroms {
            Some(c) => c,
            None => self.chroms.insert(relation_chroms(&self.inputs, relate)?),
        };
        let w = b_width(&mut self.inputs, chroms)?;
        self.b_width = Some(w);
        Ok(w)
    }

    /// Names of B file 1's columns for a TSV column-name line (D42):
    /// `b_<field>` for a VCF/BCF B, `b_1`…`b_n` otherwise.
    pub fn b_column_names(&mut self, relate: &RelateOptions) -> Result<Vec<String>, IoError> {
        let width = self.b_width(relate)?;
        Ok(match self.inputs.b.first() {
            Some(b) => b_column_names(b.as_ref(), width),
            None => numbered(width),
        })
    }
}

/// The relation's chromosome table (D51).
pub(crate) fn relation_chroms(
    inputs: &RelateInputs,
    o: &RelateOptions,
) -> Result<ChromTable, IoError> {
    let all: Vec<&dyn RecordReader> = std::iter::once(inputs.a.as_ref())
        .chain(inputs.b.iter().map(|r| r.as_ref()))
        .collect();
    chrom_table_for(o.genome.as_deref(), &all, o.chrom_order)
}

/// See [`RelationFiles::b_width`]. A BED-family B file 1 is replaced by a
/// reader that hands back the record read ahead.
pub(crate) fn b_width(
    inputs: &mut RelateInputs,
    chroms: &mut ChromTable,
) -> Result<usize, IoError> {
    let Some(first) = inputs.b.first() else {
        return Ok(3);
    };
    Ok(match (first.format(), first.header()) {
        (Format::Gff3 | Format::Gtf, _) => 9,
        (Format::Vcf | Format::Bcf, Header::Vcf(h)) => {
            if h.samples.is_empty() {
                8
            } else {
                9 + h.samples.len()
            }
        }
        (Format::Vcf | Format::Bcf, _) => 8,
        (Format::Bed | Format::BedGraph | Format::Wig, _) => {
            let b0 = inputs.b.remove(0);
            let (b0, columns) = peek_first(b0, chroms)?;
            inputs.b.insert(0, b0);
            columns.unwrap_or(3)
        }
    })
}

fn numbered(n: usize) -> Vec<String> {
    (1..=n).map(|i| format!("b_{i}")).collect()
}

fn b_column_names(b: &dyn RecordReader, width: usize) -> Vec<String> {
    match (b.format(), b.header()) {
        (Format::Vcf | Format::Bcf, Header::Vcf(h)) => {
            let mut names: Vec<String> =
                ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                    .iter()
                    .map(|n| format!("b_{n}"))
                    .collect();
            if !h.samples.is_empty() {
                names.push("b_FORMAT".to_string());
                names.extend(
                    h.samples
                        .iter()
                        .map(|s| format!("b_{}", String::from_utf8_lossy(s))),
                );
            }
            names
        }
        _ => numbered(width),
    }
}

/// How one input of a relation is read per chromosome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Source {
    /// A BGZF file with an index (D41).
    Indexed,
    /// An uncompressed BED, bedGraph or VCF regular file, split by a scan
    /// for its chromosomes' byte ranges (D76).
    Plain,
}

/// How each input (A, then each B) is read when the relation runs per
/// chromosome (D41, D76), or `None` when it streams: one thread, `--unsorted`
/// or `--sort`, or an input that is neither indexed nor a plain seekable file
/// of a format the scan handles. Standard input, `/dev/fd/N` and other
/// non-regular files, compressed files without an index, GFF and WIG stream.
/// The output is the same either way for inputs in table order; others are
/// refused alike.
pub(super) fn input_sources(
    files: &RelationFiles,
    relate: &RelateOptions,
    threads: NonZeroUsize,
) -> Option<Vec<Source>> {
    if threads.get() == 1 || relate.unsorted || relate.sort.is_some() {
        return None;
    }
    let readers =
        std::iter::once(files.inputs.a.as_ref()).chain(files.inputs.b.iter().map(|r| r.as_ref()));
    let paths = std::iter::once(&files.paths.a).chain(&files.paths.b);
    paths
        .zip(readers)
        .map(|(p, r)| {
            if all_indexed(&[p.as_str()]) {
                Some(Source::Indexed)
            } else if is_regular_file(p)
                && r.compression() == Compression::None
                && range_format(r.format())
            {
                Some(Source::Plain)
            } else {
                None
            }
        })
        .collect()
}

/// Whether a relation over these inputs runs per chromosome: see
/// [`input_sources`].
pub fn runs_per_chromosome(
    files: &RelationFiles,
    relate: &RelateOptions,
    threads: NonZeroUsize,
) -> bool {
    input_sources(files, relate, threads).is_some()
}

/// The relation step of a command: runs on one stream, or on one chunk of
/// chromosomes, and writes to `w`. A table is given when something was
/// already read with it (it must then be used).
pub(crate) type RelationStep<'a, R> = dyn Fn(
        RelateInputs,
        Option<ChromTable>,
        &mut dyn RecordWriter,
    ) -> Result<(SweepSummary, R), IoError>
    + Sync
    + 'a;

/// Run `step` over `files`: per chromosome when [`runs_per_chromosome`],
/// else once over the streams. Warns about records skipped under
/// `--skip-invalid` (every chunk's, together). Returns the summary and each
/// run's result, in chromosome order.
pub(crate) fn run_files<R: Send>(
    files: RelationFiles,
    relate: &RelateOptions,
    run: &RunOptions,
    w: &mut dyn RecordWriter,
    step: &RelationStep<R>,
) -> Result<(SweepSummary, Vec<R>), IoError> {
    let (summary, rs) = if let Some(sources) = input_sources(&files, relate, run.threads) {
        super::indexed::run_per_chromosome(files, &sources, relate, run, w, step)?
    } else {
        let RelationFiles { inputs, chroms, .. } = files;
        let (summary, r) = step(inputs, chroms, w)?;
        (summary, vec![r])
    };
    super::report_skipped(&summary);
    Ok((summary, rs))
}
