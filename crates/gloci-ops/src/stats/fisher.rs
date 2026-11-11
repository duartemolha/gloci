//! `stats fisher`: bedtools' 2x2 table of interval overlaps and klib's exact
//! test on it (replaces `bedtools fisher`).

use super::merged::{DropUnknown, MergedReader, Tally};
use crate::relate::{RelateInputs, RelateOptions, relation_chroms, run_relation_in};
use crate::skipped::{Skipped, log_skipped};
use crate::sort::{SortSettings, sort_readers};
use gloci_core::{ChromOrder, ChromTable, Genome, OverlapCriteria, StrandRule};
use gloci_engine::{ARec, Hit, Neighbours, Visitor};
use gloci_io::{
    InvalidPolicy, IoError, ReaderOptions, RecordBatch, RecordReader, UsageError, open_reader,
};
use std::sync::atomic::Ordering;

#[derive(Debug, Clone)]
pub struct FisherOptions {
    pub criteria: OverlapCriteria,
    /// `--merge-overlaps`: merge each input first (per strand under `-s`/`-S`).
    pub merge: bool,
    /// `--exclude BED`: bases of the genome to leave out of the estimate.
    pub exclude: Option<String>,
    /// `--sort`: sort A, B and the exclude file first.
    pub sort: Option<SortSettings>,
    /// The raw `-g` value, for the chromosome table.
    pub genome: Option<String>,
}

/// What the sweep counted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FisherCounts {
    pub n_a: u64,
    pub n_b: u64,
    /// Accepted (A, B) pairs.
    pub overlaps: u64,
    pub len_a: u64,
    pub len_b: u64,
    /// Sum of the `-g` lengths less the merged length of the exclude file.
    pub genome_size: u64,
}

/// The 2x2 table as bedtools builds it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FisherTable {
    pub n11: i64,
    pub n12: i64,
    pub n21: i64,
    pub n22: i64,
    /// The estimated number of possible intervals.
    pub n22_full: i64,
}

struct Pairs<'a> {
    criteria: &'a OverlapCriteria,
    pairs: u64,
}

impl Visitor for Pairs<'_> {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        _near: &Neighbours<'_>,
    ) -> Result<(), IoError> {
        for h in hits {
            if self
                .criteria
                .accepts(a.meta.span, a.meta.strand, h.span, h.strand)
            {
                self.pairs += 1;
            }
        }
        Ok(())
    }
}

/// The merged length of the intervals of `path`: overlapping and bookended
/// records count once. The file must be sorted (or `--sort`ed).
fn merged_length(path: &str, sort: Option<&SortSettings>) -> Result<u64, IoError> {
    let ro = ReaderOptions {
        threads: std::num::NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
        chrom_map: None,
        accept: gloci_io::Accept::default(),
    };
    let r = open_reader(path, &ro)?;
    let (r, mut chroms): (Box<dyn RecordReader>, ChromTable) = match sort {
        Some(s) => {
            let mut c = ChromTable::lexicographic();
            let r = sort_readers(vec![r], &mut c, s).map_err(IoError::with_sorted_input_loc)?;
            (r, c)
        }
        None => (r, ChromTable::first_appearance()),
    };
    let mut merged = MergedReader::new(r, false);
    let mut batch = RecordBatch::default();
    let mut total = 0u64;
    loop {
        batch.clear();
        let more = merged.read_batch(&mut chroms, &mut batch, 4096)?;
        total = total.saturating_add(batch.recs.iter().map(|m| m.span.end - m.span.start).sum());
        if !more {
            return Ok(total);
        }
    }
}

/// Count the records, bases and accepted pairs of A and B (merged first with
/// `merge`), and size the genome: the `-g` lengths less the merged length of
/// the exclude file. Both inputs are read to the end.
pub fn fisher_counts(
    inputs: RelateInputs,
    genome: &Genome,
    o: &FisherOptions,
) -> Result<FisherCounts, IoError> {
    let RelateInputs { a, b } = inputs;
    if b.len() != 1 {
        return Err(UsageError(format!(
            "fisher compares two files, A and one B, but {} B files were given; \
run it once per B file",
            b.len()
        ))
        .into());
    }
    let mut genome_size = 0u64;
    for i in 0..genome.len() {
        let Some(len) = genome.length(i) else {
            return Err(UsageError(format!(
                "-g has no length for {}; fisher sizes the genome from the chromosome lengths, \
so give a .genome or .fai file, a .dict, a BAM header or a built-in genome such as hg38",
                String::from_utf8_lossy(genome.name(i))
            ))
            .into());
        };
        genome_size = genome_size.saturating_add(len);
    }
    if let Some(path) = &o.exclude {
        genome_size = genome_size.saturating_sub(merged_length(path, o.sort.as_ref())?);
    }
    let relate = RelateOptions {
        criteria: o.criteria,
        window: None,
        split: false,
        split_deletions: false,
        unsorted: false,
        sort: None,
        labels: Vec::new(),
        genome: o.genome.clone(),
        chrom_order: ChromOrder::default(),
        header: false,
    };
    let inputs = RelateInputs { a, b };
    let mut chroms = relation_chroms(&inputs, &relate)?;
    let RelateInputs { a, b } = inputs;
    let by_strand = o.criteria.strand != StrandRule::Any;
    let mut skipped = Skipped::default();
    let prepare = |r: Box<dyn RecordReader>,
                   chroms: &mut ChromTable,
                   skipped: &mut Skipped|
     -> Result<Box<dyn RecordReader>, IoError> {
        match &o.sort {
            Some(s) => sort_readers(vec![r], chroms, s).map_err(IoError::with_sorted_input_loc),
            None => Ok(skipped.watch(r)),
        }
    };
    let a = prepare(a, &mut chroms, &mut skipped)?;
    let b = prepare(
        b.into_iter().next().expect("one B"),
        &mut chroms,
        &mut skipped,
    )?;
    let (name_a, name_b) = (a.name().to_string(), b.name().to_string());
    // Merged clusters of unknown strand can never match under -s/-S, so they
    // take no part (bedtools leaves them out too). Unmerged records of
    // unknown strand are still counted, as bedtools counts them.
    let wrap = |r: Box<dyn RecordReader>| -> (Box<dyn RecordReader>, _) {
        if !o.merge {
            return (r, None);
        }
        let merged = Box::new(MergedReader::new(r, by_strand));
        if by_strand {
            let (d, flag) = DropUnknown::new(merged);
            (Box::new(d), Some(flag))
        } else {
            (merged, None)
        }
    };
    let (a, flag_a) = wrap(a);
    let (b, flag_b) = wrap(b);
    let (a, a_counts) = Tally::new(a);
    let (b, b_counts) = Tally::new(b);
    let mut sum = Pairs {
        criteria: &o.criteria,
        pairs: 0,
    };
    run_relation_in(
        RelateInputs {
            a: Box::new(a),
            b: vec![Box::new(b)],
        },
        Some(chroms),
        &relate,
        None,
        &mut sum,
    )?;
    for (name, flag) in [(name_a, flag_a), (name_b, flag_b)] {
        if flag.is_some_and(|f| f.load(Ordering::Relaxed)) {
            log::warn!(
                "{name} has records with unknown strand; under -s and -S merged intervals of \
unknown strand are left out of the table (they can never match)"
            );
        }
    }
    log_skipped(&skipped.take());
    skipped.log_alignment_counts();
    Ok(FisherCounts {
        n_a: a_counts.records.load(Ordering::Relaxed),
        n_b: b_counts.records.load(Ordering::Relaxed),
        overlaps: sum.pairs,
        len_a: a_counts.bases.load(Ordering::Relaxed),
        len_b: b_counts.bases.load(Ordering::Relaxed),
        genome_size,
    })
}

fn to_i64(v: u64) -> i64 {
    i64::try_from(v).unwrap_or(i64::MAX)
}

/// `floor(G / (2 + LA/nA + LB/nB))` for nA, nB > 0: exact in u128 as
/// `G*nA*nB / (2*nA*nB + LB*nA + LA*nB)`, in f64 when that overflows.
fn estimate(c: &FisherCounts) -> i64 {
    let (g, na, nb) = (
        u128::from(c.genome_size),
        u128::from(c.n_a),
        u128::from(c.n_b),
    );
    let (la, lb) = (u128::from(c.len_a), u128::from(c.len_b));
    let exact = g
        .checked_mul(na)
        .and_then(|x| x.checked_mul(nb))
        .and_then(|num| {
            let den = 2u128
                .checked_mul(na.checked_mul(nb)?)?
                .checked_add(lb.checked_mul(na)?)?
                .checked_add(la.checked_mul(nb)?)?;
            num.checked_div(den)
        });
    match exact {
        Some(v) => i64::try_from(v).unwrap_or(i64::MAX),
        None => {
            let d = 2.0 + c.len_a as f64 / c.n_a as f64 + c.len_b as f64 / c.n_b as f64;
            (c.genome_size as f64 / d).floor() as i64
        }
    }
}

/// bedtools' table: n12 and n21 are what is left of A and B, and n22 what is
/// left of the estimated number of possible intervals. An empty A or B makes
/// bedtools' estimate NaN, which C converts to a value that loses the `max`
/// against `n11 + n12 + n21` on x86 (see differences.md#fisher-empty-estimate).
pub fn fisher_table(c: &FisherCounts) -> FisherTable {
    let n11 = to_i64(c.overlaps);
    let n12 = (to_i64(c.n_a) - n11).max(0);
    let n21 = (to_i64(c.n_b) - n11).max(0);
    let floor = n11.saturating_add(n12).saturating_add(n21);
    let n22_full = if c.n_a == 0 || c.n_b == 0 {
        floor
    } else {
        floor.max(estimate(c))
    };
    let n22 = n22_full
        .saturating_sub(n12)
        .saturating_sub(n21)
        .saturating_sub(n11)
        .max(0);
    FisherTable {
        n11,
        n12,
        n21,
        n22,
        n22_full,
    }
}

// ---- klib's kt_fisher_exact as bedtools 2.31.1 bundles it (src/fisher/kfunc.cpp) ----

fn lbinom(n: i64, k: i64) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    libm::lgamma((n + 1) as f64) - libm::lgamma((k + 1) as f64) - libm::lgamma((n - k + 1) as f64)
}

fn hypergeo(n11: i64, n1_: i64, n_1: i64, n: i64) -> f64 {
    libm::exp(lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1))
}

/// klib's `hgacc_t`: the incremental hypergeometric probability.
#[derive(Default)]
struct HgAcc {
    n11: i64,
    n1_: i64,
    n_1: i64,
    n: i64,
    p: f64,
}

impl HgAcc {
    /// `hypergeo_acc(n11, n1_, n_1, n, aux)`; the all-zero margins form
    /// (`n1_ == n_1 == n == 0`) means "only n11 changed".
    fn acc(&mut self, n11: i64, n1_: i64, n_1: i64, n: i64) -> f64 {
        if n1_ != 0 || n_1 != 0 || n != 0 {
            self.n11 = n11;
            self.n1_ = n1_;
            self.n_1 = n_1;
            self.n = n;
        } else {
            if n11 % 11 != 0 && n11 + self.n - self.n1_ - self.n_1 != 0 {
                if n11 == self.n11 + 1 {
                    self.p *= (self.n1_ - self.n11) as f64 / n11 as f64
                        * (self.n_1 - self.n11) as f64
                        / (n11 + self.n - self.n1_ - self.n_1) as f64;
                    self.n11 = n11;
                    return self.p;
                }
                if n11 == self.n11 - 1 {
                    self.p *= self.n11 as f64 / (self.n1_ - n11) as f64
                        * (self.n11 + self.n - self.n1_ - self.n_1) as f64
                        / (self.n_1 - n11) as f64;
                    self.n11 = n11;
                    return self.p;
                }
            }
            self.n11 = n11;
        }
        self.p = hypergeo(self.n11, self.n1_, self.n_1, self.n);
        self.p
    }
}

/// klib's `kt_fisher_exact`: (left, right, two-tail) p-values of the table.
pub fn kt_fisher_exact(n11: i64, n12: i64, n21: i64, n22: i64) -> (f64, f64, f64) {
    let n1_ = n11 + n12;
    let n_1 = n11 + n21;
    let n = n11 + n12 + n21 + n22;
    let max = n_1.min(n1_);
    let min = (n1_ + n_1 - n).max(0);
    if min == max {
        return (1.0, 1.0, 1.0);
    }
    let mut aux = HgAcc::default();
    let q = aux.acc(n11, n1_, n_1, n);
    if q == 0.0 {
        // The table's probability underflows, so both tails and the two-tail are 0 to a good
        // approximation; which of left and right gets the 1 depends on which side of the mode
        // n11 sits (bedtools 2.31.1 src/fisher/kfunc.cpp, compared without a division).
        return if n11.wrapping_mul(n.wrapping_add(2))
            < n_1.wrapping_add(1).wrapping_mul(n1_.wrapping_add(1))
        {
            (0.0, 1.0, 0.0)
        } else {
            (1.0, 0.0, 0.0)
        };
    }
    // left tail
    let mut p = aux.acc(min, 0, 0, 0);
    let mut left = 0.0;
    let mut i = min + 1;
    while p < 0.999_999_99 * q && i <= max {
        left += p;
        p = aux.acc(i, 0, 0, 0);
        i += 1;
    }
    i -= 1;
    if p < 1.000_000_01 * q {
        left += p;
    } else {
        i -= 1;
    }
    // right tail
    p = aux.acc(max, 0, 0, 0);
    let mut right = 0.0;
    let mut j = max - 1;
    while p < 0.999_999_99 * q && j >= 0 {
        right += p;
        p = aux.acc(j, 0, 0, 0);
        j -= 1;
    }
    j += 1;
    if p < 1.000_000_01 * q {
        right += p;
    } else {
        j += 1;
    }
    let mut two = left + right;
    if two > 1.0 {
        two = 1.0;
    }
    if (i - n11).abs() < (j - n11).abs() {
        right = 1.0 - left + q;
    } else {
        left = 1.0 - right + q;
    }
    (left, right, two)
}

fn cell(out: &mut Vec<u8>, v: i64) {
    let start = out.len();
    gloci_core::cfmt::write_i64(out, v);
    while out.len() - start < 12 {
        out.push(b' ');
    }
}

/// bedtools' report: the `#` lines, then the p-values and the odds ratio.
pub fn write_fisher(c: &FisherCounts, t: &FisherTable, out: &mut Vec<u8>) {
    use gloci_core::cfmt::{write_fixed, write_g, write_i64, write_u64};
    out.extend_from_slice(b"# Number of query intervals: ");
    write_u64(out, c.n_a);
    out.extend_from_slice(b"\n# Number of db intervals: ");
    write_u64(out, c.n_b);
    out.extend_from_slice(b"\n# Number of overlaps: ");
    write_u64(out, c.overlaps);
    out.extend_from_slice(b"\n# Number of possible intervals (estimated): ");
    write_i64(out, t.n22_full);
    out.extend_from_slice(b"\n# phyper(");
    write_i64(out, t.n11);
    out.extend_from_slice(b" - 1, ");
    write_u64(out, c.n_a);
    out.extend_from_slice(b", ");
    write_i64(out, t.n22_full);
    out.extend_from_slice(b" - ");
    write_u64(out, c.n_a);
    out.extend_from_slice(b", ");
    write_u64(out, c.n_b);
    out.extend_from_slice(b", lower.tail=F)\n# Contingency Table Of Counts\n");
    out.extend_from_slice(b"#_________________________________________\n");
    out.extend_from_slice(b"#           |  in -b       | not in -b    |\n");
    out.extend_from_slice(b"#     in -a | ");
    cell(out, t.n11);
    out.extend_from_slice(b" | ");
    cell(out, t.n12);
    out.extend_from_slice(b" |\n# not in -a | ");
    cell(out, t.n21);
    out.extend_from_slice(b" | ");
    cell(out, t.n22);
    out.extend_from_slice(b" |\n#_________________________________________\n");
    out.extend_from_slice(b"# p-values for fisher's exact test\n");
    out.extend_from_slice(b"left\tright\ttwo-tail\tratio\n");
    let (left, right, two) = kt_fisher_exact(t.n11, t.n12, t.n21, t.n22);
    for p in [left, right, two] {
        write_g(out, p, 5);
        out.push(b'\t');
    }
    let ratio = (t.n11 as f64 / t.n12 as f64) / (t.n21 as f64 / t.n22 as f64);
    if ratio.is_nan() {
        out.extend_from_slice(b"-nan");
    } else {
        write_fixed(out, ratio, 3);
    }
    out.push(b'\n');
}
