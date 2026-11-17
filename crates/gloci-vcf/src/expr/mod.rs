//! bcftools-style filter expressions (`--include` / `--exclude`, D71): a
//! supported subset of bcftools `-i`/`-e`.
//!
//! Parse once per header with [`parse`], then call [`Expr::eval`] per record.
//! Evaluation does not allocate once the [`EvalScratch`] buffers have grown to
//! the record shape; regexes are compiled at parse time.
//!
//! Wiring (see `vcf filter`): after masking, build the genotypes with
//! [`MaskedGenotypes`] (masked samples read as missing) and the [`SiteStats`]
//! from the same masking, so `FMT/` fields, `GT`, `MAF`, `AC`, `AN`,
//! `F_MISSING` and `N_MISSING` all see the masked data.

mod eval;
mod lexer;
mod parser;
mod samples;

use crate::filter::SiteView;
use crate::gt::{GenotypeSource, SampleFields};
use crate::stats::SiteStats;
use gloci_io::{UsageError, VcfHeader};
use lexer::Op;

pub use eval::MaskedGenotypes;

/// What `vcf query` does with a record under its filter (D172).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QueryPass {
    /// The record does not print.
    Drop,
    /// The record prints with every selected sample.
    All,
    /// The record prints; `pass[i]` marks the header samples that print (unselected ones are false).
    Samples,
}

/// A parsed filter expression.
#[derive(Debug)]
pub struct Expr {
    root: Node,
    needs_stats: bool,
    needs_genotypes: bool,
    needs_info: bool,
}

/// The type of an INFO field, from its header definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum InfoTy {
    Num,
    Str,
    Flag,
}

#[derive(Debug)]
pub(crate) enum Field {
    Qual,
    Pos,
    Chrom,
    Id,
    Ref,
    Alt,
    Filter,
    Type,
    NAlt,
    /// The header's sample count, used when neither stats nor genotypes are given.
    NSamples(usize),
    Info {
        key: Box<[u8]>,
        ty: InfoTy,
        idx: Option<usize>,
    },
    /// `FMT/DP`, `FMT/GQ`, `GT`: the optional index is a sample subscript
    /// (0-based, numbered over the whole header), as in bcftools.
    FmtDp(Option<usize>),
    FmtGq(Option<usize>),
    Gt(Option<usize>),
    Maf,
    Ac(Option<usize>),
    An,
    FMissing,
    NMissing,
    NPass(Box<Node>),
}

/// The right-hand side of a comparison, resolved at parse time.
pub(crate) enum Rhs {
    Num(f64),
    Str(Box<[u8]>),
    Miss,
    Re(regex::bytes::Regex),
    /// A `GT` keyword: a bit set of genotype classes.
    GtKw(u8),
    /// A `TYPE` keyword.
    TypeKw(TypeKw),
}

impl std::fmt::Debug for Rhs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Rhs::Num(n) => write!(f, "Num({n})"),
            Rhs::Str(s) => write!(f, "Str({:?})", String::from_utf8_lossy(s)),
            Rhs::Miss => write!(f, "Miss"),
            Rhs::Re(r) => write!(f, "Re({})", r.as_str()),
            Rhs::GtKw(k) => write!(f, "GtKw({k:#b})"),
            Rhs::TypeKw(k) => write!(f, "TypeKw({k:?})"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TypeKw {
    Snp,
    Mnp,
    Indel,
    Other,
    Ref,
}

#[derive(Debug)]
pub(crate) enum Node {
    Or(Vec<Node>),
    And(Vec<Node>),
    Not(Box<Node>),
    Cmp {
        field: Field,
        op: Op,
        rhs: Rhs,
    },
    /// A bare flag: true when the INFO flag is present.
    Truthy(Field),
}

/// Reusable buffers for [`Expr::eval`].
#[derive(Default)]
pub struct EvalScratch {
    fields: SampleFields,
    gt_text: Vec<u8>,
    /// Samples considered by `FMT/` fields and `N_PASS` (`None` = all).
    keep: Option<Vec<bool>>,
    /// Free per-sample masks for [`Expr::eval_query`].
    masks: Vec<Vec<bool>>,
}

impl EvalScratch {
    /// Restrict the samples `FMT/`, `GT` and `N_PASS` consider to those marked
    /// `true` (the `--samples` selection). `None` considers all. Call once per
    /// header; the selection persists across records.
    pub fn set_selection(&mut self, keep: Option<&[bool]>) {
        match (keep, &mut self.keep) {
            (None, k) => *k = None,
            (Some(src), Some(k)) => {
                k.clear();
                k.extend_from_slice(src);
            }
            (Some(src), k) => *k = Some(src.to_vec()),
        }
    }
}

/// Parse `src` against `header`, resolving field names and types.
///
/// Errors are usage errors that quote the expression, give the 1-based column
/// of the bad token and say how to fix it.
pub fn parse(src: &str, header: &VcfHeader) -> Result<Expr, UsageError> {
    parser::parse(src, header)
}

impl Expr {
    /// Does evaluation read genotypes (`FMT/`, `GT`, `N_PASS` or a computed
    /// field)? When false, pass `None` for `gt` and `stats`.
    pub fn needs_genotypes(&self) -> bool {
        self.needs_genotypes
    }

    /// Does evaluation read `stats` (`MAF`, `AC`, `AN`, `F_MISSING`,
    /// `N_MISSING`)? `N_SAMPLES` uses `stats` when given.
    pub fn needs_stats(&self) -> bool {
        self.needs_stats
    }

    /// Does evaluation read INFO values? A BCF `SiteView` carries them only
    /// when `BcfSite::info_values` is set.
    pub fn needs_info_values(&self) -> bool {
        self.needs_info
    }

    /// Evaluate on one record. `gt` must be given when
    /// [`needs_genotypes`](Expr::needs_genotypes), `stats` when
    /// [`needs_stats`](Expr::needs_stats); a missing input is an `Err`.
    pub fn eval(
        &self,
        site: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        stats: Option<&SiteStats>,
        scratch: &mut EvalScratch,
    ) -> Result<bool, String> {
        if self.needs_info && !site.info_has_values {
            return Err(
                "the expression reads INFO values but this BCF record view carries \
                 INFO keys only; build the view with BcfSite::with_info_values() \
                 (check Expr::needs_info_values()); this is a bug, please report it"
                    .to_string(),
            );
        }
        eval::eval(&self.root, site, gt, stats, scratch)
    }

    /// bcftools query's record rule (D172): `exclude` is `--exclude`. `pass` is resized to the
    /// record's sample count when the result is per sample. Same inputs as [`Expr::eval`];
    /// no allocation once `scratch` and `pass` have grown.
    pub fn eval_query(
        &self,
        exclude: bool,
        site: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        stats: Option<&SiteStats>,
        scratch: &mut EvalScratch,
        pass: &mut Vec<bool>,
    ) -> Result<QueryPass, String> {
        if self.needs_info && !site.info_has_values {
            return Err(
                "the expression reads INFO values but this BCF record view carries \
                 INFO keys only; build the view with BcfSite::with_info_values() \
                 (check Expr::needs_info_values()); this is a bug, please report it"
                    .to_string(),
            );
        }
        samples::eval_query(&self.root, exclude, site, gt, stats, scratch, pass)
    }
}
