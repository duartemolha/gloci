//! Population genetics statistics as vcftools 0.1.16 computes them, run by
//! the `vcf stats` engine ([`crate::StatsEngine`]): per-site and windowed
//! nucleotide diversity (π), Tajima's D, Weir and Cockerham's Fst, linkage
//! disequilibrium and the relatedness of pairs of samples.

mod diversity;
mod fst;
mod ld;
mod relatedness;
mod window;

pub use diversity::{SitePi, TajimaD, WindowPi};
pub use fst::{Populations, TWO_FILES as FST_TWO_FILES, WeirFst};
pub use ld::{GenoR2, HapR2, InterchromR2, LdBudget, LdWindow};
pub use relatedness::{RelatednessManichaikul, RelatednessYang};

use crate::stats::{SiteStats, no_polyploid};

/// What else `vcf popgen` allows when a genotype of ploidy above two stops it.
const POLYPLOID_HINT: &str = "or leave their samples out with --exclude-samples";

/// vcftools stops on a genotype of ploidy above two and uses a site only
/// when every other selected, unmasked sample is diploid: an error for the
/// first, false for a haploid (or lone `.`) genotype.
fn diploid(st: &SiteStats) -> Result<bool, String> {
    no_polyploid(st.polyploid, POLYPLOID_HINT)?;
    Ok(!st.haploid)
}
