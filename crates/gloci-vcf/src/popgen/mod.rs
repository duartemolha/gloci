//! Population genetics statistics as vcftools 0.1.16 computes them, run by
//! the `vcf stats` engine ([`crate::StatsEngine`]): per-site and windowed
//! nucleotide diversity (π), Tajima's D, and Weir and Cockerham's Fst.

mod diversity;
mod fst;
mod window;

pub use diversity::{SitePi, TajimaD, WindowPi};
pub use fst::{Populations, TWO_FILES as FST_TWO_FILES, WeirFst};

use crate::stats::SiteStats;

/// vcftools uses a site only when every selected, unmasked sample is diploid.
fn diploid(st: &SiteStats) -> bool {
    !st.haploid && !st.polyploid
}
