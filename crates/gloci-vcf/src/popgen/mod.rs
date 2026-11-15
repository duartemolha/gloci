//! Population genetics statistics as vcftools 0.1.16 computes them, run by
//! the `vcf stats` engine ([`crate::StatsEngine`]): per-site and windowed
//! nucleotide diversity (π) and Tajima's D.

mod diversity;
mod window;

pub use diversity::{SitePi, TajimaD, WindowPi};

use crate::stats::SiteStats;

/// vcftools uses a site only when every selected, unmasked sample is diploid.
fn diploid(st: &SiteStats) -> bool {
    !st.haploid && !st.polyploid
}
