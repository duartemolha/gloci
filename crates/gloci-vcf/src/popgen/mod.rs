//! Population genetics statistics as vcftools 0.1.16 computes them, run by
//! the `vcf stats` engine ([`crate::StatsEngine`]): per-site and windowed
//! nucleotide diversity (π) and Tajima's D.

mod diversity;
mod window;

pub use diversity::{SitePi, TajimaD, WindowPi};

use crate::filter::SiteView;
use crate::stats::SiteStats;
use std::collections::HashSet;

/// vcftools uses a site only when every selected, unmasked sample is diploid.
fn diploid(st: &SiteStats) -> bool {
    !st.haploid && !st.polyploid
}

/// The order a windowed statistic needs: each chromosome in one run and
/// positions that never go backwards within it, because rows are written as
/// the input passes them.
struct SortedInput {
    /// The option that needs the order, for the message.
    flag: &'static str,
    /// The current chromosome and the last position on it.
    chrom: Option<Vec<u8>>,
    pos: u64,
    /// The chromosomes already finished.
    done: HashSet<Vec<u8>>,
}

impl SortedInput {
    fn new(flag: &'static str) -> SortedInput {
        SortedInput {
            flag,
            chrom: None,
            pos: 0,
            done: HashSet::new(),
        }
    }

    /// Check `s` against the records before it. `Ok(true)` when it starts a
    /// chromosome (the first record included).
    fn check(&mut self, s: &SiteView) -> Result<bool, String> {
        if self.chrom.as_deref() == Some(s.chrom) {
            if s.pos < self.pos {
                return Err(self.unsorted(s, "a position goes backwards"));
            }
            self.pos = s.pos;
            return Ok(false);
        }
        if self.done.contains(s.chrom) {
            return Err(self.unsorted(s, "the chromosome reappears"));
        }
        if let Some(c) = self.chrom.replace(s.chrom.to_vec()) {
            self.done.insert(c);
        }
        self.pos = s.pos;
        Ok(true)
    }

    fn unsorted(&self, s: &SiteView, why: &str) -> String {
        format!(
            "{} needs input sorted by chromosome and position, but {why} at {}:{}; run: gloci \
             sort FILE (a streamed statistic cannot merge rows it already wrote)",
            self.flag,
            String::from_utf8_lossy(s.chrom),
            s.pos
        )
    }
}
