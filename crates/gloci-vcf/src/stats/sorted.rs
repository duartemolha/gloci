//! The input order a streamed statistic needs.

use crate::filter::SiteView;
use std::collections::HashSet;

/// The order a statistic that writes rows as the input passes them needs:
/// each chromosome in one run and, within it, an order key (the position,
/// or the bin of `--snp-density`) that never goes backwards. A chromosome
/// that reappears or a key that goes backwards is an error naming the flag,
/// the record and `gloci sort`.
pub(crate) struct SortedInput {
    /// The option that needs the order, for the message.
    flag: &'static str,
    /// What streams its rows (`statistic`, `density`), for the message.
    noun: &'static str,
    /// The current chromosome and the last key on it.
    chrom: Option<Vec<u8>>,
    key: u64,
    /// The chromosomes already finished.
    done: HashSet<Vec<u8>>,
}

impl SortedInput {
    pub(crate) fn new(flag: &'static str, noun: &'static str) -> SortedInput {
        SortedInput {
            flag,
            noun,
            chrom: None,
            key: 0,
            done: HashSet::new(),
        }
    }

    /// Check `s` against the records before it, with its position as the
    /// key. `Ok(true)` when it starts a chromosome (the first record
    /// included).
    pub(crate) fn check(&mut self, s: &SiteView) -> Result<bool, String> {
        self.check_key(s, s.pos)
    }

    /// [`SortedInput::check`] with `key` in place of the position.
    pub(crate) fn check_key(&mut self, s: &SiteView, key: u64) -> Result<bool, String> {
        if self.chrom.as_deref() == Some(s.chrom) {
            if key < self.key {
                return Err(self.unsorted(s, "a position goes backwards"));
            }
            self.key = key;
            return Ok(false);
        }
        if self.done.contains(s.chrom) {
            return Err(self.unsorted(s, "the chromosome reappears"));
        }
        if let Some(c) = self.chrom.replace(s.chrom.to_vec()) {
            self.done.insert(c);
        }
        self.key = key;
        Ok(true)
    }

    fn unsorted(&self, s: &SiteView, why: &str) -> String {
        format!(
            "{} needs input sorted by chromosome and position, but {why} at {}:{}; run: gloci \
             sort FILE (a streamed {} cannot merge rows it already wrote)",
            self.flag,
            String::from_utf8_lossy(s.chrom),
            s.pos,
            self.noun
        )
    }
}
