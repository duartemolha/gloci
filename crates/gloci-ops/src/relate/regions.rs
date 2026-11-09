//! Restricting a plain file's records to regions (D75) when it is read per
//! chromosome (D76).

use crate::bedsearch::append_matching;
use crate::view::Matcher;
use gloci_core::{ChromTable, Interval, RegionSet};
use gloci_io::{Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordReader};
use std::collections::BTreeMap;

/// The records of `inner` that overlap a region, as `view` selects them from
/// a scan: the same records with the same line numbers.
pub(super) struct RegionFilter {
    inner: Box<dyn RecordReader>,
    matcher: Matcher,
    inclusive_points: bool,
    scratch: RecordBatch,
}

impl RegionFilter {
    /// `regions` are named as the records are (after `--chrom-map`).
    pub(super) fn new(
        inner: Box<dyn RecordReader>,
        regions: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
    ) -> RegionFilter {
        let mut by_chrom: BTreeMap<Vec<u8>, Vec<Interval>> = BTreeMap::new();
        for (c, iv) in regions {
            by_chrom.entry(c.clone()).or_default().push(*iv);
        }
        let sets = by_chrom
            .into_iter()
            .map(|(c, ivs)| (c, RegionSet::new(ivs)))
            .collect();
        RegionFilter {
            inner,
            matcher: Matcher::new(sets),
            inclusive_points,
            scratch: RecordBatch::default(),
        }
    }
}

impl RecordReader for RegionFilter {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn compression(&self) -> Compression {
        self.inner.compression()
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        let before = batch.recs.len();
        loop {
            self.scratch.clear();
            let more = self.inner.read_batch(chroms, &mut self.scratch, max)?;
            let (matcher, inclusive) = (&mut self.matcher, self.inclusive_points);
            append_matching(batch, &self.scratch, |m| {
                matcher.matches(chroms, m, inclusive)
            });
            if !more || batch.recs.len() > before {
                return Ok(more);
            }
        }
    }
}
