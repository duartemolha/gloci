//! `sample`: a uniform random subset of the records (D151), the equivalent
//! of `bedtools sample`.
//!
//! Complying records (with a strand filter, those on that strand; an unknown
//! strand never complies) are numbered *j* from 1 and each gets a key from
//! stream 0 of the seed, key *j* at word 2(*j*-1). The `count` records with
//! the smallest (key, *j*) are kept: a uniform subset, written in input
//! order (bedtools writes reservoir order). The result depends on the seed
//! and input only, never on threads.

use crate::rng::Stream;
use gloci_core::{ChromTable, Strand};
use gloci_io::{IoError, RecordBatch, RecordReader, RecordWriter, next_batch_size};
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Copy)]
pub struct SampleOptions {
    /// How many records to keep.
    pub count: u64,
    pub seed: u64,
    /// Sample only the records on this strand.
    pub strand: Option<Strand>,
}

/// A kept record, ordered by (key, number) so the heap's top is the one to
/// drop first.
struct Kept {
    key: u64,
    j: u64,
    raw: Vec<u8>,
}

impl PartialEq for Kept {
    fn eq(&self, o: &Kept) -> bool {
        (self.key, self.j) == (o.key, o.j)
    }
}
impl Eq for Kept {}
impl PartialOrd for Kept {
    fn partial_cmp(&self, o: &Kept) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Kept {
    fn cmp(&self, o: &Kept) -> std::cmp::Ordering {
        (self.key, self.j).cmp(&(o.key, o.j))
    }
}

/// Write a uniform random subset of `o.count` records of `reader`, in input
/// order, and return how many were written. Fewer complying records than
/// `o.count` is an error and nothing is written. Memory is the kept records.
pub fn sample(
    mut reader: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &SampleOptions,
    w: &mut dyn RecordWriter,
) -> Result<u64, IoError> {
    let want = usize::try_from(o.count).unwrap_or(usize::MAX);
    let mut stream = Stream::new(o.seed, 0);
    let mut heap: BinaryHeap<Kept> = BinaryHeap::new();
    let mut seen = 0u64;
    let mut read = 0u64;
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = reader.read_batch(chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            read += 1;
            if o.strand.is_some_and(|s| m.strand != s) {
                continue;
            }
            seen += 1;
            let key = stream.next_u64();
            // The second word of the pair keeps key j at word 2(j-1).
            stream.next_u64();
            if want == 0 {
                continue;
            }
            if heap.len() < want {
                heap.push(Kept {
                    key,
                    j: seen,
                    raw: batch.raw(i).to_vec(),
                });
            } else if heap
                .peek()
                .is_some_and(|top| (key, seen) < (top.key, top.j))
            {
                heap.pop();
                heap.push(Kept {
                    key,
                    j: seen,
                    raw: batch.raw(i).to_vec(),
                });
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    crate::skipped::log_reader(reader.as_ref());
    if read == 0 && o.count > 0 {
        // An empty input is read as no records, as every other command does.
        log::warn!(
            "{} is empty: it has no records; no record was sampled",
            reader.name()
        );
        return Ok(0);
    }
    if seen < o.count {
        let what = match o.strand {
            Some(Strand::Plus) => " on the + strand",
            Some(Strand::Minus) => " on the - strand",
            _ => "",
        };
        return Err(IoError::Detect {
            path: reader.name().to_string(),
            msg: format!(
                "{} has only {seen} records{what}, fewer than the {} requested with -n; \
ask for at most {seen} (nothing was written)",
                reader.name(),
                o.count
            ),
        });
    }
    let mut kept = heap.into_vec();
    kept.sort_unstable_by_key(|k| k.j);
    for k in &kept {
        w.write_record(&k.raw).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })?;
    }
    Ok(kept.len() as u64)
}
