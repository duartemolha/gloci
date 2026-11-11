//! `partition` against the brute-force `naive_partition`, over several
//! inputs with nested features, points and chromosomes some inputs lack; and
//! the memory bound of the sweep.

use gloci_core::{ChromOrder, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter, chrom_table_for};
use gloci_ops::partition::{PartitionOptions, partition, partition_counting};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, NaivePartition, naive_partition};
use gloci_testkit::strategies::{CHROMS, relation_inputs};
use proptest::prelude::*;
use std::io;
use std::num::NonZeroUsize;

#[derive(Default)]
struct Sink(Vec<String>);

impl RecordWriter for Sink {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        self.0.push(String::from_utf8(raw.to_vec()).unwrap());
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), gloci_io::IoError> {
        Ok(())
    }
}

fn text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| format!("{}\n", r.bed_line())).collect()
}

fn ro() -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn partition_matches_naive(
        input in relation_inputs(),
        min in prop::sample::select(vec![None, Some(1usize), Some(2), Some(3)]),
        max in prop::sample::select(vec![None, Some(1usize), Some(2)]),
        merge_adjacent in any::<bool>(),
        membership in any::<bool>(),
        all in any::<bool>(),
    ) {
        let rank = naive::karyotype_rank(CHROMS);
        let mut files = vec![naive::sort(input.a.clone(), &rank)];
        files.extend(input.b.iter().map(|f| naive::sort(f.clone(), &rank)));
        let labels: Vec<String> = (0..files.len()).map(|i| format!("f{i}.bed")).collect();
        let readers: Vec<_> = files
            .iter()
            .zip(&labels)
            .map(|(f, l)| mem_reader(l, &text(f), &ro()))
            .collect();
        let refs: Vec<&dyn gloci_io::RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
        let mut chroms = chrom_table_for(None, &refs, ChromOrder::Karyotype).unwrap();
        prop_assume!(all || min.zip(max).is_none_or(|(a, b)| a <= b));
        let o = PartitionOptions {
            names: labels.clone(),
            min_files: if all { None } else { min },
            max_files: if all { None } else { max },
            all,
            merge_adjacent,
            membership,
            bedgraph: false,
            empty: None,
        };
        let q = NaivePartition {
            labels,
            min_files: if all { files.len() } else { min.unwrap_or(1) },
            max_files: if all { files.len() } else { max.unwrap_or(usize::MAX) },
            merge_adjacent,
            membership,
        };
        let mut sink = Sink::default();
        partition(readers, &mut chroms, None, &o, &mut sink).unwrap();
        prop_assert_eq!(sink.0, naive_partition(&files, &q, &rank));
    }
}

/// Five thousand short intervals in each of two files: at most two are open
/// at once, whatever the file size.
#[test]
fn open_intervals_stay_bounded() {
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..5000u64 {
        a += &format!("chr1\t{}\t{}\n", i * 20, i * 20 + 12);
        b += &format!("chr1\t{}\t{}\n", i * 20 + 6, i * 20 + 18);
    }
    let readers = vec![
        mem_reader("a.bed", &a, &ro()),
        mem_reader("b.bed", &b, &ro()),
    ];
    let refs: Vec<&dyn gloci_io::RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let mut chroms = chrom_table_for(None, &refs, ChromOrder::Karyotype).unwrap();
    let o = PartitionOptions {
        names: vec!["a.bed".into(), "b.bed".into()],
        membership: true,
        ..PartitionOptions::default()
    };
    let mut sink = Sink::default();
    let peak = partition_counting(readers, &mut chroms, None, &o, &mut sink).unwrap();
    assert_eq!(sink.0.len(), 5000 * 3);
    assert!(peak <= 2, "{peak} intervals open at once");
}
