//! `complement` against the brute-force `naive_complement`, with and
//! without a genome, over several inputs.

use gloci_core::OrderSource;
use gloci_core::{ChromOrder, Genome, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter, chrom_table_for};
use gloci_ops::complement::complement;
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, naive_complement};
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
    fn complement_without_genome_matches_naive(input in relation_inputs()) {
        let rank = naive::karyotype_rank(CHROMS);
        let mut files = vec![naive::sort(input.a.clone(), &rank)];
        files.extend(input.b.iter().map(|f| naive::sort(f.clone(), &rank)));
        let readers: Vec<_> = files
            .iter()
            .enumerate()
            .map(|(i, f)| mem_reader(&format!("f{i}.bed"), &text(f), &ro()))
            .collect();
        let refs: Vec<&dyn gloci_io::RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
        let mut chroms = chrom_table_for(None, &refs, ChromOrder::Karyotype).unwrap();
        let mut sink = Sink::default();
        complement(readers, &mut chroms, None, &mut sink).unwrap();
        prop_assert_eq!(sink.0, naive_complement(&files, None, &rank));
    }

    #[test]
    fn complement_with_genome_matches_naive(input in relation_inputs()) {
        // A genome order unlike karyotype order, long enough for every coordinate.
        let order = ["chrM", "chrX", "chr10", "chr2", "chr1"];
        let sizes: Vec<(String, u64)> = order.iter().map(|n| (n.to_string(), 1u64 << 42)).collect();
        let genome = Genome::new(sizes.iter().map(|(n, l)| (n.as_bytes().to_vec(), Some(*l))).collect()).unwrap();
        let rank = |c: &str| order.iter().position(|n| *n == c).unwrap() as u32;
        let mut files = vec![naive::sort(input.a.clone(), rank)];
        files.extend(input.b.iter().map(|f| naive::sort(f.clone(), rank)));
        let readers: Vec<_> = files
            .iter()
            .enumerate()
            .map(|(i, f)| mem_reader(&format!("f{i}.bed"), &text(f), &ro()))
            .collect();
        let mut chroms = gloci_core::ChromTable::closed(&genome, OrderSource::Karyotype);
        let mut sink = Sink::default();
        complement(readers, &mut chroms, Some(&genome), &mut sink).unwrap();
        prop_assert_eq!(sink.0, naive_complement(&files, Some(&sizes), rank));
    }
}
