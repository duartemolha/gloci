//! `merge` against the brute-force `naive_merge`, including negative gaps,
//! `-s` and several inputs.

use gloci_core::{ChromOrder, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter, chrom_table_for};
use gloci_ops::agg::parse_aggs;
use gloci_ops::merge::{MergeOptions, merge};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, naive_merge};
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

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn merge_matches_naive(
        input in relation_inputs(),
        gap in prop::sample::select(vec![-20i64, -2, -1, 0, 1, 2, 10, 1000]),
        per_strand in any::<bool>(),
    ) {
        let rank = naive::karyotype_rank(CHROMS);
        let mut files = vec![naive::sort(input.a.clone(), &rank)];
        files.extend(input.b.iter().map(|f| naive::sort(f.clone(), &rank)));
        let ro = ReaderOptions { accept: gloci_io::Accept::default(), chrom_map: None,
            threads: NonZeroUsize::MIN,
            format: Some(Format::Bed),
            vcf_span: VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
        };
        let readers: Vec<_> = files
            .iter()
            .enumerate()
            .map(|(i, f)| mem_reader(&format!("f{i}.bed"), &text(f), &ro))
            .collect();
        let refs: Vec<&dyn gloci_io::RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
        let mut chroms = chrom_table_for(None, &refs, ChromOrder::Karyotype).unwrap();
        let o = MergeOptions {
            max_gap: gap,
            per_strand,
            count: true,
            aggs: parse_aggs(&["4:collapse".to_string()]).unwrap(),
            empty: ".".into(),
            delim: ",".into(),
        };
        let mut sink = Sink::default();
        merge(readers, &mut chroms, &o, &mut sink).unwrap();
        prop_assert_eq!(sink.0, naive_merge(&files, &rank, gap, per_strand));
    }
}

/// With `-s`, one early `-` and one early `.` record must not hold back the
/// finished `+` clusters that follow: the reorder buffer stays small.
#[test]
fn per_strand_buffer_drains() {
    use gloci_ops::merge::merge_counting;
    let mut s = String::from("chr1\t0\t10\tm\t0\t-\nchr1\t5\t12\tu\t0\t.\n");
    for i in 0..5000u64 {
        s += &format!("chr1\t{}\t{}\tp{i}\t0\t+\n", 100 + i * 20, 110 + i * 20);
    }
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let r = mem_reader("f.bed", &s, &ro);
    let mut chroms = chrom_table_for(None, &[r.as_ref()], ChromOrder::Karyotype).unwrap();
    let o = MergeOptions {
        max_gap: 0,
        per_strand: true,
        count: false,
        aggs: vec![],
        empty: ".".into(),
        delim: ",".into(),
    };
    let mut sink = Sink::default();
    let peak = merge_counting(vec![r], &mut chroms, &o, &mut sink).unwrap();
    assert_eq!(sink.0.len(), 5002);
    assert!(peak <= 4, "reorder buffer peaked at {peak}");
}
