//! `genomecov` against the brute-force `naive_genomecov` (a depth array per
//! chromosome), in every mode with strand filters and 5'/3' ends; and the
//! memory bound of the sweep.

use gloci_core::{ChromTable, Genome, OrderSource, Strand, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::genomecov::{
    GenomecovMode, GenomecovOptions, Piece, format_g6, genomecov, genomecov_counting,
};
use gloci_testkit::TestRec;
use gloci_testkit::mem::mem_reader;
use gloci_testkit::naive::{self, NaiveCovMode, NaiveGenomecov, naive_genomecov};
use gloci_testkit::strategies::{CHROMS, relation_inputs};
use proptest::prelude::*;
use std::io;
use std::num::NonZeroUsize;

const LEN: u64 = 200;

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

fn genome() -> Genome {
    Genome::new(
        CHROMS
            .iter()
            .map(|n| (n.as_bytes().to_vec(), Some(LEN)))
            .collect(),
    )
    .unwrap()
}

fn text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| format!("{}\n", r.bed_line())).collect()
}

/// Squeeze the generated records into the small test genome.
fn shrink(recs: &[TestRec]) -> Vec<TestRec> {
    recs.iter()
        .map(|r| {
            let start = r.start % 150;
            let end = (start + (r.end - r.start) % 60).min(LEN);
            TestRec {
                start,
                end,
                ..r.clone()
            }
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn genomecov_matches_naive(
        input in relation_inputs(),
        mode in 0usize..4,
        strand in prop::sample::select(vec![None, Some('+'), Some('-')]),
        end in prop::sample::select(vec![0u8, 5, 3]),
        max in prop::sample::select(vec![None, Some(0u64), Some(1), Some(2)]),
    ) {
        let rank = naive::karyotype_rank(CHROMS);
        let recs = naive::sort(shrink(&input.a), &rank);
        let g = genome();
        let named: Vec<(&str, u64)> = CHROMS.iter().map(|n| (*n, LEN)).collect();
        let (m, nm) = [
            (GenomecovMode::Hist, NaiveCovMode::Hist),
            (GenomecovMode::BedGraph, NaiveCovMode::BedGraph),
            (GenomecovMode::BedGraphAll, NaiveCovMode::BedGraphAll),
            (GenomecovMode::PerBase, NaiveCovMode::PerBase),
        ][mode];
        let max = if m == GenomecovMode::Hist { max } else { None };
        let o = GenomecovOptions {
            mode: m,
            scale: None,
            strand: strand.map(|c| if c == '+' { Strand::Plus } else { Strand::Minus }),
            split: false,
            split_deletions: false,
            piece: match end { 5 => Piece::FivePrime, 3 => Piece::ThreePrime, _ => Piece::Whole },
            max_depth: max,
        };
        let q = NaiveGenomecov { mode: nm, strand, end, max_depth: max };
        let reader = mem_reader("in.bed", &text(&recs), &ro());
        let mut chroms = ChromTable::closed(&g, OrderSource::Karyotype);
        let mut sink = Sink::default();
        genomecov(reader, &mut chroms, &g, &o, &mut sink).unwrap();
        let want = naive_genomecov(&recs, &named, &q, |b, l| format_g6(b as f64 / l as f64));
        prop_assert_eq!(sink.0, want);
    }
}

/// Five thousand short features: a handful of events are pending at once,
/// whatever the file size, in every output mode.
#[test]
fn pending_events_stay_bounded() {
    let mut a = String::new();
    for i in 0..5000u64 {
        a += &format!("chr1\t{}\t{}\n", i * 20, i * 20 + 12);
        a += &format!("chr1\t{}\t{}\n", i * 20 + 6, i * 20 + 18);
    }
    let g = Genome::new(vec![(b"chr1".to_vec(), Some(200_000))]).unwrap();
    for mode in [
        GenomecovMode::Hist,
        GenomecovMode::BedGraph,
        GenomecovMode::PerBase,
    ] {
        let o = GenomecovOptions {
            mode,
            ..GenomecovOptions::default()
        };
        let mut chroms = ChromTable::closed(&g, OrderSource::Karyotype);
        let mut sink = Sink::default();
        let peak = genomecov_counting(
            mem_reader("in.bed", &a, &ro()),
            &mut chroms,
            &g,
            &o,
            &mut sink,
        )
        .unwrap();
        assert!(peak <= 6, "{mode:?}: {peak} events pending at once");
    }
}
