//! The `--sort` entry points of gloci-ops over in-memory readers: a relation
//! run with `RelateOptions::sort`, and `sort_together` for merge.

use gloci_core::{ChromOrder, ChromTable, OverlapCriteria, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, RecordWriter};
use gloci_ops::intersect::{IntersectOptions, Match, Mode, intersect};
use gloci_ops::merge::{MergeOptions, merge};
use gloci_ops::{RelateInputs, RelateOptions, SortSettings, sort_input_options, sort_together};
use gloci_testkit::mem::mem_reader;
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

fn user_opts() -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn settings() -> SortSettings {
    SortSettings {
        max_mem: 16 << 20,
        tmp_dir: std::env::temp_dir(),
        threads: NonZeroUsize::MIN,
        vcf_span: VcfSpan::Allele,
    }
}

const A_SHUFFLED: &str = "chr2\t10\t50\tx\nchr1\t500\t600\ty\nchr1\t10\t20\tz\nchr10\t0\t9\tw\n";
const A_SORTED: &str = "chr1\t10\t20\tz\nchr1\t500\t600\ty\nchr2\t10\t50\tx\nchr10\t0\t9\tw\n";
const B1_SHUFFLED: &str = "chr10\t5\t6\nchr1\t550\t560\nchr1\t15\t16\nchr1\t0\t1000\n";
const B1_SORTED: &str = "chr1\t0\t1000\nchr1\t15\t16\nchr1\t550\t560\nchr10\t5\t6\n";
const MERGED_SORTED: &str = "chr1\t0\t1000\nchr1\t10\t20\nchr1\t15\t16\nchr1\t500\t600\nchr1\t550\t560\nchr2\t10\t50\nchr10\t0\t9\nchr10\t5\t6\n";
const B2_SHUFFLED: &str = "chr2\t0\t100\nchr1\t599\t700\n";
const B2_SORTED: &str = "chr1\t599\t700\nchr2\t0\t100\n";

fn options(sort: Option<SortSettings>) -> IntersectOptions {
    IntersectOptions {
        mode: Mode::Pairs {
            overlap_bases: true,
            keep_unmatched: true,
        },
        matching: Match::Overlap,
        relate: RelateOptions {
            criteria: OverlapCriteria::default(),
            window: None,
            split: false,
            unsorted: false,
            sort,
            labels: vec!["B1".into(), "B2".into()],
            genome: None,
            chrom_order: ChromOrder::Karyotype,
            header: false,
        },
    }
}

fn run(a: &str, b: [&str; 2], ro: &ReaderOptions, o: &IntersectOptions) -> Vec<String> {
    let inputs = RelateInputs {
        a: mem_reader("a.bed", a, ro),
        b: vec![
            mem_reader("b1.bed", b[0], ro),
            mem_reader("b2.bed", b[1], ro),
        ],
    };
    let mut sink = Sink::default();
    intersect(inputs, o, &mut sink).unwrap();
    sink.0
}

#[test]
fn relation_with_sort_settings_equals_presorted_inputs() {
    let want = run(
        A_SORTED,
        [B1_SORTED, B2_SORTED],
        &user_opts(),
        &options(None),
    );
    assert!(want.len() > 4, "{want:?}");
    let o = options(Some(settings()));
    let got = run(
        A_SHUFFLED,
        [B1_SHUFFLED, B2_SHUFFLED],
        &sort_input_options(&user_opts()),
        &o,
    );
    assert_eq!(got, want);
}

#[test]
fn relation_sort_and_unsorted_together_is_a_usage_error() {
    let mut o = options(Some(settings()));
    o.relate.unsorted = true;
    let ro = sort_input_options(&user_opts());
    let inputs = RelateInputs {
        a: mem_reader("a.bed", A_SHUFFLED, &ro),
        b: vec![mem_reader("b1.bed", B1_SHUFFLED, &ro)],
    };
    let e = intersect(inputs, &o, &mut Sink::default()).unwrap_err();
    assert!(e.to_string().contains("cannot be combined"), "{e}");
}

#[test]
fn merge_with_sort_together_equals_presorted_input() {
    let o = MergeOptions {
        max_gap: 0,
        per_strand: false,
        count: true,
        aggs: vec![],
        empty: ".".into(),
        delim: ",".into(),
    };
    let ro = user_opts();
    let mut chroms = ChromTable::karyotype();
    let mut want = Sink::default();
    merge(
        vec![mem_reader("s.bed", MERGED_SORTED, &ro)],
        &mut chroms,
        &o,
        &mut want,
    )
    .unwrap();
    // Two shuffled inputs sorted together; the order is the sorted order.
    let ro = sort_input_options(&user_opts());
    let mut chroms = ChromTable::karyotype();
    let sorted = sort_together(
        vec![
            mem_reader("a.bed", A_SHUFFLED, &ro),
            mem_reader("b.bed", B1_SHUFFLED, &ro),
        ],
        &mut chroms,
        &settings(),
    )
    .unwrap();
    assert_eq!(sorted.len(), 1);
    let mut got = Sink::default();
    merge(sorted, &mut chroms, &o, &mut got).unwrap();
    assert!(!want.0.is_empty());
    assert_eq!(got.0, want.0);
}
