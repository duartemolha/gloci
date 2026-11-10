//! Binary search on uncompressed sorted BED (D11), and region views read
//! that way (`assume_sorted`) or by the default checked scan, against the
//! naive oracle.

use gloci_core::{ChromOrder, ChromTable, Genome, Interval, OrderSource};
use gloci_io::{InvalidPolicy, IoError, ReaderOptions, RecordBatch};
use gloci_ops::{ChromBlock, ViewOptions, find_chrom_block, open_view};
use gloci_testkit::naive;
use gloci_testkit::strategies::{ABSENT_CHROM, CHROMS, TestRec, bed_records, regions};
use proptest::prelude::*;
use std::io::Cursor;
use std::num::NonZeroUsize;

fn find(text: &str, chrom: &str, table: &ChromTable) -> ChromBlock {
    let mut c = Cursor::new(text.as_bytes());
    find_chrom_block(&mut c, text.len() as u64, chrom.as_bytes(), table).unwrap()
}

fn reader_options() -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

/// Every record line `open_view` yields, with a newline each.
fn view(path: &str, opts: &ViewOptions) -> String {
    try_view(path, opts).unwrap()
}

fn try_view(path: &str, opts: &ViewOptions) -> Result<String, IoError> {
    let (mut r, mut chroms) = open_view(path, opts, &reader_options())?;
    let mut out = String::new();
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 7)?;
        for i in 0..batch.len() {
            out.push_str(std::str::from_utf8(batch.raw(i)).unwrap());
            out.push('\n');
        }
        if !more {
            break;
        }
    }
    Ok(out)
}

/// Region strings (1-based, closed) for regions of the strategy.
fn region_strings(regs: &[(Vec<u8>, Interval)]) -> Vec<String> {
    regs.iter()
        .map(|(c, iv)| {
            let c = String::from_utf8_lossy(c);
            if iv.end == u64::MAX {
                format!("{c}:{}-", iv.start + 1)
            } else {
                format!("{c}:{}-{}", iv.start + 1, iv.end)
            }
        })
        .collect()
}

/// The records that overlap a region, once each, in file order.
fn naive_view(recs: &[TestRec], regs: &[(Vec<u8>, Interval)], inclusive: bool) -> String {
    recs.iter()
        .filter(|r| {
            regs.iter()
                .any(|(c, iv)| naive::matches(r, c, *iv, inclusive))
        })
        .map(|r| r.bed_line() + "\n")
        .collect()
}

/// Are each chromosome's records contiguous and sorted by start?
fn contiguous_sorted(recs: &[TestRec]) -> bool {
    let mut done = std::collections::HashSet::new();
    recs.windows(2).all(|w| {
        if w[0].chrom == w[1].chrom {
            w[1].start >= w[0].start
        } else {
            done.insert(w[0].chrom.clone());
            !done.contains(&w[1].chrom)
        }
    })
}

fn table(lex: bool) -> ChromTable {
    if lex {
        ChromTable::lexicographic()
    } else {
        ChromTable::karyotype()
    }
}

fn rank(lex: bool) -> Box<dyn Fn(&str) -> u32> {
    if lex {
        let mut names: Vec<&str> = CHROMS.to_vec();
        names.sort();
        Box::new(move |c: &str| names.iter().position(|n| *n == c).unwrap() as u32)
    } else {
        Box::new(naive::karyotype_rank(CHROMS))
    }
}

/// BED text of `recs`, optionally with leading header lines and with blank
/// lines and CRLF endings sprinkled in, and the byte offset of the first
/// line of each chromosome.
fn bed_text(recs: &[TestRec], noise: bool) -> (String, Vec<(String, u64)>) {
    let mut text = String::new();
    if noise {
        text.push_str("#comment\ntrack name=x\nbrowser position chr1:1-10\n");
    }
    let mut firsts: Vec<(String, u64)> = Vec::new();
    for (i, r) in recs.iter().enumerate() {
        if noise && i % 5 == 3 {
            text.push('\n');
        }
        if firsts.last().is_none_or(|(c, _)| *c != r.chrom) {
            firsts.push((r.chrom.clone(), text.len() as u64));
        }
        text.push_str(&r.bed_line());
        text.push_str(if noise && i % 3 == 1 { "\r\n" } else { "\n" });
    }
    (text, firsts)
}

/// The offset of the line that starts with `prefix`.
fn at(text: &str, prefix: &str) -> ChromBlock {
    let off = text
        .match_indices(prefix)
        .find(|&(i, _)| i == 0 || text.as_bytes()[i - 1] == b'\n')
        .expect("prefix starts a line")
        .0;
    ChromBlock::At(off as u64)
}

#[test]
fn finds_first_line_of_each_chromosome() {
    let text = "chr1\t1\t2\nchr1\t5\t9\nchr2\t1\t2\nchr10\t0\t5\nchr10\t9\t9\nchrX\t1\t2\n";
    let t = ChromTable::karyotype();
    assert_eq!(find(text, "chr1", &t), ChromBlock::At(0));
    assert_eq!(find(text, "chr2", &t), at(text, "chr2\t1"));
    assert_eq!(find(text, "chr10", &t), at(text, "chr10\t0"));
    assert_eq!(find(text, "chrX", &t), at(text, "chrX"));
    for absent in ["chr3", "chr0", "chrY", "1", "chr1_random"] {
        assert_eq!(find(text, absent, &t), ChromBlock::NotFoundByBisection);
    }
    assert_eq!(find("", "chr1", &t), ChromBlock::NotFoundByBisection);
    assert_eq!(
        find("track name=x\n\n#c\n", "chr1", &t),
        ChromBlock::NotFoundByBisection
    );
}

#[test]
fn skips_header_and_blank_lines() {
    let text = "track name=x\n#c\n\nchr1\t1\t2\n\n#late\nchr2\t1\t2\r\n\r\nchr3\t5\t6\n";
    let t = ChromTable::karyotype();
    assert_eq!(find(text, "chr1", &t), at(text, "chr1"));
    assert_eq!(find(text, "chr2", &t), at(text, "chr2"));
    assert_eq!(find(text, "chr3", &t), at(text, "chr3"));
}

#[test]
fn long_lines() {
    let pad = "x".repeat(100_000);
    let text = format!("chr1\t1\t2\t{pad}\nchr2\t1\t2\t{pad}\nchr3\t1\t2\n");
    let t = ChromTable::karyotype();
    assert_eq!(find(&text, "chr2", &t), ChromBlock::At(100_010));
    assert_eq!(find(&text, "chr3", &t), ChromBlock::At(200_020));
}

#[test]
fn table_order_decides() {
    // lexicographic file, lexicographic table
    let text = "chr1\t1\t2\nchr10\t1\t2\nchr2\t1\t2\n";
    let lex = ChromTable::lexicographic();
    assert_eq!(find(text, "chr10", &lex), at(text, "chr10"));
    assert_eq!(find(text, "chr2", &lex), at(text, "chr2"));
    // a genome table: names it cannot place stop the search
    let g = Genome::new(vec![(b"chr2".to_vec(), None), (b"chr1".to_vec(), None)]).unwrap();
    let closed = ChromTable::closed(&g, OrderSource::GenomeFile("g".into()));
    let text = "chr2\t1\t2\nchr1\t1\t2\n";
    assert_eq!(find(text, "chr1", &closed), at(text, "chr1"));
    assert_eq!(find(text, "chr2", &closed), ChromBlock::At(0));
    let text = "chr2\t1\t2\nchrUn\t1\t2\nchr1\t1\t2\n";
    assert_eq!(find(text, "chr1", &closed), ChromBlock::NotFoundByBisection);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    /// In a file sorted in the table's order, every chromosome is found at
    /// its first line and absent ones are not found.
    #[test]
    fn bisection_finds_each_chromosome_block(
        recs in bed_records(),
        lex in any::<bool>(),
        noise in any::<bool>(),
    ) {
        let recs = naive::sort(recs, rank(lex));
        let (text, firsts) = bed_text(&recs, noise);
        let t = table(lex);
        for c in CHROMS.iter().chain([&ABSENT_CHROM]) {
            let want = firsts
                .iter()
                .find(|(n, _)| n == c)
                .map_or(ChromBlock::NotFoundByBisection, |&(_, at)| ChromBlock::At(at));
            prop_assert_eq!(find(&text, c, &t), want, "chromosome {}", c);
        }
    }

    /// A region view of a plain sorted BED file with `assume_sorted`
    /// (binary search, or the verified full scan when the search misses)
    /// equals the naive query: every record overlapping a region, once, in
    /// file order.
    #[test]
    fn view_by_bisection_matches_naive(
        recs in bed_records(),
        regs in regions(4),
        file_lex in any::<bool>(),
        query_lex in any::<bool>(),
        noise in any::<bool>(),
        inclusive in any::<bool>(),
    ) {
        let recs = naive::sort(recs, rank(file_lex));
        let (text, _) = bed_text(&recs, noise);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.bed");
        std::fs::write(&path, &text).unwrap();
        let opts = ViewOptions {
            regions: region_strings(&regs),
            inclusive_points: inclusive,
            chrom_order: if query_lex { ChromOrder::Lex } else { ChromOrder::Karyotype },
            assume_sorted: true,
            ..ViewOptions::default()
        };
        prop_assert_eq!(view(&path.to_string_lossy(), &opts), naive_view(&recs, &regs, inclusive));
    }

    /// Without `assume_sorted`, a plain BED file is scanned with its order
    /// checked: the answer is the naive one exactly when each chromosome's
    /// records are contiguous and sorted, and an order error otherwise, so
    /// never a wrong answer.
    #[test]
    fn default_scan_is_right_or_refuses(
        recs in bed_records(),
        regs in regions(4),
        sorted in any::<bool>(),
        inclusive in any::<bool>(),
    ) {
        let recs = if sorted { naive::sort(recs, rank(false)) } else { recs };
        let (text, _) = bed_text(&recs, false);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("a.bed");
        std::fs::write(&path, &text).unwrap();
        let opts = ViewOptions {
            regions: region_strings(&regs),
            inclusive_points: inclusive,
            ..ViewOptions::default()
        };
        match try_view(&path.to_string_lossy(), &opts) {
            Ok(out) => {
                prop_assert!(contiguous_sorted(&recs), "unsorted input accepted");
                prop_assert_eq!(out, naive_view(&recs, &regs, inclusive));
            }
            Err(e) => {
                prop_assert!(matches!(e, IoError::Order(_)), "{}", e);
                prop_assert!(!contiguous_sorted(&recs), "sorted input refused: {}", e);
            }
        }
    }
}
