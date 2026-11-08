//! Tabix-style indexing of any BGZF TSV (`write_tsv_index`), checked by
//! querying the written index directly with noodles.

use gloci_io::{IndexKind, TsvColumns, find_index, write_tsv_index};
use gloci_testkit::fixtures;
use noodles_bgzf::io::Reader as BgzfReader;
use noodles_bgzf::io::Seek;
use noodles_core::Position;
use noodles_csi::BinningIndex;
use std::io::BufRead;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

fn nz1() -> NonZeroUsize {
    NonZeroUsize::new(1).unwrap()
}

fn tmp(name: &str) -> PathBuf {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    static N: AtomicU64 = AtomicU64::new(0);
    let dir = DIR.get_or_init(|| {
        let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join("gloci-io-index-tsv");
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    });
    dir.join(format!("{}-{name}", N.fetch_add(1, Ordering::Relaxed)))
}

fn bgzf(text: &str) -> PathBuf {
    let p = tmp("x.tsv.gz");
    std::fs::write(&p, fixtures::bgzip(text.as_bytes())).unwrap();
    p
}

fn cols(seq: usize, begin: usize, end: Option<usize>, zero_based: bool) -> TsvColumns {
    TsvColumns {
        seq,
        begin,
        end,
        zero_based,
        skip_lines: 0,
        comment: b'#',
    }
}

/// Query `chrom:start-end` (1-based, closed) through the index at `idx`,
/// returning the data lines in the chunks that fall in the region, using
/// the columns the index header records (as `tabix` does).
fn query(data: &Path, idx: &Path, chrom: &str, start: usize, end: usize) -> Vec<String> {
    let kind = find_index(data).unwrap();
    assert_eq!(kind.path, idx);
    let (index, header) = match kind.kind {
        IndexKind::Tbi => {
            let i = noodles_tabix::fs::read(idx).unwrap();
            let h = i.header().unwrap().clone();
            (Box::new(i) as Box<dyn BinningIndex>, h)
        }
        IndexKind::Csi => {
            let i = noodles_csi::fs::read(idx).unwrap();
            let h = i.header().unwrap().clone();
            (Box::new(i) as Box<dyn BinningIndex>, h)
        }
    };
    let rid = header
        .reference_sequence_names()
        .get_index_of(chrom.as_bytes())
        .unwrap();
    let region: noodles_core::region::Interval =
        (Position::try_from(start).unwrap()..=Position::try_from(end).unwrap()).into();
    let chunks = index.query(rid, region).unwrap();
    let (cs, cb, ce) = (
        header.reference_sequence_name_index(),
        header.start_position_index(),
        header.end_position_index(),
    );
    let zero = header.format().coordinate_system()
        == noodles_csi::binning_index::index::header::format::CoordinateSystem::Bed;
    let mut out = Vec::new();
    for chunk in chunks {
        let mut r = BgzfReader::new(std::fs::File::open(data).unwrap());
        r.seek_to_virtual_position(chunk.start()).unwrap();
        loop {
            if r.virtual_position() >= chunk.end() {
                break;
            }
            let mut line = String::new();
            if r.read_line(&mut line).unwrap() == 0 {
                break;
            }
            let line = line.trim_end().to_string();
            if line.starts_with(header.line_comment_prefix() as char) {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f[cs] != chrom {
                continue;
            }
            let b: u64 = f[cb].parse().unwrap();
            let (b0, e0) = match (ce, zero) {
                (Some(ce), true) => (b, f[ce].parse::<u64>().unwrap()),
                (Some(ce), false) => (b - 1, f[ce].parse::<u64>().unwrap()),
                (None, true) => (b, b + 1),
                (None, false) => (b - 1, b),
            };
            if b0 < end as u64 && e0 > start as u64 - 1 {
                out.push(line);
            }
        }
    }
    out
}

#[test]
fn tsv_index_query_matches_lines() {
    let text = "chr1\t0\t50\tx\textra\n\
chr1\t90\t99\ty\textra\n\
chr1\t99\t120\tz\textra\n\
chr1\t150\t200\tw\textra\n\
chr1\t200\t300\tv\textra\n\
chr2\t99\t200\tu\textra\n";
    let p = bgzf(text);
    let idx = write_tsv_index(&p, &p, &cols(1, 2, Some(3), true), nz1()).unwrap();
    assert!(idx.to_string_lossy().ends_with(".tbi"));
    let got = query(&p, &idx, "chr1", 100, 200);
    assert_eq!(got, ["chr1\t99\t120\tz\textra", "chr1\t150\t200\tw\textra"]);
    let tbi = noodles_tabix::fs::read(&idx).unwrap();
    let h = tbi.header().unwrap();
    assert_eq!(h.reference_sequence_name_index(), 0);
    assert_eq!(h.start_position_index(), 1);
    assert_eq!(h.end_position_index(), Some(2));
    assert_eq!(h.line_comment_prefix(), b'#');
}

#[test]
fn one_based_without_end_column() {
    let text = "chrA\t10\tfoo\nchrA\t150\tbar\nchrA\t201\tbaz\n";
    let p = bgzf(text);
    let idx = write_tsv_index(&p, &p, &cols(1, 2, None, false), nz1()).unwrap();
    // begin 150 covers [149,150): 1-based position 150
    assert_eq!(query(&p, &idx, "chrA", 150, 150), ["chrA\t150\tbar"]);
    assert_eq!(query(&p, &idx, "chrA", 149, 149), Vec::<String>::new());
    assert_eq!(query(&p, &idx, "chrA", 100, 200), ["chrA\t150\tbar"]);
}

#[test]
fn skip_lines_and_comment() {
    let text = "this line is skipped\n% a comment\nchr1\t5\t10\nchr1\t20\t30\n";
    let p = bgzf(text);
    let mut c = cols(1, 2, Some(3), true);
    c.skip_lines = 1;
    c.comment = b'%';
    let idx = write_tsv_index(&p, &p, &c, nz1()).unwrap();
    let tbi = noodles_tabix::fs::read(&idx).unwrap();
    let h = tbi.header().unwrap();
    assert_eq!(h.line_skip_count(), 1);
    assert_eq!(h.line_comment_prefix(), b'%');
    assert_eq!(query(&p, &idx, "chr1", 1, 6), ["chr1\t5\t10"]);
}

#[test]
fn csi_for_large_end() {
    let text = "chr1\t10\t20\nchr1\t100\t600000000\n";
    let p = bgzf(text);
    let idx = write_tsv_index(&p, &p, &cols(1, 2, Some(3), true), nz1()).unwrap();
    assert!(idx.to_string_lossy().ends_with(".csi"), "{idx:?}");
    assert_eq!(query(&p, &idx, "chr1", 500_000_000, 500_000_001).len(), 1);
}

#[test]
fn unsorted_names_the_line() {
    let p = bgzf("chr1\t50\t60\nchr1\t10\t20\n");
    let err = write_tsv_index(&p, &p, &cols(1, 2, Some(3), true), nz1()).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("line 2"), "{msg}");
}

#[test]
fn bad_column_is_a_data_error_naming_the_line() {
    let p = bgzf("chr1\t10\t20\nchr1\tabc\t30\n");
    let err = write_tsv_index(&p, &p, &cols(1, 2, Some(3), true), nz1()).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("line 2") && msg.contains("begin"), "{msg}");
}
