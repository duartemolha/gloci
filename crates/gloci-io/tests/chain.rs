//! Chain files for liftover (D174): the toy chain, published layouts, compression, bad input.

use gloci_core::Interval;
use gloci_io::{ChainIndex, ChainInfo, Segment};
use std::io::Write;

const TOY: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../xtask/diff/inputs/liftover/toy.chain"
);

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

fn toy() -> ChainIndex {
    ChainIndex::load(TOY).unwrap()
}

fn pieces(ix: &ChainIndex, chrom: &str, s: u64, e: u64) -> Vec<Segment> {
    let src = ix.source(chrom).unwrap();
    let mut out = vec![];
    ix.segments(&src, iv(s, e), &mut out);
    out
}

/// `(chain, source, target, reverse)` as the plan writes them.
type Row = (u32, (u64, u64), (u64, u64), bool);

fn rows(v: &[Segment]) -> Vec<Row> {
    v.iter()
        .map(|s| {
            (
                s.chain,
                (s.source.start, s.source.end),
                (s.target.start, s.target.end),
                s.reverse,
            )
        })
        .collect()
}

#[test]
fn toy_chain_shape() {
    let ix = toy();
    assert_eq!(ix.chains(), 10);
    assert_eq!(ix.blocks(), 13);
    assert_eq!(
        ix.chain(1),
        &ChainInfo {
            line: 7,
            source: "chr1".into(),
            target: "chr2".into(),
            target_size: 5000,
            reverse: true
        }
    );
    assert_eq!(ix.chain(4).line, 17);
    assert_eq!(ix.chain(9).line, 32);
}

#[test]
fn segments_follow_crossmap_order() {
    let ix = toy();
    assert_eq!(
        rows(&pieces(&ix, "chr1", 0, 1000)),
        vec![
            (0, (100, 200), (1100, 1200), false),
            (0, (200, 300), (1250, 1350), false),
            (0, (350, 600), (1350, 1600), false),
            (2, (400, 500), (0, 100), false),
            (3, (400, 450), (10, 60), false),
            (1, (700, 780), (3920, 4000), true),
            (1, (800, 900), (3800, 3900), true),
        ]
    );
    assert_eq!(
        rows(&pieces(&ix, "chr1", 150, 150)),
        vec![(0, (150, 150), (1150, 1150), false)]
    );
    assert!(pieces(&ix, "chr1", 200, 200).is_empty());
    assert!(pieces(&ix, "chr1", 600, 600).is_empty());
    let p = pieces(&ix, "chr1", 420, 420);
    assert_eq!(p.iter().map(|s| s.chain).collect::<Vec<_>>(), [0, 2, 3]);
    assert_eq!(
        p.iter().map(|s| s.target.start).collect::<Vec<_>>(),
        [1420, 20, 30]
    );
    assert_eq!(
        rows(&pieces(&ix, "chrC", 50, 150)),
        vec![
            (8, (50, 100), (50, 100), false),
            (9, (100, 150), (850, 900), true)
        ]
    );
}

#[test]
fn segments_clears_out() {
    let ix = toy();
    let src = ix.source("chr1").unwrap();
    let mut out = vec![];
    ix.segments(&src, iv(0, 1000), &mut out);
    assert_eq!(out.len(), 7);
    ix.segments(&src, iv(200, 200), &mut out);
    assert!(out.is_empty());
}

#[test]
fn source_is_exact() {
    let ix = toy();
    assert!(ix.source("chr1").is_some());
    assert!(ix.source("5").is_some());
    assert!(ix.source("chr5").is_none());
    assert!(ix.source("1").is_none());
    assert!(ix.source("chrZ").is_none());
}

fn first(ix: &ChainIndex, chrom: &str, s: u64, e: u64) -> Row {
    rows(&pieces(ix, chrom, s, e))[0]
}

#[test]
fn published_layouts_load() {
    // Ensembl: space-separated blocks, the same id under `1` and `chr1`.
    let ens = "chain 1 1 1000 + 10 110 1 1000 + 20 120 2\n60 10 10\n30\n\n\
               chain 1 chr1 1000 + 10 110 1 1000 + 20 120 2\n60 10 10\n30\n\n";
    let ix = ChainIndex::parse("e.chain", ens.as_bytes()).unwrap();
    assert_eq!(ix.chains(), 2);
    assert_eq!(first(&ix, "1", 15, 16), (0, (15, 16), (25, 26), false));
    assert_eq!(first(&ix, "chr1", 15, 16), (1, (15, 16), (25, 26), false));

    // T2T: tab-separated header, id 0 twice.
    let t2t = "chain\t100\tchrA\t1000\t+\t0\t100\tchrX\t1000\t+\t0\t100\t0\n100\n\n\
               chain\t100\tchrB\t1000\t+\t0\t100\tchrX\t1000\t+\t100\t200\t0\n100\n\n";
    let ix = ChainIndex::parse("t.chain", t2t.as_bytes()).unwrap();
    assert_eq!(ix.chains(), 2);
    assert_eq!(first(&ix, "chrB", 5, 6), (1, (5, 6), (105, 106), false));

    // A 12-field header.
    let ix = ChainIndex::parse("c", b"chain 1 a 100 + 0 50 b 100 + 10 60\n50\n").unwrap();
    assert_eq!(first(&ix, "a", 0, 5), (0, (0, 5), (10, 15), false));

    // CRLF line ends.
    let ix = ChainIndex::parse("c", b"chain 1 a 100 + 0 50 b 100 + 10 60 1\r\n50\r\n\r\n").unwrap();
    assert_eq!(first(&ix, "a", 0, 5), (0, (0, 5), (10, 15), false));

    // A `#` comment between chains.
    let text = "chain 1 a 100 + 0 50 b 100 + 10 60 1\n50\n\n# comment\n\
                chain 1 a 100 + 60 100 b 100 + 0 40 2\n40\n";
    let ix = ChainIndex::parse("c", text.as_bytes()).unwrap();
    assert_eq!((ix.chains(), ix.blocks()), (2, 2));

    // No blank line between chains.
    let text =
        "chain 1 a 100 + 0 50 b 100 + 10 60 1\n50\nchain 1 a 100 + 60 100 b 100 + 0 40 2\n40\n";
    let ix = ChainIndex::parse("c", text.as_bytes()).unwrap();
    assert_eq!((ix.chains(), ix.blocks()), (2, 2));
    assert_eq!(first(&ix, "a", 70, 71), (1, (70, 71), (10, 11), false));
}

#[test]
fn compressed_chains_load_the_same() {
    let plain = std::fs::read(TOY).unwrap();
    let dir = tempfile::tempdir().unwrap();

    let gz_path = dir.path().join("toy.chain.gz");
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&plain).unwrap();
    std::fs::write(&gz_path, gz.finish().unwrap()).unwrap();

    let bgzf_path = dir.path().join("toy.chain.bgz");
    let mut w = noodles_bgzf::io::Writer::new(std::fs::File::create(&bgzf_path).unwrap());
    w.write_all(&plain).unwrap();
    w.finish().unwrap();

    let want = toy();
    let want_rows = rows(&pieces(&want, "chr1", 0, 1000));
    for p in [gz_path, bgzf_path] {
        let ix = ChainIndex::load(p.to_str().unwrap()).unwrap();
        assert_eq!(ix.chains(), want.chains());
        assert_eq!(ix.blocks(), want.blocks());
        assert_eq!(rows(&pieces(&ix, "chr1", 0, 1000)), want_rows);
    }
}

#[test]
fn bad_chains_name_the_line() {
    let h = "chain 1 chr1 1000 + 0 100 chr2 1000 + 0 100 1\n";
    // (text, line the error names (0: not required), needles)
    let cases: Vec<(String, u64, Vec<&str>)> = vec![
        (
            "chain 1 chr1 1000 - 0 100 chr2 1000 + 0 100 1\n100\n".into(),
            1,
            vec!["source strand"],
        ),
        (
            "chain 1 chr1 1000 + 0 100 chr2 1000 ? 0 100 1\n100\n".into(),
            1,
            vec!["target strand"],
        ),
        (
            "chain 1 chr1 1000 + 0 100 chr2 1000 + 0\n100\n".into(),
            1,
            vec!["12 or 13"],
        ),
        ("100\n".into(), 1, vec!["before"]),
        (format!("{h}60 10\n"), 2, vec!["1 or 3 numbers"]),
        (format!("{h}60 10 10\n{h}"), 3, vec!["ends"]),
        (
            "chain 1 chr1 1000 + 0 650 chr2 1000 + 0 600 1\n600\n".into(),
            2,
            vec!["600", "650"],
        ),
        (
            "chain 1 chr1 1000 + 0 1100 chr2 2000 + 0 100 1\n100\n".into(),
            1,
            vec!["size"],
        ),
        (
            "chain 1 chr1 1000 + 0 10 chr2 1000 + 0 10 1\nx\n".into(),
            2,
            vec!["number"],
        ),
        (
            format!("{h}100\nchain 1 chr1 1001 + 0 100 chr2 1000 + 0 100 2\n100\n"),
            3,
            vec!["line 1"],
        ),
        (
            format!("{h}100\nchain 1 chr3 1000 + 0 100 chr2 999 + 0 100 2\n100\n"),
            3,
            vec!["line 1"],
        ),
        ("# only\n# comments\n".into(), 0, vec!["no chain"]),
        (format!("{h}60 10 10\n"), 2, vec!["ends inside"]),
    ];
    for (text, line, needles) in cases {
        let err = ChainIndex::parse("c.chain", text.as_bytes())
            .err()
            .unwrap_or_else(|| panic!("{text:?} should fail"))
            .to_string();
        for needle in needles {
            assert!(
                err.contains(needle),
                "{err:?} lacks {needle:?} for {text:?}"
            );
        }
        if line > 0 {
            assert!(
                err.contains(&format!("c.chain line {line}")),
                "{err:?} lacks line {line} for {text:?}"
            );
        }
    }
}

#[test]
fn every_fixture_chain_loads() {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../xtask/diff/inputs/liftover/"
    );
    ChainIndex::load(&format!("{dir}toy.chain")).unwrap();
    for (file, line, needle) in [
        ("bad-end.chain", "line 1", "650"),
        ("minus-source.chain", "line 1", "source strand"),
    ] {
        let msg = ChainIndex::load(&format!("{dir}{file}"))
            .err()
            .unwrap_or_else(|| panic!("{file} must not load"))
            .to_string();
        assert!(msg.contains(line) && msg.contains(needle), "{file}: {msg}");
    }
}
