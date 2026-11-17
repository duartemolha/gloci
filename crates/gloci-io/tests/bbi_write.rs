//! The bigWig and bigBed writer (D171): checked input, atomic output, the same bytes at any
//! thread count.

use bigtools::{BigBedRead, BigWigRead, Value};
use gloci_core::Genome;
use gloci_io::bbi::autosql::{
    autosql_field_count, autosql_standard_columns, bed_autosql, bed_layout,
};
use gloci_io::{
    BbiKind, BbiSizes, BbiWriteOptions, BbiWriter, IoError, OutFormat, OutputSpec, PushError,
    parse_output_spec, parse_track_output_spec,
};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::mpsc;
use std::time::Duration;

fn genome(list: &[(&str, u64)]) -> Genome {
    Genome::new(
        list.iter()
            .map(|(n, l)| (n.as_bytes().to_vec(), Some(*l)))
            .collect(),
    )
    .unwrap()
}

fn options(kind: BbiKind, sizes: &[(&str, u64)], threads: usize) -> BbiWriteOptions {
    BbiWriteOptions {
        kind,
        sizes: BbiSizes {
            genome: genome(sizes),
            source: "-g test.sizes".into(),
        },
        threads: NonZeroUsize::new(threads).unwrap(),
        autosql: None,
        autosql_name: "autosql".into(),
        defined_field_count: None,
    }
}

fn spec(path: &Path) -> OutputSpec {
    parse_track_output_spec(Some(path), None).unwrap()
}

fn writer(path: &Path, kind: BbiKind, sizes: &[(&str, u64)], threads: usize) -> BbiWriter {
    BbiWriter::create(&spec(path), options(kind, sizes, threads)).unwrap()
}

const SIZES: &[(&str, u64)] = &[("chr1", 1000), ("chr2", 1000), ("chr10", 1000)];

/// The record error pushing `lines` ends with, none of the earlier lines being refused.
fn record_error(kind: BbiKind, sizes: &[(&str, u64)], lines: &[&str]) -> String {
    let dir = tempfile::tempdir().unwrap();
    let name = if kind == BbiKind::BigWig {
        "x.out.bw"
    } else {
        "x.out.bb"
    };
    let mut w = writer(&dir.path().join(name), kind, sizes, 1);
    let (last, init) = lines.split_last().unwrap();
    for l in init {
        w.push(l.as_bytes()).unwrap();
    }
    match w.push(last.as_bytes()) {
        Err(PushError::Record(m)) => m,
        other => panic!("expected a record error for {last:?}, got {other:?}"),
    }
}

fn bw_error(lines: &[&str]) -> String {
    record_error(BbiKind::BigWig, SIZES, lines)
}

fn bb_error(lines: &[&str]) -> String {
    record_error(BbiKind::BigBed, SIZES, lines)
}

fn write_all(path: &Path, kind: BbiKind, sizes: &[(&str, u64)], threads: usize, lines: &[String]) {
    let mut w = writer(path, kind, sizes, threads);
    for l in lines {
        w.push(l.as_bytes()).unwrap();
    }
    w.finish().unwrap();
}

#[test]
fn output_spec_rules() {
    for name in ["x.bw", "x.bigwig"] {
        let e = parse_output_spec(Some(Path::new(name)), None)
            .unwrap_err()
            .0;
        assert!(
            e.contains("gloci convert") && e.contains("genomecov --bedgraph"),
            "{e}"
        );
    }
    for name in ["x.bb", "x.bigbed"] {
        let e = parse_output_spec(Some(Path::new(name)), None)
            .unwrap_err()
            .0;
        assert!(e.contains("gloci convert"), "{e}");
    }
    assert!(parse_output_spec(None, Some("bigwig")).is_err());
    let fmt = |p: Option<&str>, o: Option<&str>| {
        parse_track_output_spec(p.map(Path::new), o).unwrap().format
    };
    assert_eq!(fmt(Some("x.bw"), None), Some(OutFormat::BigWig));
    assert_eq!(fmt(Some("x.BigWig"), None), Some(OutFormat::BigWig));
    assert_eq!(fmt(None, Some("bigwig")), Some(OutFormat::BigWig));
    assert_eq!(fmt(Some("x.bb"), None), Some(OutFormat::BigBed));
    assert_eq!(fmt(Some("x.bigbed"), None), Some(OutFormat::BigBed));
    assert!(parse_track_output_spec(Some(Path::new("x.bw")), Some("bigbed")).is_err());
}

#[test]
fn refuses_non_file_outputs() {
    let sizes = || options(BbiKind::BigWig, SIZES, 1);
    let usage = |r: Result<BbiWriter, IoError>| match r {
        Err(IoError::Usage(u)) => assert!(u.0.contains("regular file"), "{}", u.0),
        other => panic!("expected a usage error, got {:?}", other.map(|_| ())),
    };
    let base = OutputSpec {
        path: None,
        format: Some(OutFormat::BigWig),
        bgzf: false,
        direct: false,
    };
    usage(BbiWriter::create(&base, sizes()));
    let dev = OutputSpec {
        path: Some("/dev/null".into()),
        direct: true,
        ..base.clone()
    };
    usage(BbiWriter::create(&dev, sizes()));
    let mem = OutputSpec {
        path: Some("mem:1".into()),
        direct: true,
        ..base
    };
    usage(BbiWriter::create(&mem, sizes()));
}

#[test]
fn bigwig_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("a.bw");
    let lines: Vec<String> = ["chr2\t0\t10\t0.5", "chr10\t0\t5\t3", "chr1\t0\t1\t-1"]
        .map(String::from)
        .to_vec();
    write_all(&out, BbiKind::BigWig, SIZES, 2, &lines);
    let mut r = BigWigRead::open_file(&out).unwrap();
    let names: Vec<String> = r.chroms().iter().map(|c| c.name.clone()).collect();
    assert_eq!(names, ["chr2", "chr10", "chr1"]);
    let get = |r: &mut BigWigRead<_>, c: &str, e: u32| -> Vec<Value> {
        r.get_interval(c, 0, e)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap()
    };
    assert_eq!(get(&mut r, "chr2", 10)[0].value, 0.5);
    assert_eq!(get(&mut r, "chr10", 5)[0].value, 3.0);
    assert_eq!(get(&mut r, "chr1", 1)[0].value, -1.0);
}

fn many_records() -> Vec<String> {
    let mut lines = Vec::new();
    for c in ["chr1", "chr2", "chr10"] {
        for i in 0..1000u64 {
            lines.push(format!(
                "{c}\t{}\t{}\t{}",
                i * 3,
                i * 3 + 2,
                (i % 17) as f32 * 0.25
            ));
        }
    }
    lines
}

#[test]
fn bytes_do_not_depend_on_threads() {
    let dir = tempfile::tempdir().unwrap();
    let sizes = [("chr1", 100_000), ("chr2", 100_000), ("chr10", 100_000)];
    let lines = many_records();
    let mut files = Vec::new();
    for t in [1, 2, 8] {
        let p = dir.path().join(format!("t{t}.bw"));
        write_all(&p, BbiKind::BigWig, &sizes, t, &lines);
        files.push(std::fs::read(&p).unwrap());
    }
    assert_eq!(files[0], files[1], "-t 1 and -t 2 differ");
    assert_eq!(files[0], files[2], "-t 1 and -t 8 differ");
    let bed: Vec<String> = lines
        .iter()
        .map(|l| l.rsplit_once('\t').unwrap().0.to_string() + "\tname\t5")
        .collect();
    let mut files = Vec::new();
    for t in [1, 2, 8] {
        let p = dir.path().join(format!("t{t}.bb"));
        write_all(&p, BbiKind::BigBed, &sizes, t, &bed);
        files.push(std::fs::read(&p).unwrap());
    }
    assert_eq!(files[0], files[1], "bigBed -t 1 and -t 2 differ");
    assert_eq!(files[0], files[2], "bigBed -t 1 and -t 8 differ");
}

#[test]
fn split_chromosome_is_an_error() {
    let e = bw_error(&["chr1\t0\t1\t1", "chr2\t0\t1\t1", "chr1\t5\t6\t1"]);
    assert!(e.contains("chr1") && e.contains("gloci sort"), "{e}");
}

#[test]
fn starts_must_ascend() {
    let e = bw_error(&["chr1\t50\t60\t1", "chr1\t10\t20\t1"]);
    assert!(e.contains("gloci sort"), "{e}");
}

#[test]
fn bigwig_overlap_is_an_error() {
    let e = bw_error(&["chr1\t0\t10\t1", "chr1\t5\t20\t2"]);
    assert!(e.contains("overlap") && e.contains("chr1:0-10"), "{e}");
}

#[test]
fn bigwig_row_shape() {
    for l in ["chr1\t0\t10", "chr1\t0\t10\t1\tx"] {
        let e = bw_error(&[l]);
        assert!(e.contains("4 columns"), "{e}");
    }
    let e = bw_error(&["chr1\t0\t10\tgeneA"]);
    assert!(e.contains("geneA") && e.contains("not a number"), "{e}");
    let e = bw_error(&["chr1\t5\t5\t1"]);
    assert!(e.contains("zero-length"), "{e}");
}

#[test]
fn value_outside_float32_is_an_error() {
    for v in ["1e39", "-1e39", "nan", "inf", "-inf"] {
        let e = bw_error(&[&format!("chr1\t0\t10\t{v}")]);
        assert!(e.contains(v), "{e}");
    }
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("f.bw");
    write_all(
        &out,
        BbiKind::BigWig,
        SIZES,
        1,
        &["chr1\t0\t10\t0.123456789".to_string()],
    );
    let mut r = BigWigRead::open_file(&out).unwrap();
    let v: Vec<Value> = r
        .get_interval("chr1", 0, 10)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(v[0].value, 0.123_456_79_f32);
}

#[test]
fn end_past_chromosome_is_an_error() {
    let sizes = [("chr1", 1000)];
    let e = record_error(BbiKind::BigWig, &sizes, &["chr1\t900\t1001\t1"]);
    assert!(e.contains("1000"), "{e}");
    let e = record_error(BbiKind::BigBed, &sizes, &["chr1\t900\t1001\tx"]);
    assert!(e.contains("1000"), "{e}");
}

#[test]
fn unknown_chromosome_names_the_sizes_source() {
    for kind in [BbiKind::BigWig, BbiKind::BigBed] {
        let line = if kind == BbiKind::BigWig {
            "chrZ\t0\t10\t1"
        } else {
            "chrZ\t0\t10\tx"
        };
        let e = record_error(kind, SIZES, &[line]);
        assert!(e.contains("chrZ") && e.contains("-g test.sizes"), "{e}");
    }
    let dir = tempfile::tempdir().unwrap();
    let g = Genome::new(vec![(b"chr1".to_vec(), None)]).unwrap();
    let mut o = options(BbiKind::BigWig, SIZES, 1);
    o.sizes.genome = g;
    let mut w = BbiWriter::create(&spec(&dir.path().join("n.bw")), o).unwrap();
    match w.push(b"chr1\t0\t10\t1") {
        Err(PushError::Record(m)) => assert!(m.contains("length"), "{m}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn chromosome_longer_than_32_bits_is_an_error() {
    let e = record_error(
        BbiKind::BigWig,
        &[("chrBig", 5_000_000_000)],
        &["chrBig\t0\t10\t1"],
    );
    assert!(e.contains("32-bit"), "{e}");
}

#[test]
fn too_many_chromosomes_is_an_error() {
    let names: Vec<(String, u64)> = (0..65_536).map(|i| (format!("s{i}"), 10)).collect();
    let sizes: Vec<(&str, u64)> = names.iter().map(|(n, l)| (n.as_str(), *l)).collect();
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(&dir.path().join("many.bw"), BbiKind::BigWig, &sizes, 2);
    for (i, (n, _)) in names.iter().enumerate() {
        match w.push(format!("{n}\t0\t5\t1").as_bytes()) {
            Ok(()) => assert!(i < 65_535, "record {i} was accepted"),
            Err(PushError::Record(m)) => {
                assert_eq!(i, 65_535);
                assert!(m.contains("65,535"), "{m}");
                return;
            }
            Err(e) => panic!("{e}"),
        }
    }
    panic!("no error");
}

fn temp_files(dir: &Path) -> Vec<String> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect()
}

#[test]
fn empty_input_is_an_error_and_leaves_no_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("e.bw");
    let w = writer(&out, BbiKind::BigWig, SIZES, 2);
    let e = w.finish().unwrap_err().to_string();
    assert!(e.contains("no records"), "{e}");
    assert!(!out.exists());
    assert_eq!(temp_files(dir.path()), Vec::<String>::new());
}

#[test]
fn drop_without_finish_removes_the_temp_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("d.bw");
    let mut w = writer(&out, BbiKind::BigWig, SIZES, 2);
    for i in 0..10 {
        w.push(format!("chr1\t{}\t{}\t1", i * 2, i * 2 + 1).as_bytes())
            .unwrap();
    }
    drop(w);
    assert!(!out.exists());
    assert_eq!(temp_files(dir.path()), Vec::<String>::new());
}

#[test]
fn bigbed_layout_rules() {
    let rows: [(&str, usize, usize); 9] = [
        ("c\t0\t1", 3, 3),
        ("c\t0\t1\tn", 4, 4),
        ("c\t0\t1\tn\tx", 4, 5),
        ("c\t0\t1\tn\t5\tx", 5, 6),
        ("c\t0\t1\tn\t5\t+\t12.5\tx\ty\tz", 6, 10),
        ("c\t0\t1\tn\t5\t+\t0", 6, 7),
        ("c\t0\t10\tn\t5\t+\t0\t10\t255,0,0", 9, 9),
        ("c\t0\t10\tn\t5\t+\t0\t10\t0\t2\t3,4,\t0,6,", 12, 12),
        ("c\t0\t10\tn\t5\t+\t0\t10\t0\t3\t3,4,\t0,6,", 9, 12),
    ];
    for (line, standard, columns) in rows {
        let l = bed_layout(line.as_bytes()).unwrap();
        assert_eq!((l.standard, l.columns), (standard, columns), "{line}");
        assert_eq!(
            autosql_field_count(&bed_autosql(&l)).unwrap(),
            columns,
            "{line}"
        );
    }
}

#[test]
fn bigbed_autosql_and_header() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("a.bb");
    write_all(
        &out,
        BbiKind::BigBed,
        SIZES,
        2,
        &["chr1\t0\t10\tn\t5\t+\t1.5\tx".to_string()],
    );
    let mut r = BigBedRead::open_file(&out).unwrap();
    assert_eq!(
        r.autosql().unwrap().unwrap(),
        "table bed\n\"BED records written by gloci\"\n(\n    string chrom; \"Chromosome\"\n    uint chromStart; \"Start, 0-based\"\n    uint chromEnd; \"End, exclusive\"\n    string name; \"Name\"\n    uint score; \"Score\"\n    char[1] strand; \"Strand: +, - or .\"\n    string field1; \"Extra column 1\"\n    string field2; \"Extra column 2\"\n)\n"
    );
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(u16::from_le_bytes([bytes[32], bytes[33]]), 8);
    assert_eq!(u16::from_le_bytes([bytes[34], bytes[35]]), 6);
}

#[test]
fn bigbed_standard_rows() {
    let l = bed_layout(b"c\t0\t10\tn\t5\t+\t0\t10\t0\t2\t3,4,\t0,6,").unwrap();
    let text = bed_autosql(&l);
    for row in [
        "    uint thickStart; \"Start of the thick part\"\n",
        "    uint thickEnd; \"End of the thick part\"\n",
        "    uint itemRgb; \"Colour, R,G,B\"\n",
        "    int blockCount; \"Number of blocks\"\n",
        "    int[blockCount] blockSizes; \"Block sizes\"\n",
        "    int[blockCount] blockStarts; \"Block starts, relative to chromStart\"\n",
    ] {
        assert!(text.contains(row), "{row}");
    }
}

#[test]
fn bigbed_records_must_keep_the_layout() {
    let e = bb_error(&["chr1\t0\t10\tn\t5\t+\t1.5\tx", "chr1\t20\t30\tn\t5\t+\t1.5"]);
    assert!(e.contains("8 columns"), "{e}");
    let e = bb_error(&["chr1\t0\t10\tn\t5", "chr1\t20\t30\tn\tx"]);
    assert!(e.contains("score"), "{e}");
}

#[test]
fn bigbed_given_autosql_must_match() {
    let dir = tempfile::tempdir().unwrap();
    let mut o = options(BbiKind::BigBed, SIZES, 1);
    o.autosql = Some(bed_autosql(&bed_layout(b"c\t0\t1\tn\t5").unwrap()));
    let mut w = BbiWriter::create(&spec(&dir.path().join("g.bb")), o).unwrap();
    match w.push(b"chr1\t0\t10\tn\t5\t+\t1.5\tx") {
        Err(PushError::Record(m)) => assert!(m.contains('5') && m.contains('8'), "{m}"),
        other => panic!("{other:?}"),
    }
}

/// An autoSql table of `fields` (`type name` each).
fn table(fields: &[&str]) -> String {
    let mut t = String::from("table t\n\"t\"\n(\n");
    for f in fields {
        t.push_str(&format!("    {f}; \"x\"\n"));
    }
    t.push_str(")\n");
    t
}

const BED6: [&str; 6] = [
    "string chrom",
    "uint chromStart",
    "uint chromEnd",
    "string name",
    "uint score",
    "char[1] strand",
];

#[test]
fn autosql_standard_columns_go_by_name_and_type() {
    let with = |extra: &[&str]| {
        let mut f = BED6.to_vec();
        f.extend(extra);
        autosql_standard_columns(&table(&f)).unwrap()
    };
    assert_eq!(with(&["uint readsA", "uint readsB", "uint replicates"]), 6);
    assert_eq!(with(&["float signalValue", "float pValue"]), 6);
    assert_eq!(
        with(&["uint thickStart", "uint thickEnd", "uint itemRgb"]),
        9
    );
    assert_eq!(
        with(&[
            "uint thickStart",
            "uint thickEnd",
            "uint reserved",
            "int blockCount",
            "int[blockCount] blockSizes",
            "int[blockCount] chromStarts"
        ]),
        12
    );
    assert_eq!(
        with(&[
            "uint thickStart",
            "uint thickEnd",
            "uint itemRgb",
            "int blockCount",
            "int[blockCount] blockSizes",
            "int[blockCount] blockStarts",
            "string more"
        ]),
        12
    );
    // thickStart without thickEnd: BED has no 7-column prefix.
    assert_eq!(with(&["uint thickStart", "uint other"]), 6);
    // blockCount without its lists: rounded down to 9.
    assert_eq!(
        with(&[
            "uint thickStart",
            "uint thickEnd",
            "uint itemRgb",
            "int blockCount"
        ]),
        9
    );
    // The type must match too.
    assert_eq!(with(&["string thickStart", "uint thickEnd"]), 6);
    let fields = [
        "string chrom",
        "uint chromStart",
        "uint chromEnd",
        "string name",
        "float score",
    ];
    assert_eq!(autosql_standard_columns(&table(&fields)).unwrap(), 4);
    let fields = [
        "string chrom",
        "uint chromStart",
        "uint chromEnd",
        "string name",
        "uint score",
        "string strand",
    ];
    assert_eq!(autosql_standard_columns(&table(&fields)).unwrap(), 5);
    // The first three columns are chrom, chromStart and chromEnd whatever their names.
    let fields = ["string chr", "uint start", "uint end", "string name"];
    assert_eq!(autosql_standard_columns(&table(&fields)).unwrap(), 4);
}

/// The definedFieldCount of the bigBed at `path`.
fn defined_field_count(path: &Path) -> u16 {
    BigBedRead::open_file(path)
        .unwrap()
        .info()
        .header
        .defined_field_count
}

#[test]
fn the_given_autosql_sets_the_standard_columns() {
    let dir = tempfile::tempdir().unwrap();
    // Integer extras would pass as thickStart, thickEnd and itemRgb.
    let line = b"chr1\t0\t10\tn\t5\t+\t12\t30\t2";
    let mut extra = BED6.to_vec();
    extra.extend(["uint readsA", "uint readsB", "uint replicates"]);
    let mut o = options(BbiKind::BigBed, SIZES, 1);
    o.autosql = Some(table(&extra));
    let out = dir.path().join("a.bb");
    let mut w = BbiWriter::create(&spec(&out), o).unwrap();
    w.push(line).unwrap();
    // Later records need only BED6: column 9 is no colour, columns 7 and 8 no thick part.
    w.push(b"chr1\t20\t30\tn\t5\t+\t12\t30\tx,y").unwrap();
    w.push(b"chr1\t20\t30\tn\t5\t+\t40\t2\t300").unwrap();
    w.finish().unwrap();
    assert_eq!(defined_field_count(&out), 6);
    // Without the autoSql the first record decides: 9.
    let plain = dir.path().join("p.bb");
    write_all(
        &plain,
        BbiKind::BigBed,
        SIZES,
        1,
        &[String::from_utf8(line.to_vec()).unwrap()],
    );
    assert_eq!(defined_field_count(&plain), 9);
}

#[test]
fn records_must_have_the_standard_columns_the_autosql_declares() {
    let dir = tempfile::tempdir().unwrap();
    let mut fields = BED6.to_vec();
    fields.extend(["uint thickStart", "uint thickEnd"]);
    let mut o = options(BbiKind::BigBed, SIZES, 1);
    o.autosql = Some(table(&fields));
    o.autosql_name = "my.as".into();
    let mut w = BbiWriter::create(&spec(&dir.path().join("d.bb")), o).unwrap();
    match w.push(b"chr1\t0\t10\tn\t5\t+\t1.5\tx") {
        Err(PushError::Record(m)) => {
            assert!(m.contains("thickStart '1.5'"), "{m}");
            assert!(m.contains("my.as declares the first 8 columns"), "{m}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_kept_defined_field_count_is_capped_at_the_first_record() {
    let dir = tempfile::tempdir().unwrap();
    for (kept, want) in [(6, 6), (12, 9), (7, 6)] {
        let mut o = options(BbiKind::BigBed, SIZES, 1);
        o.defined_field_count = Some(kept);
        let out = dir.path().join(format!("k{kept}.bb"));
        let mut w = BbiWriter::create(&spec(&out), o).unwrap();
        w.push(b"chr1\t0\t10\tn\t5\t+\t12\t30\t2").unwrap();
        w.finish().unwrap();
        assert_eq!(defined_field_count(&out), want, "kept {kept}");
    }
}

#[test]
fn bigbed_bytes() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(&dir.path().join("u.bb"), BbiKind::BigBed, SIZES, 1);
    match w.push(b"chr1\t0\t10\tx\xE9") {
        Err(PushError::Record(m)) => assert!(m.contains("UTF-8"), "{m}"),
        other => panic!("{other:?}"),
    }
}

/// Runs `body` on a thread and fails, rather than hangs, when it does not finish in time.
fn within_seconds<T: Send + 'static>(secs: u64, body: impl FnOnce() -> T + Send + 'static) -> T {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(body());
    });
    rx.recv_timeout(Duration::from_secs(secs))
        .expect("the write did not finish in time: it hangs")
}

fn guarded_write(kind: BbiKind, sizes: &'static [(&'static str, u64)], lines: &'static [&str]) {
    within_seconds(20, move || {
        let dir = tempfile::tempdir().unwrap();
        write_all(
            &dir.path().join(if kind == BbiKind::BigWig {
                "t.bw"
            } else {
                "t.bb"
            }),
            kind,
            sizes,
            2,
            &lines.iter().map(|l| l.to_string()).collect::<Vec<_>>(),
        );
    });
}

/// bigtools' writer never finishes on a file whose records all have zero length: gloci stops
/// it and says so. Without that, this test would wait forever (it fails after 10 seconds).
#[test]
fn bigbed_with_only_zero_length_records_is_an_error_not_a_hang() {
    for lines in [
        &["chr1\t5\t5\tx"][..],
        &["chr1\t5\t5\tx", "chr1\t5\t5\ty"][..],
        &["chr1\t5\t5\tx", "chr2\t7\t7\ty"][..],
    ] {
        let lines: Vec<String> = lines.iter().map(|l| l.to_string()).collect();
        let (dir, out) = {
            let dir = tempfile::tempdir().unwrap();
            let out = dir.path().join("z.bb");
            (dir, out)
        };
        let out2 = out.clone();
        let e = within_seconds(10, move || {
            let mut w = writer(&out2, BbiKind::BigBed, SIZES, 2);
            for l in &lines {
                w.push(l.as_bytes()).unwrap();
            }
            w.finish().unwrap_err().to_string()
        });
        assert!(e.contains("zero length"), "{e}");
        assert!(!out.exists());
        assert_eq!(temp_files(dir.path()), Vec::<String>::new());
    }
}

#[test]
fn bigbed_zero_length_records_beside_others_read_back() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("m.bb");
    let lines: Vec<String> = [
        "chr1\t1\t10\ta",
        "chr1\t5\t5\tb",
        "chr1\t20\t20\tc",
        "chr2\t7\t7\td",
        "chr2\t8\t9\te",
    ]
    .map(String::from)
    .to_vec();
    let (o, l) = (out.clone(), lines.clone());
    within_seconds(20, move || write_all(&o, BbiKind::BigBed, SIZES, 2, &l));
    let mut r = BigBedRead::open_file(&out).unwrap();
    let mut got = Vec::new();
    for c in ["chr1", "chr2"] {
        for e in r.get_interval(c, 0, 1000).unwrap() {
            let e = e.unwrap();
            got.push(format!("{c}\t{}\t{}\t{}", e.start, e.end, e.rest));
        }
    }
    assert_eq!(got, lines);
}

#[test]
fn bigbed_zero_length_refused_at_zero_and_at_the_chromosome_end() {
    let e = bb_error(&["chr1\t0\t0\tx"]);
    assert!(e.contains("position 0"), "{e}");
    let e = bb_error(&["chr1\t1000\t1000\tx"]);
    assert!(e.contains("1000") && e.contains("very end"), "{e}");
}

#[test]
fn tiny_chromosomes_finish() {
    guarded_write(BbiKind::BigWig, &[("chr1", 10)], &["chr1\t1\t3\t1"]);
    guarded_write(BbiKind::BigBed, &[("chr1", 10)], &["chr1\t1\t3\tx"]);
}

#[test]
fn bad_autosql_is_a_usage_error_naming_the_file_and_leaves_no_file() {
    let dir = tempfile::tempdir().unwrap();
    let mut o = options(BbiKind::BigBed, SIZES, 1);
    o.autosql = Some("table bed\n(\n  nonsense\n)".into());
    o.autosql_name = "bad.as".into();
    match BbiWriter::create(&spec(&dir.path().join("b.bb")), o) {
        Err(IoError::Usage(u)) => {
            assert!(u.0.starts_with("bad.as: "), "{}", u.0);
            assert!(!u.0.contains("Invalid"), "Display form, not Debug: {}", u.0);
        }
        other => panic!("expected a usage error, got {:?}", other.map(|_| ())),
    }
    assert_eq!(temp_files(dir.path()), Vec::<String>::new());
    assert!(!dir.path().join("b.bb").exists());
}

#[test]
fn a_refused_first_bigbed_record_does_not_set_the_layout() {
    let dir = tempfile::tempdir().unwrap();
    let mut w = writer(&dir.path().join("l.bb"), BbiKind::BigBed, SIZES, 1);
    // Eight columns on a chromosome the sizes lack: refused.
    assert!(matches!(
        w.push(b"nochr\t0\t10\tn\t5\t+\t1.5\tx"),
        Err(PushError::Record(_))
    ));
    // The layout is still open: three columns are accepted.
    w.push(b"chr1\t0\t10").unwrap();
    w.finish().unwrap();
}

#[test]
fn a_user_autosql_is_stored_as_given() {
    let dir = tempfile::tempdir().unwrap();
    let text = "table peak\n\"my peaks\"\n(\n    string chrom; \"c\"\n    uint chromStart; \"s\"\n    uint chromEnd; \"e\"\n    string name; \"n\"\n    uint score; \"sc\"\n)\n";
    let mut o = options(BbiKind::BigBed, SIZES, 2);
    o.autosql = Some(text.into());
    let out = dir.path().join("u.bb");
    let mut w = BbiWriter::create(&spec(&out), o).unwrap();
    w.push(b"chr1\t0\t10\tn\t5").unwrap();
    w.finish().unwrap();
    let mut r = BigBedRead::open_file(&out).unwrap();
    assert_eq!(r.autosql().unwrap().unwrap(), text);
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(u16::from_le_bytes([bytes[34], bytes[35]]), 5);
}

#[test]
fn drop_after_batches_were_sent_removes_the_temp_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("big.bw");
    let sizes = [("chr1", 1_000_000u64)];
    let mut w = writer(&out, BbiKind::BigWig, &sizes, 2);
    // More than QUEUE batches of 4,096: the writer thread has been fed and is mid-write.
    for i in 0..(4096 * 4 + 100) {
        w.push(format!("chr1\t{}\t{}\t1", i * 2, i * 2 + 1).as_bytes())
            .unwrap();
    }
    drop(w);
    assert!(!out.exists());
    assert_eq!(temp_files(dir.path()), Vec::<String>::new());
}

#[test]
fn bytes_do_not_depend_on_threads_across_batches() {
    let dir = tempfile::tempdir().unwrap();
    let sizes = [("chr1", 200_000u64), ("chr2", 200_000), ("chr10", 200_000)];
    let mut lines = Vec::new();
    for (c, n) in [("chr1", 8000u64), ("chr2", 7000), ("chr10", 5000)] {
        for i in 0..n {
            lines.push(format!(
                "{c}\t{}\t{}\t{}",
                i * 4,
                i * 4 + 3,
                (i % 23) as f32 * 0.5
            ));
        }
    }
    assert_eq!(lines.len(), 20_000);
    for kind in [BbiKind::BigWig, BbiKind::BigBed] {
        let ext = if kind == BbiKind::BigWig { "bw" } else { "bb" };
        let lines: Vec<String> = if kind == BbiKind::BigBed {
            lines
                .iter()
                .map(|l| l.rsplit_once('\t').unwrap().0.to_string() + "\tname")
                .collect()
        } else {
            lines.clone()
        };
        let mut files = Vec::new();
        for t in [1, 2, 8] {
            let p = dir.path().join(format!("t{t}.{ext}"));
            write_all(&p, kind, &sizes, t, &lines);
            files.push(std::fs::read(&p).unwrap());
        }
        assert_eq!(files[0], files[1], "{ext}: -t 1 and -t 2 differ");
        assert_eq!(files[0], files[2], "{ext}: -t 1 and -t 8 differ");
    }
}

#[test]
fn only_the_track_commands_list_bbi_extensions() {
    let plain = parse_output_spec(Some(Path::new("x.zzz")), None)
        .unwrap_err()
        .0;
    assert!(
        !plain.contains(".bw") && !plain.contains(".bigbed"),
        "{plain}"
    );
    let track = parse_track_output_spec(Some(Path::new("x.zzz")), None)
        .unwrap_err()
        .0;
    assert!(
        track.contains(".bw") && track.contains(".bigbed"),
        "{track}"
    );
    let flag = parse_track_output_spec(None, Some("zzz")).unwrap_err().0;
    assert!(flag.contains(".bb"), "{flag}");
}
