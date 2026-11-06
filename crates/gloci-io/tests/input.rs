use gloci_io::*;
use gloci_testkit::fixtures;
use std::io::Read;
use std::num::NonZeroUsize;

fn opts(n: usize) -> InputOptions {
    InputOptions {
        threads: NonZeroUsize::new(n).unwrap(),
        format: None,
    }
}

fn read_all(mut r: Box<dyn std::io::BufRead + Send>) -> Vec<u8> {
    let mut v = Vec::new();
    r.read_to_end(&mut v).unwrap();
    v
}

fn plain(b: &[u8]) -> Vec<u8> {
    b.to_vec()
}

#[test]
fn detects_by_content_and_keeps_peeked_bytes() {
    let dir = tempfile::tempdir().unwrap();
    for (bytes, fmt) in [
        (b"chr1\t1\t2\n".to_vec(), Format::Bed),
        (
            b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec(),
            Format::Vcf,
        ),
        (
            b"##gff-version 3\nchr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=a\n".to_vec(),
            Format::Gff3,
        ),
        (
            b"chr1\tsrc\texon\t1\t9\t.\t+\t.\tgene_id \"g\";\n".to_vec(),
            Format::Gtf,
        ),
        (
            b"track type=bedGraph\nchr1\t0\t5\t1.5\n".to_vec(),
            Format::BedGraph,
        ),
        (Vec::new(), Format::Bed),
    ] {
        for (wrap, comp) in [
            (plain as fn(&[u8]) -> Vec<u8>, Compression::None),
            (fixtures::gzip, Compression::Gzip),
            (fixtures::bgzip, Compression::Bgzf),
        ] {
            let p = fixtures::write(dir.path(), "x.txt", &wrap(&bytes));
            let i = open_input(p.to_str().unwrap(), &opts(4)).unwrap();
            assert_eq!((i.format, i.compression), (fmt, comp));
            assert_eq!(read_all(i.reader), bytes);
        }
    }
}

fn write(dir: &tempfile::TempDir, name: &str, b: &[u8]) -> String {
    fixtures::write(dir.path(), name, b)
        .to_str()
        .unwrap()
        .to_string()
}

#[test]
fn detect_errors_override_and_truncation() {
    let dir = tempfile::tempdir().unwrap();
    let e = open_input(&write(&dir, "y", b"hello world\n"), &opts(1)).unwrap_err();
    assert!(e.to_string().contains("pass --format"));
    let o = InputOptions {
        format: Some(Format::Bed),
        ..opts(1)
    };
    assert_eq!(
        open_input(&write(&dir, "z", b"hello\n"), &o)
            .unwrap()
            .format,
        Format::Bed
    );
    let mut b = fixtures::bgzip(b"chr1\t1\t2\n");
    b.truncate(b.len() - 28);
    assert!(matches!(
        open_input(&write(&dir, "t.gz", &b), &opts(1)),
        Err(IoError::Truncated { .. })
    ));
    assert!(matches!(
        open_input("/nonexistent/zzz", &opts(1)),
        Err(IoError::Open { .. })
    ));
}

#[test]
fn parallel_and_serial_decompression_match() {
    let dir = tempfile::tempdir().unwrap();
    let mut data = Vec::new();
    let mut i = 0u64;
    while data.len() < 5_000_000 {
        data.extend_from_slice(format!("chr1\t{}\t{}\n", i, i + 7).as_bytes());
        i += 3;
    }
    let big = write(&dir, "big.bed.gz", &fixtures::bgzip(&data));
    let a = read_all(open_input(&big, &opts(1)).unwrap().reader);
    let b = read_all(open_input(&big, &opts(8)).unwrap().reader);
    assert_eq!(a, data);
    assert_eq!(a, b);
}

#[test]
fn loc_display_and_format_parse() {
    let e = IoError::Parse {
        path: "a.bed".into(),
        loc: Loc::Line(1234567),
        msg: "bad".into(),
    };
    assert_eq!(e.to_string(), "a.bed line 1,234,567: bad");
    let e = IoError::Parse {
        path: "a.bcf".into(),
        loc: Loc::Record(2),
        msg: "bad".into(),
    };
    assert_eq!(e.to_string(), "a.bcf record 2: bad");
    assert_eq!("gff".parse::<Format>().unwrap(), Format::Gff3);
    assert!("nope".parse::<Format>().is_err());
}
