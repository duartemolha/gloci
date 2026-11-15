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

fn usage(e: IoError) -> String {
    match e {
        IoError::Usage(u) => u.0,
        e => panic!("not a usage error: {e}"),
    }
}

/// The inputs a command line gives `-` are given to the run's token: the second is refused by
/// both names, and another run's token never sees them.
#[test]
fn standard_input_is_given_to_one_input_of_a_run() {
    use gloci_io::cancel::RunToken;
    let run = RunToken::new();
    run.give_stdin("-a").unwrap();
    assert_eq!(
        usage(run.give_stdin("-b").unwrap_err()),
        "standard input is given for both -a and -b; give one of them as a file"
    );
    let next = RunToken::new();
    next.give_stdin("-b").unwrap();
    assert_eq!(
        usage(next.give_stdin("-b (value 2)").unwrap_err()),
        "standard input is given for both -b and -b (value 2); give one of them as a file"
    );
}

/// The reader that opens `-` second in a run is refused before it reads, through every opener;
/// the first reader takes what the command line gave. Each test claims standard input before it
/// opens `-`, so that nothing reads the test's own standard input.
#[test]
fn a_second_reader_of_standard_input_is_refused() {
    use gloci_io::cancel::{RunGuard, RunToken, claim_stdin};
    let two = "standard input is read for two inputs; give one of them as a file";
    let run = RunToken::new();
    {
        let _g = RunGuard::enter(&run);
        claim_stdin().unwrap();
        assert_eq!(usage(open_input("-", &opts(1)).unwrap_err()), two);
        assert_eq!(usage(open_raw("-").err().unwrap()), two);
        assert_eq!(usage(open_stdin().unwrap_err()), two);
        let dir = tempfile::tempdir().unwrap();
        let mut spool = InputSpool::create(dir.path()).unwrap();
        assert_eq!(usage(spool.spool("-").unwrap_err()), two);
        assert_eq!(usage(read_new_header("-").unwrap_err()), two);
    }
    // A worker of the run is the run.
    let watch = run.watch();
    std::thread::spawn(move || {
        let _g = watch.enter();
        assert_eq!(usage(claim_stdin().unwrap_err()), two);
    })
    .join()
    .unwrap();
    // The reader of what the command line gave takes it; one more is refused by that name.
    let given = RunToken::new();
    given.give_stdin("-a").unwrap();
    let _g = RunGuard::enter(&given);
    claim_stdin().unwrap();
    assert_eq!(
        usage(claim_stdin().unwrap_err()),
        "standard input is given for -a and read for another input too; give one of them as a \
file"
    );
}

/// Outside any run (no token) there is nothing to scope the claim to, and nothing is checked.
#[test]
fn standard_input_is_not_checked_outside_a_run() {
    assert!(gloci_io::cancel::current().is_none());
    gloci_io::cancel::claim_stdin().unwrap();
    gloci_io::cancel::claim_stdin().unwrap();
}
