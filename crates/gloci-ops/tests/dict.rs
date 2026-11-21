//! `dict` (D132): the normalisation, tags and errors that the samtools goldens
//! do not pin one by one.
use gloci_io::FastaStream;
use gloci_ops::dict::{DictOptions, dict, normalized_md5, parse_alt_list};
use md5::{Digest, Md5};
use std::num::NonZeroUsize;

fn run(text: &[u8], o: &DictOptions) -> Result<String, String> {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("in.fa");
    std::fs::write(&p, text).unwrap();
    let s = FastaStream::open(p.to_str().unwrap(), NonZeroUsize::new(1).unwrap()).unwrap();
    let mut out = Vec::new();
    dict(s, None, o, &mut out).map_err(|e| e.to_string())?;
    Ok(String::from_utf8(out).unwrap())
}

fn opts() -> DictOptions {
    DictOptions {
        uri: Some("x".into()),
        ..Default::default()
    }
}

fn md5_of(line: &[u8]) -> (String, u64) {
    let mut m = Md5::new();
    let n = normalized_md5(&mut m, line);
    let hex = m.finalize().iter().map(|b| format!("{b:02x}")).collect();
    (hex, n)
}

#[test]
fn dict_md5_normalisation() {
    let (a, na) = md5_of(b"ACGT");
    let (b, nb) = md5_of(b" AC\tgt\r");
    assert_eq!((a.as_str(), na), ("f1f8f4bf413b16ad135722aa4591043e", 4));
    assert_eq!((b.as_str(), nb), (a.as_str(), 4));
    let (c, nc) = md5_of(b"acgtNNRYKM*-");
    assert_eq!(nc, 12);
    assert_eq!(c, md5_of(b"ACGTNNRYKM*-").0);
    let (e, ne) = md5_of(b"");
    assert_eq!((e.as_str(), ne), ("d41d8cd98f00b204e9800998ecf8427e", 0));
    // Controls, DEL and non-ASCII bytes are dropped; `!` and `~` are kept.
    assert_eq!(md5_of(b"\x00!\x7f~\xc3\xa9\n"), md5_of(b"!~"));
    assert_eq!(md5_of(b"!~").1, 2);
}

#[test]
fn dict_tags() {
    let fa = b">chr1 d\nA\n>1\nA\n>chrM\nA\n>M\nA\n>MT\nA\n>chr\nA\n>mt\nA\n";
    let o = DictOptions {
        alias: true,
        assembly: Some("GRCh38".into()),
        species: Some("human".into()),
        alt: Some(parse_alt_list(b"MT\r\n\n1\textra\nchr1 d\n")),
        hd_line: true,
        ..opts()
    };
    let got = run(fa, &o).unwrap();
    let ans: Vec<&str> = got
        .lines()
        .skip(1)
        .map(|l| l.split('\t').find(|f| f.starts_with("AN:")).unwrap())
        .collect();
    assert_eq!(
        ans,
        [
            "AN:1",
            "AN:chr1",
            "AN:M,chrMT,MT",
            "AN:chrM,chrMT,MT",
            "AN:chrMT,chrM,M",
            "AN:",
            "AN:chrmt"
        ]
    );
    assert!(got.starts_with("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:1\tM5:"));
    // Order: M5, AH, AN, UR, AS, SP; AH only for names in the list.
    let l = got.lines().find(|l| l.contains("SN:MT\t")).unwrap();
    let tags: Vec<&str> = l.split('\t').skip(4).collect();
    assert_eq!(
        tags,
        ["AH:*", "AN:chrMT,chrM,M", "UR:x", "AS:GRCh38", "SP:human"]
    );
    assert_eq!(got.lines().filter(|l| l.contains("AH:*")).count(), 2);
}

#[test]
fn dict_default_uri_and_stdin() {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("in.fa");
    std::fs::write(&p, b">a\nAC\n").unwrap();
    let s = FastaStream::open(p.to_str().unwrap(), NonZeroUsize::new(1).unwrap()).unwrap();
    let mut out = Vec::new();
    let o = DictOptions::default();
    dict(s, Some(&p), &o, &mut out).unwrap();
    let want = format!(
        "\tUR:file://{}\n",
        std::fs::canonicalize(&p).unwrap().display()
    );
    assert!(String::from_utf8(out).unwrap().ends_with(&want));
    // No path (standard input) and no --uri: no UR tag.
    assert!(!run(b">a\nAC\n", &o).unwrap().contains("UR:"));
}

#[test]
fn dict_duplicate_names() {
    let e = run(b">a\nA\n>b\nA\n>a x\nC\n", &opts()).unwrap_err();
    assert!(e.contains("line 5"), "{e}");
    assert!(e.contains("line 1"), "{e}");
    assert!(e.contains("'a'"), "{e}");
}

#[test]
fn blank_lines_and_empty_sequences() {
    let got = run(b">a\n\nAC\n\n>b\n", &opts()).unwrap();
    assert!(got.contains("SN:a\tLN:2\t"));
    assert!(got.contains("SN:b\tLN:0\tM5:d41d8cd98f00b204e9800998ecf8427e"));
}

#[test]
fn alt_list_parsing() {
    let s = parse_alt_list(b"a\tb\nc\r\n\nd e\n");
    let mut v: Vec<_> = s.into_iter().collect();
    v.sort();
    assert_eq!(v, [b"a".to_vec(), b"c".to_vec(), b"d e".to_vec()]);
}
