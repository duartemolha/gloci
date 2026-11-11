//! FASTA indexes (`.fai`/`.gzi`, byte-equal to samtools 1.21) and random access (D108).
use gloci_core::Interval;
use gloci_io::*;
use gloci_testkit::fixtures;
use std::collections::BTreeMap;
use std::io::{Cursor, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

fn inputs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/fa")
}

/// Copy a committed fixture (without its index) into `dir`.
fn fixture(dir: &Path, name: &str) -> PathBuf {
    let dest = dir.join(name);
    std::fs::copy(inputs().join(name), &dest).unwrap();
    dest
}

fn t(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

/// A captured copy of every `warning:` the library logs.
static WARNINGS: Mutex<Vec<String>> = Mutex::new(Vec::new());

struct Capture;
impl log::Log for Capture {
    fn enabled(&self, _: &log::Metadata) -> bool {
        true
    }
    fn log(&self, r: &log::Record) {
        WARNINGS.lock().unwrap().push(r.args().to_string());
    }
    fn flush(&self) {}
}

fn capture_warnings() {
    let _ = log::set_logger(&Capture);
    log::set_max_level(log::LevelFilter::Warn);
}

#[test]
fn fai_matches_samtools() {
    let dir = tempfile::tempdir().unwrap();
    for name in ["ref.fa", "crlf.fa", "desc.fa", "empty-seq.fa", "ref.fa.gz"] {
        let p = fixture(dir.path(), name);
        for threads in [1, 2, 8] {
            let written = write_fasta_index(&p, t(threads)).unwrap();
            let fai = std::fs::read(fai_path(&p)).unwrap();
            assert_eq!(
                fai,
                std::fs::read(inputs().join(format!("{name}.fai"))).unwrap(),
                "{name} -t {threads}: .fai differs from samtools"
            );
            if name.ends_with(".gz") {
                assert_eq!(written, vec![fai_path(&p), gzi_path(&p)]);
                assert_eq!(
                    std::fs::read(gzi_path(&p)).unwrap(),
                    std::fs::read(inputs().join(format!("{name}.gzi"))).unwrap(),
                    "{name} -t {threads}: .gzi differs from htslib"
                );
            } else {
                assert_eq!(written, vec![fai_path(&p)]);
                assert!(!gzi_path(&p).exists());
            }
        }
    }
    // The fixtures differ from each other, so each comparison discriminates.
    let a = std::fs::read(inputs().join("ref.fa.fai")).unwrap();
    let b = std::fs::read(inputs().join("crlf.fa.fai")).unwrap();
    assert_ne!(a, b);
    // chr2 of ref.fa.gz spans several BGZF blocks: the .gzi has real entries.
    let gzi = std::fs::read(inputs().join("ref.fa.gz.gzi")).unwrap();
    assert!(u64::from_le_bytes(gzi[..8].try_into().unwrap()) >= 3);
}

#[test]
fn fai_edge_cases_match_samtools_rules() {
    let idx = |s: &str| build_fai(&mut Cursor::new(s.as_bytes().to_vec()), "x.fa");
    // Values taken from `samtools faidx` 1.21 on the same bytes.
    let row = |s: &str| -> String { String::from_utf8(fai_bytes(&idx(s).unwrap())).unwrap() };
    assert_eq!(row(">a\nACGT"), "a\t4\t3\t4\t5\n"); // no final newline: width counts one
    assert_eq!(row(">a\r\nACGT\r"), "a\t4\t4\t4\t6\n");
    assert_eq!(
        row(">a\nACGT\nAC\n\n\n>b\nAA\n"),
        "a\t6\t3\t4\t5\nb\t2\t16\t2\t3\n"
    );
    assert_eq!(
        row(">a\r\nACGT\r\nAC\r\n>b\nAA\n"),
        "a\t6\t4\t4\t6\nb\t2\t17\t2\t3\n"
    );
    assert_eq!(row("\n\n>a\nAC\n"), "a\t2\t5\t2\t3\n");
    assert_eq!(row(">a\n>b\nACGT\n"), "b\t4\t6\t4\t5\n"); // empty `a` is left out
    assert_eq!(row(">a b\tc\nAC\n"), "a\t2\t7\t2\t3\n");
    assert_eq!(row(">a\nA\nA\nA\n"), "a\t3\t3\t1\t2\n");
    // Rejected by samtools too.
    for bad in [
        ">a\nACGT\nACGTA\n",      // last line longer
        ">a\nACG\nACGT\n",        // short line in the middle
        ">a\nACGT\n\nAC\n",       // blank line inside
        ">a\n\nACGT\n",           // blank line after the header
        "ACGT\n>a\nAC\n",         // no header
        ">a\nACGT\n>b\n",         // last entry empty
        ">a\nACGT\n>b",           // header only at the end
        "",                       // nothing
        ">a\nACGT\r\nACGT\nAC\n", // mixed line ends in one sequence
    ] {
        assert!(idx(bad).is_err(), "{bad:?}");
    }
}

#[test]
fn fai_errors() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();

    // A uneven line: exit 1 (a data error), naming the line.
    let p = fixture(d, "irregular.fa");
    match write_fasta_index(&p, t(1)).unwrap_err() {
        IoError::Parse { loc, msg, path } => {
            assert_eq!(loc, Loc::Line(4));
            assert!(path.ends_with("irregular.fa"), "{path}");
            assert!(msg.contains("chr1"), "{msg}");
        }
        e => panic!("{e:?}"),
    }
    assert!(!fai_path(&p).exists(), "no index is published on error");

    // A duplicate name: exit 1, naming the line and the name.
    let dup = fixtures::write(d, "dup.fa", b">chr1\nACGT\n>chr2\nAC\n>chr1\nGG\n");
    match write_fasta_index(&dup, t(1)).unwrap_err() {
        IoError::Parse { loc, msg, .. } => {
            assert_eq!(loc, Loc::Line(5));
            assert!(
                msg.contains("chr1") && msg.contains("more than once"),
                "{msg}"
            );
        }
        e => panic!("{e:?}"),
    }
    // An empty name and a space inside a sequence line.
    let bad = fixtures::write(d, "noname.fa", b">\nACGT\n");
    assert!(matches!(
        write_fasta_index(&bad, t(1)).unwrap_err(),
        IoError::Parse { .. }
    ));
    let sp = fixtures::write(d, "space.fa", b">a\nAC GT\n");
    assert!(matches!(
        write_fasta_index(&sp, t(1)).unwrap_err(),
        IoError::Parse {
            loc: Loc::Line(2),
            ..
        }
    ));

    // No .fai: usage error (exit 2) with the exact command.
    let r = fixture(d, "ref.fa");
    match IndexedFasta::open(&r).unwrap_err() {
        IoError::Usage(u) => assert!(
            u.0.contains(&format!("gloci index {}", r.display())),
            "{}",
            u.0
        ),
        e => panic!("{e:?}"),
    }
    // BGZF without .gzi.
    let gz = fixture(d, "ref.fa.gz");
    std::fs::copy(inputs().join("ref.fa.gz.fai"), fai_path(&gz)).unwrap();
    match IndexedFasta::open(&gz).unwrap_err() {
        IoError::Usage(u) => assert!(
            u.0.contains(&format!("gloci index {}", gz.display())),
            "{}",
            u.0
        ),
        e => panic!("{e:?}"),
    }
    // Plain gzip: usage error with a bgzip recipe, from `open` and from `index`.
    let plain = std::fs::read(inputs().join("ref.fa")).unwrap();
    let pgz = fixtures::write(d, "plain.fa.gz", &fixtures::gzip(&plain));
    for e in [
        IndexedFasta::open(&pgz).unwrap_err(),
        write_fasta_index(&pgz, t(1)).unwrap_err(),
    ] {
        match e {
            IoError::Usage(u) => assert!(u.0.contains("bgzip"), "{}", u.0),
            e => panic!("{e:?}"),
        }
    }

    // An edited FASTA no longer matches its index.
    let ed = fixture(d, "ref.fa");
    write_fasta_index(&ed, t(1)).unwrap();
    let mut text = std::fs::read(&ed).unwrap();
    // Insert one base in the first line of chr1: every later line break moves.
    let off = std::str::from_utf8(&std::fs::read(inputs().join("ref.fa.fai")).unwrap())
        .unwrap()
        .lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .parse::<usize>()
        .unwrap();
    text.insert(off + 5, b'A');
    std::fs::write(&ed, &text).unwrap();
    let mut f = IndexedFasta::open(&ed).unwrap();
    let e = f
        .fetch(b"chr1", iv(0, 500), &mut |_| Ok(()))
        .unwrap_err()
        .to_string();
    assert!(e.contains("does not match"), "{e}");
    assert!(e.contains(&format!("gloci index {}", ed.display())), "{e}");

    // A fetch outside the sequence and an unknown name are errors, not panics.
    let ok = fixture(d, "desc.fa");
    write_fasta_index(&ok, t(1)).unwrap();
    let mut f = IndexedFasta::open(&ok).unwrap();
    assert!(f.fetch(b"chrM", iv(0, 151), &mut |_| Ok(())).is_err());
    assert!(f.fetch(b"nope", iv(0, 1), &mut |_| Ok(())).is_err());
}

#[test]
fn stale_index_warns() {
    capture_warnings();
    let dir = tempfile::tempdir().unwrap();
    let p = fixture(dir.path(), "ref.fa");
    write_fasta_index(&p, t(1)).unwrap();
    let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
    std::fs::File::options()
        .write(true)
        .open(fai_path(&p))
        .unwrap()
        .set_modified(old)
        .unwrap();
    WARNINGS.lock().unwrap().clear();
    IndexedFasta::open(&p).unwrap();
    let w = WARNINGS.lock().unwrap().join("\n");
    assert!(
        w.contains("is older than") && w.contains("gloci index"),
        "{w}"
    );
    // A current index does not warn.
    write_fasta_index(&p, t(1)).unwrap();
    WARNINGS.lock().unwrap().clear();
    IndexedFasta::open(&p).unwrap();
    assert!(WARNINGS.lock().unwrap().is_empty());
}

/// Naive reference: name -> bases, parsed without any index.
fn parse_naive(text: &[u8]) -> BTreeMap<Vec<u8>, Vec<u8>> {
    let mut m = BTreeMap::new();
    let mut cur: Option<Vec<u8>> = None;
    for line in text.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if let Some(h) = line.strip_prefix(b">") {
            let n: Vec<u8> = h
                .iter()
                .copied()
                .take_while(|b| !b" \t".contains(b))
                .collect();
            m.insert(n.clone(), Vec::new());
            cur = Some(n);
        } else if let Some(n) = &cur {
            m.get_mut(n).unwrap().extend_from_slice(line);
        }
    }
    m
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn collect(
    f: &mut IndexedFasta,
    chrom: &[u8],
    i: Interval,
    rev: bool,
) -> Result<Vec<Vec<u8>>, IoError> {
    let mut chunks = Vec::new();
    let mut cb = |c: &[u8]| {
        chunks.push(c.to_vec());
        Ok(())
    };
    if rev {
        f.fetch_rev(chrom, i, &mut cb)?;
    } else {
        f.fetch(chrom, i, &mut cb)?;
    }
    Ok(chunks)
}

fn check_against_naive(path: &Path, text: &[u8], seed: u64, n: usize) {
    write_fasta_index(path, t(2)).unwrap();
    let naive = parse_naive(text);
    let mut f = IndexedFasta::open(path).unwrap();
    let entries = f.entries().to_vec();
    let gzi = gzi_path(path);
    let blocks: Vec<u64> = if gzi.exists() {
        let b = std::fs::read(&gzi).unwrap();
        b[8..]
            .chunks(16)
            .map(|e| u64::from_le_bytes(e[8..16].try_into().unwrap()))
            .collect()
    } else {
        Vec::new()
    };
    let mut rng = Lcg(seed);
    for _ in 0..n {
        let e = &entries[rng.below(entries.len() as u64) as usize];
        // Interesting edges: line boundaries, BGZF block starts, the ends.
        let edge = |rng: &mut Lcg| -> u64 {
            let p = match rng.below(5) {
                0 => rng.below(e.len + 1),
                1 => (rng.below(e.len / e.line_bases + 1) * e.line_bases)
                    .saturating_add(rng.below(3)),
                2 if !blocks.is_empty() => {
                    let u = blocks[rng.below(blocks.len() as u64) as usize];
                    let rel = u.saturating_sub(e.offset);
                    (rel / e.line_width * e.line_bases + (rel % e.line_width).min(e.line_bases))
                        .saturating_add(rng.below(3))
                        .saturating_sub(1)
                }
                3 => e.len - rng.below(3).min(e.len),
                _ => rng.below(3).min(e.len),
            };
            p.min(e.len)
        };
        let (a, b) = (edge(&mut rng), edge(&mut rng));
        let (s, en) = if a <= b { (a, b) } else { (b, a) };
        let want = &naive[&e.name][s as usize..en as usize];
        let fwd = collect(&mut f, &e.name, iv(s, en), false).unwrap();
        assert_eq!(fwd.concat(), want, "{:?}:{s}-{en}", e.name);
        assert!(
            fwd.iter()
                .all(|c| c.len() as u64 <= FETCH_CHUNK && !c.is_empty())
        );
        let mut rev = collect(&mut f, &e.name, iv(s, en), true).unwrap();
        rev.reverse();
        assert_eq!(rev, fwd);
        if s == en {
            assert!(fwd.is_empty());
        }
    }
    // Whole sequences.
    for e in &entries {
        let got = collect(&mut f, &e.name, iv(0, e.len), false)
            .unwrap()
            .concat();
        assert_eq!(got, naive[&e.name]);
    }
}

#[test]
fn fetch_equals_naive() {
    let dir = tempfile::tempdir().unwrap();
    for (name, seed) in [
        ("ref.fa", 1),
        ("crlf.fa", 2),
        ("desc.fa", 3),
        ("ref.fa.gz", 4),
    ] {
        let p = fixture(dir.path(), name);
        let raw = std::fs::read(&p).unwrap();
        let text = if name.ends_with(".gz") {
            use std::io::Read;
            let mut v = Vec::new();
            noodles_bgzf::io::Reader::new(Cursor::new(raw))
                .read_to_end(&mut v)
                .unwrap();
            v
        } else {
            raw
        };
        check_against_naive(&p, &text, seed, 1000);
    }
}

/// A FASTA with one long sequence, wrapped at `width`.
fn long_fasta(len: usize, width: usize, crlf: bool) -> Vec<u8> {
    let mut rng = Lcg(77);
    let nl: &[u8] = if crlf { b"\r\n" } else { b"\n" };
    let mut out = b">big".to_vec();
    out.extend_from_slice(nl);
    let seq: Vec<u8> = (0..len)
        .map(|_| b"ACGTacgtN"[rng.below(9) as usize])
        .collect();
    for l in seq.chunks(width) {
        out.extend_from_slice(l);
        out.extend_from_slice(nl);
    }
    out
}

#[test]
fn big_sequences_come_in_bounded_chunks() {
    let dir = tempfile::tempdir().unwrap();
    // 3 MiB + 1000 bases: wrapped at 60, unwrapped (one 3 MiB line), CRLF, and BGZF.
    let len = 3 * (1 << 20) + 1000;
    for (name, text, bgzf) in [
        ("w60.fa", long_fasta(len, 60, false), false),
        ("flat.fa", long_fasta(len, len, false), false),
        ("crlf.fa", long_fasta(len, 70, true), false),
        ("w60.fa.gz", long_fasta(len, 60, false), true),
    ] {
        let p = if bgzf {
            fixtures::write(dir.path(), name, &fixtures::bgzip(&text))
        } else {
            fixtures::write(dir.path(), name, &text)
        };
        write_fasta_index(&p, t(2)).unwrap();
        let naive = parse_naive(&text);
        let want = &naive[&b"big"[..]];
        let mut f = IndexedFasta::open(&p).unwrap();
        let chunks = collect(&mut f, b"big", iv(0, len as u64), false).unwrap();
        assert!(chunks.len() >= 3, "{name}: {} chunks", chunks.len());
        assert!(
            chunks.iter().all(|c| c.len() as u64 <= FETCH_CHUNK),
            "{name}: a chunk is over 1 MiB"
        );
        assert_eq!(chunks.concat(), *want, "{name}");
        let mut rev = collect(&mut f, b"big", iv(0, len as u64), true).unwrap();
        rev.reverse();
        assert_eq!(rev, chunks, "{name}");
        // An interval inside the middle chunks too.
        let (s, e) = (500_000u64, 2_700_000u64);
        let mid = collect(&mut f, b"big", iv(s, e), false).unwrap().concat();
        assert_eq!(mid, want[s as usize..e as usize], "{name}");
    }
}

#[test]
fn reopen_shares_the_index_and_reads_independently() {
    let dir = tempfile::tempdir().unwrap();
    let p = fixture(dir.path(), "ref.fa.gz");
    write_fasta_index(&p, t(1)).unwrap();
    let mut a = IndexedFasta::open(&p).unwrap();
    let mut b = a.reopen().unwrap();
    assert_eq!(a.entries(), b.entries());
    assert_eq!(a.len_of(b"chr2"), Some(300_000));
    assert_eq!(a.len_of(b"chrZ"), None);
    let x = collect(&mut a, b"chr2", iv(250_000, 250_100), false).unwrap();
    let y = collect(&mut b, b"chr2", iv(10, 20), false).unwrap();
    let x2 = collect(&mut a, b"chr2", iv(250_000, 250_100), false).unwrap();
    assert_eq!(x, x2);
    assert_eq!(y[0].len(), 10);
}

#[test]
fn stream_reads_lines_headers_and_line_ends() {
    let dir = tempfile::tempdir().unwrap();
    let text = b">a desc\r\nACGT\r\nAC\r\n\r\n>b\nGG\nT";
    for (name, bytes) in [
        ("s.fa", text.to_vec()),
        ("s.fa.gz", fixtures::bgzip(text)),
        ("s.gz", fixtures::gzip(text)),
    ] {
        let p = fixtures::write(dir.path(), name, &bytes);
        for threads in [1, 2, 8] {
            let mut s = FastaStream::open(p.to_str().unwrap(), t(threads)).unwrap();
            let mut got = Vec::new();
            while let Some(l) = s.next().unwrap() {
                let (kind, data) = match l {
                    FastaLine::Header(h) => ('>', h.to_vec()),
                    FastaLine::Seq(q) => ('s', q.to_vec()),
                };
                got.push((kind, data, s.line(), s.crlf()));
            }
            let want: Vec<(char, Vec<u8>, u64, bool)> = vec![
                ('>', b"a desc".to_vec(), 1, true),
                ('s', b"ACGT".to_vec(), 2, true),
                ('s', b"AC".to_vec(), 3, true),
                ('s', b"".to_vec(), 4, true),
                ('>', b"b".to_vec(), 5, false),
                ('s', b"GG".to_vec(), 6, false),
                ('s', b"T".to_vec(), 7, false),
            ];
            assert_eq!(got, want, "{name} -t {threads}");
        }
    }
    // Not FASTA.
    let bed = fixtures::write(dir.path(), "x.bed", b"chr1\t1\t2\n");
    assert!(FastaStream::open(bed.to_str().unwrap(), t(1)).is_err());
}

#[test]
fn stream_splits_a_huge_line_into_pieces() {
    let dir = tempfile::tempdir().unwrap();
    let len = 2_500_000;
    let p = fixtures::write(dir.path(), "flat.fa", &long_fasta(len, len, true));
    let mut s = FastaStream::open(p.to_str().unwrap(), t(1)).unwrap();
    assert!(matches!(s.next().unwrap(), Some(FastaLine::Header(h)) if h == b"big"));
    let (mut total, mut pieces) = (0usize, 0);
    while let Some(l) = s.next().unwrap() {
        let FastaLine::Seq(q) = l else { panic!() };
        assert!(q.len() <= 1 << 20 && !q.contains(&b'\r'));
        total += q.len();
        pieces += 1;
    }
    assert_eq!(total, len);
    assert!(pieces >= 3);
}

#[test]
fn genome_from_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let p = fixture(dir.path(), "ref.fa");
    // Without a .fai: usage error naming the fix.
    match load_genome(&p).unwrap_err() {
        IoError::Usage(u) => assert!(u.0.contains("gloci index"), "{}", u.0),
        e => panic!("{e:?}"),
    }
    write_fasta_index(&p, t(1)).unwrap();
    let a = load_genome(&p).unwrap();
    let b = load_genome(&fai_path(&p)).unwrap();
    assert_eq!(a, b);
    let names: Vec<&[u8]> = (0..a.len()).map(|i| a.name(i)).collect();
    assert_eq!(
        names,
        [&b"chrM"[..], b"chr1", b"chr2", b"chr10", b"chrX"],
        "file order, not karyotype or lexicographic"
    );
    assert_eq!(a.length(1), Some(1000));
    // BGZF FASTA too.
    let g = fixture(dir.path(), "ref.fa.gz");
    assert!(matches!(load_genome(&g).unwrap_err(), IoError::Usage(_)));
    write_fasta_index(&g, t(1)).unwrap();
    assert_eq!(load_genome(&g).unwrap().length(2), Some(300_000));
}

fn mutate(rng: &mut Lcg, src: &[u8], max_edits: u64) -> Vec<u8> {
    let mut v = src.to_vec();
    for _ in 0..=rng.below(max_edits) {
        if v.is_empty() {
            v.push(b'>');
            continue;
        }
        let i = rng.below(v.len() as u64) as usize;
        match rng.below(5) {
            0 => v[i] = rng.below(256) as u8,
            1 => v[i] = b"\n\r>\t A"[rng.below(6) as usize],
            2 => {
                v.remove(i);
            }
            3 => v.insert(i, b"\n\r>\t 9"[rng.below(6) as usize]),
            _ => v.truncate(i),
        }
    }
    v
}

#[test]
fn fuzz_fasta() {
    let ref_fa = std::fs::read(inputs().join("ref.fa")).unwrap();
    let crlf = std::fs::read(inputs().join("crlf.fa")).unwrap();
    let ref_fai = std::fs::read(inputs().join("ref.fa.fai")).unwrap();
    let mut rng = Lcg(2024);
    // Mutated FASTA through build_fai: an error or an index, never a panic.
    for i in 0..2000 {
        let src = if i % 2 == 0 { &ref_fa } else { &crlf };
        let m = mutate(&mut rng, &src[..src.len().min(1500)], 6);
        let _ = build_fai(&mut Cursor::new(m), "fuzz.fa");
    }
    // Mutated .fai lines through open and fetch.
    let dir = tempfile::tempdir().unwrap();
    let p = fixture(dir.path(), "ref.fa");
    let fai = fai_path(&p);
    for i in 0..2000u64 {
        let mut m = mutate(&mut rng, &ref_fai, 4);
        if i % 5 == 0 {
            // Huge and zero numbers, the cases that overflow or divide.
            let n = [
                b"0".as_slice(),
                b"18446744073709551615",
                b"1",
                b"99999999999",
            ][rng.below(4) as usize];
            let at = rng.below(m.len() as u64 + 1) as usize;
            m.splice(at..at, n.iter().copied());
        }
        let mut w = std::fs::File::create(&fai).unwrap();
        w.write_all(&m).unwrap();
        drop(w);
        let Ok(mut f) = IndexedFasta::open(&p) else {
            continue;
        };
        let entries = f.entries().to_vec();
        for e in entries.iter().take(4) {
            let end = e.len.min(rng.below(3000));
            let s = rng.below(end + 1);
            let _ = f.fetch(&e.name, iv(s, end), &mut |_| Ok(()));
            let _ = f.fetch_rev(&e.name, iv(0, e.len.min(2000)), &mut |_| Ok(()));
        }
    }
    // The same through a BGZF file with a mutated .gzi.
    let g = fixture(dir.path(), "ref.fa.gz");
    write_fasta_index(&g, t(1)).unwrap();
    let gzi = std::fs::read(gzi_path(&g)).unwrap();
    for _ in 0..500 {
        std::fs::write(gzi_path(&g), mutate(&mut rng, &gzi, 4)).unwrap();
        if let Ok(mut f) = IndexedFasta::open(&g) {
            let _ = f.fetch(b"chr2", iv(100_000, 100_500), &mut |_| Ok(()));
            let _ = f.fetch(b"chrX", iv(0, 700), &mut |_| Ok(()));
        }
    }
}
