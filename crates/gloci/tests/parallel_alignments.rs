//! Per-chromosome runs on indexed BAM inputs (D101, D107, D124): the output equals the streamed
//! output at every thread count, errors name one place, and an indexed B is not drained.
mod common;
use common::*;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn bytes(args: &[String]) -> (Option<i32>, Vec<u8>, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code(),
        o.stdout,
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

/// SAM text over `contigs` (name, length) with `reads` per contig, plus unmapped reads.
fn sam_text(contigs: &[(String, u32)], reads: &[(usize, u32, u16, &str)]) -> String {
    let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, l) in contigs {
        t.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    for (i, (c, pos, flag, cigar)) in reads.iter().enumerate() {
        let len: usize = cigar.trim_end_matches('M').parse().unwrap_or(10);
        t.push_str(&format!(
            "r{i}\t{flag}\t{}\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t{}\n",
            contigs[*c].0,
            "A".repeat(len),
            "I".repeat(len)
        ));
    }
    t
}

fn bam_of(dir: &Path, name: &str, sam: &str) -> PathBuf {
    let sp = dir.join(format!("{name}.sam"));
    std::fs::write(&sp, sam).unwrap();
    let out = dir.join(name);
    gloci()
        .args(["view", &s(&sp), "-O", "bam", "-o", &s(&out)])
        .assert()
        .success();
    out
}

fn index(p: &Path) {
    gloci().args(["index", &s(p)]).assert().success();
}

fn small_contigs() -> Vec<(String, u32)> {
    ["chr1", "chr2", "chr3", "chr10", "chrX"]
        .iter()
        .map(|n| (n.to_string(), 100_000))
        .collect()
}

/// Reads on every contig, including flagged ones (unmapped, secondary, duplicate), in order.
fn small_sam() -> String {
    let mut reads = Vec::new();
    for c in 0..5usize {
        for k in 0..40u32 {
            let flag = match k % 9 {
                3 => 16,
                5 => 4,
                7 => 1024,
                _ => 0,
            };
            reads.push((c, 1 + k * 37, flag, "30M"));
        }
    }
    sam_text(&small_contigs(), &reads)
}

const BED: &str =
    "chr1\t0\t300\nchr1\t900\t1100\nchr2\t100\t900\nchr3\t0\t50\nchr10\t0\t1500\nchrX\t500\t520\n";

/// Output with the `@PG` command line removed (it may name `-t`).
fn no_cl(out: &[u8]) -> Vec<u8> {
    String::from_utf8_lossy(out)
        .lines()
        .map(|l| {
            if l.starts_with("@PG") {
                l.split('\t')
                    .filter(|f| !f.starts_with("CL:"))
                    .collect::<Vec<_>>()
                    .join("\t")
            } else {
                l.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        .into_bytes()
}

#[test]
fn bam_per_chromosome_equals_streaming() {
    let dir = tempdir();
    let bam = bam_of(dir.path(), "a.bam", &small_sam());
    index(&bam);
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let bam = s(&bam);
    let mut cases: Vec<Vec<String>> = Vec::new();
    let mk = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    for region in [None, Some("chr2:100-500")] {
        let mut extra: Vec<&str> = Vec::new();
        if let Some(r) = region {
            extra.extend(["--region", r]);
        }
        for (a, b) in [(&bam, &bed), (&bed, &bam), (&bam, &bam)] {
            for mode in ["any", "count"] {
                let mut v = mk(&["intersect", "-a", a, "-b", b, "--mode", mode]);
                if mode == "any" && a == &bam {
                    v.extend(mk(&["-O", "sam"]));
                }
                v.extend(mk(&extra));
                cases.push(v);
            }
            let mut v = mk(&["intersect", "-a", a, "-b", b, "--mode", "none"]);
            if a == &bam {
                v.extend(mk(&["-O", "sam"]));
            }
            v.extend(mk(&extra));
            cases.push(v);
            for cov in [&[][..], &["--hist"]] {
                let mut v = mk(&["coverage", "-a", a, "-b", b]);
                v.extend(mk(cov));
                v.extend(mk(&extra));
                cases.push(v);
            }
        }
    }
    for base in cases {
        let at = |t: &str| {
            let mut v = base.clone();
            v.extend(mk(&["-t", t]));
            let (code, out, err) = bytes(&v);
            (code, no_cl(&out), err)
        };
        let one = at("1");
        assert_eq!(one.0, Some(0), "{base:?}: {}", one.2);
        assert!(
            !one.1.is_empty() || base.contains(&"none".to_string()),
            "{base:?}"
        );
        for t in ["2", "8"] {
            assert_eq!(one, at(t), "{base:?} -t {t}");
        }
    }
}

#[test]
fn bam_per_chromosome_verbose_counts_and_warnings_once() {
    let dir = tempdir();
    // Reads with no place, at the end of the file, are counted from the index.
    let sam = small_sam()
        + "u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\nu2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
    let bam = bam_of(dir.path(), "a.bam", &sam);
    index(&bam);
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let run = |t: &str| {
        let args: Vec<String> = [
            "intersect",
            "-a",
            &s(&bam),
            "-b",
            &bed,
            "-O",
            "sam",
            "-v",
            "-t",
            t,
        ]
        .map(String::from)
        .to_vec();
        let (code, out, err) = bytes(&args);
        assert_eq!(code, Some(0), "{err}");
        let stats: Vec<String> = err
            .lines()
            .filter(|l| l.contains("unmapped") || l.contains("skipped"))
            .map(String::from)
            .collect();
        (no_cl(&out), stats)
    };
    let one = run("1");
    assert!(!one.1.is_empty(), "expected a skipped-reads line at -v");
    assert_eq!(one, run("2"));
    assert_eq!(one, run("8"));
}

/// An indexed BAM B is not read to its end when A ends early (the BED behaviour): a block
/// damaged far past where a tiny A ends is never read.
#[test]
fn indexed_bam_b_is_not_drained() {
    let dir = tempdir();
    let contigs = small_contigs();
    let mut reads = vec![(0usize, 100u32, 0u16, "30M")];
    for k in 0..30_000u32 {
        reads.push((4, 1 + k * 3, 0, "30M"));
    }
    let good = bam_of(dir.path(), "big.bam", &sam_text(&contigs, &reads));
    index(&good);
    let mut raw = std::fs::read(&good).unwrap();
    // Flip bytes inside a data block near the end, before the 28-byte EOF block.
    let at = raw.len() - 28 - 200;
    for b in &mut raw[at..at + 8] {
        *b ^= 0xff;
    }
    std::fs::write(&good, &raw).unwrap();
    let a = s(&plain_file(dir.path(), "a.bed", "chr1\t0\t200\n"));
    for t in ["1", "2", "8"] {
        let (code, _, err) =
            bytes(&["intersect", "-a", &a, "-b", &s(&good), "-t", t].map(String::from));
        assert_eq!(code, Some(0), "-t {t}: {err}");
    }
    // While A reaching the damaged chromosome reports it, alike at every -t.
    let a = s(&plain_file(
        dir.path(),
        "a2.bed",
        "chr1\t0\t200\nchrX\t0\t10\nchrX\t80000\t80100\n",
    ));
    let at = |t: &str| bytes(&["intersect", "-a", &a, "-b", &s(&good), "-t", t].map(String::from));
    // (Output already written before the error differs: only the code and message must not.)
    // (The stale-index warning is for runs that use the index: it is not an error line.)
    let errors = |e: &str| {
        e.lines()
            .filter(|l| l.starts_with("error"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let one = at("1");
    assert_eq!(one.0, Some(1), "{}", one.2);
    for t in ["2", "8"] {
        let r = at(t);
        assert_eq!((one.0, errors(&one.2)), (r.0, errors(&r.2)), "-t {t}");
    }
}

fn crc32(data: &[u8]) -> u32 {
    let mut c = !0u32;
    for &b in data {
        c ^= u32::from(b);
        for _ in 0..8 {
            c = (c >> 1) ^ (0xEDB8_8320 & (c & 1).wrapping_neg());
        }
    }
    !c
}

/// `payload` as BGZF made of stored (uncompressed) blocks: a block's compressed size depends
/// only on the payload's length, so a copy with one byte changed has the same block offsets,
/// and so the same index.
fn stored_bgzf(payload: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut block = |data: &[u8]| {
        let bsize = 18 + 5 + data.len() + 8 - 1;
        out.extend_from_slice(&[
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
        ]);
        out.extend_from_slice(&(bsize as u16).to_le_bytes());
        out.push(1);
        out.extend_from_slice(&(data.len() as u16).to_le_bytes());
        out.extend_from_slice(&(!(data.len() as u16)).to_le_bytes());
        out.extend_from_slice(data);
        out.extend_from_slice(&crc32(data).to_le_bytes());
        out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    };
    for c in payload.chunks(50_000) {
        block(c);
    }
    // The standard empty block that ends a BGZF file.
    out.extend_from_slice(&[
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0,
        0, 0, 0, 0,
    ]);
    out
}

/// A corrupt record on an A chromosome of an indexed B gives one message at every `-t`,
/// naming the record's place by chromosome and by file.
#[test]
fn bam_error_same_at_every_t() {
    use std::io::Read;
    let dir = tempdir();
    let mut reads = Vec::new();
    for c in 0..5usize {
        for k in 0..10u32 {
            reads.push((c, 1 + k * 100, 0u16, "30M"));
        }
    }
    let src = bam_of(dir.path(), "src.bam", &sam_text(&small_contigs(), &reads));
    let mut payload = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(&src).unwrap())
        .read_to_end(&mut payload)
        .unwrap();
    let good = dir.path().join("good.bam");
    std::fs::write(&good, stored_bgzf(&payload)).unwrap();
    index(&good);
    // Walk to the fourth record of chr2 (reference id 1) and give its CIGAR an operation
    // code BAM does not define.
    let u32_at =
        |b: &[u8], at: usize| u32::from_le_bytes(b[at..at + 4].try_into().unwrap()) as usize;
    let l_text = u32_at(&payload, 4);
    let mut at = 8 + l_text;
    let n_ref = u32_at(&payload, at);
    at += 4;
    for _ in 0..n_ref {
        at += 4 + u32_at(&payload, at) + 4;
    }
    let (mut on_chr2, mut hit) = (0, 0);
    while at < payload.len() {
        let size = u32_at(&payload, at);
        if u32_at(&payload, at + 4) == 1 {
            on_chr2 += 1;
            if on_chr2 == 4 {
                let l_name = usize::from(payload[at + 4 + 8]);
                payload[at + 4 + 32 + l_name] |= 0xF;
                hit = at;
                break;
            }
        }
        at += 4 + size;
    }
    assert!(hit > 0);
    let bad = dir.path().join("bad.bam");
    std::fs::write(&bad, stored_bgzf(&payload)).unwrap();
    let index_file = dir.path().join("bad.bam.bai");
    std::fs::copy(dir.path().join("good.bam.bai"), &index_file).unwrap();
    // The index is as new as the file.
    std::fs::File::options()
        .append(true)
        .open(&index_file)
        .unwrap()
        .set_modified(std::time::SystemTime::now())
        .unwrap();
    let a = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t0\t900\nchr2\t0\t900\nchr3\t0\t300\n",
    ));
    let at = |t: &str| {
        let (code, _, err) =
            bytes(&["intersect", "-a", &a, "-b", &s(&bad), "-t", t].map(String::from));
        let errors: Vec<String> = err
            .lines()
            .filter(|l| l.starts_with("error"))
            .map(String::from)
            .collect();
        (code, errors)
    };
    let one = at("1");
    assert_eq!(one.0, Some(1), "{one:?}");
    let msg = one.1.join("\n");
    assert!(
        msg.contains("record 4 of chr2") && msg.contains("record 14 of the file"),
        "{msg}"
    );
    for t in ["2", "8"] {
        assert_eq!(one, at(t), "-t {t}");
    }
}

#[test]
fn huge_header_per_chromosome() {
    let dir = tempdir();
    let n = 200_000usize;
    let contigs: Vec<(String, u32)> = (0..n).map(|i| (format!("c{i:06}"), 10_000)).collect();
    let mut reads = Vec::new();
    let mut bed = String::new();
    for k in 0..20_000usize {
        let c = k * 10;
        for j in 0..2u32 {
            reads.push((c, 10 + j * 100, 0u16, "20M"));
        }
        bed.push_str(&format!("{}\t0\t60\n", contigs[c].0));
    }
    let bam = bam_of(dir.path(), "huge.bam", &sam_text(&contigs, &reads));
    index(&bam);
    let bed = s(&plain_file(dir.path(), "t.bed", &bed));
    let run = |t: &str| {
        let t0 = Instant::now();
        let r = bytes(
            &[
                "intersect",
                "-a",
                &bed,
                "-b",
                &s(&bam),
                "--mode",
                "count",
                "-t",
                t,
            ]
            .map(String::from),
        );
        (r, t0.elapsed())
    };
    let (one, _) = run("1");
    assert_eq!(one.0, Some(0), "{}", one.2);
    let (eight, took) = run("8");
    assert_eq!(one, eight);
    assert!(
        took.as_secs() < 60,
        "-t 8 took {took:?}: per-chromosome work must not be quadratic in contigs"
    );
    let _ = (one, eight);
}

/// The one-time warnings of D101 (a read cut at its reference's end, a read past it) are
/// printed once per input, naming the first such read in file order, at every `-t`.
#[test]
fn bam_one_time_warnings_once_at_every_t() {
    let dir = tempdir();
    let mut reads = Vec::new();
    for c in 0..5usize {
        reads.push((c, 100, 0u16, "30M"));
        reads.push((c, 99_990, 0, "30M"));
        reads.push((c, 100_001, 0, "30M"));
    }
    let bam = bam_of(dir.path(), "edge.bam", &sam_text(&small_contigs(), &reads));
    index(&bam);
    let bed = s(&plain_file(
        dir.path(),
        "t.bed",
        "chr1\t0\t200\nchr2\t0\t200\nchr3\t0\t200\nchr10\t0\t200\nchrX\t0\t200\n",
    ));
    let warnings = |t: &str| {
        let (code, out, err) = bytes(
            &[
                "intersect",
                "-a",
                &s(&bam),
                "-b",
                &bed,
                "-O",
                "sam",
                "-t",
                t,
            ]
            .map(String::from),
        );
        assert_eq!(code, Some(0), "{err}");
        let w: Vec<String> = err
            .lines()
            .filter(|l| l.starts_with("warning"))
            .map(String::from)
            .collect();
        (no_cl(&out), w)
    };
    let one = warnings("1");
    assert_eq!(one.1.len(), 2, "{:?}", one.1);
    for t in ["2", "8"] {
        let r = warnings(t);
        assert_eq!(one.1, r.1, "-t {t}");
        assert_eq!(
            String::from_utf8_lossy(&one.0),
            String::from_utf8_lossy(&r.0),
            "-t {t}"
        );
    }
}
