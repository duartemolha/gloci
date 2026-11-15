//! Helpers shared by the behaviour test modules.

use crate::common::*;
use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::path::PathBuf;

pub(crate) fn bin() -> PathBuf {
    assert_cmd::cargo::cargo_bin("gloci")
}

/// Deterministic unsorted BED with `n` records over 12 chromosomes.
pub(crate) fn bed_text(n: usize) -> String {
    let mut s = String::with_capacity(n * 40);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let start = (x >> 33) % 5_000_000;
        writeln!(
            s,
            "chr{}\t{start}\t{}\tr{i}\t{}\t+",
            1 + (x >> 20) % 12,
            start + 50 + (x >> 50) % 500,
            i % 1000
        )
        .unwrap();
    }
    s
}

pub(crate) fn vcf_text(n: usize, sorted: bool) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    let mut x: u64 = 7;
    for i in 0..n {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
        let (chrom, pos) = if sorted {
            (1 + i * 2 / n, 1 + i * 20)
        } else {
            (1 + (x >> 40) as usize % 2, 1 + (x >> 20) as usize % 100_000)
        };
        writeln!(s, "chr{chrom}\t{pos}\t.\tA\tC\t.\t.\t.").unwrap();
    }
    s
}

/// Run with a thread option in the given spelling (`-t 1`, `-t2`, `--threads=8`).
pub(crate) fn stdout_with(args: &[&str], t: &str) -> Vec<u8> {
    let opt: Vec<String> = match t {
        "1" => vec!["-t".into(), "1".into()],
        "2" => vec!["-t2".into()],
        n => vec![format!("--threads={n}")],
    };
    let out = gloci().args(args).args(&opt).output().unwrap();
    assert_eq!(
        out.status.code(),
        Some(0),
        "{args:?} -t {t}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

pub(crate) fn index_bytes(src: &Path, t: &str) -> Vec<u8> {
    let d = tempdir();
    let f = d.path().join("in.vcf.gz");
    std::fs::copy(src, &f).unwrap();
    gloci()
        .args(["index", f.to_str().unwrap(), "-t", t])
        .assert()
        .success();
    std::fs::read(d.path().join("in.vcf.gz.tbi")).unwrap()
}

pub(crate) fn s_of(p: &Path) -> String {
    p.to_str().unwrap().to_string()
}

pub(crate) fn sv(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

/// A text VCF without the `##glociCommand` line.
pub(crate) fn without_command(bytes: &[u8]) -> Vec<u8> {
    String::from_utf8_lossy(bytes)
        .lines()
        .filter(|l| !l.starts_with("##glociCommand"))
        .flat_map(|l| [l, "\n"])
        .collect::<String>()
        .into_bytes()
}

pub(crate) fn run_stdout(args: &[&str]) -> Vec<u8> {
    let o = gloci().args(args).output().unwrap();
    assert_eq!(o.status.code(), Some(0), "{args:?}: {:?}", o.stderr);
    o.stdout
}

/// Run `args` plus `-t T` for T in 1, 2, 8 (three spellings) and require
/// byte-identical output. This does not exercise chunk boundaries:
/// `GLOCI_TEST_CHUNK_MEM_LIMIT` only shrinks the chunks of the relation
/// commands' chunk driver (intersect, subtract, closest, map, coverage),
/// which the P2 commands do not use.
pub(crate) fn same_at_every_thread_count(what: &str, args: &[&str]) -> Vec<u8> {
    let outs: Vec<Vec<u8>> = ["1", "2", "8"]
        .iter()
        .map(|t| stdout_with(args, t))
        .collect();
    assert_eq!(outs[0], outs[1], "{what}: -t 1 vs -t 2");
    assert_eq!(outs[1], outs[2], "{what}: -t 2 vs -t 8");
    outs[0].clone()
}

pub(crate) fn starts_with_error() -> predicates::str::StartsWithPredicate {
    predicates::str::starts_with("error: ")
}

/// The BCF magic after BGZF decompression.
pub(crate) fn is_bcf(bytes: &[u8]) -> bool {
    let mut all = Vec::new();
    noodles_bgzf::io::Reader::new(bytes)
        .read_to_end(&mut all)
        .is_ok()
        && all.starts_with(b"BCF\x02\x02")
}

/// Names in `dir` other than `keep`.
pub(crate) fn leftovers(dir: &Path, keep: &[&str]) -> Vec<String> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| !keep.contains(&n.as_str()))
        .collect()
}

/// One failing run of a P2 command: which output it writes decides what
/// "nothing left behind" means.
pub(crate) enum Sink {
    /// `-o DIR/out.EXT`
    File(&'static str),
    /// `--out-prefix DIR/out`
    Prefix,
}

/// Run `tpl` with `{IN}` replaced by the bad input, once per thread count:
/// it exits 1, says `msg`, and writes nothing.
pub(crate) fn assert_fails_clean(
    name: &str,
    kind: &str,
    tpl: &[&str],
    bad: &str,
    sink: &Sink,
    msg: &str,
) {
    for t in ["1", "4"] {
        let d = tempdir();
        let out_dir = d.path().join("out");
        std::fs::create_dir(&out_dir).unwrap();
        let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", bad)).collect();
        match sink {
            Sink::File(ext) => {
                args.extend(["-o".into(), s_of(&out_dir.join(format!("out.{ext}")))])
            }
            Sink::Prefix => args.extend(["--out-prefix".into(), s_of(&out_dir.join("out"))]),
        }
        args.extend(["-t".into(), t.into()]);
        let o = gloci().args(&args).output().unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert_eq!(o.status.code(), Some(1), "{name} / {kind} -t {t}: {err}");
        assert!(
            err.contains(msg),
            "{name} / {kind} -t {t}: wanted `{msg}` in: {err}"
        );
        assert!(err.starts_with("error: "), "{name} / {kind} -t {t}: {err}");
        let left: Vec<_> = std::fs::read_dir(&out_dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert!(
            left.is_empty(),
            "{name} / {kind} -t {t}: left behind {left:?}"
        );
    }
}

/// With the output directory missing the run exits 1, names the path and
/// how to fix it, and creates nothing.
pub(crate) fn assert_missing_dir(name: &str, tpl: &[&str], good: &str, sink: &Sink) {
    let d = tempdir();
    let missing = d.path().join("no/such/dir");
    let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", good)).collect();
    match sink {
        Sink::File(ext) => args.extend(["-o".into(), s_of(&missing.join(format!("out.{ext}")))]),
        Sink::Prefix => args.extend(["--out-prefix".into(), s_of(&missing.join("out"))]),
    }
    let o = gloci().args(&args).output().unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert_eq!(o.status.code(), Some(1), "{name}: {args:?}: {err}");
    assert!(
        err.contains("no/such/dir")
            && err.contains("No such file or directory")
            && err.contains("check that the directory exists"),
        "{name}: {err}"
    );
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0, "{name}");
}

/// `tpl` with `{IN}` replaced by `input` exits 0 and writes its output.
pub(crate) fn assert_succeeds(name: &str, tpl: &[&str], input: &str, sink: &Sink) {
    let d = tempdir();
    let mut args: Vec<String> = tpl.iter().map(|a| a.replace("{IN}", input)).collect();
    match sink {
        Sink::File(ext) => args.extend(["-o".into(), s_of(&d.path().join(format!("out.{ext}")))]),
        Sink::Prefix => args.extend(["--out-prefix".into(), s_of(&d.path().join("out"))]),
    }
    let o = gloci().args(&args).output().unwrap();
    assert_eq!(
        o.status.code(),
        Some(0),
        "{name}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(
        std::fs::read_dir(d.path()).unwrap().count() > 0,
        "{name}: no output"
    );
}
