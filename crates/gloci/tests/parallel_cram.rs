//! Per-chromosome runs and regions on indexed CRAM inputs (D41, D107, D109): the output equals the
//! streamed output at every thread count and equals the same reads in an indexed BAM, errors name
//! one place, and the one-time warnings and `-v` counts are printed once.
mod common;
use common::*;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn aln(name: &str) -> String {
    s(&inputs().join("aln").join(name))
}

fn reference() -> String {
    s(&inputs().join("fa/ref.fa"))
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[String]) -> (Option<i32>, Vec<u8>, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code(),
        o.stdout,
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

fn mk(v: &[&str]) -> Vec<String> {
    v.iter().map(|x| x.to_string()).collect()
}

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

/// The record lines of SAM text.
fn records(sam: &[u8]) -> Vec<u8> {
    String::from_utf8_lossy(sam)
        .lines()
        .filter(|l| !l.starts_with('@'))
        .collect::<Vec<_>>()
        .join("\n")
        .into_bytes()
}

const BED: &str =
    "chrM\t0\t100\nchr1\t0\t300\nchr1\t600\t900\nchr2\t100\t900\nchr10\t0\t500\nchrX\t100\t700\n";

/// The BAM of `cram` and its `.bai`, in `dir`.
fn bam_of(dir: &Path, cram: &str) -> String {
    let out = dir.join("copy.bam");
    gloci()
        .args(["view", cram, "--reference", &reference(), "-o", &s(&out)])
        .assert()
        .success();
    gloci().args(["index", &s(&out)]).assert().success();
    s(&out)
}

#[test]
fn cram_per_chromosome_equals_streaming() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    for file in ["many.cram", "many.multi.cram"] {
        let cram = aln(file);
        let bam = bam_of(dir.path(), &cram);
        let r = reference();
        let mut cases: Vec<Vec<String>> = Vec::new();
        for (a, b) in [(&cram, &bed), (&bed, &cram), (&cram, &cram)] {
            for mode in ["any", "count", "none"] {
                let mut v = mk(&["intersect", "-a", a, "-b", b, "--mode", mode]);
                if a == &cram && mode != "count" {
                    v.extend(mk(&["-O", "sam"]));
                }
                cases.push(v);
            }
            for cov in [&[][..], &["--hist"]] {
                let mut v = mk(&["coverage", "-a", a, "-b", b]);
                v.extend(mk(cov));
                cases.push(v);
            }
        }
        for extra in [
            &[][..],
            &["--bga"],
            &["--split-deletions"],
            &["--5p"],
            &["--per-base"],
            &["--hist"],
        ] {
            let mut v = mk(&["genomecov", &cram]);
            v.extend(mk(extra));
            cases.push(v);
        }
        cases.push(mk(&["convert", &cram]));
        for region in [
            &["chr1:100-400"][..],
            &["chr2:1-50", "chr1:390-420", "chrX"],
        ] {
            let mut v = mk(&["view", "-O", "sam", "--no-header", &cram]);
            v.extend(mk(region));
            cases.push(v);
        }
        cases.push(mk(&[
            "intersect",
            "-a",
            &cram,
            "-b",
            &bed,
            "--region",
            "chr2:100-500",
            "--mode",
            "count",
        ]));
        for base in cases {
            let at = |t: &str| {
                let mut v = base.clone();
                v.extend(mk(&["--reference", &r, "-t", t]));
                let (code, out, err) = run(&v);
                (code, no_cl(&out), err)
            };
            let one = at("1");
            assert_eq!(one.0, Some(0), "{file} {base:?}: {}", one.2);
            for t in ["2", "8"] {
                assert_eq!(one, at(t), "{file} {base:?} -t {t}");
            }
            // The same reads in a BAM: what they give, bar the header's names and the file.
            let bam_args: Vec<String> = base
                .iter()
                .map(|a| if a == &cram { bam.clone() } else { a.clone() })
                .collect();
            let mut v = bam_args;
            v.extend(mk(&["-t", "8"]));
            let (code, out, err) = run(&v);
            assert_eq!(code, Some(0), "BAM {base:?}: {err}");
            assert_eq!(
                records(&one.1),
                records(&no_cl(&out)),
                "{file} {base:?} against the BAM"
            );
        }
    }
}

/// `-t 8` on an indexed CRAM really reads per chromosome (and `-t 1` and a missing index
/// stream).
#[test]
fn cram_runs_per_chromosome_when_indexed() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let r = reference();
    let args = |cram: &str, t: &str| {
        mk(&[
            "-vv",
            "intersect",
            "-a",
            &bed,
            "-b",
            cram,
            "--reference",
            &r,
            "--mode",
            "count",
            "-t",
            t,
        ])
    };
    let (_, _, err) = run(&args(&aln("many.cram"), "8"));
    assert!(err.contains("running per chromosome"), "{err}");
    let (_, _, err) = run(&args(&aln("many.cram"), "1"));
    assert!(!err.contains("running per chromosome"), "{err}");
    // A copy with no .crai streams, with the same output.
    let bare = dir.path().join("bare.cram");
    std::fs::copy(aln("many.cram"), &bare).unwrap();
    let (code, out, err) = run(&args(&s(&bare), "8"));
    assert_eq!(code, Some(0), "{err}");
    assert!(!err.contains("running per chromosome"), "{err}");
    let (_, want, _) = run(&args(&aln("many.cram"), "8"));
    assert_eq!(out, want);
}

#[test]
fn cram_per_chromosome_verbose_counts_and_warnings_once() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let r = reference();
    for file in ["many.cram", "many.multi.cram"] {
        let go = |t: &str| {
            let args = mk(&[
                "intersect",
                "-a",
                &aln(file),
                "-b",
                &bed,
                "--reference",
                &r,
                "-O",
                "sam",
                "-v",
                "-t",
                t,
            ]);
            let (code, out, err) = run(&args);
            assert_eq!(code, Some(0), "{err}");
            let stats: Vec<String> = err
                .lines()
                .filter(|l| l.contains("unmapped") || l.contains("skipped"))
                .map(String::from)
                .collect();
            (no_cl(&out), stats)
        };
        let one = go("1");
        assert!(
            one.1.iter().any(|l| l.contains("unmapped reads left out")),
            "{file}: {:?}",
            one.1
        );
        assert_eq!(one, go("2"), "{file}");
        assert_eq!(one, go("8"), "{file}");
    }
}

/// A corrupt container on an A chromosome gives one message at every `-t`.
#[test]
fn cram_error_same_at_every_t() {
    let dir = tempdir();
    let src = aln("many.cram");
    let mut bytes = std::fs::read(&src).unwrap();
    // A byte inside the third data container of chr2 (its CRC no longer matches).
    let index = gloci_io::CramIndex::read(Path::new(&format!("{src}.crai"))).unwrap();
    let mut chr2: Vec<u64> = index
        .entries
        .iter()
        .filter(|e| e.ref_id == 2)
        .map(|e| e.container_offset)
        .collect();
    chr2.dedup();
    let at = chr2[0] as usize + 60;
    bytes[at] ^= 0x55;
    let bad = dir.path().join("bad.cram");
    std::fs::write(&bad, &bytes).unwrap();
    let crai = dir.path().join("bad.cram.crai");
    std::fs::copy(format!("{src}.crai"), &crai).unwrap();
    // The index is as new as the file.
    std::fs::File::options()
        .append(true)
        .open(&crai)
        .unwrap()
        .set_modified(std::time::SystemTime::now())
        .unwrap();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let r = reference();
    let at = |t: &str, a: &str, b: &str| {
        let (code, _, err) = run(&mk(&[
            "intersect",
            "-a",
            a,
            "-b",
            b,
            "--reference",
            &r,
            "--mode",
            "count",
            "-t",
            t,
        ]));
        let errors: Vec<String> = err
            .lines()
            .filter(|l| l.starts_with("error"))
            .map(String::from)
            .collect();
        (code, errors)
    };
    for (a, b) in [(&bed, &s(&bad)), (&s(&bad), &bed)] {
        let one = at("1", a, b);
        assert_eq!(one.0, Some(1), "{one:?}");
        assert!(
            one.1.join("\n").contains("bad.cram") && one.1.join("\n").contains("container"),
            "{one:?}"
        );
        for t in ["2", "8"] {
            assert_eq!(one, at(t, a, b), "-t {t}");
        }
    }
}

/// A wrong reference is the same error at every `-t`.
#[test]
fn cram_reference_error_same_at_every_t() {
    let dir = tempdir();
    let fa = dir.path().join("ref.fa");
    let text = std::fs::read_to_string(reference()).unwrap();
    // The first base of each line of chr2 changed: its slices' MD5 no longer match.
    let mut out = String::new();
    let mut in_chr2 = false;
    for l in text.lines() {
        if l.starts_with('>') {
            in_chr2 = l.starts_with(">chr2");
        } else if in_chr2 {
            let first = if l.starts_with('A') { 'C' } else { 'A' };
            out.push_str(&format!("{first}{}\n", &l[1..]));
            continue;
        }
        out.push_str(l);
        out.push('\n');
    }
    std::fs::write(&fa, out).unwrap();
    gloci().args(["index", &s(&fa)]).assert().success();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let at = |t: &str| {
        let (code, _, err) = run(&mk(&[
            "intersect",
            "-a",
            &bed,
            "-b",
            &aln("many.cram"),
            "--reference",
            &s(&fa),
            "--mode",
            "count",
            "-t",
            t,
        ]));
        let errors: Vec<String> = err
            .lines()
            .filter(|l| l.starts_with("error"))
            .map(String::from)
            .collect();
        (code, errors)
    };
    let one = at("1");
    assert_eq!(one.0, Some(1), "{one:?}");
    assert!(one.1.join("\n").contains("does not match"), "{one:?}");
    for t in ["2", "8"] {
        assert_eq!(one, at(t), "-t {t}");
    }
}

/// `--sort` of a CRAM sorts it into a temporary BAM and reads that.
#[test]
fn cram_sort_flag_matches() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let r = reference();
    let base = |extra: &[&str]| {
        let mut v = mk(&[
            "intersect",
            "-a",
            &bed,
            "-b",
            &aln("many.cram"),
            "--reference",
            &r,
            "--mode",
            "count",
        ]);
        v.extend(mk(extra));
        run(&v)
    };
    let want = base(&[]);
    assert_eq!(want.0, Some(0));
    assert_eq!(base(&["--sort"]).1, want.1);
}

/// CRAM 3.1 (order-1 rANS qualities, tokenised names) decodes in the chunk workers of a
/// per-chromosome run: they need the large stack the decode pool's threads have (a 2 MiB
/// default overflowed it, in release builds). The result equals the streamed one and the CRAM
/// 3.0 file of the same reads.
#[test]
fn cram_31_per_chromosome_at_every_t() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "t.bed", BED));
    let r = reference();
    let v31 = aln("many.v31.cram");
    let cases = [
        vec!["genomecov", "@"],
        vec!["genomecov", "@", "--bga", "--split-deletions"],
        vec!["intersect", "-a", "@", "-b", &bed, "--mode", "count"],
        vec!["intersect", "-a", &bed, "-b", "@", "--mode", "count"],
        vec!["coverage", "-a", &bed, "-b", "@"],
    ];
    for base in cases {
        let at = |file: &str, t: &str| {
            let mut v: Vec<String> = base
                .iter()
                .map(|a| {
                    if *a == "@" {
                        file.to_string()
                    } else {
                        a.to_string()
                    }
                })
                .collect();
            v.extend(mk(&["--reference", &r, "-t", t]));
            let (code, out, err) = run(&v);
            assert_eq!(code, Some(0), "{base:?} -t {t}: {err}");
            out
        };
        let want = at(&v31, "1");
        assert!(!want.is_empty());
        for t in ["2", "8"] {
            assert_eq!(at(&v31, t), want, "{base:?} -t {t}");
        }
        assert_eq!(at(&aln("many.cram"), "8"), want, "{base:?} against 3.0");
    }
}

/// A `.crai` made for another file, or older than its CRAM, blames the index: the message names
/// the index and says to rebuild it (D26), rather than calling the CRAM corrupt.
#[test]
fn cram_wrong_and_stale_index_blame_the_index() {
    let dir = tempdir();
    let r = reference();
    let cram = dir.path().join("a.cram");
    std::fs::copy(aln("many.cram"), &cram).unwrap();
    // The index of another CRAM: its offsets land inside this one's containers.
    std::fs::copy(aln("flags.cram.crai"), dir.path().join("a.cram.crai")).unwrap();
    for args in [
        mk(&["genomecov", &s(&cram), "--reference", &r, "-t", "4"]),
        mk(&["view", &s(&cram), "--reference", &r, "-O", "sam", "chr1"]),
    ] {
        let (code, _, err) = run(&args);
        assert_eq!(code, Some(1), "{args:?}: {err}");
        assert!(
            err.contains("a.cram.crai")
                && err.contains("out of date or belongs to another file")
                && err.contains("gloci index"),
            "{args:?}: {err}"
        );
    }
    // The right index, older than the CRAM: a warning, at -t 4 and for regions.
    std::fs::copy(aln("many.cram.crai"), dir.path().join("a.cram.crai")).unwrap();
    std::thread::sleep(std::time::Duration::from_millis(1100));
    std::fs::File::options()
        .append(true)
        .open(&cram)
        .unwrap()
        .set_modified(std::time::SystemTime::now())
        .unwrap();
    for args in [
        mk(&["genomecov", &s(&cram), "--reference", &r, "-t", "4"]),
        mk(&[
            "view",
            &s(&cram),
            "--reference",
            &r,
            "-O",
            "sam",
            "chr1:1-10",
        ]),
    ] {
        let (code, _, err) = run(&args);
        assert_eq!(code, Some(0), "{args:?}: {err}");
        assert_eq!(err.matches("a.cram.crai is older than").count(), 1, "{err}");
    }
}
