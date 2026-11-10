//! `gloci index` of BAM, `--write-index` on BAM output, regions through the index, and stale
//! and broken indexes (D19, D26, D107, D116).
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn aln(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// Run gloci with `args`, expecting success; returns stdout.
fn run(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

/// The records `view -O sam --no-header` writes for `regions`, and stderr.
fn records(bam: &Path, regions: &[&str]) -> (Vec<u8>, String) {
    let mut args = vec!["view", "-O", "sam", "--no-header", "-t", "2"];
    let b = s(bam);
    args.push(&b);
    args.extend_from_slice(regions);
    let o = gloci().args(&args).assert().success();
    let o = o.get_output();
    (
        o.stdout.clone(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

fn copy(dir: &Path, name: &str, to: &str) -> PathBuf {
    let p = dir.join(to);
    std::fs::copy(aln(name), &p).unwrap();
    p
}

#[test]
fn write_index_on_bam_output() {
    let dir = tempdir();
    let regions = ["chr1:101-110", "chr1:700-705", "chr2:51-60", "chr1:200-210"];
    let mut bais = Vec::new();
    let mut bams = Vec::new();
    for t in ["1", "2", "8"] {
        let out = dir.path().join(format!("o{t}.bam"));
        run(&[
            "view",
            &s(&aln("flags.bam")),
            "-o",
            &s(&out),
            "--write-index",
            "-t",
            t,
        ]);
        let bai = dir.path().join(format!("o{t}.bam.bai"));
        assert!(bai.exists(), "-t {t}");
        bais.push(std::fs::read(&bai).unwrap());
        bams.push(out);
    }
    // the index is the same at every -t, and the one `gloci index` writes, and samtools'
    assert!(bais.iter().all(|b| *b == bais[0]));
    let again = dir.path().join("again.bam");
    std::fs::copy(&bams[0], &again).unwrap();
    run(&["index", &s(&again)]);
    assert_eq!(
        std::fs::read(dir.path().join("again.bam.bai")).unwrap(),
        bais[0]
    );
    // reads come back through it, as a scan finds them
    let (indexed, err) = records(&bams[0], &regions);
    assert!(!err.contains("scanning"), "{err}");
    std::fs::remove_file(dir.path().join("o1.bam.bai")).unwrap();
    let (scanned, err) = records(&bams[0], &regions);
    assert!(err.contains("scanning the whole file"), "{err}");
    assert_eq!(
        String::from_utf8_lossy(&indexed),
        String::from_utf8_lossy(&scanned)
    );
    assert!(!indexed.is_empty());

    // `cat`: two halves of a sorted file, joined
    let (a, b, c) = (
        dir.path().join("a.bam"),
        dir.path().join("b.bam"),
        dir.path().join("c.bam"),
    );
    run(&["view", &s(&aln("flags.bam")), "chr1", "-o", &s(&a)]);
    run(&["view", &s(&aln("flags.bam")), "chr2", "-o", &s(&b)]);
    let mut cat_bai = Vec::new();
    for t in ["1", "8"] {
        run(&[
            "cat",
            &s(&a),
            &s(&b),
            "-o",
            &s(&c),
            "--write-index",
            "-t",
            t,
        ]);
        cat_bai.push(std::fs::read(dir.path().join("c.bam.bai")).unwrap());
    }
    assert_eq!(cat_bai[0], cat_bai[1]);
    let (indexed, err) = records(&c, &regions);
    assert!(!err.contains("scanning"), "{err}");
    std::fs::remove_file(dir.path().join("c.bam.bai")).unwrap();
    let (scanned, _) = records(&c, &regions);
    assert_eq!(
        String::from_utf8_lossy(&indexed),
        String::from_utf8_lossy(&scanned)
    );
    assert!(!indexed.is_empty());
    // the index of the joined file is the index of the file `gloci index` indexes
    run(&["index", &s(&c)]);
    assert_eq!(
        std::fs::read(dir.path().join("c.bam.bai")).unwrap(),
        cat_bai[0]
    );

    // BAM on stdout, and SAM output, cannot be indexed
    gloci()
        .args(["view", &s(&aln("flags.bam")), "-O", "bam", "--write-index"])
        .assert()
        .code(2)
        .stderr(contains("-o FILE"));
    gloci()
        .args(["cat", &s(&a), &s(&b), "--write-index"])
        .assert()
        .code(2)
        .stderr(contains("-o FILE").and(contains("needs BGZF output").not()));
    gloci()
        .args([
            "view",
            &s(&aln("flags.bam")),
            "-o",
            &s(&dir.path().join("o.sam")),
            "--write-index",
        ])
        .assert()
        .code(2)
        .stderr(contains("gloci view").and(contains("-o x.bam")));
}

#[test]
fn index_command_for_bam_and_sam() {
    let dir = tempdir();
    let bam = copy(dir.path(), "flags.bam", "x.bam");
    run(&["index", &s(&bam)]);
    assert_eq!(
        std::fs::read(dir.path().join("x.bam.bai")).unwrap(),
        std::fs::read(aln("flags.bam.bai")).unwrap()
    );
    // --csi replaces it
    run(&["index", &s(&bam), "--csi"]);
    assert!(dir.path().join("x.bam.csi").exists());
    assert!(!dir.path().join("x.bam.bai").exists());
    assert_eq!(
        records(&bam, &["chr1:700-705"]).0,
        run(&[
            "view",
            "-O",
            "sam",
            "--no-header",
            &aln("flags.sam").to_string_lossy(),
            "chr1:700-705"
        ])
    );
    // SAM has no index
    let sam = copy(dir.path(), "flags.sam", "x.sam");
    gloci()
        .args(["index", &s(&sam)])
        .assert()
        .code(2)
        .stderr(contains(format!("gloci view {} -o x.bam", s(&sam))));
    // --csi is for BAM
    let bed = bgzf_file(dir.path(), "x.bed.gz", "chr1\t1\t2\n");
    gloci()
        .args(["index", &s(&bed), "--csi"])
        .assert()
        .code(2)
        .stderr(contains("--csi is for BAM"));
    // an unsorted BAM is an error that says what to sort
    let edge = copy(dir.path(), "edge.bam", "edge.bam");
    gloci()
        .args(["index", &s(&edge)])
        .assert()
        .code(1)
        .stderr(contains("gloci sort").and(contains("panicked").not()));
    assert!(!dir.path().join("edge.bam.bai").exists());
}

/// An old `.bai` beside a BAM that gloci overwrites describes other data (D19).
#[test]
fn stale_indexes_beside_overwritten_bam_output_are_removed() {
    let dir = tempdir();
    let out = dir.path().join("x.bam");
    for stale in ["x.bam.bai", "x.bai", "x.bam.csi"] {
        std::fs::write(dir.path().join(stale), b"old").unwrap();
    }
    run(&["view", &s(&aln("flags.bam")), "-o", &s(&out)]);
    for stale in ["x.bam.bai", "x.bai", "x.bam.csi"] {
        assert!(!dir.path().join(stale).exists(), "{stale} was left behind");
    }
    // `cat` and `sort` overwrite too
    std::fs::write(dir.path().join("x.bam.bai"), b"old").unwrap();
    run(&["cat", &s(&aln("cat-a.bam")), "-o", &s(&out)]);
    assert!(!dir.path().join("x.bam.bai").exists());
    // and with --write-index the new one replaces the old
    std::fs::write(dir.path().join("x.bam.csi"), b"old").unwrap();
    run(&[
        "view",
        &s(&aln("flags.bam")),
        "-o",
        &s(&out),
        "--write-index",
    ]);
    assert!(dir.path().join("x.bam.bai").exists());
    assert!(!dir.path().join("x.bam.csi").exists());
}

#[test]
fn stale_and_broken_indexes_are_reported() {
    let dir = tempdir();
    let bam = copy(dir.path(), "flags.bam", "x.bam");
    let bai = dir.path().join("x.bam.bai");
    std::fs::copy(aln("flags.bam.bai"), &bai).unwrap();
    // a .bai older than the BAM still answers, with a warning
    let f = std::fs::File::options().write(true).open(&bai).unwrap();
    f.set_modified(std::time::SystemTime::now() - std::time::Duration::from_secs(3600))
        .unwrap();
    let (out, err) = records(&bam, &["chr1:101-110"]);
    assert!(
        err.contains("older than") && err.contains("gloci index"),
        "{err}"
    );
    assert!(!out.is_empty());
    // a corrupt one is an error that says how to fix it, not a panic or an abort
    for bad in [&b"BAI\x01\xff\xff\xff\x7f"[..], b"junk", b""] {
        std::fs::write(&bai, bad).unwrap();
        let o = gloci()
            .args(["view", "-O", "sam", &s(&bam), "chr1:101-110"])
            .assert()
            .failure();
        let err = String::from_utf8_lossy(&o.get_output().stderr).into_owned();
        assert!(
            err.contains("x.bam.bai") && !err.contains("panicked"),
            "{err}"
        );
    }
}

/// Reads through the index equal reads by scan at every thread count, for regions that repeat,
/// overlap and sit on unmapped placed reads.
#[test]
fn indexed_regions_equal_scan_at_every_thread_count() {
    let dir = tempdir();
    let bam = copy(dir.path(), "flags.bam", "x.bam");
    let scan = copy(dir.path(), "flags.bam", "scan.bam");
    let regions = [
        vec!["chr1:650-701", "chr2:51-55"],
        vec!["chr2", "chr1:101-110", "chr1:105-120"],
        vec!["chr1:701-701"],
        vec!["chr1:700-700"],
        vec!["chr1:1-100000", "nope:1-5"],
    ];
    for r in &regions {
        let mut args: Vec<String> = ["view", "-O", "sam", "--no-header"]
            .iter()
            .map(|a| a.to_string())
            .collect();
        let mut want = None;
        for t in ["1", "2", "8"] {
            // `scan.bam` has no index
            for (file, indexed) in [(&bam, true), (&scan, false)] {
                if indexed && !dir.path().join("x.bam.bai").exists() {
                    run(&["index", &s(file)]);
                }
                let mut a = args.clone();
                a.extend(["-t".into(), t.into(), s(file)]);
                a.extend(r.iter().map(|x| x.to_string()));
                let out = gloci()
                    .args(&a)
                    .assert()
                    .success()
                    .get_output()
                    .stdout
                    .clone();
                let want = want.get_or_insert_with(|| out.clone());
                assert_eq!(&out, want, "{r:?} -t {t} indexed {indexed}");
            }
        }
        args.clear();
    }
}
