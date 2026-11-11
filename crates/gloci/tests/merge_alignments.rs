//! `gloci sort -m` of alignments (D122): the header rules, the order check, ties, the cascade
//! above 128 inputs, CRAM and SAM inputs, and thread counts. The merged records and headers are
//! pinned against `samtools merge -c -p` by the `aln-merge-*` golden cases.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn fixture(name: &str) -> String {
    s(&inputs().join("aln").join(name))
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn view_sam(path: &str) -> String {
    let o = gloci().args(["view", path, "-O", "sam"]).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// The lines of SAM text, without the `@PG` lines (which carry the command line).
fn lines_no_pg(sam: &str) -> Vec<String> {
    sam.lines()
        .filter(|l| !l.starts_with("@PG\tID:gloci"))
        .map(String::from)
        .collect()
}

fn names(sam: &str) -> Vec<&str> {
    sam.lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').next().unwrap())
        .collect()
}

#[test]
fn merge_header_rules() {
    let (a, b) = (fixture("merge-a.bam"), fixture("merge-b.bam"));
    let o = gloci().args(["sort", "-m", &a, &b]).assert().success();
    let bam = o.get_output().stdout.clone();
    let dir = tempdir();
    let p = dir.path().join("m.bam");
    std::fs::write(&p, bam).unwrap();
    let sam = view_sam(&s(&p));
    let count = |l: &str| sam.lines().filter(|x| *x == l).count();
    // an identical @RG is kept once, and a new one is added
    assert_eq!(count("@RG\tID:g1\tSM:s1"), 1, "{sam}");
    assert_eq!(count("@RG\tID:g2\tSM:s2"), 1, "{sam}");
    // the same @PG ID keeps the first line
    let bwa: Vec<&str> = sam
        .lines()
        .filter(|l| l.starts_with("@PG\tID:bwa"))
        .collect();
    assert_eq!(bwa.len(), 1, "{sam}");
    assert!(bwa[0].contains("a.fq"), "{sam}");
    // @CO lines from every input, in order
    let co: Vec<&str> = sam.lines().filter(|l| l.starts_with("@CO")).collect();
    assert_eq!(co, ["@CO\tmerge a", "@CO\tmerge b"]);
    assert!(
        sam.lines()
            .next()
            .unwrap()
            .starts_with("@HD\tVN:1.6\tSO:coordinate")
    );
    assert_eq!(
        names(&sam),
        [
            "ma1", "mb1", "ma2", "mb2", "mb3", "ma3", "ma4", "mb4", "mb5", "ma5", "ma6", "mb6"
        ]
    );

    // the same @RG ID with another SM: an error naming both lines
    gloci()
        .args(["sort", "-m", &a, &fixture("merge-rgclash.bam")])
        .assert()
        .code(1)
        .stderr(
            contains("@RG\tID:g1\tSM:s1")
                .and(contains("@RG\tID:g1\tSM:other"))
                .and(contains("alignment-merge-rg")),
        );
    // different @SQ: an error naming the first difference
    gloci()
        .args(["sort", "-m", &a, &fixture("merge-c.bam")])
        .assert()
        .code(1)
        .stderr(contains("reference 2 (chr2) is 100000 long").and(contains("90000")));
    // a name-sorted input in a coordinate merge
    gloci()
        .args(["sort", "-m", &a, &fixture("merge-b.n.bam")])
        .assert()
        .code(1)
        .stderr(contains("pass -n/-N, or sort it: gloci sort"));
    // and a coordinate-sorted one in a name merge
    gloci()
        .args(["sort", "-m", "-n", &fixture("merge-a.n.bam"), &b])
        .assert()
        .code(1)
        .stderr(contains("drop -n/-N"));
    // -n and -N are distinct orders
    gloci()
        .args([
            "sort",
            "-m",
            "-N",
            &fixture("merge-a.n.bam"),
            &fixture("merge-b.n.bam"),
        ])
        .assert()
        .code(1)
        .stderr(contains("SS:queryname:natural"));
}

#[test]
fn merge_checks_order() {
    // sort.bam is not coordinate-sorted
    let sort = fixture("sort.bam");
    gloci()
        .args(["sort", "-m", &sort, &sort])
        .assert()
        .code(1)
        .stderr(
            contains("sort.bam")
                .and(contains("record 2"))
                .and(contains("gloci sort")),
        );
}

fn write_sam(dir: &Path, name: &str, records: &[(&str, u16, u32)]) -> String {
    let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n");
    for (n, flag, pos) in records {
        t.push_str(&format!(
            "{n}\t{flag}\tchr1\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
        ));
    }
    s(&plain_file(dir, name, &t))
}

#[test]
fn merge_ties_by_input() {
    let dir = tempdir();
    let x = write_sam(
        dir.path(),
        "x.sam",
        &[("x1", 0, 100), ("x2", 0, 100), ("x3", 16, 100)],
    );
    let y = write_sam(dir.path(), "y.sam", &[("y1", 0, 100), ("y2", 16, 100)]);
    let z = write_sam(dir.path(), "z.sam", &[("z2", 0, 100), ("z1", 16, 100)]);
    let merged = |files: &[&str]| {
        let o = gloci().args(["sort", "-m"]).args(files).assert().success();
        String::from_utf8(o.get_output().stdout.clone()).unwrap()
    };
    // forward strand first, then reverse; within a strand the input order, then record order
    let m = merged(&[&x, &y, &z]);
    assert_eq!(names(&m), ["x1", "x2", "y1", "z2", "x3", "y2", "z1"]);
    let m = merged(&[&z, &y, &x]);
    assert_eq!(names(&m), ["z2", "y1", "x1", "x2", "z1", "y2", "x3"]);
}

#[test]
fn merge_many_inputs() {
    let dir = tempdir();
    let mut files = Vec::new();
    let mut concat = String::from(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n",
    );
    for i in 0..200u32 {
        let mut t = String::from(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n",
        );
        let mut recs = Vec::new();
        for k in 0..3u32 {
            let (chr, pos) = if k < 2 {
                ("chr1", 1 + (i * 7 + k * 31) % 50)
            } else {
                ("chr2", 1 + (i * 3) % 20)
            };
            let flag = if (i + k) % 3 == 0 { 16 } else { 0 };
            recs.push((
                chr,
                pos,
                format!(
                    "r{i}_{k}\t{flag}\t{chr}\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                ),
            ));
        }
        recs.sort_by_key(|(c, p, _)| (*c, *p));
        for (_, _, l) in &recs {
            t.push_str(l);
            concat.push_str(l);
        }
        files.push(s(&plain_file(dir.path(), &format!("f{i:03}.sam"), &t)));
    }
    let all = s(&plain_file(dir.path(), "all.sam", &concat));
    let tmp = tempdir();
    let merged = gloci()
        .args(["sort", "-m", "--tmp-dir", &s(tmp.path())])
        .args(&files)
        .assert()
        .success();
    let merged = String::from_utf8(merged.get_output().stdout.clone()).unwrap();
    let sorted = gloci().args(["sort", &all]).assert().success();
    let sorted = String::from_utf8(sorted.get_output().stdout.clone()).unwrap();
    let recs = |t: &str| {
        t.lines()
            .filter(|l| !l.starts_with('@'))
            .map(String::from)
            .collect::<Vec<_>>()
    };
    assert_eq!(recs(&merged).len(), 600);
    assert_eq!(recs(&merged), recs(&sorted));
    // the spill directory is cleaned up
    assert_eq!(std::fs::read_dir(tmp.path()).unwrap().count(), 0);
}

#[test]
fn merge_and_sort_cram_sam_inputs() {
    let dir = tempdir();
    let reference = s(&inputs().join("fa/ref.fa"));
    let cram = fixture("flags.cram");
    // BAM and SAM twins of the CRAM, in its order
    let bam = s(&dir.path().join("flags.bam"));
    gloci()
        .args([
            "view",
            &cram,
            "--reference",
            &reference,
            "-O",
            "bam",
            "-o",
            &bam,
        ])
        .assert()
        .success();
    let sam = s(&dir.path().join("flags.sam"));
    gloci()
        .args(["view", &bam, "-O", "sam", "-o", &sam])
        .assert()
        .success();
    // sort: the CRAM gives the same records as its BAM twin, written as BAM
    let sorted_cram = s(&dir.path().join("sc.bam"));
    let sorted_bam = s(&dir.path().join("sb.bam"));
    gloci()
        .args(["sort", &cram, "--reference", &reference, "-o", &sorted_cram])
        .assert()
        .success();
    gloci()
        .args(["sort", &bam, "-o", &sorted_bam])
        .assert()
        .success();
    assert_eq!(
        lines_no_pg(&view_sam(&sorted_cram)),
        lines_no_pg(&view_sam(&sorted_bam))
    );
    // sort -m: CRAM + SAM inputs equal the BAM twins merged
    let want = s(&dir.path().join("want.bam"));
    gloci()
        .args(["sort", "-m", &bam, &bam, "-o", &want])
        .assert()
        .success();
    let got = s(&dir.path().join("got.bam"));
    gloci()
        .args([
            "sort",
            "-m",
            &cram,
            &sam,
            "--reference",
            &reference,
            "-o",
            &got,
        ])
        .assert()
        .success();
    assert_eq!(lines_no_pg(&view_sam(&got)), lines_no_pg(&view_sam(&want)));
    // a CRAM needs --reference
    gloci().args(["sort", "-m", &cram, &sam]).assert().failure();
}

#[test]
fn merge_output_format_and_threads() {
    let (a, b) = (fixture("merge-a.sam"), fixture("merge-b.sam"));
    // SAM inputs write SAM; a mix of formats writes BAM
    let o = gloci().args(["sort", "-m", &a, &b]).assert().success();
    assert!(
        String::from_utf8(o.get_output().stdout.clone())
            .unwrap()
            .starts_with("@HD")
    );
    let o = gloci()
        .args(["sort", "-m", &a, &fixture("merge-b.bam")])
        .assert()
        .success();
    assert_eq!(&o.get_output().stdout[..2], &[0x1f, 0x8b]);
    // -O sam on BAM inputs
    let o = gloci()
        .args([
            "sort",
            "-m",
            "-O",
            "sam",
            &fixture("merge-a.bam"),
            &fixture("merge-b.bam"),
        ])
        .assert()
        .success();
    assert!(o.get_output().stdout.starts_with(b"@HD"));
    // byte-identical output at 1, 2 and 8 threads; only @PG CL may differ, and it strips -t
    let dir = tempdir();
    let mut outs = Vec::new();
    for t in ["1", "2", "8"] {
        let out = s(&dir.path().join(format!("t{t}.bam")));
        gloci()
            .args([
                "sort",
                "-m",
                "-n",
                &fixture("merge-a.n.bam"),
                &fixture("merge-b.n.bam"),
            ])
            .args(["-o", &out, "-t", t])
            .assert()
            .success();
        let sam = view_sam(&out);
        let pg = sam
            .lines()
            .find(|l| l.starts_with("@PG\tID:gloci"))
            .unwrap();
        assert!(pg.contains(" sort -m -n ") && !pg.contains(" -t"), "{pg}");
        outs.push(lines_no_pg(&sam));
    }
    assert_eq!(outs[0], outs[1]);
    assert_eq!(outs[0], outs[2]);
    // --check over several inputs
    gloci()
        .args(["sort", "-m", "--check", &a, &b])
        .assert()
        .code(0);
}

/// 200 small coordinate-sorted SAM files, named `f000.sam`..; returns their paths.
fn many_sams(dir: &Path, count: u32) -> Vec<String> {
    (0..count)
        .map(|i| {
            let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n");
            for k in 0..3u32 {
                let pos = 1 + (i * 7 + k * 31) % 50 + k * 100;
                t.push_str(&format!(
                    "r{i}_{k}\t0\tchr1\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                ));
            }
            s(&plain_file(dir, &format!("f{i:03}.sam"), &t))
        })
        .collect()
}

#[test]
fn merge_pipes_above_the_fan_in() {
    let dir = tempdir();
    let files = many_sams(dir.path(), 200);
    let want = gloci().args(["sort", "-m"]).args(&files).assert().success();
    let want = String::from_utf8(want.get_output().stdout.clone()).unwrap();
    let first = std::fs::read_to_string(&files[0]).unwrap();
    // stdin first, and a process substitution in the middle and at the end
    let mut args = vec!["-m".to_string(), "-".to_string()];
    args.extend(files[1..100].iter().cloned());
    args.push(format!("<(cat {})", files[100]));
    args.extend(files[101..199].iter().cloned());
    args.push(format!("<(cat {})", files[199]));
    let tmp = tempdir();
    let o = bash_gloci(
        &format!("sort {} --tmp-dir {}", args.join(" "), s(tmp.path())),
        &first,
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let got = String::from_utf8(o.stdout).unwrap();
    // the names of the pipes differ in @PG only
    let recs = |t: &str| {
        t.lines()
            .filter(|l| !l.starts_with('@'))
            .map(String::from)
            .collect::<Vec<_>>()
    };
    assert_eq!(recs(&got).len(), 600);
    assert_eq!(recs(&got), recs(&want));
    assert_eq!(std::fs::read_dir(tmp.path()).unwrap().count(), 0);
}
