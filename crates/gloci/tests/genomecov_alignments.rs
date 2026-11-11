//! `genomecov` on alignments (D105, D117, D124): the header genome, the blocks, per-chromosome
//! runs that equal the streamed output at every thread count, the memory bound and `--sort`.
mod common;
use common::*;
use gloci_core::{ChromOrder, ChromTable, OrderSource};
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, ReaderOptions, RecordWriter, chrom_table_for, open_reader,
};
use gloci_ops::genomecov::{GenomecovMode, GenomecovOptions, genomecov_counting};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[&str]) -> (Option<i32>, String, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code(),
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

fn ok(args: &[&str]) -> String {
    let (code, out, err) = run(args);
    assert_eq!(code, Some(0), "{args:?}: {err}");
    out
}

fn sam_text(contigs: &[(&str, u32)], reads: &[(usize, u32, u16, &str)]) -> String {
    let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, l) in contigs {
        t.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    for (i, (c, pos, flag, cigar)) in reads.iter().enumerate() {
        // The sequence is as long as the operations that read it.
        let mut len = 0usize;
        let mut num = 0usize;
        for ch in cigar.chars() {
            if let Some(d) = ch.to_digit(10) {
                num = num * 10 + d as usize;
            } else {
                if "MIS=X".contains(ch) {
                    len += num;
                }
                num = 0;
            }
        }
        t.push_str(&format!(
            "r{i}\t{flag}\t{}\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t{}\n",
            contigs[*c].0,
            "A".repeat(len.max(1)),
            "I".repeat(len.max(1))
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

const CONTIGS: [(&str, u32); 5] = [
    ("chr1", 20_000),
    ("chr2", 20_000),
    ("chr3", 5_000),
    ("chr10", 20_000),
    ("chrX", 20_000),
];

/// Reads on every contig but chr3, with unmapped, reverse, duplicate, spliced and deleted ones,
/// and two reads with no place at the end.
fn reads_sam() -> String {
    let mut reads = Vec::new();
    for c in [0usize, 1, 3, 4] {
        for k in 0..60u32 {
            let flag = match k % 9 {
                3 => 16,
                5 => 4,
                7 => 1024,
                _ => 0,
            };
            let cigar = match k % 11 {
                2 => "5M300N5M",
                6 => "6M3D6M",
                _ => "30M",
            };
            reads.push((c, 1 + k * 37, flag, cigar));
        }
    }
    let mut t = sam_text(&CONTIGS, &reads);
    t.push_str("u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\nu2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
    t
}

#[test]
fn genomecov_alignment_parallel_equals_streaming() {
    let dir = tempdir();
    let bam = bam_of(dir.path(), "a.bam", &reads_sam());
    gloci().args(["index", &s(&bam)]).assert().success();
    let b = s(&bam);
    let modes: [&[&str]; 8] = [
        &[],
        &["--bga"],
        &["--bedgraph"],
        &["--per-base"],
        &["--bga", "--split"],
        &["--bga", "--split-deletions"],
        &["--per-base", "--strand", "-"],
        &["--bedgraph", "--scale", "0.5", "--min-mapq", "10"],
    ];
    for extra in modes {
        let at = |t: &str| {
            let mut v = vec!["genomecov", b.as_str()];
            v.extend_from_slice(extra);
            v.extend(["-t", t, "-v"]);
            let (code, out, err) = run(&v);
            assert_eq!(code, Some(0), "{extra:?} -t {t}: {err}");
            // The -v counts: the same lines at every thread count.
            let counts: Vec<&str> = err.lines().filter(|l| l.contains("left out")).collect();
            assert!(
                counts.iter().any(|l| l.contains("unmapped")),
                "{extra:?}: {err}"
            );
            (out, counts.join("\n"))
        };
        let one = at("1");
        assert!(!one.0.is_empty(), "{extra:?}");
        for t in ["2", "8"] {
            assert_eq!(one, at(t), "{extra:?} -t {t}");
        }
    }
    // The per-chromosome path ran at -t 2 (the debug line names it).
    let (_, _, err) = run(&["genomecov", &b, "-t", "2", "-vv"]);
    assert!(err.contains("per chromosome"), "{err}");
    let (_, _, err) = run(&["genomecov", &b, "-t", "1", "-vv"]);
    assert!(!err.contains("per chromosome"), "{err}");
}

#[test]
fn genomecov_alignment_errors_are_the_same_at_every_thread_count() {
    let dir = tempdir();
    let bam = bam_of(dir.path(), "a.bam", &reads_sam());
    gloci().args(["index", &s(&bam)]).assert().success();
    // A -g that is shorter than the reads need, and one that lacks a chromosome.
    let short = plain_file(
        dir.path(),
        "short.genome",
        "chr1\t20000\nchr2\t100\nchr3\t5000\nchr10\t20000\nchrX\t20000\n",
    );
    let lacks = plain_file(dir.path(), "lacks.genome", "chr1\t20000\nchr2\t20000\n");
    let b = s(&bam);
    for g in [&short, &lacks] {
        let at = |t: &str| {
            let (code, _, err) = run(&["genomecov", &b, "-g", &s(g), "-t", t]);
            let errs: Vec<&str> = err.lines().filter(|l| l.starts_with("error")).collect();
            (code, errs.join("\n"))
        };
        let one = at("1");
        assert_eq!(one.0, Some(1), "{one:?}");
        assert!(!one.1.is_empty());
        for t in ["2", "8"] {
            assert_eq!(one, at(t), "-t {t}");
        }
    }
    let (_, _, err) = run(&["genomecov", &b, "-g", &s(&short)]);
    assert!(err.contains("extends past the end of chr2"), "{err}");
    assert!(err.contains("record"), "{err}");
}

#[test]
fn genomecov_bam_no_genome_needed() {
    let dir = tempdir();
    let sam = sam_text(
        &[("chr1", 100), ("chr2", 50), ("chr3", 30)],
        &[(0, 11, 0, "10M"), (1, 1, 0, "5M")],
    );
    let bam = s(&bam_of(dir.path(), "a.bam", &sam));
    // The header's lengths, in header order, empty chromosomes included.
    assert_eq!(
        ok(&["genomecov", &bam, "--bga"]),
        "chr1\t0\t10\t0\nchr1\t10\t20\t1\nchr1\t20\t100\t0\n\
chr2\t0\t5\t1\nchr2\t5\t50\t0\nchr3\t0\t30\t0\n"
    );
    // The given genome wins.
    let g = s(&plain_file(
        dir.path(),
        "g.genome",
        "chr1\t200\nchr2\t50\nchr3\t30\n",
    ));
    assert_eq!(
        ok(&["genomecov", &bam, "--bga", "-g", &g]),
        "chr1\t0\t10\t0\nchr1\t10\t20\t1\nchr1\t20\t200\t0\n\
chr2\t0\t5\t1\nchr2\t5\t50\t0\nchr3\t0\t30\t0\n"
    );
    // SAM reads the same.
    let sp = dir.path().join("a.bam.sam");
    assert_eq!(
        ok(&["genomecov", &s(&sp), "--bga"]),
        ok(&["genomecov", &bam, "--bga"])
    );
    // A BED input still needs -g, and --split-deletions needs an alignment.
    let bed = s(&plain_file(dir.path(), "a.bed", "chr1\t0\t10\n"));
    let (code, _, err) = run(&["genomecov", &bed]);
    assert_eq!(code, Some(2), "{err}");
    assert!(
        err.contains("-g is required unless the input is BAM, CRAM or SAM"),
        "{err}"
    );
    let (code, _, err) = run(&["genomecov", &bed, "-g", &g, "--split-deletions"]);
    assert_eq!(code, Some(2), "{err}");
    assert!(
        err.contains("--split-deletions") && err.contains("not BAM, CRAM or SAM"),
        "{err}"
    );
}

#[test]
fn genomecov_alignment_blocks_follow_the_shared_rule() {
    let dir = tempdir();
    let sam = sam_text(
        &[("chr1", 100)],
        &[
            (0, 1, 0, "5M10N5M"),
            (0, 31, 0, "5M2D5M"),
            (0, 61, 16, "4M"),
        ],
    );
    let bam = s(&bam_of(dir.path(), "a.bam", &sam));
    let bg = |extra: &[&str]| {
        let mut v = vec!["genomecov", bam.as_str(), "--bedgraph"];
        v.extend_from_slice(extra);
        ok(&v)
    };
    // A read is its whole span, deleted bases and skipped bases included.
    assert_eq!(
        bg(&[]),
        "chr1\t0\t20\t1\nchr1\t30\t42\t1\nchr1\t60\t64\t1\n"
    );
    // --split cuts at N, --split-deletions at N and D.
    assert_eq!(
        bg(&["--split"]),
        "chr1\t0\t5\t1\nchr1\t15\t20\t1\nchr1\t30\t42\t1\nchr1\t60\t64\t1\n"
    );
    assert_eq!(
        bg(&["--split-deletions"]),
        "chr1\t0\t5\t1\nchr1\t15\t20\t1\nchr1\t30\t35\t1\nchr1\t37\t42\t1\nchr1\t60\t64\t1\n"
    );
    assert_eq!(
        bg(&["--split", "--split-deletions"]),
        bg(&["--split-deletions"])
    );
    // The strand is the read's direction; --5p and --3p take its ends.
    assert_eq!(bg(&["--strand", "-"]), "chr1\t60\t64\t1\n");
    assert_eq!(
        bg(&["--5p"]),
        "chr1\t0\t1\t1\nchr1\t30\t31\t1\nchr1\t63\t64\t1\n"
    );
    // The blocks are the ones `convert --blocks --split-deletions` writes.
    let blocks = ok(&["convert", &bam, "--blocks", "--split-deletions"]);
    let want: String = blocks
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            format!("{}\t{}\t{}\t1\n", f[0], f[1], f[2])
        })
        .collect();
    assert_eq!(bg(&["--split-deletions"]), want);
    // The two options cannot go with --5p or --3p.
    let (code, _, err) = run(&["genomecov", &bam, "--split-deletions", "--5p"]);
    assert_eq!(code, Some(2), "{err}");
}

/// A sink that drops every row.
struct Null;

impl RecordWriter for Null {
    fn write_record(&mut self, _: &[u8]) -> std::io::Result<()> {
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), gloci_io::IoError> {
        Ok(())
    }
}

/// A hundred thousand reads of 100 bases at depth 30 over a long chromosome: the events held
/// follow the depth, not the length of the chromosome or the file.
#[test]
fn genomecov_bam_events_bounded() {
    let dir = tempdir();
    let n = 100_000u32;
    let len = n * 10 / 3 + 1000;
    let mut t = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len}\n");
    let (seq, qual) = ("A".repeat(100), "I".repeat(100));
    for i in 0..n {
        t.push_str(&format!(
            "r{i}\t0\tchr1\t{}\t60\t100M\t*\t0\t0\t{seq}\t{qual}\n",
            1 + i * 10 / 3
        ));
    }
    let path = s(&plain_file(dir.path(), "deep.sam", &t));
    for mode in [
        GenomecovMode::Hist,
        GenomecovMode::BedGraph,
        GenomecovMode::PerBase,
    ] {
        let ro = ReaderOptions {
            accept: Accept {
                alignments: Some(AlignOptions::default()),
                ..Accept::default()
            },
            chrom_map: None,
            threads: NonZeroUsize::MIN,
            format: None,
            vcf_span: gloci_core::VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
        };
        let reader = open_reader(&path, &ro).unwrap();
        let mut chroms: ChromTable =
            chrom_table_for(None, &[reader.as_ref()], ChromOrder::Karyotype).unwrap();
        assert!(matches!(
            chroms.order_source(),
            OrderSource::AlignmentHeader(_)
        ));
        let genome =
            gloci_core::Genome::new(vec![(b"chr1".to_vec(), Some(u64::from(len)))]).unwrap();
        let o = GenomecovOptions {
            mode,
            ..GenomecovOptions::default()
        };
        let peak = genomecov_counting(reader, &mut chroms, &genome, &o, &mut Null).unwrap();
        assert!(peak <= 120, "{mode:?}: {peak} events pending at once");
        assert!(peak >= 30, "{mode:?}: the depth is 30, got {peak}");
    }
}

#[test]
fn genomecov_sort_flag_bam() {
    let dir = tempdir();
    let sorted_sam = reads_sam();
    let sorted = s(&bam_of(dir.path(), "sorted.bam", &sorted_sam));
    // The same reads, shuffled: the header chromosomes' reads reversed and interleaved.
    let (head, body): (Vec<&str>, Vec<&str>) = sorted_sam.lines().partition(|l| l.starts_with('@'));
    let mut body = body;
    body.reverse();
    let shuffled_sam = format!("{}\n{}\n", head.join("\n"), body.join("\n"));
    let shuffled = s(&bam_of(dir.path(), "shuffled.bam", &shuffled_sam));
    for extra in [&[][..], &["--bga"], &["--bga", "--split-deletions"]] {
        let mut want = vec!["genomecov", sorted.as_str()];
        want.extend_from_slice(extra);
        let want = ok(&want);
        let mut got = vec!["genomecov", shuffled.as_str(), "--sort", "--max-mem", "16M"];
        got.extend_from_slice(extra);
        assert_eq!(ok(&got), want, "{extra:?}");
    }
    // Without --sort the order error says how to fix it.
    let (code, _, err) = run(&["genomecov", &shuffled]);
    assert_eq!(code, Some(1), "{err}");
    assert!(err.contains("--sort"), "{err}");
    assert!(err.contains("gloci sort"), "{err}");
    assert_eq!(err.matches("--sort").count(), 1, "{err}");
}

/// A read cut at its reference's end warns once, at every thread count, also when a later read
/// ends the run with an error against a shorter `-g`.
#[test]
fn genomecov_overhang_warning_prints_once_at_every_thread_count() {
    let dir = tempdir();
    let sam = sam_text(
        &[("c1", 1000), ("c2", 1000)],
        &[(0, 979, 0, "10M"), (0, 991, 0, "20M")],
    );
    let bam = bam_of(dir.path(), "a.bam", &sam);
    gloci().args(["index", &s(&bam)]).assert().success();
    let g = s(&plain_file(dir.path(), "g.genome", "c1\t990\nc2\t1000\n"));
    for t in ["1", "2", "8"] {
        let (code, _, err) = run(&["genomecov", &s(&bam), "-g", &g, "-t", t]);
        assert_eq!(code, Some(1), "-t {t}: {err}");
        assert_eq!(
            err.matches("extends past the end of the reference").count(),
            1,
            "-t {t}: {err}"
        );
        assert!(err.contains("record 2 of c1"), "-t {t}: {err}");
    }
}

#[test]
fn sorted_cram_uses_its_reference() {
    // --sort replaces the CRAM with a sorted BAM copy; the reference was still used, and the
    // output is the streamed one.
    let inputs = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs");
    let cram = s(&inputs.join("aln/many.cram"));
    let fa = s(&inputs.join("fa/ref.fa"));
    let dir = tempdir();
    let tmp = s(dir.path());
    let (code, sorted, err) = run(&[
        "genomecov",
        &cram,
        "--reference",
        &fa,
        "--bedgraph",
        "--sort",
        "--tmp-dir",
        &tmp,
    ]);
    assert_eq!(code, Some(0), "{err}");
    assert!(!err.contains("not used"), "{err}");
    assert_eq!(
        sorted,
        ok(&["genomecov", &cram, "--reference", &fa, "--bedgraph"])
    );
}
