//! BGZF at every thread count: the compressed bytes a command writes do not depend on `-t`,
//! and an input that is not sound BGZF fails with one message and one exit status at every
//! `-t`.
mod common;
use common::*;
use gloci_testkit::fixtures;
use std::fmt::Write as _;
use std::path::Path;

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

/// Sorted BED text of about 600 KiB: ten BGZF blocks.
fn bed_text() -> String {
    let mut s = String::new();
    for i in 0..16_000u32 {
        let start = i * 40 + (i * 7919) % 13;
        writeln!(
            s,
            "chr{}\t{start}\t{}\tfeature{i}\t{}\t{}",
            1 + i / 8_000,
            start + 25,
            (i * 31) % 1000,
            if i % 3 == 0 { '-' } else { '+' }
        )
        .unwrap();
    }
    s
}

/// The lines of `text` in another order.
fn shuffled(text: &str) -> String {
    let mut lines: Vec<&str> = text.lines().collect();
    let n = lines.len();
    for i in 0..n {
        lines.swap(i, (i * 7919 + 13) % n);
    }
    lines.join("\n") + "\n"
}

/// VCF text with three samples, `records` records from `first` on.
fn vcf_text(first: u32, records: u32) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100000000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n",
    );
    for i in first..first + records {
        let gt = ["0/1\t1/1\t0/0", "0/0\t0/1\t0/1", "1/1\t0/0\t./."][(i % 3) as usize];
        writeln!(
            s,
            "chr1\t{}\tid{i}\t{}\t{}\t{}\tPASS\tDP={}\tGT\t{gt}",
            i * 5 + 1,
            ["A", "C", "G", "T"][(i % 4) as usize],
            ["C", "G", "T", "A"][(i % 4) as usize],
            30 + i % 60,
            i % 90
        )
        .unwrap();
    }
    s
}

/// Coordinate-sorted SAM text of about 1 MiB.
fn sam_text() -> String {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n");
    for i in 0..9_000u32 {
        let seq: String = (0..50)
            .map(|k| ['A', 'C', 'G', 'T'][((i * 13 + k * 7 + i / 5) % 4) as usize])
            .collect();
        let qual: String = (0..50)
            .map(|k| char::from(b'#' + ((i + k * 3) % 40) as u8))
            .collect();
        writeln!(
            s,
            "read{i}\t{}\tchr1\t{}\t{}\t50M\t*\t0\t0\t{seq}\t{qual}",
            if i % 2 == 0 { 0 } else { 16 },
            i * 100 + 1,
            i % 61
        )
        .unwrap();
    }
    s
}

/// Where each BGZF block of `bgzf` starts.
fn block_starts(bgzf: &[u8]) -> Vec<usize> {
    let mut starts = Vec::new();
    let mut at = 0;
    while at < bgzf.len() {
        starts.push(at);
        at += u16::from_le_bytes([bgzf[at + 16], bgzf[at + 17]]) as usize + 1;
    }
    assert_eq!(at, bgzf.len(), "the last block is whole");
    starts
}

/// Run `gloci ARGS -t N` in `dir` at every thread count and return the bytes of the files
/// `outs` (or of stdout when there are none), which must be the same every time and span
/// several BGZF blocks.
fn same_at_every_thread_count(dir: &Path, args: &[&str], outs: &[&str]) -> Vec<u8> {
    let mut first: Option<Vec<u8>> = None;
    for t in ["1", "2", "4", "8"] {
        for out in outs {
            let _ = std::fs::remove_file(dir.join(out));
        }
        let run = gloci()
            .current_dir(dir)
            .args(args)
            .args(["-t", t])
            .assert()
            .success();
        let bytes = if outs.is_empty() {
            run.get_output().stdout.clone()
        } else {
            outs.iter()
                .flat_map(|out| std::fs::read(dir.join(out)).unwrap())
                .collect()
        };
        match &first {
            None => {
                let data = match outs.first() {
                    Some(out) => std::fs::read(dir.join(out)).unwrap(),
                    None => bytes.clone(),
                };
                let blocks = block_starts(&data).len();
                assert!(blocks >= 4, "{args:?} writes {blocks} BGZF blocks");
                assert!(data.ends_with(&BGZF_EOF), "{args:?}");
                first = Some(bytes);
            }
            Some(first) => assert!(
                bytes == *first,
                "{args:?}: the output at -t {t} differs from the output at -t 1"
            ),
        }
    }
    first.unwrap()
}

#[test]
fn compressed_output_is_the_same_bytes_at_every_thread_count() {
    let dir = tempdir();
    let d = dir.path();
    let bed = bed_text();
    bgzf_file(d, "in.bed.gz", &bed);
    plain_file(d, "unsorted.bed", &shuffled(&bed));
    bgzf_file(d, "a.vcf.gz", &vcf_text(0, 12_000));
    bgzf_file(d, "b.vcf.gz", &vcf_text(12_000, 12_000));
    plain_file(d, "in.sam", &sam_text());

    // `view` writes its input's lines as they are, in blocks cut where noodles' writer cuts
    // them and compressed as it compresses them.
    let viewed = same_at_every_thread_count(
        d,
        &["view", "in.bed.gz", "-o", "out.bed.gz"],
        &["out.bed.gz"],
    );
    assert!(viewed == fixtures::bgzip(bed.as_bytes()));
    let piped = same_at_every_thread_count(d, &["view", "in.bed.gz", "-O", "bed.gz"], &[]);
    assert!(piped == viewed);
    same_at_every_thread_count(
        d,
        &["sort", "unsorted.bed", "-o", "out.bed.gz"],
        &["out.bed.gz"],
    );

    let cases: [(&[&str], &[&str]); 9] = [
        (
            &["sort", "unsorted.bed", "-o", "out.bed.gz", "--write-index"],
            &["out.bed.gz", "out.bed.gz.tbi"],
        ),
        (
            &[
                "view",
                "in.bed.gz",
                "-o",
                "out.bed.gz",
                "--compression-level",
                "1",
            ],
            &["out.bed.gz"],
        ),
        (&["merge", "in.bed.gz", "-o", "out.bed.gz"], &["out.bed.gz"]),
        (
            &[
                "vcf",
                "filter",
                "--min-maf",
                "0.01",
                "a.vcf.gz",
                "-o",
                "out.vcf.gz",
            ],
            &["out.vcf.gz"],
        ),
        (
            &[
                "vcf",
                "filter",
                "--min-maf",
                "0.01",
                "a.vcf.gz",
                "-o",
                "out.bcf",
            ],
            &["out.bcf"],
        ),
        (
            &["vcf", "concat", "a.vcf.gz", "b.vcf.gz", "-o", "out.vcf.gz"],
            &["out.vcf.gz"],
        ),
        (
            &["convert", "a.vcf.gz", "-o", "out.bed.gz"],
            &["out.bed.gz"],
        ),
        (&["view", "in.sam", "-o", "out.bam"], &["out.bam"]),
        (
            &["sort", "in.sam", "-o", "out.bam", "--write-index"],
            &["out.bam", "out.bam.bai"],
        ),
    ];
    for (args, outs) in cases {
        same_at_every_thread_count(d, args, outs);
    }
}

/// How a command ended: its exit status, what it wrote to stderr and stdout, and whether it
/// left its output file.
#[derive(Debug, PartialEq)]
struct Ended {
    code: Option<i32>,
    stderr: String,
    stdout: Vec<u8>,
    output_left: bool,
}

/// Run `gloci ARGS -t N` in `dir` at -t 1, 2 and 8, with `stdin` (a file of `dir`) on its
/// standard input, and return how it ended, which must be the same every time.
fn ends_the_same(dir: &Path, args: &[&str], stdin: Option<&str>, out: &str) -> Ended {
    let mut first: Option<Ended> = None;
    for t in ["1", "2", "8"] {
        let _ = std::fs::remove_file(dir.join(out));
        let mut cmd = gloci();
        cmd.current_dir(dir).args(args).args(["-t", t]);
        if let Some(file) = stdin {
            cmd.pipe_stdin(dir.join(file)).unwrap();
        }
        let o = cmd.output().unwrap();
        let ended = Ended {
            code: o.status.code(),
            stderr: String::from_utf8_lossy(&o.stderr).into_owned(),
            stdout: o.stdout,
            output_left: dir.join(out).exists(),
        };
        match &first {
            None => first = Some(ended),
            Some(first) => assert_eq!(&ended, first, "{args:?} at -t {t} and at -t 1"),
        }
    }
    first.unwrap()
}

const CORRUPT: &str = "; the file may be corrupt or the stream was cut off\n";

/// `view` of the file and of the same bytes on stdin, and `index` of the file: each fails
/// with `reason` at every thread count and leaves no output. An empty `reason` stands for
/// the inflate library's own text.
fn fails_to_read(dir: &Path, name: &str, reason: &str) {
    let index = format!("{name}.tbi");
    for (args, stdin, out, shown) in [
        (
            vec!["view", name, "-o", "out.bed.gz"],
            None,
            "out.bed.gz",
            name,
        ),
        (
            vec!["view", "-", "--format", "bed", "-o", "out.bed.gz"],
            Some(name),
            "out.bed.gz",
            "stdin",
        ),
        (vec!["index", name], None, index.as_str(), name),
    ] {
        let ended = ends_the_same(dir, &args, stdin, out);
        assert_eq!(ended.code, Some(1), "{args:?}: {}", ended.stderr);
        let head = format!("error: error reading {shown}: ");
        assert!(
            ended.stderr.starts_with(&head) && ended.stderr.ends_with(CORRUPT),
            "{args:?}: {}",
            ended.stderr
        );
        if !reason.is_empty() {
            assert_eq!(ended.stderr, format!("{head}{reason}{CORRUPT}"), "{args:?}");
        }
        assert!(!ended.output_left, "{args:?} left {out}");
        assert!(ended.stdout.is_empty());
    }
}

#[test]
fn bad_bgzf_input_fails_the_same_at_every_thread_count() {
    let dir = tempdir();
    let d = dir.path();
    let bed = bed_text();
    let good = fixtures::bgzip(bed.as_bytes());
    let starts = block_starts(&good);
    assert!(starts.len() >= 10);
    // A block in the middle, and the block after it.
    let (b, next) = (starts[5], starts[6]);
    let put = |name: &str, bytes: &[u8]| {
        fixtures::write(d, name, bytes);
    };
    let patched = |at: usize, new: &[u8]| {
        let mut bad = good.clone();
        bad[at..at + new.len()].copy_from_slice(new);
        bad
    };

    // Blocks whose content is wrong: the file ends with its end-of-file block, so it is read.
    for (name, bytes, reason) in [
        ("bad-magic.bed.gz", patched(b, b"XX"), "invalid BGZF header"),
        (
            "no-bgzf-field.bed.gz",
            patched(b + 12, b"XY"),
            "invalid BGZF header",
        ),
        (
            "bad-crc.bed.gz",
            patched(next - 8, &[!good[next - 8]]),
            "block data checksum mismatch",
        ),
        (
            "isize-70000.bed.gz",
            patched(next - 4, &70_000u32.to_le_bytes()),
            "invalid BGZF ISIZE",
        ),
        (
            "isize-large.bed.gz",
            patched(next - 4, &65_536u32.to_le_bytes()),
            "block data checksum mismatch",
        ),
        (
            "isize-small.bed.gz",
            patched(next - 4, &100u32.to_le_bytes()),
            "",
        ),
        (
            "bsize-10.bed.gz",
            patched(b + 16, &10u16.to_le_bytes()),
            "invalid frame size",
        ),
        ("bad-deflate.bed.gz", patched(b + 18, &[0xff; 64]), ""),
        (
            "cut-in-data-then-eof.bed.gz",
            [&good[..b + 500], &BGZF_EOF[..]].concat(),
            "failed to fill whole buffer",
        ),
        (
            "gzip-after-bgzf.bed.gz",
            [
                &good[..b],
                &fixtures::gzip(&bed.as_bytes()[..200_000]),
                &BGZF_EOF[..],
            ]
            .concat(),
            "invalid BGZF header",
        ),
        (
            "text-after-eof.bed.gz",
            [&good[..], &bed.as_bytes()[..200_000], &BGZF_EOF[..]].concat(),
            "invalid BGZF header",
        ),
    ] {
        put(name, &bytes);
        fails_to_read(d, name, reason);
    }

    // A file without its end-of-file block is refused before it is read; the same bytes on
    // stdin are read.
    let mut ends_between_blocks = Vec::new();
    for (name, bytes) in [
        ("cut-in-data.bed.gz", &good[..b + 500]),
        ("cut-in-header.bed.gz", &good[..b + 5]),
        ("cut-between-blocks.bed.gz", &good[..b]),
    ] {
        put(name, bytes);
        let truncated = format!(
            "error: {name}: BGZF end-of-file marker missing; the file is probably truncated\n"
        );
        for args in [vec!["view", name, "-o", "out.bed.gz"], vec!["index", name]] {
            let ended = ends_the_same(d, &args, None, "out.bed.gz");
            assert_eq!((ended.code, ended.stderr.as_str()), (Some(1), &*truncated));
            assert!(!ended.output_left);
        }
        let piped = ends_the_same(
            d,
            &["view", "-", "--format", "bed", "-o", "out.bed.gz"],
            Some(name),
            "out.bed.gz",
        );
        if name == "cut-in-data.bed.gz" {
            // The stream ends inside a block.
            assert_eq!(
                (piped.code, piped.stderr.as_str()),
                (
                    Some(1),
                    &*format!("error: error reading stdin: failed to fill whole buffer{CORRUPT}")
                )
            );
            assert!(!piped.output_left);
        } else {
            ends_between_blocks.push(piped);
        }
    }
    // A stream that ends inside a block's header reads as one that ends before that block:
    // both hold the blocks before it and nothing else.
    assert_eq!(ends_between_blocks[0], ends_between_blocks[1]);

    // Streams that hold no data, and one that is gzip but not BGZF.
    put("empty.bed.gz", b"");
    put("eof-only.bed.gz", &BGZF_EOF);
    put("plain-gzip.bed.gz", &fixtures::gzip(bed.as_bytes()));
    for name in ["empty.bed.gz", "eof-only.bed.gz"] {
        for stdin in [None, Some(name)] {
            let file = if stdin.is_some() { "-" } else { name };
            let args = ["view", file, "--format", "bed", "-o", "out.bed.gz"];
            let ended = ends_the_same(d, &args, stdin, "out.bed.gz");
            assert_eq!((ended.code, ended.stderr.as_str()), (Some(0), ""));
            assert_eq!(std::fs::read(d.join("out.bed.gz")).unwrap(), BGZF_EOF);
        }
    }
    let ended = ends_the_same(d, &["index", "empty.bed.gz"], None, "empty.bed.gz.tbi");
    assert_eq!(
        (ended.code, ended.stderr.as_str()),
        (Some(1), "error: empty.bed.gz is empty: it has no records\n")
    );
    let ended = ends_the_same(
        d,
        &["index", "eof-only.bed.gz"],
        None,
        "eof-only.bed.gz.tbi",
    );
    assert_eq!((ended.code, ended.output_left), (Some(0), true));
    let ended = ends_the_same(
        d,
        &["view", "plain-gzip.bed.gz", "-o", "out.bed.gz"],
        None,
        "out.bed.gz",
    );
    assert_eq!((ended.code, ended.stderr.as_str()), (Some(0), ""));
    assert!(std::fs::read(d.join("out.bed.gz")).unwrap() == good);
    let ended = ends_the_same(
        d,
        &["index", "plain-gzip.bed.gz"],
        None,
        "plain-gzip.bed.gz.tbi",
    );
    assert_eq!(
        (ended.code, ended.stderr.as_str()),
        (
            Some(2),
            "error: plain-gzip.bed.gz is not BGZF-compressed; recompress with: gloci view \
plain-gzip.bed.gz -o plain-gzip.bed.gz\n"
        )
    );
}

/// A bad block in a VCF, a BCF and a BAM: the `vcf` commands, `view` and the BAM index
/// builder report it the same at every thread count.
#[test]
fn a_bad_block_in_vcf_bcf_and_bam_fails_the_same_at_every_thread_count() {
    let dir = tempdir();
    let d = dir.path();
    let wrong_crc = |name: &str, good: &[u8]| {
        let starts = block_starts(good);
        assert!(starts.len() >= 8, "{name} has {} blocks", starts.len());
        let crc = starts[starts.len() / 2] - 8;
        let mut bad = good.to_vec();
        bad[crc] = !bad[crc];
        fixtures::write(d, name, &bad);
    };
    let vcf = vcf_text(0, 24_000);
    wrong_crc("bad.vcf.gz", &fixtures::bgzip(vcf.as_bytes()));
    bgzf_file(d, "good.vcf.gz", &vcf);
    gloci()
        .current_dir(d)
        .args(["view", "good.vcf.gz", "-o", "good.bcf"])
        .assert()
        .success();
    wrong_crc("bad.bcf", &std::fs::read(d.join("good.bcf")).unwrap());
    plain_file(d, "good.sam", &sam_text());
    gloci()
        .current_dir(d)
        .args(["view", "good.sam", "-o", "good.bam"])
        .assert()
        .success();
    wrong_crc("bad.bam", &std::fs::read(d.join("good.bam")).unwrap());

    let cases: [(&[&str], &str, &str); 8] = [
        (
            &[
                "vcf",
                "filter",
                "--min-maf",
                "0.01",
                "bad.vcf.gz",
                "-o",
                "out.vcf.gz",
            ],
            "out.vcf.gz",
            "bad.vcf.gz",
        ),
        (
            &["view", "bad.vcf.gz", "-o", "out.bcf"],
            "out.bcf",
            "bad.vcf.gz",
        ),
        (&["index", "bad.vcf.gz"], "bad.vcf.gz.tbi", "bad.vcf.gz"),
        (
            &[
                "vcf",
                "filter",
                "--min-maf",
                "0.01",
                "bad.bcf",
                "-o",
                "out.vcf.gz",
            ],
            "out.vcf.gz",
            "bad.bcf",
        ),
        (&["index", "bad.bcf"], "bad.bcf.csi", "bad.bcf"),
        (&["view", "bad.bam", "-o", "out.bam"], "out.bam", "bad.bam"),
        (&["index", "bad.bam"], "bad.bam.bai", "bad.bam"),
        (&["flagstat", "bad.bam"], "out", "bad.bam"),
    ];
    for (args, out, name) in cases {
        let ended = ends_the_same(d, args, None, out);
        assert_eq!(ended.code, Some(1), "{args:?}: {}", ended.stderr);
        assert_eq!(
            ended.stderr,
            format!("error: error reading {name}: block data checksum mismatch{CORRUPT}"),
            "{args:?}"
        );
        assert!(!ended.output_left, "{args:?} left {out}");
    }
}
