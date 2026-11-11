//! Tests for `gloci reheader` (D130): a new header for a BAM or CRAM file, the rest copied.
mod common;
use common::*;
use gloci_io::bgzf_copy::{inflate, read_block};
use gloci_io::{BGZF_EOF, read_bam_header_with_extent, read_cram_prelude};
use predicates::str::contains;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn samtools_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/samtools")
        .join(name)
}

fn aln_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

/// The BAM made from SAM text.
fn make_bam(dir: &Path, name: &str, sam: &str) -> PathBuf {
    let src = write(dir, &format!("{name}.sam"), sam);
    let out = dir.join(name);
    gloci()
        .args(["view", &s(&src), "-o", &s(&out)])
        .assert()
        .success();
    out
}

/// `view --no-header -O sam` of a BAM or CRAM.
fn records(p: &Path) -> String {
    let o = gloci()
        .args(["view", &s(p), "--no-header", "-O", "sam"])
        .assert()
        .success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// The header lines of a BAM or CRAM, without `@PG` lines.
fn header_lines(p: &Path) -> Vec<String> {
    let o = gloci()
        .args(["view", &s(p), "--header-only", "-O", "sam"])
        .assert()
        .success();
    String::from_utf8(o.get_output().stdout.clone())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(String::from)
        .collect()
}

/// `(stored bytes, inflated size)` of each block of a BGZF file.
fn blocks(bytes: &[u8]) -> Vec<(Vec<u8>, usize)> {
    let mut src = bytes;
    let mut v = Vec::new();
    while let Some(b) = read_block(&mut src).unwrap() {
        v.push((b.raw, b.isize));
    }
    v
}

fn inflate_file(bytes: &[u8]) -> Vec<u8> {
    let mut src = bytes;
    let mut out = Vec::new();
    while let Some(b) = read_block(&mut src).unwrap() {
        out.extend(inflate(&b).unwrap());
    }
    out
}

/// Header and extent of a BAM file's decompressed stream.
fn extent(bam: &[u8]) -> u64 {
    let data = inflate_file(bam);
    read_bam_header_with_extent(&mut data.as_slice()).unwrap().1
}

const IN_SQ: &str = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr3\tLN:3000\n";

/// A coordinate-sorted SAM of `n` reads spread over the three references.
fn many_reads(n: usize) -> String {
    let mut sam = format!("@HD\tVN:1.6\tSO:coordinate\n{IN_SQ}@RG\tID:g1\tSM:s1\n");
    let mut x = 12345u64;
    for i in 0..n {
        let mut seq = String::new();
        let mut qual = String::new();
        for _ in 0..100 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(x >> 33) as usize % 4] as char);
            qual.push((33 + ((x >> 40) % 40) as u8) as char);
        }
        let r = i * 3 / n;
        let pos = 1 + (i % (n / 3)) % 800;
        writeln!(
            sam,
            "r{i}\t0\tchr{}\t{pos}\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tRG:Z:g1",
            r + 1
        )
        .unwrap();
    }
    sam
}

const NEW_HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n@SQ\tSN:2\tLN:2000\n\
@SQ\tSN:3\tLN:3000\n@SQ\tSN:4\tLN:4000\n@RG\tID:g1\tSM:renamed\n";

#[test]
fn reheader_bam_blocks() {
    let dir = tempdir();
    let bam = make_bam(dir.path(), "in.bam", &many_reads(30_000));
    let hdr = write(dir.path(), "new.sam", NEW_HEADER);
    let out = dir.path().join("out.bam");
    gloci()
        .args(["reheader", &s(&hdr), &s(&bam), "-o", &s(&out)])
        .assert()
        .success();
    let (inb, outb) = (std::fs::read(&bam).unwrap(), std::fs::read(&out).unwrap());
    let (in_blocks, out_blocks) = (blocks(&inb), blocks(&outb));
    assert!(in_blocks.len() > 20, "{} blocks", in_blocks.len());
    // Where the old header ends in the input's data, and the new header's in the output's.
    let old_end = extent(&inb);
    let new_end = extent(&outb);
    let (mut seen, mut hit) = (0u64, None);
    for (i, (_, n)) in in_blocks.iter().enumerate() {
        if seen + *n as u64 > old_end {
            hit = Some((i, (old_end - seen) as usize));
            break;
        }
        seen += *n as u64;
    }
    let (k, cut) = hit.expect("the header is inside the input");
    assert!(
        cut > 0 || k == 0,
        "the fixture should put the header end inside a block"
    );
    // The output's data: the new header, then the input's data from the old header's end.
    let (in_data, out_data) = (inflate_file(&inb), inflate_file(&outb));
    assert_eq!(out_data[new_end as usize..], in_data[old_end as usize..]);
    // The blocks: every input block after the one holding the header's end, byte for byte,
    // then the EOF marker; the header blocks before them end exactly at the new header's end.
    let tail: Vec<&Vec<u8>> = in_blocks[k + 1..].iter().map(|(b, _)| b).collect();
    let got: Vec<&Vec<u8>> = out_blocks[out_blocks.len() - tail.len()..]
        .iter()
        .map(|(b, _)| b)
        .collect();
    assert_eq!(got, tail);
    assert_eq!(out_blocks.last().unwrap().0, BGZF_EOF);
    let mut at = 0u64;
    let mut edge = false;
    for (_, n) in &out_blocks {
        at += *n as u64;
        edge |= at == new_end;
    }
    assert!(edge, "the new header ends on a block edge");
    // The records decode the same under the new header.
    assert_eq!(records(&out).replace("\tchr", "\t"), {
        records(&bam).replace("\tchr", "\t")
    });
    let lines = header_lines(&out);
    assert!(lines.contains(&"@SQ\tSN:4\tLN:4000".to_string()));
    let o = gloci()
        .args(["view", &s(&out), "--header-only", "-O", "sam"])
        .assert()
        .success();
    let text = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(text.contains("@PG\tID:gloci\tPN:gloci"), "{text}");
}

#[test]
fn reheader_huge_header() {
    let dir = tempdir();
    const N: usize = 200_000;
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    let mut new = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for i in 0..N {
        let len = 1000 + (i % 7000);
        writeln!(sam, "@SQ\tSN:contig_{i:06}\tLN:{len}").unwrap();
        writeln!(new, "@SQ\tSN:c{i}\tLN:{len}").unwrap();
    }
    sam.push_str("@RG\tID:g1\tSM:s1\n");
    new.push_str("@RG\tID:g1\tSM:s1\n");
    let mut want = String::new();
    for i in 0..300 {
        let c = (i * 997) % N;
        let line = format!(
            "r{i}\t0\tCONTIG\t{}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:g1\n",
            1 + i
        );
        sam.push_str(&line.replace("CONTIG", &format!("contig_{c:06}")));
        want.push_str(&line.replace("CONTIG", &format!("c{c}")));
    }
    // Reads must be in reference order for a coordinate-sorted header: sort by contig.
    let mut recs: Vec<(usize, String, String)> = Vec::new();
    for (a, b) in sam
        .lines()
        .filter(|l| !l.starts_with('@'))
        .zip(want.lines())
    {
        let name = a.split('\t').nth(2).unwrap();
        recs.push((name[7..].parse().unwrap(), a.to_string(), b.to_string()));
    }
    recs.sort();
    let head: String = sam
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| format!("{l}\n"))
        .collect();
    let sam: String = head
        + &recs
            .iter()
            .map(|r| format!("{}\n", r.1))
            .collect::<String>();
    let want: String = recs.iter().map(|r| format!("{}\n", r.2)).collect();
    let bam = make_bam(dir.path(), "huge.bam", &sam);
    let hdr = write(dir.path(), "new.sam", &new);
    let out = dir.path().join("out.bam");
    let t = Instant::now();
    gloci()
        .args(["reheader", &s(&hdr), &s(&bam), "-o", &s(&out)])
        .assert()
        .success();
    assert!(t.elapsed().as_secs() < 10, "took {:?}", t.elapsed());
    let (inb, outb) = (std::fs::read(&bam).unwrap(), std::fs::read(&out).unwrap());
    let (old_end, new_end) = (extent(&inb), extent(&outb));
    assert!(old_end > 5_000_000, "header of {old_end} bytes");
    assert!(blocks(&inb).len() > 50);
    let (in_data, out_data) = (inflate_file(&inb), inflate_file(&outb));
    assert_eq!(out_data[new_end as usize..], in_data[old_end as usize..]);
    assert_eq!(records(&out), want);
    // Every input block after the header's end is in the output, verbatim.
    let in_blocks = blocks(&inb);
    let out_blocks = blocks(&outb);
    let mut seen = 0u64;
    let k = in_blocks
        .iter()
        .position(|(_, n)| {
            seen += *n as u64;
            seen > old_end
        })
        .unwrap();
    let tail: Vec<&Vec<u8>> = in_blocks[k + 1..].iter().map(|(b, _)| b).collect();
    let got: Vec<&Vec<u8>> = out_blocks[out_blocks.len() - tail.len()..]
        .iter()
        .map(|(b, _)| b)
        .collect();
    assert_eq!(got, tail);
}

/// Replace the references of a header text: names and appended lines.
fn cram_header_text(names: [&str; 2], extra: &str) -> String {
    format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{}\tLN:1000\n@SQ\tSN:{}\tLN:1234\n{extra}@RG\tID:g1\tSM:s1\n",
        names[0], names[1]
    )
}

#[test]
fn reheader_cram() {
    let dir = tempdir();
    let hdr = write(
        dir.path(),
        "new.sam",
        &cram_header_text(["1", "2"], "@SQ\tSN:3\tLN:99\n"),
    );
    // CRAM 3.0 with the reference embedded.
    let input = aln_fixture("flags.embed.cram");
    let out = dir.path().join("out.cram");
    gloci()
        .args(["reheader", &s(&hdr), &s(&input), "-o", &s(&out)])
        .assert()
        .success();
    let (inb, outb) = (std::fs::read(&input).unwrap(), std::fs::read(&out).unwrap());
    let (pi, po) = (
        read_cram_prelude(&input).unwrap(),
        read_cram_prelude(&out).unwrap(),
    );
    assert_eq!(
        (pi.major, pi.minor, pi.file_id),
        (po.major, po.minor, po.file_id)
    );
    // The data containers, EOF container included, are the input's bytes.
    assert_eq!(
        outb[po.data_start as usize..],
        inb[pi.data_start as usize..]
    );
    assert!(String::from_utf8_lossy(&po.header_text).contains("@SQ\tSN:3\tLN:99"));
    assert_eq!(
        records(&out),
        records(&input)
            .replace("\tchr1\t", "\t1\t")
            .replace("\tchr2\t", "\t2\t")
    );
    // CRAM 3.1 needs its reference; the reference's names change with the file's.
    let refdir = dir.path().join("ref");
    std::fs::create_dir(&refdir).unwrap();
    let fa = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/fa/ref.fa"),
    )
    .unwrap();
    let renamed = fa.replace(">chr1", ">1").replace(">chr2", ">2");
    let ref_fa = write(&refdir, "ref.fa", &renamed);
    gloci().args(["index", &s(&ref_fa)]).assert().success();
    let v31 = aln_fixture("flags.v31.cram");
    let out31 = dir.path().join("out31.cram");
    gloci()
        .args(["reheader", &s(&hdr), &s(&v31), "-o", &s(&out31)])
        .assert()
        .success();
    let (inb, outb) = (std::fs::read(&v31).unwrap(), std::fs::read(&out31).unwrap());
    let (pi, po) = (
        read_cram_prelude(&v31).unwrap(),
        read_cram_prelude(&out31).unwrap(),
    );
    assert_eq!((po.major, po.minor), (3, 1));
    assert_eq!(
        outb[po.data_start as usize..],
        inb[pi.data_start as usize..]
    );
    let got = gloci()
        .args([
            "view",
            &s(&out31),
            "--no-header",
            "-O",
            "sam",
            "--reference",
            &s(&ref_fa),
        ])
        .assert()
        .success();
    let got = String::from_utf8(got.get_output().stdout.clone()).unwrap();
    let orig_ref = aln_ref_fa();
    let want = gloci()
        .args([
            "view",
            &s(&v31),
            "--no-header",
            "-O",
            "sam",
            "--reference",
            &s(&orig_ref),
        ])
        .assert()
        .success();
    let want = String::from_utf8(want.get_output().stdout.clone()).unwrap();
    assert_eq!(
        got,
        want.replace("\tchr1\t", "\t1\t")
            .replace("\tchr2\t", "\t2\t")
    );
}

fn aln_ref_fa() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/fa/ref.fa")
}

#[test]
fn reheader_cram_header_is_the_fixtures_header_text_otherwise() {
    // The new header container is read back by P3a's reader: header lines, refs, @PG.
    let dir = tempdir();
    let hdr = write(dir.path(), "new.sam", &cram_header_text(["a", "b"], ""));
    let out = dir.path().join("o.cram");
    gloci()
        .args([
            "reheader",
            &s(&hdr),
            &s(&aln_fixture("flags.embed.cram")),
            "-o",
            &s(&out),
        ])
        .assert()
        .success();
    let lines = header_lines(&out);
    assert_eq!(
        lines,
        [
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:a\tLN:1000",
            "@SQ\tSN:b\tLN:1234",
            "@RG\tID:g1\tSM:s1"
        ]
    );
}

#[test]
fn reheader_inputs() {
    let dir = tempdir();
    let bam = samtools_fixture("reheader-in.bam");
    let hdr = samtools_fixture("new-header.sam");
    // stdin as IN is a usage error that names the way out.
    gloci()
        .args(["reheader", &s(&hdr), "-"])
        .write_stdin("x")
        .assert()
        .code(2)
        .stderr(contains("gloci view - -o x.bam"));
    // `-` as HEADER works.
    let text = std::fs::read(&hdr).unwrap();
    let out = dir.path().join("o.bam");
    gloci()
        .args(["reheader", "-", &s(&bam), "-o", &s(&out)])
        .write_stdin(text)
        .assert()
        .success();
    assert!(header_lines(&out).contains(&"@SQ\tSN:4\tLN:4000".to_string()));
    // A SAM IN is a usage error.
    gloci()
        .args([
            "reheader",
            &s(&hdr),
            &s(&samtools_fixture("reheader-in.sam")),
        ])
        .assert()
        .code(2)
        .stderr(contains("SAM"));
    // HEADER may be a SAM file with records: only its header lines count.
    let with_reads = write(
        dir.path(),
        "with-reads.sam",
        &format!("{NEW_HEADER}r9\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"),
    );
    let out2 = dir.path().join("o2.bam");
    gloci()
        .args(["reheader", &s(&with_reads), &s(&bam), "-o", &s(&out2)])
        .assert()
        .success();
    let want_header: Vec<String> = NEW_HEADER.lines().map(String::from).collect();
    assert_eq!(header_lines(&out2), want_header);
    assert_eq!(records(&out2), records(&out));
    // HEADER may be a BAM or CRAM file.
    let out3 = dir.path().join("o3.bam");
    gloci()
        .args([
            "reheader",
            &s(&samtools_fixture("new-header.bam")),
            &s(&bam),
            "-o",
            &s(&out3),
        ])
        .assert()
        .success();
    assert_eq!(header_lines(&out3), header_lines(&out));
    // A HEADER without header lines is refused.
    let none = write(dir.path(), "none.sam", "r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n");
    gloci()
        .args(["reheader", &s(&none), &s(&bam)])
        .assert()
        .code(1)
        .stderr(contains("no SAM header lines"));
    // A BAM without its end-of-file marker is exit 1, with nothing written.
    let cut = dir.path().join("cut.bam");
    let bytes = std::fs::read(&bam).unwrap();
    std::fs::write(&cut, &bytes[..bytes.len() - 28]).unwrap();
    let o4 = dir.path().join("o4.bam");
    gloci()
        .args(["reheader", &s(&hdr), &s(&cut), "-o", &s(&o4)])
        .assert()
        .code(1)
        .stderr(contains("end-of-file marker"));
    assert!(!o4.exists());
    // So is a CRAM without its end-of-file container.
    gloci()
        .args(["reheader", &s(&hdr), &s(&aln_fixture("flags-noeof.cram"))])
        .assert()
        .code(1);
    // A bad @SQ list: exit 1, nothing written, the message names both lines.
    let o5 = dir.path().join("o5.bam");
    gloci()
        .args([
            "reheader",
            &s(&samtools_fixture("bad-ln.sam")),
            &s(&bam),
            "-o",
            &s(&o5),
        ])
        .assert()
        .code(1)
        .stderr(contains("@SQ 2 of"))
        .stderr(contains("(SN:chr2 LN:2001)"))
        .stderr(contains("(SN:chr2 LN:2000)"))
        .stderr(contains("--chrom-map"));
    assert!(!o5.exists());
}

#[test]
fn reheader_in_place_is_atomic() {
    let dir = tempdir();
    let target = dir.path().join("x.bam");
    std::fs::copy(samtools_fixture("reheader-in.bam"), &target).unwrap();
    let before = std::fs::read(&target).unwrap();
    // A failure leaves IN as it was, and no temporary file behind.
    gloci()
        .args([
            "reheader",
            &s(&samtools_fixture("bad-order.sam")),
            &s(&target),
            "-o",
            &s(&target),
        ])
        .assert()
        .code(1);
    assert_eq!(std::fs::read(&target).unwrap(), before);
    let names: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    assert_eq!(names.len(), 1, "{names:?}");
    // Success replaces it.
    gloci()
        .args([
            "reheader",
            &s(&samtools_fixture("new-header.sam")),
            &s(&target),
            "-o",
            &s(&target),
        ])
        .assert()
        .success();
    assert_ne!(std::fs::read(&target).unwrap(), before);
    assert!(header_lines(&target).contains(&"@SQ\tSN:1\tLN:1000".to_string()));
    assert_eq!(
        records(&target)
            .replace("\t1\t", "\tchr1\t")
            .lines()
            .count(),
        7
    );
}

#[test]
fn reheader_is_deterministic_across_threads() {
    let dir = tempdir();
    let bam = samtools_fixture("reheader-in.bam");
    let cram = aln_fixture("flags.embed.cram");
    let h1 = samtools_fixture("new-header.sam");
    let h2 = write(dir.path(), "c.sam", &cram_header_text(["x", "y"], ""));
    for (hdr, input, name) in [(&h1, &bam, "o.bam"), (&h2, &cram, "o.cram")] {
        // `-t` is not part of the @PG command line, so with the same relative `-o` name (each
        // run in its own directory) the whole files are equal, @PG included.
        let mut outs = Vec::new();
        for t in ["1", "2", "8"] {
            let wd = dir.path().join(format!("t{t}-{name}"));
            std::fs::create_dir(&wd).unwrap();
            gloci()
                .current_dir(&wd)
                .args(["reheader", &s(hdr), &s(input), "-o", name, "-t", t])
                .assert()
                .success();
            outs.push(std::fs::read(wd.join(name)).unwrap());
        }
        assert_eq!(outs[0], outs[1]);
        assert_eq!(outs[0], outs[2]);
    }
}
