//! CRAM input (D109, D110): every CRAM fixture decodes to its BAM fixture's records, through the
//! interval reader and the every-record source alike; the reference is checked (missing,
//! mismatched, BGZF); truncated and corrupt files are errors, never panics; the pool gives the
//! same output at every thread count; gloci's `.crai` equals samtools'.

use gloci_core::{ChromOrder, ChromTable, Interval, Strand, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures::bgzip;
use std::io::Cursor;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn aln(name: &str) -> PathBuf {
    inputs().join("aln").join(name)
}

fn reference() -> PathBuf {
    inputs().join("fa/ref.fa")
}

/// The CRAM fixtures made from `set` (make-cram.sh).
fn variants(set: &str) -> Vec<String> {
    let kinds: &[&str] = match set {
        "pastend" => &["", ".embed"],
        "md" => &[""],
        _ => &["", ".v31", ".multi", ".embed", ".noref"],
    };
    kinds.iter().map(|k| format!("{set}{k}.cram")).collect()
}

fn ropts(reference: Option<PathBuf>, ao: AlignOptions, threads: usize) -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions { reference, ..ao }),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::new(threads).unwrap(),
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn full() -> AlignOptions {
    AlignOptions {
        full_records: true,
        ..AlignOptions::default()
    }
}

/// The BAM record `raw` as CRAM keeps it: CRAM does not tell `=` and `X` from `M`, so a run of
/// them comes back as one `M` (samtools decodes it so too).
fn as_cram_keeps(raw: &[u8]) -> Vec<u8> {
    let l_name = usize::from(raw[8]);
    let n = usize::from(u16::from_le_bytes([raw[12], raw[13]]));
    let at = 32 + l_name;
    let mut ops: Vec<(u32, u32)> = Vec::new();
    for w in raw[at..at + 4 * n].chunks_exact(4) {
        let v = u32::from_le_bytes([w[0], w[1], w[2], w[3]]);
        let (op, len) = (if matches!(v & 0xf, 7 | 8) { 0 } else { v & 0xf }, v >> 4);
        match ops.last_mut() {
            Some((o, l)) if *o == 0 && op == 0 => *l += len,
            _ => ops.push((op, len)),
        }
    }
    let mut out = raw[..at].to_vec();
    out[12..14].copy_from_slice(&(ops.len() as u16).to_le_bytes());
    for (op, len) in ops {
        out.extend_from_slice(&(len << 4 | op).to_le_bytes());
    }
    out.extend_from_slice(&raw[at + 4 * n..]);
    out
}

#[derive(Debug, PartialEq)]
struct Got {
    chrom: String,
    key: u64,
    span: Interval,
    strand: Strand,
    loc: u64,
    raw: Vec<u8>,
}

/// What the interval reader makes of a file: its records, its counts, and its format.
#[derive(Debug)]
struct Intervals {
    recs: Vec<Got>,
    counts: (u64, u64, u64),
    format: Format,
}

fn counts(d: &Diagnostics) -> (u64, u64, u64) {
    (d.unmapped, d.filtered, d.off_reference)
}

fn intervals(path: &Path, o: &ReaderOptions) -> Result<Intervals, IoError> {
    let mut r = open_reader(&path.to_string_lossy(), o)?;
    let mut table = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut recs = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut table, &mut batch, 2)?;
        for (i, m) in batch.recs.iter().enumerate() {
            recs.push(Got {
                chrom: String::from_utf8_lossy(table.name(m.chrom)).into_owned(),
                key: m.key,
                span: m.span,
                strand: m.strand,
                loc: m.loc,
                raw: batch.raw(i).to_vec(),
            });
        }
        if !more {
            break;
        }
    }
    Ok(Intervals {
        recs,
        counts: counts(r.diagnostics()),
        format: r.format(),
    })
}

/// Each record with its number, and the counts.
type Every = (Vec<(u64, Vec<u8>)>, (u64, u64, u64));

/// Every record through the every-record source, with its number, and the counts.
fn every(path: &Path, o: &ReaderOptions) -> Result<Every, IoError> {
    let mut s = open_aln_source(&path.to_string_lossy(), o)?;
    assert_eq!(s.layout(), Format::Bam);
    let mut out = Vec::new();
    let mut buf = Vec::new();
    while let Some(n) = s.next_record(&mut buf)? {
        out.push((n, std::mem::take(&mut buf)));
    }
    Ok((out, counts(s.diagnostics())))
}

/// The CRAM `cram` gives the records, places and counts of the BAM `bam`, through both read
/// paths, with full and slim records and with the read filters.
fn assert_same_as_bam(cram: &Path, bam: &Path, reference: Option<PathBuf>) {
    let filters = AlignOptions {
        min_mapq: 10,
        exclude_flags: 0x704,
        ..AlignOptions::default()
    };
    for ao in [full(), AlignOptions::default(), filters] {
        let mut want = intervals(bam, &ropts(None, ao.clone(), 1)).unwrap();
        for g in &mut want.recs {
            g.raw = as_cram_keeps(&g.raw);
        }
        let got = intervals(cram, &ropts(reference.clone(), ao.clone(), 1))
            .unwrap_or_else(|e| panic!("{}: {e}", cram.display()));
        assert_eq!(got.format, Format::Cram);
        assert_eq!(got.recs, want.recs, "{} ({ao:?})", cram.display());
        assert_eq!(
            got.counts,
            want.counts,
            "{} counts ({ao:?})",
            cram.display()
        );
        let mut want = every(bam, &ropts(None, ao.clone(), 1)).unwrap();
        for (_, raw) in &mut want.0 {
            *raw = as_cram_keeps(raw);
        }
        let got = every(cram, &ropts(reference.clone(), ao.clone(), 1)).unwrap();
        assert_eq!(got, want, "{} every record ({ao:?})", cram.display());
    }
}

#[test]
fn cram_equals_bam() {
    for set in ["cigars", "flags", "tags", "pastend"] {
        for v in variants(set) {
            assert_same_as_bam(&aln(&v), &aln(&format!("{set}.bam")), Some(reference()));
        }
    }
    // Unmapped reads come through the every-record source (D120).
    let (recs, _) = every(&aln("flags.cram"), &ropts(Some(reference()), full(), 1)).unwrap();
    assert_eq!(recs.len(), 14);
}

#[test]
fn cram_no_ref_needed() {
    for set in ["cigars", "flags", "tags"] {
        for v in [".noref", ".embed"] {
            let cram = aln(&format!("{set}{v}.cram"));
            assert_same_as_bam(&cram, &aln(&format!("{set}.bam")), None);
        }
    }
    // Reads past the end of chr1 decode against `N` from the embedded bases too.
    assert_same_as_bam(&aln("pastend.embed.cram"), &aln("pastend.bam"), None);
}

#[test]
fn cram_missing_reference() {
    for v in [
        "flags.cram",
        "flags.v31.cram",
        "flags.multi.cram",
        "tags.cram",
    ] {
        let first = if v.starts_with("tags") {
            "chrM"
        } else {
            "chr1"
        };
        for t in [1, 2, 8] {
            let e = every(&aln(v), &ropts(None, full(), t))
                .unwrap_err()
                .to_string();
            assert!(
                e.contains(&format!("needs the reference sequence of {first}"))
                    && e.contains("--reference FASTA")
                    && e.contains("M5 ")
                    && e.contains(v),
                "{v} -t {t}: {e}"
            );
            let e = intervals(&aln(v), &ropts(None, full(), t))
                .unwrap_err()
                .to_string();
            assert!(e.contains("needs the reference sequence"), "{v}: {e}");
        }
    }
    // A FASTA without the sequence says so.
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("other.fa");
    std::fs::write(&fa, ">chrQ\nACGT\n").unwrap();
    write_fasta_index(&fa, NonZeroUsize::MIN).unwrap();
    let e = every(&aln("flags.cram"), &ropts(Some(fa), full(), 1))
        .unwrap_err()
        .to_string();
    assert!(
        e.contains("needs the reference sequence of chr1")
            && e.contains("other.fa has no sequence chr1"),
        "{e}"
    );
}

/// ref.fa with `edit` applied to its text, indexed, in a new directory.
fn edited_reference(edit: impl Fn(&str) -> String) -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let text = std::fs::read_to_string(reference()).unwrap();
    let fa = dir.path().join("ref.fa");
    std::fs::write(&fa, edit(&text)).unwrap();
    write_fasta_index(&fa, NonZeroUsize::MIN).unwrap();
    (dir, fa)
}

/// ref.fa with chr1 one base shorter.
fn shorter_chr1(text: &str) -> String {
    let mut out = String::new();
    for rec in text
        .split_inclusive('\n')
        .collect::<String>()
        .split('>')
        .skip(1)
    {
        let (head, seq) = rec.split_once('\n').unwrap();
        let mut bases: String = seq.split('\n').collect();
        if head == "chr1" {
            bases.pop();
        }
        out.push('>');
        out.push_str(head);
        out.push('\n');
        for chunk in bases.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(chunk).unwrap());
            out.push('\n');
        }
    }
    out
}

#[test]
fn cram_mismatched_reference() {
    // chr1:101 (the first base the reads of flags.cram cover) changed: the slice checksum fails.
    let (_d, fa) = edited_reference(|t| {
        let at = t.find(">chr1\n").unwrap() + 6;
        // 100 bases in: line 2 (60 bases per line, 61 bytes), column 40.
        let i = at + 61 + 40;
        let mut b = t.as_bytes().to_vec();
        b[i] = if b[i] == b'A' { b'C' } else { b'A' };
        String::from_utf8(b).unwrap()
    });
    for t in [1, 2, 8] {
        let e = every(&aln("flags.cram"), &ropts(Some(fa.clone()), full(), t))
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("the sequence chr1 in")
                && e.contains("does not match the one")
                && e.contains("(chr1:101-701)"),
            "-t {t}: {e}"
        );
    }
    // A sequence of another length: refused before any record is read.
    let (_d, fa) = edited_reference(shorter_chr1);
    let r = open_reader(
        &aln("flags.cram").to_string_lossy(),
        &ropts(Some(fa), full(), 1),
    );
    let e = r.err().expect("refused at open").to_string();
    assert!(
        e.contains("flags.cram was not made with") && e.contains("chr1 is 1,000 bases long"),
        "{e}"
    );
}

#[test]
fn cram_reference_bgzf() {
    let dir = tempfile::tempdir().unwrap();
    let gz = dir.path().join("ref.fa.gz");
    std::fs::write(&gz, bgzip(&std::fs::read(reference()).unwrap())).unwrap();
    write_fasta_index(&gz, NonZeroUsize::MIN).unwrap();
    assert_same_as_bam(&aln("flags.cram"), &aln("flags.bam"), Some(gz.clone()));
    std::fs::remove_file(dir.path().join("ref.fa.gz.gzi")).unwrap();
    let e = every(&aln("flags.cram"), &ropts(Some(gz), full(), 1)).unwrap_err();
    assert!(matches!(e, IoError::Usage(_)), "{e:?}");
    assert!(e.to_string().contains("has no index"), "{e}");
}

/// The byte offsets where `cram`'s containers start, and its length.
fn container_starts(cram: &Path) -> (Vec<u64>, u64) {
    let mut starts: Vec<u64> = CramIndex::read(&crai_path(cram))
        .unwrap()
        .entries
        .iter()
        .map(|e| e.container_offset)
        .collect();
    starts.dedup();
    (starts, std::fs::metadata(cram).unwrap().len())
}

#[test]
fn cram_truncated() {
    let dir = tempfile::tempdir().unwrap();
    for v in ["flags.cram", "tags.multi.cram"] {
        let bytes = std::fs::read(aln(v)).unwrap();
        let (starts, len) = container_starts(&aln(v));
        let eof = len - 38;
        let mut cuts: Vec<u64> = vec![10, 26, 40, starts[0] - 1, eof, eof + 20, len - 1];
        for w in starts.windows(2) {
            cuts.extend([w[0], w[0] + 5, w[0] + 30, (w[0] + w[1]) / 2, w[1] - 1]);
        }
        cuts.extend(starts.iter().copied());
        for cut in cuts {
            let path = dir.path().join(format!("cut-{cut}-{v}"));
            std::fs::write(&path, &bytes[..cut as usize]).unwrap();
            let name = path.display().to_string();
            for t in [1, 4] {
                let e = every(&path, &ropts(Some(reference()), full(), t))
                    .expect_err("a truncated CRAM is an error")
                    .to_string();
                assert!(e.contains(&name), "cut at {cut}: {e}");
                // As a stream (no up-front check of the end), too.
                let mut s = open_aln_source_from(
                    "piped.cram",
                    Box::new(Cursor::new(bytes[..cut as usize].to_vec())),
                    &ropts(Some(reference()), full(), t),
                );
                let e = match &mut s {
                    Err(e) => e.to_string(),
                    Ok(s) => {
                        let mut buf = Vec::new();
                        loop {
                            match s.next_record(&mut buf) {
                                Ok(Some(_)) => buf.clear(),
                                Ok(None) => {
                                    panic!("cut at {cut}: a truncated stream read to its end")
                                }
                                Err(e) => break e.to_string(),
                            }
                        }
                    }
                };
                assert!(e.contains("piped.cram"), "cut at {cut}: {e}");
            }
        }
    }
}

#[test]
fn cram_corrupt_bytes_are_errors_not_panics() {
    let path = aln("tags.cram");
    let bytes = std::fs::read(&path).unwrap();
    let (starts, len) = container_starts(&path);
    let o = ropts(Some(reference()), full(), 1);
    for i in 0..bytes.len() {
        let mut b = bytes.clone();
        b[i] ^= 0x5a;
        let r = open_aln_source_from("x.cram", Box::new(Cursor::new(b)), &o).and_then(|mut s| {
            let mut buf = Vec::new();
            while s.next_record(&mut buf)?.is_some() {
                buf.clear();
            }
            Ok(())
        });
        match r {
            Err(e) => assert!(e.to_string().contains("x.cram"), "byte {i}: {e}"),
            // The file id, the header container's padding and the end-of-file container's
            // block are read by no one; every byte of the file definition's version and of a
            // data container is checked (a checksum covers each container header and block).
            Ok(()) => assert!(
                (6..26).contains(&i)
                    || (26..starts[0]).contains(&(i as u64))
                    || i as u64 >= len - 15,
                "byte {i}: a corrupt data container was read"
            ),
        }
    }
    // Not CRAM 3, and not CRAM at all past its magic.
    for (head, want) in [
        (&b"CRAM\x02\x01"[..], "CRAM 2.1, which gloci does not read"),
        (&b"CRAM\x03\x05"[..], "CRAM 3.5, which gloci does not read"),
        (
            &b"CRAM\x03\x00garbagegarbagegarbagegarbagegarbage"[..],
            "x.cram",
        ),
    ] {
        let mut b = head.to_vec();
        b.resize(b.len().max(64), b'#');
        let e = open_aln_source_from("x.cram", Box::new(Cursor::new(b)), &o)
            .err()
            .expect("refused")
            .to_string();
        assert!(e.contains(want), "{e}");
    }
}

#[test]
fn cram_pool_deterministic() {
    for set in ["cigars", "flags", "tags", "pastend"] {
        for v in variants(set) {
            let path = aln(&v);
            let one = every(&path, &ropts(Some(reference()), full(), 1)).unwrap();
            let ione =
                intervals(&path, &ropts(Some(reference()), AlignOptions::default(), 1)).unwrap();
            for t in [2, 8] {
                let got = every(&path, &ropts(Some(reference()), full(), t)).unwrap();
                assert_eq!(got, one, "{v} -t {t}");
                let igot = intervals(&path, &ropts(Some(reference()), AlignOptions::default(), t))
                    .unwrap();
                assert_eq!(igot.recs, ione.recs, "{v} -t {t}");
                assert_eq!(igot.counts, ione.counts, "{v} -t {t}");
            }
        }
    }
}

#[test]
fn crai_matches_samtools() {
    for set in ["cigars", "flags", "tags", "pastend", "md"] {
        for v in variants(set) {
            let path = aln(&v);
            let want = CramIndex::read(&crai_path(&path)).unwrap();
            assert!(!want.entries.is_empty());
            for t in [1, 2, 8] {
                let got = CramIndex::build(&path, None, NonZeroUsize::new(t).unwrap()).unwrap();
                assert_eq!(got.entries, want.entries, "{v} -t {t}");
            }
        }
    }
}

#[test]
fn crai_written_and_found() {
    let dir = tempfile::tempdir().unwrap();
    let cram = dir.path().join("x.cram");
    std::fs::copy(aln("tags.multi.cram"), &cram).unwrap();
    assert_eq!(find_crai(&cram), None);
    let p = write_crai(&cram, None, NonZeroUsize::MIN).unwrap();
    assert_eq!(p, dir.path().join("x.cram.crai"));
    assert_eq!(
        CramIndex::read(&p).unwrap(),
        CramIndex::read(&crai_path(&aln("tags.multi.cram"))).unwrap()
    );
    assert_eq!(find_crai(&cram), Some((p.clone(), false)));
    std::fs::rename(&p, dir.path().join("x.crai")).unwrap();
    assert_eq!(
        find_crai(&cram).map(|f| f.0),
        Some(dir.path().join("x.crai"))
    );
    // A reference of other lengths is refused; indexing itself needs none.
    let (_d, fa) = edited_reference(shorter_chr1);
    let e = write_crai(&cram, Some(&fa), NonZeroUsize::MIN).unwrap_err();
    assert!(e.to_string().contains("was not made with"), "{e}");
}

#[test]
fn genome_from_cram_header() {
    // No reference is needed for -g x.cram, whatever the file was written with.
    for v in ["tags.noref.cram", "tags.cram"] {
        let g = load_genome(&aln(v)).unwrap();
        let got: Vec<(String, Option<u64>)> = (0..g.len())
            .map(|i| (String::from_utf8_lossy(g.name(i)).into_owned(), g.length(i)))
            .collect();
        assert_eq!(
            got,
            [
                ("chrM", 150),
                ("chr1", 1000),
                ("chr2", 1234),
                ("chr10", 600),
                ("chrX", 777)
            ]
            .map(|(n, l)| (n.to_string(), Some(l)))
        );
    }
    let flags = aln("flags.cram").display().to_string();
    let t = chrom_table_for(Some(&flags), &[], ChromOrder::Karyotype).unwrap();
    assert_eq!(t.name(t.get(b"chr2").unwrap()), b"chr2");
}
