//! Fixtures of the alignment and FASTA command matrices.

use crate::common::*;
use crate::helpers::*;
use gloci_testkit::fixtures::{bgzip, write};
use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

pub(crate) fn p3a_inputs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

pub(crate) fn p3a_aln(name: &str) -> String {
    s_of(&p3a_inputs().join("aln").join(name))
}

pub(crate) fn p3a_fa(name: &str) -> String {
    s_of(&p3a_inputs().join("fa").join(name))
}

/// SAM text, record lines only, with `@PG` command lines dropped (they name `-t`).
/// BAM input (gzip magic) is decoded through `gloci view` first.
pub(crate) fn p3a_norm(bytes: &[u8]) -> Vec<u8> {
    let text = if bytes.starts_with(&[0x1f, 0x8b]) {
        let d = tempdir();
        let f = write(d.path(), "x.bam", bytes);
        gloci()
            .args(["view", "-O", "sam", &s_of(&f)])
            .output()
            .unwrap()
            .stdout
    } else {
        bytes.to_vec()
    };
    String::from_utf8_lossy(&text)
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .flat_map(|l| [l, "\n"])
        .collect::<String>()
        .into_bytes()
}

/// The three alignment formats over the same 1210 reads, plus a FASTA and BED to go with them.
pub(crate) struct P3a {
    pub(crate) dir: tempfile::TempDir,
    pub(crate) cram: String,
    pub(crate) bam: String,
    pub(crate) sam: String,
    pub(crate) fa: String,
    pub(crate) fa_gz: String,
    pub(crate) bed: String,
    pub(crate) masks: String,
}

impl P3a {
    pub(crate) fn new() -> P3a {
        let dir = tempdir();
        let d = dir.path();
        let cram = d.join("many.cram");
        std::fs::copy(p3a_aln("many.cram"), &cram).unwrap();
        std::fs::copy(p3a_aln("many.cram.crai"), d.join("many.cram.crai")).unwrap();
        let fa = d.join("ref.fa");
        std::fs::copy(p3a_fa("ref.fa"), &fa).unwrap();
        std::fs::copy(p3a_fa("ref.fa.fai"), d.join("ref.fa.fai")).unwrap();
        let fa_gz = d.join("ref.fa.gz");
        for ext in ["", ".fai", ".gzi"] {
            std::fs::copy(
                p3a_fa(&format!("ref.fa.gz{ext}")),
                d.join(format!("ref.fa.gz{ext}")),
            )
            .unwrap();
        }
        let sam = d.join("many.sam");
        let bam = d.join("many.bam");
        for out in [&sam, &bam] {
            gloci()
                .args(["view", &s_of(&cram), "--reference", &s_of(&fa), "-o"])
                .arg(out)
                .assert()
                .success();
        }
        gloci().args(["index", &s_of(&bam)]).assert().success();
        // strictly sorted copies, for `sort -m` (many.cram is only position-sorted)
        for (src, name) in [(&bam, "sorted.bam"), (&sam, "sorted.sam")] {
            gloci()
                .args(["sort", &s_of(src), "-o", &s_of(&d.join(name))])
                .assert()
                .success();
        }
        let bed = plain_file(
            d,
            "t.bed",
            "chrM\t0\t100\nchr1\t0\t300\nchr1\t600\t900\nchr2\t100\t900\nchr10\t0\t500\nchrX\t100\t700\n",
        );
        let masks = plain_file(
            d,
            "masks.bed",
            "chrM\t5\t20\nchr1\t0\t30\nchr1\t100\t140\nchr2\t500\t700\nchr10\t10\t20\n",
        );
        P3a {
            cram: s_of(&cram),
            bam: s_of(&bam),
            sam: s_of(&sam),
            fa: s_of(&fa),
            fa_gz: s_of(&fa_gz),
            bed: s_of(&bed),
            masks: s_of(&masks),
            dir,
        }
    }

    pub(crate) fn path(&self, name: &str) -> String {
        s_of(&self.dir.path().join(name))
    }

    /// `(kind, input, extra args)` for each alignment format.
    pub(crate) fn kinds(&self) -> Vec<(&'static str, String, Vec<String>)> {
        vec![
            ("bam", self.bam.clone(), vec![]),
            ("sam", self.sam.clone(), vec![]),
            (
                "cram",
                self.cram.clone(),
                vec!["--reference".into(), self.fa.clone()],
            ),
        ]
    }

    /// Commands over one alignment input, as `(name, args)`. `-t` is added by the caller.
    pub(crate) fn alignment_cmds(&self, kind: &str, input: &str) -> Vec<(String, Vec<String>)> {
        let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
        let bed = self.bed.as_str();
        let mut c: Vec<(String, Vec<String>)> = vec![
            ("convert".into(), v(&["convert", input])),
            ("convert --bed12".into(), v(&["convert", input, "--bed12"])),
            ("convert --fastq".into(), v(&["convert", input, "--fastq"])),
            ("view".into(), v(&["view", "-O", "sam", input])),
            (
                "view --min-mapq".into(),
                v(&["view", "-O", "sam", "--min-mapq", "30", input]),
            ),
            (
                "intersect A".into(),
                v(&["intersect", "-a", input, "-b", bed]),
            ),
            (
                "intersect B".into(),
                v(&["intersect", "-a", bed, "-b", input]),
            ),
            (
                "intersect counts".into(),
                v(&["intersect", "-a", bed, "-b", input, "--mode", "count"]),
            ),
            ("coverage".into(), v(&["coverage", "-a", bed, "-b", input])),
            ("genomecov".into(), v(&["genomecov", input])),
            (
                "genomecov bg".into(),
                v(&["genomecov", input, "--bedgraph"]),
            ),
            ("sort".into(), v(&["sort", "-O", "sam", input])),
            ("sort -n".into(), v(&["sort", "-n", "-O", "sam", input])),
            ("sort -N".into(), v(&["sort", "-N", "-O", "sam", input])),
            (
                "sort small budget".into(),
                v(&["sort", "--max-mem", "16M", "-O", "sam", input]),
            ),
        ];
        if kind != "cram" {
            let sorted = input.replace("many.", "sorted.");
            c.push((
                "sort -m".into(),
                v(&["sort", "-m", "-O", "sam", &sorted, &sorted]),
            ));
        }
        if kind != "sam" {
            c.push((
                "view region".into(),
                v(&["view", "-O", "sam", input, "chr1:1-500", "chr2:100-900"]),
            ));
        }
        if kind == "bam" {
            let sorted = input.replace("many.", "sorted.");
            c.push(("cat".into(), v(&["cat", &sorted, &sorted])));
        }
        c
    }

    pub(crate) fn fasta_cmds(&self) -> Vec<(String, Vec<String>)> {
        let v = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
        let bed12 = p3a_fa("gf12.bed");
        let gf = p3a_fa("gf.bed");
        let nuc = p3a_fa("nuc.bed");
        let mut c = Vec::new();
        for (tag, fa) in [("plain", &self.fa), ("bgzf", &self.fa_gz)] {
            c.push((
                format!("getfasta {tag}"),
                v(&["getfasta", "--fasta", fa, &gf, "-s", "--name"]),
            ));
            c.push((
                format!("getfasta --split {tag}"),
                v(&[
                    "getfasta",
                    "--fasta",
                    fa,
                    &bed12,
                    "--split",
                    "-s",
                    "--bed-out",
                ]),
            ));
            c.push((
                format!("getfasta --tab {tag}"),
                v(&["getfasta", "--fasta", fa, &gf, "--tab"]),
            ));
            c.push((
                format!("nuc {tag}"),
                v(&["nuc", "--fasta", fa, &nuc, "-s", "--seq", "--pattern", "AC"]),
            ));
            c.push((
                format!("maskfasta {tag}"),
                v(&["maskfasta", "--fasta", fa, "-b", &self.masks]),
            ));
            c.push((
                format!("maskfasta --soft {tag}"),
                v(&["maskfasta", "--fasta", fa, "-b", &self.masks, "--soft"]),
            ));
        }
        c
    }
}

pub(crate) fn p3a_run_t(args: &[String], t: &str) -> std::process::Output {
    let opt: Vec<String> = match t {
        "1" => vec!["-t".into(), "1".into()],
        "2" => vec!["-t2".into()],
        n => vec![format!("--threads={n}")],
    };
    gloci().args(args).args(&opt).output().unwrap()
}

/// A SAM (and its BAM, indexed when sorted) of `n` 100 bp reads over chr1 and chr2 (10 Mb each),
/// in coordinate order or shuffled.
pub(crate) fn p3a_big_alignments(dir: &Path, n: usize, shuffled: bool) -> (String, String) {
    let mut x: u64 = 0xC0FFEE;
    let mut next = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x >> 20
    };
    let mut reads: Vec<(usize, u64, u64)> = (0..n)
        .map(|_| ((next() % 2) as usize, 1 + next() % 9_999_000, next()))
        .collect();
    if !shuffled {
        reads.sort_by_key(|r| (r.0, r.1, r.2 & 1));
    }
    let mut sam = format!(
        "@HD\tVN:1.6\tSO:{}\n@SQ\tSN:chr1\tLN:10000000\n@SQ\tSN:chr2\tLN:10000000\n",
        if shuffled { "unsorted" } else { "coordinate" }
    );
    for (i, (c, pos, r)) in reads.iter().enumerate() {
        let seq: String = (0..100)
            .map(|k| b"ACGT"[((r >> (k % 30)) as usize + k) & 3] as char)
            .collect();
        writeln!(
            sam,
            "read{i}\t{}\tchr{}\t{pos}\t60\t100M\t*\t0\t0\t{seq}\t{}",
            if r & 1 == 0 { 0 } else { 16 },
            c + 1,
            "I".repeat(100)
        )
        .unwrap();
    }
    let name = if shuffled { "shuf" } else { "big" };
    let sam_path = plain_file(dir, &format!("{name}.sam"), &sam);
    let bam_path = dir.join(format!("{name}.bam"));
    gloci()
        .args(["view", &s_of(&sam_path), "-o", &s_of(&bam_path)])
        .assert()
        .success();
    if !shuffled {
        gloci().args(["index", &s_of(&bam_path)]).assert().success();
    }
    (s_of(&sam_path), s_of(&bam_path))
}

/// A 3 Mb FASTA (`big.fa`, indexed) and BED of 6000 intervals over it, for output that outgrows a
/// pipe.
pub(crate) fn p3a_big_fasta(dir: &Path) -> (String, String) {
    let mut x: u64 = 99;
    let mut fa = String::from(">big\n");
    let bases = b"ACGT";
    for line in 0..50_000 {
        for _ in 0..60 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            fa.push(bases[(x >> 60) as usize & 3] as char);
        }
        fa.push('\n');
        let _ = line;
    }
    let f = plain_file(dir, "big.fa", &fa);
    gloci().args(["index", &s_of(&f)]).assert().success();
    let mut bed = String::new();
    for i in 0..6000 {
        writeln!(bed, "big\t{}\t{}\tn{i}\t0\t+", i * 400, i * 400 + 380).unwrap();
    }
    let b = plain_file(dir, "big.bed", &bed);
    (s_of(&f), s_of(&b))
}

/// Records of SAM text (no header lines).
pub(crate) fn p3a_records(sam: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(&p3a_norm(sam))
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(String::from)
        .collect()
}

/// A BAM cut inside its last record but ending with a valid BGZF end-of-file block (so only the
/// record itself is short).
pub(crate) fn p3a_cut_record_bam(dir: &Path) -> String {
    let mut raw = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(p3a_aln("cat-big.bam")).unwrap())
        .read_to_end(&mut raw)
        .unwrap();
    raw.truncate(raw.len() - 7);
    s_of(&write(dir, "cut-record.bam", &bgzip(&raw)))
}

/// Run gloci with its file size limited to `limit` bytes (writes past it fail with EFBIG, as on a
/// full disk, instead of killing the process with SIGXFSZ).
#[cfg(unix)]
pub(crate) fn p3a_run_size_limited(args: &[String], limit: u64) -> std::process::Output {
    use std::os::unix::process::CommandExt;
    let mut cmd = Command::new(bin());
    cmd.args(args);
    unsafe {
        cmd.pre_exec(move || {
            libc::signal(libc::SIGXFSZ, libc::SIG_IGN);
            let rlim = libc::rlimit {
                rlim_cur: limit,
                rlim_max: limit,
            };
            if libc::setrlimit(libc::RLIMIT_FSIZE, &rlim) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            Ok(())
        });
    }
    cmd.output().unwrap()
}
