//! Fixtures of the statistics, randomness, pairs and VCF command matrices.

use crate::common::*;
use crate::fixtures_alignments::*;
use crate::fixtures_bed_vcf::*;
use crate::helpers::*;
use gloci_testkit::fixtures::{bgzip, vcf_to_bcf, write};
use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use std::process::Stdio;

pub(crate) fn p3b_xtask(rel: &str) -> String {
    s_of(&p3a_inputs().join(rel))
}

/// Sorted BED6: `per` intervals on each of chr1..chr`chroms`, `step` apart from `off`.
pub(crate) fn p3b_bed(chroms: usize, per: usize, step: usize, off: usize, w: usize) -> String {
    let mut t = String::new();
    for c in 1..=chroms {
        for i in 0..per {
            let st = off + i * step + (i % 3) * 7;
            writeln!(
                t,
                "chr{c}\t{st}\t{}\tf{c}_{i}\t{}\t{}",
                st + w + i % 45,
                i % 100,
                if i % 2 == 0 { '+' } else { '-' }
            )
            .unwrap();
        }
    }
    t
}

/// BEDPE pairs in no particular order: ends 500 bases apart on one chromosome, some on two.
pub(crate) fn p3b_bedpe(n: usize, salt: usize) -> String {
    let mut t = String::new();
    for i in 0..n {
        let c = 1 + (i + salt) % 4;
        let s = 100 + ((i * 7919 + salt) % 3000) * 300;
        let c2 = if i % 17 == 0 { 1 + (c % 4) } else { c };
        writeln!(
            t,
            "chr{c}\t{s}\t{}\tchr{c2}\t{}\t{}\tp{salt}_{i}\t{}\t{}\t{}",
            s + 50 + i % 40,
            s + 500,
            s + 560,
            i % 100,
            if i % 2 == 0 { '+' } else { '-' },
            if i % 3 == 0 { '+' } else { '-' }
        )
        .unwrap();
    }
    t
}

pub(crate) fn p3b_vcf(chroms: std::ops::RangeInclusive<usize>, per: usize) -> String {
    p2_vcf(chroms, per).replace("length=10000>", "length=100000000>")
}

/// `p3b_vcf` with a second alternate allele on every third record.
pub(crate) fn p3b_multi_vcf(per: usize) -> String {
    p3b_vcf(1..=3, per)
        .lines()
        .enumerate()
        .flat_map(|(k, l)| {
            let l = if !l.starts_with('#') && k % 3 == 0 {
                l.replacen("\tA\tC\t", "\tA\tC,G\t", 1)
            } else {
                l.to_string()
            };
            [l, "\n".to_string()]
        })
        .collect()
}

pub(crate) fn p3b_rename_samples(vcf: &str) -> String {
    vcf.replace("\ts1\ts2\ts3\n", "\tt1\tt2\tt3\n")
}

/// A 4-chromosome FASTA, indexed, long enough for `p3b_vcf` records up to 3000 per chromosome.
pub(crate) fn p3b_fasta(dir: &Path) -> String {
    let mut x: u64 = 99;
    let mut t = String::new();
    for c in 1..=4 {
        writeln!(t, ">chr{c}").unwrap();
        let seq: Vec<u8> = (0..120_000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(x >> 40) as usize & 3]
            })
            .collect();
        for l in seq.chunks(60) {
            t.push_str(std::str::from_utf8(l).unwrap());
            t.push('\n');
        }
    }
    let p = s_of(&plain_file(dir, "ref.fa", &t));
    gloci().args(["index", &p]).assert().success();
    p
}

/// A name-sorted SAM of `n` proper pairs on chr1.
pub(crate) fn p3b_pairs_sam(n: usize) -> String {
    let mut t = String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:10000000\n");
    let seq = "ACGT".repeat(12) + "AC";
    let qual = "I".repeat(50);
    for i in 0..n {
        let pos = 1 + i * 100;
        writeln!(
            t,
            "p{i:07}\t99\tchr1\t{pos}\t60\t50M\t=\t{}\t250\t{seq}\t{qual}",
            pos + 200
        )
        .unwrap();
        writeln!(
            t,
            "p{i:07}\t147\tchr1\t{}\t60\t50M\t=\t{pos}\t-250\t{seq}\t{qual}",
            pos + 200
        )
        .unwrap();
    }
    t
}

pub(crate) struct P3b {
    /// Keeps the fixture directory alive.
    #[allow(dead_code)]
    pub(crate) fx: P2,
    pub(crate) genome: String,
    pub(crate) a: String,
    pub(crate) b: String,
    pub(crate) unsorted_a: String,
    pub(crate) pe_a: String,
    pub(crate) pe_b: String,
    pub(crate) vcf: String,
    pub(crate) vcf_gz: String,
    pub(crate) vcf_bcf: String,
    vcf2: String,
    pub(crate) multi: String,
    pub(crate) renamed: String,
    pub(crate) fasta: String,
    pub(crate) pairs_sam: String,
}

impl P3b {
    pub(crate) fn new() -> P3b {
        let fx = P2::new();
        let genome = fx.plain(
            "g.genome",
            "chr1\t1000000\nchr2\t1000000\nchr3\t1000000\nchr4\t1000000\n",
        );
        let a_txt = p3b_bed(4, 3000, 300, 100, 60);
        let a = fx.plain("a.bed", &a_txt);
        let b = fx.plain("b.bed", &p3b_bed(4, 2000, 440, 130, 90));
        let mut lines: Vec<&str> = a_txt.lines().collect();
        lines.reverse();
        let unsorted_a = fx.plain("ua.bed", &(lines.join("\n") + "\n"));
        let pe_a = fx.plain("a.bedpe", &p3b_bedpe(6000, 0));
        let pe_b = fx.plain("b.bedpe", &p3b_bedpe(300, 0));
        let vtxt = p3b_vcf(1..=4, 3000);
        let vcf = fx.plain("v.vcf", &vtxt);
        let vcf_gz = fx.indexed("v.vcf.gz", &vtxt);
        let vcf_bcf = fx.bcf("v.bcf", &vtxt);
        let vcf2 = fx.plain("v2.vcf", &p3b_vcf(2..=4, 3000));
        let multi = fx.plain("multi.vcf", &p3b_multi_vcf(3000));
        let renamed = fx.plain("t.vcf", &p3b_rename_samples(&p3b_vcf(2..=4, 3000)));
        let fasta = p3b_fasta(fx.dir.path());
        let pairs_sam = fx.plain("pairs.sam", &p3b_pairs_sam(6000));
        P3b {
            fx,
            genome,
            a,
            b,
            unsorted_a,
            pe_a,
            pe_b,
            vcf,
            vcf_gz,
            vcf_bcf,
            vcf2,
            multi,
            renamed,
            fasta,
            pairs_sam,
        }
    }

    /// Commands that write to standard output, with a name each.
    pub(crate) fn stdout_cmds(&self) -> Vec<(String, Vec<String>)> {
        let (a, b, g, ua) = (&self.a, &self.b, &self.genome, &self.unsorted_a);
        let (pa, pb) = (&self.pe_a, &self.pe_b);
        let mut c: Vec<(String, Vec<String>)> = Vec::new();
        let mut add = |n: &str, args: Vec<String>| c.push((n.to_string(), args));
        add("jaccard", sv(&["stats", "jaccard", "-a", a, "-b", b]));
        add(
            "jaccard -f -r -s",
            sv(&[
                "stats", "jaccard", "-a", a, "-b", b, "-f", "0.4", "-r", "-s",
            ]),
        );
        add(
            "jaccard --sort",
            sv(&["stats", "jaccard", "-a", ua, "-b", b, "--sort"]),
        );
        add(
            "fisher",
            sv(&["stats", "fisher", "-a", a, "-b", b, "-g", g]),
        );
        add(
            "fisher -m --sort",
            sv(&[
                "stats", "fisher", "-a", ua, "-b", b, "-g", g, "-m", "--sort",
            ]),
        );
        add(
            "fisher --exclude",
            sv(&["stats", "fisher", "-a", a, "-b", b, "-g", g, "--exclude", b]),
        );
        add("reldist", sv(&["stats", "reldist", "-a", a, "-b", b]));
        add(
            "reldist --detail",
            sv(&["stats", "reldist", "-a", a, "-b", b, "--detail"]),
        );
        add(
            "reldist --sort",
            sv(&["stats", "reldist", "-a", ua, "-b", b, "--sort"]),
        );
        add(
            "reldist --unsorted",
            sv(&["stats", "reldist", "-a", ua, "-b", b, "--unsorted"]),
        );
        add(
            "random",
            sv(&["random", "-g", g, "-n", "30000", "--seed", "11"]),
        );
        add(
            "random -l",
            sv(&["random", "-g", g, "-n", "5000", "-l", "250", "--seed", "12"]),
        );
        add(
            "random bed.gz",
            sv(&[
                "random", "-g", g, "-n", "5000", "--seed", "13", "-O", "bed.gz",
            ]),
        );
        add("sample", sv(&["sample", "-n", "500", "--seed", "3", a]));
        add(
            "sample --strand",
            sv(&["sample", "-n", "500", "--seed", "4", "--strand", "+", a]),
        );
        add(
            "sample vcf",
            sv(&["sample", "-n", "300", "--seed", "5", &self.vcf]),
        );
        add("shuffle", sv(&["shuffle", "-g", g, "--seed", "9", a]));
        add(
            "shuffle --same-chrom",
            sv(&["shuffle", "-g", g, "--seed", "9", "--same-chrom", a]),
        );
        add(
            "shuffle --chrom-first",
            sv(&["shuffle", "-g", g, "--seed", "9", "--chrom-first", a]),
        );
        add(
            "shuffle --exclude --no-overlapping",
            sv(&[
                "shuffle",
                "-g",
                g,
                "--seed",
                "9",
                "--exclude",
                b,
                "--no-overlapping",
                a,
            ]),
        );
        add(
            "shuffle bedpe",
            sv(&["shuffle", "-g", g, "--seed", "9", "--format", "bedpe", pa]),
        );
        for ty in [
            "either", "both", "xor", "neither", "notboth", "ispan", "ospan", "notispan", "notospan",
        ] {
            add(
                &format!("pair-intersect {ty}"),
                sv(&["pair-intersect", "-a", pa, "-b", b, "--type", ty]),
            );
        }
        add(
            "pair-intersect -f -r -s",
            sv(&["pair-intersect", "-a", pa, "-b", b, "-f", "0.3", "-r", "-s"]),
        );
        add(
            "pair-intersect vcf B",
            sv(&["pair-intersect", "-a", pa, "-b", &self.vcf]),
        );
        add(
            "pair-intersect bedpe.gz",
            sv(&["pair-intersect", "-a", pa, "-b", b, "-O", "bedpe.gz"]),
        );
        for ty in ["both", "either", "neither", "notboth"] {
            add(
                &format!("pair-pair {ty}"),
                sv(&["pair-pair", "-a", pa, "-b", pb, "--type", ty]),
            );
        }
        add(
            "pair-pair -w",
            sv(&["pair-pair", "-a", pa, "-b", pb, "-w", "400"]),
        );
        add(
            "pair-pair -f -r -s",
            sv(&["pair-pair", "-a", pa, "-b", pb, "-f", "0.3", "-r", "-s"]),
        );
        add(
            "convert --bedpe",
            sv(&["convert", "--bedpe", &self.pairs_sam]),
        );
        add(
            "convert --bedpe fixture sam",
            sv(&["convert", "--bedpe", &p3b_xtask("aln/pairs.sam")]),
        );
        add(
            "convert --bedpe fixture bam",
            sv(&["convert", "--bedpe", &p3b_xtask("aln/pairs.bam")]),
        );
        add(
            "convert --bam bedpe",
            sv(&["convert", "--bam", "--format", "bedpe", pa, "-g", g]),
        );
        add("convert --bam", sv(&["convert", "--bam", a, "-g", g]));
        add(
            "convert --fastq fixture pairs",
            sv(&[
                "convert",
                "--fastq",
                &p3b_xtask("aln/fastq-pairs.bam"),
                "-o",
                "/dev/stdout",
                "--fastq2",
                "/dev/stdout",
            ]),
        );
        add(
            "convert --fastq --tags fixture",
            sv(&[
                "convert",
                "--fastq",
                "--tags",
                &p3b_xtask("aln/fastq-tags.bam"),
                "-o",
                "/dev/stdout",
                "--fastq2",
                "/dev/stdout",
            ]),
        );
        add("vcf split", sv(&["vcf", "split", &self.multi]));
        add(
            "vcf split --types",
            sv(&["vcf", "split", &self.multi, "--types", "snps"]),
        );
        add(
            "vcf norm",
            sv(&["vcf", "norm", &self.vcf, "-f", &self.fasta, "-c", "s"]),
        );
        add(
            "vcf norm --split",
            sv(&[
                "vcf",
                "norm",
                &self.multi,
                "-f",
                &self.fasta,
                "-c",
                "s",
                "--split",
            ]),
        );
        add(
            "vcf annotate",
            sv(&[
                "vcf",
                "annotate",
                &self.vcf,
                "-a",
                &self.vcf2,
                "-c",
                "ID,INFO/DP",
                "-m",
                "+HIT",
            ]),
        );
        add(
            "vcf merge",
            sv(&["vcf", "merge", &self.vcf2, &self.renamed, "-m", "none"]),
        );
        add(
            "vcf merge both",
            sv(&["vcf", "merge", &self.vcf2, &self.renamed, "-m", "both"]),
        );
        for (tag, f) in [
            ("plain", &self.vcf),
            ("bgzf", &self.vcf_gz),
            ("bcf", &self.vcf_bcf),
        ] {
            add(
                &format!("vcf filter {tag}"),
                sv(&["vcf", "filter", f, "--min-qual", "30", "--pass-only"]),
            );
            add(
                &format!("vcf extract {tag}"),
                sv(&["vcf", "extract", f, "--info", "DP"]),
            );
            add(
                &format!("vcf extract --format {tag}"),
                sv(&["vcf", "extract", f, "--format", "GT"]),
            );
        }
        c
    }
}

/// The files `vcf compare`, `vcf stats` and `vcf popgen` write, at one thread count.
pub(crate) fn p3b_prefix_files(args: &[String], exts: &[&str], t: &str) -> Vec<Vec<u8>> {
    let d = tempdir();
    let prefix = s_of(&d.path().join("o"));
    let o = gloci()
        .args(args)
        .args(["--out-prefix", &prefix, "-t", t])
        .output()
        .unwrap();
    assert_eq!(
        o.status.code(),
        Some(0),
        "{args:?} -t {t}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    exts.iter()
        .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap_or_else(|_| panic!("no .{e}")))
        .collect()
}

pub(crate) fn p3b_prefix_cmds(fx: &P3b) -> Vec<(String, Vec<String>, Vec<&'static str>)> {
    let cmp = ["diff.sites", "diff.indv", "diff.discordance_matrix"];
    let stats = ["frq", "lqual", "imiss", "het", "singletons", "snpden"];
    let pop1 = fx.fx.plain("fst-pop1.txt", "s1\n");
    let pop2 = fx.fx.plain("fst-pop2.txt", "s2\ns3\n");
    let mut c = vec![
        (
            "vcf compare".to_string(),
            sv(&[
                "vcf",
                "compare",
                &fx.vcf,
                &fx.vcf2,
                "--site-discordance",
                "--sample-discordance",
                "--discordance-matrix",
            ]),
            cmp.to_vec(),
        ),
        (
            "vcf compare bcf/bgzf --min-gq".to_string(),
            sv(&[
                "vcf",
                "compare",
                &fx.vcf_bcf,
                &fx.vcf_gz,
                "--min-dp",
                "3",
                "--site-discordance",
                "--sample-discordance",
                "--discordance-matrix",
            ]),
            cmp.to_vec(),
        ),
    ];
    for (tag, f) in [
        ("plain", &fx.vcf),
        ("bgzf", &fx.vcf_gz),
        ("bcf", &fx.vcf_bcf),
    ] {
        c.push((
            format!("vcf stats {tag}"),
            sv(&[
                "vcf",
                "stats",
                f,
                "--freq",
                "--site-quality",
                "--missing-indv",
                "--het",
                "--singletons",
                "--snp-density",
                "1000",
            ]),
            stats.to_vec(),
        ));
        c.push((
            format!("vcf popgen {tag}"),
            sv(&[
                "vcf",
                "popgen",
                f,
                "--site-pi",
                "--window-pi",
                "1000",
                "--window-pi-step",
                "300",
                "--tajima-d",
                "700",
            ]),
            ["sites.pi", "windowed.pi", "Tajima.D"].to_vec(),
        ));
        c.push((
            format!("vcf popgen fst {tag}"),
            sv(&[
                "vcf",
                "popgen",
                f,
                "--fst-pop",
                &pop1,
                "--fst-pop",
                &pop2,
                "--fst-window-size",
                "1000",
                "--fst-window-step",
                "300",
            ]),
            ["weir.fst", "weir.fst.summary", "windowed.weir.fst"].to_vec(),
        ));
        c.push((
            format!("vcf popgen ld {tag}"),
            sv(&[
                "vcf",
                "popgen",
                f,
                "--geno-r2",
                "--hap-r2",
                "--ld-window",
                "10",
                "--min-dp",
                "3",
            ]),
            ["geno.ld", "hap.ld"].to_vec(),
        ));
    }
    c
}

/// Output larger than a pipe, closed after 64 bytes: exit 0 and silent, as for the P2 commands.
pub(crate) fn p3b_assert_broken_pipe(name: &str, args: &[String], min_bytes: usize) {
    let full = gloci().args(args).output().unwrap();
    assert_eq!(
        full.status.code(),
        Some(0),
        "{name}: {}",
        String::from_utf8_lossy(&full.stderr)
    );
    assert!(
        full.stdout.len() > min_bytes,
        "{name}: only {} bytes of output",
        full.stdout.len()
    );
    for t in ["1", "4"] {
        let mut child = Command::new(bin())
            .args(args)
            .args(["-t", t])
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut first = [0u8; 64];
        child
            .stdout
            .as_mut()
            .unwrap()
            .read_exact(&mut first)
            .unwrap();
        drop(child.stdout.take());
        let out = child.wait_with_output().unwrap();
        assert_eq!(
            out.status.code(),
            Some(0),
            "{name} -t {t}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            out.stderr.is_empty(),
            "{name} -t {t}: stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
}

/// Lines of a text output that are not `##glociCommand`.
pub(crate) fn p3b_text(bytes: &[u8]) -> String {
    String::from_utf8(without_command(bytes)).unwrap()
}

pub(crate) fn p3b_header_lines(bytes: &[u8]) -> Vec<String> {
    p3b_text(bytes)
        .lines()
        .filter(|l| l.starts_with('#'))
        .map(String::from)
        .collect()
}

/// Bad VCF inputs of the P3b matrices: sorted rows with a swapped pair, a bad POS, and
/// truncated BGZF / BCF files; `renamed` gives them the sample names `t1..t3` (for merge).
pub(crate) struct P3bBadVcf {
    pub(crate) good: String,
    pub(crate) unsorted: String,
    pub(crate) malformed: String,
    pub(crate) trunc_gz: String,
    pub(crate) nofooter_gz: String,
    pub(crate) trunc_bcf: String,
}

pub(crate) fn p3b_bad_vcf(fx: &P2, renamed: bool) -> P3bBadVcf {
    let tag = if renamed { "r" } else { "s" };
    let adj = |t: String| if renamed { p3b_rename_samples(&t) } else { t };
    let text = adj(p3b_vcf(1..=2, 40));
    let rows: Vec<&str> = text.lines().collect();
    let hdr = rows.iter().filter(|l| l.starts_with('#')).count();
    let body: Vec<&str> = rows[hdr..].to_vec();
    let head = rows[..hdr].join("\n");
    let mut swapped = body.clone();
    swapped.swap(5, 6);
    let bad_pos = {
        let mut f: Vec<&str> = body[6].split('\t').collect();
        f[1] = "XX";
        f.join("\t")
    };
    let mut malformed = body.clone();
    malformed[6] = &bad_pos;
    let gz = bgzip(text.as_bytes());
    let bcf = vcf_to_bcf(text.as_bytes());
    P3bBadVcf {
        good: fx.plain(&format!("good_{tag}.vcf"), &text),
        unsorted: fx.plain(
            &format!("uns_{tag}.vcf"),
            &format!("{head}\n{}\n", swapped.join("\n")),
        ),
        malformed: fx.plain(
            &format!("mal_{tag}.vcf"),
            &format!("{head}\n{}\n", malformed.join("\n")),
        ),
        trunc_gz: s_of(&write(
            fx.dir.path(),
            &format!("t_{tag}.vcf.gz"),
            &gz[..gz.len() / 2],
        )),
        nofooter_gz: s_of(&write(
            fx.dir.path(),
            &format!("nf_{tag}.vcf.gz"),
            &gz[..gz.len() - 28],
        )),
        trunc_bcf: s_of(&write(
            fx.dir.path(),
            &format!("t_{tag}.bcf"),
            &bcf[..bcf.len() / 2],
        )),
    }
}
