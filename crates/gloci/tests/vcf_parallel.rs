//! The parallel selection stage (D157): every `vcf` command gives the same
//! bytes and the same error at every thread count, in bounded memory.
mod common;
use common::*;
use gloci_io::{Accept, InvalidPolicy, ReaderOptions};
use gloci_ops::{VcfSelect, inflight_peak, open_selected, reset_inflight_peak, run_selection};
use gloci_vcf::{GtFields, RewritePlan};
use std::fmt::Write as _;
use std::num::NonZero;
use std::path::Path;

const STATS: [&str; 15] = [
    "--freq",
    "--counts",
    "--depth",
    "--site-depth",
    "--site-mean-depth",
    "--geno-depth",
    "--site-quality",
    "--missing-site",
    "--missing-indv",
    "--het",
    "--hardy",
    "--singletons",
    "--filter-summary",
    "--tstv-summary",
    "--indel-hist",
];

const STAT_EXT: [&str; 15] = [
    "frq",
    "frq.count",
    "idepth",
    "ldepth",
    "ldepth.mean",
    "gdepth",
    "lqual",
    "lmiss",
    "imiss",
    "het",
    "hwe",
    "singletons",
    "FILTER.summary",
    "TsTv.summary",
    "indel.hist",
];

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn read(p: &str) -> Vec<u8> {
    std::fs::read(p).unwrap_or_else(|e| panic!("{p}: {e}"))
}

/// Assert two outputs are equal, naming the first line that differs.
fn same(a: &[u8], b: &[u8], what: &str) {
    if a != b {
        let (ta, tb) = (String::from_utf8_lossy(a), String::from_utf8_lossy(b));
        let at = ta.lines().zip(tb.lines()).position(|(x, y)| x != y);
        match at {
            Some(i) => panic!(
                "{what}: line {} differs:\n{}\n{}",
                i + 1,
                ta.lines()
                    .nth(i)
                    .unwrap_or("")
                    .chars()
                    .take(300)
                    .collect::<String>(),
                tb.lines()
                    .nth(i)
                    .unwrap_or("")
                    .chars()
                    .take(300)
                    .collect::<String>()
            ),
            None => panic!("{what}: {} bytes against {}", a.len(), b.len()),
        }
    }
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % n
    }
}

const HEAD: &str = "##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
##FILTER=<ID=LowQ,Description=\"Low quality\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">
";

/// A seeded cohort: `samples` samples, `records` records over two chromosomes.
/// `bad` lists record numbers (1-based) whose first sample has a broken GT.
fn cohort(samples: usize, records: usize, bad: &[usize]) -> String {
    let mut rng = Lcg(12345);
    let mut t = String::from(HEAD);
    t.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..samples {
        write!(t, "\tS{i}").unwrap();
    }
    t.push('\n');
    let mut pos = 0u64;
    for r in 1..=records {
        let chrom = if r <= records / 2 { "chr1" } else { "chr2" };
        if r == records / 2 + 1 {
            pos = 0;
        }
        pos += 5 + rng.next(40);
        let (refb, alt) = match rng.next(6) {
            0 => ("A", "AT"),
            1 => ("C", "T,G"),
            2 => ("G", "A"),
            3 => ("AT", "A"),
            _ => ("T", "C"),
        };
        let n_alt = alt.split(',').count() as u64;
        let filter = if rng.next(5) == 0 { "LowQ" } else { "PASS" };
        write!(
            t,
            "{chrom}\t{pos}\trs{r}\t{refb}\t{alt}\t{}\t{filter}\tDP={}\tGT:DP:GQ:FT",
            20 + rng.next(40),
            5 + rng.next(30)
        )
        .unwrap();
        for i in 0..samples {
            if i == 0 && bad.contains(&r) {
                t.push_str("\t0/x:3:3:PASS");
                continue;
            }
            let a = rng.next(n_alt + 1);
            let b = rng.next(n_alt + 1);
            let sep = if rng.next(3) == 0 { '|' } else { '/' };
            let gt = match rng.next(12) {
                0 => "./.".to_string(),
                1 => format!("{a}"),
                _ => format!("{a}{sep}{b}"),
            };
            let ft = if rng.next(9) == 0 { "LowQ" } else { "PASS" };
            write!(t, "\t{gt}:{}:{}:{ft}", 2 + rng.next(30), 5 + rng.next(60)).unwrap();
        }
        t.push('\n');
    }
    t
}

fn run(args: &[&str]) {
    gloci().args(args).assert().success();
}

fn run_in(dir: &Path, args: &[&str]) {
    gloci().current_dir(dir).args(args).assert().success();
}

#[test]
fn parallel_selection_equals_sequential() {
    let d = tempdir();
    let (samples, records) = (600, 1200);
    let vcf = s(&plain_file(
        d.path(),
        "c.vcf",
        &cohort(samples, records, &[]),
    ));
    let bcf = s(&d.path().join("c.bcf"));
    run(&["vcf", "filter", &vcf, "-o", &bcf]);
    let names: Vec<String> = (0..samples).step_by(2).map(|i| format!("S{i}")).collect();
    let names = names.join(",");
    let sel: Vec<&str> = vec![
        "--samples",
        &names,
        "--min-dp",
        "6",
        "--include",
        "INFO/DP>9",
        "--thin",
        "50",
        "--region",
        "chr1:1-9000,chr2:1-9000",
    ];
    for input in [&vcf, &bcf] {
        let tag = if input == &vcf { "vcf" } else { "bcf" };
        let ext = if input == &vcf { "vcf" } else { "bcf" };
        let mut filtered = Vec::new();
        let mut extracted = Vec::new();
        let mut stats = Vec::new();
        for t in ["1", "2", "8"] {
            // The recorded command line names the output, so each thread
            // count writes the same relative names in its own directory.
            let dir = d.path().join(format!("{tag}-t{t}"));
            std::fs::create_dir_all(&dir).unwrap();
            // filter
            let out = format!("out.{ext}");
            let mut a = vec!["vcf", "filter", input, "-o", &out, "--threads", t];
            a.extend(&sel);
            run_in(&dir, &a);
            filtered.push(read(&s(&dir.join(&out))));
            // extract: a table and the 012 matrix
            let tab = "out.tsv".to_string();
            let mut a = vec![
                "vcf",
                "extract",
                input,
                "--format",
                "GT,DP,GQ",
                "--info",
                "DP",
                "-o",
                &tab,
                "--threads",
                t,
            ];
            a.extend(&sel);
            // --info and --format together are one table at a time
            a.drain(5..7);
            run_in(&dir, &a);
            let m = "m".to_string();
            let mut a = vec![
                "vcf",
                "extract",
                input,
                "--genotype-matrix",
                "012",
                "--out-prefix",
                &m,
                "--threads",
                t,
            ];
            a.extend(&sel);
            run_in(&dir, &a);
            let rd = |f: &str| read(&s(&dir.join(f)));
            extracted.push((
                rd(&tab),
                rd(&format!("{m}.012")),
                rd(&format!("{m}.012.pos")),
            ));
            // stats: all 15, plus the summaries
            let p = "st".to_string();
            let mut a = vec!["vcf", "stats", input, "--out-prefix", &p, "--threads", t];
            a.extend(STATS);
            a.extend(["--summary", "--json", "--min-gq", "10"]);
            a.extend(&sel);
            run_in(&dir, &a);
            stats.push(s(&dir.join(&p)));
        }
        assert!(
            filtered[0].len() > 1000,
            "{tag}: the selection kept nothing"
        );
        for i in 1..3 {
            same(&filtered[0], &filtered[i], &format!("{tag} filter"));
            assert_eq!(extracted[0], extracted[i], "{tag} extract");
            for e in STAT_EXT.iter().chain(&["summary.txt", "summary.json"]) {
                assert_eq!(
                    read(&format!("{}.{e}", stats[0])),
                    read(&format!("{}.{e}", stats[i])),
                    "{tag} stats .{e}"
                );
            }
        }
        assert!(extracted[0].0.len() > 100 && extracted[0].1.len() > 100);
    }
}

#[test]
fn parallel_selection_first_error() {
    let d = tempdir();
    // Two broken records; the earlier one is the one to report.
    let vcf = s(&plain_file(
        d.path(),
        "bad.vcf",
        &cohort(600, 1200, &[800, 1100]),
    ));
    let mut msgs = Vec::new();
    for t in ["1", "8"] {
        let out = d.path().join(format!("out{t}.vcf"));
        let o = gloci()
            .args([
                "vcf",
                "filter",
                &vcf,
                "--min-dp",
                "1",
                "-o",
                &s(&out),
                "--threads",
                t,
            ])
            .assert()
            .failure();
        assert!(!out.exists(), "-t {t}: an output file was published");
        msgs.push(String::from_utf8(o.get_output().stderr.clone()).unwrap());
        let prefix = s(&d.path().join(format!("st{t}")));
        let o = gloci()
            .args([
                "vcf",
                "stats",
                &vcf,
                "--out-prefix",
                &prefix,
                "--freq",
                "--threads",
                t,
            ])
            .assert()
            .failure();
        msgs.push(String::from_utf8(o.get_output().stderr.clone()).unwrap());
        assert!(
            !Path::new(&format!("{prefix}.frq")).exists(),
            "-t {t}: a statistics file was published"
        );
    }
    assert!(msgs[0].contains("record 800"), "{}", msgs[0]);
    assert_eq!(msgs[0], msgs[2], "filter");
    assert_eq!(msgs[1], msgs[3], "stats");
}

#[test]
fn parallel_selection_bounded() {
    let d = tempdir();
    let samples = 20_000;
    let mut t = String::from(HEAD);
    t.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..samples {
        write!(t, "\tS{i}").unwrap();
    }
    t.push('\n');
    let records = 200;
    for r in 1..=records {
        write!(t, "chr1\t{r}\t.\tA\tC\t50\tPASS\t.\tGT").unwrap();
        for i in 0..samples {
            t.push_str(["\t0/0", "\t0/1", "\t1/1"][(i + r) % 3]);
        }
        t.push('\n');
    }
    let path = s(&plain_file(d.path(), "wide.vcf", &t));
    let ro = ReaderOptions {
        accept: Accept::default(),
        chrom_map: None,
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
        threads: NonZero::new(8).unwrap(),
    };
    let sel = VcfSelect::default();
    let input = open_selected(&path, &ro, &sel).unwrap();
    let (reader, mut chroms) = input.into_parts();
    let mut stage =
        gloci_ops::SelectStage::new("test", reader.as_ref(), &sel, &RewritePlan::default(), true)
            .unwrap();
    stage.share_decoded(GtFields::GT_ONLY);
    stage.set_threads(NonZero::new(8).unwrap());
    reset_inflight_peak();
    let mut seen = 0u64;
    let mut genotypes = 0u64;
    let n = run_selection(reader, &mut chroms, &mut stage, |hit| {
        seen += 1;
        assert_eq!(hit.rec_no, seen);
        genotypes += hit.genotypes().unwrap().n_samples() as u64;
        Ok(())
    })
    .unwrap();
    assert_eq!(n, records as u64);
    assert_eq!(seen, records as u64);
    assert_eq!(genotypes, records as u64 * samples as u64);
    let peak = inflight_peak();
    let workers = 7;
    assert!(peak > 0, "nothing was held in flight");
    assert!(
        peak <= 2 * workers * 524_288,
        "{peak} genotypes in flight; the bound is {}",
        2 * workers * 524_288
    );
}
