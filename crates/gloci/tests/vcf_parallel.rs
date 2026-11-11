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
    cohort_with(samples, records, bad, true)
}

/// [`cohort`], with `haploid` calls (a lone allele) or without: `--het`
/// skips every site that has a haploid sample.
fn cohort_with(samples: usize, records: usize, bad: &[usize], haploid: bool) -> String {
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
                1 if haploid => format!("{a}"),
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

fn lines_of(b: &[u8], prefix_skip: &str) -> usize {
    String::from_utf8_lossy(b)
        .lines()
        .filter(|l| !l.starts_with(prefix_skip))
        .count()
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
    let with_sel: Vec<&str> = vec![
        "--samples",
        &names,
        "--min-dp",
        "6",
        "--include",
        "INFO/DP>9",
        "--thin",
        "50",
        "--region",
        "chr1:1-9000",
        "--region",
        "chr2:1-9000",
    ];
    for (sel_name, sel) in [("sel", &with_sel), ("none", &Vec::new())] {
        for input in [&vcf, &bcf] {
            let tag = format!("{sel_name}-{}", if input == &vcf { "vcf" } else { "bcf" });
            let ext = if input == &vcf { "vcf" } else { "bcf" };
            let mut filtered = Vec::new();
            let mut extracted = Vec::new();
            let mut stats = Vec::new();
            for t in ["1", "2", "8"] {
                // The recorded command line names the output, so each thread
                // count writes the same relative names in its own directory.
                let dir = d.path().join(format!("{tag}-t{t}"));
                std::fs::create_dir_all(&dir).unwrap();
                let out = format!("out.{ext}");
                let mut a = vec!["vcf", "filter", input, "-o", &out, "--threads", t];
                a.extend(sel);
                run_in(&dir, &a);
                filtered.push(read(&s(&dir.join(&out))));
                let tab = "out.tsv".to_string();
                let mut a = vec![
                    "vcf",
                    "extract",
                    input,
                    "--format",
                    "GT,DP,GQ",
                    "-o",
                    &tab,
                    "--threads",
                    t,
                ];
                a.extend(sel);
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
                a.extend(sel);
                run_in(&dir, &a);
                let rd = |f: &str| read(&s(&dir.join(f)));
                extracted.push((
                    rd(&tab),
                    rd(&format!("{m}.012")),
                    rd(&format!("{m}.012.pos")),
                ));
                let p = "st".to_string();
                let mut a = vec!["vcf", "stats", input, "--out-prefix", &p, "--threads", t];
                a.extend(STATS);
                a.extend(["--summary", "--json", "--min-gq", "10"]);
                a.extend(sel);
                run_in(&dir, &a);
                stats.push(s(&dir.join(&p)));
            }
            // The runs must have something to compare.
            if ext == "vcf" {
                assert!(
                    lines_of(&filtered[0], "#") > 100,
                    "{tag}: filter kept little"
                );
            } else {
                assert!(filtered[0].len() > 10_000, "{tag}: filter kept little");
            }
            assert!(lines_of(&extracted[0].0, "#") > 100, "{tag}: extract table");
            assert!(lines_of(&extracted[0].2, "#") > 100, "{tag}: .012.pos");
            assert!(
                lines_of(&read(&format!("{}.lmiss", stats[0])), "CHR") > 100,
                "{tag}: .lmiss rows"
            );
            for i in 1..3 {
                same(&filtered[0], &filtered[i], &format!("{tag} filter"));
                assert_eq!(extracted[0], extracted[i], "{tag} extract");
                for e in STAT_EXT.iter().chain(&["summary.txt", "summary.json"]) {
                    same(
                        &read(&format!("{}.{e}", stats[0])),
                        &read(&format!("{}.{e}", stats[i])),
                        &format!("{tag} stats .{e}"),
                    );
                }
            }
        }
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

fn in_process_ro() -> ReaderOptions {
    ReaderOptions {
        accept: Accept::default(),
        chrom_map: None,
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
        threads: NonZero::new(8).unwrap(),
    }
}

#[test]
fn parallel_selection_bounded() {
    let _g = IN_PROCESS.lock().unwrap_or_else(|e| e.into_inner());
    let d = tempdir();
    let samples = 20_000;
    // 2,000 records x 20,000 samples = 40M genotypes, five times the bound.
    let records = 2_000;
    let mut t = String::from(HEAD);
    t.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..samples {
        write!(t, "\tS{i}").unwrap();
    }
    t.push('\n');
    for r in 1..=records {
        write!(t, "chr1\t{r}\t.\tA\tC\t50\tPASS\t.\tGT").unwrap();
        for i in 0..samples {
            t.push_str(["\t0/0", "\t0/1", "\t1/1"][(i + r) % 3]);
        }
        t.push('\n');
    }
    let path = s(&plain_file(d.path(), "wide.vcf", &t));
    drop(t);
    let sel = VcfSelect::default();
    let input = open_selected(&path, &in_process_ro(), &sel).unwrap();
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
        // A stalled consumer lets the workers run ahead as far as the stage
        // allows: wait until a full complement of batches is in flight (13
        // of 26 x 20,000), then give a stage without a bound time to exceed it.
        if seen == 1 {
            let start = std::time::Instant::now();
            while inflight_peak() < 13 * 26 * 20_000 && start.elapsed().as_secs() < 120 {
                std::thread::sleep(std::time::Duration::from_millis(100));
            }
            std::thread::sleep(std::time::Duration::from_secs(3));
        }
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
    let batch = 26 * samples; // max(1, 524288 / 20000) records
    let bound = 2 * workers * 524_288;
    assert!(
        peak > 5 * batch,
        "{peak} genotypes in flight: no read-ahead"
    );
    assert!(
        peak <= bound,
        "{peak} genotypes in flight; the bound is {bound}"
    );
    assert!(
        records * samples > 3 * bound,
        "the input must exceed the bound"
    );
}

/// The in-process tests share the process-wide test hooks and counters.
static IN_PROCESS: std::sync::Mutex<()> = std::sync::Mutex::new(());

#[test]
fn parallel_selection_reports_a_worker_panic() {
    let _g = IN_PROCESS.lock().unwrap_or_else(|e| e.into_inner());
    let d = tempdir();
    // Several batches (about 150 records each), so the driver is still
    // dispatching when the worker dies; at -t 2 the only worker is gone.
    let path = s(&plain_file(d.path(), "p.vcf", &cohort(2000, 800, &[])));
    for (threads, seq) in [(2usize, 0u64), (2, 1), (4, 1)] {
        let (tx, rx) = std::sync::mpsc::channel();
        let path = path.clone();
        std::thread::spawn(move || {
            let sel = VcfSelect::default();
            let input = open_selected(&path, &in_process_ro(), &sel).unwrap();
            let (reader, mut chroms) = input.into_parts();
            let mut stage = gloci_ops::SelectStage::new(
                "test",
                reader.as_ref(),
                &sel,
                &RewritePlan::default(),
                true,
            )
            .unwrap();
            stage.set_threads(NonZero::new(threads).unwrap());
            gloci_ops::set_worker_panic_at(Some(seq));
            let r = run_selection(reader, &mut chroms, &mut stage, |_| Ok(()));
            gloci_ops::set_worker_panic_at(None);
            let _ = tx.send(r.map_err(|e| e.to_string()));
        });
        let r = rx
            .recv_timeout(std::time::Duration::from_secs(120))
            .expect("the selection hung after a worker panic");
        let e = r.expect_err("a worker panic must be an error");
        assert!(
            e.contains("a selection worker panicked on the batch starting at record"),
            "-t {threads}, batch {seq}: {e}"
        );
    }
}

/// Both errors lie in one batch: the earlier record's wins at every -t.
#[test]
fn parallel_selection_earliest_record_error_in_a_batch() {
    let d = tempdir();
    let text = cohort(40, 1200, &[]);
    // Record 400: sample 0 names allele 9 of a site with at most 2 ALTs.
    // Record 450: out of order (position 1, after larger ones).
    let mut out = String::new();
    let mut n = 0;
    for l in text.lines() {
        n += if l.starts_with('#') { 0 } else { 1 };
        let mut l = l.to_string();
        if n == 400 {
            let cols: Vec<&str> = l.split('\t').collect();
            let mut c: Vec<String> = cols.iter().map(|x| x.to_string()).collect();
            c[9] = "9/9:3:3:PASS".into();
            l = c.join("\t");
        }
        if n == 450 {
            let mut c: Vec<String> = l.split('\t').map(|x| x.to_string()).collect();
            c[1] = "1".into();
            l = c.join("\t");
        }
        out.push_str(&l);
        out.push('\n');
    }
    let vcf = s(&plain_file(d.path(), "e.vcf", &out));
    let bcf = s(&d.path().join("e.bcf"));
    run(&["vcf", "filter", &vcf, "-o", &bcf]);
    for input in [&vcf, &bcf] {
        let mut msgs = Vec::new();
        for t in ["1", "8"] {
            let o = gloci()
                .args([
                    "vcf",
                    "filter",
                    input,
                    "--min-dp",
                    "1",
                    "--thin",
                    "50",
                    "--threads",
                    t,
                ])
                .assert()
                .failure();
            msgs.push(String::from_utf8(o.get_output().stderr.clone()).unwrap());
        }
        for m in &msgs {
            assert!(
                m.contains("record 400") && m.contains("allele 9"),
                "{input}: {m}"
            );
        }
    }
}

#[test]
fn worker_tallies_equal_sequential() {
    // Diploid data with missing and masked genotypes, so --het has sites;
    // --depth, --missing-indv and --het run on the workers' counters, apart
    // from the --thin run, which keeps them on the driver.
    let d = tempdir();
    let (samples, records) = (300, 1500);
    let vcf = s(&plain_file(
        d.path(),
        "t.vcf",
        &cohort_with(samples, records, &[], false),
    ));
    let names: Vec<String> = (0..samples)
        .filter(|i| i % 5 != 0)
        .map(|i| format!("S{i}"))
        .collect();
    let names = names.join(",");
    let runs: [(&str, Vec<&str>); 4] = [
        ("plain", vec![]),
        (
            "masked",
            vec![
                "--min-dp",
                "8",
                "--min-gq",
                "20",
                "--mask-filtered-genotypes",
            ],
        ),
        ("samples", vec!["--samples", &names, "--max-missing", "0.5"]),
        (
            "thin",
            vec!["--thin", "40", "--region", "chr1:1-20000", "--min-dp", "5"],
        ),
    ];
    // With --geno-depth the records carry their samples; without it nothing
    // reads them and the workers pack none (the bench's path).
    let sets: [(&str, &[&str], &[&str]); 2] = [
        (
            "gd",
            &[
                "--depth",
                "--missing-indv",
                "--het",
                "--freq",
                "--geno-depth",
            ],
            &["idepth", "imiss", "het", "frq", "gdepth"],
        ),
        (
            "nogd",
            &["--depth", "--missing-indv", "--het", "--freq"],
            &["idepth", "imiss", "het", "frq"],
        ),
    ];
    for ((tag, sel), (set, flags, exts)) in
        runs.iter().flat_map(|r| sets.iter().map(move |x| (r, x)))
    {
        let tag = format!("{tag}-{set}");
        let mut outs = Vec::new();
        for t in ["1", "2", "8"] {
            let dir = d.path().join(format!("{tag}-t{t}"));
            std::fs::create_dir_all(&dir).unwrap();
            let mut a = vec!["vcf", "stats", &vcf, "--out-prefix", "o", "--threads", t];
            a.extend(flags.iter());
            a.extend(sel);
            run_in(&dir, &a);
            outs.push(dir);
        }
        for e in exts.iter() {
            let base = read(&s(&outs[0].join(format!("o.{e}"))));
            for o in &outs[1..] {
                same(
                    &base,
                    &read(&s(&o.join(format!("o.{e}")))),
                    &format!("{tag} .{e}"),
                );
            }
        }
        let het = read(&s(&outs[0].join("o.het")));
        assert!(lines_of(&het, "INDV") > 50, "{tag}: --het had no rows");
        let imiss = read(&s(&outs[0].join("o.imiss")));
        assert!(
            lines_of(&imiss, "INDV") > 50,
            "{tag}: --missing-indv had no rows"
        );
    }
}

/// `-t 8` is a cap on BGZF decompression threads plus selection workers
/// together, not on each (M9): polled from `/proc/PID/status`, the process has
/// at most 8 such threads besides the fixed ones, and the output is the `-t 1` output.
#[cfg(target_os = "linux")]
#[test]
fn threads_cap_decompression_plus_selection_workers() {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
    let d = tempdir();
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n",
    );
    for i in 0..300_000u64 {
        let _ = writeln!(
            text,
            "chr1\t{}\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\t1/1",
            100 + i * 3
        );
    }
    let f = bgzf_file(d.path(), "in.vcf.gz", &text)
        .to_string_lossy()
        .into_owned();
    let one = d.path().join("one.vcf");
    gloci()
        .args(["vcf", "filter", &f, "-t", "1", "-O", "vcf", "-o"])
        .arg(&one)
        .assert()
        .success();
    let eight = d.path().join("eight.vcf");
    let mut child = std::process::Command::new(assert_cmd::cargo::cargo_bin("gloci"))
        .args(["vcf", "filter", &f, "-t", "8", "-O", "vcf", "-o"])
        .arg(&eight)
        .spawn()
        .unwrap();
    let pid = child.id();
    let peak = Arc::new(AtomicU64::new(0));
    let stop = Arc::new(AtomicBool::new(false));
    let poller = {
        let (peak, stop) = (peak.clone(), stop.clone());
        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                if let Ok(s) = std::fs::read_to_string(format!("/proc/{pid}/status"))
                    && let Some(l) = s.lines().find_map(|l| l.strip_prefix("Threads:"))
                    && let Ok(t) = l.trim().parse::<u64>()
                {
                    peak.fetch_max(t, Ordering::Relaxed);
                }
                std::thread::sleep(std::time::Duration::from_micros(200));
            }
        })
    };
    assert!(child.wait().unwrap().success());
    stop.store(true, Ordering::Relaxed);
    poller.join().unwrap();
    let peak = peak.load(Ordering::Relaxed);
    // 3 fixed (main, signal handler, the command's thread), 8 for the idle
    // global rayon pool, and at most 8 for decompression plus selection
    // workers; one spare. Before the cap this was 27.
    assert!(peak <= 3 + 8 + 8 + 1, "peak {peak} threads at -t 8");
    // Records are the same (the ##glociCommand line differs by -t).
    let body = |p: &Path| {
        std::fs::read_to_string(p)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with("##glociCommand"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&one), body(&eight));
}
