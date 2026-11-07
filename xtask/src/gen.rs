//! Deterministic synthetic BED and VCF generators for tests and benchmarks.
//!
//! Output is written as it is generated. Unsorted BED output holds only the
//! current record in memory; `--sorted` holds one chromosome at a time.
use anyhow::{Result, bail, ensure};
use clap::ValueEnum;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::io::Write;

/// Every synthetic chromosome has this length.
pub const CHROM_LEN: u64 = 250_000_000;

/// Shape of the generated BED intervals.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Profile {
    /// Random positions, log-uniform lengths.
    Uniform,
    /// Many intervals centred on a few shared points, so they nest deeply.
    DeepNesting,
    /// Intervals spanning almost the whole chromosome.
    WholeChrom,
    /// Zero-length intervals (start == end).
    Points,
}

#[derive(Clone, Debug)]
pub struct BedOpts {
    pub records: u64,
    pub chroms: u32,
    pub seed: u64,
    pub profile: Profile,
    pub len_min: u64,
    pub len_max: u64,
    pub sorted: bool,
}

impl Default for BedOpts {
    fn default() -> Self {
        BedOpts {
            records: 1000,
            chroms: 24,
            seed: 1,
            profile: Profile::Uniform,
            len_min: 50,
            len_max: 10_000,
            sorted: false,
        }
    }
}

/// Number of shared centres per chromosome in the deep-nesting profile.
const CENTRES: u64 = 100;

struct Rec {
    start: u64,
    end: u64,
    score: u32,
    strand: u8,
}

fn log_uniform(rng: &mut ChaCha8Rng, min: u64, max: u64) -> u64 {
    if min >= max {
        return min;
    }
    let (a, b) = ((min.max(1) as f64).ln(), (max as f64).ln());
    let x = rng.random_range(a..=b).exp() as u64;
    x.clamp(min, max)
}

fn make_rec(rng: &mut ChaCha8Rng, o: &BedOpts) -> Rec {
    let (start, end) = match o.profile {
        Profile::Uniform => {
            let len = log_uniform(rng, o.len_min, o.len_max).min(CHROM_LEN);
            let s = rng.random_range(0..=CHROM_LEN - len);
            (s, s + len)
        }
        Profile::DeepNesting => {
            let k = rng.random_range(0..CENTRES);
            let centre = (k + 1) * (CHROM_LEN / (CENTRES + 1));
            let len = log_uniform(rng, o.len_min, o.len_max.saturating_mul(100))
                .min(CHROM_LEN / (CENTRES + 1));
            (centre - len / 2, centre - len / 2 + len)
        }
        Profile::WholeChrom => {
            let s = rng.random_range(0..1_000_000);
            let e = CHROM_LEN - rng.random_range(0..1_000_000);
            (s, e)
        }
        Profile::Points => {
            let p = rng.random_range(0..=CHROM_LEN);
            (p, p)
        }
    };
    Rec {
        start,
        end,
        score: rng.random_range(0..=1000),
        strand: if rng.random_bool(0.5) { b'+' } else { b'-' },
    }
}

fn write_rec(w: &mut impl Write, chrom: u32, r: &Rec, idx: u64) -> Result<()> {
    writeln!(
        w,
        "chr{chrom}\t{}\t{}\tr{idx}\t{}\t{}",
        r.start, r.end, r.score, r.strand as char
    )?;
    Ok(())
}

/// Write a BED6 file. Chromosomes are `chr1..chrK`.
pub fn bed(o: &BedOpts, out: &mut impl Write) -> Result<()> {
    ensure!(o.chroms >= 1, "--chroms must be at least 1");
    ensure!(
        o.len_min <= o.len_max,
        "--len-min must not exceed --len-max"
    );
    let mut rng = ChaCha8Rng::seed_from_u64(o.seed);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    if !o.sorted {
        for i in 0..o.records {
            let c = rng.random_range(1..=o.chroms);
            let r = make_rec(&mut rng, o);
            write_rec(&mut w, c, &r, i)?;
        }
        return Ok(w.flush()?);
    }
    // Sorted: chromosomes in byte order, one at a time.
    let mut names: Vec<u32> = (1..=o.chroms).collect();
    names.sort_by_key(|c| format!("chr{c}"));
    let k = u64::from(o.chroms);
    let mut idx = 0;
    for (n, c) in names.iter().enumerate() {
        let count = o.records / k + u64::from((n as u64) < o.records % k);
        let mut recs: Vec<Rec> = (0..count).map(|_| make_rec(&mut rng, o)).collect();
        recs.sort_by_key(|r| (r.start, r.end));
        for r in &recs {
            write_rec(&mut w, *c, r, idx)?;
            idx += 1;
        }
    }
    Ok(w.flush()?)
}

#[derive(Clone, Debug)]
pub struct VcfOpts {
    pub records: u64,
    pub samples: u32,
    pub seed: u64,
    pub chroms: u32,
}

const BASES: [u8; 4] = *b"ACGT";

fn base(rng: &mut ChaCha8Rng) -> u8 {
    BASES[rng.random_range(0..4)]
}

fn other_base(rng: &mut ChaCha8Rng, b: u8) -> u8 {
    loop {
        let x = base(rng);
        if x != b {
            return x;
        }
    }
}

fn seq(rng: &mut ChaCha8Rng, n: usize) -> String {
    (0..n).map(|_| base(rng) as char).collect()
}

/// Write a sorted VCF (chromosomes in header order, increasing POS) mixing
/// SNVs, MNVs, insertions, deletions, multi-allelic sites and `<DEL>` with END.
pub fn vcf(o: &VcfOpts, out: &mut impl Write) -> Result<()> {
    ensure!(o.chroms >= 1, "--chroms must be at least 1");
    let mut rng = ChaCha8Rng::seed_from_u64(o.seed);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    writeln!(w, "##fileformat=VCFv4.3")?;
    for c in 1..=o.chroms {
        writeln!(w, "##contig=<ID=chr{c},length={CHROM_LEN}>")?;
    }
    writeln!(w, "##ALT=<ID=DEL,Description=\"Deletion\">")?;
    writeln!(
        w,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for s in 1..=o.samples {
        write!(w, "\tS{s}")?;
    }
    writeln!(w)?;
    let k = u64::from(o.chroms);
    for c in 1..=o.chroms {
        let n = o.records / k + u64::from(c == o.chroms) * (o.records % k);
        let max_gap = (CHROM_LEN / n.max(1)).saturating_sub(1).max(1);
        let mut pos = 0u64;
        for _ in 0..n {
            pos += rng.random_range(1..=max_gap);
            let rb = base(&mut rng);
            let mut info = String::from(".");
            let (r, alt, multi) = match rng.random_range(0..100) {
                0..40 => (
                    (rb as char).to_string(),
                    (other_base(&mut rng, rb) as char).to_string(),
                    false,
                ),
                40..50 => {
                    let r = {
                        let n = rng.random_range(2..=3);
                        seq(&mut rng, n)
                    };
                    let a: String = r.bytes().map(|b| other_base(&mut rng, b) as char).collect();
                    (r, a, false)
                }
                50..65 => {
                    let ins = {
                        let n = rng.random_range(1..=5);
                        seq(&mut rng, n)
                    };
                    let r = (rb as char).to_string();
                    (r.clone(), format!("{r}{ins}"), false)
                }
                65..80 => {
                    let del = {
                        let n = rng.random_range(1..=5);
                        seq(&mut rng, n)
                    };
                    let r = (rb as char).to_string();
                    (format!("{r}{del}"), r, false)
                }
                80..92 => {
                    let mut alts = vec![other_base(&mut rng, rb)];
                    loop {
                        let b = other_base(&mut rng, rb);
                        if !alts.contains(&b) {
                            alts.push(b);
                            break;
                        }
                    }
                    let a: Vec<String> = alts.iter().map(|b| (*b as char).to_string()).collect();
                    ((rb as char).to_string(), a.join(","), true)
                }
                _ => {
                    let len = log_uniform(&mut rng, 50, 100_000);
                    info = format!("SVTYPE=DEL;END={}", (pos + len).min(CHROM_LEN));
                    ((rb as char).to_string(), "<DEL>".to_string(), false)
                }
            };
            if info == "." {
                info = format!("DP={}", rng.random_range(5..200));
            }
            let qual = rng.random_range(10..60);
            write!(w, "chr{c}\t{pos}\t.\t{r}\t{alt}\t{qual}\tPASS\t{info}\tGT")?;
            let max_gt = if multi { 2 } else { 1 };
            for _ in 0..o.samples {
                if rng.random_range(0..50) == 0 {
                    write!(w, "\t./.")?;
                } else {
                    let (a, b) = (rng.random_range(0..=max_gt), rng.random_range(0..=max_gt));
                    write!(w, "\t{}/{}", a.min(b), a.max(b))?;
                }
            }
            writeln!(w)?;
        }
    }
    Ok(w.flush()?)
}

/// Parse the generator arguments of a `cases.toml` `setup` string or the
/// command line: `bed --records N ...` or `vcf --records N ...`.
pub fn run_args(args: &[String], out: &mut impl Write) -> Result<()> {
    let (kind, rest) = match args.split_first() {
        Some((k, r)) => (k.as_str(), r),
        None => bail!("gen needs a kind: bed or vcf"),
    };
    let mut bed_o = BedOpts::default();
    let mut vcf_o = VcfOpts {
        records: 1000,
        samples: 3,
        seed: 1,
        chroms: 5,
    };
    let mut it = rest.iter();
    while let Some(flag) = it.next() {
        let mut val = || -> Result<&String> {
            it.next()
                .ok_or_else(|| anyhow::anyhow!("{flag} needs a value"))
        };
        match (kind, flag.as_str()) {
            ("bed", "--records") => bed_o.records = val()?.parse()?,
            ("bed", "--chroms") => bed_o.chroms = val()?.parse()?,
            ("bed", "--seed") => bed_o.seed = val()?.parse()?,
            ("bed", "--len-min") => bed_o.len_min = val()?.parse()?,
            ("bed", "--len-max") => bed_o.len_max = val()?.parse()?,
            ("bed", "--sorted") => bed_o.sorted = true,
            ("bed", "--profile") => {
                let v = val()?;
                bed_o.profile = Profile::from_str(v, true).map_err(|e| anyhow::anyhow!("{e}"))?
            }
            ("vcf", "--records") => vcf_o.records = val()?.parse()?,
            ("vcf", "--samples") => vcf_o.samples = val()?.parse()?,
            ("vcf", "--seed") => vcf_o.seed = val()?.parse()?,
            ("vcf", "--chroms") => vcf_o.chroms = val()?.parse()?,
            _ => bail!("unknown gen argument {flag} for {kind}"),
        }
    }
    match kind {
        "bed" => bed(&bed_o, out),
        "vcf" => vcf(&vcf_o, out),
        _ => bail!("unknown gen kind {kind} (bed or vcf)"),
    }
}
