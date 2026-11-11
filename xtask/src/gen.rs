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
pub struct PairOpts {
    pub records: u64,
    pub depth: f64,
    pub chroms: u32,
    pub seed: u64,
    pub len_min: u64,
    pub len_max: u64,
    pub sorted: bool,
}

impl Default for PairOpts {
    fn default() -> Self {
        PairOpts {
            records: 1000,
            depth: 1.0,
            chroms: 24,
            seed: 1,
            len_min: 50,
            len_max: 10_000,
            sorted: true,
        }
    }
}

/// Write A/B BED pairs with controlled overlap depth.
pub fn pair(o: &PairOpts, out_a: &mut impl Write, out_b: &mut impl Write) -> Result<()> {
    ensure!(o.chroms >= 1, "--chroms must be at least 1");
    ensure!(
        o.len_min <= o.len_max,
        "--len-min must not exceed --len-max"
    );
    ensure!(o.depth >= 0.0, "--depth must be non-negative");

    let mut rng = ChaCha8Rng::seed_from_u64(o.seed);
    let mut wa = std::io::BufWriter::with_capacity(1 << 20, out_a);
    let mut wb = std::io::BufWriter::with_capacity(1 << 20, out_b);

    let k = u64::from(o.chroms);
    let mut names: Vec<u32> = (1..=o.chroms).collect();
    if o.sorted {
        names.sort_by_key(|c| format!("chr{c}"));
    }

    let mut idx_a = 0u64;
    let mut idx_b = 0u64;

    for (n, c) in names.iter().enumerate() {
        let n_records = o.records / k + u64::from((n as u64) < o.records % k);
        if n_records == 0 {
            continue;
        }

        let mut recs_a = Vec::with_capacity(n_records as usize);
        let mut recs_b = Vec::with_capacity(n_records as usize);

        if o.depth == 0.0 {
            let stride = (CHROM_LEN as f64 - 100.0) / (n_records as f64);
            for i in 0..n_records {
                let slot_start = 50 + (i as f64 * stride).round() as u64;
                let slot_end = 50 + ((i + 1) as f64 * stride).round() as u64;
                let slot_end = slot_end.min(CHROM_LEN - 10);
                let mid = (slot_start + slot_end) / 2;

                let w_a = mid.saturating_sub(slot_start).max(2);
                let max_la = (w_a / 2).min(o.len_max).max(1);
                let min_la = o.len_min.min(max_la).max(1);
                let len_a = log_uniform(&mut rng, min_la, max_la);
                let s_a = slot_start + rng.random_range(0..=(w_a.saturating_sub(len_a)));
                let e_a = s_a + len_a;
                recs_a.push(Rec {
                    start: s_a,
                    end: e_a,
                    score: rng.random_range(0..=1000),
                    strand: if rng.random_bool(0.5) { b'+' } else { b'-' },
                });

                let w_b = slot_end.saturating_sub(mid).max(2);
                let max_lb = (w_b / 2).min(o.len_max).max(1);
                let min_lb = o.len_min.min(max_lb).max(1);
                let len_b = log_uniform(&mut rng, min_lb, max_lb);
                let s_b = mid + rng.random_range(0..=(w_b.saturating_sub(len_b)));
                let e_b = s_b + len_b;
                recs_b.push(Rec {
                    start: s_b,
                    end: e_b,
                    score: rng.random_range(0..=1000),
                    strand: if rng.random_bool(0.5) { b'+' } else { b'-' },
                });
            }
        } else {
            let margin = ((CHROM_LEN as f64 / (n_records as f64 + 1.0)) * (o.depth / 2.0)).round()
                as u64
                + 50;
            let margin = margin.min(CHROM_LEN / 4);
            let usable = CHROM_LEN.saturating_sub(margin * 2).max(1);
            let stride = usable as f64 / (n_records as f64);
            let l = (stride * (o.depth / 2.0)).max(1.0);

            for i in 0..n_records {
                let len_a = ((l * rng.random_range(0.8..1.2)).round() as u64).max(1);
                let c_a = margin as f64 + (i as f64 + rng.random_range(0.0..1.0)) * stride;
                let s_a = ((c_a - len_a as f64 / 2.0).round() as u64)
                    .max(1)
                    .min(CHROM_LEN - len_a - 1);
                let e_a = (s_a + len_a).min(CHROM_LEN);
                recs_a.push(Rec {
                    start: s_a,
                    end: e_a,
                    score: rng.random_range(0..=1000),
                    strand: if rng.random_bool(0.5) { b'+' } else { b'-' },
                });

                let len_b = ((l * rng.random_range(0.8..1.2)).round() as u64).max(1);
                let c_b = margin as f64 + (i as f64 + rng.random_range(0.0..1.0)) * stride;
                let s_b = ((c_b - len_b as f64 / 2.0).round() as u64)
                    .max(1)
                    .min(CHROM_LEN - len_b - 1);
                let e_b = (s_b + len_b).min(CHROM_LEN);
                recs_b.push(Rec {
                    start: s_b,
                    end: e_b,
                    score: rng.random_range(0..=1000),
                    strand: if rng.random_bool(0.5) { b'+' } else { b'-' },
                });
            }
        }

        if o.sorted {
            recs_a.sort_by_key(|r| (r.start, r.end));
            recs_b.sort_by_key(|r| (r.start, r.end));
        }

        for r in &recs_a {
            write_rec(&mut wa, *c, r, idx_a)?;
            idx_a += 1;
        }
        for r in &recs_b {
            write_rec(&mut wb, *c, r, idx_b)?;
            idx_b += 1;
        }
    }

    wa.flush()?;
    wb.flush()?;
    Ok(())
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
/// The header says VCFv4.2, the newest version vcftools 0.1.16 reads (nothing
/// here needs 4.3).
pub fn vcf(o: &VcfOpts, out: &mut impl Write) -> Result<()> {
    ensure!(o.chroms >= 1, "--chroms must be at least 1");
    let mut rng = ChaCha8Rng::seed_from_u64(o.seed);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    writeln!(w, "##fileformat=VCFv4.2")?;
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
    writeln!(w, "##FILTER=<ID=LowQual,Description=\"Low quality\">")?;
    writeln!(w, "##FILTER=<ID=q10,Description=\"Quality below 10\">")?;
    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype Filter\">"
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
            let qual_str = match rng.random_range(0..100) {
                0..5 => ".".to_string(),
                _ => rng.random_range(10..60).to_string(),
            };
            let filter_str = match rng.random_range(0..100) {
                0..75 => "PASS",
                75..85 => "LowQual",
                85..95 => "q10",
                _ => ".",
            };
            write!(
                w,
                "chr{c}\t{pos}\t.\t{r}\t{alt}\t{qual_str}\t{filter_str}\t{info}\tGT:DP:GQ:FT"
            )?;
            let max_gt = if multi { 2 } else { 1 };
            for _ in 0..o.samples {
                if rng.random_range(0..20) == 0 {
                    write!(w, "\t./.:.:.:.")?;
                } else {
                    let (a, b) = (rng.random_range(0..=max_gt), rng.random_range(0..=max_gt));
                    let dp = rng.random_range(5..100);
                    let gq = rng.random_range(10..99);
                    let ft = if dp < 10 { "LowDP" } else { "PASS" };
                    write!(w, "\t{}/{}:{dp}:{gq}:{ft}", a.min(b), a.max(b))?;
                }
            }
            writeln!(w)?;
        }
    }
    Ok(w.flush()?)
}

// ---- SAM and FASTA ----

/// A synthetic chromosome never shrinks below this, so a read plus its
/// longest intron always fits.
const MIN_SAM_CHROM: u64 = 50_000;
/// Longest intron (CIGAR N) and deletion (CIGAR D) of a generated read.
const MAX_INTRON: u32 = 5_000;
const MAX_DEL: u32 = 5;
/// Fraction of reads with a deletion (the rest, after the spliced ones, are
/// one plain match block).
const DEL_FRAC: f64 = 0.03;

#[derive(Clone, Debug)]
pub struct SamOpts {
    pub records: u64,
    pub read_len: u32,
    /// Mean depth of the aligned bases over the reference.
    pub depth: f64,
    pub chroms: u32,
    /// Fraction of reads with a CIGAR N (a spliced read).
    pub spliced_frac: f64,
    /// Records in random order instead of coordinate order.
    pub shuffled: bool,
    pub seed: u64,
    /// Salt of the reads' random stream, so several parts of one reference
    /// hold different reads. The reference depends on `seed` alone.
    pub part: u64,
}

impl Default for SamOpts {
    fn default() -> Self {
        SamOpts {
            records: 1000,
            read_len: 100,
            depth: 30.0,
            chroms: 24,
            spliced_frac: 0.1,
            shuffled: false,
            seed: 1,
            part: 0,
        }
    }
}

fn check_sam(o: &SamOpts) -> Result<()> {
    ensure!(o.chroms >= 1, "--chroms must be at least 1");
    ensure!(o.read_len >= 20, "--read-len must be at least 20");
    ensure!(o.read_len <= 10_000, "--read-len must be at most 10000");
    ensure!(o.depth > 0.0, "--depth must be positive");
    ensure!(
        (0.0..=1.0).contains(&o.spliced_frac),
        "--spliced-frac must be between 0 and 1"
    );
    Ok(())
}

/// Length of every chromosome of a SAM: `records * read_len / (depth * chroms)`.
pub fn sam_chrom_len(o: &SamOpts) -> u64 {
    let total = o.records as f64 * f64::from(o.read_len) / o.depth;
    ((total / f64::from(o.chroms)).ceil() as u64).max(MIN_SAM_CHROM)
}

/// The A/C/G/T sequence of chromosome `idx` (0-based) for `seed`.
fn ref_chrom(seed: u64, idx: u32, len: u64) -> Vec<u8> {
    let mut rng =
        ChaCha8Rng::seed_from_u64(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (u64::from(idx) + 1));
    let mut v = Vec::with_capacity(len as usize);
    while (v.len() as u64) < len {
        let mut x: u64 = rng.random();
        for _ in 0..32 {
            if (v.len() as u64) >= len {
                break;
            }
            v.push(BASES[(x & 3) as usize]);
            x >>= 2;
        }
    }
    v
}

fn write_fasta_seq(w: &mut impl Write, name: &str, seq: &[u8], width: usize) -> Result<()> {
    writeln!(w, ">{name}")?;
    for line in seq.chunks(width) {
        w.write_all(line)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// The reference FASTA the reads of `sam(o)` align to (`chr1..chrK`, 60
/// columns).
pub fn sam_reference(o: &SamOpts, out: &mut impl Write) -> Result<()> {
    check_sam(o)?;
    let len = sam_chrom_len(o);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    for c in 0..o.chroms {
        let seq = ref_chrom(o.seed, c, len);
        write_fasta_seq(&mut w, &format!("chr{}", c + 1), &seq, 60)?;
    }
    Ok(w.flush()?)
}

#[derive(Clone, Copy)]
struct Read {
    start: u32,
    /// First match block; the whole read when `gap` is 0.
    first: u32,
    /// Reference bases skipped after the first block.
    gap: u32,
    /// The gap is an intron (N), not a deletion (D).
    intron: bool,
    reverse: bool,
    mapq: u8,
}

fn make_read(rng: &mut ChaCha8Rng, o: &SamOpts, chrom_len: u64) -> Read {
    let u: f64 = rng.random();
    let (first, gap, intron) = if u < o.spliced_frac {
        (
            rng.random_range(10..=o.read_len - 10),
            rng.random_range(100..=MAX_INTRON),
            true,
        )
    } else if u < o.spliced_frac + DEL_FRAC {
        (
            rng.random_range(10..=o.read_len - 10),
            rng.random_range(1..=MAX_DEL),
            false,
        )
    } else {
        (o.read_len, 0, false)
    };
    let span = u64::from(o.read_len) + u64::from(gap);
    Read {
        start: rng.random_range(0..=chrom_len - span) as u32,
        first,
        gap,
        intron,
        reverse: rng.random_bool(0.5),
        mapq: if rng.random_bool(0.9) {
            60
        } else {
            rng.random_range(0..60)
        },
    }
}

fn write_read(
    w: &mut impl Write,
    line: &mut Vec<u8>,
    o: &SamOpts,
    chrom: u32,
    r: &Read,
    idx: u64,
    seq: &[u8],
) -> Result<()> {
    line.clear();
    let s = r.start as usize;
    let cigar = if r.gap == 0 {
        format!("{}M", o.read_len)
    } else {
        format!(
            "{}M{}{}{}M",
            r.first,
            r.gap,
            if r.intron { 'N' } else { 'D' },
            o.read_len - r.first
        )
    };
    write!(
        line,
        "r{idx}\t{}\tchr{chrom}\t{}\t{}\t{cigar}\t*\t0\t0\t",
        if r.reverse { 16 } else { 0 },
        s + 1,
        r.mapq
    )?;
    line.extend_from_slice(&seq[s..s + r.first as usize]);
    if r.gap > 0 {
        let rest = s + (r.first + r.gap) as usize;
        line.extend_from_slice(&seq[rest..rest + (o.read_len - r.first) as usize]);
    }
    line.push(b'\t');
    let salt = (idx % 21) as u32;
    for i in 0..o.read_len {
        line.push(33 + 20 + ((i * 7 + salt) % 21) as u8);
    }
    line.push(b'\n');
    w.write_all(line)?;
    Ok(())
}

/// Write single-end alignments as SAM text, to the reference of
/// [`sam_reference`]: every read matches it exactly at its position (its
/// bases skip the intron or deletion). Coordinate order by default (header
/// `SO:coordinate`), or `--shuffled` (`SO:unsorted`); either way the file is
/// produced a chromosome (sorted) or a record (shuffled) at a time.
pub fn sam(o: &SamOpts, out: &mut impl Write) -> Result<()> {
    check_sam(o)?;
    let len = sam_chrom_len(o);
    let refs: Vec<Vec<u8>> = (0..o.chroms).map(|c| ref_chrom(o.seed, c, len)).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(o.seed ^ o.part.wrapping_mul(0xD1B5_4A32_D192_ED03));
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    writeln!(
        w,
        "@HD\tVN:1.6\tSO:{}",
        if o.shuffled { "unsorted" } else { "coordinate" }
    )?;
    for c in 1..=o.chroms {
        writeln!(w, "@SQ\tSN:chr{c}\tLN:{len}")?;
    }
    let mut line = Vec::with_capacity(1024);
    if o.shuffled {
        for i in 0..o.records {
            let c = rng.random_range(0..o.chroms);
            let r = make_read(&mut rng, o, len);
            write_read(&mut w, &mut line, o, c + 1, &r, i, &refs[c as usize])?;
        }
        return Ok(w.flush()?);
    }
    let k = u64::from(o.chroms);
    let mut idx = 0;
    for c in 0..o.chroms {
        let count = o.records / k + u64::from(u64::from(c) < o.records % k);
        let mut reads: Vec<Read> = (0..count).map(|_| make_read(&mut rng, o, len)).collect();
        // samtools order: by position, a forward read before a reverse one.
        reads.sort_by_key(|r| (r.start, r.reverse));
        for r in &reads {
            write_read(&mut w, &mut line, o, c + 1, r, idx, &refs[c as usize])?;
            idx += 1;
        }
    }
    Ok(w.flush()?)
}

/// Write `count` sorted BED4 target windows (200 to 2000 bases) spread over
/// the chromosomes of `sam(o)`, in header order.
pub fn sam_targets(o: &SamOpts, count: u64, out: &mut impl Write) -> Result<()> {
    check_sam(o)?;
    let len = sam_chrom_len(o);
    let mut rng = ChaCha8Rng::seed_from_u64(o.seed ^ 0x007A_46E7);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    let k = u64::from(o.chroms);
    let mut idx = 0;
    for c in 0..o.chroms {
        let n = count / k + u64::from(u64::from(c) < count % k);
        if n == 0 {
            continue;
        }
        let stride = len / n;
        for i in 0..n {
            let wlen = rng.random_range(200..=2000).min(stride.max(2) / 2).max(1);
            let start = i * stride + rng.random_range(0..=(stride / 2));
            let end = (start + wlen).min(len);
            writeln!(w, "chr{}\t{start}\t{end}\tt{idx}", c + 1)?;
            idx += 1;
        }
    }
    Ok(w.flush()?)
}

#[derive(Clone, Debug)]
pub struct FastaOpts {
    /// Intervals of the BED written beside the FASTA.
    pub records: u64,
    pub chroms: u32,
    pub chrom_len: u64,
    pub line_width: usize,
    /// Fraction of the bases soft-masked (lowercase).
    pub soft_frac: f64,
    /// Fraction of the bases that are N.
    pub n_frac: f64,
    pub seed: u64,
    /// Bounds of the interval lengths of the BED.
    pub len_min: u64,
    pub len_max: u64,
}

impl Default for FastaOpts {
    fn default() -> Self {
        FastaOpts {
            records: 1000,
            chroms: 24,
            chrom_len: 1_000_000,
            line_width: 60,
            soft_frac: 0.05,
            n_frac: 0.001,
            seed: 1,
            len_min: 50,
            len_max: 10_000,
        }
    }
}

fn check_fasta(o: &FastaOpts) -> Result<()> {
    ensure!(o.chroms >= 1, "--chroms must be at least 1");
    ensure!(o.chrom_len >= 20_000, "--chrom-len must be at least 20000");
    ensure!(o.line_width >= 1, "--line-width must be at least 1");
    ensure!(
        1 <= o.len_min && o.len_min <= o.len_max && o.len_max <= o.chrom_len,
        "--len-min and --len-max must satisfy 1 <= min <= max <= --chrom-len"
    );
    for (f, n) in [(o.soft_frac, "--soft-frac"), (o.n_frac, "--n-frac")] {
        ensure!((0.0..1.0).contains(&f), "{n} must be in [0, 1)");
    }
    Ok(())
}

/// Apply runs of `run_min..=run_max` bases to `seq` covering about `frac` of
/// it, each run changed by `f`.
fn apply_runs(
    rng: &mut ChaCha8Rng,
    seq: &mut [u8],
    frac: f64,
    run_min: u64,
    run_max: u64,
    f: impl Fn(&mut u8),
) {
    if frac <= 0.0 {
        return;
    }
    let mean_run = (run_min + run_max) as f64 / 2.0;
    let mean_gap = (mean_run * (1.0 - frac) / frac).max(1.0);
    let mut pos = 0u64;
    let n = seq.len() as u64;
    loop {
        pos += rng.random_range(1..=(2.0 * mean_gap) as u64 + 1);
        if pos >= n {
            return;
        }
        let end = (pos + rng.random_range(run_min..=run_max)).min(n);
        seq[pos as usize..end as usize].iter_mut().for_each(&f);
        pos = end;
    }
}

/// Write a FASTA of `chr1..chrK`: random A/C/G/T with soft-masked and N runs.
/// Its A/C/G/T are those of [`sam_reference`] for the same seed and length.
pub fn fasta(o: &FastaOpts, out: &mut impl Write) -> Result<()> {
    check_fasta(o)?;
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    for c in 0..o.chroms {
        let mut seq = ref_chrom(o.seed, c, o.chrom_len);
        let mut rng = ChaCha8Rng::seed_from_u64(o.seed ^ (u64::from(c) + 1).wrapping_mul(0xA24B));
        apply_runs(&mut rng, &mut seq, o.soft_frac, 50, 2000, |b| {
            *b = b.to_ascii_lowercase()
        });
        apply_runs(&mut rng, &mut seq, o.n_frac, 100, 2000, |b| *b = b'N');
        write_fasta_seq(&mut w, &format!("chr{}", c + 1), &seq, o.line_width)?;
    }
    Ok(w.flush()?)
}

/// Write `records` BED6 intervals (`len_min` to `len_max` bases, log-uniform) inside the
/// chromosomes of `fasta(o)`, sorted in chromosome order.
pub fn fasta_intervals(o: &FastaOpts, out: &mut impl Write) -> Result<()> {
    check_fasta(o)?;
    let mut rng = ChaCha8Rng::seed_from_u64(o.seed ^ 0xB3D);
    let mut w = std::io::BufWriter::with_capacity(1 << 20, out);
    let k = u64::from(o.chroms);
    let mut idx = 0;
    for c in 0..o.chroms {
        let n = o.records / k + u64::from(u64::from(c) < o.records % k);
        let mut v: Vec<(u64, u64)> = (0..n)
            .map(|_| {
                let len = log_uniform(&mut rng, o.len_min, o.len_max);
                let s = rng.random_range(0..=o.chrom_len - len);
                (s, s + len)
            })
            .collect();
        v.sort_unstable();
        for (s, e) in v {
            let strand = if rng.random_bool(0.5) { '+' } else { '-' };
            writeln!(w, "chr{}\t{s}\t{e}\ti{idx}\t0\t{strand}", c + 1)?;
            idx += 1;
        }
    }
    Ok(w.flush()?)
}

/// Parse the generator arguments of a `cases.toml` `setup` string or the
/// command line: `bed --records N ...`, `vcf --records N ...`, `pair --records N ...`,
/// `sam --records N ...` or `fasta --records N ...`.
pub fn run_args(args: &[String], out: &mut impl Write) -> Result<()> {
    let (kind, rest) = match args.split_first() {
        Some((k, r)) => (k.as_str(), r),
        None => bail!("gen needs a kind: bed, vcf, pair, sam or fasta"),
    };
    let mut bed_o = BedOpts::default();
    let mut vcf_o = VcfOpts {
        records: 1000,
        samples: 3,
        seed: 1,
        chroms: 5,
    };
    let mut pair_o = PairOpts::default();
    let mut sam_o = SamOpts::default();
    let mut fasta_o = FastaOpts::default();
    let mut ref_out: Option<String> = None;
    let mut targets_out: Option<String> = None;
    let mut targets: Option<u64> = None;
    let mut bed_out: Option<String> = None;
    let mut out_a_path: Option<String> = None;
    let mut out_b_path: Option<String> = None;

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
            ("pair", "--records") => pair_o.records = val()?.parse()?,
            ("pair", "--depth") => pair_o.depth = val()?.parse()?,
            ("pair", "--chroms") => pair_o.chroms = val()?.parse()?,
            ("pair", "--seed") => pair_o.seed = val()?.parse()?,
            ("pair", "--len-min") => pair_o.len_min = val()?.parse()?,
            ("pair", "--len-max") => pair_o.len_max = val()?.parse()?,
            ("pair", "--sorted") => pair_o.sorted = true,
            ("pair", "--out-a") => out_a_path = Some(val()?.clone()),
            ("pair", "--out-b") => out_b_path = Some(val()?.clone()),
            ("pair", "--prefix") => {
                let p = val()?;
                out_a_path = Some(format!("{p}_a.bed"));
                out_b_path = Some(format!("{p}_b.bed"));
            }
            ("sam", "--records") => sam_o.records = val()?.parse()?,
            ("sam", "--read-len") => sam_o.read_len = val()?.parse()?,
            ("sam", "--depth") => sam_o.depth = val()?.parse()?,
            ("sam", "--chroms") => sam_o.chroms = val()?.parse()?,
            ("sam", "--spliced-frac") => sam_o.spliced_frac = val()?.parse()?,
            ("sam", "--shuffled") => sam_o.shuffled = true,
            ("sam", "--seed") => sam_o.seed = val()?.parse()?,
            ("sam", "--part") => sam_o.part = val()?.parse()?,
            ("sam", "--ref-out") => ref_out = Some(val()?.clone()),
            ("sam", "--targets-out") => targets_out = Some(val()?.clone()),
            ("sam", "--targets") => targets = Some(val()?.parse()?),
            ("fasta", "--records") => fasta_o.records = val()?.parse()?,
            ("fasta", "--chroms") => fasta_o.chroms = val()?.parse()?,
            ("fasta", "--chrom-len") => fasta_o.chrom_len = val()?.parse()?,
            ("fasta", "--line-width") => fasta_o.line_width = val()?.parse()?,
            ("fasta", "--soft-frac") => fasta_o.soft_frac = val()?.parse()?,
            ("fasta", "--n-frac") => fasta_o.n_frac = val()?.parse()?,
            ("fasta", "--seed") => fasta_o.seed = val()?.parse()?,
            ("fasta", "--len-min") => fasta_o.len_min = val()?.parse()?,
            ("fasta", "--len-max") => fasta_o.len_max = val()?.parse()?,
            ("fasta", "--bed-out") => bed_out = Some(val()?.clone()),
            _ => bail!("unknown gen argument {flag} for {kind}"),
        }
    }
    match kind {
        "bed" => bed(&bed_o, out),
        "vcf" => vcf(&vcf_o, out),
        "pair" => {
            if let (Some(pa), Some(pb)) = (out_a_path, out_b_path) {
                let mut fa = std::fs::File::create(pa)?;
                let mut fb = std::fs::File::create(pb)?;
                pair(&pair_o, &mut fa, &mut fb)
            } else {
                bail!("pair needs --out-a and --out-b (or --prefix)");
            }
        }
        "sam" => {
            if let Some(p) = ref_out {
                sam_reference(&sam_o, &mut std::fs::File::create(p)?)?;
            }
            if let Some(p) = targets_out {
                let n = targets.unwrap_or((sam_o.records / 50).max(1));
                sam_targets(&sam_o, n, &mut std::fs::File::create(p)?)?;
            }
            sam(&sam_o, out)
        }
        "fasta" => {
            if let Some(p) = bed_out {
                fasta_intervals(&fasta_o, &mut std::fs::File::create(p)?)?;
            }
            fasta(&fasta_o, out)
        }
        _ => bail!("unknown gen kind {kind} (bed, vcf, pair, sam or fasta)"),
    }
}
