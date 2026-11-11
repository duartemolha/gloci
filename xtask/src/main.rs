use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "xtask", about = "Developer tasks for gloci")]
struct Cli {
    #[command(subcommand)]
    task: Task,
}

#[derive(Subcommand)]
enum Task {
    /// Generate man pages in target/assets/man
    Man,
    /// Generate shell completions in target/assets/completions
    Completions,
    /// Generate documentation artifacts (SKILL.md)
    Docs {
        /// Write the SKILL.md agent usage guide to PATH
        #[arg(long, value_name = "PATH", default_value = "SKILL.md")]
        skill: PathBuf,
    },
    /// Generate synthetic data
    Gen {
        #[command(subcommand)]
        kind: Gen,
    },
    /// Compare gloci with the pinned reference tools (needs Docker)
    Diff {
        #[arg(long, default_value = xtask::diff::DEFAULT_IMAGE)]
        image: String,
        #[arg(long, default_value = "xtask/diff/cases.toml")]
        cases: PathBuf,
        /// Run only this case
        #[arg(long)]
        only: Option<String>,
        /// Rewrite crates/gloci/tests/golden from the reference tools
        #[arg(long)]
        update_golden: bool,
        /// gloci binary to test (default: cargo build for the host)
        #[arg(long)]
        gloci: Option<PathBuf>,
    },
    /// Maintain the built-in genomes in crates/gloci-io/data/genomes
    Genomes {
        #[command(subcommand)]
        action: GenomesAction,
    },
    /// Run the benchmark suite
    Bench {
        #[arg(long, default_value = "xtask/bench/cases.toml")]
        cases: PathBuf,
        #[arg(long, default_value = "bench-results")]
        out: PathBuf,
        /// Evaluate each case's targets and exit 1 on a miss
        #[arg(long)]
        check_targets: bool,
        /// gloci binary to benchmark (default: build --release for the host)
        #[arg(long)]
        gloci: Option<PathBuf>,
        /// Run hyperfine and the reference tools on the host, not in the container
        #[arg(long)]
        no_docker: bool,
        #[arg(long, default_value = xtask::bench::DEFAULT_IMAGE)]
        image: String,
        /// Run only this case
        #[arg(long, visible_alias = "case")]
        only: Option<String>,
        /// Override --records of every case (reduced local runs)
        #[arg(long)]
        records: Option<u64>,
        /// Where generated inputs are cached
        #[arg(long, default_value = "target/bench-data")]
        data_dir: PathBuf,
    },
}

#[derive(Subcommand)]
enum GenomesAction {
    /// Download every source in SOURCES.toml and rewrite the .genome files
    Refresh {
        /// Refresh only this genome id
        #[arg(long)]
        only: Option<String>,
        #[arg(long, default_value = xtask::genomes::DEFAULT_DIR)]
        dir: PathBuf,
    },
}

#[derive(Subcommand)]
enum Gen {
    Bed {
        #[arg(long)]
        records: u64,
        #[arg(long, default_value_t = 24)]
        chroms: u32,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[arg(long, value_enum, default_value = "uniform")]
        profile: xtask::synth::Profile,
        #[arg(long, default_value_t = 50)]
        len_min: u64,
        #[arg(long, default_value_t = 10_000)]
        len_max: u64,
        /// Write records in sorted (byte-order chromosome) order
        #[arg(long)]
        sorted: bool,
        #[arg(short, long)]
        output: PathBuf,
    },
    Vcf {
        #[arg(long)]
        records: u64,
        #[arg(long, default_value_t = 3)]
        samples: u32,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[arg(long, default_value_t = 5)]
        chroms: u32,
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Single-end alignments (SAM) whose reads match a generated reference
    Sam {
        #[arg(long)]
        records: u64,
        #[arg(long, default_value_t = 100)]
        read_len: u32,
        /// Mean depth of the aligned bases (sets the chromosome length)
        #[arg(long, default_value_t = 30.0)]
        depth: f64,
        #[arg(long, default_value_t = 24)]
        chroms: u32,
        /// Fraction of reads with an intron (CIGAR N)
        #[arg(long, default_value_t = 0.1)]
        spliced_frac: f64,
        /// Records in random order instead of coordinate order
        #[arg(long)]
        shuffled: bool,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Salt of the reads' randomness: parts of one reference differ
        #[arg(long, default_value_t = 0)]
        part: u64,
        /// Also write the reference FASTA the reads match
        #[arg(long)]
        ref_out: Option<PathBuf>,
        /// Also write a sorted BED of target windows
        #[arg(long)]
        targets_out: Option<PathBuf>,
        /// Number of target windows (default records / 50)
        #[arg(long)]
        targets: Option<u64>,
        #[arg(short, long)]
        output: PathBuf,
    },
    /// A FASTA with soft-masked and N runs, and a BED of intervals inside it
    Fasta {
        /// Intervals of the BED (--bed-out)
        #[arg(long)]
        records: u64,
        #[arg(long, default_value_t = 24)]
        chroms: u32,
        #[arg(long, default_value_t = 1_000_000)]
        chrom_len: u64,
        #[arg(long, default_value_t = 60)]
        line_width: usize,
        #[arg(long, default_value_t = 0.05)]
        soft_frac: f64,
        #[arg(long, default_value_t = 0.001)]
        n_frac: f64,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Shortest interval of the BED
        #[arg(long, default_value_t = 50)]
        len_min: u64,
        /// Longest interval of the BED
        #[arg(long, default_value_t = 10_000)]
        len_max: u64,
        #[arg(long)]
        bed_out: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
    },
    Pair {
        #[arg(long)]
        records: u64,
        #[arg(long, default_value_t = 1.0)]
        depth: f64,
        #[arg(long, default_value_t = 24)]
        chroms: u32,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[arg(long, default_value_t = 50)]
        len_min: u64,
        #[arg(long, default_value_t = 10_000)]
        len_max: u64,
        #[arg(long)]
        sorted: bool,
        #[arg(long)]
        out_a: Option<PathBuf>,
        #[arg(long)]
        out_b: Option<PathBuf>,
        #[arg(long)]
        prefix: Option<PathBuf>,
    },
}

fn run(cli: Cli) -> Result<bool> {
    match cli.task {
        Task::Man => xtask::docs::generate_man(Path::new("target/assets/man"))?,
        Task::Completions => {
            xtask::docs::generate_completions(Path::new("target/assets/completions"))?
        }
        Task::Docs { skill } => xtask::docs::generate_skill(&skill)?,
        Task::Gen { kind } => match kind {
            Gen::Bed {
                records,
                chroms,
                seed,
                profile,
                len_min,
                len_max,
                sorted,
                output,
            } => {
                let o = xtask::synth::BedOpts {
                    records,
                    chroms,
                    seed,
                    profile,
                    len_min,
                    len_max,
                    sorted,
                };
                xtask::synth::bed(&o, &mut std::fs::File::create(output)?)?
            }
            Gen::Vcf {
                records,
                samples,
                seed,
                chroms,
                output,
            } => {
                let o = xtask::synth::VcfOpts {
                    records,
                    samples,
                    seed,
                    chroms,
                };
                xtask::synth::vcf(&o, &mut std::fs::File::create(output)?)?
            }
            Gen::Sam {
                records,
                read_len,
                depth,
                chroms,
                spliced_frac,
                shuffled,
                seed,
                part,
                ref_out,
                targets_out,
                targets,
                output,
            } => {
                let o = xtask::synth::SamOpts {
                    records,
                    read_len,
                    depth,
                    chroms,
                    spliced_frac,
                    shuffled,
                    seed,
                    part,
                };
                if let Some(p) = ref_out {
                    xtask::synth::sam_reference(&o, &mut std::fs::File::create(p)?)?;
                }
                if let Some(p) = targets_out {
                    let n = targets.unwrap_or((records / 50).max(1));
                    xtask::synth::sam_targets(&o, n, &mut std::fs::File::create(p)?)?;
                }
                xtask::synth::sam(&o, &mut std::fs::File::create(output)?)?
            }
            Gen::Fasta {
                records,
                chroms,
                chrom_len,
                line_width,
                soft_frac,
                n_frac,
                seed,
                len_min,
                len_max,
                bed_out,
                output,
            } => {
                let o = xtask::synth::FastaOpts {
                    records,
                    chroms,
                    chrom_len,
                    line_width,
                    soft_frac,
                    n_frac,
                    seed,
                    len_min,
                    len_max,
                };
                if let Some(p) = bed_out {
                    xtask::synth::fasta_intervals(&o, &mut std::fs::File::create(p)?)?;
                }
                xtask::synth::fasta(&o, &mut std::fs::File::create(output)?)?
            }
            Gen::Pair {
                records,
                depth,
                chroms,
                seed,
                len_min,
                len_max,
                sorted,
                out_a,
                out_b,
                prefix,
            } => {
                let (pa, pb) = match (out_a, out_b, prefix) {
                    (Some(a), Some(b), _) => (a, b),
                    (_, _, Some(p)) => {
                        let stem = p.to_string_lossy();
                        (
                            PathBuf::from(format!("{stem}_a.bed")),
                            PathBuf::from(format!("{stem}_b.bed")),
                        )
                    }
                    _ => anyhow::bail!("pair needs --out-a and --out-b (or --prefix)"),
                };
                let o = xtask::synth::PairOpts {
                    records,
                    depth,
                    chroms,
                    seed,
                    len_min,
                    len_max,
                    sorted,
                };
                xtask::synth::pair(
                    &o,
                    &mut std::fs::File::create(pa)?,
                    &mut std::fs::File::create(pb)?,
                )?
            }
        },
        Task::Genomes { action } => match action {
            GenomesAction::Refresh { only, dir } => {
                return xtask::genomes::run(&xtask::genomes::Args { dir, only });
            }
        },
        Task::Diff {
            image,
            cases,
            only,
            update_golden,
            gloci,
        } => {
            return xtask::diff::run(&xtask::diff::Args {
                image,
                cases,
                only,
                update_golden,
                gloci,
            });
        }
        Task::Bench {
            cases,
            out,
            check_targets,
            gloci,
            no_docker,
            image,
            only,
            records,
            data_dir,
        } => {
            return xtask::bench::run(&xtask::bench::Args {
                cases,
                out,
                check_targets,
                gloci,
                no_docker,
                image,
                only,
                records,
                data_dir,
            });
        }
    }
    Ok(true)
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
