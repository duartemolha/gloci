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
    /// Maintain the built-in genomes in data/genomes
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
        #[arg(long)]
        case: Option<String>,
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
}

fn run(cli: Cli) -> Result<bool> {
    match cli.task {
        Task::Man => xtask::docs::generate_man(Path::new("target/assets/man"))?,
        Task::Completions => {
            xtask::docs::generate_completions(Path::new("target/assets/completions"))?
        }
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
            case,
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
                only: case,
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
