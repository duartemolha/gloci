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
