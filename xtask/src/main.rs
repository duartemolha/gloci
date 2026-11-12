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
    /// Paper datasets (D95): never run by CI
    Paper {
        #[command(subcommand)]
        action: PaperAction,
    },
    /// Maintain the built-in genomes in crates/gloci-io/data/genomes
    Genomes {
        #[command(subcommand)]
        action: GenomesAction,
    },
    /// Check the release metadata (version, CHANGELOG, CITATION.cff, LICENSE); exit 1 on a problem
    ReleaseCheck {
        /// The pushed tag, which must be v<workspace version>
        #[arg(long)]
        tag: Option<String>,
    },
    /// Print a version's CHANGELOG section (the release notes)
    ReleaseNotes {
        /// The version, without the leading v
        #[arg(long)]
        version: String,
        /// Write to this file instead of stdout
        #[arg(long)]
        out: Option<PathBuf>,
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
enum PaperAction {
    /// Download the paper datasets and build the derived inputs
    Fetch {
        /// Write the sha256 of each derived input that has none
        #[arg(long)]
        record: bool,
        /// Only this dataset or derived input
        #[arg(long)]
        only: Option<String>,
        /// Build derived inputs on the host, not in the bench image
        #[arg(long)]
        no_docker: bool,
        #[arg(long, default_value = "gloci-bench:p1")]
        image: String,
        #[arg(long, default_value = "target/bench-data")]
        data_dir: PathBuf,
    },
    /// Regenerate paper/fetch-datasets.sh from datasets.toml
    FetchScript {
        /// Write here instead of paper/fetch-datasets.sh
        #[arg(long)]
        out: Option<PathBuf>,
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
        /// With --out-b: B's called genotypes are redrawn with this probability
        #[arg(long, default_value_t = 0.02)]
        discord: f64,
        /// Also write a second VCF of the same sites whose genotypes differ (vcf compare)
        #[arg(long)]
        out_b: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
    },
    /// BEDPE pairs, and optionally a second file to test them against
    Bedpe {
        #[arg(long)]
        pairs: u64,
        #[arg(long, default_value_t = 24)]
        chroms: u32,
        /// Fraction of the pairs with the ends on two chromosomes
        #[arg(long, default_value_t = 0.1)]
        inter_frac: f64,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[arg(long, default_value_t = 50)]
        len_min: u64,
        #[arg(long, default_value_t = 1000)]
        len_max: u64,
        /// Write to --out-b the first N pairs of A's stream, each end moved by up to 40 bases
        #[arg(long, conflicts_with = "b_bed")]
        b_pairs: Option<u64>,
        /// Write this many random BED intervals to --out-b
        #[arg(long)]
        b_bed: Option<u64>,
        #[arg(long)]
        out_b: Option<PathBuf>,
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
        /// Records are read pairs (--records counts both reads); needs --by-name
        #[arg(long)]
        paired: bool,
        /// Name-grouped order, the mates of a pair next to each other (SO:queryname)
        #[arg(long)]
        by_name: bool,
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
                discord,
                out_b,
                output,
            } => {
                let o = xtask::synth::VcfOpts {
                    records,
                    samples,
                    seed,
                    chroms,
                };
                let mut a = std::fs::File::create(output)?;
                match out_b {
                    Some(b) => {
                        xtask::synth::vcf_pair(&o, discord, &mut a, &mut std::fs::File::create(b)?)?
                    }
                    None => xtask::synth::vcf(&o, &mut a)?,
                }
            }
            Gen::Bedpe {
                pairs,
                chroms,
                inter_frac,
                seed,
                len_min,
                len_max,
                b_pairs,
                b_bed,
                out_b,
                output,
            } => {
                let o = xtask::synth::BedpeOpts {
                    pairs,
                    chroms,
                    inter_frac,
                    seed,
                    len_min,
                    len_max,
                };
                xtask::synth::bedpe(&o, &mut std::fs::File::create(output)?)?;
                match (out_b, b_pairs, b_bed) {
                    (Some(p), Some(n), _) => {
                        let ob = xtask::synth::BedpeOpts { pairs: n, ..o };
                        xtask::synth::bedpe_jittered(
                            &ob,
                            xtask::synth::B_JITTER,
                            &mut std::fs::File::create(p)?,
                        )?
                    }
                    (Some(p), None, Some(n)) => {
                        xtask::synth::bedpe_bed(&o, n, &mut std::fs::File::create(p)?)?
                    }
                    (None, None, None) => {}
                    _ => anyhow::bail!("--out-b goes with exactly one of --b-pairs and --b-bed"),
                }
            }
            Gen::Sam {
                records,
                read_len,
                depth,
                chroms,
                spliced_frac,
                shuffled,
                paired,
                by_name,
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
                    paired,
                    by_name,
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
        Task::Paper { action } => match action {
            PaperAction::Fetch {
                record,
                only,
                no_docker,
                image,
                data_dir,
            } => {
                return xtask::paper::fetch::run(&xtask::paper::fetch::FetchArgs {
                    record,
                    only,
                    image: (!no_docker).then_some(image),
                    data_dir,
                });
            }
            PaperAction::FetchScript { out } => {
                xtask::paper::fetch::write_script(out.as_deref())?;
            }
        },
        Task::Genomes { action } => match action {
            GenomesAction::Refresh { only, dir } => {
                return xtask::genomes::run(&xtask::genomes::Args { dir, only });
            }
        },
        Task::ReleaseCheck { tag } => {
            let root = xtask::diff::workspace_root();
            let problems = xtask::release::release_check(&root, tag.as_deref())?;
            for p in &problems {
                eprintln!("release-check: {p}");
            }
            return Ok(problems.is_empty());
        }
        Task::ReleaseNotes { version, out } => {
            let root = xtask::diff::workspace_root();
            let text = std::fs::read_to_string(root.join("CHANGELOG.md"))?;
            let notes = xtask::release::release_notes(&text, &version)?;
            match out {
                Some(p) => std::fs::write(p, format!("{notes}\n"))?,
                None => println!("{notes}"),
            }
        }
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
