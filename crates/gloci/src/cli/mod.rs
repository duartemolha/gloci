//! clap definitions: global flags, shared option groups and the subcommands.
use clap::{Parser, Subcommand};

mod alignments;
mod fasta;
mod files;
mod intervals;
mod relate;
mod sampling;
mod shared;
mod vcf;
mod vcf_edit;

pub use alignments::*;
pub use fasta::*;
pub use files::*;
pub use intervals::*;
pub use relate::*;
pub use sampling::*;
pub use shared::*;
pub use vcf::*;
pub use vcf_edit::*;

/// The first line of the help of every command outside the stable core
/// ([`crate::CORE_COMMANDS`]); `every_command_is_core_or_marked_experimental` keeps the two in step.
pub const EXPERIMENTAL: &str = "Experimental: checked against the tool it replaces, but its options and output may \
still change between releases.";

#[derive(Parser)]
#[command(
    name = "gloci",
    version,
    about = "Genomic interval toolkit",
    arg_required_else_help = true,
    after_help = "Stable core: genomes, index, sort, view, intersect, merge, subtract, closest, coverage.\n\
Every other command is experimental: checked against the tool it replaces, but its options and\n\
output may still change between releases."
)]
pub struct Cli {
    /// Only print errors.
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    pub quiet: bool,
    /// Print more progress information (repeat for debug output).
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    Genomes(GenomesArgs),
    Index(IndexArgs),
    Sort(SortArgs),
    View(ViewArgs),
    #[command(before_help = EXPERIMENTAL)]
    Cat(CatArgs),
    #[command(before_help = EXPERIMENTAL)]
    Check(CheckArgs),
    Intersect(IntersectArgs),
    Merge(MergeArgs),
    #[command(before_help = EXPERIMENTAL)]
    Cluster(ClusterArgs),
    #[command(before_help = EXPERIMENTAL)]
    Spacing(SpacingArgs),
    #[command(before_help = EXPERIMENTAL)]
    Complement(ComplementArgs),
    #[command(before_help = EXPERIMENTAL)]
    Partition(PartitionArgs),
    #[command(before_help = EXPERIMENTAL)]
    Genomecov(GenomecovArgs),
    #[command(before_help = EXPERIMENTAL)]
    Slop(SlopArgs),
    #[command(before_help = EXPERIMENTAL)]
    Flank(FlankArgs),
    #[command(before_help = EXPERIMENTAL)]
    Shift(ShiftArgs),
    #[command(before_help = EXPERIMENTAL)]
    Windows(WindowsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Convert(ConvertArgs),
    #[command(before_help = EXPERIMENTAL)]
    Getfasta(GetfastaArgs),
    #[command(before_help = EXPERIMENTAL)]
    Nuc(NucArgs),
    Subtract(SubtractArgs),
    Closest(ClosestArgs),
    #[command(before_help = EXPERIMENTAL)]
    Map(MapArgs),
    Coverage(CoverageArgs),
    #[command(before_help = EXPERIMENTAL)]
    Groupby(GroupbyArgs),
    #[command(before_help = EXPERIMENTAL)]
    Maskfasta(MaskfastaArgs),
    #[command(before_help = EXPERIMENTAL)]
    Vcf(Box<VcfArgs>),
    // P3b Task 2
    #[command(before_help = EXPERIMENTAL)]
    Stats(StatsArgs),
    // P3b Task 7
    #[command(before_help = EXPERIMENTAL)]
    PairIntersect(PairIntersectArgs),
    // P3b Task 8
    #[command(before_help = EXPERIMENTAL)]
    PairPair(PairPairArgs),
    // P3b Task 5
    #[command(before_help = EXPERIMENTAL)]
    Random(RandomArgs),
    #[command(before_help = EXPERIMENTAL)]
    Sample(SampleArgs),
    // P3b Task 6
    #[command(before_help = EXPERIMENTAL)]
    Shuffle(ShuffleArgs),
    #[command(before_help = EXPERIMENTAL)]
    Faidx(FaidxArgs),
    #[command(before_help = EXPERIMENTAL)]
    Dict(DictArgs),
    #[command(before_help = EXPERIMENTAL)]
    Flags(FlagsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Idxstats(IdxstatsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Flagstat(FlagstatArgs),
    #[command(before_help = EXPERIMENTAL)]
    Head(HeadArgs),
    #[command(before_help = EXPERIMENTAL)]
    Samples(SamplesArgs),
    #[command(before_help = EXPERIMENTAL)]
    Depth(DepthArgs),
    #[command(before_help = EXPERIMENTAL)]
    Bedcov(BedcovArgs),
    #[command(before_help = EXPERIMENTAL)]
    Reheader(ReheaderArgs),
}
