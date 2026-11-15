//! clap definitions: global flags, shared option groups and the subcommands.
use clap::{Command, FromArgMatches, Parser, Subcommand};

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
    arg_required_else_help = true
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
    /// The inputs the command line gives `-` (standard input) to, in the order given, as the
    /// user can find them (`-a`, `-b (value 2)`, `A`): read from the parse by
    /// [`try_parse_from`], not an argument. A run gives them to its token before it starts
    /// ([`crate::embed::execute`]), which refuses a second.
    #[arg(skip)]
    pub stdin_inputs: Vec<String>,
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
    #[command(before_help = EXPERIMENTAL)]
    Igv(IgvArgs),
    #[command(before_help = EXPERIMENTAL)]
    Links(LinksArgs),
    #[command(before_help = EXPERIMENTAL)]
    Tag(TagArgs),
}

/// The groups of the commands outside the stable core, in the order `gloci --help` lists them:
/// the title of each group and its commands. This is the one place that says which group an
/// experimental command is in; the stable group is [`crate::CORE_COMMANDS`].
/// `every_command_is_in_one_group` fails when a command is in none or in two.
pub const EXPERIMENTAL_GROUPS: &[(&str, &[&str])] = &[
    (
        "Intervals",
        &[
            "cluster",
            "spacing",
            "complement",
            "partition",
            "genomecov",
            "slop",
            "flank",
            "shift",
            "windows",
            "map",
            "groupby",
            "convert",
            "igv",
            "links",
        ],
    ),
    (
        "Statistics and sampling",
        &[
            "stats",
            "random",
            "sample",
            "shuffle",
            "pair-intersect",
            "pair-pair",
        ],
    ),
    ("VCF", &["vcf"]),
    (
        "Alignments",
        &[
            "cat", "flags", "idxstats", "flagstat", "samples", "depth", "bedcov", "reheader", "tag",
        ],
    ),
    ("FASTA", &["getfasta", "maskfasta", "nuc", "faidx", "dict"]),
    ("Any file", &["check", "head"]),
];

/// The heading and commands of each group of the top-level help, stable core first.
pub fn command_groups() -> Vec<(String, Vec<&'static str>)> {
    let mut groups = vec![("Stable core".to_string(), crate::CORE_COMMANDS.to_vec())];
    for (title, names) in EXPERIMENTAL_GROUPS {
        groups.push((format!("{title} (experimental)"), names.to_vec()));
    }
    groups
}

/// The line under the command list.
const TIERS: &str =
    "Stable core: changes only with a version bump.\nExperimental: may change between releases.";

/// The clap command tree. clap 4.6 has one heading for all subcommands, so the root command's
/// help template lists them in [`command_groups`] with their one-line summaries; the usage,
/// the options and the version are clap's own.
pub fn command() -> Command {
    let root = <Cli as clap::CommandFactory>::command();
    let styles = root.get_styles().clone();
    let (head, reset) = (
        styles.get_header().render().to_string(),
        styles.get_header().render_reset().to_string(),
    );
    let (lit, lit_reset) = (
        styles.get_literal().render().to_string(),
        styles.get_literal().render_reset().to_string(),
    );
    let width = root
        .get_subcommands()
        .map(|c| c.get_name().len())
        .max()
        .unwrap_or(0);
    let mut list = String::new();
    for (title, names) in command_groups() {
        list.push_str(&format!("{head}{title}:{reset}\n"));
        for name in names {
            let about = root
                .find_subcommand(name)
                .and_then(|c| c.get_about())
                .map(|a| a.to_string())
                .unwrap_or_default();
            list.push_str(&format!("  {lit}{name:<width$}{lit_reset}  {about}\n"));
        }
        list.push('\n');
    }
    // `{` and `}` are template syntax; no command name or summary has one.
    let template = format!(
        "{{about-with-newline}}\n{{usage-heading}} {{usage}}\n\n{list}{head}Options:{reset}\n{{options}}\n\n{TIERS}\n"
    );
    root.help_template(template)
}

/// Parse an argv (`gloci` first) with [`command`], so that help prints with the grouped list.
pub fn try_parse_from<I, T>(args: I) -> Result<Cli, clap::Error>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let mut cmd = command();
    let matches = cmd.try_get_matches_from_mut(args)?;
    let mut cli = Cli::from_arg_matches(&matches).map_err(|e| e.format(&mut cmd))?;
    cli.stdin_inputs = crate::embed::stdin_inputs(&cmd, &matches);
    Ok(cli)
}
