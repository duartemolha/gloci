use clap::builder::TypedValueParser as _;
use clap::{Args, ValueEnum};
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// Parse a size such as `512M` or `2G`: an integer with a required K/M/G/T
/// suffix in binary units (powers of 1024). `10` and `1.5G` are errors.
pub fn parse_size(s: &str) -> Result<u64, String> {
    let bad = || {
        format!("invalid size {s:?}: use an integer with a K, M, G or T suffix, e.g. 512M or 2G")
    };
    let last = s.chars().last().ok_or_else(bad)?;
    let shift = match last.to_ascii_uppercase() {
        'K' => 10,
        'M' => 20,
        'G' => 30,
        'T' => 40,
        _ => return Err(bad()),
    };
    let digits = &s[..s.len() - 1];
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad());
    }
    let n: u64 = digits.parse().map_err(|_| bad())?;
    n.checked_mul(1 << shift)
        .ok_or_else(|| format!("size {s:?} is too large"))
}

fn parse_threads(s: &str) -> Result<NonZeroUsize, String> {
    s.parse::<NonZeroUsize>()
        .map_err(|_| format!("invalid thread count {s:?}: use an integer of at least 1"))
}

pub(super) fn parse_level(s: &str) -> Result<u8, String> {
    match s.parse::<u8>() {
        Ok(n) if n <= 12 => Ok(n),
        _ => Err(format!("invalid compression level {s:?}: use 0-12")),
    }
}

/// Output destination and encoding.
#[derive(Args, Debug, Clone)]
pub struct OutputArgs {
    /// Write output to FILE (format from its extension) instead of stdout. A
    /// device or pipe such as /dev/null or >(cmd) is written to directly.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Output format when writing to stdout.
    #[arg(short = 'O', long, value_name = "FMT")]
    pub output_format: Option<String>,
    /// Also write a .tbi/.csi index next to the output (.bai or .csi for BAM output).
    #[arg(long)]
    pub write_index: bool,
    /// BGZF compression level, 0-12.
    #[arg(long, value_name = "N", default_value = "6", value_parser = parse_level)]
    pub compression_level: u8,
}

/// Thread count.
#[derive(Args, Debug, Clone)]
pub struct ThreadArgs {
    /// Number of threads (default: all available cores).
    #[arg(short = 't', long, value_name = "N", env = "GLOCI_THREADS", value_parser = parse_threads)]
    pub threads: Option<NonZeroUsize>,
}

#[derive(Args, Debug, Clone)]
pub struct GenomeArg {
    /// Chromosome order and sizes: a built-in genome such as hg38 (see
    /// gloci genomes list), or a .fai, .genome, .dict, VCF/BCF or BAM header file.
    #[arg(short, long, value_name = "NAME|FILE")]
    pub genome: Option<String>,
}

/// The help of `--chrom-order` on a command with `-g`.
pub const CHROM_ORDER_HELP: &str = "Chromosome order when neither -g nor the header of a \
VCF/BCF or alignment input sets one: lex (byte order, as sort -k1,1, sort-bed and bedtools sort; \
the default) or karyotype (chr1..chr22, chrX, chrY, chrM)";

/// The help of `--chrom-order` on a command without `-g` (see [`help_without_g`]).
pub const CHROM_ORDER_HELP_NO_G: &str = "Chromosome order when no input header lists contigs: lex \
(byte order, as sort -k1,1, sort-bed and bedtools sort; the default) or karyotype (chr1..chr22, \
chrX, chrY, chrM)";

#[derive(Args, Debug, Clone)]
pub struct ChromOrderArg {
    #[arg(
        long,
        help = CHROM_ORDER_HELP,
        value_name = "ORDER",
        default_value = "lex",
        value_parser = clap::builder::PossibleValuesParser::new(["karyotype", "lex"])
            .try_map(|s| s.parse::<gloci_core::ChromOrder>())
    )]
    pub chrom_order: gloci_core::ChromOrder,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum VcfSpan {
    Allele,
    Ref,
}

#[derive(Args, Debug, Clone)]
pub struct VcfSpanArg {
    /// Which footprint a VCF record covers: the allele span or the REF span.
    #[arg(long, value_enum, default_value = "allele", value_name = "MODE")]
    pub vcf_span: VcfSpan,
}

/// The help of `--chrom-map` on a command with `-g`.
pub const CHROM_MAP_HELP: &str = "Rename chromosomes on input with a two-column table (old name, \
new name; tab-separated, `#` comments allowed). Unlisted names are kept. A file that has both an \
old name and its new name is refused. Order checks, -g and the output all see the new names; -g \
itself is never renamed";

/// The help of `--chrom-map` on a command without `-g` (see [`help_without_g`]).
pub const CHROM_MAP_HELP_NO_G: &str = "Rename chromosomes on input with a two-column table (old \
name, new name; tab-separated, `#` comments allowed). Unlisted names are kept. A file that has \
both an old name and its new name is refused. Order checks and the output see the new names";

/// `--chrom-map` (D62), on every command that reads intervals or variants.
#[derive(Args, Debug, Clone)]
pub struct ChromMapArg {
    #[arg(long, value_name = "FILE", help = CHROM_MAP_HELP)]
    pub chrom_map: Option<String>,
}

/// `cmd` and its subcommands, with the shared help of `--chrom-order` and `--chrom-map`,
/// which names `-g`, said without it on a command that has no `-g`. A command whose help for
/// them is its own keeps it.
/// (The genome option is found by its long name: `vcf merge` has a hidden `-g/--gvcf`.)
pub fn help_without_g(cmd: clap::Command) -> clap::Command {
    let has_g = cmd.get_arguments().any(|a| a.get_long() == Some("genome"));
    let cmd = if has_g {
        cmd
    } else {
        cmd.mut_args(|a| {
            let help = a.get_help().map(|h| h.to_string());
            match (a.get_id().as_str(), help.as_deref()) {
                ("chrom_order", Some(CHROM_ORDER_HELP)) => a.help(CHROM_ORDER_HELP_NO_G),
                ("chrom_map", Some(CHROM_MAP_HELP)) => a.help(CHROM_MAP_HELP_NO_G),
                _ => a,
            }
        })
    };
    cmd.mut_subcommands(help_without_g)
}

impl ChromMapArg {
    /// Read the table, when one was given.
    pub fn load(&self) -> Result<Option<std::sync::Arc<gloci_io::ChromMap>>, gloci_io::IoError> {
        self.chrom_map
            .as_deref()
            .map(|p| gloci_io::ChromMap::load(p).map(std::sync::Arc::new))
            .transpose()
    }
}

/// The alignment filters (D101), flattened into each command that reads alignments.
#[derive(Args, Debug, Clone, Default)]
pub struct AlignArgs {
    /// Skip alignments with a mapping quality below N (BAM, CRAM and SAM input).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_mapq: u8,
    /// Keep only alignments that have all of these flags: a number (0x2 or 2)
    /// or names such as paired,proper_pair (as samtools view -f).
    #[arg(long, value_name = "FLAGS")]
    pub require_flags: Option<String>,
    /// Drop alignments that have any of these flags: a number or names such as
    /// unmap,secondary,dup (as samtools view -F).
    #[arg(long, value_name = "FLAGS")]
    pub exclude_flags: Option<String>,
    /// The FASTA a CRAM input was made with (bedtools --cram-ref), with its
    /// .fai (and .gzi when BGZF). Not needed for a CRAM written without a
    /// reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
}

impl AlignArgs {
    /// The reader options; `full_records` keeps whole records (for `view` and `sort`).
    pub fn options(
        &self,
        full_records: bool,
    ) -> Result<gloci_io::AlignOptions, gloci_io::UsageError> {
        let flags = |name: &str, v: &Option<String>| match v {
            None => Ok(0),
            Some(s) => {
                gloci_io::parse_flags(s).map_err(|e| gloci_io::UsageError(format!("{name}: {e}")))
            }
        };
        Ok(gloci_io::AlignOptions {
            min_mapq: self.min_mapq,
            require_flags: flags("--require-flags", &self.require_flags)?,
            exclude_flags: flags("--exclude-flags", &self.exclude_flags)?,
            full_records,
            reference: self.reference.clone(),
            score_tag: None,
        })
    }
}

#[derive(Args, Debug, Clone)]
pub struct InvalidArg {
    /// Skip invalid records with a warning instead of failing.
    #[arg(long)]
    pub skip_invalid: bool,
}

/// `--sort` and temporary directory options for relation commands.
#[derive(Args, Debug, Clone)]
pub struct RelationSortArg {
    /// Sort the inputs first, for inputs that are unsorted or too large to
    /// sort in memory (external sort bounded by --max-mem). The sorted
    /// records are written to a temporary file in --tmp-dir and read back,
    /// so besides the spill runs (about the size of the input, compressed)
    /// the directory needs one uncompressed sorted copy: allow roughly twice
    /// the input size. Everything is removed afterwards.
    #[arg(long)]
    pub sort: bool,
    /// Memory budget for --sort, in binary units (K, M, G, T); at least 16M.
    #[arg(long, value_name = "SIZE", default_value = "2G", value_parser = parse_size, requires = "sort")]
    pub max_mem: u64,
    /// Directory for temporary files: --sort's, and the spills of
    /// per-chromosome runs on indexed and plain files; ignored when nothing spills. It
    /// must exist and be writable (default: $GLOCI_TMPDIR if set and
    /// non-empty, else the system temp dir).
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct TmpDirArg {
    /// Directory for temporary spill files; it must exist and be writable
    /// (default: $GLOCI_TMPDIR if set and non-empty, else the system temp dir).
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct RegionArgs {
    /// Restrict to a region such as chr1:1,000-2,000 (repeatable).
    #[arg(long, value_name = "REGION")]
    pub region: Vec<String>,
    /// Read regions from FILE (BED: 0-based, half-open); see --region.
    #[arg(long, value_name = "FILE")]
    pub regions_file: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct HeaderArg {
    /// Keep header lines in the output.
    #[arg(long)]
    pub header: bool,
}

#[derive(Args, Debug, Clone)]
pub struct FormatArg {
    /// Input format, overriding detection: bed, bedgraph, wig, gff, gtf, vcf
    /// or bcf; bigwig or bigbed (from a file); bam, cram or sam where the
    /// command reads alignments; psl or rmsk for convert.
    #[arg(long, value_name = "FMT")]
    pub format: Option<String>,
}

/// Parse a fraction such as `0.5` or `50%` (spec 4.3).
pub(super) fn parse_fraction(s: &str) -> Result<gloci_core::Fraction, String> {
    gloci_core::Fraction::parse(s).map_err(|e| e.to_string())
}

/// The A input of a relation command.
#[derive(Args, Debug, Clone)]
pub struct AArg {
    /// A: the file whose records are reported (`-` for stdin).
    #[arg(short = 'a', long = "a", value_name = "FILE", required = true)]
    pub a: String,
    /// A's format, overriding detection: bed, bedgraph, wig, gff, gtf, vcf or
    /// bcf; bigwig or bigbed (from a file); bam, cram or sam where the
    /// command reads alignments.
    #[arg(long, value_name = "FMT")]
    pub a_format: Option<String>,
}

/// Overlap thresholds, strand rule and block mode (spec 4.3).
#[derive(Args, Debug, Clone)]
pub struct OverlapArgs {
    /// Minimum fraction of A covered, e.g. 0.5 or 50% (4.3: fractions).
    #[arg(short = 'f', long, value_name = "FRAC", value_parser = parse_fraction)]
    pub min_frac_a: Option<gloci_core::Fraction>,
    /// Minimum fraction of B covered, e.g. 0.5 or 50% (4.3: fractions).
    #[arg(short = 'F', long, value_name = "FRAC", value_parser = parse_fraction)]
    pub min_frac_b: Option<gloci_core::Fraction>,
    /// Require both fractions to be met (4.3: -r).
    #[arg(short = 'r', long)]
    pub reciprocal: bool,
    /// Require either fraction to be met (4.3: -e).
    #[arg(short = 'e', long)]
    pub either: bool,
    /// Minimum number of overlapping bases (4.3: --min-bases).
    #[arg(long, value_name = "N")]
    pub min_bases: Option<u64>,
    /// Only overlaps on the same strand; unknown strands never match (4.3: strand).
    #[arg(short = 's', long = "same-strand", conflicts_with = "opposite_strand")]
    pub same_strand: bool,
    /// Only overlaps on opposite strands; unknown strands never match (4.3: strand).
    #[arg(short = 'S', long = "opposite-strand")]
    pub opposite_strand: bool,
    /// Compare BED12 blocks (exons) and the aligned segments of BAM/SAM reads (cut at
    /// CIGAR N) instead of whole spans (4.3: BED12 blocks).
    #[arg(long)]
    pub split: bool,
    /// A point on the edge of an interval overlaps it (4.3: points).
    #[arg(long)]
    pub inclusive_points: bool,
}

impl OverlapArgs {
    pub fn criteria(&self) -> gloci_core::OverlapCriteria {
        use gloci_core::StrandRule;
        gloci_core::OverlapCriteria {
            min_frac_a: self.min_frac_a,
            min_frac_b: self.min_frac_b,
            reciprocal: self.reciprocal,
            either: self.either,
            min_bases: self.min_bases,
            strand: match (self.same_strand, self.opposite_strand) {
                (true, _) => StrandRule::Same,
                (_, true) => StrandRule::Opposite,
                _ => StrandRule::Any,
            },
            inclusive_points: self.inclusive_points,
        }
    }
}

/// B padding (spec 4.3: windows). Cannot be combined with thresholds (D37).
#[derive(Args, Debug, Clone)]
#[command(group = clap::ArgGroup::new("pad").multiple(true).args(["window", "left", "right"]))]
pub struct WindowArgs {
    /// Pad B by N on both sides before testing (4.3: windows).
    #[arg(
        short = 'w',
        long = "window",
        value_name = "N",
        conflicts_with_all = WINDOW_CONFLICTS
    )]
    pub window: Option<u64>,
    /// Pad B by N on the left (4.3: windows).
    #[arg(long, value_name = "N", conflicts_with_all = WINDOW_CONFLICTS)]
    pub left: Option<u64>,
    /// Pad B by N on the right (4.3: windows).
    #[arg(long, value_name = "N", conflicts_with_all = WINDOW_CONFLICTS)]
    pub right: Option<u64>,
    /// Left and right follow A's strand (4.3: windows).
    #[arg(long, requires = "pad", conflicts_with_all = WINDOW_CONFLICTS)]
    pub strand_window: bool,
}

/// The threshold flags a window cannot be combined with (D37).
const WINDOW_CONFLICTS: [&str; 5] = [
    "min_frac_a",
    "min_frac_b",
    "reciprocal",
    "either",
    "min_bases",
];

impl WindowArgs {
    /// The window, when any window flag is given. `--left`/`--right`
    /// override `-w` for their side.
    pub fn window(&self) -> Option<gloci_core::Window> {
        if self.window.is_none() && self.left.is_none() && self.right.is_none() {
            return None;
        }
        let w = self.window.unwrap_or(0);
        Some(gloci_core::Window {
            left: self.left.unwrap_or(w),
            right: self.right.unwrap_or(w),
            strand_aware: self.strand_window,
        })
    }
}

/// Sort-free mode for B (spec 8.2).
#[derive(Args, Debug, Clone)]
pub struct UnsortedArg {
    /// Accept inputs in any order and keep A's order in the output; every -b
    /// file is held in memory. Cannot be combined with --sort.
    #[arg(long, conflicts_with = "sort")]
    pub unsorted: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{BArgs, SortArg};
    use clap::Parser;
    use gloci_core::ChromOrder;

    #[derive(Parser)]
    struct T {
        #[command(flatten)]
        o: ChromOrderArg,
    }

    #[test]
    fn align_args_options() {
        let a = AlignArgs {
            min_mapq: 20,
            require_flags: Some("paired".into()),
            exclude_flags: Some("0x904".into()),
            reference: None,
        };
        let o = a.options(true).unwrap();
        assert_eq!(
            (o.min_mapq, o.require_flags, o.exclude_flags),
            (20, 1, 0x904)
        );
        assert!(o.full_records);
        let bad = AlignArgs {
            exclude_flags: Some("bogus".into()),
            ..AlignArgs::default()
        };
        assert!(
            bad.options(false)
                .unwrap_err()
                .0
                .contains("--exclude-flags")
        );
    }

    #[test]
    fn chrom_order_arg() {
        let order = |a: &[&str]| T::try_parse_from(a).map(|t| t.o.chrom_order);
        assert_eq!(order(&["t"]).unwrap(), ChromOrder::Lex);
        assert_eq!(
            order(&["t", "--chrom-order", "karyotype"]).unwrap(),
            ChromOrder::Karyotype
        );
        assert!(order(&["t", "--chrom-order", "natural"]).is_err());
    }

    #[derive(Parser)]
    struct R {
        #[command(flatten)]
        a: AArg,
        #[command(flatten)]
        b: BArgs,
        #[command(flatten)]
        o: OverlapArgs,
        #[command(flatten)]
        w: WindowArgs,
        #[command(flatten)]
        u: UnsortedArg,
        #[command(flatten)]
        sort: SortArg,
    }

    fn parse(args: &[&str]) -> Result<R, clap::Error> {
        let mut v = vec!["t", "-a", "a.bed", "-b", "b.bed"];
        v.extend_from_slice(args);
        R::try_parse_from(v)
    }

    #[test]
    fn relation_groups_parse() {
        let r = R::try_parse_from([
            "t",
            "-a",
            "a.bed",
            "--a-format",
            "bed",
            "-b",
            "x.bed",
            "y.bed",
            "-b",
            "z.bed",
            "--names",
            "x,y,z",
            "--b-format",
            "bed",
            "-f",
            "50%",
            "-F",
            "0.25",
            "-r",
            "-e",
            "--min-bases",
            "3",
            "-s",
            "--split",
            "--inclusive-points",
            "--unsorted",
        ])
        .unwrap();
        assert_eq!(r.a.a, "a.bed");
        assert_eq!(r.b.b, ["x.bed", "y.bed", "z.bed"]);
        assert_eq!(r.b.labels().unwrap(), ["x", "y", "z"]);
        let c = r.o.criteria();
        assert_eq!(
            c.min_frac_a,
            Some(gloci_core::Fraction::parse("0.5").unwrap())
        );
        assert_eq!(c.min_bases, Some(3));
        assert_eq!(c.strand, gloci_core::StrandRule::Same);
        assert!(r.o.split && r.u.unsorted && r.w.window().is_none());
        assert!(parse(&["-f", "1.5"]).is_err());
        assert!(
            R::try_parse_from(["t", "-a", "a.bed"]).is_err(),
            "-b is required"
        );
    }

    #[test]
    fn names_must_match_b_count() {
        assert!(parse(&["--names", "x,y"]).unwrap().b.labels().is_err());
        let r = parse(&[]).unwrap();
        assert_eq!(r.b.labels().unwrap(), ["b.bed"]);
    }

    #[test]
    fn window_conflicts_with_thresholds() {
        use clap::error::ErrorKind::ArgumentConflict;
        for w in [["-w", "5"], ["--left", "5"], ["--right", "5"]] {
            for t in [
                &["-f", "0.5"][..],
                &["-F", "0.5"],
                &["-r"],
                &["-e"],
                &["--min-bases", "1"],
            ] {
                let args: Vec<&str> = w.iter().chain(t).copied().collect();
                let e = parse(&args)
                    .err()
                    .unwrap_or_else(|| panic!("{args:?} parsed"));
                assert_eq!(e.kind(), ArgumentConflict, "{args:?}");
                assert_eq!(e.exit_code(), 2);
            }
        }
        assert!(parse(&["--strand-window", "-r"]).is_err());
        assert!(parse(&["--strand-window"]).is_err(), "needs a padding flag");
        let r = parse(&["-w", "10", "--left", "2", "--strand-window", "-s"]).unwrap();
        let w = r.w.window().unwrap();
        assert_eq!((w.left, w.right, w.strand_aware), (2, 10, true));
    }

    #[test]
    fn same_and_opposite_strand_conflict() {
        assert!(parse(&["-s", "-S"]).is_err());
        assert!(parse(&["-S"]).is_ok());
    }
}
