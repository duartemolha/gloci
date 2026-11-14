use super::shared::{ChromMapArg, ChromOrderArg, InvalidArg, OutputArgs, ThreadArgs};
use clap::{Args, ValueEnum};
use std::path::PathBuf;

/// Write VCF/BCF files one after another (bcftools concat, vcf-concat).
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf concat chr1.vcf.gz chr2.vcf.gz -o all.vcf.gz
  gloci vcf concat chr*.bcf -o all.bcf
  gloci vcf concat --naive chr1.vcf.gz chr2.vcf.gz -o all.vcf.gz")]
pub struct VcfConcatArgs {
    /// Input VCF or BCF files, in output order (`-` for stdin). The samples
    /// must be identical and in the same order, the INFO, FILTER, FORMAT,
    /// ALT and contig definitions must match, and the last record of each
    /// file must sort at or before the first record of the next (exit 1
    /// otherwise; gloci sort fixes the order). The output has the first
    /// file's header. No ##glociCommand line is added, because keeping the
    /// first file's header unchanged lets --naive and checked output stay
    /// identical.
    #[arg(value_name = "FILES", required = true, num_args = 1..)]
    pub files: Vec<String>,
    /// Copy the BGZF blocks of the inputs without recompressing them. Needs
    /// BGZF-compressed VCF files with byte-identical headers (exit 2
    /// otherwise; BCF inputs are not supported) and writes BGZF VCF. The
    /// order is checked across files only, from each file's first and last
    /// record. The output is the same as without --naive, after decompression
    /// (bcftools concat -n).
    #[arg(long)]
    pub naive: bool,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    /// Applies to every input: contig names in each header and the records
    /// are renamed before the headers are compared, the order is checked and
    /// the output is written. Not allowed with --naive, which copies blocks
    /// unchanged.
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

// P3b Task 20

/// Which multi-allelic records `vcf split` splits.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum SplitTypesArg {
    Snps,
    Indels,
    Both,
    Any,
}

/// Split multi-allelic records into biallelic ones (bcftools norm -m-).
///
/// A record with N ALT alleles becomes N records, one per ALT, in ALT
/// order; ID, QUAL and FILTER are copied. In each, the kept ALT is
/// allele 1 in GT and every other allele is 0 (missing stays missing;
/// phase and ploidy are kept). INFO and FORMAT fields follow their
/// declared Number: A keeps that ALT's value, R keeps REF's and that
/// ALT's, G keeps that ALT's genotype values (haploid and diploid);
/// other fields are copied. INFO and FORMAT keys the header does not
/// declare are dropped from split records. Alleles are not trimmed or
/// left-aligned. A record with one ALT is written as it is. A field
/// with the wrong number of values for its Number is an error (exit 1).
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf split calls.vcf.gz -o calls.split.vcf.gz
  gloci vcf split calls.bcf --types snps -o snps-split.bcf
  gloci vcf split calls.vcf.gz | gloci subtract -a - -b known.vcf.gz --match alleles")]
pub struct VcfSplitArgs {
    /// Input VCF or BCF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Split only records with an ALT of this type: snps (a record with a
    /// SNP ALT), indels (an insertion or deletion ALT), both or any (every
    /// multi-allelic record; bcftools splits both and any alike). A split
    /// record is split into all its ALTs; other records are written as they
    /// are.
    #[arg(long, value_enum, default_value = "any", value_name = "TYPES")]
    pub types: SplitTypesArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

// P3b Task 21

/// `--check-ref` modes (bcftools `-c`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum CheckRefArg {
    /// Fail on the first disagreement (the default).
    #[value(name = "e")]
    E,
    /// Warn and write the record as it is.
    #[value(name = "w")]
    W,
    /// Drop the record.
    #[value(name = "x")]
    X,
    /// Set REF from the FASTA (an ALT that equals it is swapped into REF).
    #[value(name = "s")]
    S,
}

/// `--rm-dup` types (bcftools `-d`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum RmDupArg {
    Snps,
    Indels,
    Both,
    All,
    Exact,
}

/// Left-align and normalise records against a reference (bcftools norm -f).
///
/// Every indel, and every record whose REF and ALT share a first or last
/// base, is trimmed and moved to its leftmost position, as `bcftools
/// norm -f REF` does without -a: the bases all alleles end with are
/// removed (the alleles grow to the left from the reference when one
/// would be empty), then the bases they all start with. Records that
/// move are re-sorted within --window bases and written in position
/// order; INFO/END follows the new REF. REF is checked against the FASTA
/// first (--check-ref). The FASTA needs a .fai (and a .gzi for BGZF):
/// run `gloci index REF` once; gloci never builds one itself. Only
/// windows of the reference are read, never whole chromosomes. A record
/// with any symbolic ALT (`<DEL>`, `<DUP>`, ...) is not trimmed or moved.
/// --skip-invalid also skips records that fail --check-ref e. `bcftools
/// norm -a`, `-m+`, `--old-rec-tag` and `-N` are not available.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf norm calls.vcf.gz -f ref.fa -o calls.norm.vcf.gz
  gloci vcf norm calls.bcf --fasta ref.fa --split -o calls.norm.bcf
  gloci vcf norm calls.vcf.gz -f ref.fa --check-ref x --rm-dup exact | gloci vcf stats -
  gloci vcf norm calls.vcf.gz -f ref.fa --check-ref s -o fixed.vcf.gz")]
pub struct VcfNormArgs {
    /// Input VCF or BCF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Reference FASTA, with its .fai (and .gzi if BGZF-compressed); create
    /// them with `gloci index REF`.
    #[arg(short = 'f', long, visible_alias = "fasta-ref", value_name = "REF")]
    pub fasta: PathBuf,
    /// Split multi-allelic records first, as `bcftools norm -m-TYPES -f`:
    /// snps, indels, both or any (the default; bcftools splits both and any
    /// alike). The same split as `gloci vcf split`.
    #[arg(
        long,
        value_enum,
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "any",
        value_name = "TYPES"
    )]
    pub split: Option<SplitTypesArg>,
    /// What to do with a record whose REF is not the FASTA's (or has a base
    /// other than ACGTN, or an ALT equal to REF): e (fail, the default), w
    /// (warn, write the record as it is), x (drop it) or s (set REF from the
    /// FASTA; an ALT that equals the FASTA is swapped into REF and GT is
    /// recoded).
    #[arg(
        short = 'c',
        long,
        value_enum,
        default_value = "e",
        value_name = "MODE"
    )]
    pub check_ref: CheckRefArg,
    /// Remove duplicate records after normalising, keeping the first: snps
    /// (a SNP or MNP record at the position of a kept one), indels, both
    /// (each type on its own), all (any record at a kept record's position)
    /// or exact (the same REF and ALT alleles).
    #[arg(short = 'd', long, value_enum, value_name = "TYPES")]
    pub rm_dup: Option<RmDupArg>,
    /// Records that move left are re-sorted within this many bases (bcftools
    /// --site-win). A window smaller than the largest shift can order records
    /// differently from bcftools; the default 1000 matches it.
    #[arg(short = 'w', long, default_value = "1000", value_name = "BASES", value_parser = clap::value_parser!(u64).range(1..))]
    pub window: u64,
    // bcftools norm options gloci does not offer: each is a usage error that
    // names the option.
    #[arg(short = 'a', long = "atomize", hide = true)]
    pub atomize: bool,
    #[arg(
        short = 'm',
        long = "multiallelics",
        hide = true,
        value_name = "TYPE",
        allow_hyphen_values = true
    )]
    pub multiallelics: Option<String>,
    #[arg(long = "old-rec-tag", hide = true, value_name = "STR")]
    pub old_rec_tag: Option<String>,
    #[arg(short = 'N', long = "do-not-normalize", hide = true)]
    pub do_not_normalize: bool,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

// P3b Task 22

/// Transfer, remove and add fields of a VCF/BCF (bcftools annotate).
///
/// --annotations SRC names a VCF/BCF, a tab-delimited file (plain or
/// BGZF) or a BED file; --columns says what to take from it. FILE and SRC
/// must both be sorted in one chromosome order (FILE's header contigs,
/// else --chrom-order); both are read as a stream, no index is needed
/// and neither is held in memory. A record of a VCF SRC annotates the
/// FILE record at the same CHROM and POS that shares an allele with it
/// (bcftools' default --pair-logic); a tab file with REF and ALT columns
/// matches the same way, one with CHROM,POS matches the records whose
/// REF span covers the position, and a BED file or CHROM,FROM,TO columns
/// match overlapping records (--min-overlap restricts them). --columns
/// takes ID, QUAL, FILTER, INFO/TAG (or TAG), NEW:=OLD, INFO (every tag
/// of a VCF SRC) and ^INFO/TAG lists; +TAG fills only missing values;
/// `.TAG` carries a tab file's missing values over; `-` skips a tab
/// column. A tag the FILE header lacks is declared from a VCF SRC, or
/// must come from --header-lines. --remove drops ID, QUAL, FILTER, INFO,
/// FORMAT (GT stays), INFO/TAG, FORMAT/TAG or FILTER/NAME, and their
/// header lines; ^ keeps only the listed tags of its kind. --set-id
/// writes IDs from %CHROM, %POS, %REF, %FIRST_ALT, %ID and %INFO/TAG
/// (a + prefix fills only missing IDs). --mark-sites +TAG flags the
/// records that have a match, -TAG those that do not. INFO values are
/// copied as written, not reformatted; Number=A and R values follow
/// their allele, not their position. Not available: -s/-S, -i/-e,
/// --rename-chrs, --merge-logic, --single-overlaps, --rename-annots,
/// FORMAT transfers.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf annotate calls.vcf.gz -a dbsnp.vcf.gz -c ID,INFO/RS -o calls.annot.vcf.gz
  gloci vcf annotate calls.vcf.gz -a genes.bed.gz -c CHROM,FROM,TO,GENE --header-lines genes.hdr
  gloci vcf annotate calls.vcf.gz -x INFO,FORMAT/PL --mark-sites +CLEAN
  gloci vcf annotate calls.vcf.gz -I '+%CHROM:%POS:%REF:%FIRST_ALT' -o ids.vcf.gz")]
pub struct VcfAnnotateArgs {
    /// Input VCF or BCF (`-` for stdin), sorted.
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The annotation source: a VCF/BCF, a tab-delimited file (plain or
    /// BGZF) or a BED file, sorted in FILE's chromosome order. A tab file
    /// needs --columns; a BED file is a tab file named *.bed, *.bed.gz or
    /// *.bed.bgz (0-based starts).
    #[arg(short = 'a', long, value_name = "SRC")]
    pub annotations: Option<String>,
    /// What to take from --annotations, comma-separated: ID, QUAL, FILTER,
    /// INFO/TAG or TAG, NEW:=OLD (a VCF source's OLD stored as NEW), INFO
    /// (every tag of a VCF source), ^INFO/TAG (INFO without it), +TAG (only
    /// fill missing values), .TAG (carry a tab file's missing values) and
    /// `-` (skip a tab column). A tab file's list also names its columns:
    /// CHROM, POS (or FROM and TO) and optionally REF and ALT, in file order,
    /// e.g. CHROM,POS,REF,ALT,ID,INFO/AF.
    #[arg(short = 'c', long, value_name = "LIST", allow_hyphen_values = true)]
    pub columns: Option<String>,
    /// Remove, comma-separated: ID, QUAL, FILTER, INFO, FORMAT (GT is kept),
    /// INFO/TAG, FORMAT/TAG, FILTER/NAME, or ##key for a header line. A
    /// leading ^ (^INFO/TAG, ^FORMAT/TAG, ^FILTER/NAME) keeps only the
    /// listed tags of that kind and removes the rest. Their header lines go
    /// too. Removals run before anything is transferred.
    #[arg(short = 'x', long, value_name = "LIST")]
    pub remove: Option<String>,
    /// Append the header lines of FILE (## lines such as ##INFO=<...>) to the
    /// header; they can declare a tag --columns then uses.
    #[arg(long, value_name = "FILE")]
    pub header_lines: Option<PathBuf>,
    /// Write IDs from a format of %CHROM, %POS, %REF, %FIRST_ALT, %ID and
    /// %INFO/TAG; a leading + sets only missing IDs, e.g.
    /// '+%CHROM:%POS:%REF:%FIRST_ALT'. A tag name is letters, digits, _ and
    /// . only: write %POS:%REF, not %POS_%REF.
    #[arg(
        short = 'I',
        long,
        value_name = "[+]FORMAT",
        allow_hyphen_values = true
    )]
    pub set_id: Option<String>,
    /// Add INFO/TAG as a flag to the records that have a match in
    /// --annotations (+TAG, or TAG) or that have none (-TAG); without
    /// --annotations every record gets it.
    #[arg(short = 'm', long, value_name = "[+-]TAG", allow_hyphen_values = true)]
    pub mark_sites: Option<String>,
    /// For a BED or FROM,TO annotation file: the fraction of the annotation
    /// (ANN), of the record (:VCF) or of both (ANN:VCF) that must overlap.
    #[arg(long, value_name = "ANN[:VCF]")]
    pub min_overlap: Option<String>,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    // bcftools annotate options gloci does not offer: each is a usage error
    // that names the option.
    #[arg(long = "rename-chrs", hide = true, value_name = "FILE")]
    pub rename_chrs: Option<String>,
    #[arg(short = 's', long = "samples", hide = true, value_name = "LIST")]
    pub samples: Option<String>,
    #[arg(short = 'S', long = "samples-file", hide = true, value_name = "FILE")]
    pub samples_file: Option<String>,
    #[arg(short = 'i', long = "include", hide = true, value_name = "EXPR")]
    pub include: Option<String>,
    #[arg(short = 'e', long = "exclude", hide = true, value_name = "EXPR")]
    pub exclude: Option<String>,
    #[arg(short = 'l', long = "merge-logic", hide = true, value_name = "LOGIC")]
    pub merge_logic: Option<String>,
    #[arg(long = "single-overlaps", hide = true)]
    pub single_overlaps: bool,
    #[arg(long = "rename-annots", hide = true, value_name = "FILE")]
    pub rename_annots: Option<String>,
    #[arg(long = "pair-logic", hide = true, value_name = "LOGIC")]
    pub pair_logic: Option<String>,
    #[arg(short = 'C', long = "columns-file", hide = true, value_name = "FILE")]
    pub columns_file: Option<String>,
    #[arg(short = 'H', long = "header-line", hide = true, value_name = "STR")]
    pub header_line: Option<String>,
    #[arg(short = 'k', long = "keep-sites", hide = true)]
    pub keep_sites: bool,
    #[arg(short = 'r', long = "regions", hide = true, value_name = "REGION")]
    pub regions: Option<String>,
    #[arg(short = 'R', long = "regions-file", hide = true, value_name = "FILE")]
    pub regions_file: Option<String>,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

// P3b Task 23

/// `bcftools merge -m`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum MergeModeArg {
    None,
    Snps,
    Indels,
    Both,
    All,
    Id,
}

/// `bcftools merge -F`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum FilterLogicArg {
    /// Records are merged with all their filters (`+`).
    #[value(name = "+")]
    All,
    /// A PASS in any input removes the filters (`x`).
    #[value(name = "x")]
    Pass,
}

/// Combine the samples of several VCF/BCF files (bcftools merge).
///
/// The inputs must be sorted in one chromosome order (the first input's
/// ##contig lines, else --chrom-order); they are read as a stream, so no
/// index is needed and only the records of one position are held. Records
/// at the same CHROM and POS are merged as --merge says: ALT alleles are
/// united in input order (a shorter REF is extended to the longest), and
/// GT and Number=A, R and G FORMAT values are re-indexed to the merged
/// alleles. A sample absent from a record is ./. (or 0/0 with -0) in the
/// record's ploidy. QUAL is the largest, FILTER the union without PASS,
/// ID the distinct IDs joined with ;. Sample names must be distinct
/// (--force-samples renames a repeat to N:NAME). The header is the first
/// input's with the other inputs' new lines added. The inputs are open
/// together: merging more files than `ulimit -n` allows fails with the
/// operating system's error; raise the limit (ulimit -n 4096). bcftools
/// merge -g, -L, -r, -R, -T, --print-header and --use-header are not
/// available.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf merge a.vcf.gz b.vcf.gz c.vcf.gz -o cohort.vcf.gz
  gloci vcf merge --file-list shards.txt -m none --force-samples -o cohort.bcf
  gloci vcf merge a.bcf b.bcf -i DP:sum,AF:avg -0 | gloci vcf stats - --out-prefix qc
  gloci vcf merge a.vcf.gz b.vcf.gz -f PASS,. -m id -o by-id.vcf.gz")]
pub struct VcfMergeArgs {
    /// Input VCF or BCF files (two or more; `-` is not allowed), each
    /// sorted in the same chromosome order. Or use --file-list.
    #[arg(value_name = "FILE")]
    pub files: Vec<String>,
    /// Read the input file names from FILE, one per line (bcftools -l).
    #[arg(short = 'l', long, value_name = "FILE")]
    pub file_list: Option<PathBuf>,
    /// Which records at one position are merged into one (bcftools -m):
    /// none (only records with the same alleles), snps, indels, both (SNPs
    /// together and indels together; the default), all (any records) or id
    /// (records with the same ID).
    #[arg(
        short = 'm',
        long,
        value_enum,
        default_value = "both",
        value_name = "MODE"
    )]
    pub merge: MergeModeArg,
    /// Rename a repeated sample name to N:NAME (N is the number of its
    /// file) instead of failing.
    #[arg(long)]
    pub force_samples: bool,
    /// How INFO tags are merged, TAG:METHOD with METHOD one of sum, avg,
    /// min, max, join, comma-separated; `-` for none. The default is
    /// DP:sum,DP4:sum (and AC:sum,AN:sum when there are no samples; with
    /// samples AC and AN are recomputed from the merged genotypes). A
    /// tag with no rule takes the first record's value; Number=A and R tags
    /// are re-indexed to the merged alleles.
    #[arg(short = 'i', long, value_name = "RULES", allow_hyphen_values = true)]
    pub info_rules: Option<String>,
    /// Fill a sample absent from a record with 0/0 instead of ./. (bcftools -0).
    #[arg(short = '0', long)]
    pub missing_to_ref: bool,
    /// Skip, before merging, the records whose FILTER has none of these
    /// names (comma-separated; `.` is a missing FILTER, PASS is PASS).
    #[arg(short = 'f', long, value_name = "LIST")]
    pub apply_filters: Option<String>,
    /// What to do with the FILTER of merged records: + (all filters kept,
    /// PASS dropped when another filter is set; the default) or x (PASS in
    /// any input makes the result PASS).
    #[arg(
        short = 'F',
        long,
        value_enum,
        default_value = "+",
        value_name = "LOGIC"
    )]
    pub filter_logic: FilterLogicArg,
    /// Accepted and ignored: gloci merges sorted inputs as a stream and
    /// never needs an index (bcftools merge needs one unless --no-index).
    #[arg(long, alias = "force-no-index", hide = true)]
    pub no_index: bool,
    /// Accepted: merge a single input (bcftools merge --force-single).
    #[arg(long, hide = true)]
    pub force_single: bool,
    // bcftools merge options gloci does not offer: each is a usage error
    // that names the option.
    #[arg(short = 'g', long = "gvcf", hide = true, value_name = "REF")]
    pub gvcf: Option<String>,
    #[arg(short = 'L', long = "local-alleles", hide = true, value_name = "N")]
    pub local_alleles: Option<String>,
    #[arg(long = "print-header", hide = true)]
    pub print_header: bool,
    #[arg(long = "use-header", hide = true, value_name = "FILE")]
    pub use_header: Option<String>,
    #[arg(short = 'r', long = "regions", hide = true, value_name = "REGION")]
    pub regions: Option<String>,
    #[arg(short = 'R', long = "regions-file", hide = true, value_name = "FILE")]
    pub regions_file: Option<String>,
    #[arg(short = 'T', long = "targets-file", hide = true, value_name = "FILE")]
    pub targets_file: Option<String>,
    #[arg(long = "merge-logic", hide = true, value_name = "LOGIC")]
    pub merge_logic: Option<String>,
    #[arg(short = 'M', long = "missing-rules", hide = true, value_name = "RULES")]
    pub missing_rules: Option<String>,
    #[arg(long = "regions-overlap", hide = true, value_name = "0|1|2")]
    pub regions_overlap: Option<String>,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}
