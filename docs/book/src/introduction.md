# gloci

gloci is a single fast binary for genomic interval, VCF, alignment and sequence work. It replaces the
everyday parts of bedtools, bedops, vcftools, tabix and bgzip, and common samtools commands, and it is
built for multi-GB VCF, BED and BAM files. The same commands are available from Python as lazy,
streaming pipelines.

Why another tool:

- **One set of rules.** Coordinates are 0-based half-open `u64` everywhere, a VCF record is the bases
  its alleles change, and a zero-length interval is a point. [Semantics](semantics.md) is the contract;
  a rule never changes without that page changing first.
- **Same answers as the tools it replaces, except where documented.** Output is checked against pinned
  bedtools, bedops, vcftools, htslib and bcftools releases through committed golden files. Every
  intentional difference has a reason and a test: see [Differences](differences.md).
- **Streaming and parallel.** Sorted inputs are swept in a single pass; output is byte-identical at every
  `--threads` value, and files are written atomically.
- **Errors that help.** A failure says what went wrong, where, and how to fix it. Exit codes are 0 (success), 1
  (data or runtime error) and 2 (usage error).
- **One static binary.** No runtime dependencies on Linux (musl) or macOS, and no htslib.

Where to go next:

- New to gloci: [Install](install.md), then the [Quickstart](quickstart.md).
- Coming from bedtools, bedops, vcftools or samtools: [Migrating](migration.md) maps each command and
  flag.
- Looking for one command: the [command reference](reference/sort.md) has the usage text and examples of
  every command.
- Using Python: [gloci for Python](python.md) and the [Python API](python-api.md).
