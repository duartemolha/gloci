"""polars-bio: merge overlapping and bookended records. Usage: IN_A OUT."""

import sys

import polars_bio as pb

a, out = sys.argv[1:3]
# polars-bio counts a gap in positions: min_dist=1 merges bookended records, as bedtools does.
df = pb.merge(pb.scan_bed(a, use_zero_based=True), min_dist=1, output_type="polars.DataFrame")
df.select(["chrom", "start", "end"]).write_csv(out, separator="\t", include_header=False)
