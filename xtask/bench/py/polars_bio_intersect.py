"""polars-bio: A records overlapping any B record. Usage: IN_A IN_B OUT."""

import sys

import polars_bio as pb

a, b, out = sys.argv[1:4]
df = pb.overlap(
    pb.scan_bed(a, use_zero_based=True),
    pb.scan_bed(b, use_zero_based=True),
    overlap_output="left",
    distinct_output=True,
    output_type="polars.DataFrame",
)
df.write_csv(out, separator="\t", include_header=False)
