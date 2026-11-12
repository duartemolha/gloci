"""PyRanges v1: merge overlapping and bookended records. Usage: IN_A OUT."""

import sys
from pathlib import Path

import pyranges1 as pr

a, out = sys.argv[1:3]
pr.read_bed(Path(a)).merge_overlaps(use_strand=False, slack=1).to_bed(out)
