"""PyRanges v1: A records overlapping any B record. Usage: IN_A IN_B OUT."""

import sys
from pathlib import Path

import pyranges1 as pr

a, b, out = sys.argv[1:4]
pr.read_bed(Path(a)).overlap(pr.read_bed(Path(b)), strand_behavior="ignore").to_bed(out)
