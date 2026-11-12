"""pybedtools: A records overlapping any B record. Usage: IN_A IN_B OUT."""

import os
import sys

import pybedtools

a, b, out = sys.argv[1:4]
# Its temporary files go next to the output, not into a small /tmp.
pybedtools.set_tempdir(os.path.dirname(os.path.abspath(out)))
pybedtools.BedTool(a).intersect(pybedtools.BedTool(b), u=True, sorted=True).saveas(out)
pybedtools.cleanup()
