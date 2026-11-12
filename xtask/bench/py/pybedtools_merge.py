"""pybedtools: merge overlapping and bookended records. Usage: IN_A OUT."""

import os
import sys

import pybedtools

a, out = sys.argv[1:3]
pybedtools.set_tempdir(os.path.dirname(os.path.abspath(out)))
pybedtools.BedTool(a).merge().saveas(out)
pybedtools.cleanup()
