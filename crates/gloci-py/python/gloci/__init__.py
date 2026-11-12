"""gloci: genomic intervals and VCF as lazy, streaming pipelines, run in-process.

    import gloci
    gloci.open("calls.vcf.gz").intersect(b="targets.bed").merge(max_gap=50).save("out.bed")

Each step is a gloci command (`gloci intersect`, `gloci merge`), run in the same code as the
command line, so the output equals the shell pipeline that `Pipeline.cli()` prints.
"""

from ._native import __version__
from ._pipeline import Pipeline, open
from .errors import DataError, GlociError, UnsortedError, UsageError

__all__ = [
    "DataError",
    "GlociError",
    "Pipeline",
    "UnsortedError",
    "UsageError",
    "__version__",
    "open",
]
