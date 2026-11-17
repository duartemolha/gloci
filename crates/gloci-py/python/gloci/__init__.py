"""gloci: genomic intervals and VCF as lazy, streaming pipelines, run in-process.

    import gloci
    gloci.open("calls.vcf.gz").intersect(b="targets.bed").merge(max_gap=50).save("out.bed")

Each step is a gloci command (`gloci intersect`, `gloci merge`), run in the same code as the
command line, so the output equals the shell pipeline that `Pipeline.cli()` prints. A source is
a path, a file object, an iterable of lines or tuples, an Arrow table (pyarrow, polars,
pandas) or another pipeline; `gloci.windows()` and `gloci.random()` start a pipeline with no
input, and `vcf_stats()`, `vcf_extract_matrix()`, `vcf_compare()`, `vcf_popgen()` and `index()`
write files at once. A pipeline's output comes back as bytes, a count, `gloci.Record`s (iteration), or an
Arrow table (`to_arrow_stream()`, `to_arrow()`, `to_polars()`, `to_pandas()`).
"""

from ._functions import (
    depth_files,
    index,
    random,
    vcf_compare,
    vcf_extract_matrix,
    vcf_popgen,
    vcf_stats,
    windows,
)
from ._native import ArrowStream, Record, RecordIterator, __version__
from ._pipeline import Pipeline, open
from ._types import Source, SupportsArrowStream
from .errors import DataError, GlociError, UnsortedError, UsageError

__all__ = [
    "ArrowStream",
    "DataError",
    "GlociError",
    "Pipeline",
    "Record",
    "RecordIterator",
    "Source",
    "SupportsArrowStream",
    "UnsortedError",
    "UsageError",
    "__version__",
    "depth_files",
    "index",
    "open",
    "random",
    "vcf_compare",
    "vcf_extract_matrix",
    "vcf_popgen",
    "vcf_stats",
    "windows",
]
