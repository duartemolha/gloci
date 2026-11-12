"""The types of the public API: what `gloci.Source` stands for, as the type stub declares it."""

from __future__ import annotations

import os
from collections.abc import Iterable
from typing import IO, Protocol, Union

from ._native import Record
from ._pipeline import Pipeline

__all__ = ["Source", "SupportsArrowStream"]


class SupportsArrowStream(Protocol):
    """An object that exports an Arrow C stream: a pyarrow Table or RecordBatchReader, a polars
    or pandas DataFrame, or a `gloci.ArrowStream`."""

    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...


#: What a pipeline reads: a path, another Pipeline, a file object, an iterable of lines or
#: tuples of values, or an Arrow table or stream.
Source = Union[
    str,
    os.PathLike[str],
    Pipeline,
    IO[str],
    IO[bytes],
    Iterable[Union[str, bytes, tuple[object, ...], Record]],
    SupportsArrowStream,
]
