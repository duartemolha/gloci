"""The exceptions gloci raises (spec 8.5, D90).

`str(error)` is the command line's message without its `error: ` prefix. Exceptions raised by
Python code that gloci calls (a source's iterator, a log handler) propagate unchanged.
"""

from __future__ import annotations

__all__ = ["DataError", "GlociError", "UnsortedError", "UsageError"]


class GlociError(Exception):
    """The base class of gloci's errors. Raised on its own for an internal error, which is a
    bug: its message says so."""


class UsageError(GlociError, ValueError):
    """A call gloci cannot run as given: an unknown keyword argument, a value of the wrong
    kind, a missing option. Raised when the step is called, or by a terminal method before
    anything runs; messages name keyword arguments (`max_gap=`), not flags."""


class DataError(GlociError):
    """An input or output that cannot be used: a malformed record, a file that cannot be
    opened, a truncated file, or data the command cannot handle.

    Attributes:
        path: The input or output the error is about, as it was given (an in-memory input by
            its display name, such as `<intersect output>`), or None when it names none.
        line: The 1-based line of the error (in a binary file, its record), or None.
    """

    def __init__(self, message: str, path: str | None = None, line: int | None = None) -> None:
        super().__init__(message)
        self.path = path
        self.line = line

    def __reduce__(self):  # type: ignore[no-untyped-def]
        return (type(self), (str(self), self.path, self.line))


class UnsortedError(DataError):
    """An input that is not sorted. `path` names it and `line` is the later record's line.

    Attributes:
        records: Both records, as (earlier, later), each its text.
        lines: Their lines (in a binary file, their records), as (earlier, later).
        hint: How to fix it: sort the input with `gloci sort`, or pass `sort=True`.
    """

    def __init__(
        self,
        message: str,
        path: str | None = None,
        line: int | None = None,
        records: tuple[str, str] | None = None,
        lines: tuple[int, int] | None = None,
        hint: str | None = None,
    ) -> None:
        super().__init__(message, path, line)
        self.records = records
        self.lines = lines
        self.hint = hint

    def __reduce__(self):  # type: ignore[no-untyped-def]
        args = (str(self), self.path, self.line, self.records, self.lines, self.hint)
        return (type(self), args)
