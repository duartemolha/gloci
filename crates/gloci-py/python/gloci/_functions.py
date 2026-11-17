"""D87's module functions: the sources `windows()` and `random()`, which start a pipeline with
no input, and the functions that write files at once, as the command line does
(`vcf_stats()`, `vcf_extract_matrix()`, `vcf_compare()`, `vcf_popgen()`, `index()`).

Keyword arguments are the command's flags with `-` written as `_` (D86), checked before
anything runs. The functions run through the same driver as a pipeline: the GIL is released,
Ctrl-C cancels, and errors are gloci's exceptions.
"""

from __future__ import annotations

import os
import pathlib
from typing import Any

from . import _native
from ._pipeline import (
    Path,
    Pipeline,
    Source,
    _Command,
    _kwarg,
    _placeholders,
    _Plan,
    classify_source,
)
from ._signature import SIGNATURES, make_doc, make_signature
from .errors import UsageError

__all__ = [
    "index",
    "random",
    "vcf_compare",
    "vcf_extract_matrix",
    "vcf_popgen",
    "vcf_stats",
    "depth_files",
    "windows",
]

#: Each module function's command, by name (`_native.module_functions()`).
_FUNCTIONS: dict[str, dict[str, Any]] = {
    f["function"]: dict(f, method=f["function"]) for f in _native.module_functions()
}
if set(_FUNCTIONS) != set(__all__):
    raise RuntimeError(
        f"gloci's module functions {sorted(_FUNCTIONS)} are not those of gloci-py {__all__}"
    )


def _source(name: str, kwargs: dict[str, Any]) -> Pipeline:
    info = _FUNCTIONS[name]
    kw = {k: _kwarg(f"{name}()", info["path"], k, v) for k, v in kwargs.items()}
    # The user's values, checked at once (D85).
    _native.source_argv(info["path"], _placeholders(kw), tokens=False)
    return Pipeline._make(_Command((info, kw), f"<{name}>"), None, ())


def windows(**kwargs: Any) -> Pipeline:
    # The docstring is made from the command's spec (below).
    return _source("windows", kwargs)


def random(**kwargs: Any) -> Pipeline:
    # The docstring is made from the command's spec (below).
    return _source("random", kwargs)


def _prefix(where: str, out: object) -> str:
    if not isinstance(out, (str, os.PathLike)):
        raise TypeError(f"{where}: out= takes a path prefix (str or os.PathLike)")
    return os.fsdecode(os.fspath(out))


def _run(name: str, sources: list[Source], kwargs: dict[str, Any]) -> list[pathlib.Path]:
    """Run the module function `name` on `sources`, and return the files it wrote, sorted: the
    output files the run itself committed, as the command named them, whatever else changed in
    the directory meanwhile."""
    info = _FUNCTIONS[name]
    where = f"{name}()"
    kw = {k: _kwarg(where, info["path"], k, v) for k, v in kwargs.items()}
    # The user's values, checked before anything runs (D85).
    _native.function_argv(info["path"], ["-"] * len(sources), _placeholders(kw), tokens=False)
    plan = _Plan()
    inputs = [plan.input(s) for s in sources]
    real_kw, shown_kw = plan.kwargs(kw)
    argv = _native.function_argv(info["path"], [r for r, _ in inputs], real_kw)
    shown = _native.function_argv(info["path"], [s for _, s in inputs], shown_kw, check=False)
    plan.add(argv, shown)
    written: list[str] = plan.run_plan("save")
    return sorted({pathlib.Path(p) for p in written})


def _with_prefix(where: str, out: object, kwargs: dict[str, Any]) -> dict[str, Any]:
    if "out_prefix" in kwargs:
        raise UsageError(f"{where}: out= and out_prefix= both set the prefix; pass out= only")
    return {**kwargs, "out_prefix": _prefix(where, out)}


def vcf_stats(src: object, *, out: str | os.PathLike[str], **kwargs: Any) -> list[pathlib.Path]:
    # The docstring is made from the command's spec (below).
    kw = _with_prefix("vcf_stats()", out, kwargs)
    return _run("vcf_stats", [classify_source(src, "vcf_stats()")], kw)


def depth_files(src: object, *, out: str | os.PathLike[str], **kwargs: Any) -> list[pathlib.Path]:
    # The docstring is made from the command's spec (below).
    kw = _with_prefix("depth_files()", out, kwargs)
    return _run("depth_files", [classify_source(src, "depth_files()")], kw)


def vcf_extract_matrix(
    src: object, *, out: str | os.PathLike[str], **kwargs: Any
) -> list[pathlib.Path]:
    # The docstring is made from the command's spec (below).
    kw = _with_prefix("vcf_extract_matrix()", out, kwargs)
    kw.setdefault("genotype_matrix", "012")
    src = classify_source(src, "vcf_extract_matrix()")
    return _run("vcf_extract_matrix", [src], kw)


def vcf_compare(
    a: object, b: object, *, out: str | os.PathLike[str], **kwargs: Any
) -> list[pathlib.Path]:
    # The docstring is made from the command's spec (below).
    kw = _with_prefix("vcf_compare()", out, kwargs)
    sources = [classify_source(a, "vcf_compare()"), classify_source(b, "vcf_compare()")]
    return _run("vcf_compare", sources, kw)


def vcf_popgen(src: object, *, out: str | os.PathLike[str], **kwargs: Any) -> list[pathlib.Path]:
    # The docstring is made from the command's spec (below).
    kw = _with_prefix("vcf_popgen()", out, kwargs)
    return _run("vcf_popgen", [classify_source(src, "vcf_popgen()")], kw)


def index(path: str | os.PathLike[str], **kwargs: Any) -> list[pathlib.Path]:
    # The docstring is made from the command's spec (below).
    src = classify_source(path, "index()")
    if not isinstance(src, Path):
        shown = "a Pipeline" if isinstance(path, Pipeline) else src.display
        raise UsageError(
            f"index needs a file, not {shown}: an index is written next to the file it indexes, "
            "so save() the data to a file and index that"
        )
    return _run("index", [src], kwargs)


# The signatures and docstrings come from the command specs, as the type stub's do.
for _name, _info in _FUNCTIONS.items():
    _f = globals()[_name]
    _f.__doc__ = make_doc(tuple(_info["path"]), SIGNATURES[_name])
    _f.__signature__ = make_signature(tuple(_info["path"]), SIGNATURES[_name], method=False)
del _name, _info, _f
