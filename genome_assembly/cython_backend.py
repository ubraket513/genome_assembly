"""Optional Cython accelerator bridge."""

from __future__ import annotations

from importlib import import_module
from importlib.util import find_spec
from typing import Iterable

_REQUIRED_SYMBOLS = ("count_kmers", "build_edges", "compact_contigs")


def cython_available() -> bool:
    """Return True when the compiled Cython backend is importable and current."""

    if find_spec("genome_assembly._cython_backend") is None:
        return False
    try:
        module = import_module("genome_assembly._cython_backend")
    except ImportError:
        return False
    return all(hasattr(module, symbol) for symbol in _REQUIRED_SYMBOLS)


def require_cython_backend():
    """Import the compiled Cython backend or raise an actionable error."""

    try:
        module = import_module("genome_assembly._cython_backend")
    except ImportError as exc:
        raise RuntimeError(
            "Cython backend requested, but genome_assembly._cython_backend is not built. "
            "Build it with: python -m pip install '.[perf]' && python setup.py build_ext --inplace"
        ) from exc

    missing = [symbol for symbol in _REQUIRED_SYMBOLS if not hasattr(module, symbol)]
    if missing:
        raise RuntimeError(
            "Cython backend requested, but genome_assembly._cython_backend is outdated "
            f"(missing: {', '.join(missing)}). Rebuild it with: python setup.py build_ext --inplace"
        )

    return module


def count_kmers(
    reads: Iterable[str],
    k: int,
    *,
    skip_ambiguous: bool = True,
) -> list[tuple[str, int]]:
    """Count k-mers with the optional Cython extension."""

    if k < 1:
        raise ValueError("k must be >= 1")
    module = require_cython_backend()
    return module.count_kmers(list(reads), k, skip_ambiguous)


def build_edges(
    reads: Iterable[str],
    node_k: int,
    *,
    min_abundance: int = 1,
    skip_ambiguous: bool = True,
) -> tuple[int, list[tuple[str, str, str, int]]]:
    """Build a deterministic de Bruijn edge table with Cython."""

    if node_k < 1:
        raise ValueError("node_k must be >= 1")
    if min_abundance < 1:
        raise ValueError("min_abundance must be >= 1")
    module = require_cython_backend()
    raw_edge_count, edges = module.build_edges(
        list(reads),
        node_k,
        min_abundance,
        skip_ambiguous,
    )
    return raw_edge_count, edges


def compact_contigs(
    node_k: int,
    edge_rows: Iterable[tuple[str, str, str, int]],
    *,
    min_length: int = 0,
) -> list[tuple[str, float, int]]:
    """Compact de Bruijn graph edges into contig rows with Cython."""

    if node_k < 1:
        raise ValueError("node_k must be >= 1")
    if min_length < 0:
        raise ValueError("min_length must be >= 0")
    module = require_cython_backend()
    return module.compact_contigs(node_k, list(edge_rows), min_length)
