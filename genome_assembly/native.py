"""Optional native backend bridge."""

from __future__ import annotations

from importlib import import_module
from importlib.util import find_spec
from typing import Iterable

_REQUIRED_SYMBOLS = ("count_kmers", "build_edges")


def native_available() -> bool:
    """Return True when the current Rust extension module is importable."""

    if find_spec("genome_assembly_native") is None:
        return False
    try:
        native = import_module("genome_assembly_native")
    except ImportError:
        return False
    return all(hasattr(native, symbol) for symbol in _REQUIRED_SYMBOLS)


def require_native():
    """Import the native extension or raise an actionable error."""

    try:
        native = import_module("genome_assembly_native")
    except ImportError as exc:
        raise RuntimeError(
            "Native backend requested, but genome_assembly_native is not installed. "
            "Build it with: cd native/genome_assembly_native && maturin develop --release"
        ) from exc

    missing = [symbol for symbol in _REQUIRED_SYMBOLS if not hasattr(native, symbol)]
    if missing:
        raise RuntimeError(
            "Native backend requested, but genome_assembly_native is outdated "
            f"(missing: {', '.join(missing)}). Rebuild it with: "
            "cd native/genome_assembly_native && maturin develop --release"
        )

    return native


def count_kmers(
    reads: Iterable[str],
    k: int,
    *,
    skip_ambiguous: bool = True,
) -> list[tuple[str, int]]:
    """Count k-mers with the optional Rust extension."""

    if k < 1:
        raise ValueError("k must be >= 1")
    native = require_native()
    return native.count_kmers(list(reads), k, skip_ambiguous)


def build_edges(
    reads: Iterable[str],
    node_k: int,
    *,
    min_abundance: int = 1,
    skip_ambiguous: bool = True,
) -> tuple[int, list[tuple[str, str, str, int]]]:
    """Build a deterministic de Bruijn edge table with the optional Rust extension."""

    if node_k < 1:
        raise ValueError("node_k must be >= 1")
    if min_abundance < 1:
        raise ValueError("min_abundance must be >= 1")
    native = require_native()
    raw_edge_count, edges = native.build_edges(
        list(reads),
        node_k,
        min_abundance,
        skip_ambiguous,
    )
    return raw_edge_count, edges
