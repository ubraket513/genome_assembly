"""Assembly quality metrics."""

from __future__ import annotations

from typing import Iterable


def _lengths(contigs: Iterable[str]) -> list[int]:
    return sorted((len(contig) for contig in contigs), reverse=True)


def nx(contigs: Iterable[str], x: int | float, *, reference_length: int | None = None) -> int:
    """Return the Nx contig length.

    If reference_length is omitted, the total assembled length is used.
    """

    if not 0 < float(x) <= 100:
        raise ValueError("x must be in the interval (0, 100]")

    contig_list = list(contigs)
    lengths = _lengths(contig_list)
    if not lengths:
        return 0

    total = reference_length if reference_length is not None else sum(lengths)
    if total <= 0:
        return 0

    threshold = total * float(x) / 100.0
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def n50(contigs: Iterable[str], *, reference_length: int | None = None) -> int:
    return nx(contigs, 50, reference_length=reference_length)


def assembly_stats(contigs: Iterable[str], *, reference_length: int | None = None) -> dict[str, int | float | None]:
    contig_list = list(contigs)
    lengths = _lengths(contig_list)
    total_bp = sum(lengths)
    denominator = reference_length if reference_length is not None else total_bp
    gc_bases = 0

    for contig in contig_list:
        sequence = contig.upper()
        gc_bases += sequence.count("G") + sequence.count("C")

    return {
        "contigs": len(lengths),
        "total_bp": total_bp,
        "max_contig": lengths[0] if lengths else 0,
        "n50": nx(contig_list, 50),
        "n90": nx(contig_list, 90),
        "ng50": nx(contig_list, 50, reference_length=reference_length) if reference_length else None,
        "ng90": nx(contig_list, 90, reference_length=reference_length) if reference_length else None,
        "gc_percent": round((gc_bases / total_bp) * 100, 4) if total_bp else 0.0,
        "reference_length": reference_length,
        "coverage_percent": round((total_bp / denominator) * 100, 4) if denominator else None,
    }
