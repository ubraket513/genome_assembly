"""Reproducible benchmark helpers for assembly backends."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import platform
import sys
import time
import tracemalloc
from typing import Iterable

from .assemble import assemble_short_reads
from .config import AssemblyConfig
from .cython_backend import cython_available
from .io import read_fasta
from .native import native_available
from .simulate import simulate_reads

DEFAULT_BACKENDS = ("python", "cython", "native")


def parse_backend_list(value: str | Iterable[str]) -> list[str]:
    """Parse backend selections from comma-separated text or an iterable."""

    if isinstance(value, str):
        backends = [part.strip() for part in value.split(",")]
    else:
        backends = [str(part).strip() for part in value]

    parsed = [backend for backend in backends if backend]
    if not parsed:
        raise ValueError("at least one backend is required")
    return parsed


def _backend_unavailable_reason(backend: str) -> str | None:
    if backend == "python":
        return None
    if backend == "cython" and not cython_available():
        return "Cython extension is not built"
    if backend == "native" and not native_available():
        return "Rust native extension is not installed or is outdated"
    return None


def _platform_metadata() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
    }


def run_benchmark(
    reference: str | Path,
    *,
    backends: Iterable[str] = DEFAULT_BACKENDS,
    k: int = 31,
    read_length: int = 150,
    coverage: float = 5.0,
    seed: int | None = 7,
    min_abundance: int = 1,
    min_contig_length: int = 0,
    threads: int = 1,
    command: list[str] | None = None,
) -> dict[str, object]:
    """Benchmark selected backends on deterministic simulated reads."""

    selected_backends = parse_backend_list(backends)
    records = read_fasta(reference)
    reference_record = records[0]
    simulated = simulate_reads(
        reference_record.sequence,
        read_length=read_length,
        coverage=coverage,
        seed=seed,
    )

    runs: list[dict[str, object]] = []
    for backend in selected_backends:
        config = AssemblyConfig(
            k=k,
            min_abundance=min_abundance,
            min_contig_length=min_contig_length,
            backend=backend,
            threads=threads,
        )
        try:
            config.validate()
        except ValueError as exc:
            runs.append({"backend": backend, "status": "error", "error": str(exc)})
            continue

        unavailable_reason = _backend_unavailable_reason(backend)
        if unavailable_reason is not None:
            runs.append({"backend": backend, "status": "skipped", "reason": unavailable_reason})
            continue

        tracemalloc.start()
        start = time.perf_counter()
        try:
            result = assemble_short_reads(simulated.reads, config)
            _, peak_bytes = tracemalloc.get_traced_memory()
            wall_time_seconds = time.perf_counter() - start
        except Exception as exc:  # noqa: BLE001 - benchmark report should capture backend failures.
            _, peak_bytes = tracemalloc.get_traced_memory()
            wall_time_seconds = time.perf_counter() - start
            runs.append(
                {
                    "backend": backend,
                    "status": "error",
                    "wall_time_seconds": wall_time_seconds,
                    "peak_traced_bytes": peak_bytes,
                    "error": str(exc),
                }
            )
        else:
            runs.append(
                {
                    "backend": backend,
                    "status": "ok",
                    "wall_time_seconds": wall_time_seconds,
                    "peak_traced_bytes": peak_bytes,
                    "graph": asdict(result.summary),
                    "stats": result.stats(reference_length=len(reference_record.sequence)),
                }
            )
        finally:
            tracemalloc.stop()

    mean_coverage = sum(simulated.coverage) / len(simulated.coverage)
    return {
        "benchmark": {
            "command": command,
            "reference": str(reference),
            "reference_name": reference_record.name,
            "reference_length": len(reference_record.sequence),
            "read_count": len(simulated.reads),
            "read_length": simulated.read_length,
            "target_coverage": simulated.target_coverage,
            "mean_coverage": mean_coverage,
            "seed": simulated.seed,
            "backends": selected_backends,
            "config": {
                "k": k,
                "min_abundance": min_abundance,
                "min_contig_length": min_contig_length,
                "threads": threads,
            },
        },
        "platform": _platform_metadata(),
        "runs": runs,
    }
