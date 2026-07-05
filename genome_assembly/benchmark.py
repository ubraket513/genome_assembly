"""Reproducible benchmark helpers for assembly backends."""

from __future__ import annotations

from dataclasses import asdict
import gc
from pathlib import Path
import platform
import sys
import threading
import time
import tracemalloc
from typing import Iterable

from .assemble import assemble_short_reads
from .config import AssemblyConfig
from .cython_backend import cython_available
from .io import read_fasta
from .native import native_available
from .simulate import generate_random_genome, simulate_reads

DEFAULT_BACKENDS = ("python", "cython", "native")


def _current_rss_bytes() -> int | None:
    """Resident set size of this process in bytes, or None if unavailable.

    Reads Linux /proc, which captures native (Rust) allocations that Python's
    tracemalloc cannot see. Returns None on platforms without /proc.
    """

    try:
        with open("/proc/self/statm", "r", encoding="ascii") as handle:
            resident_pages = int(handle.read().split()[1])
    except (OSError, IndexError, ValueError):
        return None
    import resource

    return resident_pages * resource.getpagesize()


class _PeakRSSSampler:
    """Poll process RSS on a background thread and keep the maximum seen.

    Peak RSS is the honest memory metric for comparing native backends, since
    tracemalloc only tracks the Python heap.
    """

    def __init__(self, interval: float = 0.01) -> None:
        self.interval = interval
        self.baseline = _current_rss_bytes()
        self.peak = self.baseline
        self._supported = self.baseline is not None
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def __enter__(self) -> "_PeakRSSSampler":
        if self._supported:
            self._thread = threading.Thread(target=self._run, daemon=True)
            self._thread.start()
        return self

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            current = _current_rss_bytes()
            if current is not None and (self.peak is None or current > self.peak):
                self.peak = current

    def __exit__(self, *exc: object) -> None:
        if self._thread is not None:
            self._stop.set()
            self._thread.join()
        current = _current_rss_bytes()  # final reading in case the run was short
        if current is not None and (self.peak is None or current > self.peak):
            self.peak = current

    def result(self) -> dict[str, int | None]:
        delta = None
        if self.peak is not None and self.baseline is not None:
            delta = max(self.peak - self.baseline, 0)
        return {
            "peak_rss_bytes": self.peak,
            "baseline_rss_bytes": self.baseline,
            "rss_delta_bytes": delta,
        }


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
    reference: str | Path | None = None,
    *,
    synthetic_genome_size: int | None = None,
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
    """Benchmark selected backends on deterministic simulated reads.

    Provide either ``reference`` (a FASTA path) or ``synthetic_genome_size`` to
    generate a deterministic random genome of that many bases in memory. The
    synthetic path enables bacterial-scale and larger tiers without downloads.
    """

    selected_backends = parse_backend_list(backends)
    if (reference is None) == (synthetic_genome_size is None):
        raise ValueError("provide exactly one of reference or synthetic_genome_size")

    if synthetic_genome_size is not None:
        reference_name = f"synthetic_{synthetic_genome_size}bp"
        reference_label = reference_name
        reference_sequence = generate_random_genome(synthetic_genome_size, seed=seed)
    else:
        reference_record = read_fasta(reference)[0]
        reference_name = reference_record.name
        reference_label = str(reference)
        reference_sequence = reference_record.sequence

    simulated = simulate_reads(
        reference_sequence,
        read_length=read_length,
        coverage=coverage,
        seed=seed,
    )
    reference_length = len(reference_sequence)

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

        gc.collect()  # release the prior backend's objects so RSS delta is per-run
        tracemalloc.start()
        sampler = _PeakRSSSampler()
        start = time.perf_counter()
        try:
            with sampler:
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
                    **sampler.result(),
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
                    **sampler.result(),
                    "graph": asdict(result.summary),
                    "stats": result.stats(reference_length=reference_length),
                }
            )
        finally:
            tracemalloc.stop()

    mean_coverage = sum(simulated.coverage) / len(simulated.coverage)
    return {
        "benchmark": {
            "command": command,
            "reference": reference_label,
            "reference_name": reference_name,
            "reference_length": reference_length,
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
