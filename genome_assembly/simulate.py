"""Read simulation utilities."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .kmers import normalize_sequence


@dataclass(frozen=True)
class SimulatedReads:
    reads: list[str]
    coverage: list[int]
    read_length: int
    target_coverage: float
    seed: int | None


def simulate_reads(
    sequence: str,
    *,
    read_length: int = 150,
    coverage: float = 30.0,
    seed: int | None = None,
) -> SimulatedReads:
    """Simulate fixed-length reads at approximately the requested depth."""

    sequence = normalize_sequence(sequence)
    if not sequence:
        raise ValueError("sequence must not be empty")
    if read_length < 1:
        raise ValueError("read_length must be >= 1")
    if read_length > len(sequence):
        raise ValueError("read_length cannot exceed sequence length")
    if coverage <= 0:
        raise ValueError("coverage must be > 0")

    rng = random.Random(seed)
    read_count = max(1, math.ceil(len(sequence) * coverage / read_length))
    per_base = [0] * len(sequence)
    reads: list[str] = []
    max_start = len(sequence) - read_length

    for _ in range(read_count):
        start = rng.randint(0, max_start)
        end = start + read_length
        reads.append(sequence[start:end])
        for index in range(start, end):
            per_base[index] += 1

    return SimulatedReads(
        reads=reads,
        coverage=per_base,
        read_length=read_length,
        target_coverage=coverage,
        seed=seed,
    )
