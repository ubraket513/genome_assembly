"""Minimizer and syncmer sketching primitives (pure-Python oracle).

These are the building blocks for minimizer-space assembly (mdBG) and for
HPC-style minimizer bucketing, where each k-mer is deterministically routed to
an independent stream by the hash of its minimizer.

This module is the correctness oracle. The Rust native crate mirrors it for
speed; `test_core.py` asserts parity when the extension is built.
"""

from __future__ import annotations

from .kmers import normalize_sequence

_FNV64_OFFSET = 0xCBF29CE484222325
_FNV64_PRIME = 0x100000001B3
_MASK64 = 0xFFFFFFFFFFFFFFFF


def hash64(text: str) -> int:
    """Deterministic 64-bit FNV-1a hash of an ASCII string.

    Python's built-in ``hash`` is per-process randomized; sketching needs a
    stable order, so we fix the hash here and in the Rust mirror.
    """

    h = _FNV64_OFFSET
    for byte in text.encode("ascii", "ignore"):
        h = ((h ^ byte) * _FNV64_PRIME) & _MASK64
    return h


def minimizers(sequence: str, w: int, m: int) -> list[tuple[int, str]]:
    """Return windowed (w, m)-minimizers as ``(position, m-mer)`` pairs.

    For each window of ``w`` consecutive m-mers, the lowest-hash m-mer is
    selected (leftmost on ties). Consecutive windows that pick the same position
    are de-duplicated, so the result is the ordered set of distinct minimizers.
    """

    if w < 1:
        raise ValueError("w must be >= 1")
    if m < 1:
        raise ValueError("m must be >= 1")

    sequence = normalize_sequence(sequence)
    mmer_count = len(sequence) - m + 1
    if mmer_count < 1:
        return []

    hashes = [hash64(sequence[i : i + m]) for i in range(mmer_count)]

    selected: list[tuple[int, str]] = []
    last_position = -1
    window_count = mmer_count - w + 1
    if window_count < 1:
        # Sequence shorter than one full window: take the global minimum.
        window_count = 1
        w = mmer_count

    for start in range(window_count):
        best = start
        for offset in range(start + 1, start + w):
            if hashes[offset] < hashes[best]:
                best = offset
        if best != last_position:
            selected.append((best, sequence[best : best + m]))
            last_position = best
    return selected


def syncmers(sequence: str, k: int, s: int, t: int = 0) -> list[tuple[int, str]]:
    """Return open ``(k, s)``-syncmers as ``(position, k-mer)`` pairs.

    A k-mer is an open syncmer when the lowest-hash s-mer inside it starts at
    offset ``t``. Open syncmers have a bounded, window-guaranteed spacing that
    minimizers lack, which is why they are preferred for modern sketching.
    """

    if k < 1:
        raise ValueError("k must be >= 1")
    if s < 1 or s > k:
        raise ValueError("s must satisfy 1 <= s <= k")
    inner = k - s + 1
    if t < 0 or t >= inner:
        raise ValueError("t must satisfy 0 <= t < k - s + 1")

    sequence = normalize_sequence(sequence)
    if len(sequence) < k:
        return []

    result: list[tuple[int, str]] = []
    for start in range(len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        best_offset = 0
        best_hash = hash64(kmer[0:s])
        for offset in range(1, inner):
            candidate = hash64(kmer[offset : offset + s])
            if candidate < best_hash:
                best_hash = candidate
                best_offset = offset
        if best_offset == t:
            result.append((start, kmer))
    return result


def minimizer_bucket(kmer: str, m: int, num_buckets: int) -> int:
    """Route a k-mer to one of ``num_buckets`` streams by its minimizer hash.

    Identical k-mers always route to the same bucket, and consecutive k-mers
    that share a minimizer land together, so buckets are independent, lock-free
    partitions for parallel or distributed counting.
    """

    if num_buckets < 1:
        raise ValueError("num_buckets must be >= 1")
    kmer = normalize_sequence(kmer)
    if m < 1 or m > len(kmer):
        raise ValueError("m must satisfy 1 <= m <= len(kmer)")
    best = min(hash64(kmer[i : i + m]) for i in range(len(kmer) - m + 1))
    return best % num_buckets
