"""K-mer utilities."""

from __future__ import annotations

DNA_ALPHABET = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def normalize_sequence(sequence: str) -> str:
    """Return an uppercase sequence with whitespace removed."""

    return "".join(str(sequence).split()).upper()


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""

    return normalize_sequence(sequence).translate(_COMPLEMENT)[::-1]


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically canonical orientation of a k-mer."""

    kmer = normalize_sequence(kmer)
    return min(kmer, reverse_complement(kmer))


def is_unambiguous_dna(sequence: str) -> bool:
    """Return True when a sequence contains only A, C, G, and T."""

    return set(normalize_sequence(sequence)).issubset(DNA_ALPHABET)


def iter_kmers(sequence: str, k: int, *, skip_ambiguous: bool = True):
    """Yield k-mers from a DNA sequence."""

    if k < 1:
        raise ValueError("k must be >= 1")

    sequence = normalize_sequence(sequence)
    if len(sequence) < k:
        return

    for start in range(0, len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        if skip_ambiguous and not is_unambiguous_dna(kmer):
            continue
        yield kmer
