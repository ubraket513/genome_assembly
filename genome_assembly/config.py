"""Configuration objects for assembly workflows."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AssemblyConfig:
    """Settings for short-read de Bruijn graph assembly."""

    k: int = 31
    min_abundance: int = 1
    min_contig_length: int = 0
    skip_ambiguous: bool = True
    backend: str = "python"
    threads: int = 1

    def validate(self) -> None:
        if self.k < 1:
            raise ValueError("k must be >= 1")
        if self.min_abundance < 1:
            raise ValueError("min_abundance must be >= 1")
        if self.min_contig_length < 0:
            raise ValueError("min_contig_length must be >= 0")
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.backend not in {"python", "cython", "native"}:
            raise ValueError("backend must be 'python', 'cython', or 'native'")
