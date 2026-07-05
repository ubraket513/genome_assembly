# ADR-001: Use Rust as the First Native Assembly Core

## Status

Superseded by ADR-003

## Context

The current project began as Python prototype code for read simulation and
de Bruijn graph assembly. Production use needs faster k-mer processing, lower
memory overhead, multicore execution, and a clean path toward HPC/GPU backends.

## Decision

Use a Python public API and CLI with a Rust native backend as the first
acceleration target. Keep C++/CUDA as an optional backend for GPU-heavy kernels.

This was first superseded by ADR-002, which prioritized Cython and
PyPy-compatible Python paths. ADR-003 restores Rust as the preferred production
native core while retaining Cython as a tactical accelerator.

## Rationale

- Rust gives native performance with stronger memory safety than C/C++.
- PyO3/maturin provides a practical Python packaging path.
- Rayon-style parallelism is a good fit for k-mer counting, minimizer bucketing,
  and graph compaction.
- C++/CUDA remains available for specialized GPU work without forcing all core
  development into a higher-risk codebase from day one.

## Trade-offs

- Rust has fewer mature bioinformatics assembly libraries than C++.
- CUDA integration is less direct than in C++.
- The first pure-Python backend must remain tested so behavior stays clear while
  the native backend is introduced.
