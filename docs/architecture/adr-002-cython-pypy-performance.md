# ADR-002: Prioritize Cython and PyPy-Compatible Performance Paths

## Status

Superseded by ADR-003

## Context

The project needs performance improvements while staying approachable as a
Python library and CLI. The previous direction made Rust the first native core.
The current priority is to optimize within the Python ecosystem before pushing
more work into external native stacks.

## Decision

Use Cython as the first compiled accelerator for CPython and keep the pure-Python
backend compatible with PyPy where practical. Rust remains an experimental
secondary backend for future multicore/HPC kernels.

## Rationale

- Cython keeps implementation close to the Python correctness oracle.
- Cython can accelerate graph construction and contig compaction without changing
  the public API.
- PyPy may speed up pure-Python loops through JIT compilation, but PyPy is not a
  good target for CPython C extension speedups.
- Keeping `backend="python"` as the default gives a clean PyPy-compatible path.

## Trade-offs

- Cython extensions target CPython, not PyPy portability.
- Cython is easier to adopt now but may not match Rust/C++ for long-term HPC and
  CUDA-heavy kernels.
- The project must maintain parity tests across Python, Cython, and optional
  Rust paths.

## Consequences

- `backend="cython"` becomes the preferred near-term accelerator.
- `backend="python"` remains the portable baseline and PyPy candidate.
- `backend="native"` remains available but experimental.
- Benchmarking should compare CPython/Python, CPython/Cython, PyPy/Python, and
  Rust/native when available.

## Supersession Note

ADR-003 keeps the Cython backend but changes the strategic priority: Rust is the
preferred production native core, Cython is a near-term CPython accelerator, and
C++ is reserved for benchmark-justified specialist kernels or integrations.
