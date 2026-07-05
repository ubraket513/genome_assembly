# ADR-004: Defer the C++ and GPU Lanes Until Benchmarks Justify Them

## Status

Accepted

## Context

The roadmap reserves two specialist lanes:

- A C++/pybind11 lane (Stage 9 / M6) for binding existing C++ bioinformatics
  libraries or writing SIMD/CUDA-heavy kernels.
- A GPU research spike (Stage 11 / M8) for kernels such as k-mer sorting,
  minimizer bucketing, or read mapping.

Both were always gated on benchmark evidence. The benchmark harness now exists
(`ga benchmark`) and the native Rust backend already parallelizes k-mer counting
across threads with measured speedup (about 1.7x at 4 threads on 66k synthetic
reads). No profiling run to date has identified a bottleneck that Rust or Cython
cannot address on the CPU.

## Decision

Do not open the C++ or GPU lanes for the initial release. Keep Rust/PyO3 as the
production native core and Cython as the tactical CPython accelerator, per
ADR-003.

Revisit each lane only when its documented entry criteria are met:

- **C++ (Stage 9):** an established C++ library needs binding, or a benchmark
  shows a specific bottleneck that Rust/Cython measurably cannot meet.
- **GPU (Stage 11):** profiling identifies a regular, high-arithmetic-intensity
  kernel (k-mer sort/count or minimizer bucketing are the likeliest), and a
  prototype beats the CPU native backend on wall time after accounting for host
  to device transfer, memory, and reproducibility.

## Rationale

- The current bottlenecks are unmeasured at genome scale, so committing to C++
  build/ABI complexity or a CUDA dependency now would be speculative — exactly
  the overcommitment ADR-003 warns against.
- Irregular graph traversal (compaction, cleaning) is a poor GPU fit and is not
  a measured hot path; moving it to a GPU would add complexity for no proven
  gain.
- Deferring keeps a fresh clone buildable with only a Python toolchain (plus an
  optional Rust toolchain), which matters for the first release.

## Consequences

- The initial release ships Python, Cython, and Rust native backends only.
- The C++ and GPU sections of the roadmap remain open with explicit,
  benchmark-driven entry criteria rather than a scheduled implementation.
- When a lane is opened, it enters as an isolated spike with parity tests and a
  follow-up ADR recording the benchmark evidence that justified it.
