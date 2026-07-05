# ADR-003: Use Rust as the Primary Native Core

## Status

Accepted

## Context

The project now has three viable performance paths:

- Pure Python, which remains the correctness oracle and easiest install path.
- Cython, which is close to Python and useful for near-term CPython speedups.
- Rust/PyO3, which already exists in the repository and can grow into a safer
  production native core for graph algorithms, multicore execution, and wheel
  packaging.

C++ is also a serious option for assembly work because many bioinformatics and
HPC libraries are written in C++, and C++ has the most direct CUDA ecosystem.
However, adopting C++ as the default core would increase memory-safety and build
complexity before benchmarks prove that trade-off is necessary.

## Decision

Use Rust as the preferred production native core. Keep Cython as a tactical
CPython accelerator and Python-adjacent prototyping lane. Evaluate C++ through
isolated pybind11 or CMake-backed modules only when one of these is true:

- We need to bind an established C++ bioinformatics library.
- A SIMD/CUDA-heavy kernel is clearly easier or faster in C++.
- Benchmarks show Rust/Cython cannot meet a specific performance target.

## Considered Options

### Rust/PyO3 Primary

Pros:

- Strong memory safety for graph and sequence-processing code.
- Good Python extension story through PyO3 and maturin.
- Good path to deterministic multicore CPU work through Rayon-style designs.
- Already scaffolded in this repository.

Cons:

- Smaller existing assembler ecosystem than C++.
- CUDA integration is less direct than C++.
- Requires Rust toolchain for contributors working on native kernels.

### Cython Primary

Pros:

- Fastest path from Python prototype to compiled CPython extension.
- Keeps algorithm code close to the Python oracle.
- Good for incremental hot-path acceleration.

Cons:

- CPython-oriented, not PyPy-friendly.
- Less attractive for long-term HPC kernels, memory ownership, and threading.
- Can drift into Python-object-heavy code if not carefully profiled.

### C++/pybind11 Primary

Pros:

- Mature HPC and CUDA ecosystem.
- Direct access to many existing bioinformatics libraries.
- Strong control over low-level memory layout and SIMD.

Cons:

- Higher memory-safety risk.
- More complex cross-platform builds and ABI management.
- Easier to overcommit before real bottlenecks are measured.

## Rationale

Rust is the best default for new native code because it balances speed,
maintainability, memory safety, and Python packaging. Cython remains valuable as
a short feedback-loop accelerator. C++ should be pulled in deliberately, not as
the default, when existing libraries or GPU kernels justify the cost.

## Consequences

- `backend="native"` means Rust/PyO3 unless a future backend name says
  otherwise.
- `backend="cython"` remains supported for CPython acceleration and parity
  comparisons.
- C++ work should enter as a documented spike with benchmark targets, not as an
  unbounded rewrite.
- The benchmark harness must compare Python, Cython, and Rust first; C++ enters
  only after those baselines exist.

## Implementation Notes

- Finish Rust contig compaction parity before adding multicore execution.
- Keep pure-Python behavior as oracle tests for every native kernel.
- Use maturin for Rust extension development and wheel builds.
- Use pybind11/CMake only for a future `cpp` backend or specialist extension.
