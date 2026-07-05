# Genome Assembly Roadmap

## Target

Build a production-grade genome assembly library and CLI that starts with
short-read de Bruijn assembly and grows into a high-performance CPU/HPC/GPU
toolkit.

The intended product shape is:

- Python API and friendly CLI for users, notebooks, reports, and orchestration.
- Rust native core for production hot paths, multicore CPU work, and future HPC
  kernels.
- Cython accelerator for tactical CPython hot paths and rapid native prototypes.
- PyPy-compatible pure-Python backend where practical, without treating PyPy as
  a compiled-backend target.
- Optional C++/CUDA lane for existing C++ libraries, SIMD-heavy code, or
  GPU-heavy kernels once Python, Cython, and Rust baselines are proven.
- Interoperability with FASTA, FASTQ, GFA, QUAST, and established assemblers.

## Current State

Completed:

- Installable Python package via `pyproject.toml`.
- Dependency-light FASTA/FASTQ I/O.
- Deterministic read simulation.
- Pure-Python de Bruijn graph construction from `(k + 1)`-mers.
- Maximal non-branching path compaction into contigs.
- Assembly metrics including `N50/N90` and reference-aware `NG50/NG90`.
- Typer CLI: `simulate`, `assemble`, and `stats`.
- Initial tests and CLI smoke flow.
- Rust/PyO3 native backend for k-mer counting, edge-table construction, and
  contig compaction.
- Optional native backend bridge with Python parity tests.
- Historical ADR-002 for the Cython/PyPy lane.
- ADR-003 selecting Rust as the primary production native core while keeping
  Cython and C++ as scoped lanes.
- Cython backend for k-mer counting, graph edge-table construction, and contig
  compaction.
- CLI `--backend` selection for `python`, `cython`, and `native`.
- CLI `benchmark` command for reproducible backend comparisons on simulated
  reads.
- Parallel native k-mer counting honoring `--threads`, with deterministic
  per-thread map merging and Rust/Python parity tests.
- Graph cleaning: conservative tip clipping and simple bubble popping with
  reportable removed-edge counts, shared across all backends.

Known boundaries:

- The current assembler is correctness-first and not yet optimized for large
  genomes.
- Rust is the preferred production native core; Cython remains a tactical
  CPython accelerator.
- C++ is not a default core. Use it only for benchmark-justified specialist
  kernels or external C++/CUDA integrations.
- PyPy has not been validated locally yet and should run only the pure-Python
  backend.
- Paired-end logic, minimizer partitioning, multi-k assembly, and GPU kernels
  are not yet implemented. Graph cleaning covers conservative tip clipping and
  simple bubble popping only.

## Milestones

### M0: Python Product Baseline

Status: complete.

Keep the Python implementation as the correctness oracle. Every accelerator must
match this behavior before replacing it.

Success criteria:

- Python tests pass with no compiled extensions installed.
- API behavior stays stable for `backend="python"`.

### M1: Cython CPython Accelerator

Status: complete.

Use Cython for k-mer counting, graph edge-table construction, and contig
compaction while preserving the public Python API.

Success criteria:

- Python tests pass without the Cython extension installed.
- Cython parity tests pass when the extension is built.
- CLI supports `--backend cython`.
- Missing Cython extension produces a clear actionable error.

### M2: Rust Native Core

Status: complete for k-mer counting, edge-table construction, and contig
compaction.

Move graph hot paths into Rust while preserving Python API behavior.

Success criteria:

- `cargo test` passes for the Rust crate.
- Native and Python graph summaries match on synthetic fixtures.
- Native and Python contig sequences match on small fixtures.
- Native backend is ready for SARS-CoV-2 and bacterial-scale benchmark
  comparison.
- The Python graph API remains unchanged.

### M3: Benchmark Harness

Status: complete for deterministic simulated-read benchmarks.

Add reproducible benchmark commands and fixture tiers.

Benchmark tiers:

- Tiny synthetic graphs for correctness.
- SARS-CoV-2 fixture for smoke performance.
- Bacterial-scale simulated reads.
- Optional public short-read datasets for larger validation.

Metrics:

- Wall time.
- Peak traced Python memory initially; process RSS should be added for native
  memory accounting later.
- Threads used.
- Contig count, total bp, N50, NG50, GC percent.

### M4: Multicore CPU/HPC

Status: complete for parallel native k-mer counting; chunked I/O and minimizer
partitioning still open.

Parallel native k-mer counting now runs on a Rayon thread pool sized by
`--threads`, with deterministic per-thread map merging. Chunked/streaming I/O,
minimizer partitioning, and thread-aware memory controls remain future work and
are deferred until a streaming read path exists. Cython stays single-process.

Success criteria:

- `--threads` affects native work. Done.
- Deterministic output across thread counts. Done (Rust + Python parity tests).
- Benchmark report shows speedup and memory behavior. Partial: `ga benchmark`
  honors `--threads` and records wall time; peak memory is still Python-traced
  only, so native RSS accounting stays open (see M3).

### M5: PyPy Compatibility Lane

Status: blocked on toolchain (`pypy3` not installed locally). The pure-Python
backend has no CPython-only dependency, so it is expected to run on PyPy; this
is unvalidated until a PyPy interpreter is available.

Validate the pure-Python backend under PyPy and document unsupported paths.

Success criteria:

- PyPy can run core tests for `backend="python"`.
- Docs state that Cython/Rust extension backends are CPython-oriented.
- Any PyPy-specific test exclusions are explicit.

### M6: C++/CUDA Evaluation Lane

Status: deferred (ADR-004). No benchmark bottleneck justifies C++ yet; re-entry
criteria are documented.

Evaluate C++ only when it solves a specific problem better than Rust or Cython.

Entry criteria:

- Existing C++ bioinformatics library needs binding.
- CUDA/SIMD kernel is clearly better served by C++ tooling.
- Benchmark harness shows a bottleneck Rust/Cython cannot meet.

Success criteria:

- A pybind11 or CMake-backed spike has isolated scope and parity tests.
- The decision to keep or drop C++ is documented with benchmark data.

### M7: Assembly Quality Features

Status: initial graph cleaning complete; multi-k and paired-end still open.

Add graph cleaning and better short-read assembly behavior.

Feature order:

- Done: reportable removed-edge counts (`tips_removed`, `bubble_edges_removed`).
- Done: tip clipping (`--tip-length`, conservative, off by default).
- Done: simple bubble popping (`--bubble-length`, keeps highest-coverage path).
- Open: multi-k assembly exploration.
- Open: paired-end design hooks.

### M8: GPU Exploration

Status: deferred (ADR-004). Revisit once profiling identifies a regular,
GPU-friendly kernel that beats the CPU native backend end to end.

Prototype GPU only where it is likely to pay off.

Candidate kernels:

- K-mer sorting/counting.
- Minimizer bucketing.
- Read mapping/seed extension.
- Consensus or partial-order alignment.

Avoid moving irregular graph traversal to GPU until profiling proves it is worth
the complexity.

### M9: Distribution and Integration

Status: complete for the initial release.

Prepare release-quality packaging.

Targets:

- Done: pure-Python wheel builds with `python -m build` and installs without a
  compiler.
- Done: native wheels build via maturin for Linux/macOS/Windows in the release
  workflow.
- Partial: CLI has useful, actionable errors and reproducible reports. Typer
  provides shell completion out of the box; a documented completion recipe is
  future polish.
- Done: documentation includes install, examples, troubleshooting, backend
  status, and the staged plan.

## Strategic Edge

Do not try to beat SPAdes, MEGAHIT, hifiasm, Flye, or Verkko on every dimension
immediately. The best edge is a production-grade library and CLI that is easier
to embed, easier to inspect, benchmark-friendly, and eventually accelerated by a
safe native core.
