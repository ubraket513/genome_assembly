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

Known boundaries:

- The current assembler is correctness-first and not yet optimized for large
  genomes.
- Rust is the preferred production native core; Cython remains a tactical
  CPython accelerator.
- C++ is not a default core. Use it only for benchmark-justified specialist
  kernels or external C++/CUDA integrations.
- PyPy has not been validated locally yet and should run only the pure-Python
  backend.
- Graph cleaning, paired-end logic, minimizer partitioning, and GPU kernels are
  not yet implemented.

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

Status: next.

Introduce parallel Rust execution with chunked I/O, minimizer partitioning, and
thread-aware memory controls. Keep Cython single-process unless profiling shows
clear value.

Success criteria:

- `--threads` affects native work.
- Deterministic output across thread counts.
- Benchmark report shows speedup and memory behavior.

### M5: PyPy Compatibility Lane

Status: optional after benchmark harness.

Validate the pure-Python backend under PyPy and document unsupported paths.

Success criteria:

- PyPy can run core tests for `backend="python"`.
- Docs state that Cython/Rust extension backends are CPython-oriented.
- Any PyPy-specific test exclusions are explicit.

### M6: C++/CUDA Evaluation Lane

Status: gated by benchmarks.

Evaluate C++ only when it solves a specific problem better than Rust or Cython.

Entry criteria:

- Existing C++ bioinformatics library needs binding.
- CUDA/SIMD kernel is clearly better served by C++ tooling.
- Benchmark harness shows a bottleneck Rust/Cython cannot meet.

Success criteria:

- A pybind11 or CMake-backed spike has isolated scope and parity tests.
- The decision to keep or drop C++ is documented with benchmark data.

### M7: Assembly Quality Features

Add graph cleaning and better short-read assembly behavior.

Feature order:

- Low-abundance edge filtering improvements.
- Tip clipping.
- Simple bubble detection.
- Multi-k assembly exploration.
- Paired-end design hooks.

### M8: GPU Exploration

Prototype GPU only where it is likely to pay off.

Candidate kernels:

- K-mer sorting/counting.
- Minimizer bucketing.
- Read mapping/seed extension.
- Consensus or partial-order alignment.

Avoid moving irregular graph traversal to GPU until profiling proves it is worth
the complexity.

### M9: Distribution and Integration

Prepare release-quality packaging.

Targets:

- Pure-Python wheel remains easy to install.
- Native wheels are built for common Linux/macOS platforms.
- CLI supports shell completion, useful errors, and reproducible reports.
- Documentation includes install, examples, performance notes, and handoff
  instructions.

## Strategic Edge

Do not try to beat SPAdes, MEGAHIT, hifiasm, Flye, or Verkko on every dimension
immediately. The best edge is a production-grade library and CLI that is easier
to embed, easier to inspect, benchmark-friendly, and eventually accelerated by a
safe native core.
