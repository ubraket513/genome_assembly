# Genome Assembly Roadmap

## Target

Build a production-grade genome assembly library and CLI that starts with
short-read de Bruijn assembly and grows into a high-performance CPU/HPC/GPU
toolkit.

The intended product shape is:

- Python API and friendly CLI for users, notebooks, reports, and orchestration.
- Rust native core for k-mer counting, graph construction, compaction, and
  multicore CPU execution.
- Optional C++/CUDA backend for GPU-heavy kernels once CPU baselines are proven.
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
- ADR selecting Rust as the first native core.
- Rust/PyO3 native scaffold for k-mer counting and edge-table construction.
- Optional native backend bridge with Python parity tests.

Known boundaries:

- The current assembler is correctness-first and not yet optimized for large
  genomes.
- The native backend is scaffolded incrementally and should remain optional
  until parity and packaging are mature.
- Graph cleaning, paired-end logic, minimizer partitioning, and GPU kernels are
  not yet implemented.

## Milestones

### M0: Python Product Baseline

Status: complete.

Keep the Python implementation as the correctness oracle. Every native kernel
must match this behavior before replacing it.

### M1: Native Backend Scaffold

Status: complete.

Add a standalone Rust/PyO3 crate exposing deterministic k-mer counting. Wire the
Python `backend="native"` setting to use native k-mer counts when the extension
is installed, while preserving the default pure-Python install.

Success criteria:

- `cargo test` passes for the Rust crate.
- Python tests pass with no native extension installed.
- A missing native extension produces a clear actionable error.
- API behavior stays stable for `backend="python"`.

### M2: Native Graph Construction

Status: complete for edge-table construction; compaction remains Python.

Move graph edge creation into Rust and return a compact, deterministic edge
table to Python.

Success criteria:

- Native and Python graph summaries match on synthetic fixtures.
- Native edge counting is measurably faster on SARS-CoV-2 and bacterial-scale
  synthetic reads.
- The Python graph API remains unchanged.

### M3: Native Contig Compaction

Status: next.

Move maximal non-branching path compaction into Rust.

Success criteria:

- Contig sequences match Python behavior on linear, branching, circular, tip,
  and bubble fixtures.
- Native assembly produces the same FASTA/GFA outputs modulo ordering rules
  documented in tests.

### M4: Benchmark Harness

Add reproducible benchmark commands and fixture tiers.

Benchmark tiers:

- Tiny synthetic graphs for correctness.
- SARS-CoV-2 fixture for smoke performance.
- Bacterial-scale simulated reads.
- Optional public short-read datasets for larger validation.

Metrics:

- Wall time.
- Peak memory.
- Threads used.
- Contig count, total bp, N50, NG50, GC percent.

### M5: Multicore CPU/HPC

Introduce parallel Rust execution with chunked I/O, minimizer partitioning, and
thread-aware memory controls.

Success criteria:

- `--threads` affects native work.
- Deterministic output across thread counts.
- Benchmark report shows speedup and memory behavior.

### M6: Assembly Quality Features

Add graph cleaning and better short-read assembly behavior.

Feature order:

- Low-abundance edge filtering improvements.
- Tip clipping.
- Simple bubble detection.
- Multi-k assembly exploration.
- Paired-end design hooks.

### M7: GPU Exploration

Prototype GPU only where it is likely to pay off.

Candidate kernels:

- K-mer sorting/counting.
- Minimizer bucketing.
- Read mapping/seed extension.
- Consensus or partial-order alignment.

Avoid moving irregular graph traversal to GPU until profiling proves it is worth
the complexity.

### M8: Distribution and Integration

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
