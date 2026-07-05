# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While the version is `0.x`, minor releases may include breaking API changes and
patch releases are reserved for backward-compatible fixes.

## [Unreleased]

### Added

- Parallel native k-mer counting on a Rayon thread pool sized by `--threads`,
  with deterministic per-thread map merging. Thread count changes runtime, not
  output.
- Graph cleaning: conservative tip clipping (`--tip-length`) and simple bubble
  popping (`--bubble-length`), both disabled by default. Removed-edge counts are
  reported in `GraphSummary` and the `summary.json` graph block.
- `GraphSummary.tips_removed` and `GraphSummary.bubble_edges_removed` fields.
- Continuous integration for Python tests, Rust tests, clippy, and formatting.
- `LICENSE` (MIT) and this changelog.

### Changed

- `--threads` now drives the native backend instead of being a reserved no-op.

## [0.1.0]

### Added

- Installable Python package with a Typer CLI (`ga`): `simulate`, `assemble`,
  `stats`, `benchmark`.
- Dependency-light FASTA/FASTQ I/O and deterministic read simulation.
- Pure-Python de Bruijn graph construction and maximal non-branching path
  compaction.
- Assembly metrics including N50/N90 and reference-aware NG50/NG90.
- Cython accelerator backend (`backend="cython"`).
- Rust/PyO3 native backend (`backend="native"`) for k-mer counting, edge-table
  construction, and contig compaction, with Python parity tests.
