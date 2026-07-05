# Genome Assembly Implementation Plan

This document is written as a handoff for Claude or another coding agent. Follow
the stages in order. Do not skip validation commands.

## Ground Rules

- Preserve the pure-Python backend as the correctness oracle.
- Keep `backend="python"` as the default.
- Only make `backend="native"` use compiled code when the extension is actually
  installed.
- Add tests before relying on a native replacement.
- Keep generated artifacts out of git unless they are intentional lockfiles.
- Prefer deterministic ordering for all outputs.

## Stage 1: Native Backend Scaffold

Status: complete.

Goal: create the first Rust/PyO3 native module and connect it to Python without
breaking default installs.

Implementation tasks:

- Add `native/genome_assembly_native/Cargo.toml`.
- Add `native/genome_assembly_native/src/lib.rs`.
- Expose `count_kmers(reads, k, skip_ambiguous=True)`.
- Add Rust unit tests for linear reads, repeated k-mers, ambiguous bases, and
  invalid `k`.
- Add `genome_assembly/native.py` with:
  - `native_available()`
  - `require_native()`
  - `count_kmers()`
- Modify graph construction so `AssemblyConfig(backend="native")` uses native
  edge-table construction when available.
- Keep a clear error when the native module is missing.

Validation commands:

```bash
python -m unittest discover -s tests -v
python -m py_compile simul.py de_bruijn.py genome_assembly/*.py tests/test_core.py
cargo test --manifest-path native/genome_assembly_native/Cargo.toml
```

Acceptance criteria:

- Python tests pass without the native extension installed.
- Rust tests pass.
- `backend="python"` behavior is unchanged.
- `backend="native"` fails with a helpful install/build message if the extension
  has not been installed.

## Stage 2: Optional Native Install Flow

Status: complete.

Goal: document and validate a local developer flow for installing the native
extension.

Implementation tasks:

- Add native install instructions to `README.md`.
- Add a short `native/genome_assembly_native/README.md`.
- Validate one local flow:

```bash
python -m pip install maturin
source .venv/bin/activate
unset CONDA_PREFIX  # needed when this shell was launched from Conda
cd native/genome_assembly_native
maturin develop --release
cd ../..
python - <<'PY'
from genome_assembly.native import native_available, count_kmers
print(native_available())
print(count_kmers(["ACGTACGT"], 3))
PY
```

Acceptance criteria:

- Native install docs match the actual command path.
- Native count output matches Python `Counter(iter_kmers(...))` ordering.

## Stage 3: Native Edge Table Integration

Status: complete.

Goal: move graph edge counting and edge table creation fully into Rust.

Implementation tasks:

- Extend Rust output to include `(prefix, suffix, sequence, count)`.
- Add Python parity tests comparing native and Python `GraphSummary`.
- Keep deterministic sort by edge sequence.
- Benchmark edge creation on the bundled SARS-CoV-2 fixture.

Acceptance criteria:

- Native and Python summaries match exactly.
- Native backend produces identical contig sequences on small fixtures.

## Stage 4: Native Contig Compaction

Status: next.

Goal: move maximal non-branching path compaction into Rust.

Implementation tasks:

- Implement Rust adjacency maps and one-in-one-out traversal.
- Add circular graph handling.
- Add parity tests for:
  - linear genome
  - overlapping reads
  - branch
  - cycle
  - tip
  - bubble
- Return contigs with name, sequence, mean abundance, and edge count.

Acceptance criteria:

- Native contigs match Python contigs on all fixtures.
- Python objects and CLI output remain unchanged.

## Stage 5: Benchmark Harness

Goal: make performance claims reproducible.

Implementation tasks:

- Add `ga benchmark`.
- Add benchmark fixture generation.
- Record JSON output with command, config, wall time, contig stats, and platform
  metadata.
- Add docs explaining benchmark tiers and limitations.

Acceptance criteria:

- Benchmarks can compare `backend="python"` and `backend="native"`.
- Results are written in machine-readable JSON.

## Stage 6: Multicore CPU/HPC

Goal: make `--threads` meaningful for native work.

Implementation tasks:

- Add Rayon to the Rust crate.
- Parallelize read chunk processing.
- Merge per-thread k-mer maps deterministically.
- Add tests proving identical output for `threads=1` and `threads>1`.
- Add memory budget configuration before processing large files.

Acceptance criteria:

- Thread count changes runtime, not output.
- Benchmarks include thread count and speedup.

## Stage 7: Graph Cleaning

Goal: improve short-read assembly quality.

Implementation tasks:

- Implement tip clipping.
- Implement simple bubble detection.
- Improve abundance filtering with reportable removed-edge counts.
- Add CLI flags for cleaning behavior, defaulting conservatively.

Acceptance criteria:

- Tests show error k-mers are removed without breaking clean linear assemblies.
- Summary reports graph-cleaning counts.

## Stage 8: GPU Research Spike

Goal: decide whether GPU acceleration is worth implementing.

Implementation tasks:

- Profile native CPU backend first.
- Identify the top two bottlenecks.
- Prototype only one GPU-friendly kernel, likely k-mer sorting/counting or
  minimizer bucketing.
- Compare against CPU native backend on wall time, transfer overhead, memory, and
  reproducibility.

Acceptance criteria:

- Written decision: continue CUDA backend, defer it, or drop it.
- No production GPU dependency unless the benchmark justifies it.

## Stage 9: Release Hardening

Goal: make the project usable by others.

Implementation tasks:

- Add CI for Python tests, Rust tests, and linting.
- Add wheel build workflow.
- Add examples and troubleshooting docs.
- Add changelog and versioning policy.
- Add license file.

Acceptance criteria:

- Fresh clone can run all documented commands.
- Release artifacts are reproducible.
- Docs state which backends are stable and which are experimental.
