# Genome Assembly Implementation Plan

This document is written as a handoff for Claude or another coding agent. Follow
the stages in order. Do not skip validation commands.

## Ground Rules

- Preserve the pure-Python backend as the correctness oracle.
- Keep `backend="python"` as the default.
- Prefer `backend="native"` for production native work. This means Rust/PyO3.
- Keep `backend="cython"` as a tactical CPython accelerator.
- Consider C++ only for isolated, benchmark-justified pybind11/CMake spikes.
- Use PyPy only with `backend="python"` unless a backend is explicitly validated
  under PyPy.
- Only make compiled backends use compiled code when the extension is actually
  installed.
- Add tests before relying on a native replacement.
- Keep generated artifacts out of git unless they are intentional lockfiles.
- Prefer deterministic ordering for all outputs.

## Stage 1: Python Product Baseline

Status: complete.

Goal: keep the pure-Python package installable, tested, and suitable as the
correctness oracle.

Implementation tasks:

- Maintain FASTA/FASTQ I/O, simulation, graph construction, contig compaction,
  CLI, and tests.
- Keep `backend="python"` dependency-light.

Validation commands:

```bash
python -m unittest discover -s tests -v
python -m py_compile simul.py de_bruijn.py genome_assembly/*.py tests/test_core.py
```

Acceptance criteria:

- Python tests pass without compiled extensions installed.
- `backend="python"` behavior is unchanged.

## Stage 2: Cython Backend

Status: complete.

Goal: keep Cython available as the tactical CPython acceleration path.

Implementation tasks:

- Add `genome_assembly/_cython_backend.pyx`.
- Add `genome_assembly/cython_backend.py`.
- Add `setup.py` with `cythonize()`.
- Add `perf` optional dependencies for `cython` and `setuptools`.
- Expose:
  - `count_kmers()`
  - `build_edges()`
  - `compact_contigs()`
- Add `AssemblyConfig(backend="cython")`.
- Add CLI `--backend cython`.
- Add parity tests against `backend="python"`.

Validation commands:

```bash
python -m pip install -e '.[perf]'
python setup.py build_ext --inplace
python -m unittest discover -s tests -v
ga assemble reads.fastq --k 31 --backend cython --outdir assembly_out
```

Acceptance criteria:

- Tests pass without Cython built, skipping Cython parity.
- Tests pass with Cython built, enabling Cython parity.
- Missing Cython extension has a clear build message.

## Stage 3: Rust Native Install Flow

Status: complete.

Goal: document and validate a local developer flow for installing the native
extension.

Implementation tasks:

- Maintain native install instructions in `README.md`.
- Maintain `native/genome_assembly_native/README.md`.
- Validate one local flow:

```bash
python -m pip install -e '.[native]'
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

## Stage 4: Rust Native Edge Table Integration

Status: complete.

Goal: move graph edge counting and edge table creation fully into Rust.

Implementation tasks:

- Extend Rust output to include `(prefix, suffix, sequence, count)`.
- Add Python parity tests comparing native and Python `GraphSummary`.
- Keep deterministic sort by edge sequence.
- Prepare edge creation for benchmark comparison on the bundled SARS-CoV-2
  fixture.

Acceptance criteria:

- Native and Python summaries match exactly.
- Native backend produces identical contig sequences on small fixtures.

## Stage 5: Rust Native Contig Compaction

Status: complete.

Goal: move maximal non-branching path compaction into Rust.

Implementation tasks:

- Implement Rust adjacency maps and one-in-one-out traversal.
- Add circular graph handling.
- Expose `compact_contigs()` through PyO3.
- Add Python native bridge support.
- Dispatch `backend="native"` compaction from `DeBruijnGraph`.
- Add Rust unit tests for linear, circular, short-filtered, and invalid-input
  compaction.
- Keep Python parity tests comparing native and Python backend outputs.

Validation commands:

```bash
cargo test --manifest-path native/genome_assembly_native/Cargo.toml
cd native/genome_assembly_native
maturin develop --release
cd ../..
python -m unittest discover -s tests -v
```

Acceptance criteria:

- Native contigs match Python contigs on all fixtures.
- Python objects and CLI output remain unchanged.

## Stage 6: Benchmark Harness

Status: complete for deterministic simulated-read benchmarks.

Goal: make performance claims reproducible.

Implementation tasks:

- Add `ga benchmark`.
- Add benchmark fixture generation.
- Record JSON output with command, config, wall time, traced Python memory,
  contig stats, and platform metadata.
- Add docs explaining benchmark tiers and limitations.

Validation commands:

```bash
ga benchmark genomic.fna --backends python,cython,native --coverage 5 --output benchmark.json
python -m json.tool benchmark.json >/dev/null
```

Acceptance criteria:

- Benchmarks can compare `backend="python"`, `backend="cython"`, and
  `backend="native"` when each backend is available.
- Results are written in machine-readable JSON.

## Stage 7: Multicore CPU/HPC

Status: complete for parallel native k-mer counting.

Goal: make `--threads` meaningful for native work.

Implementation tasks:

- Done: parallelize native read-chunk k-mer counting with Rayon, honoring a
  per-call thread pool sized by `threads`.
- Done: merge per-thread k-mer maps deterministically. Integer merge is
  associative/commutative and the final row order is sorted, so any thread count
  yields identical output.
- Done: thread `config.threads` through `native.build_edges` and
  `DeBruijnGraph.from_reads`; `ga assemble --threads` / `ga benchmark --threads`
  now drive the native backend.
- Done: Rust test `parallel_counts_match_single_thread` and Python test
  `test_native_multithread_output_matches_single_thread` prove `threads=1`
  equals `threads>1`.
- Deferred: memory-budget configuration. Reads are fully materialized into a
  Python list today, so a budget knob without streaming ingestion is theater.
  Add it together with chunked/streaming FASTQ I/O, not before.

Notes:

- Contig compaction stays single-threaded. It is an irregular graph traversal
  and not a measured bottleneck; revisit only if benchmarks flag it.
- `backend="cython"` stays single-process per ADR-003.

Validation commands:

```bash
cargo test --manifest-path native/genome_assembly_native/Cargo.toml
cd native/genome_assembly_native && maturin develop --release && cd ../..
python -m unittest discover -s tests -v
```

Acceptance criteria:

- Thread count changes runtime, not output. Verified on ~66k synthetic reads at
  k=21: about 1.7x wall-time reduction at 4 threads with identical k-mer counts.
- `--threads` is honored by `ga assemble` and `ga benchmark` for the native
  backend. Comparing two `ga benchmark` runs at different `--threads` gives a
  speedup measurement without a bespoke sweep command.

## Stage 8: PyPy Compatibility Lane

Status: blocked on toolchain. `pypy3` is not installed in the current
environment, so the pure-Python backend has not been validated under PyPy yet.
The only runtime dependency is Typer (pure-Python Click), which has no CPython-C
requirement, so the `backend="python"` path is expected to run on PyPy; the
Cython and Rust extension backends are CPython-oriented and are not advertised
as PyPy accelerators. Run the validation commands below once `pypy3` is
available to confirm and lift this block.

Goal: validate whether PyPy provides useful speedups for the pure-Python backend.

Implementation tasks:

- Install or locate `pypy3`.
- Run core tests under PyPy without compiled extensions.
- Document unsupported extension backends under PyPy.
- Add a PyPy benchmark row once benchmark harness exists.

Validation commands:

```bash
pypy3 -m pip install -e .
pypy3 -m unittest discover -s tests -v
pypy3 - <<'PY'
from genome_assembly import AssemblyConfig, assemble_short_reads
print(assemble_short_reads(["ACGTTA", "GTTACC"], AssemblyConfig(k=3)).stats())
PY
```

Acceptance criteria:

- Pure-Python backend works on PyPy, or blockers are documented.
- Cython and Rust extension backends are not advertised as PyPy accelerators.

## Stage 9: C++/CUDA Evaluation Lane

Status: decided — deferred. No benchmark bottleneck currently justifies the C++
build/ABI cost. Recorded in `docs/architecture/adr-004-defer-cpp-and-gpu-lanes.md`
with explicit re-entry criteria.

Goal: decide whether C++ should enter as a specialist backend.

Implementation tasks:

- Identify one bottleneck from benchmark output.
- Decide whether Rust, Cython, or C++ is the right target for that bottleneck.
- If C++ is justified, add a narrow pybind11/CMake spike under a separate
  backend name instead of replacing Rust.
- Add parity tests before exposing the C++ path to the CLI.

Acceptance criteria:

- C++ enters only with a written benchmark reason.
- The spike can be removed without affecting Python, Cython, or Rust backends.

## Stage 10: Graph Cleaning

Status: complete.

Goal: improve short-read assembly quality.

Implementation tasks:

- Done: tip clipping. `DeBruijnGraph._clip_tips` removes short dead-end branches
  anchored at a junction, iterating to convergence. A pure linear path is never
  a tip, so clean assemblies are preserved at any threshold.
- Done: simple bubble popping. `DeBruijnGraph._pop_bubbles` collapses parallel
  simple paths between the same two nodes, keeping the highest-coverage path
  (deterministic tie-break by sequence).
- Done: reportable removed-edge counts. `GraphSummary` gains `tips_removed` and
  `bubble_edges_removed`; the abundance-filter count is snapshotted at build
  time so cleaning does not inflate it.
- Done: conservative CLI flags. `--tip-length` and `--bubble-length` both
  default to 0 (disabled), so default assembly output is unchanged.

Design notes:

- Cleaning runs on the shared Python edge list inside `assemble_short_reads`,
  after edge construction and before compaction, so python/cython/native all
  benefit from one implementation. `test_cleaning_is_backend_agnostic` guards
  this.
- Tip clipping is length-thresholded only; keep the threshold small (about 2k)
  so the longer arm of a real branch is not clipped.

Validation commands:

```bash
python -m unittest discover -s tests -v
ga assemble reads.fastq --k 31 --tip-length 62 --bubble-length 130 --outdir assembly_out
```

Acceptance criteria:

- Tests show error k-mers are removed without breaking clean linear assemblies.
  Covered by `GraphCleaningTests`.
- Summary reports graph-cleaning counts. `summary.json` graph block carries
  `tips_removed` and `bubble_edges_removed`.

## Stage 11: GPU Research Spike

Status: decided — deferred. No profiled bottleneck yet warrants a GPU kernel,
and irregular graph traversal is a poor GPU fit. Recorded in
`docs/architecture/adr-004-defer-cpp-and-gpu-lanes.md` with re-entry criteria.

Goal: decide whether GPU acceleration is worth implementing.

Implementation tasks:

- Profile Python, Cython, and native CPU backends first.
- Identify the top two bottlenecks.
- Prototype only one GPU-friendly kernel, likely k-mer sorting/counting or
  minimizer bucketing.
- Compare against CPU native backend on wall time, transfer overhead, memory, and
  reproducibility.

Acceptance criteria:

- Written decision: continue CUDA backend, defer it, or drop it.
- No production GPU dependency unless the benchmark justifies it.

## Stage 12: Release Hardening

Status: complete.

Goal: make the project usable by others.

Implementation tasks:

- Done: CI (`.github/workflows/ci.yml`) runs Python tests on 3.10-3.12, Rust
  tests, clippy, `cargo fmt --check`, and a full-parity job that builds the
  Cython and native extensions.
- Done: release workflow (`.github/workflows/release.yml`) builds the
  pure-Python sdist/wheel and native wheels for Linux/macOS/Windows as
  artifacts. Publishing is left as a manual, human-authorized step.
- Done: examples in `README.md` and `docs/TROUBLESHOOTING.md`.
- Done: `CHANGELOG.md` with a stated `0.x` versioning policy.
- Done: `LICENSE` (MIT), matching `pyproject.toml`.

Acceptance criteria:

- Fresh clone can run all documented commands. CI exercises the documented
  install and test flows on a clean runner.
- Release artifacts are reproducible. `python -m build` produces the
  pure-Python wheel/sdist (verified locally); native wheels come from maturin.
- Docs state which backends are stable and which are experimental. See the
  backend status table in `docs/TROUBLESHOOTING.md`.
