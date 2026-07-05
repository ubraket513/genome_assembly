# Troubleshooting and Backend Status

## Backend stability

| Backend             | Status       | Install                         | Notes |
|---------------------|--------------|---------------------------------|-------|
| `python` (default)  | Stable       | `pip install -e .`              | Correctness oracle. Pure Python, no build step, PyPy-friendly. |
| `cython`            | Stable       | `pip install -e '.[perf]'` then `python setup.py build_ext --inplace` | CPython only. Tactical accelerator. |
| `native` (Rust)     | Stable       | `pip install -e '.[native]'` then `maturin develop --release` in the crate | Production native core. Only backend that uses `--threads`. |
| C++ / GPU           | Not shipped  | —                               | Deferred; see `docs/architecture/adr-004-defer-cpp-and-gpu-lanes.md`. |

If a compiled backend is requested but not built, the CLI and API raise an
actionable error telling you exactly which command to run. The pure-Python
backend always works, and the test suite skips compiled-backend parity checks
when the extensions are absent.

## Common issues

### `Native backend requested, but genome_assembly_native is not installed`

Build the Rust extension:

```bash
python -m pip install -e '.[native]'
cd native/genome_assembly_native && maturin develop --release && cd ../..
```

### `maturin develop` fails with a virtualenv/prefix error under Conda

If the shell was launched from Conda, `maturin` can pick up the wrong prefix.
Activate the project venv and clear `CONDA_PREFIX` first:

```bash
source .venv/bin/activate
unset CONDA_PREFIX
cd native/genome_assembly_native && maturin develop --release && cd ../..
```

### `Cython backend requested, but ... is not built`

```bash
python -m pip install -e '.[perf]'
python setup.py build_ext --inplace
```

### `No graph edges remain after filtering`

`k` is too large for the read length, or `--min-abundance` is too high for the
coverage. Lower `k`, lower `--min-abundance`, or provide longer/deeper reads.

### `--threads` does not change runtime

`--threads` only affects `backend="native"`. The `python` and `cython` backends
run single-process. It also has little effect on tiny inputs where the parallel
setup cost dominates; the speedup shows on larger read sets.

### Graph cleaning removed too much

`--tip-length` clips any short dead-end arm off a junction, so a value larger
than the real sequence length between branches will erase real contigs. Keep it
small, around `2 * k`. Both `--tip-length` and `--bubble-length` default to 0
(disabled).

## Reproducibility

All outputs are deterministically ordered. Simulation is seeded (`--seed`), and
every backend and thread count produces byte-identical contigs for the same
inputs. `ga benchmark` writes machine-readable JSON with platform metadata,
wall time, and traced Python memory (native RSS accounting is future work).
