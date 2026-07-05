# genome_assembly_native

Rust/PyO3 native backend scaffold for `genome-assembly`.

The crate is intentionally separate from the default Python package build so the
pure-Python install remains simple while native acceleration matures.

## Test

```bash
cargo test --manifest-path native/genome_assembly_native/Cargo.toml
```

## Local Python Install

From the repository root:

```bash
python -m pip install maturin
source .venv/bin/activate
unset CONDA_PREFIX  # needed when this shell was launched from Conda
cd native/genome_assembly_native
maturin develop --release
cd ../..
```

Then:

```python
from genome_assembly.native import native_available, count_kmers

print(native_available())
print(count_kmers(["ACGTACGT"], 3))
```

The installed extension also supports the internal graph-building primitive:

```python
from genome_assembly.native import build_edges

print(build_edges(["ACGTTA", "GTTACC"], 3))
```
