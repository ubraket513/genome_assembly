# Genome Assembly

This repository is being shaped from a small de Bruijn graph prototype into a
production-oriented genome assembly library and CLI.

The first implementation slice focuses on a correct, dependency-light short-read
assembly core:

- FASTA/FASTQ parsing and writing
- deterministic read simulation
- de Bruijn graph construction from `(k + 1)`-mers
- maximal non-branching path compaction into contigs
- N50/Nx assembly metrics
- a Typer-powered CLI
- backend selection for Python, Cython, and Rust native acceleration

## Install

```bash
python -m pip install -e .
```

## CLI

Simulate reads from the bundled SARS-CoV-2 reference:

```bash
ga simulate genomic.fna reads.fastq --read-length 150 --coverage 30 --seed 7
```

Assemble reads:

```bash
ga assemble reads.fastq --k 31 --min-abundance 1 --outdir assembly_out --emit-gfa
```

Compute contig statistics:

```bash
ga stats assembly_out/contigs.fasta --reference genomic.fna
```

Benchmark available backends on deterministic simulated reads:

```bash
ga benchmark genomic.fna --backends python,cython,native --coverage 5 --output benchmark.json
```

## Python API

```python
from genome_assembly import AssemblyConfig, assemble_short_reads

reads = ["ACGTTGA", "GTTGACC"]
result = assemble_short_reads(reads, AssemblyConfig(k=3))

for contig in result.contigs:
    print(contig.name, contig.sequence)
```

## Performance Backends

The public API is intentionally separated from the execution backend.

- `backend="python"` is the default correctness oracle and easiest install.
- `backend="cython"` is the tactical CPython accelerator.
- `backend="native"` is the preferred production native path, implemented with
  Rust/PyO3.

Build the Cython accelerator:

```bash
python -m pip install -e '.[perf]'
python setup.py build_ext --inplace
```

Use it from the CLI:

```bash
ga assemble reads.fastq --k 31 --backend cython --outdir assembly_out
```

Use it from Python:

```python
from genome_assembly import AssemblyConfig, assemble_short_reads

result = assemble_short_reads(reads, AssemblyConfig(k=31, backend="cython"))
```

PyPy should use `backend="python"` unless a compiled backend is explicitly
validated for PyPy.

The Rust backend lives in `native/genome_assembly_native` and is the primary
native core for k-mer counting, graph edge-table construction, and contig
compaction.

```bash
cargo test --manifest-path native/genome_assembly_native/Cargo.toml
```

For local native-backend experimentation:

```bash
python -m pip install -e '.[native]'
source .venv/bin/activate
unset CONDA_PREFIX  # needed when this shell was launched from Conda
cd native/genome_assembly_native
maturin develop --release
cd ../..
```

Then use:

```python
from genome_assembly import AssemblyConfig, assemble_short_reads

result = assemble_short_reads(reads, AssemblyConfig(k=31, backend="native"))
```

C++ remains a candidate for future specialist work through pybind11/CMake,
especially for existing C++ bioinformatics libraries or CUDA-heavy kernels. It
is not the default core until benchmarks justify that added build and memory
safety cost.

See `docs/ROADMAP.md` for the product outline and
`docs/IMPLEMENTATION_PLAN.md` for the staged engineering plan.
