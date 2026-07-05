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
- a native backend seam for future Rust/C++/CUDA acceleration

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

## Python API

```python
from genome_assembly import AssemblyConfig, assemble_short_reads

reads = ["ACGTTGA", "GTTGACC"]
result = assemble_short_reads(reads, AssemblyConfig(k=3))

for contig in result.contigs:
    print(contig.name, contig.sequence)
```

## Native Acceleration Roadmap

The public API is intentionally separated from the execution backend. The
current backend is pure Python for correctness and testability. The next backend
should move k-mer counting, graph construction, and unitig compaction into Rust,
with optional C++/CUDA kernels later for GPU-heavy stages such as sorting,
minimizer bucketing, mapping, and consensus.

The Rust backend scaffold lives in `native/genome_assembly_native`.

```bash
cargo test --manifest-path native/genome_assembly_native/Cargo.toml
```

For local native-backend experimentation:

```bash
python -m pip install maturin
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

See `docs/ROADMAP.md` for the product outline and
`docs/IMPLEMENTATION_PLAN.md` for the staged engineering plan.
