"""Production-oriented genome assembly library."""

from .assemble import AssemblyResult, assemble_short_reads
from .benchmark import run_benchmark
from .config import AssemblyConfig
from .cython_backend import cython_available
from .graph import Contig, DeBruijnGraph, GraphSummary
from .io import FastaRecord, FastqRecord, read_fasta, read_fastq, write_fasta, write_fastq
from .metrics import assembly_stats, n50, nx
from .native import native_available
from .simulate import SimulatedReads, generate_random_genome, simulate_reads
from .sketch import minimizer_bucket, minimizers, syncmers

__version__ = "0.1.0"

__all__ = [
    "AssemblyConfig",
    "AssemblyResult",
    "Contig",
    "cython_available",
    "DeBruijnGraph",
    "FastaRecord",
    "FastqRecord",
    "GraphSummary",
    "SimulatedReads",
    "assemble_short_reads",
    "assembly_stats",
    "generate_random_genome",
    "minimizer_bucket",
    "minimizers",
    "n50",
    "native_available",
    "nx",
    "syncmers",
    "read_fasta",
    "read_fastq",
    "run_benchmark",
    "simulate_reads",
    "write_fasta",
    "write_fastq",
]
