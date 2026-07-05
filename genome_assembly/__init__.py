"""Production-oriented genome assembly library."""

from .assemble import AssemblyResult, assemble_short_reads
from .config import AssemblyConfig
from .graph import Contig, DeBruijnGraph, GraphSummary
from .io import FastaRecord, FastqRecord, read_fasta, read_fastq, write_fasta, write_fastq
from .metrics import assembly_stats, n50, nx
from .simulate import SimulatedReads, simulate_reads

__version__ = "0.1.0"

__all__ = [
    "AssemblyConfig",
    "AssemblyResult",
    "Contig",
    "DeBruijnGraph",
    "FastaRecord",
    "FastqRecord",
    "GraphSummary",
    "SimulatedReads",
    "assemble_short_reads",
    "assembly_stats",
    "n50",
    "nx",
    "read_fasta",
    "read_fastq",
    "simulate_reads",
    "write_fasta",
    "write_fastq",
]
