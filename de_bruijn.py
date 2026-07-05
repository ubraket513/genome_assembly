"""Compatibility wrapper around the production package graph core."""

from genome_assembly.config import AssemblyConfig
from genome_assembly.graph import DeBruijnGraph
from genome_assembly.kmers import iter_kmers


def get_kmers(read, k: int):
    """Return k-mers for one read or a collection of reads."""

    if isinstance(read, str):
        return list(iter_kmers(read, k))
    kmers = []
    for item in read:
        sequence = "".join(item) if not isinstance(item, str) else item
        kmers.extend(iter_kmers(sequence, k))
    return kmers


def create(reads, k: int):
    """Create a dependency-light DeBruijnGraph from reads."""

    normalized_reads = ["".join(read) if not isinstance(read, str) else read for read in reads]
    return DeBruijnGraph.from_reads(normalized_reads, AssemblyConfig(k=k))


def is_node_1_to_1(graph, node):
    """Check whether a graph node has exactly one incoming and outgoing edge."""

    return graph.in_degree(node) == 1 and graph.out_degree(node) == 1


def generate_contigs(graph):
    """Return contig sequences from maximal non-branching paths."""

    return [contig.sequence for contig in graph.compact_contigs()]


def compress(graph):
    """Return compacted contigs for compatibility with the old prototype."""

    return graph.compact_contigs()
