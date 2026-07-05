"""Dependency-light de Bruijn graph primitives."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

from .config import AssemblyConfig
from .kmers import iter_kmers


@dataclass(frozen=True)
class Edge:
    prefix: str
    suffix: str
    sequence: str
    count: int


@dataclass(frozen=True)
class Contig:
    name: str
    sequence: str
    mean_abundance: float
    edge_count: int

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class GraphSummary:
    k: int
    nodes: int
    edges: int
    raw_edges: int
    filtered_edges: int
    min_abundance: int


class DeBruijnGraph:
    """A compact in-memory directed de Bruijn graph."""

    def __init__(self, k: int, edges: list[Edge], raw_edge_count: int, min_abundance: int) -> None:
        self.k = k
        self.edges = edges
        self.raw_edge_count = raw_edge_count
        self.min_abundance = min_abundance
        self._out_edges: dict[str, list[Edge]] = defaultdict(list)
        self._in_edges: dict[str, list[Edge]] = defaultdict(list)
        self._nodes: set[str] = set()

        for edge in edges:
            self._out_edges[edge.prefix].append(edge)
            self._in_edges[edge.suffix].append(edge)
            self._nodes.add(edge.prefix)
            self._nodes.add(edge.suffix)

    @classmethod
    def from_reads(cls, reads: list[str], config: AssemblyConfig) -> "DeBruijnGraph":
        config.validate()
        edge_size = config.k + 1
        counts: Counter[str] = Counter()

        for read in reads:
            counts.update(iter_kmers(read, edge_size, skip_ambiguous=config.skip_ambiguous))

        raw_edge_count = sum(counts.values())
        edges = [
            Edge(kmer[:-1], kmer[1:], kmer, count)
            for kmer, count in sorted(counts.items())
            if count >= config.min_abundance
        ]

        if not edges:
            raise ValueError(
                "No graph edges remain after filtering; lower k/min_abundance or provide longer reads"
            )

        return cls(config.k, edges, raw_edge_count, config.min_abundance)

    @property
    def nodes(self) -> set[str]:
        return set(self._nodes)

    def in_degree(self, node: str) -> int:
        return len(self._in_edges.get(node, []))

    def out_degree(self, node: str) -> int:
        return len(self._out_edges.get(node, []))

    def out_edges(self, node: str) -> list[Edge]:
        return list(self._out_edges.get(node, []))

    def summary(self) -> GraphSummary:
        retained_edge_observations = sum(edge.count for edge in self.edges)
        return GraphSummary(
            k=self.k,
            nodes=len(self._nodes),
            edges=len(self.edges),
            raw_edges=self.raw_edge_count,
            filtered_edges=max(self.raw_edge_count - retained_edge_observations, 0),
            min_abundance=self.min_abundance,
        )

    def _is_one_in_one_out(self, node: str) -> bool:
        return self.in_degree(node) == 1 and self.out_degree(node) == 1

    def compact_contigs(self, *, min_length: int = 0) -> list[Contig]:
        """Return maximal non-branching paths as contigs."""

        contigs: list[Contig] = []
        visited: set[tuple[str, str]] = set()

        def edge_key(edge: Edge) -> tuple[str, str]:
            return edge.prefix, edge.suffix

        def walk(start_edge: Edge) -> tuple[str, list[int], int]:
            visited.add(edge_key(start_edge))
            sequence = start_edge.sequence
            counts = [start_edge.count]
            current = start_edge.suffix
            edge_count = 1

            while self._is_one_in_one_out(current):
                next_edges = self.out_edges(current)
                if not next_edges:
                    break
                next_edge = next_edges[0]
                if edge_key(next_edge) in visited:
                    break
                visited.add(edge_key(next_edge))
                sequence += next_edge.sequence[-1]
                counts.append(next_edge.count)
                current = next_edge.suffix
                edge_count += 1

            return sequence, counts, edge_count

        for node in sorted(self._nodes):
            if self._is_one_in_one_out(node):
                continue
            for edge in self.out_edges(node):
                if edge_key(edge) in visited:
                    continue
                sequence, counts, edge_count = walk(edge)
                if len(sequence) >= min_length:
                    contigs.append(
                        Contig(
                            name=f"contig_{len(contigs) + 1}",
                            sequence=sequence,
                            mean_abundance=sum(counts) / len(counts),
                            edge_count=edge_count,
                        )
                    )

        for edge in self.edges:
            if edge_key(edge) in visited:
                continue
            sequence, counts, edge_count = walk(edge)
            if len(sequence) >= min_length:
                contigs.append(
                    Contig(
                        name=f"contig_{len(contigs) + 1}",
                        sequence=sequence,
                        mean_abundance=sum(counts) / len(counts),
                        edge_count=edge_count,
                    )
                )

        contigs.sort(key=lambda contig: (-contig.length, contig.name))
        return [
            Contig(
                name=f"contig_{index}",
                sequence=contig.sequence,
                mean_abundance=contig.mean_abundance,
                edge_count=contig.edge_count,
            )
            for index, contig in enumerate(contigs, start=1)
        ]

    def to_gfa(self) -> str:
        """Serialize the edge graph in simple GFA 1.0 form."""

        lines = ["H\tVN:Z:1.0"]
        for node in sorted(self._nodes):
            lines.append(f"S\t{node}\t{node}")
        for edge in self.edges:
            lines.append(f"L\t{edge.prefix}\t+\t{edge.suffix}\t+\t{self.k - 1}M\tKC:i:{edge.count}")
        return "\n".join(lines) + "\n"
