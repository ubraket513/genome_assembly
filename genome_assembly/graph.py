"""Dependency-light de Bruijn graph primitives."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

from .config import AssemblyConfig
from .cython_backend import build_edges as cython_build_edges
from .cython_backend import compact_contigs as cython_compact_contigs
from .kmers import iter_kmers
from .native import build_edges as native_build_edges
from .native import compact_contigs as native_compact_contigs


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
    tips_removed: int = 0
    bubble_edges_removed: int = 0


class DeBruijnGraph:
    """A compact in-memory directed de Bruijn graph."""

    def __init__(
        self,
        k: int,
        edges: list[Edge],
        raw_edge_count: int,
        min_abundance: int,
        backend: str = "python",
        threads: int = 1,
    ) -> None:
        self.k = k
        self.edges = edges
        self.raw_edge_count = raw_edge_count
        self.min_abundance = min_abundance
        self.backend = backend
        self.threads = threads
        # Snapshot of retained edge observations at build time. Kept separate so
        # graph cleaning does not inflate the abundance-filter count in summary().
        self._built_edge_observations = sum(edge.count for edge in edges)
        self._tips_removed = 0
        self._bubble_edges_removed = 0
        self._reindex()

    def _reindex(self) -> None:
        self._out_edges: dict[str, list[Edge]] = defaultdict(list)
        self._in_edges: dict[str, list[Edge]] = defaultdict(list)
        self._nodes: set[str] = set()
        for edge in self.edges:
            self._out_edges[edge.prefix].append(edge)
            self._in_edges[edge.suffix].append(edge)
            self._nodes.add(edge.prefix)
            self._nodes.add(edge.suffix)

    @classmethod
    def from_reads(cls, reads: list[str], config: AssemblyConfig) -> "DeBruijnGraph":
        config.validate()
        edge_size = config.k + 1

        if config.backend == "cython":
            raw_edge_count, edge_rows = cython_build_edges(
                reads,
                config.k,
                min_abundance=config.min_abundance,
                skip_ambiguous=config.skip_ambiguous,
            )
            edges = [
                Edge(prefix, suffix, sequence, count)
                for prefix, suffix, sequence, count in edge_rows
            ]
        elif config.backend == "native":
            raw_edge_count, edge_rows = native_build_edges(
                reads,
                config.k,
                min_abundance=config.min_abundance,
                skip_ambiguous=config.skip_ambiguous,
                threads=config.threads,
            )
            edges = [
                Edge(prefix, suffix, sequence, count)
                for prefix, suffix, sequence, count in edge_rows
            ]
        else:
            python_counts: Counter[str] = Counter()
            for read in reads:
                python_counts.update(iter_kmers(read, edge_size, skip_ambiguous=config.skip_ambiguous))

            raw_edge_count = sum(python_counts.values())
            edges = [
                Edge(kmer[:-1], kmer[1:], kmer, count)
                for kmer, count in sorted(python_counts.items())
                if count >= config.min_abundance
            ]

        if not edges:
            raise ValueError(
                "No graph edges remain after filtering; lower k/min_abundance or provide longer reads"
            )

        return cls(
            config.k, edges, raw_edge_count, config.min_abundance, config.backend, config.threads
        )

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
        return GraphSummary(
            k=self.k,
            nodes=len(self._nodes),
            edges=len(self.edges),
            raw_edges=self.raw_edge_count,
            filtered_edges=max(self.raw_edge_count - self._built_edge_observations, 0),
            min_abundance=self.min_abundance,
            tips_removed=self._tips_removed,
            bubble_edges_removed=self._bubble_edges_removed,
        )

    def _is_one_in_one_out(self, node: str) -> bool:
        return self.in_degree(node) == 1 and self.out_degree(node) == 1

    @staticmethod
    def _edge_key(edge: Edge) -> tuple[str, str]:
        return edge.prefix, edge.suffix

    def _chain_sequence(self, chain: list[Edge]) -> str:
        sequence = chain[0].sequence
        for edge in chain[1:]:
            sequence += edge.sequence[-1]
        return sequence

    def _remove_edges(self, keys: set[tuple[str, str]]) -> int:
        if not keys:
            return 0
        self.edges = [edge for edge in self.edges if self._edge_key(edge) not in keys]
        self._reindex()
        return len(keys)

    def clean(self, *, tip_length: int = 0, bubble_length: int = 0) -> None:
        """Clip short dead-end tips and pop simple bubbles in place.

        Both passes are no-ops at 0, so default assembly behavior is unchanged.
        Cleaning runs on the shared Python edge list, so every backend benefits.
        """

        if tip_length > 0:
            self._tips_removed += self._clip_tips(tip_length)
        if bubble_length > 0:
            self._bubble_edges_removed += self._pop_bubbles(bubble_length)

    def _iter_unitigs(self) -> list[list[Edge]]:
        """Yield maximal non-branching edge chains (unitigs)."""

        visited: set[tuple[str, str]] = set()
        unitigs: list[list[Edge]] = []

        def walk(start_edge: Edge) -> list[Edge]:
            chain = [start_edge]
            visited.add(self._edge_key(start_edge))
            current = start_edge.suffix
            while self._is_one_in_one_out(current):
                next_edges = self.out_edges(current)
                if not next_edges:
                    break
                next_edge = next_edges[0]
                if self._edge_key(next_edge) in visited:
                    break
                visited.add(self._edge_key(next_edge))
                chain.append(next_edge)
                current = next_edge.suffix
            return chain

        for node in sorted(self._nodes):
            if self._is_one_in_one_out(node):
                continue
            for edge in self.out_edges(node):
                if self._edge_key(edge) not in visited:
                    unitigs.append(walk(edge))
        for edge in self.edges:  # remaining pure cycles
            if self._edge_key(edge) not in visited:
                unitigs.append(walk(edge))
        return unitigs

    def _clip_tips(self, max_tip_bp: int) -> int:
        """Remove short dead-end branches anchored at a junction, to convergence."""

        removed = 0
        while True:
            remove_keys: set[tuple[str, str]] = set()
            for chain in self._iter_unitigs():
                start_node = chain[0].prefix
                end_node = chain[-1].suffix
                length_bp = self.k + len(chain)
                if length_bp >= max_tip_bp:
                    continue
                is_source_tip = self.in_degree(start_node) == 0 and self.in_degree(end_node) >= 2
                is_sink_tip = self.out_degree(end_node) == 0 and self.out_degree(start_node) >= 2
                if is_source_tip or is_sink_tip:
                    remove_keys.update(self._edge_key(edge) for edge in chain)
            if not remove_keys:
                return removed
            removed += self._remove_edges(remove_keys)

    def _walk_simple_path(self, start_edge: Edge, max_bp: int) -> list[Edge] | None:
        """Walk a simple path through 1-in-1-out interiors up to max_bp; None if it loops or overruns."""

        chain = [start_edge]
        seen = {self._edge_key(start_edge)}
        current = start_edge.suffix
        while self._is_one_in_one_out(current):
            next_edges = self.out_edges(current)
            if not next_edges:
                break
            next_edge = next_edges[0]
            key = self._edge_key(next_edge)
            if key in seen:
                return None
            seen.add(key)
            chain.append(next_edge)
            current = next_edge.suffix
            if self.k + len(chain) > max_bp:
                return None
        return chain

    def _pop_bubbles(self, max_bubble_bp: int) -> int:
        """Collapse parallel simple paths between the same two nodes, keeping the highest-coverage path."""

        removed = 0
        while True:
            popped = self._pop_one_bubble(max_bubble_bp)
            if popped == 0:
                return removed
            removed += popped

    def _pop_one_bubble(self, max_bubble_bp: int) -> int:
        for node in sorted(self._nodes):
            if self.out_degree(node) < 2:
                continue
            paths_by_end: dict[str, list[list[Edge]]] = defaultdict(list)
            for edge in self.out_edges(node):
                chain = self._walk_simple_path(edge, max_bubble_bp)
                if chain is not None:
                    paths_by_end[chain[-1].suffix].append(chain)
            for end_node, group in paths_by_end.items():
                if end_node == node or len(group) < 2:
                    continue
                group.sort(
                    key=lambda chain: (
                        -sum(edge.count for edge in chain) / len(chain),
                        self._chain_sequence(chain),
                    )
                )
                keep_keys = {self._edge_key(edge) for edge in group[0]}
                remove_keys: set[tuple[str, str]] = set()
                for chain in group[1:]:
                    remove_keys.update(self._edge_key(edge) for edge in chain)
                remove_keys -= keep_keys
                if remove_keys:
                    return self._remove_edges(remove_keys)
        return 0

    def compact_contigs(self, *, min_length: int = 0) -> list[Contig]:
        """Return maximal non-branching paths as contigs."""

        if self.backend == "cython":
            contig_rows = cython_compact_contigs(
                self.k,
                [(edge.prefix, edge.suffix, edge.sequence, edge.count) for edge in self.edges],
                min_length=min_length,
            )
            return [
                Contig(
                    name=f"contig_{index}",
                    sequence=sequence,
                    mean_abundance=mean_abundance,
                    edge_count=edge_count,
                )
                for index, (sequence, mean_abundance, edge_count) in enumerate(contig_rows, start=1)
            ]
        if self.backend == "native":
            contig_rows = native_compact_contigs(
                self.k,
                [(edge.prefix, edge.suffix, edge.sequence, edge.count) for edge in self.edges],
                min_length=min_length,
                threads=self.threads,
            )
            return [
                Contig(
                    name=f"contig_{index}",
                    sequence=sequence,
                    mean_abundance=mean_abundance,
                    edge_count=edge_count,
                )
                for index, (sequence, mean_abundance, edge_count) in enumerate(contig_rows, start=1)
            ]

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

        contigs.sort(key=lambda contig: (-contig.length, contig.sequence))
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
