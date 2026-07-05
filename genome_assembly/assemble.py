"""High-level assembly workflows."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .config import AssemblyConfig
from .graph import Contig, DeBruijnGraph, GraphSummary
from .metrics import assembly_stats


@dataclass(frozen=True)
class AssemblyResult:
    config: AssemblyConfig
    graph: DeBruijnGraph
    contigs: list[Contig]

    @property
    def summary(self) -> GraphSummary:
        return self.graph.summary()

    def stats(self, *, reference_length: int | None = None) -> dict[str, int | float | None]:
        return assembly_stats([contig.sequence for contig in self.contigs], reference_length=reference_length)

    def to_report(self, *, reference_length: int | None = None) -> dict[str, object]:
        return {
            "config": asdict(self.config),
            "graph": asdict(self.summary),
            "stats": self.stats(reference_length=reference_length),
        }


def assemble_short_reads(reads: list[str], config: AssemblyConfig | None = None) -> AssemblyResult:
    """Assemble short reads with a de Bruijn graph."""

    config = config or AssemblyConfig()
    config.validate()
    if config.backend == "native":
        raise NotImplementedError("The native backend is not implemented yet; use backend='python'")

    graph = DeBruijnGraph.from_reads(reads, config)
    contigs = graph.compact_contigs(min_length=config.min_contig_length)
    return AssemblyResult(config=config, graph=graph, contigs=contigs)
