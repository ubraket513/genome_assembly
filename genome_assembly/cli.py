"""Command-line interface for genome assembly workflows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer

from . import __version__
from .assemble import assemble_short_reads
from .config import AssemblyConfig
from .io import FastaRecord, FastqRecord, read_fasta, read_sequences, write_fasta, write_fastq
from .metrics import assembly_stats
from .simulate import simulate_reads

app = typer.Typer(
    help="Genome assembly library and CLI.",
    no_args_is_help=True,
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"genome-assembly {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool | None,
        typer.Option("--version", callback=_version_callback, is_eager=True, help="Show version and exit."),
    ] = None,
) -> None:
    """Genome assembly workflows."""


@app.command()
def simulate(
    reference: Annotated[Path, typer.Argument(exists=True, readable=True, help="Reference FASTA file.")],
    output: Annotated[Path, typer.Argument(help="Output FASTQ path.")],
    read_length: Annotated[int, typer.Option("--read-length", "-l", min=1)] = 150,
    coverage: Annotated[float, typer.Option("--coverage", "-c", min=0.0001)] = 30.0,
    seed: Annotated[int | None, typer.Option("--seed", help="Random seed for reproducibility.")] = 7,
) -> None:
    """Simulate fixed-length reads from the first FASTA record."""

    records = read_fasta(reference)
    simulated = simulate_reads(records[0].sequence, read_length=read_length, coverage=coverage, seed=seed)
    fastq_records = [
        FastqRecord(f"read_{index}", read, "I" * len(read))
        for index, read in enumerate(simulated.reads, start=1)
    ]
    write_fastq(fastq_records, output)
    mean_coverage = sum(simulated.coverage) / len(simulated.coverage)
    typer.echo(f"Wrote {len(simulated.reads)} reads to {output} (mean coverage {mean_coverage:.2f}x)")


@app.command()
def assemble(
    reads: Annotated[Path, typer.Argument(exists=True, readable=True, help="Input FASTA or FASTQ reads.")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory.")] = Path("assembly_out"),
    k: Annotated[int, typer.Option("--k", "-k", min=1, help="Node k-mer size.")] = 31,
    min_abundance: Annotated[
        int,
        typer.Option("--min-abundance", "-m", min=1, help="Drop edges below this abundance."),
    ] = 1,
    min_contig_length: Annotated[
        int,
        typer.Option("--min-contig-length", min=0, help="Suppress contigs shorter than this length."),
    ] = 0,
    threads: Annotated[int, typer.Option("--threads", "-t", min=1, help="Reserved for native backends.")] = 1,
    emit_gfa: Annotated[bool, typer.Option("--emit-gfa", help="Write graph.gfa alongside contigs.")] = False,
) -> None:
    """Assemble short reads into contigs."""

    outdir.mkdir(parents=True, exist_ok=True)
    sequences = read_sequences(reads)
    config = AssemblyConfig(
        k=k,
        min_abundance=min_abundance,
        min_contig_length=min_contig_length,
        threads=threads,
    )
    result = assemble_short_reads(sequences, config)

    contig_records = [
        FastaRecord(
            contig.name,
            contig.sequence,
            f"length={contig.length} mean_abundance={contig.mean_abundance:.3f} edges={contig.edge_count}",
        )
        for contig in result.contigs
    ]
    write_fasta(contig_records, outdir / "contigs.fasta")

    if emit_gfa:
        (outdir / "graph.gfa").write_text(result.graph.to_gfa(), encoding="utf-8")

    report = result.to_report()
    (outdir / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    stats = report["stats"]
    typer.echo(
        f"Wrote {stats['contigs']} contigs to {outdir / 'contigs.fasta'} "
        f"(N50={stats['n50']}, total={stats['total_bp']} bp)"
    )


@app.command()
def stats(
    contigs: Annotated[Path, typer.Argument(exists=True, readable=True, help="Contigs FASTA file.")],
    reference: Annotated[
        Path | None,
        typer.Option("--reference", "-r", exists=True, readable=True, help="Optional reference FASTA."),
    ] = None,
) -> None:
    """Compute assembly statistics for contigs."""

    contig_records = read_fasta(contigs)
    reference_length = None
    if reference is not None:
        reference_length = sum(len(record.sequence) for record in read_fasta(reference))
    report = assembly_stats([record.sequence for record in contig_records], reference_length=reference_length)
    typer.echo(json.dumps(report, indent=2))


if __name__ == "__main__":
    app()
