"""Command-line interface for genome assembly workflows."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from textwrap import dedent
from typing import Annotated

import typer

from . import __version__
from .assemble import assemble_short_reads
from .benchmark import parse_backend_list, run_benchmark
from .config import AssemblyConfig
from .io import FastaRecord, FastqRecord, read_fasta, read_sequences, write_fasta, write_fastq
from .metrics import assembly_stats
from .simulate import simulate_reads

_EPILOG = dedent(
    """\
    Examples:
      ga simulate genome.fna reads.fastq --coverage 30
      ga assemble reads.fastq --k 31 --outdir out
      ga stats out/contigs.fasta --reference genome.fna

    New here? Run:  ga quickstart
    """
)

app = typer.Typer(
    help="Genome assembly library and CLI. Turn reads into contigs in three steps.",
    no_args_is_help=True,
    epilog=_EPILOG,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _fail(message: str) -> None:
    """Print a clean one-line error and exit, instead of a traceback."""

    typer.secho(f"Error: {message}", fg=typer.colors.RED, err=True)
    raise typer.Exit(code=1)


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
def shell() -> None:
    """Launch the interactive oxidas shell (agentic REPL)."""

    from .agent.repl import main as run_shell

    run_shell()


@app.command()
def quickstart() -> None:
    """Print a copy-paste walkthrough for first-time users."""

    typer.echo(
        dedent(
            """\
            Genome assembly in three steps
            ==============================

            You need one input: a reads file (FASTA or FASTQ). If you only have a
            reference genome, simulate reads from it first.

            1) (optional) Simulate reads from a reference genome:
                 ga simulate genome.fna reads.fastq --coverage 30

            2) Assemble the reads into contigs:
                 ga assemble reads.fastq --k 31 --outdir out
               Outputs: out/contigs.fasta and out/summary.json

            3) Check how good the assembly is:
                 ga stats out/contigs.fasta --reference genome.fna

            Helpful extras
            --------------
              --backend native      faster assembly (needs the built Rust extension)
              --threads 4           use 4 cores (native backend only)
              --tip-length 62       clean short error branches (~2*k is a safe value)
              ga benchmark ...      compare backend speeds on simulated reads

            Every command has its own help, for example:
                 ga assemble --help
            """
        )
    )


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
    typer.secho(f"Next: ga assemble {output} --outdir assembly_out", fg=typer.colors.GREEN)


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
    backend: Annotated[
        str,
        typer.Option("--backend", help="Execution backend: python, cython, or native."),
    ] = "python",
    threads: Annotated[int, typer.Option("--threads", "-t", min=1, help="Worker threads for the native backend (ignored by python/cython).")] = 1,
    tip_length: Annotated[
        int,
        typer.Option("--tip-length", min=0, help="Clip dead-end tips shorter than this many bp (0 disables)."),
    ] = 0,
    bubble_length: Annotated[
        int,
        typer.Option("--bubble-length", min=0, help="Pop simple bubbles up to this many bp, keeping the highest-coverage path (0 disables)."),
    ] = 0,
    emit_gfa: Annotated[bool, typer.Option("--emit-gfa", help="Write graph.gfa alongside contigs.")] = False,
) -> None:
    """Assemble short reads into contigs."""

    outdir.mkdir(parents=True, exist_ok=True)
    sequences = read_sequences(reads)
    if not sequences:
        _fail(f"No reads found in {reads}. Is it a valid FASTA/FASTQ file?")
    config = AssemblyConfig(
        k=k,
        min_abundance=min_abundance,
        min_contig_length=min_contig_length,
        backend=backend,
        threads=threads,
        tip_length=tip_length,
        bubble_length=bubble_length,
    )
    try:
        result = assemble_short_reads(sequences, config)
    except ValueError as exc:
        # e.g. k too large for the reads, or min-abundance too high for coverage.
        _fail(f"{exc}")
    except RuntimeError as exc:
        # e.g. a compiled backend was requested but is not built; message is actionable.
        _fail(f"{exc}")

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
    graph = report["graph"]
    if graph["tips_removed"] or graph["bubble_edges_removed"]:
        typer.echo(
            f"Graph cleaning removed {graph['tips_removed']} tip edges "
            f"and {graph['bubble_edges_removed']} bubble edges"
        )
    typer.secho(f"Next: ga stats {outdir / 'contigs.fasta'}", fg=typer.colors.GREEN)


@app.command()
def mdbg(
    reads: Annotated[Path, typer.Argument(exists=True, readable=True, help="Long/accurate reads (FASTA or FASTQ).")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory.")] = Path("mdbg_out"),
    window: Annotated[int, typer.Option("--window", "-w", min=1, help="Minimizer window (w).")] = 5,
    minimizer_length: Annotated[int, typer.Option("--minimizer-length", "-l", min=1, help="Minimizer length (m).")] = 8,
    kmin: Annotated[int, typer.Option("--kmin", min=2, help="Minimizers per k-min-mer (k).")] = 3,
    min_contig_length: Annotated[int, typer.Option("--min-contig-length", min=0, help="Suppress short contigs.")] = 0,
) -> None:
    """Assemble long, accurate reads in minimizer space (mdBG, native only)."""

    from .native import mdbg_assemble, native_available

    if not native_available():
        _fail(
            "mdBG needs the native Rust extension. Build it with: "
            "cd native/genome_assembly_native && maturin develop --release"
        )
    sequences = read_sequences(reads)
    if not sequences:
        _fail(f"No reads found in {reads}.")
    try:
        contigs = mdbg_assemble(sequences, window, minimizer_length, kmin, min_length=min_contig_length)
    except (ValueError, RuntimeError) as exc:
        _fail(str(exc))

    outdir.mkdir(parents=True, exist_ok=True)
    records = [
        FastaRecord(f"contig_{index}", sequence, f"length={len(sequence)} mean_abundance={abundance:.3f}")
        for index, (sequence, abundance) in enumerate(contigs, start=1)
    ]
    write_fasta(records, outdir / "contigs.fasta")
    total = sum(len(sequence) for sequence, _ in contigs)
    typer.echo(f"Wrote {len(contigs)} mdBG contigs to {outdir / 'contigs.fasta'} (total {total} bp)")
    if contigs:
        typer.secho(f"Next: ga stats {outdir / 'contigs.fasta'}", fg=typer.colors.GREEN)


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


@app.command()
def benchmark(
    reference: Annotated[
        Path | None,
        typer.Argument(exists=True, readable=True, help="Reference FASTA file (omit if using --genome-size)."),
    ] = None,
    genome_size: Annotated[
        int | None,
        typer.Option("--genome-size", min=1, help="Benchmark on a synthetic random genome of this many bp instead of a reference file. Enables bacterial-scale tiers without downloads."),
    ] = None,
    output: Annotated[Path, typer.Option("--output", "-o", help="Benchmark JSON output path.")] = Path(
        "benchmark.json"
    ),
    backends: Annotated[
        str,
        typer.Option("--backends", help="Comma-separated backends: python,cython,native."),
    ] = "python,cython,native",
    k: Annotated[int, typer.Option("--k", "-k", min=1, help="Node k-mer size.")] = 31,
    read_length: Annotated[int, typer.Option("--read-length", "-l", min=1)] = 150,
    coverage: Annotated[float, typer.Option("--coverage", "-c", min=0.0001)] = 5.0,
    seed: Annotated[int | None, typer.Option("--seed", help="Random seed for reproducibility.")] = 7,
    min_abundance: Annotated[
        int,
        typer.Option("--min-abundance", "-m", min=1, help="Drop edges below this abundance."),
    ] = 1,
    min_contig_length: Annotated[
        int,
        typer.Option("--min-contig-length", min=0, help="Suppress contigs shorter than this length."),
    ] = 0,
    threads: Annotated[int, typer.Option("--threads", "-t", min=1, help="Worker threads for the native backend (ignored by python/cython).")] = 1,
) -> None:
    """Run deterministic backend benchmarks on simulated reads."""

    if (reference is None) == (genome_size is None):
        _fail("provide exactly one of a REFERENCE file argument or --genome-size")

    report = run_benchmark(
        reference,
        synthetic_genome_size=genome_size,
        backends=parse_backend_list(backends),
        k=k,
        read_length=read_length,
        coverage=coverage,
        seed=seed,
        min_abundance=min_abundance,
        min_contig_length=min_contig_length,
        threads=threads,
        command=sys.argv,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    runs = report["runs"]
    ok_count = sum(1 for run in runs if run["status"] == "ok")
    typer.echo(f"Wrote benchmark report to {output} ({ok_count}/{len(runs)} backends completed)")


if __name__ == "__main__":
    app()
