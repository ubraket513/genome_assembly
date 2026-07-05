"""Execute resolved plans, offloading blocking work off the UI thread.

`run_plan` is the blocking dispatcher into the core library (this is what would
call the heavy Rust/PyO3 code). `execute_plan` runs it inside a thread-pool
executor so the asyncio UI loop keeps animating a transient Rich spinner, then
hands back a Markdown summary. `run_plan`/`apply_summary` are Rich-free and can
be driven headlessly in tests.
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

from ..assemble import assemble_short_reads
from ..benchmark import run_benchmark
from ..config import AssemblyConfig
from ..io import FastaRecord, FastqRecord, read_fasta, read_sequences, write_fasta, write_fastq
from ..metrics import assembly_stats
from ..simulate import simulate_reads
from .intents import Plan
from .state import SessionState


def run_plan(plan: Plan) -> dict:
    """Blocking dispatch into the core library. Runs inside a worker thread."""

    kind = plan.kind
    p = plan.params

    if kind == "simulate":
        records = read_fasta(p["reference"])
        simulated = simulate_reads(records[0].sequence, coverage=p["coverage"], seed=7)
        fastq = [
            FastqRecord(f"read_{i}", read, "I" * len(read))
            for i, read in enumerate(simulated.reads, start=1)
        ]
        write_fastq(fastq, p["output"])
        return {
            "reads": len(simulated.reads),
            "output": str(p["output"]),
            "mean_coverage": sum(simulated.coverage) / len(simulated.coverage),
        }

    if kind == "assemble":
        outdir: Path = p["outdir"]
        outdir.mkdir(parents=True, exist_ok=True)
        sequences = read_sequences(p["reads"])
        config = AssemblyConfig(
            k=p["k"],
            backend=p["backend"],
            threads=p["threads"],
            tip_length=p["tip_length"],
            bubble_length=p["bubble_length"],
        )
        result = assemble_short_reads(sequences, config)
        contigs_path = outdir / "contigs.fasta"
        write_fasta(
            [FastaRecord(c.name, c.sequence, f"length={c.length}") for c in result.contigs],
            contigs_path,
        )
        report = result.to_report()
        (outdir / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        return {
            "outdir": str(outdir),
            "contigs_path": str(contigs_path),
            "stats": report["stats"],
            "graph": report["graph"],
        }

    if kind == "stats":
        contig_records = read_fasta(p["contigs"])
        reference_length = None
        if p["reference"] is not None:
            reference_length = sum(len(r.sequence) for r in read_fasta(p["reference"]))
        report = assembly_stats(
            [r.sequence for r in contig_records], reference_length=reference_length
        )
        return {"stats": report}

    if kind == "benchmark":
        report = run_benchmark(
            p["reference"],
            synthetic_genome_size=p["genome_size"],
            backends=p["backends"],
            coverage=p["coverage"],
            threads=p["threads"],
        )
        return {"report": report}

    raise ValueError(f"unknown plan kind: {kind}")


def apply_summary(plan: Plan, result: dict, state: SessionState) -> str:
    """Update session memory and return a Markdown summary of the result."""

    kind = plan.kind
    if kind == "simulate":
        state.last_reads = Path(result["output"])
        state.last_reference = plan.params["reference"]
        return (
            f"Wrote **{result['reads']}** reads to `{result['output']}` "
            f"(mean coverage {result['mean_coverage']:.1f}x).\n\n"
            f"Next: `assemble {result['output']}`"
        )

    if kind == "assemble":
        state.last_outdir = Path(result["outdir"])
        state.last_contigs = Path(result["contigs_path"])
        state.last_reads = plan.params["reads"]
        stats = result["stats"]
        graph = result["graph"]
        lines = [
            f"**{stats['contigs']}** contigs · N50 **{stats['n50']}** · "
            f"total **{stats['total_bp']} bp** → `{result['contigs_path']}`"
        ]
        if graph["tips_removed"] or graph["bubble_edges_removed"]:
            lines.append(
                f"Cleaning removed {graph['tips_removed']} tip edges and "
                f"{graph['bubble_edges_removed']} bubble edges."
            )
        lines.append("Next: `run stats on that output`")
        return "\n\n".join(lines)

    if kind == "stats":
        stats = result["stats"]
        return (
            f"N50 **{stats['n50']}** · NG50 **{stats['ng50']}** · "
            f"contigs **{stats['contigs']}** · coverage **{stats['coverage_percent']}%**"
        )

    if kind == "benchmark":
        rows = result["report"]["runs"]
        lines = ["| backend | status | wall (s) | peak RSS (MB) |", "| --- | --- | --- | --- |"]
        for run in rows:
            if run["status"] == "ok":
                rss = (run.get("peak_rss_bytes") or 0) / 1e6
                lines.append(
                    f"| {run['backend']} | ok | {run['wall_time_seconds']:.3f} | {rss:.0f} |"
                )
            else:
                lines.append(f"| {run['backend']} | {run['status']} | - | - |")
        return "\n".join(lines)

    return "Done."


async def execute_plan(plan: Plan, state: SessionState, console) -> str:
    """Offload the blocking run to a worker thread; keep the UI thread free.

    The Rich status spinner is transient: it animates on the UI thread while the
    native work runs in the executor, then vanishes, leaving the Markdown
    summary behind.
    """

    loop = asyncio.get_running_loop()
    with console.status(f"[bold cyan]{plan.running_message}…", spinner="dots"):
        result = await loop.run_in_executor(None, run_plan, plan)
    return apply_summary(plan, result, state)
