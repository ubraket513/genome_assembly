"""Rule-based intent parsing and plan building for the agentic REPL.

Honest scope: this is a deterministic keyword/regex parser, not a language
model. It maps a line of input to one of the toolkit's operations and fills in
parameters from the text and from session memory. Everything here is pure and
unit-tested; the REPL and the executor build on top of it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

from .state import SessionState

_REFERENCE_SUFFIXES = (".fna", ".fa", ".fasta", ".fasta.gz", ".fna.gz")
_READS_SUFFIXES = (".fastq", ".fq", ".fasta", ".fastq.gz", ".fq.gz")
_BACKENDS = ("native", "cython", "python")


@dataclass
class Intent:
    """A recognized operation plus the raw parameters pulled from the text."""

    kind: str
    params: dict[str, object] = field(default_factory=dict)


@dataclass
class Plan:
    """A ready-to-run, resolved operation shown to the user before execution."""

    kind: str
    params: dict[str, object]
    description: str
    running_message: str
    needs_confirm: bool = True
    error: str | None = None


def _first_path(text: str, suffixes: tuple[str, ...]) -> Path | None:
    for token in re.split(r"\s+", text):
        cleaned = token.strip("\"'`,")
        if cleaned.lower().endswith(suffixes):
            return Path(cleaned)
    return None


def _number_after(text: str, keywords: tuple[str, ...]) -> float | None:
    """Find a number tied to a keyword, in either order ("threads 4" or "4 threads")."""

    for keyword in keywords:
        after = re.search(rf"{keyword}\s*[=:]?\s*(\d+(?:\.\d+)?)", text)
        if after:
            return float(after.group(1))
        before = re.search(rf"(\d+(?:\.\d+)?)\s*{keyword}\b", text)
        if before:
            return float(before.group(1))
    return None


def _coverage(text: str) -> float | None:
    inline = re.search(r"(\d+(?:\.\d+)?)\s*x\b", text)
    if inline:
        return float(inline.group(1))
    return _number_after(text, ("coverage", "depth", "cov"))


def _references_self(text: str) -> bool:
    return bool(re.search(r"\b(that|it|those|previous|last)\b", text))


def parse_intent(text: str) -> Intent | None:
    """Classify a line of input. Returns None when nothing is recognized."""

    lowered = text.strip().lower()
    if not lowered:
        return None
    if lowered in {"exit", "quit", "bye", ":q"} or lowered.startswith("exit"):
        return Intent("exit")
    if lowered in {"help", "?", "commands"} or lowered.startswith("help"):
        return Intent("help")

    backend = next((b for b in _BACKENDS if re.search(rf"\b{b}\b", lowered)), None)
    threads = _number_after(lowered, ("threads", "cores", "thread"))
    genome_size = _number_after(lowered, ("genome", "size", "bp", "bases"))

    if re.search(r"\b(benchmark|compare|bench|speed)\b", lowered):
        return Intent(
            "benchmark",
            {
                "reference": _first_path(text, _REFERENCE_SUFFIXES),
                "genome_size": int(genome_size) if genome_size else None,
                "backend": backend,
                "coverage": _coverage(lowered),
                "threads": int(threads) if threads else None,
                "self_ref": _references_self(lowered),
            },
        )
    if re.search(r"\b(simulate|generate|make)\b", lowered) and "read" in lowered:
        return Intent(
            "simulate",
            {
                "reference": _first_path(text, _REFERENCE_SUFFIXES),
                "coverage": _coverage(lowered),
                "self_ref": _references_self(lowered),
            },
        )
    if re.search(r"\b(stats|statistics|evaluate|quality|score|assess)\b", lowered):
        return Intent(
            "stats",
            {
                "contigs": _first_path(text, (".fasta", ".fa", ".fna")),
                "reference": _first_path(text, _REFERENCE_SUFFIXES),
                "self_ref": _references_self(lowered),
            },
        )
    if re.search(r"\b(assemble|assembly|build|contigs?)\b", lowered):
        return Intent(
            "assemble",
            {
                "reads": _first_path(text, _READS_SUFFIXES),
                "k": _number_after(lowered, ("k", "kmer", "k-mer")),
                "backend": backend,
                "threads": int(threads) if threads else None,
                "tip_length": _number_after(lowered, ("tip", "tip-length")),
                "bubble_length": _number_after(lowered, ("bubble", "bubble-length")),
                "self_ref": _references_self(lowered),
            },
        )
    return None


def build_plan(intent: Intent, state: SessionState) -> Plan:
    """Resolve an intent against session memory into a runnable, described plan."""

    if intent.kind in {"exit", "help"}:
        return Plan(intent.kind, {}, "", "", needs_confirm=False)

    if intent.kind == "simulate":
        reference = intent.params.get("reference") or state.last_reference
        if reference is None:
            return _unresolved("simulate", "which reference FASTA should I simulate reads from?")
        coverage = intent.params.get("coverage") or 30.0
        output = Path("reads.fastq")
        params = {"reference": Path(reference), "coverage": float(coverage), "output": output}
        return Plan(
            "simulate",
            params,
            f"Simulate reads from **{reference}** at **{coverage:g}x** coverage into `{output}`.",
            "Simulating reads",
        )

    if intent.kind == "assemble":
        reads = intent.params.get("reads") or state.last_reads
        if reads is None:
            return _unresolved("assemble", "which reads file should I assemble?")
        k = int(intent.params.get("k") or 31)
        backend = intent.params.get("backend") or "python"
        threads = int(intent.params.get("threads") or 1)
        tip = int(intent.params.get("tip_length") or 0)
        bubble = int(intent.params.get("bubble_length") or 0)
        outdir = Path("assembly_out")
        params = {
            "reads": Path(reads),
            "k": k,
            "backend": backend,
            "threads": threads,
            "tip_length": tip,
            "bubble_length": bubble,
            "outdir": outdir,
        }
        extras = f", backend **{backend}**, threads **{threads}**"
        if tip or bubble:
            extras += f", cleaning (tip {tip}, bubble {bubble})"
        return Plan(
            "assemble",
            params,
            f"Assemble **{reads}** with k=**{k}**{extras} into `{outdir}/`.",
            "Assembling contigs",
        )

    if intent.kind == "stats":
        contigs = intent.params.get("contigs")
        if contigs is None and (intent.params.get("self_ref") or state.last_contigs):
            contigs = state.last_contigs
        if contigs is None:
            return _unresolved("stats", "which contigs FASTA should I evaluate?")
        reference = intent.params.get("reference") or state.last_reference
        params = {"contigs": Path(contigs), "reference": Path(reference) if reference else None}
        ref_text = f" against **{reference}**" if reference else ""
        return Plan(
            "stats",
            params,
            f"Compute assembly statistics for **{contigs}**{ref_text}.",
            "Scoring assembly",
        )

    if intent.kind == "benchmark":
        reference = intent.params.get("reference") or state.last_reference
        genome_size = intent.params.get("genome_size")
        if reference is None and genome_size is None:
            genome_size = 100_000  # sensible default synthetic tier
        backends = intent.params.get("backend")
        backend_list = [backends] if backends else ["python", "cython", "native"]
        coverage = float(intent.params.get("coverage") or 5.0)
        threads = int(intent.params.get("threads") or 1)
        params = {
            "reference": Path(reference) if reference else None,
            "genome_size": genome_size,
            "backends": backend_list,
            "coverage": coverage,
            "threads": threads,
        }
        source = f"**{reference}**" if reference else f"a synthetic **{genome_size} bp** genome"
        return Plan(
            "benchmark",
            params,
            f"Benchmark {', '.join(backend_list)} on {source} "
            f"at {coverage:g}x coverage, threads **{threads}**.",
            "Benchmarking backends",
        )

    return _unresolved(intent.kind, "I could not turn that into an operation.")


def _unresolved(kind: str, question: str) -> Plan:
    return Plan(kind, {}, "", "", needs_confirm=False, error=question)
