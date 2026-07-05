"""FASTA and FASTQ I/O helpers with no runtime dependency on BioPython."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

from .kmers import normalize_sequence


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class FastqRecord:
    name: str
    sequence: str
    quality: str
    description: str = ""


def _open_text(path: str | Path):
    return Path(path).open("rt", encoding="utf-8")


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read FASTA records from disk."""

    records: list[FastaRecord] = []
    name: str | None = None
    description = ""
    chunks: list[str] = []

    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(FastaRecord(name, normalize_sequence("".join(chunks)), description))
                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                name = parts[0] if parts else ""
                description = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        records.append(FastaRecord(name, normalize_sequence("".join(chunks)), description))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def read_fastq(path: str | Path) -> list[FastqRecord]:
    """Read FASTQ records from disk."""

    records: list[FastqRecord] = []
    with _open_text(path) as handle:
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            line_number += 4

            if not quality:
                raise ValueError(f"Incomplete FASTQ record ending near line {line_number} in {path}")
            if not header.startswith("@"):
                raise ValueError(f"Expected FASTQ header '@' near line {line_number - 3} in {path}")
            if not plus.startswith("+"):
                raise ValueError(f"Expected FASTQ separator '+' near line {line_number - 1} in {path}")

            header_text = header[1:].strip()
            parts = header_text.split(maxsplit=1)
            name = parts[0] if parts else ""
            description = parts[1] if len(parts) > 1 else ""
            seq = normalize_sequence(sequence)
            qual = quality.rstrip("\n\r")
            if len(seq) != len(qual):
                raise ValueError(f"Sequence and quality length mismatch for FASTQ record {name!r}")
            records.append(FastqRecord(name, seq, qual, description))

    if not records:
        raise ValueError(f"No FASTQ records found in {path}")

    return records


def read_sequences(path: str | Path) -> list[str]:
    """Read sequences from FASTA or FASTQ based on the first byte."""

    path = Path(path)
    with _open_text(path) as handle:
        first = handle.read(1)

    if first == ">":
        return [record.sequence for record in read_fasta(path)]
    if first == "@":
        return [record.sequence for record in read_fastq(path)]
    raise ValueError(f"Cannot infer FASTA/FASTQ format for {path}")


def wrap_sequence(sequence: str, width: int = 80) -> Iterator[str]:
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def write_fasta(records: Iterable[FastaRecord | tuple[str, str]], path: str | Path) -> None:
    """Write FASTA records to disk."""

    with Path(path).open("wt", encoding="utf-8") as handle:
        for record in records:
            if isinstance(record, FastaRecord):
                name = record.name
                description = f" {record.description}" if record.description else ""
                sequence = record.sequence
            else:
                name, sequence = record
                description = ""
            handle.write(f">{name}{description}\n")
            for line in wrap_sequence(normalize_sequence(sequence)):
                handle.write(f"{line}\n")


def write_fastq(records: Iterable[FastqRecord | tuple[str, str, str]], path: str | Path) -> None:
    """Write FASTQ records to disk."""

    with Path(path).open("wt", encoding="utf-8") as handle:
        for record in records:
            if isinstance(record, FastqRecord):
                name = record.name
                description = f" {record.description}" if record.description else ""
                sequence = record.sequence
                quality = record.quality
            else:
                name, sequence, quality = record
                description = ""
            sequence = normalize_sequence(sequence)
            if len(sequence) != len(quality):
                raise ValueError(f"Sequence and quality length mismatch for record {name!r}")
            handle.write(f"@{name}{description}\n{sequence}\n+\n{quality}\n")
