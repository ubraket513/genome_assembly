"""Finite-state-machine states and session memory for the agentic REPL."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class Phase(str, Enum):
    """States of the read-eval-print loop.

    IDLE -> PARSING -> CONFIRM -> RUNNING -> SUMMARY -> IDLE, with ERROR as a
    recoverable side exit that returns to IDLE.
    """

    IDLE = "idle"
    PARSING = "parsing"
    CONFIRM = "confirm"
    RUNNING = "running"
    SUMMARY = "summary"
    ERROR = "error"
    EXIT = "exit"


@dataclass
class SessionState:
    """What the session remembers between turns.

    This is what lets a user say "now run stats on that output" without
    repeating file paths: the last assembly's outputs are recorded here.
    """

    phase: Phase = Phase.IDLE
    last_reads: Path | None = None
    last_reference: Path | None = None
    last_outdir: Path | None = None
    last_contigs: Path | None = None
    history: list[str] = field(default_factory=list)

    def note(self, line: str) -> None:
        self.history.append(line)
