"""Agentic REPL for the genome assembly toolkit (the `oxidas` shell).

The heavy lifting lives in the core library; this package only adds an
interactive front end: a finite-state loop, a rule-based intent parser, and a
Rich/prompt-toolkit UI that offloads blocking native work off the UI thread.
"""

from .state import Phase, SessionState

__all__ = ["Phase", "SessionState"]
