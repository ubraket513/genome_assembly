"""The `oxidas` interactive shell: an async REPL over the toolkit.

Design choice that keeps the terminal sane: prompt-toolkit owns the *input*
phase and Rich owns the *output/progress* phase — they never drive the screen at
the same time. A single asyncio loop is the driver; blocking native work is
offloaded to a thread pool by `runner.execute_plan`.
"""

from __future__ import annotations

import asyncio

from . import __doc__ as _pkg_doc  # noqa: F401  (keeps package import side effects explicit)
from .intents import build_plan, parse_intent
from .runner import execute_plan
from .state import Phase, SessionState

_MISSING_DEPS_HINT = (
    "The oxidas shell needs extra packages. Install them with:\n"
    "    python -m pip install 'genome-assembly[agent]'"
)

_COMMANDS = [
    "simulate",
    "assemble",
    "stats",
    "benchmark",
    "help",
    "exit",
    "native",
    "cython",
    "python",
    "coverage",
    "threads",
]

_HELP = """\
# oxidas — interactive genome assembly

Type natural language or shorthand. Examples:

- `simulate reads from genome.fna at 30x`
- `assemble reads.fastq with k 31 using native backend on 4 threads`
- `now run stats on that output`
- `benchmark native on a 500000 bp genome`

Each step is shown as a plan and asks for confirmation before running.
Type `help` for this message, `exit` to leave.
"""


def _load_ui():
    """Import the optional UI stack, or return None with an actionable hint."""

    try:
        from prompt_toolkit import PromptSession
        from prompt_toolkit.completion import WordCompleter
        from prompt_toolkit.history import InMemoryHistory
        from rich.console import Console
        from rich.markdown import Markdown
        from rich.panel import Panel
    except ImportError:
        return None
    return {
        "PromptSession": PromptSession,
        "WordCompleter": WordCompleter,
        "InMemoryHistory": InMemoryHistory,
        "Console": Console,
        "Markdown": Markdown,
        "Panel": Panel,
    }


async def _confirm(session, prompt_html) -> bool:
    answer = (await session.prompt_async(prompt_html)).strip().lower()
    return answer in {"", "y", "yes"}


async def run_repl(ui: dict, state: SessionState | None = None) -> None:
    console = ui["Console"]()
    Markdown = ui["Markdown"]
    Panel = ui["Panel"]
    session = ui["PromptSession"](
        history=ui["InMemoryHistory"](),
        completer=ui["WordCompleter"](_COMMANDS, ignore_case=True),
    )
    state = state or SessionState()

    console.print(
        Panel.fit(
            "[bold]oxidas[/bold] — agentic genome assembly shell\n"
            "Type a request in plain English, or [bold]help[/bold] / [bold]exit[/bold].",
            border_style="cyan",
        )
    )

    from prompt_toolkit.formatted_text import HTML

    while True:
        state.phase = Phase.IDLE
        try:
            text = await session.prompt_async(HTML("<ansicyan><b>oxidas</b></ansicyan> › "))
        except (EOFError, KeyboardInterrupt):
            break

        state.note(text)
        state.phase = Phase.PARSING
        intent = parse_intent(text)
        if intent is None:
            console.print("[yellow]Not sure what to run.[/yellow] Type [bold]help[/bold] for examples.")
            continue
        if intent.kind == "exit":
            break
        if intent.kind == "help":
            console.print(Markdown(_HELP))
            continue

        plan = build_plan(intent, state)
        if plan.error:
            state.phase = Phase.ERROR
            console.print(f"[yellow]{plan.error}[/yellow]")
            continue

        state.phase = Phase.CONFIRM
        console.print(Markdown(f"**Plan:** {plan.description}"))
        if plan.needs_confirm and not await _confirm(session, HTML("Run this? [<b>Y</b>/n] ")):
            console.print("[dim]Skipped.[/dim]")
            continue

        state.phase = Phase.RUNNING
        try:
            summary = await execute_plan(plan, state, console)
        except Exception as exc:  # noqa: BLE001 - surface any failure without crashing the shell
            state.phase = Phase.ERROR
            console.print(f"[red]Error:[/red] {exc}")
            continue

        state.phase = Phase.SUMMARY
        console.print(Markdown(summary))

    console.print("[dim]bye[/dim]")


def main() -> None:
    """Entry point for the `oxidas` script."""

    ui = _load_ui()
    if ui is None:
        print(_MISSING_DEPS_HINT)
        raise SystemExit(1)
    try:
        asyncio.run(run_repl(ui))
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
