"""plantvelo.commands
~~~~~~~~~~~~~~~~~~~~
CLI entry point: ``plantvelo``
"""
import click
from importlib import import_module
import logging
import sys
from collections import OrderedDict

import plantvelo._version


class NaturalOrderGroup(click.Group):
    """Preserve command order without importing counting or merge engines."""

    lazy_commands = OrderedDict([
        ("run", "Count molecules with prior-aware or velocyto Default logic"),
        ("merge", "Strictly validate and combine two or more LOOMFILES"),
        ("irprior", "Extract intron-retention evidence from one long-read BAM"),
    ])

    def list_commands(self, ctx):
        return list(self.commands) + [
            name for name in self.lazy_commands if name not in self.commands
        ]

    def get_command(self, ctx, name):
        command = super().get_command(ctx, name)
        if command is not None or name not in self.lazy_commands:
            return command
        module = import_module("plantvelo.commands." + name)
        command = getattr(module, name)
        import_module("plantvelo.commands").__dict__[name] = command
        return command

    def format_commands(self, ctx, formatter):
        rows = []
        for name in self.list_commands(ctx):
            if name in self.commands:
                command = self.commands[name]
                if command.hidden:
                    continue
                help_text = command.get_short_help_str()
            else:
                help_text = self.lazy_commands[name]
            rows.append((name, help_text))
        if rows:
            with formatter.section("Commands"):
                formatter.write_dl(rows)


@click.version_option(version=plantvelo._version.__version__)
@click.group(
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings=dict(max_content_width=300, terminal_width=300),
)
def cli() -> None:
    """PlantVelocity – RNA velocity for plant single-cell data.

    Counts molecules with prior-aware PlantVelo or velocyto Default logic,
    combines outputs that share one classification schema, and extracts
    single-sample intron-retention evidence from long-read BAM files.
    """
    logging.basicConfig(
        stream=sys.stdout,
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG,
    )
