"""CLI adapter for the single packaged long-read extraction engine."""
import click


@click.command(
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "help_option_names": [],
    },
    short_help="Extract intron-retention evidence from one long-read BAM",
)
@click.pass_context
def irprior(ctx):
    """Forward arguments to the extractor's shared argparse interface."""
    from plantvelo.irprior import build_arg_parser, parse_and_run

    parser = build_arg_parser()
    parser.prog = ctx.command_path
    parse_and_run(parser, ctx.args)
