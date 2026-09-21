"""CLI wrapper for strict PlantVelo Loom merging."""
import click

from plantvelo.merge import MergeValidationError, merge_loom_files


@click.command()
@click.argument(
    "loomfiles",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, readable=True),
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Combined Loom output path.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Atomically replace an existing output after validation.",
)
def merge(loomfiles, output, force):
    """Strictly validate and combine two or more LOOMFILES."""
    if len(loomfiles) < 2:
        raise click.UsageError("at least two input Loom files are required")
    try:
        summary = merge_loom_files(loomfiles, output, force=force)
    except MergeValidationError as error:
        raise click.ClickException(str(error)) from error
    click.echo(
        "Merged {} genes x {} cells: {}".format(
            summary.gene_count,
            summary.cell_count,
            output,
        )
    )
