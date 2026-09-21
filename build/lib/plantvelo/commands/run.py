"""The public ``plantvelo run`` command."""
from typing import Optional, Tuple

import click

from plantvelo.commands._run import _run
from plantvelo.species import discover_species


@click.command(
    short_help="Count molecules with prior-aware or velocyto Default logic"
)
@click.argument(
    "bamfile",
    nargs=-1,
    required=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
)
@click.argument(
    "gtffile",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "--ir-prior",
    required=False,
    help="CSV/TSV high-confidence intron-retention prior.",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "--species",
    required=False,
    default=None,
    type=click.Choice(discover_species(), case_sensitive=True),
    help="Built-in species IR prior selected from packaged TSV resources.",
)
@click.option(
    "--ir-mode",
    type=click.Choice(("prior-aware", "off"), case_sensitive=True),
    default="prior-aware",
    show_default=True,
    help="Enable prior-aware IR counting or use velocyto Default counting.",
)
@click.option(
    "--bcfile",
    "-b",
    help=(
        "Valid barcodes file, to filter the bam. If omitted, all cell "
        "barcodes are considered. Entries must match the BAM CB tag."
    ),
    default=None,
    type=click.Path(resolve_path=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--outputfolder",
    "-o",
    help="Output folder, created if it does not exist.",
    default=None,
    type=click.Path(exists=False),
)
@click.option(
    "--sample-name",
    help=(
        "Output filename stem and Seurat-compatible cell barcode prefix: "
        "<sample_name>_<barcode>-1."
    ),
    default=None,
    type=str,
)
@click.option(
    "--metadatatable",
    "-s",
    help="CSV table containing sample metadata.",
    default=None,
    type=click.Path(resolve_path=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--mask",
    "-m",
    help="GTF file containing intervals to mask (repeat regions).",
    default=None,
    type=click.Path(resolve_path=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--onefilepercell",
    "-c",
    help="Interpret each BAM file as one independent cell.",
    default=False,
    is_flag=True,
)
@click.option(
    "--without-umi",
    "-U",
    help="Count reads instead of UMI-collapsed molecules.",
    default=False,
    is_flag=True,
)
@click.option(
    "--umi-extension",
    "-u",
    help="Extend UMI by `chr`, `Gene`, or `[N]bp` (default: no).",
    default="no",
)
@click.option(
    "--multimap",
    "-M",
    help="Consider non-unique mappings (not recommended).",
    default=False,
    is_flag=True,
)
@click.option(
    "--samtools-threads",
    "-@",
    help="Number of threads for samtools sort (default: 16).",
    default=16,
)
@click.option(
    "--samtools-memory",
    help="MB of memory per samtools thread (default: 2048).",
    default=2048,
)
@click.option(
    "--dtype",
    "-t",
    help="Numeric dtype for loom layers (default: uint32).",
    default="uint32",
)
@click.option(
    "--dump",
    "-d",
    help="Dump molecular mapping diagnostics every N cells (default: 0).",
    default="0",
)
@click.option(
    "--verbose",
    "-v",
    help="Verbosity: -v warnings, -vv info, -vvv debug.",
    count=True,
    default=1,
)
def run(
    bamfile: Tuple[str, ...],
    gtffile: str,
    ir_prior: Optional[str],
    species: Optional[str],
    ir_mode: str,
    bcfile: str,
    outputfolder: str,
    sample_name: str,
    metadatatable: str,
    mask: str,
    onefilepercell: bool,
    without_umi: bool,
    umi_extension: str,
    multimap: bool,
    samtools_threads: int,
    samtools_memory: int,
    dtype: str,
    dump: str,
    verbose: int,
    additional_ca: dict = {},
) -> None:
    """Count BAMFILE molecules against GTFFILE."""

    if ir_mode == "off" and (species is not None or ir_prior is not None):
        raise click.UsageError(
            "--species and --ir-prior must not be provided when --ir-mode=off"
        )
    if ir_mode == "prior-aware" and species is not None and ir_prior is not None:
        raise click.UsageError("--species and --ir-prior are mutually exclusive")
    if ir_mode == "prior-aware" and species is None and ir_prior is None:
        raise click.UsageError("one of --species or --ir-prior is required")

    return _run(
        bamfile=bamfile,
        gtffile=gtffile,
        ir_mode=ir_mode,
        ir_prior=ir_prior,
        species=species,
        bcfile=bcfile,
        outputfolder=outputfolder,
        sampleid=sample_name,
        sample_name=sample_name,
        metadatatable=metadatatable,
        repmask=mask,
        onefilepercell=onefilepercell,
        without_umi=without_umi,
        umi_extension=umi_extension,
        multimap=multimap,
        test=False,
        samtools_threads=samtools_threads,
        samtools_memory=samtools_memory,
        dump=dump,
        loom_numeric_dtype=dtype,
        verbose=verbose,
        additional_ca=additional_ca,
    )
