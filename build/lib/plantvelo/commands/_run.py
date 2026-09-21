"""Core prior-aware pipeline, based on velocyto's BAM counting engine."""
import sys
import os
import gzip
import numpy as np
import subprocess
import multiprocessing
import itertools
import logging
from typing import Any, Dict, Optional, Tuple

import velocyto as vcy

from plantvelo.ir_prior import sha256_file
from plantvelo.species import resolve_prior_source
from plantvelo.output import (
    build_file_attrs,
    build_velocyto_file_attrs,
    create_loom_atomic,
    prepare_loom_layers,
)
from plantvelo.pipeline import load_annotations_for_mode
from plantvelo.qc import (
    write_classification_qc,
    write_prior_qc,
    write_unmatched_prior,
)
from plantvelo.run_mode import resolve_run_mode, validate_logic_contract


def _run(
    *,
    bamfile: Tuple[str, ...],
    gtffile: str,
    ir_mode: str,
    ir_prior: Optional[str],
    species: Optional[str] = None,
    bcfile: str,
    outputfolder: str,
    sampleid: str,
    sample_name: Optional[str] = None,
    metadatatable: str,
    repmask: str,
    onefilepercell: bool,
    without_umi: bool,
    umi_extension: str,
    multimap: bool,
    test: bool,
    samtools_threads: int,
    samtools_memory: int,
    loom_numeric_dtype: str,
    dump: bool,
    verbose: int,
    additional_ca: dict = {},
) -> None:
    """Run prior-aware molecule counting and write loom plus QC files.

    Parameters
    ----------
    sample_name : str, optional
        When provided, used as the output loom filename and as prefix for
        CellIDs in Seurat-compatible format (<sample_name>_<barcode>-1).
        When None, CellIDs use the legacy format (<sampleid>:<barcode><gem_grp>)
        and the filename is derived from the BAM file.
    """

    ########################
    #    Resolve Inputs    #
    ########################

    logging.basicConfig(
        stream=sys.stdout,
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=[logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG][verbose],
    )

    ir_registry = set()
    mode_profile = resolve_run_mode(
        ir_mode,
        ir_prior,
        ir_registry,
        species=species,
    )
    prior_source = None
    if mode_profile.requires_ir_prior:
        prior_source = resolve_prior_source(species, ir_prior)

    if isinstance(bamfile, tuple) and len(bamfile) > 1 and bamfile[-1][-4:] in [".bam", ".sam"]:
        multi = True
    elif isinstance(bamfile, tuple) and len(bamfile) == 1:
        multi = False
    else:
        raise IOError(f"Something went wrong in the argument parsing. You passed as bamfile: {bamfile}")

    if onefilepercell and multi:
        if bcfile is not None:
            raise ValueError("Inputs incompatibility. --bcfile/-b option was used together with --onefilepercell/-c option.")
        logging.warning("Each bam file will be interpreted as a DIFFERENT cell")
    elif not onefilepercell and multi:
        logging.warning("Several input files but --onefilepercell is False. Each bam file will be interpreted as containing a SET of cells!!!")

    if sampleid is None:
        assert metadatatable is None, "--metadatatable was specified but cannot fetch sample metadata without valid sampleid"
        if multi:
            logging.warning("When using multiple files you may want to use --sample-name option to specify the name of the output file")
        if multi and not onefilepercell:
            full_name = "_".join([os.path.basename(bamfile[i]).split(".")[0] for i in range(len(bamfile))])
            if len(full_name) > 50:
                sampleid = f'multi_input_{os.path.basename(bamfile[0]).split(".")[0]}'
            else:
                sampleid = f'multi_input_{full_name}'
        elif multi and onefilepercell:
            sampleid = f'onefilepercell_{os.path.basename(bamfile[0]).split(".")[0]}'
        else:
            sampleid = os.path.basename(bamfile[0]).split(".")[0]
        logging.info(f"No SAMPLEID specified, the sample will be called {sampleid}")

    if outputfolder is None:
        outputfolder = os.path.join(os.path.split(bamfile[0])[0], "plantvelo")
        logging.info(f"No OUTPUTFOLDER specified, find output files inside {outputfolder}")
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    ########################
    #     Barcodes         #
    ########################

    if bcfile is None:
        logging.debug("Cell barcodes will be determined while reading the .bam file")
        valid_bcset = None
    else:
        valid_bcs_list = (gzip.open(bcfile).read().decode() if bcfile.endswith(".gz") else open(bcfile).read()).rstrip().split()
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])
        if len(set(bc.split('-')[0] for bc in valid_bcs_list)) == 1:
            gem_grp = f"-{valid_bcs_list[0].split('-')[-1]}"
        else:
            gem_grp = "x"
        valid_bcset = set(bc.split('-')[0] for bc in valid_bcs_list)
        logging.info(f"Read {len(valid_bcs_list)} cell barcodes from {bcfile}")
        logging.debug(f"Example barcode: {valid_bcs_list[0].split('-')[0]} cell_id: {valid_cellid_list[0]}")

    if metadatatable:
        try:
            sample_metadata = vcy.MetadataCollection(metadatatable)
            sample = sample_metadata.where("SampleID", sampleid)
            if len(sample) == 0:
                logging.error(f"Sample ID {sampleid} not found in sample sheet")
                sample = {}
            elif len(sample) > 1:
                logging.error(f"Sample ID {sampleid} has multiple lines in sample sheet")
                sys.exit(1)
            else:
                sample = sample[0].dict
        except (NameError, TypeError):
            logging.warning("SAMPLEFILE was not specified. add -s SAMPLEFILE to add metadata.")
            sample = {}
    else:
        sample = {}

    ########################
    #     Start Analysis   #
    ########################

    if without_umi:
        if umi_extension != "no":
            logging.warning("--umi-extension was specified but incompatible with --without-umi, it will be ignored!")
        umi_extension = "without_umi"

    exincounter = vcy.ExInCounter(
        sampleid=sampleid,
        logic=mode_profile.logic_factory,
        valid_bcset=valid_bcset,
        umi_extension=umi_extension,
        onefilepercell=onefilepercell,
        dump_option=dump,
        outputfolder=outputfolder,
    )
    logic_obj = exincounter.logic
    validate_logic_contract(mode_profile, logic_obj)

    if mode_profile.requires_ir_prior:
        logging.info(
            "Load annotation from %s and IR prior from %s",
            gtffile,
            prior_source.logical_path,
        )
    else:
        logging.info("Load annotation from %s without an IR prior", gtffile)
    annotations_by_chrm_strand, prior_match = load_annotations_for_mode(
        exincounter,
        mode_profile,
        prior_source,
        gtffile,
        ir_registry,
    )
    transcript_models = list(
        itertools.chain.from_iterable(
            transcripts.values()
            for transcripts in annotations_by_chrm_strand.values()
        )
    )
    feature_count = sum(1 for model in transcript_models for _ in model)
    logging.debug(
        "Generated %d features from %d transcript models",
        feature_count,
        len(transcript_models),
    )
    if prior_match is not None and prior_match.unmatched:
        logging.warning(
            "%d high-confidence IR prior records did not match the GTF; "
            "see the unmatched QC file",
            prior_match.unmatched_unique_introns,
        )

    try:
        mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.warning("Cannot determine available memory; assuming 32 GB")
        mb_available = 32000

    threads_to_use = min(samtools_threads, multiprocessing.cpu_count())
    mb_to_use = int(min(samtools_memory, mb_available / (len(bamfile) * threads_to_use)))
    compression = vcy.BAM_COMPRESSION

    if onefilepercell and without_umi:
        tagname = "NOTAG"
    elif onefilepercell:
        tagname = "NOTAG"
        exincounter.peek_umi_only(bamfile[0])
    else:
        exincounter.peek(bamfile[0])
        tagname = exincounter.cellbarcode_str

    if multi and onefilepercell:
        bamfile_cellsorted = list(bamfile)
    elif onefilepercell:
        bamfile_cellsorted = [bamfile[0]]
    else:
        bamfile_cellsorted = [
            f"{os.path.join(os.path.dirname(bmf), 'cellsorted_' + os.path.basename(bmf))}"
            for bmf in bamfile
        ]

    sorting_process: Dict[int, Any] = {}
    check_end_process = False
    for ni, bmf_cellsorted in enumerate(bamfile_cellsorted):
        command = f"samtools sort -l {compression} -m {mb_to_use}M -t {tagname} -O BAM -@ {threads_to_use} -o {bmf_cellsorted} {bamfile[ni]}"
        if os.path.exists(bmf_cellsorted):
            logging.warning(f"The file {bmf_cellsorted} already exists. Sorting step will be skipped.")
        else:
            sorting_process[ni] = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
            logging.info(f"Sorting {bamfile[ni]} → {bmf_cellsorted}")
            check_end_process = True

    if repmask is not None:
        logging.info(f"Load repeat mask from {repmask}")
        exincounter.read_repeats(repmask)

    exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)

    if check_end_process:
        logging.info("Waiting for bam sorting to finish…")
        for k in sorting_process.keys():
            returncode = sorting_process[k].wait()
            if returncode == 0:
                logging.info(f"bam file #{k} sorted")
            else:
                raise MemoryError(
                    f"bam file #{k} could not be sorted by cells. "
                    "Try upgrading samtools (>= 1.6) or increase --samtools-memory."
                )

    logging.debug("Start molecule counting!")
    results = exincounter.count(bamfile_cellsorted, multimap=multimap)
    dict_list_arrays, cell_bcs_order = results

    ########################
    #         Output       #
    ########################

    if not exincounter.filter_mode:
        valid_bcset = exincounter.valid_bcset
        valid_bcs_list = list(valid_bcset)
        gem_grp = ""
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])

    if sample_name is not None:
        # Reformat barcodes for Seurat compatibility:
        # {sampleid}:{barcode}{gem_grp}  →  {sample_name}_{barcode}-1
        # gem_grp "x" (mixed gem groups) is normalised to "-1"; a numeric
        # suffix such as "-1" is kept as-is; empty suffix gets "-1".
        _gem_suffix = gem_grp if (gem_grp and gem_grp != "x") else "-1"
        cell_ids = np.array([f"{sample_name}_{v_bc}{_gem_suffix}" for v_bc in cell_bcs_order])
    else:
        cell_ids = np.array([f"{sampleid}:{v_bc}{gem_grp}" for v_bc in cell_bcs_order])
    ca = {"CellID": cell_ids}
    ca.update(additional_ca)
    for key, value in sample.items():
        ca[key] = np.full(len(cell_bcs_order), value)

    loom_name = sample_name if sample_name is not None else sampleid
    outfile = os.path.join(outputfolder, f"{loom_name}.loom")
    logging.debug(f"Generating output file {outfile}")

    atr_table = (
        ("Gene",       "genename", str),
        ("Accession",  "geneid",   str),
        ("Chromosome", "chrom",    str),
        ("Strand",     "strand",   str),
        ("Start",      "start",    int),
        ("End",        "end",      int),
    )

    ra: Dict[str, np.ndarray] = {}
    for name_col_attr, name_obj_attr, dtyp in atr_table:
        tmp_array = np.zeros((len(exincounter.genes),), dtype=object)
        for gene_id, gene_info in exincounter.genes.items():
            tmp_array[exincounter.geneid2ix[gene_id]] = getattr(gene_info, name_obj_attr)
        ra[name_col_attr] = tmp_array.astype(dtyp)

    layers: Dict[str, np.ndarray] = {}
    for layer_name in logic_obj.layers:
        layers[layer_name] = np.concatenate(dict_list_arrays[layer_name], axis=1)
        del dict_list_arrays[layer_name]

    import plantvelo
    loom_layers = prepare_loom_layers(
        layers,
        loom_numeric_dtype,
        mode_profile.schema,
    )
    if mode_profile.name == "prior-aware":
        logic_obj.classification_qc.validate_conservation()
        file_attrs = build_file_attrs(
            version=plantvelo.__version__,
            ir_prior_path=prior_match.prior.path,
            ir_prior_sha256=prior_match.prior.sha256,
            matched_ir_introns=prior_match.matched_unique_introns,
            matched_ir_genes=prior_match.matched_genes,
            gtf_path=gtffile,
            gtf_sha256=sha256_file(gtffile),
            source_kind=prior_match.prior.source_kind,
            species=prior_match.prior.species or "none",
            filter_description=prior_match.prior.filter_description,
        )
    else:
        file_attrs = build_velocyto_file_attrs(
            version=plantvelo.__version__,
            velocyto_version=vcy.__version__,
            logic_requested="Default",
            logic_resolved=getattr(
                logic_obj,
                "name",
                type(logic_obj).__name__,
            ),
            gtf_path=gtffile,
            gtf_sha256=sha256_file(gtffile),
        )
    create_loom_atomic(
        outfile,
        loom_layers,
        ra,
        ca,
        file_attrs,
    )

    if mode_profile.writes_prior_qc:
        qc_prefix = os.path.join(outputfolder, loom_name)
        write_prior_qc(f"{qc_prefix}.ir_prior_qc.tsv", prior_match, gtffile)
        write_unmatched_prior(
            f"{qc_prefix}.ir_prior_unmatched.tsv",
            prior_match,
        )
        write_classification_qc(
            f"{qc_prefix}.classification_qc.tsv",
            logic_obj.classification_qc,
        )
    logging.info(f"Loom file written: {outfile}")
    logging.debug("Terminated successfully!")
