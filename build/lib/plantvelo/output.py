"""Loom layer validation, provenance attributes, and atomic creation."""
from __future__ import annotations

import os
from collections import OrderedDict
from pathlib import Path
from typing import Callable, Mapping, Optional

import numpy as np

from plantvelo.loom_schema import (
    PRIOR_AWARE_SCHEMA_ID,
    VELOCYTO_DEFAULT_SCHEMA_ID,
    LoomSchemaProfile,
)

def prepare_loom_layers(
    matrices: Mapping[str, np.ndarray],
    numeric_dtype: str,
    profile: LoomSchemaProfile,
) -> Mapping[str, np.ndarray]:
    if set(matrices) != set(profile.layers):
        raise ValueError(
            "loom matrices do not match schema '{}': expected {}".format(
                profile.schema_id,
                ", ".join(profile.layers),
            )
        )
    shapes = {np.asarray(matrices[name]).shape for name in profile.layers}
    if len(shapes) != 1:
        raise ValueError("loom layer dimensions are inconsistent")

    layers = OrderedDict()
    total = np.zeros(next(iter(shapes)), dtype=np.float32)
    for name in profile.layers:
        total += np.asarray(matrices[name], dtype=np.float32)
    layers[""] = total
    for name in profile.layers:
        layers[name] = np.asarray(matrices[name]).astype(
            numeric_dtype, order="C", copy=False
        )
    return layers


def build_file_attrs(
    *,
    version: str,
    ir_prior_path: str,
    ir_prior_sha256: str,
    matched_ir_introns: int,
    matched_ir_genes: int,
    gtf_path: str,
    gtf_sha256: str,
    source_kind: str = "custom",
    species: str = "none",
    filter_description: str = "IR_class=high_confidence_IR;min_protocols=1",
) -> Mapping[str, object]:
    return {
        "plantvelo_version": version,
        "classification_schema_version": PRIOR_AWARE_SCHEMA_ID,
        "counting_mode": "prior-aware",
        "classification_precedence": "U>R>S",
        "ir_prior_source": source_kind,
        "species": species,
        "ir_prior_path": (
            str(ir_prior_path)
            if source_kind == "builtin"
            else str(Path(ir_prior_path).resolve())
        ),
        "ir_prior_sha256": ir_prior_sha256,
        "ir_prior_filter": filter_description,
        "matched_ir_introns": matched_ir_introns,
        "matched_ir_genes": matched_ir_genes,
        "gtf_path": str(Path(gtf_path).resolve()),
        "gtf_sha256": gtf_sha256,
    }


def build_velocyto_file_attrs(
    *,
    version: str,
    velocyto_version: str,
    logic_requested: str,
    logic_resolved: str,
    gtf_path: str,
    gtf_sha256: str,
) -> Mapping[str, object]:
    return {
        "plantvelo_version": version,
        "classification_schema_version": VELOCYTO_DEFAULT_SCHEMA_ID,
        "counting_mode": "off",
        "velocyto_version": velocyto_version,
        "velocyto_logic": logic_requested,
        "velocyto_logic_resolved": logic_resolved,
        "gtf_path": str(Path(gtf_path).resolve()),
        "gtf_sha256": gtf_sha256,
    }


def create_loom_atomic(
    output_path: str,
    layers: Mapping[str, np.ndarray],
    row_attrs: Mapping[str, np.ndarray],
    col_attrs: Mapping[str, np.ndarray],
    file_attrs: Mapping[str, object],
    *,
    create: Optional[Callable] = None,
) -> None:
    if create is None:
        import loompy

        create = loompy.create

    temporary_path = output_path + ".tmp"
    if os.path.exists(temporary_path):
        os.unlink(temporary_path)
    try:
        create(
            temporary_path,
            layers,
            row_attrs,
            col_attrs,
            file_attrs=file_attrs,
        )
        os.replace(temporary_path, output_path)
    except Exception:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise
