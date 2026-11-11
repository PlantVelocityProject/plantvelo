"""Loom classification schemas shared by output and merge validation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Tuple


PRIOR_AWARE_SCHEMA_ID = "prior-aware-v1"
VELOCYTO_DEFAULT_SCHEMA_ID = "velocyto-default-v1"


class LoomSchemaError(ValueError):
    """Raised when Loom schema metadata is missing or unsupported."""


@dataclass(frozen=True)
class LoomSchemaProfile:
    schema_id: str
    counting_mode: str
    layers: Tuple[str, ...]
    required_semantic_attrs: Tuple[str, ...]


COMMON_SEMANTIC_ATTRS = (
    "plantvelo_version",
    "classification_schema_version",
    "counting_mode",
    "gtf_sha256",
)

PRIOR_AWARE_PROFILE = LoomSchemaProfile(
    schema_id=PRIOR_AWARE_SCHEMA_ID,
    counting_mode="prior-aware",
    layers=("spliced", "unspliced", "retained", "ambiguous"),
    required_semantic_attrs=COMMON_SEMANTIC_ATTRS
    + (
        "classification_precedence",
        "ir_prior_sha256",
        "ir_prior_filter",
        "matched_ir_introns",
        "matched_ir_genes",
    ),
)

VELOCYTO_DEFAULT_PROFILE = LoomSchemaProfile(
    schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
    counting_mode="off",
    layers=("spliced", "unspliced", "ambiguous"),
    required_semantic_attrs=COMMON_SEMANTIC_ATTRS
    + (
        "velocyto_version",
        "velocyto_logic",
        "velocyto_logic_resolved",
    ),
)

SCHEMA_PROFILES = {
    PRIOR_AWARE_PROFILE.schema_id: PRIOR_AWARE_PROFILE,
    VELOCYTO_DEFAULT_PROFILE.schema_id: VELOCYTO_DEFAULT_PROFILE,
}


def get_loom_schema(schema_id: object) -> LoomSchemaProfile:
    key = str(schema_id)
    try:
        return SCHEMA_PROFILES[key]
    except KeyError as error:
        raise LoomSchemaError(
            "unknown classification schema: {}".format(key)
        ) from error


def normalize_semantic_attrs(
    profile: LoomSchemaProfile,
    attrs: Mapping[str, object],
) -> Mapping[str, object]:
    normalized = dict(attrs)
    if (
        profile.schema_id == PRIOR_AWARE_SCHEMA_ID
        and "counting_mode" not in normalized
    ):
        normalized["counting_mode"] = "prior-aware"

    for name in profile.required_semantic_attrs:
        if name not in normalized:
            raise LoomSchemaError(
                "missing semantic attribute '{}' for schema '{}'".format(
                    name,
                    profile.schema_id,
                )
            )
    if normalized["classification_schema_version"] != profile.schema_id:
        raise LoomSchemaError("classification schema metadata is inconsistent")
    if normalized["counting_mode"] != profile.counting_mode:
        raise LoomSchemaError("counting mode metadata is inconsistent")
    return normalized
