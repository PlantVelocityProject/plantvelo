"""Small orchestration boundaries shared by the command pipeline."""
from __future__ import annotations

from typing import MutableSet, Tuple

from plantvelo.classification import IntronKey
from plantvelo.ir_prior import PriorMatchResult, load_ir_prior, match_ir_prior
from plantvelo.species import PriorSource


def bind_ir_prior(
    exincounter,
    prior_source,
    gtf_path: str,
    registry: MutableSet[IntronKey],
) -> Tuple[object, PriorMatchResult]:
    """Load annotations and bind the required prior to a shared registry."""

    if registry:
        raise ValueError("IR registry must be empty before annotation binding")
    if isinstance(prior_source, str):
        prior_source = PriorSource(
            path=prior_source,
            logical_path=prior_source,
            source_kind="custom",
            species=None,
            loader=load_ir_prior,
        )
    prior = prior_source.loader(prior_source.path, prior_source.species) \
        if prior_source.source_kind == "builtin" else prior_source.loader(prior_source.path)
    if prior_source.logical_path != prior.path:
        from dataclasses import replace

        prior = replace(prior, path=prior_source.logical_path)
    annotations = exincounter.read_transcriptmodels(gtf_path)
    match = match_ir_prior(prior, gtf_path, annotations)
    registry.update(match.runtime_introns)
    return annotations, match


def load_annotations_for_mode(
    exincounter,
    profile,
    prior_source,
    gtf_path: str,
    registry: MutableSet[IntronKey],
):
    """Load GTF annotations and bind an IR prior only when required."""

    if profile.requires_ir_prior:
        return bind_ir_prior(
            exincounter,
            prior_source,
            gtf_path,
            registry,
        )
    if registry:
        raise ValueError("IR registry must stay empty when --ir-mode=off")
    return exincounter.read_transcriptmodels(gtf_path), None
