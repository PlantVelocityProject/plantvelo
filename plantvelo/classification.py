"""Pure molecule-level evidence and classification rules."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Set


@dataclass(frozen=True, order=True)
class IntronKey:
    """A genomic intron in normalized 0-based, half-open coordinates."""

    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str


@dataclass
class MoleculeEvidence:
    """All evidence collected for one cell-barcode, UMI, and gene."""

    ordinary_unspliced: bool = False
    hc_retained: Set[IntronKey] = field(default_factory=set)
    hc_spliced: Set[IntronKey] = field(default_factory=set)
    exonic: bool = False
    model_conflict: bool = False
    multigene: bool = False
    invalid_alignment: bool = False


class Classification(str, Enum):
    SPLICED = "spliced"
    UNSPLICED = "unspliced"
    RETAINED = "retained"
    AMBIGUOUS = "ambiguous"
    DISCARDED_MULTIGENE = "discarded_multigene"
    DISCARDED_UNANNOTATED = "discarded_unannotated"
    DISCARDED_INVALID_ALIGNMENT = "discarded_invalid_alignment"


def classify_molecule(evidence: MoleculeEvidence) -> Classification:
    """Assign one mutually exclusive state to a molecule's merged evidence."""

    if evidence.multigene:
        return Classification.DISCARDED_MULTIGENE
    if evidence.invalid_alignment:
        return Classification.DISCARDED_INVALID_ALIGNMENT
    if evidence.model_conflict:
        return Classification.AMBIGUOUS
    if evidence.hc_retained & evidence.hc_spliced:
        return Classification.AMBIGUOUS
    if evidence.ordinary_unspliced:
        return Classification.UNSPLICED
    if evidence.hc_retained:
        return Classification.RETAINED
    if evidence.exonic or evidence.hc_spliced:
        return Classification.SPLICED
    return Classification.AMBIGUOUS
