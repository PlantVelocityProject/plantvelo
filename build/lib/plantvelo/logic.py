"""Prior-aware molecule classification adapter for velocyto."""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, Mapping, Optional, Set, Tuple

from velocyto.logic import Logic

from plantvelo.classification import (
    Classification,
    IntronKey,
    MoleculeEvidence,
    classify_molecule,
)
from plantvelo.qc import ClassificationQC


def _feature_intron_key(feature) -> IntronKey:
    transcript = feature.transcript_model
    chromstrand = transcript.chromstrand
    if len(chromstrand) < 2 or chromstrand[-1] not in ("+", "-"):
        raise ValueError("invalid transcript chromstrand")
    return IntronKey(
        transcript.geneid,
        chromstrand[:-1],
        feature.start - 1,
        feature.end,
        chromstrand[-1],
    )


def _segment_tuple(segment) -> Tuple[int, int]:
    if not isinstance(segment, (tuple, list)) or len(segment) != 2:
        raise ValueError("alignment segment must contain start and end")
    start, end = segment
    if not isinstance(start, int) or not isinstance(end, int) or start > end:
        raise ValueError("alignment segment coordinates are invalid")
    return start, end


def _has_exact_splice(
    segments: Set[Tuple[int, int]], intron: IntronKey
) -> bool:
    for left in segments:
        if left[1] != intron.start:
            continue
        for right in segments:
            if right[0] == intron.end + 1 and left[1] < right[0]:
                return True
    return False


def extract_molecule_evidence(
    molitem, ir_introns: Set[IntronKey]
) -> Tuple[MoleculeEvidence, Optional[str]]:
    """Translate one velocyto Molitem into merged gene-level evidence."""

    evidence = MoleculeEvidence()
    mappings = getattr(molitem, "mappings_record", None)
    if not mappings:
        return evidence, None

    gene_ids = {getattr(transcript, "geneid", None) for transcript in mappings}
    if None in gene_ids:
        evidence.invalid_alignment = True
        return evidence, None
    if len(gene_ids) != 1:
        evidence.multigene = True
        return evidence, None
    gene_id = next(iter(gene_ids))

    model_signatures_by_segment = defaultdict(list)

    for transcript, segment_matches in mappings.items():
        matches_by_segment = defaultdict(list)
        spliced_segments: Set[Tuple[int, int]] = set()

        for segment_match in segment_matches:
            try:
                segment = _segment_tuple(segment_match.segment)
                feature = segment_match.feature
                kind = feature.kind
                if kind not in (ord("e"), ord("i")):
                    raise ValueError("unsupported feature type")
            except (AttributeError, TypeError, ValueError):
                evidence.invalid_alignment = True
                continue
            matches_by_segment[segment].append(segment_match)
            if bool(getattr(segment_match, "is_spliced", False)):
                spliced_segments.add(segment)

        for segment, matches in matches_by_segment.items():
            interpretation = set()
            intron_matches = [
                match for match in matches if match.feature.kind == ord("i")
            ]
            if intron_matches:
                has_ordinary = False
                has_retained = False
                for match in intron_matches:
                    try:
                        key = _feature_intron_key(match.feature)
                    except (AttributeError, TypeError, ValueError):
                        evidence.invalid_alignment = True
                        continue
                    if key in ir_introns:
                        evidence.hc_retained.add(key)
                        has_retained = True
                    else:
                        evidence.ordinary_unspliced = True
                        has_ordinary = True
                if has_ordinary:
                    interpretation.add("ordinary_unspliced")
                if has_retained:
                    interpretation.add("hc_retained")
            elif any(match.feature.kind == ord("e") for match in matches):
                evidence.exonic = True
                interpretation.add("exonic")
            else:
                evidence.invalid_alignment = True
            if interpretation:
                model_signatures_by_segment[segment].append(
                    frozenset(interpretation)
                )

        try:
            transcript_features = tuple(transcript)
        except TypeError:
            evidence.invalid_alignment = True
            continue
        for feature in transcript_features:
            if getattr(feature, "kind", None) != ord("i"):
                continue
            try:
                key = _feature_intron_key(feature)
            except (AttributeError, TypeError, ValueError):
                evidence.invalid_alignment = True
                continue
            if key in ir_introns and _has_exact_splice(spliced_segments, key):
                evidence.hc_spliced.add(key)

    evidence.model_conflict = any(
        len(set(signatures)) > 1
        for signatures in model_signatures_by_segment.values()
    )
    return evidence, gene_id


class PriorAwareLogic(Logic):
    """The only PlantVelo 0.2 classification logic."""

    def __init__(self, ir_introns: Set[IntronKey]) -> None:
        self.name = "PriorAwareLogic"
        # Keep the shared set by reference so the pipeline can bind it after
        # velocyto has loaded transcript models.
        self.ir_introns = ir_introns
        self.classification_qc = ClassificationQC()

    @property
    def layers(self):
        return ["spliced", "unspliced", "retained", "ambiguous"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return False

    @property
    def accept_discordant(self) -> bool:
        return False

    @staticmethod
    def _discard_code(outcome: Classification) -> int:
        if outcome == Classification.DISCARDED_MULTIGENE:
            return 1
        if outcome == Classification.DISCARDED_UNANNOTATED:
            return 2
        return 4

    def count(
        self,
        molitem,
        cell_bcidx: int,
        dict_layers_columns: Mapping[str, object],
        geneid2ix: Dict[str, int],
    ) -> int:
        mappings = getattr(molitem, "mappings_record", None)
        if not mappings:
            outcome = Classification.DISCARDED_UNANNOTATED
            self.classification_qc.record(outcome)
            return self._discard_code(outcome)

        evidence, gene_id = extract_molecule_evidence(molitem, self.ir_introns)
        outcome = classify_molecule(evidence)
        self.classification_qc.record(outcome)

        if outcome.value in self.layers:
            if gene_id is None:
                raise ValueError("classified molecule does not have a unique gene")
            gene_index = geneid2ix[gene_id]
            dict_layers_columns[outcome.value][gene_index, cell_bcidx] += 1
            return 0
        return self._discard_code(outcome)


__all__ = ["PriorAwareLogic", "extract_molecule_evidence"]
