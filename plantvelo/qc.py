"""Quality-control counters and deterministic TSV output."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple

from plantvelo.classification import Classification
from plantvelo.ir_prior import PriorMatchResult


CLASSIFICATION_METRICS = (
    "spliced",
    "unspliced",
    "retained",
    "ambiguous",
    "discarded_multigene",
    "discarded_unannotated",
    "discarded_invalid_alignment",
)


@dataclass
class ClassificationQC:
    processed_molecules: int = 0
    spliced: int = 0
    unspliced: int = 0
    retained: int = 0
    ambiguous: int = 0
    discarded_multigene: int = 0
    discarded_unannotated: int = 0
    discarded_invalid_alignment: int = 0

    def record(self, outcome: Classification) -> None:
        name = outcome.value
        if name not in CLASSIFICATION_METRICS:
            raise ValueError("unsupported classification outcome: {}".format(name))
        setattr(self, name, getattr(self, name) + 1)
        self.processed_molecules += 1

    def validate_conservation(self) -> None:
        classified = sum(getattr(self, name) for name in CLASSIFICATION_METRICS)
        if classified != self.processed_molecules:
            raise ValueError(
                "classification conservation failed: processed={}, accounted={}".format(
                    self.processed_molecules, classified
                )
            )

    def rows(self) -> Tuple[Tuple[str, object], ...]:
        return (("processed_molecules", self.processed_molecules),) + tuple(
            (name, getattr(self, name)) for name in CLASSIFICATION_METRICS
        )


def _write_metrics(path: str, rows: Iterable[Tuple[str, object]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(rows)


def write_classification_qc(path: str, qc: ClassificationQC) -> None:
    qc.validate_conservation()
    _write_metrics(path, qc.rows())


def write_prior_qc(
    path: str, match: PriorMatchResult, gtf_annotation_path: str
) -> None:
    prior = match.prior
    rows = (
        ("input_rows", prior.input_rows),
        ("high_confidence_rows", prior.high_confidence_rows),
        ("duplicate_rows_collapsed", prior.duplicate_rows_collapsed),
        ("matched_unique_introns", match.matched_unique_introns),
        ("matched_genes", match.matched_genes),
        ("unmatched_introns", match.unmatched_unique_introns),
        ("gtf_annotation_path", str(Path(gtf_annotation_path).resolve())),
        ("ir_prior_path", prior.path),
        ("ir_prior_sha256", prior.sha256),
        ("ir_prior_source", getattr(prior, "source_kind", "custom")),
        ("species", getattr(prior, "species", None) or "none"),
        (
            "ir_prior_filter",
            getattr(prior, "filter_description", "IR_class=high_confidence_IR;min_protocols=1"),
        ),
    )
    _write_metrics(path, rows)


def write_unmatched_prior(path: str, match: PriorMatchResult) -> None:
    fieldnames = tuple(match.prior.fieldnames) + ("reason",)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for unmatched in match.unmatched:
            row = dict(unmatched.record.values)
            row["reason"] = unmatched.reason
            writer.writerow(row)
