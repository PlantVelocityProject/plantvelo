"""Loading and exact GTF matching for high-confidence IR priors."""
from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Set, Tuple

from plantvelo.classification import IntronKey


REQUIRED_FIELDS = (
    "gene_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "IR_class",
)
HIGH_CONFIDENCE_CLASS = "high_confidence_IR"


class IRPriorError(ValueError):
    """Raised when an IR prior cannot be validated or matched safely."""


@dataclass(frozen=True)
class PriorRecord:
    key: IntronKey
    values: Mapping[str, str]


@dataclass(frozen=True)
class PriorLoadResult:
    path: str
    sha256: str
    fieldnames: Tuple[str, ...]
    input_rows: int
    high_confidence_rows: int
    duplicate_rows_collapsed: int
    records: Tuple[PriorRecord, ...]
    source_kind: str = "custom"
    species: Optional[str] = None
    filter_description: str = "IR_class=high_confidence_IR"


@dataclass(frozen=True)
class UnmatchedPrior:
    record: PriorRecord
    reason: str


@dataclass(frozen=True)
class PriorMatchResult:
    prior: PriorLoadResult
    runtime_introns: Set[IntronKey]
    matched_records: Tuple[PriorRecord, ...]
    unmatched: Tuple[UnmatchedPrior, ...]

    @property
    def matched_unique_introns(self) -> int:
        return len(self.runtime_introns)

    @property
    def matched_genes(self) -> int:
        return len({key.gene_id for key in self.runtime_introns})

    @property
    def unmatched_unique_introns(self) -> int:
        return len({item.record.key for item in self.unmatched})


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _delimiter(path: str) -> str:
    name = str(path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".csv"):
        return ","
    if name.endswith(".tsv"):
        return "\t"
    raise IRPriorError("IR prior must be .csv, .tsv, .csv.gz, or .tsv.gz")


def _open_text(path: str):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _normalize_interval(start_text: str, end_text: str, row_number: int) -> Tuple[int, int]:
    try:
        start = int(start_text)
        end = int(end_text)
    except (TypeError, ValueError):
        raise IRPriorError(
            "row {} start/end must be positive integers with start <= end".format(
                row_number
            )
        )
    if start < 1 or end < 1 or start > end:
        raise IRPriorError(
            "row {} start/end must be positive integers with start <= end".format(
                row_number
            )
        )
    return start - 1, end


def _load_prior_rows(
    path: str,
    required_fields: Tuple[str, ...],
    include_row: Callable[[Mapping[str, str]], bool],
    source_kind: str,
    species: Optional[str],
    filter_description: str,
) -> PriorLoadResult:
    """Load rows under a source-specific schema into the shared result model."""

    delimiter = _delimiter(path)
    input_rows = 0
    high_confidence_rows = 0
    duplicate_rows = 0
    records: List[PriorRecord] = []
    seen_rows: Set[Tuple[str, ...]] = set()

    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = tuple(reader.fieldnames or ())
        missing = [field for field in required_fields if field not in fieldnames]
        if missing:
            raise IRPriorError(
                "IR prior missing required fields: {}".format(", ".join(missing))
            )
        if source_kind == "builtin" and set(fieldnames) != set(required_fields):
            raise IRPriorError(
                "built-in IR prior must contain exactly: {}".format(
                    ", ".join(required_fields)
                )
            )

        for row_number, raw_row in enumerate(reader, start=2):
            input_rows += 1
            row = {field: (raw_row.get(field) or "").strip() for field in fieldnames}
            if not include_row(row):
                continue
            high_confidence_rows += 1

            start, end = _normalize_interval(
                row["start"], row["end"], row_number
            )
            if row["strand"] not in ("+", "-"):
                raise IRPriorError(
                    "row {} strand must be '+' or '-'".format(row_number)
                )
            if not row["gene_id"] or not row["chromosome"]:
                raise IRPriorError(
                    "row {} gene_id and chromosome must not be empty".format(
                        row_number
                    )
                )

            fingerprint = tuple(row[field] for field in fieldnames)
            if fingerprint in seen_rows:
                duplicate_rows += 1
                continue
            seen_rows.add(fingerprint)
            records.append(
                PriorRecord(
                    key=IntronKey(
                        row["gene_id"],
                        row["chromosome"],
                        start,
                        end,
                        row["strand"],
                    ),
                    values=row,
                )
            )

    if high_confidence_rows == 0:
        raise IRPriorError("IR prior contains no high-confidence IR records")

    coordinate_strands: Dict[Tuple[str, int, int], Set[str]] = defaultdict(set)
    for record in records:
        coordinate = (
            record.key.chromosome,
            record.key.start,
            record.key.end,
        )
        coordinate_strands[coordinate].add(record.key.strand)
    conflicts = [
        key for key, strands in coordinate_strands.items() if len(strands) > 1
    ]
    if conflicts:
        raise IRPriorError(
            "IR prior has a coordinate with conflicting strands: {}".format(
                conflicts[0]
            )
        )

    records.sort(key=lambda record: record.key)
    return PriorLoadResult(
        path=str(Path(path).resolve()),
        sha256=sha256_file(path),
        fieldnames=fieldnames,
        input_rows=input_rows,
        high_confidence_rows=high_confidence_rows,
        duplicate_rows_collapsed=duplicate_rows,
        records=tuple(records),
        source_kind=source_kind,
        species=species,
        filter_description=filter_description,
    )


def load_ir_prior(path: str) -> PriorLoadResult:
    """Load, filter, validate, normalize, and deduplicate a custom IR prior."""

    return _load_prior_rows(
        path=path,
        required_fields=REQUIRED_FIELDS,
        include_row=lambda row: row["IR_class"] == HIGH_CONFIDENCE_CLASS,
        source_kind="custom",
        species=None,
        filter_description="IR_class=high_confidence_IR;min_protocols=1",
    )


BUILTIN_FIELDS = (
    "gene_id",
    "chromosome",
    "start",
    "end",
    "strand",
)


def load_builtin_ir_prior(path: str, species: str) -> PriorLoadResult:
    """Load a bundled five-field prior; every row is high-confidence IR."""

    return _load_prior_rows(
        path=path,
        required_fields=BUILTIN_FIELDS,
        include_row=lambda row: True,
        source_kind="builtin",
        species=species,
        filter_description="all records treated as high-confidence IR",
    )


_GENE_ID = re.compile(r'gene_id\s+"([^"]+)"')
_TRANSCRIPT_ID = re.compile(r'transcript_id\s+"([^"]+)"')


@dataclass(frozen=True)
class _GTFIndex:
    introns: Set[IntronKey]
    chromosomes: Set[str]
    genes: Set[str]


def _raw_gtf_index(gtf_path: str) -> _GTFIndex:
    transcript_exons: Dict[
        Tuple[str, str, str, str], List[Tuple[int, int]]
    ] = defaultdict(list)
    chromosomes: Set[str] = set()
    genes: Set[str] = set()

    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise IRPriorError(
                    "GTF line {} does not have 9 fields".format(line_number)
                )
            chromosome, _, feature_type, start_text, end_text, _, strand, _, attrs = fields
            if feature_type != "exon":
                continue
            gene_match = _GENE_ID.search(attrs)
            transcript_match = _TRANSCRIPT_ID.search(attrs)
            if gene_match is None or transcript_match is None:
                raise IRPriorError(
                    "GTF exon line {} lacks gene_id or transcript_id".format(
                        line_number
                    )
                )
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                raise IRPriorError(
                    "GTF line {} has invalid coordinates".format(line_number)
                )
            transcript_exons[
                (chromosome, strand, gene_match.group(1), transcript_match.group(1))
            ].append((start, end))
            chromosomes.add(chromosome)
            genes.add(gene_match.group(1))

    introns: Set[IntronKey] = set()
    for (chromosome, strand, gene_id, _), exons in transcript_exons.items():
        ordered = sorted(exons)
        for left, right in zip(ordered, ordered[1:]):
            intron_start = left[1]
            intron_end = right[0] - 1
            if intron_start < intron_end:
                introns.add(
                    IntronKey(
                        gene_id,
                        chromosome,
                        intron_start,
                        intron_end,
                        strand,
                    )
                )
    return _GTFIndex(introns=introns, chromosomes=chromosomes, genes=genes)


def _velocyto_chromosome(chromosome: str) -> str:
    return chromosome[3:] if "chr" in chromosome[:4] else chromosome


def _runtime_key(key: IntronKey) -> IntronKey:
    return IntronKey(
        key.gene_id,
        _velocyto_chromosome(key.chromosome),
        key.start,
        key.end,
        key.strand,
    )


def _runtime_introns(annotations: Mapping[str, Mapping[str, object]]) -> Set[IntronKey]:
    introns: Set[IntronKey] = set()
    for transcripts in annotations.values():
        for transcript in transcripts.values():
            chromstrand = transcript.chromstrand
            for feature in transcript:
                if feature.kind == ord("i"):
                    introns.add(
                        IntronKey(
                            transcript.geneid,
                            chromstrand[:-1],
                            feature.start - 1,
                            feature.end,
                            chromstrand[-1],
                        )
                    )
    return introns


def _unmatched_reason(key: IntronKey, gtf_index: _GTFIndex) -> str:
    if key.chromosome not in gtf_index.chromosomes:
        return "chromosome_not_found"
    if key.gene_id not in gtf_index.genes:
        return "gene_not_found"
    coordinate = (key.gene_id, key.chromosome, key.start, key.end)
    if coordinate in {
        (item.gene_id, item.chromosome, item.start, item.end)
        for item in gtf_index.introns
    }:
        return "strand_mismatch"
    return "intron_coordinate_not_found"


def match_ir_prior(
    prior: PriorLoadResult,
    gtf_path: str,
    annotations: Mapping[str, Mapping[str, object]],
) -> PriorMatchResult:
    """Match a validated prior to raw GTF introns and loaded features."""

    gtf_index = _raw_gtf_index(gtf_path)
    available_runtime = _runtime_introns(annotations)
    runtime_introns: Set[IntronKey] = set()
    matched_records: List[PriorRecord] = []
    unmatched: List[UnmatchedPrior] = []

    for record in prior.records:
        if record.key not in gtf_index.introns:
            unmatched.append(
                UnmatchedPrior(record, _unmatched_reason(record.key, gtf_index))
            )
            continue
        candidate = _runtime_key(record.key)
        if candidate not in available_runtime:
            unmatched.append(UnmatchedPrior(record, "runtime_feature_not_found"))
            continue
        runtime_introns.add(candidate)
        matched_records.append(record)

    if not runtime_introns:
        raise IRPriorError("matched 0 high-confidence IR introns to the GTF")

    return PriorMatchResult(
        prior=prior,
        runtime_introns=runtime_introns,
        matched_records=tuple(matched_records),
        unmatched=tuple(unmatched),
    )
