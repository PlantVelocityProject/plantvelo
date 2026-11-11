"""Strict validation and merging of PlantVelo Loom outputs."""
from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import FrozenSet, Mapping, Sequence, Tuple

import numpy as np

from plantvelo.loom_schema import (
    LoomSchemaError,
    get_loom_schema,
    normalize_semantic_attrs,
)

REQUIRED_ROW_ATTRS = (
    "Accession",
    "Gene",
    "Chromosome",
    "Start",
    "End",
    "Strand",
)
PROVENANCE_ATTRS = frozenset(
    {
        "CreationDate",
        "gtf_path",
        "ir_prior_path",
        "merge_key",
        "merge_input_count",
        "merge_inputs",
        "merge_created_at",
    }
)


class MergeValidationError(ValueError):
    """Raised when Loom inputs cannot be merged without semantic loss."""


@dataclass(frozen=True)
class MergeSummary:
    """Validated shape, identity, and metadata shared by merge inputs."""

    input_paths: Tuple[str, ...]
    schema_id: str
    layers: Tuple[str, ...]
    gene_count: int
    cell_count: int
    accessions: Tuple[str, ...]
    cell_ids: FrozenSet[str]
    row_attrs: Tuple[str, ...]
    col_attrs: Tuple[str, ...]
    col_attr_dtypes: Tuple[Tuple[str, str], ...]
    layer_dtypes: Tuple[Tuple[str, str], ...]
    file_attrs: Mapping[str, object]


def _normalize_scalar(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _normalize_values(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind in {"O", "S", "U"}:
        return np.asarray([_normalize_scalar(value) for value in array])
    return array


def _array_values_equal(left, right) -> bool:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    if np.array_equal(left_array, right_array):
        return True
    if left_array.shape != right_array.shape:
        return False
    if left_array.dtype.kind in {"f", "c"} and right_array.dtype.kind in {
        "f",
        "c",
    }:
        return np.array_equal(left_array, right_array, equal_nan=True)
    return False


def _find_duplicate(values: Sequence[str]):
    seen = set()
    for value in values:
        if value in seen:
            return value
        seen.add(value)
    return None


def _validate_loom_spec(path: Path) -> None:
    import loompy

    validator = loompy.LoomValidator()
    try:
        valid = validator.validate(str(path), strictness="speconly")
    except Exception as error:
        raise MergeValidationError(
            "{}: invalid Loom file: {}".format(path, error)
        ) from error
    if not valid:
        details = "; ".join(validator.errors) or "unknown specification error"
        raise MergeValidationError(
            "{}: invalid Loom file: {}".format(path, details)
        )


def _validate_layer_conservation(
    ds,
    path: Path,
    batch_size: int,
    layers: Sequence[str],
) -> None:
    if batch_size <= 0:
        raise ValueError("batch_size must be greater than zero")
    for start in range(0, ds.shape[1], batch_size):
        stop = min(start + batch_size, ds.shape[1])
        expected = np.zeros(
            (ds.shape[0], stop - start),
            dtype=ds.layers[""].dtype,
        )
        for layer in layers:
            expected += np.asarray(
                ds.layers[layer][:, start:stop],
                dtype=expected.dtype,
            )
        observed = ds[:, start:stop]
        if not np.array_equal(observed, expected):
            raise MergeValidationError(
                "{}: default matrix does not equal the sum of schema layers "
                "for columns {}:{}".format(path, start, stop)
            )


def _read_file_attrs(ds) -> Mapping[str, object]:
    attrs = {}
    for name in ds.attrs.keys():
        value = ds.attrs[name]
        if isinstance(value, np.ndarray):
            if value.ndim == 0:
                attrs[name] = _normalize_scalar(value.item())
            else:
                attrs[name] = tuple(
                    _normalize_scalar(item) for item in value.tolist()
                )
        else:
            attrs[name] = _normalize_scalar(value)
    return attrs


def _equal_attr_values(left, right) -> bool:
    if isinstance(left, tuple) or isinstance(right, tuple):
        return _array_values_equal(left, right)
    if left == right:
        return True
    try:
        return bool(np.isnan(left) and np.isnan(right))
    except TypeError:
        return False


def validate_loom_inputs(
    input_paths: Sequence[object],
    *,
    batch_size: int = 1000,
) -> MergeSummary:
    """Validate PlantVelo Loom files and summarize their merge dimensions."""
    import loompy

    if not input_paths:
        raise MergeValidationError("the input file list is empty")

    paths = tuple(Path(path).expanduser().resolve() for path in input_paths)
    first_accessions = None
    first_row_attrs = None
    first_col_attrs = None
    first_col_attr_dtypes = None
    first_layer_dtypes = None
    first_file_attrs = None
    first_row_values = None
    first_comparable_attrs = None
    first_profile = None
    all_cell_ids = set()
    cell_owners = {}
    cell_count = 0

    for path in paths:
        _validate_loom_spec(path)
        try:
            ds = loompy.connect(str(path), mode="r")
        except Exception as error:
            raise MergeValidationError(
                "{}: cannot open Loom file: {}".format(path, error)
            ) from error

        try:
            raw_file_attrs = _read_file_attrs(ds)
            try:
                profile = get_loom_schema(
                    raw_file_attrs.get("classification_schema_version", "")
                )
                file_attrs = normalize_semantic_attrs(
                    profile,
                    raw_file_attrs,
                )
            except LoomSchemaError as error:
                raise MergeValidationError(
                    "{}: {}".format(path, error)
                ) from error

            if first_profile is None:
                first_profile = profile
            elif profile.schema_id != first_profile.schema_id:
                raise MergeValidationError(
                    "{}: cannot mix classification schemas '{}' and '{}'".format(
                        path,
                        first_profile.schema_id,
                        profile.schema_id,
                    )
                )

            named_layers = set(ds.layers.keys()) - {""}
            required_layers = set(profile.layers)
            missing_layers = required_layers - named_layers
            if missing_layers:
                layer = sorted(missing_layers)[0]
                raise MergeValidationError(
                    "{}: missing required layer '{}'".format(path, layer)
                )
            extra_layers = named_layers - required_layers
            if extra_layers:
                layer = sorted(extra_layers)[0]
                raise MergeValidationError(
                    "{}: unexpected layer '{}'".format(path, layer)
                )

            row_attr_names = tuple(sorted(ds.ra.keys()))
            missing_row_attrs = set(REQUIRED_ROW_ATTRS) - set(row_attr_names)
            if missing_row_attrs:
                name = sorted(missing_row_attrs)[0]
                raise MergeValidationError(
                    "{}: missing required row attribute '{}'".format(path, name)
                )
            col_attr_names = tuple(sorted(ds.ca.keys()))
            if "CellID" not in col_attr_names:
                raise MergeValidationError(
                    "{}: missing required column attribute 'CellID'".format(path)
                )
            col_attr_dtypes = tuple(
                (name, str(ds.ca[name].dtype)) for name in col_attr_names
            )

            accessions = tuple(
                str(value)
                for value in _normalize_values(ds.ra["Accession"])
            )
            duplicate_accession = _find_duplicate(accessions)
            if duplicate_accession is not None:
                raise MergeValidationError(
                    "{}: duplicate Accession '{}'".format(
                        path, duplicate_accession
                    )
                )

            cell_ids = tuple(
                str(value) for value in _normalize_values(ds.ca["CellID"])
            )
            duplicate_cell = _find_duplicate(cell_ids)
            if duplicate_cell is not None:
                raise MergeValidationError(
                    "{}: duplicate CellID '{}'".format(path, duplicate_cell)
                )

            _validate_layer_conservation(
                ds,
                path,
                batch_size,
                profile.layers,
            )

            layer_dtypes = tuple(
                (name, str(ds.layers[name].dtype))
                for name in ("",) + profile.layers
            )
            comparable_attrs = {
                name: value
                for name, value in file_attrs.items()
                if name not in PROVENANCE_ATTRS
            }
            row_values = {
                name: _normalize_values(ds.ra[name]) for name in row_attr_names
            }

            if first_accessions is None:
                first_accessions = accessions
                first_row_attrs = row_attr_names
                first_col_attrs = col_attr_names
                first_col_attr_dtypes = col_attr_dtypes
                first_layer_dtypes = layer_dtypes
                first_file_attrs = file_attrs
                first_row_values = row_values
                first_comparable_attrs = comparable_attrs
            else:
                if row_attr_names != first_row_attrs:
                    raise MergeValidationError(
                        "{}: row attribute names differ from the first input".format(
                            path
                        )
                    )
                if col_attr_names != first_col_attrs:
                    raise MergeValidationError(
                        "{}: column attribute names differ from the first input".format(
                            path
                        )
                    )
                if col_attr_dtypes != first_col_attr_dtypes:
                    raise MergeValidationError(
                        "{}: column attribute dtypes differ from the first "
                        "input".format(path)
                    )
                if layer_dtypes != first_layer_dtypes:
                    raise MergeValidationError(
                        "{}: layer dtypes differ from the first input".format(path)
                    )
                if set(accessions) != set(first_accessions):
                    raise MergeValidationError(
                        "{}: Accession set differs from the first input".format(path)
                    )

                index_by_accession = {
                    accession: index
                    for index, accession in enumerate(accessions)
                }
                indexer = np.asarray(
                    [index_by_accession[value] for value in first_accessions]
                )
                for name in first_row_attrs:
                    aligned = row_values[name][indexer]
                    if not _array_values_equal(
                        aligned, first_row_values[name]
                    ):
                        raise MergeValidationError(
                            "{}: row attribute '{}' differs after Accession "
                            "alignment".format(path, name)
                        )

                if set(comparable_attrs) != set(first_comparable_attrs):
                    raise MergeValidationError(
                        "{}: file attribute names differ from the first input".format(
                            path
                        )
                    )
                for name in sorted(first_comparable_attrs):
                    if not _equal_attr_values(
                        comparable_attrs[name], first_comparable_attrs[name]
                    ):
                        raise MergeValidationError(
                            "{}: file attribute '{}' differs from the first "
                            "input".format(path, name)
                        )

            for cell_id in cell_ids:
                if cell_id in cell_owners:
                    raise MergeValidationError(
                        "{}: duplicate CellID across inputs '{}' (first seen "
                        "in {})".format(path, cell_id, cell_owners[cell_id])
                    )
                cell_owners[cell_id] = path

            all_cell_ids.update(cell_ids)
            cell_count += len(cell_ids)
        finally:
            ds.close()

    return MergeSummary(
        input_paths=tuple(str(path) for path in paths),
        schema_id=first_profile.schema_id,
        layers=first_profile.layers,
        gene_count=len(first_accessions),
        cell_count=cell_count,
        accessions=first_accessions,
        cell_ids=frozenset(all_cell_ids),
        row_attrs=first_row_attrs,
        col_attrs=first_col_attrs,
        col_attr_dtypes=first_col_attr_dtypes,
        layer_dtypes=first_layer_dtypes,
        file_attrs=first_file_attrs,
    )


def _validate_merged_output(
    output_path: Path,
    expected: MergeSummary,
    expected_merge_attrs: Mapping[str, object],
    batch_size: int,
) -> MergeSummary:
    actual = validate_loom_inputs([output_path], batch_size=batch_size)
    checks = (
        (actual.schema_id == expected.schema_id, "classification schema"),
        (actual.layers == expected.layers, "layer names"),
        (actual.gene_count == expected.gene_count, "gene count"),
        (actual.cell_count == expected.cell_count, "cell count"),
        (actual.accessions == expected.accessions, "Accession order"),
        (actual.cell_ids == expected.cell_ids, "CellID set"),
        (actual.row_attrs == expected.row_attrs, "row attribute names"),
        (actual.col_attrs == expected.col_attrs, "column attribute names"),
        (
            actual.col_attr_dtypes == expected.col_attr_dtypes,
            "column attribute dtypes",
        ),
        (actual.layer_dtypes == expected.layer_dtypes, "layer dtypes"),
    )
    for condition, description in checks:
        if not condition:
            raise MergeValidationError(
                "{}: merged output {} differs from validated inputs".format(
                    output_path, description
                )
            )

    expected_semantic_attrs = {
        name: value
        for name, value in expected.file_attrs.items()
        if name not in PROVENANCE_ATTRS
    }
    actual_semantic_attrs = {
        name: value
        for name, value in actual.file_attrs.items()
        if name not in PROVENANCE_ATTRS
    }
    if set(actual_semantic_attrs) != set(expected_semantic_attrs):
        raise MergeValidationError(
            "{}: merged output file attribute names differ from validated "
            "inputs".format(output_path)
        )
    for name, expected_value in expected_semantic_attrs.items():
        if not _equal_attr_values(actual_semantic_attrs[name], expected_value):
            raise MergeValidationError(
                "{}: merged output file attribute '{}' differs from validated "
                "inputs".format(output_path, name)
            )

    for name, expected_value in expected_merge_attrs.items():
        if name not in actual.file_attrs or not _equal_attr_values(
            actual.file_attrs[name], expected_value
        ):
            raise MergeValidationError(
                "{}: merged output merge provenance '{}' is invalid".format(
                    output_path, name
                )
            )
    _validate_merged_attribute_values(output_path, expected)
    return actual


def _validate_merged_attribute_values(
    output_path: Path,
    expected: MergeSummary,
) -> None:
    import loompy

    with loompy.connect(str(output_path), mode="r") as output_ds:
        output_row_values = {
            name: _normalize_values(output_ds.ra[name])
            for name in expected.row_attrs
        }
        output_col_values = {
            name: _normalize_values(output_ds.ca[name])
            for name in expected.col_attrs
        }

    with loompy.connect(expected.input_paths[0], mode="r") as first_ds:
        for name in expected.row_attrs:
            reference = _normalize_values(first_ds.ra[name])
            if not _array_values_equal(output_row_values[name], reference):
                raise MergeValidationError(
                    "{}: merged output row attribute '{}' differs from "
                    "validated inputs".format(output_path, name)
                )

    offset = 0
    for input_path in expected.input_paths:
        with loompy.connect(input_path, mode="r") as input_ds:
            stop = offset + input_ds.shape[1]
            for name in expected.col_attrs:
                reference = _normalize_values(input_ds.ca[name])
                observed = output_col_values[name][offset:stop]
                if not _array_values_equal(observed, reference):
                    raise MergeValidationError(
                        "{}: merged output column attribute '{}' differs from "
                        "validated inputs".format(output_path, name)
                    )
            offset = stop


def _absolute_path(path: object) -> Path:
    return Path(os.path.abspath(os.path.expanduser(os.fspath(path))))


def _publish_merged_output(
    temporary_output: Path,
    output: Path,
    *,
    force: bool,
) -> None:
    if force:
        os.replace(str(temporary_output), str(output))
        return
    try:
        os.link(str(temporary_output), str(output))
    except FileExistsError as error:
        raise MergeValidationError(
            "output already exists; pass --force to replace it atomically: "
            "{}".format(output)
        ) from error


def merge_loom_files(
    input_paths: Sequence[object],
    output_path: object,
    *,
    force: bool = False,
    batch_size: int = 1000,
) -> MergeSummary:
    """Strictly validate, combine, revalidate, and publish Loom inputs."""
    import loompy

    if len(input_paths) < 2:
        raise MergeValidationError("at least two input Loom files are required")

    resolved_inputs = tuple(
        Path(path).expanduser().resolve() for path in input_paths
    )
    if len(set(resolved_inputs)) != len(resolved_inputs):
        raise MergeValidationError("duplicate input path after resolution")
    for path in resolved_inputs:
        if not path.is_file():
            raise MergeValidationError("{}: input is not a file".format(path))

    output = _absolute_path(output_path)
    output_identity = output.resolve(strict=False)
    if output_identity in set(resolved_inputs):
        raise MergeValidationError("output path is also an input Loom file")
    if not output.parent.exists():
        raise MergeValidationError(
            "output parent directory does not exist: {}".format(output.parent)
        )
    if not output.parent.is_dir():
        raise MergeValidationError(
            "output parent is not a directory: {}".format(output.parent)
        )
    if output.exists() and not force:
        raise MergeValidationError(
            "output already exists; pass --force to replace it atomically: "
            "{}".format(output)
        )

    summary = validate_loom_inputs(resolved_inputs, batch_size=batch_size)
    now = datetime.now(timezone.utc)
    merge_attrs = {
        "CreationDate": now.strftime("%Y%m%dT%H%M%S.%fZ"),
        "merge_key": "Accession",
        "merge_input_count": len(resolved_inputs),
        "merge_inputs": json.dumps([str(path) for path in resolved_inputs]),
        "merge_created_at": now.isoformat().replace("+00:00", "Z"),
    }

    with tempfile.TemporaryDirectory(
        prefix=".plantvelo-merge-",
        dir=str(output.parent),
    ) as tempdir:
        temporary_output = Path(tempdir) / output.name
        loompy.combine(
            files=[str(path) for path in resolved_inputs],
            output_file=str(temporary_output),
            key="Accession",
            file_attrs=merge_attrs,
            batch_size=batch_size,
        )
        _validate_merged_output(
            temporary_output,
            summary,
            merge_attrs,
            batch_size,
        )
        _publish_merged_output(temporary_output, output, force=force)

    return summary
