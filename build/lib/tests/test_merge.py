import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import loompy
import numpy as np

from plantvelo.loom_schema import (
    PRIOR_AWARE_SCHEMA_ID,
    VELOCYTO_DEFAULT_SCHEMA_ID,
    get_loom_schema,
)
from plantvelo.merge import (
    MergeValidationError,
    _read_file_attrs,
    merge_loom_files,
    validate_loom_inputs,
)


PRIOR_AWARE_ATTRS = {
    "plantvelo_version": "0.2.0",
    "classification_schema_version": PRIOR_AWARE_SCHEMA_ID,
    "counting_mode": "prior-aware",
    "classification_precedence": "U>R>S",
    "ir_prior_sha256": "prior-sha",
    "ir_prior_filter": "IR_class=high_confidence_IR;min_protocols=1",
    "matched_ir_introns": 2,
    "matched_ir_genes": 2,
    "gtf_sha256": "gtf-sha",
    "gtf_path": "/reference/genes.gtf",
    "ir_prior_path": "/reference/prior.tsv",
}

VELOCYTO_DEFAULT_ATTRS = {
    "plantvelo_version": "0.2.0",
    "classification_schema_version": VELOCYTO_DEFAULT_SCHEMA_ID,
    "counting_mode": "off",
    "velocyto_version": "0.17.17",
    "velocyto_logic": "Default",
    "velocyto_logic_resolved": "Permissive10X",
    "gtf_sha256": "gtf-sha",
    "gtf_path": "/reference/genes.gtf",
}


def create_plantvelo_loom(
    path,
    *,
    schema_id=PRIOR_AWARE_SCHEMA_ID,
    accessions=("AT1G00010", "AT1G00020"),
    cell_ids=("sample_AAAC-1",),
    omit_layer=None,
    extra_layers=None,
    corrupt_default=False,
    layer_dtype=np.uint32,
    row_attr_overrides=None,
    extra_row_attrs=None,
    extra_col_attrs=None,
    file_attr_overrides=None,
    omit_file_attrs=(),
):
    profile = get_loom_schema(schema_id)
    n_genes = len(accessions)
    n_cells = len(cell_ids)
    shape = (n_genes, n_cells)
    gene_number = {
        "AT1G00010": 10,
        "AT1G00020": 20,
    }
    base = np.asarray(
        [gene_number.get(accession, 10) for accession in accessions],
        dtype=layer_dtype,
    )[:, np.newaxis]
    base = np.repeat(base, n_cells, axis=1)
    all_matrices = {
        "spliced": base,
        "unspliced": base + 1,
        "retained": base + 2,
        "ambiguous": base + 3,
    }
    matrices = {name: all_matrices[name] for name in profile.layers}
    if omit_layer is not None:
        del matrices[omit_layer]
    if extra_layers:
        matrices.update(extra_layers)

    default = np.zeros(shape, dtype=np.float32)
    for matrix in matrices.values():
        default += matrix
    if corrupt_default:
        default[0, 0] += 1

    layers = {"": default}
    layers.update(matrices)
    row_attrs = {
        "Accession": np.asarray(accessions),
        "Gene": np.asarray(accessions),
        "Chromosome": np.asarray(["1"] * n_genes),
        "Start": np.asarray(
            [gene_number.get(accession, 10) for accession in accessions]
        ),
        "End": np.asarray(
            [gene_number.get(accession, 10) + 5 for accession in accessions]
        ),
        "Strand": np.asarray(["+"] * n_genes),
    }
    if row_attr_overrides:
        row_attrs.update(row_attr_overrides)
    if extra_row_attrs:
        row_attrs.update(extra_row_attrs)
    col_attrs = {"CellID": np.asarray(cell_ids)}
    if extra_col_attrs:
        col_attrs.update(extra_col_attrs)
    file_attrs = dict(
        PRIOR_AWARE_ATTRS
        if schema_id == PRIOR_AWARE_SCHEMA_ID
        else VELOCYTO_DEFAULT_ATTRS
    )
    if file_attr_overrides:
        file_attrs.update(file_attr_overrides)
    for name in omit_file_attrs:
        file_attrs.pop(name, None)
    loompy.create(
        str(path),
        layers,
        row_attrs,
        col_attrs,
        file_attrs=file_attrs,
    )
    return Path(path)


class MergeValidationTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def loom(self, name, **kwargs):
        return create_plantvelo_loom(self.root / name, **kwargs)

    def test_normalizes_zero_dimensional_array_file_attribute(self):
        class Dataset:
            attrs = {"scalar": np.asarray("value")}

        self.assertEqual(_read_file_attrs(Dataset()), {"scalar": "value"})

    def test_valid_inputs_return_summary(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom("second.loom", cell_ids=("s2_AAAC-1",))

        summary = validate_loom_inputs([first, second], batch_size=1)

        self.assertEqual(summary.gene_count, 2)
        self.assertEqual(summary.cell_count, 2)
        self.assertEqual(
            summary.cell_ids,
            frozenset({"s1_AAAC-1", "s2_AAAC-1"}),
        )

    def test_accepts_two_velocyto_default_inputs(self):
        first = self.loom(
            "first.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s1_AAAC-1",),
        )
        second = self.loom(
            "second.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s2_AAAC-1",),
        )

        summary = validate_loom_inputs([first, second], batch_size=1)

        self.assertEqual(summary.schema_id, VELOCYTO_DEFAULT_SCHEMA_ID)
        self.assertEqual(
            summary.layers,
            ("spliced", "unspliced", "ambiguous"),
        )

    def test_rejects_mixed_schema_inputs(self):
        first = self.loom("prior.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "off.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s2_AAAC-1",),
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "prior-aware-v1.*velocyto-default-v1",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_external_three_layer_loom_without_schema(self):
        invalid = self.loom(
            "external.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            omit_file_attrs=("classification_schema_version",),
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "unknown classification schema",
        ):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_velocyto_version_mismatch(self):
        first = self.loom(
            "first.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s1_AAAC-1",),
        )
        second = self.loom(
            "second.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s2_AAAC-1",),
            file_attr_overrides={"velocyto_version": "0.18.0"},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "velocyto_version",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_velocyto_logic_mismatch(self):
        first = self.loom(
            "first.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s1_AAAC-1",),
        )
        second = self.loom(
            "second.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            cell_ids=("s2_AAAC-1",),
            file_attr_overrides={
                "velocyto_logic_resolved": "Intermediate10X"
            },
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "velocyto_logic_resolved",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_accepts_legacy_prior_aware_without_counting_mode(self):
        legacy = self.loom(
            "legacy.loom",
            omit_file_attrs=("counting_mode",),
        )

        summary = validate_loom_inputs([legacy], batch_size=1)

        self.assertEqual(summary.schema_id, PRIOR_AWARE_SCHEMA_ID)

    def test_rejects_bad_velocyto_default_total(self):
        invalid = self.loom(
            "bad-off-total.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            corrupt_default=True,
        )

        with self.assertRaisesRegex(MergeValidationError, "default matrix"):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_invalid_loom_file(self):
        invalid = self.root / "invalid.loom"
        invalid.write_text("not hdf5", encoding="utf-8")

        with self.assertRaisesRegex(MergeValidationError, "invalid.loom"):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_missing_required_layer(self):
        invalid = self.loom("missing.loom", omit_layer="retained")

        with self.assertRaisesRegex(
            MergeValidationError,
            "missing required layer 'retained'",
        ):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_unexpected_layer(self):
        invalid = self.loom(
            "extra.loom",
            extra_layers={
                "legacy_retained": np.zeros((2, 1), dtype=np.uint32)
            },
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "unexpected layer 'legacy_retained'",
        ):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_duplicate_accession(self):
        invalid = self.loom(
            "duplicate-genes.loom",
            accessions=("AT1G00010", "AT1G00010"),
        )

        with self.assertRaisesRegex(MergeValidationError, "duplicate Accession"):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_duplicate_cellid_within_file(self):
        invalid = self.loom(
            "duplicate-cells.loom",
            cell_ids=("sample_AAAC-1", "sample_AAAC-1"),
        )

        with self.assertRaisesRegex(MergeValidationError, "duplicate CellID"):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_default_matrix_that_breaks_layer_conservation(self):
        invalid = self.loom("bad-total.loom", corrupt_default=True)

        with self.assertRaisesRegex(MergeValidationError, "default matrix"):
            validate_loom_inputs([invalid], batch_size=1)

    def test_rejects_duplicate_cellid_across_files(self):
        first = self.loom("first.loom", cell_ids=("same_AAAC-1",))
        second = self.loom("second.loom", cell_ids=("same_AAAC-1",))

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*duplicate CellID across inputs",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_mismatched_accession_set(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            accessions=("AT1G00010", "AT1G99999"),
            cell_ids=("s2_AAAC-1",),
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*Accession set",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_accepts_different_gene_order_and_keeps_first_order(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            accessions=("AT1G00020", "AT1G00010"),
            cell_ids=("s2_AAAC-1",),
        )

        summary = validate_loom_inputs([first, second], batch_size=1)

        self.assertEqual(summary.accessions, ("AT1G00010", "AT1G00020"))

    def test_rejects_row_attribute_value_mismatch_after_alignment(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            row_attr_overrides={"Start": np.asarray([999, 20])},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*row attribute 'Start'",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_row_attribute_name_mismatch(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            extra_row_attrs={"Extra": np.asarray([1, 2])},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*row attribute names",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_column_attribute_name_mismatch(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            extra_col_attrs={"Sample": np.asarray(["s2"])},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*column attribute names",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_column_attribute_dtype_mismatch(self):
        first = self.loom(
            "first.loom",
            cell_ids=("s1_AAAC-1",),
            extra_col_attrs={"Batch": np.asarray([1], dtype=np.int32)},
        )
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            extra_col_attrs={"Batch": np.asarray([1.5], dtype=np.float32)},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*column attribute dtypes",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_named_layer_dtype_mismatch(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            layer_dtype=np.uint16,
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*layer dtypes",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_semantic_attribute_mismatch(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            file_attr_overrides={"ir_prior_sha256": "different-prior"},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*ir_prior_sha256",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_rejects_unknown_file_attribute_mismatch(self):
        first = self.loom("first.loom", cell_ids=("s1_AAAC-1",))
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            file_attr_overrides={"future_semantics": "new-value"},
        )

        with self.assertRaisesRegex(
            MergeValidationError,
            "second.loom.*file attribute names",
        ):
            validate_loom_inputs([first, second], batch_size=1)

    def test_accepts_different_nonsemantic_provenance(self):
        first = self.loom(
            "first.loom",
            cell_ids=("s1_AAAC-1",),
            file_attr_overrides={
                "CreationDate": "2026-01-01T00:00:00Z",
                "gtf_path": "/first/genes.gtf",
                "ir_prior_path": "/first/prior.tsv",
            },
        )
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            file_attr_overrides={
                "CreationDate": "2026-02-01T00:00:00Z",
                "gtf_path": "/second/genes.gtf",
                "ir_prior_path": "/second/prior.tsv",
            },
        )

        summary = validate_loom_inputs([first, second], batch_size=1)

        self.assertEqual(summary.cell_count, 2)


class MergeExecutionTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def loom(self, name, **kwargs):
        return create_plantvelo_loom(self.root / name, **kwargs)

    def inputs(self):
        first = self.loom(
            "first.loom",
            cell_ids=("s1_AAAC-1",),
            extra_col_attrs={"Sample": np.asarray(["s1"])},
        )
        second = self.loom(
            "second.loom",
            accessions=("AT1G00020", "AT1G00010"),
            cell_ids=("s2_AAAC-1",),
            extra_col_attrs={"Sample": np.asarray(["s2"])},
        )
        return first, second

    def test_combines_two_sampled_velocyto_default_looms(self):
        first = self.loom(
            "sample-1-small.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            accessions=("AT1G00010", "AT1G00020"),
            cell_ids=("s1_AAAC-1", "s1_AAAG-1"),
            extra_col_attrs={
                "Sample": np.asarray(["s1", "s1"]),
            },
        )
        second = self.loom(
            "sample-2-small.loom",
            schema_id=VELOCYTO_DEFAULT_SCHEMA_ID,
            accessions=("AT1G00020", "AT1G00010"),
            cell_ids=("s2_AAAC-1", "s2_AAAG-1"),
            extra_col_attrs={
                "Sample": np.asarray(["s2", "s2"]),
            },
        )
        output = self.root / "combined-small-off.loom"

        summary = merge_loom_files(
            [first, second],
            output,
            batch_size=1,
        )

        self.assertEqual(summary.schema_id, VELOCYTO_DEFAULT_SCHEMA_ID)
        self.assertEqual(summary.gene_count, 2)
        self.assertEqual(summary.cell_count, 4)
        with loompy.connect(str(output), mode="r") as ds:
            self.assertEqual(
                set(ds.layers.keys()),
                {"", "spliced", "unspliced", "ambiguous"},
            )
            self.assertEqual(
                list(ds.ra["Accession"]),
                ["AT1G00010", "AT1G00020"],
            )
            self.assertEqual(
                set(ds.ca["CellID"]),
                {
                    "s1_AAAC-1",
                    "s1_AAAG-1",
                    "s2_AAAC-1",
                    "s2_AAAG-1",
                },
            )
            np.testing.assert_array_equal(
                ds[:, :],
                ds.layers["spliced"][:, :]
                + ds.layers["unspliced"][:, :]
                + ds.layers["ambiguous"][:, :],
            )
            self.assertEqual(
                ds.attrs["classification_schema_version"],
                VELOCYTO_DEFAULT_SCHEMA_ID,
            )
            self.assertEqual(ds.attrs["merge_input_count"], 2)

    def test_combines_real_looms_and_revalidates_output(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"

        summary = merge_loom_files(
            [first, second],
            output,
            batch_size=1,
        )

        self.assertEqual(summary.gene_count, 2)
        self.assertEqual(summary.cell_count, 2)
        with loompy.connect(str(output), mode="r") as ds:
            self.assertEqual(ds.shape, (2, 2))
            self.assertEqual(
                list(ds.ra["Accession"]),
                ["AT1G00010", "AT1G00020"],
            )
            self.assertEqual(
                set(ds.ca["CellID"]),
                {"s1_AAAC-1", "s2_AAAC-1"},
            )
            self.assertEqual(set(ds.ca["Sample"]), {"s1", "s2"})
            self.assertEqual(
                set(ds.layers.keys()),
                {"", "spliced", "unspliced", "retained", "ambiguous"},
            )
            np.testing.assert_array_equal(
                ds.layers["spliced"][:, :],
                np.asarray([[10, 10], [20, 20]], dtype=np.uint32),
            )
            self.assertEqual(ds.attrs["merge_key"], "Accession")
            self.assertEqual(ds.attrs["merge_input_count"], 2)
            self.assertEqual(
                json.loads(ds.attrs["merge_inputs"]),
                [str(first.resolve()), str(second.resolve())],
            )
            self.assertIn("merge_created_at", ds.attrs)
            self.assertIn("CreationDate", ds.attrs)

    def test_preserves_nan_column_attribute_values(self):
        first = self.loom(
            "first.loom",
            cell_ids=("s1_AAAC-1",),
            extra_col_attrs={"Score": np.asarray([np.nan])},
        )
        second = self.loom(
            "second.loom",
            cell_ids=("s2_AAAC-1",),
            extra_col_attrs={"Score": np.asarray([np.nan])},
        )
        output = self.root / "combined.loom"

        merge_loom_files([first, second], output, batch_size=1)

        with loompy.connect(str(output), mode="r") as ds:
            np.testing.assert_array_equal(
                ds.ca["Score"],
                np.asarray([np.nan, np.nan]),
            )

    def test_rejects_fewer_than_two_inputs(self):
        first = self.loom("first.loom")

        with self.assertRaisesRegex(MergeValidationError, "at least two"):
            merge_loom_files([first], self.root / "combined.loom")

    def test_rejects_duplicate_resolved_input_paths(self):
        first = self.loom("first.loom")

        with self.assertRaisesRegex(MergeValidationError, "duplicate input path"):
            merge_loom_files(
                [first, self.root / "." / "first.loom"],
                self.root / "combined.loom",
            )

    def test_rejects_output_that_is_also_an_input(self):
        first, second = self.inputs()

        with self.assertRaisesRegex(MergeValidationError, "also an input"):
            merge_loom_files([first, second], first, force=True)

    def test_rejects_missing_output_parent(self):
        first, second = self.inputs()

        with self.assertRaisesRegex(MergeValidationError, "parent.*does not exist"):
            merge_loom_files(
                [first, second],
                self.root / "missing" / "combined.loom",
            )

    def test_rejects_existing_output_without_force(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")

        with self.assertRaisesRegex(MergeValidationError, "already exists"):
            merge_loom_files([first, second], output)

        self.assertEqual(output.read_bytes(), b"old output")

    def test_no_force_does_not_replace_output_created_during_merge(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        from plantvelo.merge import _validate_merged_output

        def validate_then_create_destination(*args, **kwargs):
            result = _validate_merged_output(*args, **kwargs)
            output.write_bytes(b"concurrent output")
            return result

        with patch(
            "plantvelo.merge._validate_merged_output",
            side_effect=validate_then_create_destination,
        ):
            with self.assertRaisesRegex(MergeValidationError, "already exists"):
                merge_loom_files([first, second], output, batch_size=1)

        self.assertEqual(output.read_bytes(), b"concurrent output")

    def test_force_atomically_replaces_existing_output(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")

        merge_loom_files([first, second], output, force=True, batch_size=1)

        validator = loompy.LoomValidator()
        self.assertTrue(validator.validate(str(output), strictness="speconly"))

    def test_combine_failure_preserves_existing_output(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")

        with patch("loompy.combine", side_effect=RuntimeError("combine failed")):
            with self.assertRaisesRegex(RuntimeError, "combine failed"):
                merge_loom_files([first, second], output, force=True)

        self.assertEqual(output.read_bytes(), b"old output")

    def test_postvalidation_failure_preserves_existing_output(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")

        with patch(
            "plantvelo.merge._validate_merged_output",
            side_effect=MergeValidationError("postvalidation failed"),
        ):
            with self.assertRaisesRegex(
                MergeValidationError,
                "postvalidation failed",
            ):
                merge_loom_files([first, second], output, force=True)

        self.assertEqual(output.read_bytes(), b"old output")

    def test_postvalidation_rejects_changed_semantic_attribute(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")
        real_combine = loompy.combine

        def combine_with_changed_semantics(*args, **kwargs):
            real_combine(*args, **kwargs)
            with loompy.connect(kwargs["output_file"]) as ds:
                ds.attrs["ir_prior_sha256"] = "changed-after-combine"

        with patch("loompy.combine", side_effect=combine_with_changed_semantics):
            with self.assertRaisesRegex(
                MergeValidationError,
                "file attribute 'ir_prior_sha256'",
            ):
                merge_loom_files(
                    [first, second], output, force=True, batch_size=1
                )

        self.assertEqual(output.read_bytes(), b"old output")

    def test_postvalidation_rejects_changed_merge_provenance(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")
        real_combine = loompy.combine

        def combine_with_changed_provenance(*args, **kwargs):
            real_combine(*args, **kwargs)
            with loompy.connect(kwargs["output_file"]) as ds:
                ds.attrs["merge_key"] = "Gene"

        with patch("loompy.combine", side_effect=combine_with_changed_provenance):
            with self.assertRaisesRegex(
                MergeValidationError,
                "merge provenance 'merge_key'",
            ):
                merge_loom_files(
                    [first, second], output, force=True, batch_size=1
                )

        self.assertEqual(output.read_bytes(), b"old output")

    def test_postvalidation_rejects_changed_row_attribute_values(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")
        real_combine = loompy.combine

        def combine_with_changed_row_attribute(*args, **kwargs):
            real_combine(*args, **kwargs)
            with loompy.connect(kwargs["output_file"]) as ds:
                ds.ra["Gene"] = np.asarray(["wrong-1", "wrong-2"])

        with patch(
            "loompy.combine",
            side_effect=combine_with_changed_row_attribute,
        ):
            with self.assertRaisesRegex(
                MergeValidationError,
                "row attribute 'Gene'",
            ):
                merge_loom_files(
                    [first, second], output, force=True, batch_size=1
                )

        self.assertEqual(output.read_bytes(), b"old output")

    def test_postvalidation_rejects_changed_column_attribute_values(self):
        first, second = self.inputs()
        output = self.root / "combined.loom"
        output.write_bytes(b"old output")
        real_combine = loompy.combine

        def combine_with_changed_column_attribute(*args, **kwargs):
            real_combine(*args, **kwargs)
            with loompy.connect(kwargs["output_file"]) as ds:
                ds.ca["Sample"] = np.asarray(["wrong-1", "wrong-2"])

        with patch(
            "loompy.combine",
            side_effect=combine_with_changed_column_attribute,
        ):
            with self.assertRaisesRegex(
                MergeValidationError,
                "column attribute 'Sample'",
            ):
                merge_loom_files(
                    [first, second], output, force=True, batch_size=1
                )

        self.assertEqual(output.read_bytes(), b"old output")


if __name__ == "__main__":
    unittest.main()
