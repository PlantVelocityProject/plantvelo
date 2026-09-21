import os
import tempfile
import unittest

import numpy as np

from plantvelo.loom_schema import (
    PRIOR_AWARE_PROFILE,
    VELOCYTO_DEFAULT_PROFILE,
)
from plantvelo.output import (
    build_file_attrs,
    build_velocyto_file_attrs,
    create_loom_atomic,
    prepare_loom_layers,
)


class OutputTests(unittest.TestCase):
    def matrices(self):
        return {
            "spliced": np.array([[1, 0]], dtype=np.uint32),
            "unspliced": np.array([[0, 2]], dtype=np.uint32),
            "retained": np.array([[3, 0]], dtype=np.uint32),
            "ambiguous": np.array([[0, 4]], dtype=np.uint32),
        }

    def test_prepares_prior_aware_layers_and_total(self):
        layers = prepare_loom_layers(
            self.matrices(),
            "uint32",
            PRIOR_AWARE_PROFILE,
        )
        self.assertEqual(
            tuple(layers),
            ("", "spliced", "unspliced", "retained", "ambiguous"),
        )
        np.testing.assert_array_equal(layers[""], np.array([[4, 6]]))
        self.assertEqual(layers[""].dtype, np.float32)
        self.assertEqual(layers["retained"].dtype, np.uint32)

    def test_prepares_velocyto_default_layers_and_total(self):
        matrices = self.matrices()
        del matrices["retained"]

        layers = prepare_loom_layers(
            matrices,
            "uint32",
            VELOCYTO_DEFAULT_PROFILE,
        )

        self.assertEqual(
            tuple(layers),
            ("", "spliced", "unspliced", "ambiguous"),
        )
        np.testing.assert_array_equal(layers[""], np.array([[1, 6]]))

    def test_rejects_missing_layer(self):
        matrices = self.matrices()
        del matrices["retained"]
        with self.assertRaisesRegex(ValueError, "prior-aware-v1"):
            prepare_loom_layers(
                matrices,
                "uint32",
                PRIOR_AWARE_PROFILE,
            )

    def test_rejects_layers_from_a_different_schema(self):
        with self.assertRaisesRegex(ValueError, "velocyto-default-v1"):
            prepare_loom_layers(
                self.matrices(),
                "uint32",
                VELOCYTO_DEFAULT_PROFILE,
            )

    def test_rejects_mismatched_layer_dimensions(self):
        matrices = self.matrices()
        matrices["ambiguous"] = np.zeros((2, 2))
        with self.assertRaisesRegex(ValueError, "dimensions"):
            prepare_loom_layers(
                matrices,
                "uint32",
                PRIOR_AWARE_PROFILE,
            )

    def test_builds_required_provenance_attributes(self):
        attrs = build_file_attrs(
            version="0.2.0",
            ir_prior_path="/data/prior.tsv",
            ir_prior_sha256="prior-sha",
            matched_ir_introns=2,
            matched_ir_genes=1,
            gtf_path="/data/genes.gtf",
            gtf_sha256="gtf-sha",
        )
        self.assertEqual(attrs["plantvelo_version"], "0.2.0")
        self.assertEqual(attrs["counting_mode"], "prior-aware")
        self.assertEqual(
            attrs["classification_schema_version"], "prior-aware-v1"
        )
        self.assertEqual(attrs["classification_precedence"], "U>R>S")
        self.assertEqual(
            attrs["ir_prior_filter"],
            "IR_class=high_confidence_IR;min_protocols=1",
        )
        self.assertEqual(attrs["matched_ir_introns"], 2)
        self.assertEqual(attrs["gtf_sha256"], "gtf-sha")

    def test_builds_velocyto_default_provenance(self):
        attrs = build_velocyto_file_attrs(
            version="0.2.0",
            velocyto_version="0.17.17",
            logic_requested="Default",
            logic_resolved="Permissive10X",
            gtf_path="/data/genes.gtf",
            gtf_sha256="gtf-sha",
        )

        self.assertEqual(attrs["counting_mode"], "off")
        self.assertEqual(
            attrs["classification_schema_version"],
            "velocyto-default-v1",
        )
        self.assertEqual(attrs["velocyto_version"], "0.17.17")
        self.assertEqual(attrs["velocyto_logic"], "Default")
        self.assertEqual(
            attrs["velocyto_logic_resolved"],
            "Permissive10X",
        )

    def test_builds_builtin_prior_provenance(self):
        attrs = build_file_attrs(
            version="0.2.0",
            ir_prior_path="builtin://osa.tsv",
            ir_prior_sha256="prior-sha",
            matched_ir_introns=2,
            matched_ir_genes=1,
            gtf_path="/data/genes.gtf",
            gtf_sha256="gtf-sha",
            source_kind="builtin",
            species="osa",
            filter_description="all records treated as high-confidence IR",
        )
        self.assertEqual(attrs["ir_prior_path"], "builtin://osa.tsv")
        self.assertEqual(attrs["ir_prior_source"], "builtin")
        self.assertEqual(attrs["species"], "osa")
        self.assertEqual(
            attrs["ir_prior_filter"],
            "all records treated as high-confidence IR",
        )

    def test_atomic_create_replaces_only_after_success(self):
        with tempfile.TemporaryDirectory() as tempdir:
            output = os.path.join(tempdir, "sample.loom")
            calls = []

            def create(path, layers, row_attrs, col_attrs, file_attrs):
                calls.append(path)
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write("complete")

            create_loom_atomic(output, {}, {}, {}, {}, create=create)
            self.assertEqual(calls, [output + ".tmp"])
            with open(output, "r", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "complete")
            self.assertFalse(os.path.exists(output + ".tmp"))

    def test_atomic_create_removes_partial_temporary_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            output = os.path.join(tempdir, "sample.loom")

            def create(path, layers, row_attrs, col_attrs, file_attrs):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write("partial")
                raise RuntimeError("loom failed")

            with self.assertRaisesRegex(RuntimeError, "loom failed"):
                create_loom_atomic(output, {}, {}, {}, {}, create=create)
            self.assertFalse(os.path.exists(output))
            self.assertFalse(os.path.exists(output + ".tmp"))


if __name__ == "__main__":
    unittest.main()
