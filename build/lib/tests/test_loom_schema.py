import unittest

from plantvelo.loom_schema import (
    PRIOR_AWARE_SCHEMA_ID,
    VELOCYTO_DEFAULT_SCHEMA_ID,
    LoomSchemaError,
    get_loom_schema,
    normalize_semantic_attrs,
)


class LoomSchemaTests(unittest.TestCase):
    def test_prior_aware_profile_has_four_layers(self):
        profile = get_loom_schema(PRIOR_AWARE_SCHEMA_ID)

        self.assertEqual(
            profile.layers,
            ("spliced", "unspliced", "retained", "ambiguous"),
        )
        self.assertEqual(profile.counting_mode, "prior-aware")

    def test_velocyto_default_profile_has_three_layers(self):
        profile = get_loom_schema(VELOCYTO_DEFAULT_SCHEMA_ID)

        self.assertEqual(
            profile.layers,
            ("spliced", "unspliced", "ambiguous"),
        )
        self.assertEqual(profile.counting_mode, "off")

    def test_rejects_unknown_schema(self):
        with self.assertRaisesRegex(
            LoomSchemaError, "unknown classification schema"
        ):
            get_loom_schema("unknown-v1")

    def test_normalizes_legacy_prior_aware_counting_mode(self):
        profile = get_loom_schema(PRIOR_AWARE_SCHEMA_ID)
        attrs = {
            "plantvelo_version": "0.2.0",
            "classification_schema_version": PRIOR_AWARE_SCHEMA_ID,
            "classification_precedence": "U>R>S",
            "ir_prior_sha256": "prior-sha",
            "ir_prior_filter": (
                "IR_class=high_confidence_IR;min_protocols=1"
            ),
            "matched_ir_introns": 2,
            "matched_ir_genes": 1,
            "gtf_sha256": "gtf-sha",
        }

        normalized = normalize_semantic_attrs(profile, attrs)

        self.assertEqual(normalized["counting_mode"], "prior-aware")

    def test_off_schema_requires_plantvelo_provenance(self):
        profile = get_loom_schema(VELOCYTO_DEFAULT_SCHEMA_ID)

        with self.assertRaisesRegex(
            LoomSchemaError, "missing semantic attribute"
        ):
            normalize_semantic_attrs(
                profile,
                {
                    "classification_schema_version": (
                        VELOCYTO_DEFAULT_SCHEMA_ID
                    ),
                    "counting_mode": "off",
                },
            )


if __name__ == "__main__":
    unittest.main()
