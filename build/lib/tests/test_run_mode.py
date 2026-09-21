import unittest
from types import SimpleNamespace

from plantvelo.loom_schema import (
    PRIOR_AWARE_SCHEMA_ID,
    VELOCYTO_DEFAULT_SCHEMA_ID,
)
from plantvelo.run_mode import (
    RunModeError,
    resolve_run_mode,
    validate_logic_contract,
)


class FakeDefaultLogic:
    name = "Permissive10X"

    @property
    def layers(self):
        return ["spliced", "unspliced", "ambiguous"]


class FakePriorAwareLogic:
    pass


class RunModeTests(unittest.TestCase):
    def test_resolves_prior_aware_profile(self):
        registry = set()

        profile = resolve_run_mode(
            "prior-aware",
            "/data/prior.tsv",
            registry,
            prior_logic_class=FakePriorAwareLogic,
        )

        self.assertEqual(profile.schema.schema_id, PRIOR_AWARE_SCHEMA_ID)
        self.assertTrue(profile.requires_ir_prior)

    def test_resolves_off_to_velocyto_default(self):
        velocyto_module = SimpleNamespace(Default=FakeDefaultLogic)

        profile = resolve_run_mode(
            "off",
            None,
            set(),
            velocyto_module=velocyto_module,
        )

        self.assertIs(profile.logic_factory, FakeDefaultLogic)
        self.assertEqual(
            profile.schema.schema_id, VELOCYTO_DEFAULT_SCHEMA_ID
        )
        self.assertFalse(profile.requires_ir_prior)

    def test_direct_api_prior_aware_requires_prior(self):
        with self.assertRaisesRegex(RunModeError, "required"):
            resolve_run_mode(
                "prior-aware",
                None,
                set(),
                prior_logic_class=FakePriorAwareLogic,
            )

    def test_direct_api_off_rejects_prior(self):
        with self.assertRaisesRegex(RunModeError, "must not"):
            resolve_run_mode(
                "off",
                "/data/prior.tsv",
                set(),
                velocyto_module=SimpleNamespace(Default=FakeDefaultLogic),
            )

    def test_validates_actual_logic_layers(self):
        profile = resolve_run_mode(
            "off",
            None,
            set(),
            velocyto_module=SimpleNamespace(Default=FakeDefaultLogic),
        )

        validate_logic_contract(profile, FakeDefaultLogic())

    def test_rejects_changed_default_layers(self):
        profile = resolve_run_mode(
            "off",
            None,
            set(),
            velocyto_module=SimpleNamespace(Default=FakeDefaultLogic),
        )
        incompatible = SimpleNamespace(
            name="Permissive10X",
            layers=["spliced", "unspliced"],
        )

        with self.assertRaisesRegex(RunModeError, "unexpected layers"):
            validate_logic_contract(profile, incompatible)


if __name__ == "__main__":
    unittest.main()
