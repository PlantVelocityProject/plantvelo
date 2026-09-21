import importlib.util
import inspect
import os
import sys
import tempfile
import types
import unittest

from click.testing import CliRunner


class CLITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calls = []
        fake_run_module = types.ModuleType("plantvelo.commands._run")

        def fake_run(**kwargs):
            cls.calls.append(kwargs)

        fake_run_module._run = fake_run
        cls.previous_run_module = sys.modules.get("plantvelo.commands._run")
        sys.modules["plantvelo.commands._run"] = fake_run_module

        source = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "plantvelo",
            "commands",
            "run.py",
        )
        spec = importlib.util.spec_from_file_location(
            "plantvelo_run_under_test", source
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        cls.module = module
        cls.command = module.run

    @classmethod
    def tearDownClass(cls):
        if cls.previous_run_module is None:
            sys.modules.pop("plantvelo.commands._run", None)
        else:
            sys.modules["plantvelo.commands._run"] = cls.previous_run_module

    def setUp(self):
        self.calls.clear()

    def invoke_with_inputs(self, extra_arguments):
        with tempfile.TemporaryDirectory() as tempdir:
            paths = {}
            for name in ("prior.csv", "reads.bam", "genes.gtf"):
                path = os.path.join(tempdir, name)
                with open(path, "w", encoding="utf-8"):
                    pass
                paths[name] = path
            result = CliRunner().invoke(
                self.command,
                list(extra_arguments(paths))
                + [paths["reads.bam"], paths["genes.gtf"]],
            )
        return result, paths

    def test_help_lists_ir_mode_and_removes_legacy_options(self):
        result = CliRunner().invoke(self.command, ["--help"])
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("--ir-mode", result.output)
        self.assertIn("prior-aware", result.output)
        self.assertIn("off", result.output)
        self.assertIn("--ir-prior", result.output)
        self.assertNotIn("--logic", result.output)
        self.assertNotIn("--ir-flanking", result.output)

    def test_help_keeps_general_velocyto_options(self):
        result = CliRunner().invoke(self.command, ["--help"])
        for option in (
            "--bcfile",
            "--outputfolder",
            "--mask",
            "--without-umi",
            "--umi-extension",
            "--multimap",
            "--samtools-threads",
            "--dtype",
        ):
            self.assertIn(option, result.output)

    def test_callback_signature_has_prior_but_no_legacy_arguments(self):
        parameters = inspect.signature(self.command.callback).parameters
        self.assertIn("ir_mode", parameters)
        self.assertIn("ir_prior", parameters)
        self.assertNotIn("logic", parameters)
        self.assertNotIn("ir_flanking", parameters)

    def test_prior_aware_requires_species_or_prior(self):
        result, _ = self.invoke_with_inputs(lambda paths: [])

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("one of --species or --ir-prior is required", result.output)
        self.assertEqual(self.calls, [])

    def test_species_forwards_species(self):
        result, _ = self.invoke_with_inputs(
            lambda paths: ["--species", "osa"]
        )

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(self.calls[0]["species"], "osa")

    def test_species_and_prior_are_mutually_exclusive(self):
        result, paths = self.invoke_with_inputs(
            lambda paths: [
                "--species",
                "osa",
                "--ir-prior",
                paths["prior.csv"],
            ]
        )

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("mutually exclusive", result.output)

    def test_off_rejects_species(self):
        result, _ = self.invoke_with_inputs(
            lambda paths: ["--ir-mode", "off", "--species", "osa"]
        )

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("must not", result.output)

    def test_prior_aware_forwards_prior_and_mode(self):
        result, paths = self.invoke_with_inputs(
            lambda paths: ["--ir-prior", paths["prior.csv"]]
        )

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(self.calls[0]["ir_mode"], "prior-aware")
        self.assertEqual(self.calls[0]["ir_prior"], paths["prior.csv"])
        self.assertNotIn("logic", self.calls[0])
        self.assertNotIn("ir_flanking", self.calls[0])

    def test_off_allows_missing_prior(self):
        result, _ = self.invoke_with_inputs(
            lambda paths: ["--ir-mode", "off"]
        )

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(self.calls[0]["ir_mode"], "off")
        self.assertIsNone(self.calls[0]["ir_prior"])

    def test_off_rejects_prior(self):
        result, paths = self.invoke_with_inputs(
            lambda paths: [
                "--ir-mode",
                "off",
                "--ir-prior",
                paths["prior.csv"],
            ]
        )

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("must not be provided", result.output)
        self.assertEqual(self.calls, [])


if __name__ == "__main__":
    unittest.main()
