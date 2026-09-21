import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from click.testing import CliRunner

from plantvelo.commands.merge import merge
from plantvelo.commands.plantvelo import cli
from plantvelo.merge import MergeValidationError


class MergeCommandTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.first = self.root / "first.loom"
        self.second = self.root / "second.loom"
        self.first.touch()
        self.second.touch()
        self.output = self.root / "combined.loom"
        self.runner = CliRunner()

    def tearDown(self):
        self.tempdir.cleanup()

    def test_help_lists_inputs_output_and_force(self):
        result = self.runner.invoke(merge, ["--help"])

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("LOOMFILES", result.output)
        self.assertIn("--output", result.output)
        self.assertIn("--force", result.output)

    def test_rejects_one_input_as_usage_error(self):
        result = self.runner.invoke(
            merge,
            ["--output", str(self.output), str(self.first)],
        )

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("at least two", result.output)

    def test_forwards_inputs_and_output_to_core(self):
        summary = SimpleNamespace(gene_count=2, cell_count=3)
        with patch(
            "plantvelo.commands.merge.merge_loom_files",
            return_value=summary,
        ) as core:
            result = self.runner.invoke(
                merge,
                [
                    "--output",
                    str(self.output),
                    str(self.first),
                    str(self.second),
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        core.assert_called_once_with(
            (str(self.first), str(self.second)),
            str(self.output),
            force=False,
        )
        self.assertIn("2 genes x 3 cells", result.output)
        self.assertIn(str(self.output), result.output)

    def test_forwards_force(self):
        summary = SimpleNamespace(gene_count=2, cell_count=3)
        with patch(
            "plantvelo.commands.merge.merge_loom_files",
            return_value=summary,
        ) as core:
            result = self.runner.invoke(
                merge,
                [
                    "--force",
                    "--output",
                    str(self.output),
                    str(self.first),
                    str(self.second),
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        core.assert_called_once_with(
            (str(self.first), str(self.second)),
            str(self.output),
            force=True,
        )

    def test_presents_validation_error_without_traceback(self):
        with patch(
            "plantvelo.commands.merge.merge_loom_files",
            side_effect=MergeValidationError("strict validation failed"),
        ):
            result = self.runner.invoke(
                merge,
                [
                    "--output",
                    str(self.output),
                    str(self.first),
                    str(self.second),
                ],
            )

        self.assertEqual(result.exit_code, 1, result.output)
        self.assertIn("Error: strict validation failed", result.output)
        self.assertNotIn("Traceback", result.output)

    def test_top_level_help_lists_run_before_merge(self):
        result = self.runner.invoke(cli, ["--help"])

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("prior-aware or velocyto Default", result.output)
        self.assertNotIn("required high-confidence", result.output)
        run_position = result.output.index("run")
        merge_position = result.output.index("merge")
        self.assertLess(run_position, merge_position)


if __name__ == "__main__":
    unittest.main()
