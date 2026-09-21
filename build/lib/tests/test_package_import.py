import os
import subprocess
import sys
import unittest


ROOT = os.path.dirname(os.path.dirname(__file__))


class PackageImportTests(unittest.TestCase):
    def test_commands_package_import_is_lazy(self):
        code = (
            "import sys; "
            "sys.path.insert(0, {!r}); "
            "import plantvelo.commands; "
            "print('imported')"
        ).format(ROOT)
        environment = dict(os.environ)
        environment.pop("PYTHONPATH", None)
        result = subprocess.run(
            [sys.executable, "-c", code],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "imported")


if __name__ == "__main__":
    unittest.main()
