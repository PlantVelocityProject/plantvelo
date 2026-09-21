import unittest

from plantvelo.species import resolve_prior_source


class PriorSourceTests(unittest.TestCase):
    def test_resolves_custom_source(self):
        source = resolve_prior_source(None, "/tmp/prior.tsv")
        self.assertEqual(source.source_kind, "custom")
        self.assertIsNone(source.species)

    def test_resolves_builtin_source(self):
        source = resolve_prior_source("osa", None)
        self.assertEqual(source.source_kind, "builtin")
        self.assertEqual(source.species, "osa")
        self.assertEqual(source.path, source.resource.path)


if __name__ == "__main__":
    unittest.main()
