import os
import unittest

from plantvelo.species import (
    SpeciesResourceError,
    discover_species,
    resolve_builtin_species,
    validate_species_name,
)


class SpeciesTests(unittest.TestCase):
    def test_discovers_species_from_tsv_stems(self):
        self.assertEqual(discover_species(), ("ath", "gmx", "osa", "sly"))

    def test_resolves_builtin_species_with_stable_logical_path(self):
        source = resolve_builtin_species("osa")
        self.assertEqual(source.name, "osa")
        self.assertEqual(source.filename, "osa.tsv")
        self.assertTrue(os.path.isfile(source.path))
        self.assertEqual(source.logical_path, "builtin://osa.tsv")

    def test_rejects_unsafe_species_name(self):
        with self.assertRaises(SpeciesResourceError):
            validate_species_name("../osa")

    def test_rejects_unknown_species(self):
        with self.assertRaises(SpeciesResourceError):
            resolve_builtin_species("unknown")


if __name__ == "__main__":
    unittest.main()
