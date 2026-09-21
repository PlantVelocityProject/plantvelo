import unittest

from plantvelo.classification import (
    Classification,
    IntronKey,
    MoleculeEvidence,
    classify_molecule,
)


IR1 = IntronKey("G1", "Chr1", 100, 200, "+")
IR2 = IntronKey("G1", "Chr1", 300, 400, "+")


class ClassificationTests(unittest.TestCase):
    def test_exon_only_is_spliced(self):
        evidence = MoleculeEvidence(exonic=True)
        self.assertEqual(classify_molecule(evidence), Classification.SPLICED)

    def test_ordinary_intron_is_unspliced(self):
        evidence = MoleculeEvidence(ordinary_unspliced=True)
        self.assertEqual(classify_molecule(evidence), Classification.UNSPLICED)

    def test_high_confidence_retention_is_retained(self):
        evidence = MoleculeEvidence(hc_retained={IR1})
        self.assertEqual(classify_molecule(evidence), Classification.RETAINED)

    def test_unspliced_precedes_retained_and_spliced(self):
        evidence = MoleculeEvidence(
            ordinary_unspliced=True,
            hc_retained={IR1},
            exonic=True,
        )
        self.assertEqual(classify_molecule(evidence), Classification.UNSPLICED)

    def test_retained_precedes_exonic(self):
        evidence = MoleculeEvidence(hc_retained={IR1}, exonic=True)
        self.assertEqual(classify_molecule(evidence), Classification.RETAINED)

    def test_same_high_confidence_intron_conflict_is_ambiguous(self):
        evidence = MoleculeEvidence(hc_retained={IR1}, hc_spliced={IR1})
        self.assertEqual(classify_molecule(evidence), Classification.AMBIGUOUS)

    def test_different_retained_and_spliced_introns_is_retained(self):
        evidence = MoleculeEvidence(hc_retained={IR1}, hc_spliced={IR2})
        self.assertEqual(classify_molecule(evidence), Classification.RETAINED)

    def test_model_conflict_is_ambiguous_before_state_precedence(self):
        evidence = MoleculeEvidence(
            ordinary_unspliced=True,
            hc_retained={IR1},
            model_conflict=True,
        )
        self.assertEqual(classify_molecule(evidence), Classification.AMBIGUOUS)

    def test_multigene_is_discarded(self):
        evidence = MoleculeEvidence(multigene=True)
        self.assertEqual(
            classify_molecule(evidence), Classification.DISCARDED_MULTIGENE
        )

    def test_invalid_alignment_is_discarded(self):
        evidence = MoleculeEvidence(invalid_alignment=True)
        self.assertEqual(
            classify_molecule(evidence),
            Classification.DISCARDED_INVALID_ALIGNMENT,
        )

    def test_no_evidence_is_ambiguous(self):
        self.assertEqual(
            classify_molecule(MoleculeEvidence()), Classification.AMBIGUOUS
        )


if __name__ == "__main__":
    unittest.main()
