import unittest

import numpy as np

import plantvelo.logic as logic_module
from plantvelo.classification import (
    Classification,
    IntronKey,
)

PriorAwareLogic = logic_module.PriorAwareLogic
extract_molecule_evidence = logic_module.extract_molecule_evidence


class FakeTranscript:
    def __init__(self, geneid="G1", chromstrand="Chr1+"):
        self.geneid = geneid
        self.chromstrand = chromstrand
        self.list_features = []

    def __iter__(self):
        return iter(self.list_features)


class FakeFeature:
    def __init__(self, transcript, start, end, kind):
        self.transcript_model = transcript
        self.start = start
        self.end = end
        self.kind = ord(kind)
        transcript.list_features.append(self)


class FakeSegmentMatch:
    def __init__(self, segment, feature, is_spliced=False):
        self.segment = segment
        self.feature = feature
        self.is_spliced = is_spliced

    @property
    def maps_to_intron(self):
        return self.feature.kind == ord("i")

    @property
    def maps_to_exon(self):
        return self.feature.kind == ord("e")


class FakeMolecule:
    def __init__(self, mappings_record):
        self.mappings_record = mappings_record


def intron_key(feature):
    transcript = feature.transcript_model
    return IntronKey(
        transcript.geneid,
        transcript.chromstrand[:-1],
        feature.start - 1,
        feature.end,
        transcript.chromstrand[-1],
    )


class LogicTests(unittest.TestCase):
    def setUp(self):
        self.layers = {
            name: np.zeros((2, 1), dtype=np.uint32)
            for name in ("spliced", "unspliced", "retained", "ambiguous")
        }

    def make_model(self, geneid="G1"):
        transcript = FakeTranscript(geneid)
        exon1 = FakeFeature(transcript, 1, 100, "e")
        intron1 = FakeFeature(transcript, 101, 200, "i")
        exon2 = FakeFeature(transcript, 201, 300, "e")
        intron2 = FakeFeature(transcript, 301, 400, "i")
        exon3 = FakeFeature(transcript, 401, 500, "e")
        return transcript, exon1, intron1, exon2, intron2, exon3

    def test_exposes_only_new_layers_and_no_legacy_classes(self):
        logic = PriorAwareLogic(set())
        self.assertEqual(
            logic.layers, ["spliced", "unspliced", "retained", "ambiguous"]
        )
        self.assertFalse(hasattr(logic_module, "PlantPermissive10X"))
        self.assertFalse(hasattr(logic_module, "PlantValidated10X"))

    def test_ordinary_intron_internal_segment_is_unspliced(self):
        transcript, _, intron, _, _, _ = self.make_model()
        molecule = FakeMolecule(
            {transcript: [FakeSegmentMatch((120, 150), intron)]}
        )
        evidence, gene_id = extract_molecule_evidence(molecule, set())
        self.assertEqual(gene_id, "G1")
        self.assertTrue(evidence.ordinary_unspliced)

    def test_ordinary_intron_boundary_segment_is_unspliced(self):
        transcript, exon, intron, _, _, _ = self.make_model()
        segment = (90, 110)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch(segment, exon),
                    FakeSegmentMatch(segment, intron),
                ]
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, set())
        self.assertTrue(evidence.ordinary_unspliced)

    def test_high_confidence_internal_segment_is_retained(self):
        transcript, _, intron, _, _, _ = self.make_model()
        key = intron_key(intron)
        molecule = FakeMolecule(
            {transcript: [FakeSegmentMatch((120, 150), intron)]}
        )
        evidence, _ = extract_molecule_evidence(molecule, {key})
        self.assertEqual(evidence.hc_retained, {key})
        self.assertFalse(evidence.ordinary_unspliced)

    def test_high_confidence_boundary_segment_is_retained(self):
        transcript, exon, intron, _, _, _ = self.make_model()
        key = intron_key(intron)
        segment = (90, 110)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch(segment, exon),
                    FakeSegmentMatch(segment, intron),
                ]
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, {key})
        self.assertEqual(evidence.hc_retained, {key})

    def test_exact_spliced_junction_over_high_confidence_ir_is_spliced(self):
        transcript, exon1, intron, exon2, _, _ = self.make_model()
        key = intron_key(intron)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch((50, 100), exon1, is_spliced=True),
                    FakeSegmentMatch((201, 250), exon2, is_spliced=True),
                ]
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, {key})
        self.assertEqual(evidence.hc_spliced, {key})
        self.assertTrue(evidence.exonic)

    def test_non_exact_splice_does_not_mark_high_confidence_ir(self):
        transcript, exon1, intron, exon2, _, _ = self.make_model()
        key = intron_key(intron)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch((50, 99), exon1, is_spliced=True),
                    FakeSegmentMatch((201, 250), exon2, is_spliced=True),
                ]
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, {key})
        self.assertEqual(evidence.hc_spliced, set())

    def test_same_segment_with_different_model_states_is_conflict(self):
        exon_model = FakeTranscript("G1")
        exon = FakeFeature(exon_model, 1, 300, "e")
        intron_model = FakeTranscript("G1")
        intron = FakeFeature(intron_model, 101, 200, "i")
        segment = (120, 150)
        molecule = FakeMolecule(
            {
                exon_model: [FakeSegmentMatch(segment, exon)],
                intron_model: [FakeSegmentMatch(segment, intron)],
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, set())
        self.assertTrue(evidence.model_conflict)

    def test_different_read_segments_support_u_and_r_without_model_conflict(self):
        transcript, _, ordinary, _, retained, _ = self.make_model()
        retained_key = intron_key(retained)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch((120, 150), ordinary),
                    FakeSegmentMatch((320, 350), retained),
                ]
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, {retained_key})
        self.assertFalse(evidence.model_conflict)
        self.assertTrue(evidence.ordinary_unspliced)
        self.assertEqual(evidence.hc_retained, {retained_key})

    def test_one_segment_can_support_ordinary_and_retained_introns(self):
        transcript, _, ordinary, _, retained, _ = self.make_model()
        retained_key = intron_key(retained)
        segment = (90, 410)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch(segment, ordinary),
                    FakeSegmentMatch(segment, retained),
                ]
            }
        )
        evidence, _ = extract_molecule_evidence(molecule, {retained_key})
        self.assertFalse(evidence.model_conflict)
        self.assertTrue(evidence.ordinary_unspliced)
        self.assertEqual(evidence.hc_retained, {retained_key})
        self.assertEqual(
            logic_module.classify_molecule(evidence), Classification.UNSPLICED
        )

    def test_same_ir_retained_and_spliced_counts_ambiguous_once(self):
        transcript, exon1, intron, exon2, _, _ = self.make_model()
        key = intron_key(intron)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch((120, 150), intron),
                    FakeSegmentMatch((50, 100), exon1, is_spliced=True),
                    FakeSegmentMatch((201, 250), exon2, is_spliced=True),
                ]
            }
        )
        logic = PriorAwareLogic({key})
        result = logic.count(molecule, 0, self.layers, {"G1": 0})
        self.assertEqual(result, 0)
        self.assertEqual(self.layers["ambiguous"][0, 0], 1)
        self.assertEqual(logic.classification_qc.ambiguous, 1)

    def test_u_precedes_r_when_counting(self):
        transcript, _, ordinary, _, retained, _ = self.make_model()
        retained_key = intron_key(retained)
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch((120, 150), ordinary),
                    FakeSegmentMatch((320, 350), retained),
                ]
            }
        )
        logic = PriorAwareLogic({retained_key})
        logic.count(molecule, 0, self.layers, {"G1": 0})
        self.assertEqual(self.layers["unspliced"][0, 0], 1)
        self.assertEqual(self.layers["retained"][0, 0], 0)

    def test_registry_can_be_populated_after_logic_construction(self):
        registry = set()
        transcript, _, intron, _, _, _ = self.make_model()
        key = intron_key(intron)
        logic = PriorAwareLogic(registry)
        registry.add(key)
        molecule = FakeMolecule(
            {transcript: [FakeSegmentMatch((120, 150), intron)]}
        )
        logic.count(molecule, 0, self.layers, {"G1": 0})
        self.assertEqual(self.layers["retained"][0, 0], 1)

    def test_multigene_is_discarded_and_recorded(self):
        model1, exon1, _, _, _, _ = self.make_model("G1")
        model2, exon2, _, _, _, _ = self.make_model("G2")
        molecule = FakeMolecule(
            {
                model1: [FakeSegmentMatch((10, 40), exon1)],
                model2: [FakeSegmentMatch((10, 40), exon2)],
            }
        )
        logic = PriorAwareLogic(set())
        result = logic.count(molecule, 0, self.layers, {"G1": 0, "G2": 1})
        self.assertEqual(result, 1)
        self.assertEqual(logic.classification_qc.discarded_multigene, 1)
        self.assertEqual(sum(layer.sum() for layer in self.layers.values()), 0)

    def test_empty_mapping_is_discarded_unannotated(self):
        logic = PriorAwareLogic(set())
        result = logic.count(FakeMolecule({}), 0, self.layers, {"G1": 0})
        self.assertEqual(result, 2)
        self.assertEqual(logic.classification_qc.discarded_unannotated, 1)

    def test_malformed_segment_is_discarded_invalid_alignment(self):
        transcript, exon, _, _, _, _ = self.make_model()
        molecule = FakeMolecule(
            {transcript: [FakeSegmentMatch((100,), exon)]}
        )
        logic = PriorAwareLogic(set())
        result = logic.count(molecule, 0, self.layers, {"G1": 0})
        self.assertEqual(result, 4)
        self.assertEqual(logic.classification_qc.discarded_invalid_alignment, 1)


if __name__ == "__main__":
    unittest.main()
