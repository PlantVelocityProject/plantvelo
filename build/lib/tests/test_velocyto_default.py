import unittest

import numpy as np
import velocyto as vcy


class FakeTranscript:
    def __init__(self, geneid):
        self.geneid = geneid


class FakeFeature:
    def __init__(
        self,
        *,
        validated=False,
        start_overlap=False,
        end_overlap=False,
        upstream=None,
        downstream=None,
    ):
        self.is_validated = validated
        self.is_last_3prime = False
        self._start_overlap = start_overlap
        self._end_overlap = end_overlap
        self._upstream = upstream
        self._downstream = downstream

    def start_overlaps_with_part_of(self, segment):
        return self._start_overlap

    def end_overlaps_with_part_of(self, segment):
        return self._end_overlap

    def get_upstream_exon(self):
        return self._upstream

    def get_downstream_exon(self):
        return self._downstream


class FakeSegmentMatch:
    def __init__(
        self,
        feature,
        *,
        maps_to_exon=False,
        maps_to_intron=False,
        is_spliced=False,
    ):
        self.feature = feature
        self.segment = (10, 20)
        self.maps_to_exon = maps_to_exon
        self.maps_to_intron = maps_to_intron
        self.is_spliced = is_spliced


class FakeMolecule:
    def __init__(self, mappings_record):
        self.mappings_record = mappings_record


class VelocytoDefaultContractTests(unittest.TestCase):
    def setUp(self):
        self.logic = vcy.Default()
        self.layers = {
            name: np.zeros((2, 1), dtype=np.uint32)
            for name in ("spliced", "unspliced", "ambiguous")
        }
        self.geneid2ix = {"G1": 0, "G2": 1}

    def test_default_resolves_to_permissive_10x(self):
        self.assertEqual(self.logic.name, "Permissive10X")
        self.assertEqual(
            tuple(self.logic.layers),
            ("spliced", "unspliced", "ambiguous"),
        )

    def test_exon_only_counts_spliced(self):
        transcript = FakeTranscript("G1")
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch(FakeFeature(), maps_to_exon=True)
                ]
            }
        )

        result = self.logic.count(
            molecule, 0, self.layers, self.geneid2ix
        )

        self.assertEqual(result, 0)
        self.assertEqual(self.layers["spliced"][0, 0], 1)

    def test_intron_only_counts_unspliced(self):
        transcript = FakeTranscript("G1")
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch(FakeFeature(), maps_to_intron=True)
                ]
            }
        )

        result = self.logic.count(
            molecule, 0, self.layers, self.geneid2ix
        )

        self.assertEqual(result, 0)
        self.assertEqual(self.layers["unspliced"][0, 0], 1)

    def test_exon_intron_spanning_counts_unspliced(self):
        downstream_exon = FakeFeature(start_overlap=True)
        intron = FakeFeature(
            validated=True,
            end_overlap=True,
            downstream=downstream_exon,
        )
        transcript = FakeTranscript("G1")
        molecule = FakeMolecule(
            {
                transcript: [
                    FakeSegmentMatch(intron, maps_to_intron=True)
                ]
            }
        )

        result = self.logic.count(
            molecule, 0, self.layers, self.geneid2ix
        )

        self.assertEqual(result, 0)
        self.assertEqual(self.layers["unspliced"][0, 0], 1)

    def test_intron_only_and_exon_only_models_count_ambiguous(self):
        intron_transcript = FakeTranscript("G1")
        exon_transcript = FakeTranscript("G1")
        molecule = FakeMolecule(
            {
                intron_transcript: [
                    FakeSegmentMatch(FakeFeature(), maps_to_intron=True)
                ],
                exon_transcript: [
                    FakeSegmentMatch(FakeFeature(), maps_to_exon=True)
                ],
            }
        )

        result = self.logic.count(
            molecule, 0, self.layers, self.geneid2ix
        )

        self.assertEqual(result, 0)
        self.assertEqual(self.layers["ambiguous"][0, 0], 1)

    def test_unannotated_molecule_is_discarded(self):
        result = self.logic.count(
            FakeMolecule({}), 0, self.layers, self.geneid2ix
        )

        self.assertEqual(result, 2)

    def test_multigene_molecule_is_discarded(self):
        first = FakeTranscript("G1")
        second = FakeTranscript("G2")
        molecule = FakeMolecule(
            {
                first: [
                    FakeSegmentMatch(FakeFeature(), maps_to_exon=True)
                ],
                second: [
                    FakeSegmentMatch(FakeFeature(), maps_to_exon=True)
                ],
            }
        )

        result = self.logic.count(
            molecule, 0, self.layers, self.geneid2ix
        )

        self.assertEqual(result, 3)


if __name__ == "__main__":
    unittest.main()
