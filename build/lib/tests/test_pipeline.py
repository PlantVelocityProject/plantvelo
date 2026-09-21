import ast
import csv
import os
import tempfile
import unittest
from collections import OrderedDict
from types import SimpleNamespace

from plantvelo.classification import IntronKey
from plantvelo.pipeline import bind_ir_prior, load_annotations_for_mode


class FakeFeature:
    kind = ord("i")

    def __init__(self):
        self.start = 101
        self.end = 200


class FakeTranscript:
    geneid = "G1"
    chromstrand = "Chr1+"

    def __iter__(self):
        return iter((FakeFeature(),))


class FakeCounter:
    def __init__(self, annotations):
        self.annotations = annotations
        self.loaded_path = None

    def read_transcriptmodels(self, path):
        self.loaded_path = path
        return self.annotations


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.prior_path = os.path.join(self.tempdir.name, "prior.csv")
        with open(self.prior_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "gene_id",
                    "chromosome",
                    "start",
                    "end",
                    "strand",
                    "IR_class",
                ),
            )
            writer.writeheader()
            writer.writerow(
                {
                    "gene_id": "G1",
                    "chromosome": "Chr1",
                    "start": "101",
                    "end": "200",
                    "strand": "+",
                    "IR_class": "high_confidence_IR",
                }
            )
        self.gtf_path = os.path.join(self.tempdir.name, "genes.gtf")
        with open(self.gtf_path, "w", encoding="utf-8") as handle:
            handle.write(
                'Chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "G1"; '
                'transcript_id "T1"; exon_number "1";\n'
            )
            handle.write(
                'Chr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id "G1"; '
                'transcript_id "T1"; exon_number "2";\n'
            )

    def tearDown(self):
        self.tempdir.cleanup()

    def test_loads_annotations_and_populates_shared_registry(self):
        annotations = {
            "Chr1+": OrderedDict((("T1", FakeTranscript()),))
        }
        counter = FakeCounter(annotations)
        registry = set()
        returned_annotations, match = bind_ir_prior(
            counter, self.prior_path, self.gtf_path, registry
        )
        self.assertIs(returned_annotations, annotations)
        self.assertEqual(counter.loaded_path, self.gtf_path)
        self.assertEqual(
            registry, {IntronKey("G1", "Chr1", 100, 200, "+")}
        )
        self.assertEqual(match.matched_unique_introns, 1)

    def test_rejects_rebinding_nonempty_registry(self):
        counter = FakeCounter({})
        with self.assertRaisesRegex(ValueError, "empty"):
            bind_ir_prior(counter, self.prior_path, self.gtf_path, {object()})

    def test_off_mode_loads_gtf_without_prior_binding(self):
        annotations = {
            "Chr1+": OrderedDict((("T1", FakeTranscript()),))
        }
        counter = FakeCounter(annotations)
        profile = SimpleNamespace(requires_ir_prior=False)

        returned_annotations, match = load_annotations_for_mode(
            counter,
            profile,
            None,
            self.gtf_path,
            set(),
        )

        self.assertIs(returned_annotations, annotations)
        self.assertIsNone(match)
        self.assertEqual(counter.loaded_path, self.gtf_path)

    def test_prior_mode_delegates_to_prior_binding(self):
        annotations = {
            "Chr1+": OrderedDict((("T1", FakeTranscript()),))
        }
        counter = FakeCounter(annotations)
        registry = set()
        profile = SimpleNamespace(requires_ir_prior=True)

        _, match = load_annotations_for_mode(
            counter,
            profile,
            self.prior_path,
            self.gtf_path,
            registry,
        )

        self.assertEqual(match.matched_unique_introns, 1)
        self.assertEqual(
            registry,
            {IntronKey("G1", "Chr1", 100, 200, "+")},
        )

    def test_run_core_signature_requires_prior_and_has_no_legacy_parameters(self):
        path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "plantvelo",
            "commands",
            "_run.py",
        )
        with open(path, "r", encoding="utf-8") as handle:
            tree = ast.parse(handle.read())
        function = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "_run"
        )
        arguments = [argument.arg for argument in function.args.kwonlyargs]
        self.assertIn("ir_mode", arguments)
        self.assertIn("ir_prior", arguments)
        self.assertNotIn("logic", arguments)
        self.assertNotIn("ir_flanking", arguments)


if __name__ == "__main__":
    unittest.main()
