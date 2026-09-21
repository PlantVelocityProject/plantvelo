import csv
import os
import tempfile
import unittest
from types import SimpleNamespace

from plantvelo.classification import Classification, IntronKey
from plantvelo.ir_prior import PriorRecord, UnmatchedPrior
from plantvelo.qc import (
    ClassificationQC,
    write_classification_qc,
    write_prior_qc,
    write_unmatched_prior,
)


class QCTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def path(self, name):
        return os.path.join(self.tempdir.name, name)

    def read_tsv(self, path):
        with open(path, "r", encoding="utf-8", newline="") as handle:
            return list(csv.reader(handle, delimiter="\t"))

    def test_records_every_classification_category(self):
        qc = ClassificationQC()
        for outcome in Classification:
            qc.record(outcome)
        qc.validate_conservation()
        self.assertEqual(qc.processed_molecules, len(Classification))
        for outcome in Classification:
            self.assertEqual(getattr(qc, outcome.value), 1)

    def test_classification_conservation_failure_is_fatal(self):
        qc = ClassificationQC(processed_molecules=2, spliced=1)
        with self.assertRaisesRegex(ValueError, "conservation"):
            qc.validate_conservation()

    def test_classification_qc_writer_has_stable_required_metrics(self):
        qc = ClassificationQC()
        qc.record(Classification.SPLICED)
        path = self.path("classification.tsv")
        write_classification_qc(path, qc)
        rows = self.read_tsv(path)
        self.assertEqual(rows[0], ["metric", "value"])
        self.assertEqual(
            [row[0] for row in rows[1:]],
            [
                "processed_molecules",
                "spliced",
                "unspliced",
                "retained",
                "ambiguous",
                "discarded_multigene",
                "discarded_unannotated",
                "discarded_invalid_alignment",
            ],
        )

    def test_prior_qc_writer_includes_provenance_and_counts(self):
        prior = SimpleNamespace(
            input_rows=4,
            high_confidence_rows=3,
            duplicate_rows_collapsed=1,
            path="/data/prior.tsv",
            sha256="abc123",
        )
        matched = SimpleNamespace(
            prior=prior,
            matched_unique_introns=1,
            matched_genes=1,
            unmatched=(object(),),
            unmatched_unique_introns=1,
        )
        path = self.path("prior_qc.tsv")
        write_prior_qc(path, matched, "/data/genes.gtf")
        rows = dict(self.read_tsv(path)[1:])
        self.assertEqual(rows["input_rows"], "4")
        self.assertEqual(rows["unmatched_introns"], "1")
        self.assertEqual(rows["gtf_annotation_path"], "/data/genes.gtf")
        self.assertEqual(rows["ir_prior_sha256"], "abc123")

    def test_prior_qc_writer_includes_source_and_species(self):
        prior = SimpleNamespace(
            input_rows=1,
            high_confidence_rows=1,
            duplicate_rows_collapsed=0,
            path="builtin://osa.tsv",
            sha256="abc123",
            source_kind="builtin",
            species="osa",
            filter_description="all records treated as high-confidence IR",
        )
        matched = SimpleNamespace(
            prior=prior,
            matched_unique_introns=1,
            matched_genes=1,
            unmatched=(),
            unmatched_unique_introns=0,
        )
        path = self.path("builtin_prior_qc.tsv")
        write_prior_qc(path, matched, "/data/genes.gtf")
        rows = dict(self.read_tsv(path)[1:])
        self.assertEqual(rows["ir_prior_source"], "builtin")
        self.assertEqual(rows["species"], "osa")
        self.assertEqual(
            rows["ir_prior_filter"],
            "all records treated as high-confidence IR",
        )

    def test_unmatched_writer_preserves_prior_columns_and_reason(self):
        record = PriorRecord(
            IntronKey("G1", "Chr1", 100, 200, "+"),
            {
                "gene_id": "G1",
                "chromosome": "Chr1",
                "start": "101",
                "end": "200",
                "strand": "+",
                "IR_class": "high_confidence_IR",
            },
        )
        prior = SimpleNamespace(
            fieldnames=(
                "gene_id",
                "chromosome",
                "start",
                "end",
                "strand",
                "IR_class",
            )
        )
        matched = SimpleNamespace(
            prior=prior,
            unmatched=(UnmatchedPrior(record, "chromosome_not_found"),),
        )
        path = self.path("unmatched.tsv")
        write_unmatched_prior(path, matched)
        rows = self.read_tsv(path)
        self.assertEqual(rows[0][-1], "reason")
        self.assertEqual(rows[1][-1], "chromosome_not_found")
        self.assertEqual(rows[1][0], "G1")


if __name__ == "__main__":
    unittest.main()
