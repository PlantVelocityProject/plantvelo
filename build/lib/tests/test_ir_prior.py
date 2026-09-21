import csv
import gzip
import os
import tempfile
import unittest
from collections import OrderedDict

from plantvelo.classification import IntronKey
from plantvelo.ir_prior import (
    IRPriorError,
    load_builtin_ir_prior,
    load_ir_prior,
    match_ir_prior,
)


FIELDS = [
    "gene_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "IR_class",
    "intron_id",
]


class FakeFeature:
    def __init__(self, start, end, kind=ord("i")):
        self.start = start
        self.end = end
        self.kind = kind


class FakeTranscript:
    def __init__(self, geneid, chromstrand, features):
        self.geneid = geneid
        self.chromstrand = chromstrand
        self.list_features = features

    def __iter__(self):
        return iter(self.list_features)


class IRPriorTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def write_prior(self, rows, suffix=".csv", fields=FIELDS):
        path = os.path.join(self.tempdir.name, "prior" + suffix)
        delimiter = "\t" if ".tsv" in suffix else ","
        opener = gzip.open if suffix.endswith(".gz") else open
        kwargs = {"mode": "wt", "newline": "", "encoding": "utf-8"}
        with opener(path, **kwargs) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def write_builtin_prior(self, rows):
        return self.write_prior(
            rows,
            suffix=".tsv",
            fields=("gene_id", "chromosome", "start", "end", "strand"),
        )

    def write_gtf(self, chromosome="Chr1", gene_id="G1", transcripts=("T1",)):
        path = os.path.join(self.tempdir.name, "annotation.gtf")
        with open(path, "w", encoding="utf-8") as handle:
            for transcript_id in transcripts:
                attrs = (
                    'gene_id "{}"; transcript_id "{}"; exon_number "1";'
                ).format(gene_id, transcript_id)
                handle.write(
                    "{}\ttest\texon\t1\t100\t.\t+\t.\t{}\n".format(
                        chromosome, attrs
                    )
                )
                attrs = (
                    'gene_id "{}"; transcript_id "{}"; exon_number "2";'
                ).format(gene_id, transcript_id)
                handle.write(
                    "{}\ttest\texon\t201\t300\t.\t+\t.\t{}\n".format(
                        chromosome, attrs
                    )
                )
        return path

    @staticmethod
    def row(**overrides):
        row = {
            "gene_id": "G1",
            "chromosome": "Chr1",
            "start": "101",
            "end": "200",
            "strand": "+",
            "IR_class": "high_confidence_IR",
            "intron_id": "I1",
        }
        row.update(overrides)
        return row

    @staticmethod
    def annotations(chromosome="Chr1", copies=1):
        transcripts = OrderedDict()
        for index in range(copies):
            transcripts["T{}".format(index)] = FakeTranscript(
                "G1", chromosome + "+", [FakeFeature(101, 200)]
            )
        return {chromosome + "+": transcripts}

    def test_loads_csv_and_normalizes_coordinates(self):
        result = load_ir_prior(self.write_prior([self.row()]))
        self.assertEqual(
            result.records[0].key, IntronKey("G1", "Chr1", 100, 200, "+")
        )
        self.assertEqual(result.input_rows, 1)
        self.assertEqual(result.high_confidence_rows, 1)

    def test_loads_tsv_gzip(self):
        result = load_ir_prior(self.write_prior([self.row()], ".tsv.gz"))
        self.assertEqual(len(result.records), 1)
        self.assertEqual(len(result.sha256), 64)

    def test_loads_builtin_five_field_tsv(self):
        row = {field: self.row()[field] for field in FIELDS[:5]}
        result = load_builtin_ir_prior(
            self.write_builtin_prior([row]),
            species="osa",
        )
        self.assertEqual(
            result.records[0].key, IntronKey("G1", "Chr1", 100, 200, "+")
        )
        self.assertEqual(result.input_rows, 1)
        self.assertEqual(result.high_confidence_rows, 1)
        self.assertEqual(result.source_kind, "builtin")
        self.assertEqual(result.species, "osa")
        self.assertEqual(
            result.filter_description,
            "all records treated as high-confidence IR",
        )

    def test_builtin_loader_requires_exact_five_fields(self):
        fields = ("gene_id", "chromosome", "start", "end")
        row = {field: self.row()[field] for field in fields}
        with self.assertRaisesRegex(IRPriorError, "missing required fields"):
            load_builtin_ir_prior(
                self.write_prior([row], suffix=".tsv", fields=fields),
                species="osa",
            )

    def test_rejects_missing_required_field(self):
        fields = [field for field in FIELDS if field != "strand"]
        row = {field: value for field, value in self.row().items() if field in fields}
        path = self.write_prior([row], fields=fields)
        with self.assertRaisesRegex(IRPriorError, "missing required fields.*strand"):
            load_ir_prior(path)

    def test_filters_before_validating_non_high_confidence_rows(self):
        rows = [
            self.row(IR_class="low_confidence_IR", start="bad"),
            self.row(),
        ]
        result = load_ir_prior(self.write_prior(rows))
        self.assertEqual(result.input_rows, 2)
        self.assertEqual(result.high_confidence_rows, 1)

    def test_rejects_invalid_coordinates(self):
        path = self.write_prior([self.row(start="201", end="200")])
        with self.assertRaisesRegex(IRPriorError, "positive integers.*start <= end"):
            load_ir_prior(path)

    def test_rejects_invalid_strand(self):
        path = self.write_prior([self.row(strand=".")])
        with self.assertRaisesRegex(IRPriorError, "strand"):
            load_ir_prior(path)

    def test_rejects_prior_without_high_confidence_ir(self):
        path = self.write_prior([self.row(IR_class="other")])
        with self.assertRaisesRegex(IRPriorError, "no high-confidence IR"):
            load_ir_prior(path)

    def test_collapses_exact_duplicate_rows(self):
        row = self.row()
        result = load_ir_prior(self.write_prior([row, dict(row)]))
        self.assertEqual(len(result.records), 1)
        self.assertEqual(result.duplicate_rows_collapsed, 1)

    def test_allows_coordinate_gene_conflicts_on_same_strand(self):
        rows = [self.row(), self.row(gene_id="G2", intron_id="I2")]
        result = load_ir_prior(self.write_prior(rows))
        self.assertEqual(len(result.records), 2)
        self.assertEqual(
            {record.key.gene_id for record in result.records}, {"G1", "G2"}
        )

    def test_rejects_coordinate_strand_conflicts(self):
        rows = [self.row(), self.row(strand="-", intron_id="I2")]
        with self.assertRaisesRegex(IRPriorError, "conflicting strands"):
            load_ir_prior(self.write_prior(rows))

    def test_matches_raw_gtf_then_binds_runtime_feature(self):
        loaded = load_ir_prior(self.write_prior([self.row()]))
        matched = match_ir_prior(
            loaded,
            self.write_gtf(),
            self.annotations(chromosome="Chr1"),
        )
        self.assertEqual(
            matched.runtime_introns,
            {IntronKey("G1", "Chr1", 100, 200, "+")},
        )
        self.assertEqual(matched.matched_unique_introns, 1)
        self.assertEqual(matched.matched_genes, 1)

    def test_does_not_normalize_chromosome_before_matching(self):
        loaded = load_ir_prior(
            self.write_prior([self.row(chromosome="Chr1")])
        )
        with self.assertRaisesRegex(IRPriorError, "matched 0"):
            match_ir_prior(
                loaded,
                self.write_gtf(chromosome="1"),
                self.annotations(chromosome="1"),
            )

    def test_duplicate_transcript_features_count_once(self):
        loaded = load_ir_prior(self.write_prior([self.row()]))
        matched = match_ir_prior(
            loaded,
            self.write_gtf(transcripts=("T1", "T2")),
            self.annotations(chromosome="Chr1", copies=2),
        )
        self.assertEqual(matched.matched_unique_introns, 1)
        self.assertEqual(len(matched.runtime_introns), 1)

    def test_partial_unmatched_records_include_reason(self):
        rows = [self.row(), self.row(chromosome="Chr2", intron_id="I2")]
        loaded = load_ir_prior(self.write_prior(rows))
        matched = match_ir_prior(
            loaded,
            self.write_gtf(),
            self.annotations(chromosome="Chr1"),
        )
        self.assertEqual(len(matched.unmatched), 1)
        self.assertEqual(matched.unmatched[0].reason, "chromosome_not_found")

    def test_gene_with_no_matching_intron_is_not_reported_as_missing(self):
        gtf_path = self.write_gtf()
        with open(gtf_path, "a", encoding="utf-8") as handle:
            handle.write(
                'Chr1\ttest\texon\t401\t600\t.\t+\t.\tgene_id "G2"; '
                'transcript_id "T2"; exon_number "1";\n'
            )
        rows = [
            self.row(),
            self.row(
                gene_id="G2",
                start="451",
                end="500",
                intron_id="I2",
            ),
        ]
        loaded = load_ir_prior(self.write_prior(rows))
        matched = match_ir_prior(
            loaded,
            gtf_path,
            self.annotations(chromosome="Chr1"),
        )
        self.assertEqual(len(matched.unmatched), 1)
        self.assertEqual(
            matched.unmatched[0].reason, "intron_coordinate_not_found"
        )

    def test_unmatched_unique_introns_collapse_metadata_variants(self):
        rows = [
            self.row(),
            self.row(chromosome="Chr2", intron_id="I2"),
            self.row(chromosome="Chr2", intron_id="I3"),
        ]
        loaded = load_ir_prior(self.write_prior(rows))
        matched = match_ir_prior(
            loaded,
            self.write_gtf(),
            self.annotations(chromosome="Chr1"),
        )
        self.assertEqual(len(matched.unmatched), 2)
        self.assertEqual(matched.unmatched_unique_introns, 1)


if __name__ == "__main__":
    unittest.main()
