#!/usr/bin/env python3

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import warnings
import pysam

if __package__:
    from ._version import __version__
else:
    from _version import __version__

STRAND_MODES = ("unstranded", "tag", "flag", "tag_then_flag")


@dataclass
class Intron:
    gene_id: str
    intron_id: str
    chrom: str
    start0: int
    end0: int
    strand: str
    up_start0: int = None
    up_end0: int = None
    down_start0: int = None
    down_end0: int = None
    anchors: list = field(default_factory=list)
    exon_overlap_warning: str = "no"

    def __post_init__(self):
        if not self.anchors:
            self.anchors = [
                (
                    ".",
                    self.up_start0,
                    self.up_end0,
                    self.down_start0,
                    self.down_end0,
                )
            ]

        first_anchor = self.anchors[0]
        if self.up_start0 is None:
            self.up_start0 = first_anchor[1]
        if self.up_end0 is None:
            self.up_end0 = first_anchor[2]
        if self.down_start0 is None:
            self.down_start0 = first_anchor[3]
        if self.down_end0 is None:
            self.down_end0 = first_anchor[4]

    @property
    def supporting_transcripts(self):
        return ",".join(sorted({anchor[0] for anchor in self.anchors}))

    @property
    def anchor_count(self):
        return len(self.anchors)

    @property
    def intron_uid(self):
        return f"{self.gene_id}|{self.chrom}:{self.start0 + 1}-{self.end0}:{self.strand}"


@dataclass
class ReadEvidence:
    blocks: list
    splice_intervals: list


def parse_gtf_attributes(attr_string):
    """
    Parse GTF/GFF-like attributes.
    Mainly designed for GTF:
        gene_id "AT1G01010"; transcript_id "AT1G01010.1";
    """
    attrs = {}
    attr_string = attr_string.strip().rstrip(";")

    for part in attr_string.split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part and " " not in part.split("=")[0]:
            key, value = part.split("=", 1)
        else:
            items = part.split(None, 1)
            if len(items) != 2:
                continue
            key, value = items

        value = value.strip().strip('"')
        attrs[key] = value

    return attrs


def load_introns_from_gtf(gtf_file):
    """
    Build introns and their transcript-specific exon anchors from GTF exons.

    GTF coordinates:
        1-based inclusive

    Internal coordinates:
        0-based half-open

    Return:
        List of introns.
    """
    exons_by_tx = defaultdict(list)
    exons_by_gene = defaultdict(list)

    with open(gtf_file) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields

            if feature != "exon":
                continue

            start0 = int(start) - 1
            end0 = int(end)
            if start0 < 0 or end0 <= start0:
                raise ValueError(
                    f"Invalid exon coordinates in {gtf_file}:{line_number}"
                )

            attr_dict = parse_gtf_attributes(attrs)

            gene_id = attr_dict.get("gene_id")
            tx_id = attr_dict.get("transcript_id")

            if gene_id is None or tx_id is None:
                continue

            exons_by_tx[tx_id].append(
                (chrom, start0, end0, strand, gene_id)
            )
            exons_by_gene[gene_id].append((chrom, start0, end0, strand, tx_id))

    intron_dict = {}

    for tx_id, exons in exons_by_tx.items():
        if len(exons) < 2:
            continue

        exons = sorted(exons, key=lambda x: x[1])

        for i in range(len(exons) - 1):
            chrom1, s1, e1, strand1, gene1 = exons[i]
            chrom2, s2, e2, strand2, gene2 = exons[i + 1]

            if chrom1 != chrom2 or strand1 != strand2 or gene1 != gene2:
                continue

            intron_start0 = e1
            intron_end0 = s2

            if intron_start0 >= intron_end0:
                continue

            key = (gene1, chrom1, intron_start0, intron_end0, strand1)

            if key not in intron_dict:
                intron_dict[key] = {
                    "gene_id": gene1,
                    "chrom": chrom1,
                    "start0": intron_start0,
                    "end0": intron_end0,
                    "strand": strand1,
                    "anchors": [],
                }
            intron_dict[key]["anchors"].append((tx_id, s1, e1, s2, e2))

    for item in intron_dict.values():
        item["exon_overlap_warning"] = "no"
        for exon_chrom, exon_start0, exon_end0, exon_strand, tx_id in exons_by_gene[
            item["gene_id"]
        ]:
            if exon_chrom != item["chrom"] or exon_strand != item["strand"]:
                continue
            if overlap_len(item["start0"], item["end0"], exon_start0, exon_end0) > 0:
                item["exon_overlap_warning"] = "yes"
                break

    introns_by_gene = defaultdict(list)

    for key, item in intron_dict.items():
        introns_by_gene[item["gene_id"]].append(item)

    introns = []
    for gene_id, items in introns_by_gene.items():
        strand = items[0]["strand"]

        if strand == "-":
            items = sorted(items, key=lambda x: x["start0"], reverse=True)
        else:
            items = sorted(items, key=lambda x: x["start0"])

        for idx, item in enumerate(items, start=1):
            intron_id = f"{gene_id}.intron{idx}"
            introns.append(
                Intron(
                    gene_id=item["gene_id"],
                    intron_id=intron_id,
                    chrom=item["chrom"],
                    start0=item["start0"],
                    end0=item["end0"],
                    strand=item["strand"],
                    anchors=item["anchors"],
                    exon_overlap_warning=item["exon_overlap_warning"],
                )
            )

    return introns


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def get_splice_intervals(read):
    """
    Return N-cigar intervals as 0-based half-open reference intervals.
    """
    intervals = []

    if read.cigartuples is None:
        return intervals

    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        # CIGAR ops:
        # M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8
        if op == 3:
            intervals.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (0, 2, 7, 8):
            ref_pos += length
        else:
            pass

    return intervals


def get_read_evidence(read):
    # pysam blocks cover M/=/X reference positions. D and N leave gaps,
    # while I can split adjacent blocks without consuming reference positions.
    return ReadEvidence(
        blocks=read.get_blocks(),
        splice_intervals=get_splice_intervals(read),
    )


def blocks_overlap(blocks, start0, end0):
    total = 0
    for b_start, b_end in blocks:
        total += overlap_len(b_start, b_end, start0, end0)
    return total


def validate_intron_coverage(value):
    if not 0 < value <= 1:
        raise ValueError("min_intron_coverage must be greater than 0 and at most 1")
    return value


def validate_integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def parse_nonnegative_integer(value):
    try:
        return validate_integer(int(value), "value")
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("must be an integer >= 0") from error


def parse_positive_integer(value):
    try:
        return validate_integer(int(value), "value", minimum=1)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("must be an integer > 0") from error


def validate_evidence_parameters(min_anchor, min_intron_coverage, splice_tolerance):
    validate_integer(min_anchor, "min_anchor", minimum=1)
    validate_intron_coverage(min_intron_coverage)
    validate_integer(splice_tolerance, "splice_tolerance")


def validate_count_parameters(min_mapq, min_anchor, min_intron_coverage,
                              splice_tolerance, strand_mode):
    validate_evidence_parameters(min_anchor, min_intron_coverage, splice_tolerance)
    validate_integer(min_mapq, "min_mapq")
    if strand_mode not in STRAND_MODES:
        raise ValueError(f"strand_mode must be one of: {', '.join(STRAND_MODES)}")


def parse_intron_coverage(value):
    try:
        return validate_intron_coverage(float(value))
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            "--min-intron-coverage must be a fraction greater than 0 and at most 1"
        ) from error


def classify_read_for_intron(
    read,
    intron,
    min_anchor=8,
    min_intron_coverage=0.90,
    splice_tolerance=5,
    evidence=None,
):
    """
    Classify a single read against a single intron.

    Retention requires both exon anchors, no CIGAR N overlapping the path
    between those anchors, and M/=/X coverage of at least min_intron_coverage
    of the intron. D does not count as coverage; I consumes no reference bases.

    Return one of:
        retained
        spliced
        conflict
        uninformative
    """
    validate_evidence_parameters(min_anchor, min_intron_coverage, splice_tolerance)
    if evidence is None:
        evidence = get_read_evidence(read)

    # 1. Spliced evidence: CIGAR N matches this intron.
    spliced = False
    for n_start, n_end in evidence.splice_intervals:
        if (
            abs(n_start - intron.start0) <= splice_tolerance
            and abs(n_end - intron.end0) <= splice_tolerance
            and blocks_overlap(evidence.blocks, n_start - min_anchor, n_start)
            >= min_anchor
            and blocks_overlap(evidence.blocks, n_end, n_end + min_anchor)
            >= min_anchor
        ):
            spliced = True
            break

    # 2. Retained evidence: sufficient coverage along an unspliced E-I-E path.
    intron_ov = blocks_overlap(evidence.blocks, intron.start0, intron.end0)
    intron_length = intron.end0 - intron.start0

    retained = False
    if intron_length > 0 and intron_ov / intron_length >= min_intron_coverage:
        for tx_id, up_start0, up_end0, down_start0, down_end0 in intron.anchors:
            up_anchor_start = max(up_start0, up_end0 - min_anchor)
            up_anchor_end = up_end0

            down_anchor_start = down_start0
            down_anchor_end = min(down_end0, down_start0 + min_anchor)

            if any(
                overlap_len(n_start, n_end, up_anchor_start, down_anchor_end) > 0
                for n_start, n_end in evidence.splice_intervals
            ):
                continue

            up_ov = blocks_overlap(evidence.blocks, up_anchor_start, up_anchor_end)
            down_ov = blocks_overlap(
                evidence.blocks, down_anchor_start, down_anchor_end
            )

            if up_ov >= min_anchor and down_ov >= min_anchor:
                retained = True
                break

    if retained and spliced:
        return "conflict"
    elif retained:
        return "retained"
    elif spliced:
        return "spliced"
    else:
        return "uninformative"


def normalize_protocol(x):
    original = x.strip()
    x = original.lower()
    if x in ["isoseq", "iso-seq", "pacbio", "pacbio_isoseq"]:
        return "Iso-seq"
    if x in ["ont", "nanopore", "ont_cdna", "ont_drna", "direct_rna"]:
        return "ONT"
    if x in ["flep", "flep-seq", "flepseq"]:
        return "FLEP-seq"
    return original


def pass_read_filter(read, min_mapq):
    if read.is_unmapped:
        return False
    if read.is_secondary:
        return False
    if read.is_supplementary:
        return False
    if read.mapping_quality < min_mapq:
        return False
    return True


def infer_transcript_strand(read, strand_mode):
    """Use ts/XS values as reference-strand labels, with optional flag fallback.

    This preserves existing inference. Tag meanings and flag orientation must
    be checked against the aligner and library protocol; no orientation is
    inferred automatically from the sequencing technology name.
    """
    if strand_mode not in STRAND_MODES:
        raise ValueError(f"strand_mode must be one of: {', '.join(STRAND_MODES)}")
    if strand_mode == "unstranded":
        return None
    if strand_mode in {"tag", "tag_then_flag"}:
        for tag in ("ts", "XS"):
            if read.has_tag(tag):
                strand = read.get_tag(tag)
                if strand in {"+", "-"}:
                    return strand
        if strand_mode == "tag":
            return None
    return "-" if read.is_reverse else "+"


def read_reference_span_overlaps(read, start0, end0):
    return read.reference_start < end0 and read.reference_end > start0


def retained_read_is_specific(
    read,
    intron,
    introns_in_gene,
    min_anchor=8,
    min_intron_coverage=0.90,
    splice_tolerance=5,
    evidence=None,
    _statuses=None,
):
    """
    Check specificity after the target intron has been classified as retained.

    Every other intron overlapping the read reference span must be spliced,
    including partially overlapping introns. Uncovered introns are ignored.
    A read overlapping no other introns passes; this does not establish
    specificity across the entire gene.
    """
    validate_evidence_parameters(min_anchor, min_intron_coverage, splice_tolerance)
    if evidence is None:
        evidence = get_read_evidence(read)

    for other in introns_in_gene:
        if other.intron_id == intron.intron_id:
            continue

        if not read_reference_span_overlaps(read, other.start0, other.end0):
            continue

        other_status = None if _statuses is None else _statuses.get(other.intron_id)
        if other_status is None:
            other_status = classify_read_for_intron(
                read,
                other,
                min_anchor=min_anchor,
                min_intron_coverage=min_intron_coverage,
                splice_tolerance=splice_tolerance,
                evidence=evidence,
            )

        if other_status != "spliced":
            return False

    return True


class IntronCounts(dict):
    """Per-intron counts, with a summary of input checks and processing units."""

    def __init__(self, introns, summary):
        super().__init__((intron.intron_id, {"retained": 0, "spliced": 0, "specific": 0})
                         for intron in introns)
        self.summary = summary


def check_bam_introns(bam, introns):
    """Return analyzable introns and skip metadata; fail on invalid inputs."""
    if not introns:
        raise ValueError("Annotation contains no valid introns")
    try:
        bam.check_index()
    except (ValueError, OSError) as error:
        raise ValueError("BAM index is missing or unreadable; index the coordinate-sorted BAM first") from error
    sort_order = bam.header.to_dict().get("HD", {}).get("SO")
    if sort_order not in (None, "unknown", "coordinate"):
        raise ValueError(f"BAM must be coordinate sorted (header SO={sort_order})")
    lengths = dict(zip(bam.references, bam.lengths))
    skipped = defaultdict(int)
    valid = []
    gene_locations = {}
    for intron in introns:
        if intron.chrom not in lengths:
            skipped[intron.chrom] += 1
            continue
        chrom_length = lengths[intron.chrom]
        if not 0 <= intron.start0 < intron.end0 <= chrom_length:
            raise ValueError(f"Invalid intron coordinates or reference bounds: {intron.intron_uid}")
        for tx_id, up_start, up_end, down_start, down_end in intron.anchors:
            if not (0 <= up_start < up_end == intron.start0
                    and intron.end0 == down_start < down_end <= chrom_length):
                raise ValueError(f"Invalid exon anchor coordinates or reference bounds: {intron.intron_uid}")
        location = (intron.chrom, intron.strand)
        previous = gene_locations.setdefault(intron.gene_id, location)
        if previous != location:
            raise ValueError(f"Gene {intron.gene_id} has introns on multiple reference sequences or strands")
        valid.append(intron)
    if not valid:
        raise ValueError("GTF introns and BAM have no common reference sequence names")
    skipped = dict(sorted(skipped.items()))
    messages = []
    if skipped:
        details = ", ".join(f"{chrom}: {count} introns" for chrom, count in skipped.items())
        message = f"Skipping GTF reference sequences absent from BAM ({details}); skipped introns are omitted from output"
        warnings.warn(message, UserWarning, stacklevel=2)
        messages.append(message)
    return valid, {
        "annotation": {"input_introns": len(introns), "analyzed_introns": len(valid),
                       "skipped_introns": sum(skipped.values()), "skipped_chromosomes": skipped},
        "warnings": messages,
    }


def count_introns(
    introns,
    bam_file,
    min_mapq=20,
    min_anchor=8,
    min_intron_coverage=0.90,
    splice_tolerance=5,
    strand_mode="unstranded",
):
    """Count each read name per intron, fetching each gene region once.

    Results contain zero counts for valid uncovered introns and omit introns
    on missing references. The attached summary labels all strand counts as
    read-intron pairs, never as unique reads across the sample.
    """
    validate_count_parameters(min_mapq, min_anchor, min_intron_coverage,
                              splice_tolerance, strand_mode)
    if not introns:
        raise ValueError("Annotation contains no valid introns")
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        valid, summary = check_bam_introns(bam, introns)
        stats = IntronCounts(valid, summary)
        strand_stats = {"mode": strand_mode, "unit": "read_intron_pairs",
                        "pass": 0, "failed": 0, "unknown": 0,
                        "unknown_annotation_strand": 0}
        processing = {"gene_regions": 0, "evidence_extractions": 0,
                      "intron_classifications": 0}
        summary["strand_filter"] = strand_stats
        summary["processing"] = processing
        introns_by_gene = defaultdict(list)
        for intron in valid:
            introns_by_gene[intron.gene_id].append(intron)

        for gene_introns in introns_by_gene.values():
            chrom = gene_introns[0].chrom
            # Preserve the original per-intron candidate windows while
            # fetching their enclosing gene region only once.
            windows = {i.intron_id: (max(0, i.start0 - 50), i.end0 + 50)
                       for i in gene_introns}
            fetch_start = min(start for start, end in windows.values())
            fetch_end = max(end for start, end in windows.values())
            names = {i.intron_id: {"retained": set(), "spliced": set(), "specific": set()}
                     for i in gene_introns}
            processing["gene_regions"] += 1
            for read in bam.fetch(chrom, fetch_start, fetch_end):
                if not pass_read_filter(read, min_mapq):
                    continue
                candidates = [i for i in gene_introns
                              if read_reference_span_overlaps(read, *windows[i.intron_id])]
                read_strand = infer_transcript_strand(read, strand_mode)
                accepted = []
                for intron in candidates:
                    if strand_mode != "unstranded":
                        if intron.strand not in {"+", "-"}:
                            strand_stats["unknown_annotation_strand"] += 1
                        elif read_strand is None:
                            strand_stats["unknown"] += 1
                            continue
                        elif read_strand != intron.strand:
                            strand_stats["failed"] += 1
                            continue
                        else:
                            strand_stats["pass"] += 1
                    accepted.append(intron)
                if not accepted:
                    continue
                evidence = get_read_evidence(read)
                processing["evidence_extractions"] += 1
                statuses = {}
                for intron in accepted:
                    statuses[intron.intron_id] = classify_read_for_intron(
                        read, intron, min_anchor=min_anchor,
                        min_intron_coverage=min_intron_coverage,
                        splice_tolerance=splice_tolerance, evidence=evidence,
                    )
                    processing["intron_classifications"] += 1
                for intron in accepted:
                    status = statuses[intron.intron_id]
                    if status not in {"retained", "spliced"}:
                        continue
                    names[intron.intron_id][status].add(read.query_name)
                    if status == "retained" and retained_read_is_specific(
                        read, intron, gene_introns, min_anchor=min_anchor,
                        min_intron_coverage=min_intron_coverage,
                        splice_tolerance=splice_tolerance,
                        evidence=evidence, _statuses=statuses,
                    ):
                        names[intron.intron_id]["specific"].add(read.query_name)
            for intron in gene_introns:
                groups = names[intron.intron_id]
                conflicts = groups["retained"] & groups["spliced"]
                stats[intron.intron_id] = {
                    key: len(read_names - conflicts) for key, read_names in groups.items()
                }
    if strand_mode != "unstranded":
        for key in ("pass", "failed", "unknown"):
            print(f"N_strand_{key}_pairs\t{strand_stats[key]}", file=sys.stderr)
    return stats


def IRprior(
    gtf,
    bam,
    sample_id,
    protocol,
    species,
    out,
    min_mapq=20,
    min_anchor=8,
    min_intron_coverage=0.90,
    splice_tolerance=5,
    strand_mode="unstranded",
    _load_introns_from_gtf=None,
    _count_introns=None,
):
    started_at = datetime.now(timezone.utc).isoformat()
    validate_count_parameters(min_mapq, min_anchor, min_intron_coverage,
                              splice_tolerance, strand_mode)
    for name, value in (("sample_id", sample_id), ("protocol", protocol), ("species", species)):
        if not isinstance(value, str) or not value.strip() or any(c in value for c in "\t\r\n"):
            raise ValueError(f"{name} must be a non-empty value without tabs or newlines")
    output_paths = (Path(out).resolve(), Path(str(out) + ".run.json").resolve())
    if any(path in (Path(gtf).resolve(), Path(bam).resolve()) for path in output_paths):
        raise ValueError("Output paths must not overwrite the GTF or BAM input")
    protocol = normalize_protocol(protocol)

    if _load_introns_from_gtf is None:
        _load_introns_from_gtf = load_introns_from_gtf
    if _count_introns is None:
        _count_introns = count_introns

    introns = _load_introns_from_gtf(gtf)
    if not introns:
        raise ValueError("Annotation contains no valid introns")

    stats = _count_introns(
        introns=introns,
        bam_file=bam,
        min_mapq=min_mapq,
        min_anchor=min_anchor,
        min_intron_coverage=min_intron_coverage,
        splice_tolerance=splice_tolerance,
        strand_mode=strand_mode,
    )

    for destination in (Path(out), Path(str(out) + ".run.json")):
        if destination.is_dir():
            raise IsADirectoryError(f"Output path is a directory: {destination}")
    # Prepare both files before replacing existing results. Parsing, counting,
    # or report-writing failures must not truncate a previous successful TSV.
    with tempfile.TemporaryDirectory(prefix=".plantvelo-irprior-", dir=Path(out).parent) as temporary:
        staged_tsv = Path(temporary) / "result.tsv"
        staged_report = Path(temporary) / "run.json"
        with open(staged_tsv, "w") as out_handle:
            header = [
                "sample_id",
                "gene_id",
                "intron_id",
                "intron_uid",
                "chromosome",
                "start",
                "end",
                "strand",
                "supporting_transcripts",
                "anchor_count",
                "exon_overlap_warning",
                "N_retained",
                "N_spliced",
                "N_specific",
                "N_informative",
            ]
            header.append("protocol_support")
            header.append("species")

            out_handle.write("\t".join(header) + "\n")

            for intron in introns:
                if intron.intron_id not in stats:
                    continue
                counts = stats[intron.intron_id]
                n_retained = counts["retained"]
                n_spliced = counts["spliced"]
                n_specific = counts["specific"]
                n_informative = n_retained + n_spliced

                row = [
                    sample_id,
                    intron.gene_id,
                    intron.intron_id,
                    intron.intron_uid,
                    intron.chrom,
                    str(intron.start0 + 1),
                    str(intron.end0),
                    intron.strand,
                    intron.supporting_transcripts,
                    str(intron.anchor_count),
                    intron.exon_overlap_warning,
                    str(n_retained),
                    str(n_spliced),
                    str(n_specific),
                    str(n_informative),
                ]
                row.append(protocol)
                row.append(species)

                out_handle.write("\t".join(row) + "\n")

        report = {
            "software": {"name": "plantvelo", "version": __version__,
                         "pysam_version": pysam.__version__, "python_version": platform.python_version()},
            "started_at": started_at,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "parameters": {"gtf": str(Path(gtf).resolve()), "bam": str(Path(bam).resolve()),
                           "out": str(Path(out).resolve()), "sample_id": sample_id,
                           "protocol": protocol, "species": species, "min_mapq": min_mapq,
                           "min_anchor": min_anchor, "min_intron_coverage": min_intron_coverage,
                           "splice_tolerance": splice_tolerance, "strand_mode": strand_mode},
            "counts": {"unit": "read_intron_pairs",
                       "retained": sum(c["retained"] for c in stats.values()),
                       "spliced": sum(c["spliced"] for c in stats.values()),
                       "specific": sum(c["specific"] for c in stats.values()),
                       "informative": sum(c["retained"] + c["spliced"] for c in stats.values())},
        }
        report.update(getattr(stats, "summary", {}))
        with open(staged_report, "w") as handle:
            json.dump(report, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        os.replace(staged_tsv, out)
        os.replace(staged_report, str(out) + ".run.json")


def add_irprior_arguments(parser):
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF file.")
    parser.add_argument("--bam", required=True, help="Input coordinate-sorted BAM file.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--protocol", required=True, help="Sequencing protocol name.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--out", required=True, help="Output TSV file.")

    parser.add_argument("--min-mapq", type=parse_nonnegative_integer, default=20,
                        help="Minimum mapping quality (default: %(default)s).")
    parser.add_argument("--min-anchor", type=parse_positive_integer, default=8,
                        help="Aligned reference bases required on each side for retained and spliced evidence (default: %(default)s).")
    parser.add_argument(
        "--min-intron-coverage",
        type=parse_intron_coverage,
        default=0.90,
        help=(
            "Minimum fraction of intron reference bases covered by M/=/X "
            "for retention (0 < fraction <= 1; default: %(default)s). "
            "Both exon anchors and no CIGAR N along their connecting path "
            "are also required; deletions do not count as coverage."
        ),
    )
    parser.add_argument("--splice-tolerance", type=parse_nonnegative_integer, default=5,
                        help="Allowed difference at each splice boundary in bp (default: %(default)s).")
    parser.add_argument(
        "--strand-mode",
        choices=STRAND_MODES,
        default="unstranded",
        help=(
            "unstranded disables strand filtering (default); tag uses ts/XS; "
            "flag uses the alignment direction; tag_then_flag uses tags with "
            "flag fallback. Tag and flag meanings must match the library/aligner."
        ),
    )
    return parser


def run_irprior_from_args(args):
    IRprior(
        gtf=args.gtf,
        bam=args.bam,
        sample_id=args.sample_id,
        protocol=args.protocol,
        species=args.species,
        out=args.out,
        min_mapq=args.min_mapq,
        min_anchor=args.min_anchor,
        min_intron_coverage=args.min_intron_coverage,
        splice_tolerance=args.splice_tolerance,
        strand_mode=args.strand_mode,
    )


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Calculate single-sample intron retention evidence from one long-read BAM file."
    )
    return add_irprior_arguments(parser)


def parse_and_run(parser, argv=None):
    args = parser.parse_args(argv)
    try:
        run_irprior_from_args(args)
    except (ValueError, OSError) as error:
        parser.error(str(error))


def main(argv=None):
    parse_and_run(build_arg_parser(), argv)


if __name__ == "__main__":
    main()
