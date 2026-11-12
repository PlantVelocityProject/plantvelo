# PlantVelo

**Plant RNA velocity counting with high-confidence intron-retention priors**

Version 1.0.0 | Python >= 3.8 | BSD License

PlantVelo extends [velocyto.py](https://github.com/velocyto-team/velocyto.py)
with a retained state supported by high-confidence intron-retention (IR)
priors. It reads single-cell RNA-seq BAM files and a matching GTF annotation
to produce Loom 3.0 count matrices.

## Features

- **Prior-aware counting:** spliced (S), unspliced (U), retained (R), and ambiguous (A) counts.
- **Velocyto baseline:** original velocyto Default S/U/A counting for comparison.
- **Long-read evidence:** per-intron retention and splicing counts from a long-read BAM.
- **QC and merging:** prior-matching and classification reports, plus merging of compatible Loom files.

## Installation

```bash
conda create -n plantvelo python=3.10
conda activate plantvelo
conda install bioconda::velocyto.py
python -m pip install git+https://github.com/PlantVelocityProject/plantvelo.git
plantvelo --version
```

Loom output and merging require `loompy >= 3.0`.

## Usage

### Prior-aware counting

Prior-aware mode is the default and requires exactly one of `--species` or
`--ir-prior`. For example, use the built-in rice prior:

```bash
plantvelo run \
    --species osa \
    --bcfile filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outputfolder plantvelo_output \
    --sample-name sample01 \
    sample.bam annotation.gtf
```

For a custom prior, replace `--species osa` with `--ir-prior intron_prior.tsv`.
Barcodes in `--bcfile` must match the BAM `CB` tags.
Add `--mask repeats.gtf` to mask repeat regions.

### Velocyto baseline

```bash
plantvelo run \
    --ir-mode off \
    --bcfile filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outputfolder velocyto_baseline \
    --sample-name sample01 \
    sample.bam annotation.gtf
```

Baseline mode does not accept prior options (`--species`, `--ir-prior`, or
`--prior-dir`). See `plantvelo run --help` for all counting options.

## IR Priors

[PlantIRPrior](https://plantvelocityproject.github.io/PlantIRprior/)
provides species-specific, high-confidence IR priors derived from long-read
transcriptomic evidence. Visit the website to browse IR events, inspect
source samples and reference versions, and download versioned priors.

To use a downloaded prior, select the **PlantVelo prior** TSV matching your
species, reference genome assembly and annotation, then supply it with
`--ir-prior path/to/high_confidence_IR.tsv`.

### Built-in priors

These priors are included for offline use:

| Code | Species | Reference / annotation |
|---|---|---|
| `ath` | *Arabidopsis thaliana* | TAIR10 / Araport11 |
| `osa` | *Oryza sativa* (rice) | IRGSP-1.0 / RAP-DB |
| `sly` | *Solanum lycopersicum* (tomato) | SL4.0 / ITAG4.0 |
| `zma` | *Zea mays* (maize) | Zm-B73-REFERENCE-NAM-5.0 / Zm00001eb.1 |
| `gmx` | *Glycine max* (soybean) | Glycine_max_v2.1 / Wm82.a2.v1 |

### Custom priors

`--ir-prior` accepts CSV or TSV files, optionally gzip-compressed, containing
exactly five columns:

```text
gene_id    chromosome    start    end    strand
```

All input records are treated as high-confidence IR. Coordinates must be
**1-based closed intron intervals**. Gene IDs, chromosome names, intron coordinates, and strands must match the GTF exactly.

## Long-Read IR Evidence

Extract evidence from a coordinate-sorted, indexed long-read BAM and matching
GTF. The output directory must already exist.

```bash
plantvelo irprior \
    --gtf annotation.gtf \
    --bam sample.sorted.bam \
    --sample-id sample01 \
    --protocol nanopore \
    --species ath \
    --out sample01.irprior.tsv
```

The TSV contains `N_retained`, `N_spliced`, `N_specific`, and `N_informative`
counts per intron. The accompanying `<out>.run.json` records parameters,
versions, and run statistics. To turn this raw evidence into a custom prior,
select high-confidence IR introns and export the five-column prior format.

See `plantvelo irprior --help` for evidence thresholds and strand filtering
(default: `unstranded`).

## Output

Prior-aware counting produces:

```text
<output>/<sample>.loom
<output>/<sample>.ir_prior_qc.tsv
<output>/<sample>.ir_prior_unmatched.tsv
<output>/<sample>.classification_qc.tsv
```

The Loom layers are `spliced`, `unspliced`, `retained`, and `ambiguous`;
the default matrix is their sum. By default, molecules are grouped by cell
barcode, UMI, and gene, then classified with priority `U > R > S`.
Only introns in the high-confidence prior contribute retained evidence.
Unresolved states are ambiguous; molecules without a unique gene assignment
are excluded.

Baseline mode writes only the Loom file, with `spliced`, `unspliced`, and
`ambiguous` layers.

## Merge Loom Files

```bash
plantvelo merge \
    --output combined.plantvelo.loom \
    sample1.loom sample2.loom sample3.loom
```

Inputs must share the classification schema, GTF/prior hashes, gene set,
and compatible counting metadata and attributes. Prior-aware and baseline
outputs cannot be mixed; cell IDs must be globally unique. Use a distinct
`--sample-name` for each sample when counting.

QC TSV files are not merged. Add `--force` to replace an existing output
after validation.

## Support

Open a [GitHub issue](https://github.com/PlantVelocityProject/plantvelo/issues) for bugs or feature requests, or contact [jdluttzxr@stu.xmu.edu.cn](mailto:jdluttzxr@stu.xmu.edu.cn).
Pull requests are welcome.

## Citation

Please cite the velocyto publication:

> La Manno G, Soldatov R, Zeisel A, et al. RNA velocity of single cells.
> Nature. 2018;560:494–498.

If you use PlantVelocity for downstream analysis, please also cite:

> [Authors]. [PlantVelocity article title]. [Journal]. [Year]. DOI: [DOI].

Citation details will be added upon publication.
