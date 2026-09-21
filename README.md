# PlantVelo

**Plant RNA velocity counting with high-confidence intron-retention priors**

> Version 0.3.0 | Python >= 3.7 | BSD License

PlantVelo is a Python command-line tool for generating RNA velocity count
matrices from single-cell RNA-seq BAM files. It extends
[velocyto.py](https://github.com/velocyto-team/velocyto.py) with an optional
high-confidence intron-retention (IR) state for plant transcriptomes.

---

## Overview

PlantVelo supports two counting modes:

- **Prior-aware**: produces spliced (S), unspliced (U), retained (R), and
  ambiguous (A) layers using a high-confidence IR prior.
- **Velocyto baseline**: produces the original velocyto Default spliced,
  unspliced, and ambiguous layers for comparison.

Counts are written as Loom 3.0 files for downstream analysis with
[PlantVelocity](https://github.com/plantvelocity/PlantVelocity) or other
RNA velocity workflows.

## Key Features

- High-confidence IR-aware molecule classification
- UMI-aware counting from BAM and GTF files
- Standard velocyto counting mode for benchmarking
- QC reports for prior matching and molecule classification
- Strict merging of compatible PlantVelo loom files

## Installation

Create a conda environment and install velocyto first:

```bash
conda create -n plantvelo python=3.10
conda activate plantvelo
conda install -c bioconda velocyto.py
pip install git+https://github.com/PlantVelocityProject/plantvelo.git
```

Verify the installation:

```bash
plantvelo --version
```

PlantVelo requires `loompy >= 3.0` to create and merge Loom 3.0 files.

## IR Prior

Prior-aware mode (the default) requires exactly one of `--species` or
`--ir-prior`. These options are mutually exclusive and are both prohibited
with `--ir-mode off`.

### Built-in priors

Built-in priors are installed with the package and are available offline:

| Code | Species |
|---|---|
| `ath` | Arabidopsis thaliana |
| `osa` | Oryza sativa (rice) |
| `sly` | Solanum lycopersicum (tomato) |
| `zma` | Zea mays (maize) |
| `gmx` | Glycine max (soybean) |

Packaged `data/<species>.tsv` files use a five-column format:
`gene_id`, `chromosome`, `start`, `end`, and `strand`.
All records are treated as high-confidence IR. Available species are
discovered from the TSV filenames.

### Custom priors

The `--ir-prior` option remains supported and is not deprecated.
Custom files use a six-column format and may be CSV or TSV, optionally
gzip-compressed (`.csv`, `.tsv`, `.csv.gz`, or `.tsv.gz`):

```text
gene_id
chromosome
start
end
strand
IR_class
```

For custom priors, only rows with `IR_class == "high_confidence_IR"` are used. Coordinates must
be 1-based closed intervals for both built-in and custom priors and must match the GTF gene IDs, chromosome names,
coordinates, and strands.

## Usage

### Built-in prior (rice)

```bash
plantvelo run \
    --ir-mode prior-aware \
    --species osa \
    --bcfile filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outputfolder plantvelo_output \
    --sample-name sample01 \
    --mask repeat_masker.gtf \
    sample.bam annotation.gtf
```

### Custom prior

```bash
plantvelo run \
    --ir-prior intron_prior.tsv \
    --bcfile filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outputfolder plantvelo_output \
    --sample-name sample01 \
    sample.bam annotation.gtf
```

### Velocyto baseline

```bash
plantvelo run \
    --ir-mode off \
    --bcfile filtered_feature_bc_matrix/barcodes.tsv.gz \
    --outputfolder velocyto_baseline \
    --sample-name sample01 \
    --mask repeats.gtf \
    sample.bam annotation.gtf
```

The repeat mask (`--mask`) is optional. Use `plantvelo run --help` to view all available counting options.

## Merge Loom Files

Merge two or more PlantVelo loom files with the same classification schema:

```bash
plantvelo merge \
    --output combined.plantvelo.loom \
    sample1.loom sample2.loom sample3.loom
```

Use `--force` to replace an existing output after validation. Prior-aware and
velocyto baseline loom files cannot be mixed in one merge.
Inputs must also have compatible GTF/prior hashes, counting metadata, genes,
and attributes. Use a distinct `--sample-name` for each sample to keep
cell IDs globally unique. QC TSV files are not merged.

## Output

Prior-aware mode produces:

```text
<output>/<sample>.loom
<output>/<sample>.ir_prior_qc.tsv
<output>/<sample>.ir_prior_unmatched.tsv
<output>/<sample>.classification_qc.tsv
```

The loom file contains `spliced`, `unspliced`, `retained`, and `ambiguous`
layers. Velocyto baseline mode produces only `<output>/<sample>.loom`, with
`spliced`, `unspliced`, and `ambiguous` layers.

## Citation

PlantVelo builds on velocyto:

> La Manno G, Soldatov R, Zeisel A, et al. RNA velocity of single cells.
> Nature. 2018;560:494-498.

If you use PlantVelo with PlantVelocity, please also cite the corresponding
PlantVelocity publication when available.

## Contributing / Issues

- **Bug reports and feature requests:**
  [PlantVelo issues](https://github.com/PlantVelocityProject/plantvelo/issues)
- **Contact:** [jdluttzxr@stu.xmu.edu.cn](mailto:jdluttzxr@stu.xmu.edu.cn)
- **Pull requests:** Contributions are welcome.
