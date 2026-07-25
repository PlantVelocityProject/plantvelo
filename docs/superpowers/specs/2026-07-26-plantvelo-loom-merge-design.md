# PlantVelo Strict Loom Merge Design

## Purpose

PlantVelo must provide a supported way to combine multiple per-sample Loom
outputs into one standard Loom file for downstream analysis and publication.
The command must validate every input before merging, preserve all PlantVelo
layers and metadata that define counting semantics, align genes through the
unique `Accession` row attribute, and never publish a partial output.

The public interface is:

```bash
plantvelo merge \
    --output combined.loom \
    [--force] \
    sample1.loom sample2.loom sample3.loom
```

Only explicit input paths are accepted. Directory discovery and glob expansion
are left to the shell so that the exact input list remains visible and
reproducible.

## Architecture

The feature has separate command and domain layers:

- `plantvelo/commands/merge.py` defines the Click interface and translates
  domain failures into concise CLI errors.
- `plantvelo/merge.py` contains validation and merge behavior independent of
  Click. Its main public functions are `validate_loom_inputs()` and
  `merge_loom_files()`.
- `plantvelo/commands/plantvelo.py` registers `merge` after `run`, preserving
  natural help ordering.
- `plantvelo/commands/__init__.py` lazily exports the new command consistently
  with the existing `run` command.

The core implementation continues to depend on `loompy`; it does not introduce
AnnData, Scanpy, or direct HDF5 concatenation.

## CLI Contract

`plantvelo merge` requires at least two positional input files and one
`--output/-o` path. It rejects:

- missing or non-file inputs;
- the same input path more than once, including aliases that resolve to the
  same path;
- an output path that is also one of the inputs;
- an existing output unless `--force` is supplied;
- an output whose parent directory does not exist.

`--force` permits atomic replacement only. It never deletes or truncates the
existing output before the new temporary Loom has passed post-merge validation.

## Input Validation

Validation completes for all inputs before `loompy.combine` is called. Each
failure identifies the input file and the failed field or invariant.

### Loom specification

Each input must pass:

```python
loompy.LoomValidator().validate(path, strictness="speconly")
```

PlantVelo does not require Loom convention validation because its retained and
ambiguous layers and provenance attributes intentionally extend common RNA
velocity conventions.

### PlantVelo structure

Each input must contain exactly these named layers:

```text
spliced
unspliced
retained
ambiguous
```

The default matrix must be present. Required row attributes are `Accession`,
`Gene`, `Chromosome`, `Start`, `End`, and `Strand`; `CellID` is required as a
column attribute. All row and column attribute arrays must have the dimensions
required by the Loom specification.

The default matrix must equal the sum of the four mutually exclusive named
layers. This invariant is checked in column batches so validation never loads
an entire large Loom matrix into memory.

### Identity and compatibility

Within every input, `Accession` and `CellID` must each be unique. Across all
inputs, complete `CellID` values must also be globally unique. The merge command
does not infer a sample prefix from the filename because a valid Loom may have
been renamed or may use the standard PlantVelo ID generated without
`--sample-name`. Missing sample prefixes are nevertheless detected whenever
they cause duplicate complete CellIDs.

The current project data demonstrates why this rule is necessary: all 105,849
complete CellIDs are unique, while 1,368 bare Cell Ranger barcodes occur in
more than one sample.

All inputs must have the same unique `Accession` set. Row order may differ;
`loompy.combine(key="Accession")` performs the required alignment. After
alignment by `Accession`, every row attribute, including extra row attributes,
must have identical names and values across inputs. Column attribute name sets
must also be identical so all per-cell metadata can be concatenated without
discarding fields.

Named-layer sets and layer dtypes must match. These file attributes define
counting semantics and must be equal:

```text
plantvelo_version
classification_schema_version
classification_precedence
ir_prior_sha256
ir_prior_filter
matched_ir_introns
matched_ir_genes
gtf_sha256
```

After excluding non-semantic provenance, file-attribute name sets must be
identical and every value must match. The allowlist is `CreationDate`,
`gtf_path`, `ir_prior_path`, `merge_key`, `merge_input_count`, `merge_inputs`,
and `merge_created_at`. A missing, extra, or different attribute outside that
allowlist fails closed so newly introduced semantic metadata cannot be silently
ignored.

## Merge And Publication

After preflight validation, the core creates a private temporary directory
inside the output directory and calls:

```python
loompy.combine(
    files=input_paths,
    output_file=temporary_output,
    key="Accession",
    batch_size=1000,
)
```

Using the output directory ensures that final publication with `os.replace()`
stays on one filesystem. The combined file inherits the compatible PlantVelo
semantic metadata and sets `merge_key` to `Accession`, `merge_input_count` to
the number of inputs, `merge_inputs` to a JSON array of resolved input paths,
and `merge_created_at` to an ISO 8601 UTC timestamp. `CreationDate` describes
creation of the combined artifact rather than the first input.

The command does not infer or add a sample column. Existing column attributes
are retained and concatenated exactly; sample interpretation remains owned by
the producing workflow.

Before publication, the temporary result must pass the same Loom and PlantVelo
structural checks. It must additionally have:

- the validated input gene count;
- the sum of all input cell counts;
- the same globally unique CellID set as the inputs;
- the expected named layers and row/column attributes;
- a default matrix equal to the sum of the four named layers.

Only then does `os.replace()` publish the result. Any exception removes the
temporary directory and leaves a pre-existing output unchanged.

## Error Handling

The domain layer raises a dedicated validation exception containing actionable
context. The Click layer presents it as a normal command error without a Python
traceback. Representative messages include:

```text
BL-1h.loom: missing required layer 'retained'
BRZ-rep2.loom: ir_prior_sha256 differs from the first input
duplicate CellID across WT.loom and BL-0_5h.loom: AAAC...-1
output already exists; pass --force to replace it atomically
```

No validation failure is downgraded to a warning in strict mode.

## Testing

Tests use small real Loom fixtures created with `loompy.create` and the existing
`unittest` style. They cover:

- CLI help, registration, positional inputs, `--output`, and `--force`;
- successful merge with all four PlantVelo layers;
- correct alignment when one input has a different gene order;
- rejection of an invalid Loom file;
- rejection of missing or extra layers;
- rejection of duplicate or mismatched `Accession` values;
- rejection of row-attribute and semantic-provenance differences;
- rejection of duplicate CellIDs within one input and across inputs;
- preservation and concatenation of compatible column attributes;
- chunked default-matrix conservation checks;
- refusal to overwrite by default;
- atomic replacement with `--force`;
- preservation of an existing output when combine or post-validation fails;
- expected output dimensions, CellID set, layers, and merge provenance.

The complete existing test suite runs after the focused merge tests. README
documentation includes the command syntax, strict compatibility requirements,
shell glob example, overwrite behavior, and a compact Python API example.

## Non-Goals

This feature does not:

- discover Loom files recursively;
- merge gene unions or intersections;
- rename CellIDs or infer sample names;
- convert Loom to AnnData or another format;
- combine the per-run QC TSV files;
- reinterpret `retained` or `ambiguous` as standard scVelo layers;
- support a permissive compatibility mode.
