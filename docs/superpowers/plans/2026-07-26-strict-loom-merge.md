# Strict Loom Merge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `plantvelo merge` command that strictly validates explicit PlantVelo Loom inputs, combines them with `loompy.combine(key="Accession")`, and atomically publishes a revalidated standard Loom output.

**Architecture:** Put Click parsing in `plantvelo/commands/merge.py` and all reusable validation/merge behavior in `plantvelo/merge.py`. The core validates Loom structure, PlantVelo semantics, gene identity, complete CellID uniqueness, and layer conservation in bounded batches before merging; it repeats output validation before `os.replace()`.

**Tech Stack:** Python 3.7+, Click, NumPy, loompy 3.x, standard-library `dataclasses`, `datetime`, `json`, `pathlib`, `tempfile`, `unittest`.

---

## File Map

- Create `plantvelo/merge.py`: domain exception, normalized Loom inspection, cross-file compatibility, chunked conservation, atomic combine, post-validation.
- Create `plantvelo/commands/merge.py`: Click command and error presentation.
- Create `tests/test_merge.py`: real Loom fixtures and core behavior tests.
- Create `tests/test_merge_cli.py`: command help, forwarding, errors, and top-level registration.
- Modify `plantvelo/commands/plantvelo.py`: register `merge` after `run`.
- Modify `plantvelo/commands/__init__.py`: lazy export of `merge`.
- Modify `README.md`: CLI, strict rules, overwrite behavior, and Python API.

### Task 1: Build Real Loom Test Fixtures And Establish Validation API

**Files:**
- Create: `tests/test_merge.py`
- Create: `plantvelo/merge.py`

- [ ] **Step 1: Write failing tests for one valid file and structural failures**

Create a fixture helper using `loompy.create` with two genes, configurable CellIDs, all four named layers, the default layer sum, required row attributes, and the semantic attributes from the design. Add tests that call `validate_loom_inputs()` and assert a valid summary, then reject a plain-text file, a missing layer, duplicate `Accession`, duplicate within-file `CellID`, and a corrupted default matrix.

The intended public API is:

```python
summary = validate_loom_inputs([first, second], batch_size=1)
assert summary.gene_count == 2
assert summary.cell_count == 2
assert summary.cell_ids == frozenset({"s1_AAAC-1", "s2_AAAC-1"})
```

Failures use:

```python
with self.assertRaisesRegex(MergeValidationError, "missing required layer"):
    validate_loom_inputs([bad, good], batch_size=1)
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run:

```bash
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/python \
  -m unittest -v tests.test_merge
```

Expected: import failure for `plantvelo.merge` because the production module does not exist.

- [ ] **Step 3: Implement minimal single-file inspection and summary types**

Create `plantvelo/merge.py` with:

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import FrozenSet, Mapping, Sequence, Tuple

import numpy as np

REQUIRED_LAYERS = ("spliced", "unspliced", "retained", "ambiguous")
REQUIRED_ROW_ATTRS = ("Accession", "Gene", "Chromosome", "Start", "End", "Strand")
REQUIRED_SEMANTIC_ATTRS = (
    "plantvelo_version",
    "classification_schema_version",
    "classification_precedence",
    "ir_prior_sha256",
    "ir_prior_filter",
    "matched_ir_introns",
    "matched_ir_genes",
    "gtf_sha256",
)
PROVENANCE_ATTRS = frozenset({
    "CreationDate", "gtf_path", "ir_prior_path", "merge_key",
    "merge_input_count", "merge_inputs", "merge_created_at",
})

class MergeValidationError(ValueError):
    pass

@dataclass(frozen=True)
class MergeSummary:
    input_paths: Tuple[str, ...]
    gene_count: int
    cell_count: int
    accessions: Tuple[str, ...]
    cell_ids: FrozenSet[str]
    row_attrs: Tuple[str, ...]
    col_attrs: Tuple[str, ...]
    layer_dtypes: Tuple[Tuple[str, str], ...]
    file_attrs: Mapping[str, object]
```

Implement byte/string normalization, `LoomValidator(strictness="speconly")`, exact layer checks, required attributes, uniqueness, and `_validate_layer_conservation()` using `range(0, n_cells, batch_size)` and slices. Never read an entire matrix in one call.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run the Task 1 command. Expected: all initial validation tests pass.

- [ ] **Step 5: Commit Task 1**

```bash
git add plantvelo/merge.py tests/test_merge.py
git commit -m "feat: validate PlantVelo loom inputs"
```

### Task 2: Enforce Strict Cross-File Compatibility

**Files:**
- Modify: `tests/test_merge.py`
- Modify: `plantvelo/merge.py`

- [ ] **Step 1: Add failing compatibility tests**

Add distinct tests for:

```text
duplicate complete CellID across files
mismatched Accession set
different input gene order with otherwise identical row metadata
mismatched row-attribute value after Accession alignment
mismatched row-attribute name set
mismatched column-attribute name set
mismatched named-layer dtype
mismatched required semantic attribute
extra unknown file attribute in only one input
different CreationDate/GTF path/IR path being accepted
```

For the reordered-gene case, assert the validation summary keeps the first
input order and succeeds. For every mismatch, assert that the exception names
the second file and field.

- [ ] **Step 2: Run focused tests and verify RED**

Run the Task 1 command. Expected: new compatibility tests fail because the
first implementation only validates files independently.

- [ ] **Step 3: Implement compatibility comparison keyed by Accession**

Store the first input's normalized row attributes, column-attribute names,
layer dtypes, and non-provenance file attributes as the reference. For each
later input:

```python
index_by_accession = {value: index for index, value in enumerate(accessions)}
indexer = np.asarray([index_by_accession[value] for value in reference_accessions])
aligned = current_values[indexer]
```

Require `np.array_equal(aligned, reference_values)` for every row attribute.
Compare file-attribute name sets and values only after removing
`PROVENANCE_ATTRS`. Track the owner of every complete CellID and fail on a
second owner. Return one `MergeSummary` containing the total cells and global
CellID set.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run the Task 1 command. Expected: every validation and compatibility test passes.

- [ ] **Step 5: Commit Task 2**

```bash
git add plantvelo/merge.py tests/test_merge.py
git commit -m "feat: enforce strict loom compatibility"
```

### Task 3: Add Atomic Combine And Post-Validation

**Files:**
- Modify: `tests/test_merge.py`
- Modify: `plantvelo/merge.py`

- [ ] **Step 1: Add failing end-to-end merge tests**

Use real fixtures to test that `merge_loom_files()`:

- invokes real `loompy.combine` behavior and aligns a reversed second gene axis;
- produces `(reference_genes, sum_cells)` and exactly the global CellID set;
- preserves all four layers and compatible column attributes;
- writes `merge_key`, `merge_input_count`, JSON `merge_inputs`,
  `merge_created_at`, and a new `CreationDate`;
- rejects fewer than two inputs, duplicate resolved input paths, output-as-input,
  a missing output parent, and an existing output without `force=True`;
- replaces an existing valid output with `force=True`;
- preserves an existing output if `loompy.combine` is patched to raise;
- preserves an existing output if post-validation is patched to raise.

- [ ] **Step 2: Run focused tests and verify RED**

Run the Task 1 command. Expected: failures because `merge_loom_files()` is absent.

- [ ] **Step 3: Implement path validation and atomic combine**

Implement:

```python
def merge_loom_files(
    input_paths: Sequence[str],
    output_path: str,
    *,
    force: bool = False,
    batch_size: int = 1000,
) -> MergeSummary:
```

Resolve all paths, enforce the CLI-independent path contract, and validate all
inputs. Create `TemporaryDirectory(prefix=".plantvelo-merge-", dir=parent)`;
write the temporary Loom inside it using:

```python
loompy.combine(
    files=list(summary.input_paths),
    output_file=str(temporary_output),
    key="Accession",
    file_attrs=merge_attrs,
    batch_size=batch_size,
)
```

Reinspect the single temporary output, assert dimensions, Accessions, CellID
set, attributes, and layer conservation against the preflight summary, then
call `os.replace(str(temporary_output), str(output))`. Context cleanup removes
all partial files. Build timestamps with timezone-aware UTC and serialize input
paths with `json.dumps()`.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run the Task 1 command. Expected: all core merge tests pass.

- [ ] **Step 5: Commit Task 3**

```bash
git add plantvelo/merge.py tests/test_merge.py
git commit -m "feat: combine loom files atomically"
```

### Task 4: Expose The Click Command

**Files:**
- Create: `tests/test_merge_cli.py`
- Create: `plantvelo/commands/merge.py`
- Modify: `plantvelo/commands/plantvelo.py`
- Modify: `plantvelo/commands/__init__.py`

- [ ] **Step 1: Write failing CLI tests**

Use `CliRunner` and `unittest.mock.patch` to assert:

```text
plantvelo merge --help lists INPUTS, --output, and --force
one input produces a Click usage error
two inputs and an output are forwarded to merge_loom_files
force=True is forwarded only when --force is present
MergeValidationError becomes Error: ... without traceback
top-level plantvelo --help lists run before merge
```

- [ ] **Step 2: Run CLI tests and verify RED**

```bash
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/python \
  -m unittest -v tests.test_merge_cli
```

Expected: import/registration failures because the command is absent.

- [ ] **Step 3: Implement and register the command**

Create a Click command with `@click.argument("loomfiles", nargs=-1,
type=click.Path(exists=True, dir_okay=False, readable=True))`, required
`--output/-o`, and boolean `--force`. Raise `click.UsageError` when fewer than
two files are supplied. Catch only `MergeValidationError` and raise
`click.ClickException(str(error))`. On success print the output path and merged
gene/cell counts.

Import and add the command after `run` in `commands/plantvelo.py`; add `merge`
to the lazy export allowlist and `__all__` in `commands/__init__.py`.

- [ ] **Step 4: Run CLI tests and full tests**

Run the Task 4 command, then:

```bash
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/python \
  -m unittest discover -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit Task 4**

```bash
git add plantvelo/commands/merge.py plantvelo/commands/plantvelo.py \
  plantvelo/commands/__init__.py tests/test_merge_cli.py
git commit -m "feat: add plantvelo merge command"
```

### Task 5: Document And Verify The Release Surface

**Files:**
- Modify: `README.md`
- Modify: `tests/test_merge.py` only if an actual integration defect requires a regression test

- [ ] **Step 1: Add README usage and guarantees**

Document:

```bash
plantvelo merge -o combined.loom sample1.loom sample2.loom
plantvelo merge -o combined.loom --force 22.plantvelo/*/*.loom
```

Explain strict semantic/hash/layer/gene checks, global complete CellID
uniqueness, `Accession` alignment, atomic overwrite behavior, output provenance,
and the equivalent Python call:

```python
from plantvelo.merge import merge_loom_files
merge_loom_files(inputs, "combined.loom", force=False)
```

- [ ] **Step 2: Run package and command smoke checks**

```bash
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/plantvelo --help
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/plantvelo merge --help
```

Expected: `run` and `merge` are listed, and merge options match the contract.

- [ ] **Step 3: Run focused and complete verification**

```bash
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/python \
  -m unittest -v tests.test_merge tests.test_merge_cli
MPLCONFIGDIR=/tmp/matplotlib \
  /home/Zhuxiangrong/software/anaconda3/envs/plantvelo/bin/python \
  -m unittest discover -v
git diff --check
```

Expected: all tests pass and `git diff --check` reports no errors.

- [ ] **Step 4: Commit documentation**

```bash
git add README.md
git commit -m "docs: explain strict loom merging"
```

- [ ] **Step 5: Review final scope**

Confirm the diff changes only merge core/CLI/tests/docs plus the approved plan,
does not alter `run` counting semantics, does not add dependencies, and does
not modify or commit unrelated pre-existing worktree changes.
