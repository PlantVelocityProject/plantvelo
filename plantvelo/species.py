"""Data-driven discovery of bundled species IR priors."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Tuple

import data as bundled_data


SPECIES_NAME = re.compile(r"^[a-z0-9][a-z0-9_-]*$")


class SpeciesResourceError(ValueError):
    """Raised when a bundled species resource is invalid or unavailable."""


@dataclass(frozen=True)
class BuiltinSpecies:
    name: str
    filename: str
    path: str
    logical_path: str


@dataclass(frozen=True)
class PriorSource:
    path: str
    logical_path: str
    source_kind: str
    species: Optional[str]
    resource: Optional[BuiltinSpecies] = None
    loader: Optional[Callable] = None


def _resource_dir() -> Path:
    return Path(bundled_data.__file__).resolve().parent


def validate_species_name(name: str) -> str:
    if not isinstance(name, str) or SPECIES_NAME.fullmatch(name) is None:
        raise SpeciesResourceError(
            "invalid built-in species name: {}".format(name)
        )
    return name


def discover_species() -> Tuple[str, ...]:
    names = []
    for path in sorted(_resource_dir().glob("*.tsv")):
        if not path.is_file():
            continue
        names.append(validate_species_name(path.stem))
    if len(names) != len(set(names)):
        raise SpeciesResourceError("duplicate built-in species resource")
    return tuple(names)


def resolve_builtin_species(name: str) -> BuiltinSpecies:
    name = validate_species_name(name)
    if name not in discover_species():
        raise SpeciesResourceError("unknown built-in species: {}".format(name))
    filename = "{}.tsv".format(name)
    path = _resource_dir() / filename
    return BuiltinSpecies(
        name=name,
        filename=filename,
        path=str(path),
        logical_path="builtin://{}".format(filename),
    )


def resolve_prior_source(species: Optional[str], ir_prior: Optional[str]) -> PriorSource:
    """Resolve either the selected bundled species or a custom prior path."""

    if species is not None and ir_prior is not None:
        raise SpeciesResourceError("--species and --ir-prior are mutually exclusive")
    if species is not None:
        resource = resolve_builtin_species(species)
        from plantvelo.ir_prior import load_builtin_ir_prior

        return PriorSource(
            path=resource.path,
            logical_path=resource.logical_path,
            source_kind="builtin",
            species=resource.name,
            resource=resource,
            loader=load_builtin_ir_prior,
        )
    if ir_prior is not None:
        from plantvelo.ir_prior import load_ir_prior

        return PriorSource(
            path=ir_prior,
            logical_path=str(Path(ir_prior).resolve()),
            source_kind="custom",
            species=None,
            loader=load_ir_prior,
        )
    raise SpeciesResourceError("one of --species or --ir-prior is required")
