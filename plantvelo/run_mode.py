"""Runtime selection for prior-aware and velocyto Default counting."""
from __future__ import annotations

from dataclasses import dataclass
from functools import partial
from typing import Callable, Optional, Set

from plantvelo.classification import IntronKey
from plantvelo.loom_schema import (
    PRIOR_AWARE_PROFILE,
    VELOCYTO_DEFAULT_PROFILE,
    LoomSchemaProfile,
)


class RunModeError(ValueError):
    """Raised when the selected counting mode is invalid or incompatible."""


@dataclass(frozen=True)
class RunModeProfile:
    name: str
    logic_factory: Callable
    schema: LoomSchemaProfile
    requires_ir_prior: bool
    writes_prior_qc: bool


def resolve_run_mode(
    ir_mode: str,
    ir_prior: Optional[str],
    ir_registry: Set[IntronKey],
    *,
    species: Optional[str] = None,
    velocyto_module=None,
    prior_logic_class=None,
) -> RunModeProfile:
    if ir_mode == "prior-aware":
        if ir_prior is not None and species is not None:
            raise RunModeError("--species and --ir-prior are mutually exclusive")
        if ir_prior is None and species is None:
            raise RunModeError(
                "one of --species or --ir-prior is required"
            )
        if prior_logic_class is None:
            from plantvelo.logic import PriorAwareLogic

            prior_logic_class = PriorAwareLogic
        return RunModeProfile(
            name="prior-aware",
            logic_factory=partial(prior_logic_class, ir_registry),
            schema=PRIOR_AWARE_PROFILE,
            requires_ir_prior=True,
            writes_prior_qc=True,
        )

    if ir_mode == "off":
        if ir_prior is not None or species is not None:
            raise RunModeError(
                "--species and --ir-prior must not be provided when --ir-mode=off"
            )
        if velocyto_module is None:
            import velocyto as velocyto_module
        logic_factory = getattr(velocyto_module, "Default", None)
        if logic_factory is None:
            raise RunModeError("velocyto.Default is not available")
        return RunModeProfile(
            name="off",
            logic_factory=logic_factory,
            schema=VELOCYTO_DEFAULT_PROFILE,
            requires_ir_prior=False,
            writes_prior_qc=False,
        )

    raise RunModeError("unsupported IR mode: {}".format(ir_mode))


def validate_logic_contract(profile: RunModeProfile, logic_obj) -> None:
    actual_layers = tuple(logic_obj.layers)
    if actual_layers != profile.schema.layers:
        raise RunModeError(
            "logic '{}' has unexpected layers: {}".format(
                getattr(logic_obj, "name", type(logic_obj).__name__),
                ", ".join(actual_layers),
            )
        )
