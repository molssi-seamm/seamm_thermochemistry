"""Shared atomic reference-energy database and formation-energy arithmetic for SEAMM.

See ``~/Sites/reference-energy/2026-07-24_reference-energy/`` for the design
rationale. In short: a raw total energy from Gaussian/ORCA/Psi4/VASP has an
arbitrary, code-dependent zero that means nothing to a non-expert user and
cannot be compared across codes. This package holds the one shared table of
atomic reference energies, and the arithmetic to turn a system's electronic
energy into a physically meaningful energy or enthalpy of formation.
"""

from .db import ThermoDB, DEFAULT_DB_PATH  # noqa: F401
from .formation import (  # noqa: F401
    atomization_energy,
    formation_energy,
    formation_enthalpy,
    formation_gibbs_energy,
    MissingReferenceData,
)
from .derive import derive_dfH0_0K  # noqa: F401
from .report import format_report  # noqa: F401
from ._version import __version__  # noqa: F401

__all__ = [
    "ThermoDB",
    "DEFAULT_DB_PATH",
    "atomization_energy",
    "formation_energy",
    "formation_enthalpy",
    "formation_gibbs_energy",
    "MissingReferenceData",
    "derive_dfH0_0K",
    "format_report",
    "__version__",
]
