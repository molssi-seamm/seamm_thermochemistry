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
from .atct import ATcT, ATcTEntry, ATcTError, AmbiguousEntry, MissingEntry  # noqa: F401
from .atct import reaction_uncertainty  # noqa: F401
from .fragment_store import FragmentStore, FragmentRecord, harvest_job  # noqa: F401
from .fragment_store import load_species_map  # noqa: F401
from .fragment_store import (  # noqa: F401
    check_structure,
    harvest_configuration,
    store_is_writable,
)
from .geometry_store import GeometryRecord, EnergyRecord, PairingError  # noqa: F401
from .geometry_store import split_model_chemistry  # noqa: F401
from .ladder import Ladder, build_ladder, LadderError  # noqa: F401
from .zpe import zpe_scale_factor  # noqa: F401
from .frequency_corrections import frequency_correction  # noqa: F401
from .fragment_store import official_store_path  # noqa: F401
from .layered_store import LayeredStore  # noqa: F401
from .store_tools import (  # noqa: F401
    build_official,
    export_local,
    layer_conflicts,
    merge_stores,
    migrate_store,
)

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
    "ATcT",
    "FragmentStore",
    "FragmentRecord",
    "harvest_job",
    "harvest_configuration",
    "check_structure",
    "store_is_writable",
    "load_species_map",
    "GeometryRecord",
    "EnergyRecord",
    "PairingError",
    "split_model_chemistry",
    "Ladder",
    "build_ladder",
    "LadderError",
    "zpe_scale_factor",
    "frequency_correction",
    "LayeredStore",
    "official_store_path",
    "build_official",
    "export_local",
    "layer_conflicts",
    "merge_stores",
    "migrate_store",
    "ATcTEntry",
    "ATcTError",
    "AmbiguousEntry",
    "MissingEntry",
    "reaction_uncertainty",
    "__version__",
]
