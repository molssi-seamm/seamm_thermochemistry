"""The shared arithmetic: atomization and formation energies from a ThermoDB.

Every code (Gaussian, Psi4, ORCA, VASP, ...) that reports a formation energy
is doing the same subtraction:

    atomization = sum(n_X * E_ref(X)) - E(system)

and, when an experimental anchor is available,

    formation = sum(n_X * DfH(X)) - atomization

They differ only in the choice of E_ref (an isolated gas atom, or an atom of
the element's standard-state phase -- `ref_type` on ThermoDB) and whether an
experimental anchor is added at all. See
``~/Sites/reference-energy/2026-07-24_reference-energy/`` for the full
design rationale, including why the arbitrary code-dependent zero of a raw
total energy is the problem this solves, and the "no thermochemistry ->
energy of formation, not enthalpy" case this module is aimed at.

Sign convention, checked against the existing gaussian_step / psi4_step /
vasp_step implementations:

- ``atomization_energy`` matches Gaussian/Psi4's "E atomization" (>0 for a
  bound system: separated atoms sit higher in energy than the molecule).
- ``formation_energy(..., anchor=True)`` matches the existing
  ``DfH_at - H_atomization`` construction (Gaussian/Psi4's enthalpy of
  formation), generalized to also give a ZPE-free "energy of formation"
  when the 0 K anchor is used and `system_energy` excludes ZPE.
- ``formation_energy(..., anchor=False)`` matches VASP's existing
  ``DfE0 = E(system) - sum(n_X * mu_X)`` (materials formation-energy
  convention, no experimental data required) -- note this is
  ``-atomization_energy``, not `atomization_energy` itself; the two
  conventions report the quantity with opposite sign, which is why they
  are two functions rather than one flag.
"""

from collections import Counter

from seamm_util import Q_

__all__ = ["atomization_energy", "formation_energy", "MissingReferenceData"]


class MissingReferenceData(KeyError):
    """Raised when a composition needs reference data the ThermoDB doesn't have."""


def _as_counter(composition):
    return composition if isinstance(composition, Counter) else Counter(composition)


def atomization_energy(
    composition,
    system_energy,
    db,
    code,
    method,
    *,
    ref_type="atom",
    settings="",
    units="kJ/mol",
):
    """sum(n_X * E_ref(X)) - E(system): the energy to separate all atoms.

    This is exactly Gaussian/Psi4's "Atomization Energy" (electronic-only,
    no ZPE, if `system_energy` is itself ZPE-free) and VASP's "cohesive
    energy". Positive for a bound system.

    Parameters
    ----------
    composition : dict or collections.Counter
        {element_symbol: count}, e.g. {"C": 2, "H": 6}.
    system_energy : float
        The system's electronic energy, in `units`.
    db : ThermoDB
        Open reference-energy database.
    code, method : str
        Identify which tabulated atom energies to use.
    ref_type : {"atom", "element_phase"}
        See ``ThermoDB.add_atom_energy``.
    settings : str
        Must match the `settings` used when the atom energies were stored.
    units : str
        A pint-parseable unit string for `system_energy` and the return
        value, e.g. "kJ/mol", "eV", "E_h".
    """
    composition = _as_counter(composition)
    ref_energies = db.get_reference_energies(
        code, method, ref_type=ref_type, settings=settings, units=units
    )
    missing = sorted(el for el in composition if el not in ref_energies)
    if missing:
        raise MissingReferenceData(
            f"No {code}/{method} ({ref_type}) reference energy for: {', '.join(missing)}"
        )
    reference_sum = sum(n * ref_energies[el] for el, n in composition.items())
    return reference_sum - system_energy


def formation_energy(
    composition,
    system_energy,
    db,
    code,
    method,
    *,
    ref_type="atom",
    settings="",
    anchor=True,
    anchor_at_0K=True,
    units="kJ/mol",
):
    """The formation-referenced energy of a system.

    Parameters
    ----------
    composition, system_energy, db, code, method, ref_type, settings, units
        As in `atomization_energy`.
    anchor : bool
        If True (default), add the experimental atomic heat of formation so
        the result is a true energy/enthalpy of formation -- the physically
        meaningful, cross-code-comparable quantity this package exists for.
        If False, return the bare, no-experimental-data-needed quantity:
        VASP's existing ``DfE0`` convention (``ref_type="element_phase"``
        is the natural pairing, but any ref_type is accepted).
    anchor_at_0K : bool
        Use the 0 K experimental anchor (recommended: pairs with a
        ZPE-free `system_energy` to give an "energy of formation", degrades
        gracefully to DfH(0K) if ZPE is added back in, and to DfH(298K) if
        the molecule's thermal correction is added) or the 298 K anchor
        (matches the existing enthalpy-of-formation code paths directly).
        Ignored if `anchor` is False.

    Returns
    -------
    float
        Energy/enthalpy of formation (`anchor=True`) or the bare
        atomization-referenced energy in VASP's DfE0 sign convention
        (`anchor=False`), in `units`.
    """
    atomization = atomization_energy(
        composition,
        system_energy,
        db,
        code,
        method,
        ref_type=ref_type,
        settings=settings,
        units=units,
    )

    if not anchor:
        # VASP's DfE0 = E(system) - sum(n_X * mu_X) = -atomization_energy(...)
        return -atomization

    composition = _as_counter(composition)
    total_anchor = 0.0
    missing = []
    for el, n in composition.items():
        dfH0 = db.dfH0(el, at_0K=anchor_at_0K)
        if dfH0 is None:
            missing.append(el)
            continue
        total_anchor += n * Q_(dfH0, "kJ/mol").m_as(units)
    if missing:
        anchor_label = "0 K" if anchor_at_0K else "298 K"
        raise MissingReferenceData(
            f"No experimental {anchor_label} heat of formation for: {', '.join(sorted(missing))}"
        )

    return total_anchor - atomization
