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

__all__ = [
    "atomization_energy",
    "formation_energy",
    "formation_enthalpy",
    "formation_gibbs_energy",
    "MissingReferenceData",
]

# The temperature at which the tabulated standard-state data apply, K.
_T_REF = 298.15

# Elements whose standard state is a diatomic gas X2. The database tabulates
# `h298_minus_h0_std_state` per formula unit of the standard state, so for these
# it is per X2 and must be halved to give the per-atom value (see derive.py).
_DIATOMIC_STANDARD_STATE = {"H", "N", "O", "F", "Cl", "Br"}


def _std_state_thermal_enthalpy(db, symbol, temperature):
    """H(T) - H(0) of the element's standard state, per atom, in kJ/mol.

    Exact at 298.15 K, where the database tabulates it (JANAF). At other
    temperatures it is scaled linearly in T. For the gaseous standard states
    (H2, N2, O2, F2, Cl2, and the noble gases) this is within about 1 kJ/mol
    per atom up to 1000 K, since their heat capacity is nearly constant; for
    solids such as graphite the error grows with |T - 298.15 K| and reaches
    several kJ/mol per atom by 1000 K. Returns None if the database lacks the
    value.
    """
    el = db.get_element(symbol=symbol)
    if el is None:
        return None
    h = el["h298_minus_h0_std_state"]
    if h is None:
        return None
    standard_state = el["standard_state"] or ""
    if standard_state.startswith("1/2 ") or symbol in _DIATOMIC_STANDARD_STATE:
        h /= 2
    return h * temperature / _T_REF


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
            f"No {code}/{method} ({ref_type}) reference energy for: "
            f"{', '.join(missing)}"
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
            f"No experimental {anchor_label} heat of formation for: "
            f"{', '.join(sorted(missing))}"
        )

    return total_anchor - atomization


def formation_enthalpy(
    composition,
    H_system,
    temperature,
    db,
    code,
    method,
    *,
    ref_type="atom",
    settings="",
    units="kJ/mol",
):
    """DfH(T): enthalpy of formation from the elements, at temperature T.

    Generalizes ``formation_energy(..., anchor=True, anchor_at_0K=True)``
    from the electronic-only, 0 K quantity to a real temperature. The
    Hess cycle is elements (standard states, T) -> gas atoms (T) -> molecule
    (T). The gas atoms' own H(T)-H(0) appears in both halves of the cycle
    and cancels exactly, leaving the standard-state elements' thermal
    enthalpy (the form used by G4 and similar composite methods)::

        DfH(T) = sum(n_X * [DfH0(X, 0 K) - (H_std(X, T) - H_std(X, 0))])
                 - [sum(n_X * E_ref(X)) - H_system]

    where H_std is the enthalpy of the element's standard state, per atom
    (`ThermoDB.element.h298_minus_h0_std_state`, halved for the diatomic
    gases). Equivalently, DfH(T) = sum(n_X * DfH0(X, T)) - [atomization
    enthalpy at T]; at 298.15 K the two forms use the same tabulated data
    and agree to rounding.

    Approximation: the standard-state H(T)-H(0) is tabulated at 298.15 K
    only and is scaled linearly in T elsewhere -- see
    `_std_state_thermal_enthalpy` for the size of the error. At 298.15 K the
    result is exact.

    History: until 2026.10.28 this function subtracted (5/2)RT per atom,
    the gas atom's own thermal enthalpy, instead of the standard-state
    term. That made DfH(298 K) too negative by 1.96 kJ/mol per H, 5.15 per
    C and 1.86 per O (5.8 kJ/mol for water).

    Parameters
    ----------
    composition, db, code, method, ref_type, settings, units
        As in `atomization_energy`.
    H_system : float
        The system's total enthalpy H(T) -- electronic + ZPE + thermal
        correction -- on the SAME absolute energy scale as the electronic
        energy stored for the reference atoms (not a thermal correction
        alone), in `units`.
    temperature : float
        Temperature, in K.

    Returns
    -------
    float
        DfH(T), in `units`.
    """
    composition = _as_counter(composition)
    ref_energies = db.get_reference_energies(
        code, method, ref_type=ref_type, settings=settings, units=units
    )
    missing = sorted(el for el in composition if el not in ref_energies)
    if missing:
        raise MissingReferenceData(
            f"No {code}/{method} ({ref_type}) reference energy for: "
            f"{', '.join(missing)}"
        )
    ref_sum = sum(n * ref_energies[el] for el, n in composition.items())

    anchor = 0.0
    missing_anchor = []
    for el, n in composition.items():
        dfH0 = db.dfH0(el, at_0K=True)
        h_std = _std_state_thermal_enthalpy(db, el, temperature)
        if dfH0 is None or h_std is None:
            missing_anchor.append(el)
            continue
        anchor += n * Q_(dfH0 - h_std, "kJ/mol").m_as(units)
    if missing_anchor:
        raise MissingReferenceData(
            "No experimental 0 K heat of formation and/or standard-state "
            f"H(298 K)-H(0 K) for: {', '.join(sorted(missing_anchor))}"
        )

    return anchor - ref_sum + H_system


def formation_gibbs_energy(
    composition,
    G_system,
    temperature,
    db,
    code,
    method,
    *,
    ref_type="atom",
    settings="",
    units="kJ/mol",
):
    """DfG(T): Gibbs energy of formation from the elements, at T.

    Same Hess-cycle construction as `formation_enthalpy`, extended to G.
    Requires each element's `s298_std_state` (the *standard-state phase's*
    molar entropy, per atom -- see `ThermoDB.s298_std_state`) in addition
    to the usual 0 K atom energies, heats of formation and standard-state
    H(298 K)-H(0 K); raises `MissingReferenceData` for any element missing
    any of them.

    Derivation: on the computed energy scale the standard-state element's
    Gibbs energy per atom is::

        G_std(X, T) = E_ref(X) - DfH0(X, 0K) + [H_std(X, T) - H_std(X, 0)]
                      - T * S_std(X, T)

    (the gas atom at 0 K sits DfH0(X, 0K) above the standard state, and
    has no ZPE), so::

        DfG(T) = G_system - sum(n_X * G_std(X, T))
               = sum(n_X * [DfH0(X, 0K) - (H_std(X, T) - H_std(X, 0))])
                 + T * sum(n_X * S298_std_state(X)) / 1000
                 - [sum(n_X * E_ref(X)) - G_system]

    No atomic entropy is needed: the gas atom is only an intermediate.

    Approximations: `S298_std_state` is used at its tabulated 298.15 K
    value regardless of T, and the standard-state H(T)-H(0) is scaled
    linearly from its 298.15 K value (see `_std_state_thermal_enthalpy`).
    At T=0 this reduces exactly to
    ``formation_energy(..., anchor=True, anchor_at_0K=True)`` evaluated
    with `G_system` in place of the electronic energy.

    History: until 2026.10.28 this function subtracted (5/2)RT per atom
    instead of the standard-state thermal enthalpy; see
    `formation_enthalpy`.

    Parameters
    ----------
    composition, db, code, method, ref_type, settings, units
        As in `atomization_energy`.
    G_system : float
        The system's total Gibbs free energy G(T), on the same absolute
        scale as the electronic energy, in `units`.
    temperature : float
        Temperature, in K.

    Returns
    -------
    float
        DfG(T), in `units`.
    """
    composition = _as_counter(composition)
    ref_energies = db.get_reference_energies(
        code, method, ref_type=ref_type, settings=settings, units=units
    )
    missing = sorted(el for el in composition if el not in ref_energies)
    if missing:
        raise MissingReferenceData(
            f"No {code}/{method} ({ref_type}) reference energy for: "
            f"{', '.join(missing)}"
        )
    ref_sum = sum(n * ref_energies[el] for el, n in composition.items())

    anchor = 0.0
    std_state_term = 0.0
    missing_anchor = []
    for el, n in composition.items():
        dfH0 = db.dfH0(el, at_0K=True)
        h_std = _std_state_thermal_enthalpy(db, el, temperature)
        s_std = db.s298_std_state(el)
        if dfH0 is None or h_std is None or s_std is None:
            missing_anchor.append(el)
            continue
        anchor += n * Q_(dfH0 - h_std, "kJ/mol").m_as(units)
        std_state_term += n * Q_(temperature * s_std / 1000.0, "kJ/mol").m_as(units)
    if missing_anchor:
        raise MissingReferenceData(
            "No experimental 0 K heat of formation, standard-state "
            "H(298 K)-H(0 K) and/or standard-state entropy for: "
            f"{', '.join(sorted(missing_anchor))}"
        )

    return anchor + std_state_term - ref_sum + G_system
