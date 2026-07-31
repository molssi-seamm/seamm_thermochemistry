# -*- coding: utf-8 -*-

"""A detailed, citable text report of a formation-energy calculation.

Every code that reports a formation-referenced energy used to build (and
duplicate) its own ASCII table of "here is exactly what atom energies went
into this number" -- see e.g. the ~200-line ``calculate_enthalpy_of_formation``
this package's ``formation.py`` replaced the *arithmetic* of. This module is
the shared replacement for the *reporting* half: one place that turns a
composition + a ThermoDB + a set of system energies into the explicit,
per-atom breakdown table plus references, so a chemist can see exactly which
database, which atomic energy numbers, and which experimental citations
produced the headline number -- not just trust it.
"""

from datetime import datetime, timezone
from pathlib import Path

from tabulate import tabulate

from .formation import (
    formation_energy,
    formation_enthalpy,
    formation_gibbs_energy,
    MissingReferenceData,
)

__all__ = ["format_report"]

_DELTA = "\N{GREEK CAPITAL LETTER DELTA}"
_DEGREE = "\N{DEGREE SIGN}"


def _db_provenance(db):
    """A short, citable description of which database produced a report.

    A published database carries a Zenodo DOI (see `ThermoDB.doi`) -- the
    permanent, citable identifier for exactly this data snapshot, and by
    far the most useful thing to cite for reproducibility, so it takes
    priority whenever present. A working/unpublished database (the normal
    case during development) has no DOI yet, so this falls back to the
    file path and its last-modified date instead.
    """
    doi = db.doi()
    if doi:
        return f"{db.path} (DOI: https://doi.org/{doi})"
    try:
        mtime = datetime.fromtimestamp(Path(db.path).stat().st_mtime, tz=timezone.utc)
        mtime_str = mtime.strftime("%Y-%m-%d")
    except OSError:
        mtime_str = "unknown"
    return f"{db.path} (unpublished working copy, last updated {mtime_str})"


def _atom_table(composition, db, code, method, ref_type, settings, units):
    """Per-atom breakdown: symbol, count, E_ref, n*E_ref, and provenance
    (source/computed_date, when the atom_energy row carries them).

    Raises `MissingReferenceData` for any element with no tabulated
    (code, method, ref_type, settings) energy -- mirrors
    `atomization_energy`'s own check, so callers can rely on this table
    build failing exactly when the arithmetic itself would.
    """
    ref_energies = db.get_reference_energies(
        code, method, ref_type=ref_type, settings=settings, units=units
    )
    missing = sorted(el for el in composition if el not in ref_energies)
    if missing:
        raise MissingReferenceData(
            f"No {code}/{method} ({ref_type}) reference energy for: "
            f"{', '.join(missing)}"
        )

    rows = []
    ref_sum = 0.0
    for el, n in sorted(composition.items()):
        e = ref_energies[el]
        ref_sum += n * e
        detail = db.get_atom_energy_row(
            el, code, method, ref_type=ref_type, settings=settings
        )
        source = ""
        if detail:
            bits = [b for b in (detail.get("source"), detail.get("computed_date")) if b]
            source = " / ".join(bits)
        rows.append([f"{el}(g)", n, f"{e:.4f}", f"{n * e:.4f}", source])
    tmp = tabulate(
        rows,
        headers=["Atom", "n", f"E_ref ({units})", f"n*E_ref ({units})", "Source"],
        tablefmt="rounded_outline",
        colalign=("center", "center", "decimal", "decimal", "left"),
        disable_numparse=True,
    )
    return tmp, ref_sum


def _anchor_table(composition, db):
    """Per-atom experimental 0 K heat-of-formation anchor + citation."""
    rows = []
    for el, n in sorted(composition.items()):
        dfH0 = db.dfH0(el, at_0K=True)
        el_row = db.get_element(symbol=el) or {}
        citation = el_row.get("reference") or ""
        note = el_row.get("reference_note") or ""
        cite = f"{note}: {citation}" if note and citation else (note or citation)
        rows.append([f"{el}(g)", n, "" if dfH0 is None else f"{dfH0:.3f}", cite])
    return tabulate(
        rows,
        headers=["Atom", "n", "DfH0(0K) (kJ/mol)", "Reference"],
        tablefmt="rounded_outline",
        colalign=("center", "center", "decimal", "left"),
        disable_numparse=True,
    )


def _std_state_table(composition, db):
    """Per-atom standard-state-phase entropy + citation, for the Gibbs
    section."""
    rows = []
    for el, n in sorted(composition.items()):
        s_std = db.s298_std_state(el)
        el_row = db.get_element(symbol=el) or {}
        cite = el_row.get("s298_std_state_reference") or ""
        rows.append(
            [
                f"{el}(g)",
                n,
                "" if s_std is None else f"{s_std:.3f}",
                el_row.get("standard_state") or "",
                cite,
            ]
        )
    return tabulate(
        rows,
        headers=[
            "Atom",
            "n",
            "S298_std_state (J/mol K)",
            "Standard state",
            "Reference",
        ],
        tablefmt="rounded_outline",
        colalign=("center", "center", "decimal", "center", "left"),
        disable_numparse=True,
    )


def format_report(
    composition,
    db,
    code,
    method,
    *,
    ref_type="atom",
    settings="",
    units="kJ/mol",
    name=None,
    level_label=None,
    system_energy=None,
    system_enthalpy=None,
    system_gibbs_energy=None,
    temperature=None,
):
    """Build the detailed formation-energy report text.

    Parameters
    ----------
    composition, db, code, method, ref_type, settings, units
        As in `atomization_energy`.
    name : str, optional
        A human-readable name for the system (falls back to a generic
        label).
    level_label : str, optional
        A human-readable level of theory (falls back to ``f"{code}/{method}"``).
    system_energy : float, optional
        The system's electronic energy (0 K, ZPE-free). If given, the
        Atomization Energy and Energy-of-Formation (DfE0) sections are
        included.
    system_enthalpy : float, optional
        The system's total enthalpy H(T). If given (together with
        `temperature`), the Enthalpy-of-Formation (DfH(T)) section is
        included.
    system_gibbs_energy : float, optional
        The system's total Gibbs free energy G(T). If given (together with
        `temperature`), the Gibbs-Energy-of-Formation (DfG(T)) section is
        included.
    temperature : float, optional
        Temperature, in K, for the DfH(T)/DfG(T) sections.

    Returns
    -------
    str
        The formatted report. Never raises: any missing reference data is
        reported as a plain-text note in the relevant section instead of
        stopping the whole report (a molecule that is missing, say, the
        standard-state entropy for one of its elements should still get an
        atomization-energy and DfE0 section).
    """
    composition = dict(composition)
    name = name or "the system"
    level_label = level_label or f"{code}/{method}"

    title = f"Thermochemistry of {name} with {level_label}"
    lines = [title, "=" * len(title), ""]
    lines.append(f"Reference database: {_db_provenance(db)}")
    lines.append("")

    if system_energy is None:
        return "\n".join(lines)

    # ---- Atomization energy -------------------------------------------------
    try:
        table, ref_sum = _atom_table(
            composition, db, code, method, ref_type, settings, units
        )
        atomization = ref_sum - system_energy
    except MissingReferenceData as e:
        lines.append(f"Cannot calculate the atomization energy: {e}")
        return "\n".join(lines)

    lines.append("Atomization Energy")
    lines.append("-------------------")
    lines.append(
        f"The atomization energy, {_DELTA}atE, is the electronic energy to "
        "separate the system into gas-phase atoms, computed from the "
        f"per-atom reference energies below ({units}):"
    )
    lines.append("")
    lines.append(table)
    lines.append("")
    lines.append(f"    {_DELTA}atE = {atomization:.4f} {units}")
    lines.append("")

    # ---- Energy of formation, 0 K --------------------------------------------
    try:
        dfe0 = formation_energy(
            composition,
            system_energy,
            db,
            code,
            method,
            ref_type=ref_type,
            settings=settings,
            anchor=True,
            anchor_at_0K=True,
            units=units,
        )
        lines.append("Energy of Formation (0 K, electronic-only)")
        lines.append("-------------------------------------------")
        lines.append(
            f"{_DELTA}fE references the system to the elements in their "
            "standard states at 0 K, via each atom's experimental heat of "
            f"formation ({units}):"
        )
        lines.append("")
        lines.append(_anchor_table(composition, db))
        lines.append("")
        lines.append(f"    {_DELTA}fE = {dfe0:.4f} {units}")
        lines.append("")
    except MissingReferenceData as e:
        lines.append(f"Cannot calculate the energy of formation: {e}")
        lines.append("")

    # ---- Enthalpy of formation, T --------------------------------------------
    if system_enthalpy is not None and temperature is not None:
        try:
            dfh_t = formation_enthalpy(
                composition,
                system_enthalpy,
                temperature,
                db,
                code,
                method,
                ref_type=ref_type,
                settings=settings,
                units=units,
            )
            lines.append(
                f"Enthalpy of Formation, {_DELTA}fH{_DEGREE}({temperature:.2f} K)"
            )
            lines.append("-" * (len(lines[-1])))
            lines.append(
                f"{_DELTA}fH(T) extends {_DELTA}fE to the system's total enthalpy "
                f"H(T) (electronic + ZPE + thermal correction), treating each "
                "reference atom as a monatomic ideal gas -- see "
                "`formation_enthalpy` for the exact construction and its "
                "approximations."
            )
            lines.append("")
            lines.append(
                f"    {_DELTA}fH{_DEGREE}({temperature:.2f} K) = {dfh_t:.4f} {units}"
            )
            lines.append("")
        except MissingReferenceData as e:
            lines.append(f"Cannot calculate the enthalpy of formation: {e}")
            lines.append("")

    # ---- Gibbs energy of formation, T ----------------------------------------
    if system_gibbs_energy is not None and temperature is not None:
        try:
            dfg_t = formation_gibbs_energy(
                composition,
                system_gibbs_energy,
                temperature,
                db,
                code,
                method,
                ref_type=ref_type,
                settings=settings,
                units=units,
            )
            lines.append(
                f"Gibbs Energy of Formation, {_DELTA}fG{_DEGREE}({temperature:.2f} K)"
            )
            lines.append("-" * (len(lines[-1])))
            lines.append(
                f"{_DELTA}fG(T) additionally needs the entropy of each element's "
                "standard-state phase (below) -- see `formation_gibbs_energy` "
                "for the exact construction and its approximations."
            )
            lines.append("")
            lines.append(_std_state_table(composition, db))
            lines.append("")
            lines.append(
                f"    {_DELTA}fG{_DEGREE}({temperature:.2f} K) = {dfg_t:.4f} {units}"
            )
            lines.append("")
        except MissingReferenceData as e:
            lines.append(f"Cannot calculate the Gibbs energy of formation: {e}")
            lines.append("")

    return "\n".join(lines)
