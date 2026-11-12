# -*- coding: utf-8 -*-

"""Reading hand-curated reference values into a reference library.

The curated files (one CSV per source: printed JANAF, CODATA key values,
spectroscopic D0, W4-17, derived values) share common columns::

    smiles, formula, charge, multiplicity, name, phase, dfH0_kJmol,
    dfH298_kJmol, unc_kJmol, unc_level, status, citation, notes

plus columns of their own. This is a maintainer tool: it builds a library
file, which is then published; users install the published file.

Each row is added as a measurement of its source, then the value of each
species and source is chosen, the overrides applied and the derivations
checked. Rows the rules refuse are reported, not imported.
"""

import csv
import json
import math
from pathlib import Path

from .references import ReferenceError

#: The files, their source and the source's kind, rank and year
CURATED_SOURCES = {
    "janaf_printed.csv": (
        "JANAF-1998",
        "evaluated",
        10,
        1998,
        "NIST-JANAF Thermochemical Tables, 4th ed., M. W. Chase, J. Phys. Chem. "
        "Ref. Data Monograph 9 (1998)",
    ),
    "derived.csv": ("derived", "derived", 11, None, "Values derived from evaluations"),
    "codata.csv": (
        "CODATA-1989",
        "evaluated",
        20,
        1989,
        "J. D. Cox, D. D. Wagman and V. A. Medvedev, CODATA Key Values for "
        "Thermodynamics (Hemisphere, 1989)",
    ),
    "d0.csv": ("D0", "measured", 30, None, "Spectroscopic dissociation energies"),
    "w4_17.csv": (
        "W4-17",
        "literature",
        40,
        2017,
        "A. Karton, N. Sylvetsky and J. M. L. Martin, J. Comput. Chem. 38, 2063 "
        "(2017)",
    ),
}

#: Ground-state spin multiplicities of the atoms, for rows that leave it blank
ATOM_MULTIPLICITY = {
    "H": 2, "Li": 2, "Be": 1, "B": 2, "C": 3, "N": 4, "O": 3, "F": 2, "Na": 2,
    "Mg": 1, "Al": 2, "Si": 3, "P": 4, "S": 3, "Cl": 2, "K": 2, "Ca": 1, "Ga": 2,
    "Ge": 3, "As": 4, "Se": 3, "Br": 2,
}  # fmt: skip


def _float(text):
    text = (text or "").strip()
    return float(text) if text else None


def _status(row):
    status = (row.get("status") or "").strip().lower()
    notes = (row.get("notes") or "").strip()
    if (
        "held" in status
        or "do-not-use" in status
        or notes.upper().startswith(("HELD", "DO NOT USE"))
    ):
        return "held"
    if status.startswith("estimate") or row.get("janaf_estimate") == "True":
        return "estimate"
    if status.startswith("derived"):
        return "derived"
    if status.startswith("computed"):
        return "computed"
    if status in ("evaluated", "spectroscopic d0", "measured"):
        return "measured"
    raise ReferenceError(f"Unknown status {row.get('status')!r}")


def _uncertainty(row):
    """(uncertainty, level) in the library's vocabulary."""
    u = _float(row.get("unc_kJmol"))
    text = (row.get("unc_level") or "").lower()
    method = _float(row.get("method_err_1sigma_kJmol"))
    if method is not None:
        # W4-17: the atoms at ~2 sigma, plus twice the method's 1 sigma
        return math.hypot(u or 0.0, 2.0 * method), "mixed"
    if text.startswith("estimate"):
        return u, "estimate"
    if "unstated" in text:
        return u, "unstated"
    if "2 sigma" in text or "95" in text:
        # D0 rows combine the CODATA atoms (~2 sigma) with D0 as stated
        return u, ("mixed" if "+" in text else "2 sigma")
    if "1 sigma" in text:
        return u, "1 sigma"
    if u is None:
        return None, "none"
    return u, "unstated"


def _derivation(row):
    """The derivation of a derived row, from ``partners`` and
    ``drH298_kJmol``: partners are ``SMILES:SOURCE:COEFFICIENT`` separated
    by semicolons, with ``element`` as the source of an element in its
    reference state; the target has coefficient +1."""
    partners = (row.get("partners") or "").strip()
    drh = _float(row.get("drH298_kJmol"))
    if not partners or drh is None:
        return None
    reaction = [{"target": True, "coefficient": 1}]
    for item in partners.split(";"):
        smiles, source, coefficient = (x.strip() for x in item.rsplit(":", 2))
        entry = {"smiles": smiles, "coefficient": float(coefficient)}
        if source == "element":
            entry["element"] = True
        else:
            entry["source"] = source
        reaction.append(entry)
    return reaction, drh, _float(row.get("drH298_unc_kJmol"))


def import_curated(library, directory, *, files=None):
    """Read the curated CSVs in `directory` into `library`.

    Returns a list of (file, smiles, outcome, detail) for every row, where
    outcome is "added", "refused" or "held", and the problems found by
    :meth:`ReferenceLibrary.check`.
    """
    directory = Path(directory).expanduser()
    report = []
    touched = set()
    derivations = []
    for filename, (source, kind, rank, year, citation) in CURATED_SOURCES.items():
        if files is not None and filename not in files:
            continue
        path = directory / filename
        if not path.exists():
            continue
        library.add_source(source, kind, rank, citation=citation, retrieved=str(path))
        with open(path, newline="", encoding="utf-8") as fd:
            for row in csv.DictReader(fd):
                smiles = row["smiles"].strip()
                try:
                    charge = int(row.get("charge") or 0)
                    multiplicity = row.get("multiplicity", "").strip()
                    if multiplicity:
                        multiplicity = int(multiplicity)
                    else:
                        multiplicity = ATOM_MULTIPLICITY[smiles.strip("[]")]
                    status = _status(row)
                    uncertainty, unc_level = _uncertainty(row)
                    if status == "estimate":
                        unc_level = "estimate"
                    species = library.add_species(
                        smiles,
                        charge,
                        multiplicity,
                        formula=row.get("formula") or None,
                        name=row.get("name") or None,
                    )
                    h298_h0 = _float(row.get("H298_H0_mol")) or _float(
                        row.get("H298_H0_atom")
                    )
                    measurement = library.add_measurement(
                        species,
                        source,
                        citation=row["citation"],
                        dfh_298k=_float(row.get("dfH298_kJmol")),
                        dfh_0k=_float(row.get("dfH0_kJmol")),
                        uncertainty=uncertainty,
                        unc_level=unc_level,
                        status=status,
                        year=year,
                        h298_h0=h298_h0,
                        standard_state=row.get("ref_state") or "0.1 MPa",
                        ion_convention=row.get("ion_convention") or None,
                        comment=row.get("unc_level") or None,
                        notes=row.get("notes") or None,
                        retrieved=f"{filename}",
                        method=row.get("protocol") or row.get("conversion") or None,
                    )
                except (ReferenceError, KeyError, ValueError) as e:
                    report.append((filename, smiles, "refused", str(e)))
                    continue
                touched.add((species, source))
                if status == "derived" or row.get("partners"):
                    derivation = _derivation(row)
                    if derivation is not None:
                        derivations.append((measurement, derivation, row))
                report.append(
                    (filename, smiles, "held" if status == "held" else "added", status)
                )
    for measurement, (reaction, drh, drh_u), row in derivations:
        library.add_derivation(
            measurement,
            reaction,
            drh,
            drh_uncertainty=drh_u,
            citation=row.get("citation"),
            notes=json.dumps({"inputs": row.get("inputs")}),
        )
    for species, source in sorted(touched):
        library.choose(species, source)
    library.apply_overrides()
    return report, library.check()
