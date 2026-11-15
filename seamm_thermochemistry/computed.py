# -*- coding: utf-8 -*-

"""Computed reference values: the "computed" layer of the reference library.

Enthalpies of formation computed by SEAMM with the high-level ("tier-A")
recipes, for reference species that neither ATcT nor the experimental library
has. They are kept apart from the experimental library (Paul, 2026-11-15): in
their own file, ``references-computed.db``, beside ``references-official.db``,
whose layer is "computed" and whose one source, :data:`SOURCE`, is of kind
"computed".

How they are used (the resolver, :mod:`.resolver`):

* Only with "Supplementary references" = "compilations and computed".
* After every experimental value: a species with an experimental value never
  takes a computed one.
* A result that uses one is labelled "computed references".

Each row of the CSV read by :func:`import_computed` is one species. Columns:

``smiles``, ``formula``, ``charge``, ``multiplicity``, ``name``, ``phase`` (g)
    The species.
``dfH0_kJmol``, ``dfH298_kJmol``
    ΔfH° at 0 K and 298.15 K (at least one).
``unc_kJmol``, ``unc_level``, ``basis``
    The uncertainty, at 95 % (``unc_level`` must say "95 %"), and how it was
    set (from the recipe, before any CBH result exists).
``status``
    "computed (tier-A)".
``recipe``, ``job``, ``iteration``, ``model_chemistry``
    The recipe, and the SEAMM job (its id only) and loop iteration that
    computed it.
``geometry``, ``conformer``, ``H298_H0_kJmol``, ``thermal``
    The geometry level and conformer, the species' own H298 - H0 and how it
    was computed (e.g. harmonic, ZPE scale factor).
``atomization``, ``atom_references``
    The route to ΔfH and the atoms' reference values (e.g. ATcT 1.228).
``computed_on``, ``citation``, ``notes``
    The date (YYYY-MM-DD), a citation (e.g. "SEAMM tier-A, job 5738
    iteration 3, <recipe>") and anything else.

Importing is idempotent: a row whose species, job and iteration are already in
the file is passed over. A row from another job for a species already there
replaces it: the earlier value is kept, held, with a note saying which job
replaced it.

The file may be published one day, so a value that looks like a local path is
refused.
"""

import csv
import json
import re
from pathlib import Path

from .references import COMPUTED_RANK, ReferenceError, references_directory

#: The one source of the computed layer
SOURCE = "SEAMM-tierA"
#: Its rank: after every experimental source, and after W4-17
SOURCE_RANK = COMPUTED_RANK + 10
SOURCE_CITATION = (
    "Enthalpies of formation computed with SEAMM's high-level (tier-A) recipes; "
    "each value cites its recipe and SEAMM job"
)
#: The layer recorded in the file's meta
LAYER = "computed"

STATUS = "computed (tier-A)"
REQUIRED = ("smiles", "recipe", "job", "unc_kJmol", "unc_level", "computed_on")
#: The provenance kept with each value, in its notes (as JSON)
PROVENANCE = (
    "recipe",
    "job",
    "iteration",
    "model_chemistry",
    "geometry",
    "conformer",
    "thermal",
    "atomization",
    "atom_references",
    "basis",
    "computed_on",
    "notes",
)

_PATH = re.compile(r"(^|[\s'\"(=])(~/|/Users/|/home/|/D\d/|/projects/|[A-Za-z]:\\)")


def computed_references_path(path=None, official=None):
    """The computed layer, ``references-computed.db``: beside the official
    library `official` when given, else beside the ATcT cache."""
    if path is not None:
        return Path(path).expanduser()
    if official is not None:
        return Path(official).expanduser().parent / "references-computed.db"
    return references_directory() / "references-computed.db"


def _float(text):
    text = (text or "").strip()
    return float(text) if text else None


def _check(row):
    """Refuse a row that lacks what the layer needs or carries a local path."""
    for key in REQUIRED:
        if not (row.get(key) or "").strip():
            raise ReferenceError(f"no {key}")
    if (row.get("status") or "").strip() != STATUS:
        raise ReferenceError(f"status must be {STATUS!r}, not {row.get('status')!r}")
    if "95" not in row["unc_level"]:
        raise ReferenceError(
            f"the uncertainty must be at 95 %, not {row['unc_level']!r}"
        )
    if (row.get("phase") or "g").strip() != "g":
        raise ReferenceError("only gas-phase values")
    if (
        _float(row.get("dfH298_kJmol")) is None
        and _float(row.get("dfH0_kJmol")) is None
    ):
        raise ReferenceError("no enthalpy of formation given")
    if not re.fullmatch(r"\d+", row["job"].strip()):
        raise ReferenceError(f"job must be a SEAMM job id, not {row['job']!r}")
    if not re.fullmatch(r"\d{4}-\d\d-\d\d", row["computed_on"].strip()):
        raise ReferenceError(
            f"computed_on must be YYYY-MM-DD, not {row['computed_on']!r}"
        )
    for key, value in row.items():
        if value and _PATH.search(value):
            raise ReferenceError(f"{key} looks like a local path: {value!r}")


def import_computed(library, path):
    """Read a CSV of computed values (see the module docstring) into
    `library`, the computed layer.

    Returns a list of (smiles, outcome, detail): outcome is "added",
    "replaced" (another job's value is now held), "present" (this job's value
    is already there) or "refused".
    """
    if library.layer != LAYER:
        raise ReferenceError(
            f"{library.path} is the {library.layer!r} layer; computed values go "
            f"only in the {LAYER!r} layer ({computed_references_path().name})."
        )
    library.add_source(SOURCE, "computed", SOURCE_RANK, citation=SOURCE_CITATION)
    db = library._writer()
    report = []
    with open(Path(path).expanduser(), newline="", encoding="utf-8") as fd:
        for row in csv.DictReader(fd):
            row = {k: (v or "").strip() for k, v in row.items() if k}
            smiles = row.get("smiles", "")
            try:
                _check(row)
                charge = int(row.get("charge") or 0)
                multiplicity = int(row.get("multiplicity") or 1)
                species = library.add_species(
                    smiles,
                    charge,
                    multiplicity,
                    formula=row.get("formula") or None,
                    name=row.get("name") or None,
                )
            except (ReferenceError, ValueError) as e:
                report.append((smiles, "refused", str(e)))
                continue
            job = row["job"]
            iteration = row.get("iteration") or ""
            source_ref = f"SEAMM job {job}" + (
                f", iteration {iteration}" if iteration else ""
            )
            earlier = [
                m
                for m in library.measurements(species, SOURCE)
                if m["status"] == "computed"
            ]
            if any(m["source_ref"] == source_ref for m in earlier):
                report.append((smiles, "present", source_ref))
                continue
            provenance = {k: row[k] for k in PROVENANCE if row.get(k)}
            try:
                library.add_measurement(
                    species,
                    SOURCE,
                    citation=row.get("citation") or f"SEAMM tier-A, {source_ref}",
                    dfh_298k=_float(row.get("dfH298_kJmol")),
                    dfh_0k=_float(row.get("dfH0_kJmol")),
                    uncertainty=_float(row["unc_kJmol"]),
                    unc_level="95 %",
                    status="computed",
                    year=int(row["computed_on"][:4]),
                    h298_h0=_float(row.get("H298_H0_kJmol")),
                    standard_state="0.1 MPa",
                    ion_convention="stationary-electron" if charge else None,
                    source_ref=source_ref,
                    method=row["recipe"],
                    comment=row.get("basis") or None,
                    notes=json.dumps(provenance, ensure_ascii=False),
                )
            except ReferenceError as e:
                report.append((smiles, "refused", str(e)))
                continue
            for m in earlier:
                notes = json.loads(m["notes"] or "{}")
                notes["replaced by"] = source_ref
                db.execute(
                    "UPDATE measurement SET status = 'held', notes = ? WHERE id = ?",
                    (json.dumps(notes, ensure_ascii=False), m["id"]),
                )
            db.commit()
            library.choose(species, SOURCE)
            report.append((smiles, "replaced" if earlier else "added", source_ref))
    return report
