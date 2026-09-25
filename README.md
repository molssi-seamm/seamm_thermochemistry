seamm_thermochemistry
======================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/molssi-seamm/seamm_thermochemistry/workflows/CI/badge.svg)](https://github.com/molssi-seamm/seamm_thermochemistry/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/molssi-seamm/seamm_thermochemistry/branch/main/graph/badge.svg)](https://codecov.io/gh/molssi-seamm/seamm_thermochemistry/branch/main)

Shared atomic reference-energy database and formation-energy arithmetic for
SEAMM.

### Why

A raw total energy from Gaussian, ORCA, Psi4, or VASP has an arbitrary,
code-dependent zero -- meaningless to a non-expert SEAMM user, and not
comparable across codes. Re-expressing it as an **energy/enthalpy of
formation** (relative to the elements in their standard states) fixes both
problems at once. See the design doc for the full rationale:
`~/Sites/reference-energy/2026-07-24_reference-energy/`.

Today, `gaussian_step`, `psi4_step`, and `vasp_step` each carry an
independent, mostly-duplicated copy of this logic and its data (a
~5000-column, mostly-empty CSV per molecular code; a separate workbook for
VASP). This package is the single shared replacement:

- **`db.py`** -- a SQLite-backed `ThermoDB` with two tables (`element`: the
  experimental reference data every code needs; `atom_energy`: one row per
  computed atomic reference energy, keyed by element/code/method/ref_type/
  settings, with room for provenance) and a small helper API
  (`add_element`, `add_atom_energy`, `get_reference_energies`, `missing`,
  `dump_*_csv`, ...). Zero third-party dependencies beyond `seamm_util`.
- **`formation.py`** -- `atomization_energy()` and `formation_energy()`,
  the arithmetic every plugin's `calculate_enthalpy_of_formation` currently
  reimplements, generalized to also produce a ZPE-free **energy** of
  formation when no harmonic thermochemistry has been run.
- **`importers.py`** -- one-off loaders from the three legacy master files
  (Paul's experimental-data workbook, the VASP element-energy workbook, and
  the gaussian_step wide CSV) into a `ThermoDB`. Needs the `import` extra
  (`pandas`, `openpyxl`). (psi4_step's CSV is not imported: it is a copy of
  gaussian_step's, i.e. Gaussian numbers, not Psi4 results.)

### Two reference conventions, one schema

`ref_type` on `atom_energy` distinguishes:

- `"atom"` -- isolated gas-phase atom (Gaussian/Psi4/ORCA's convention, and
  the target for VASP once the atom-in-a-box calcs are wired in). Pairs
  with an experimental anchor (`element.dfH0_0K` / `dfH0_298K`) to give a
  true, cross-code-comparable energy/enthalpy of formation.
- `"element_phase"` -- energy per atom of the element's standard-state
  phase (bulk metal, graphite, O2(g), ...). VASP's existing convention
  (`element_energies.csv`'s plain `<method>@<encut>` columns) -- no
  experimental anchor needed, and useful as a fallback reference for
  elements (e.g. Mn) where the free atom is a poor DFT target.

`formation_energy(..., anchor=True)` matches the existing
`gaussian_step`/`psi4_step` enthalpy-of-formation arithmetic exactly (and
gives an energy of formation, not enthalpy, when `system_energy` excludes
ZPE and the 0 K anchor is used). `formation_energy(..., anchor=False)`
matches `vasp_step`'s existing `DfE0` exactly. Both are exercised in
`tests/test_formation.py` against hand-worked numbers.

### Status

The reference database is published on Zenodo (a DOI per version) and
fetched with `seamm-thermochemistry-installer install` -- not bundled in
the Python package. `gaussian_step`'s `calculate_energy_of_formation`
already consumes it in production. Current coverage: the full Gaussian
composite-method/basis grid, VASP (PBE family, both the isolated-atom and
standard-state-phase conventions), and ORCA (several DFT methods across the
full def2 basis family), all vetted and imported via
`seamm-thermochemistry-import-orca` / the `importers` module. There are no
Psi4 atom energies yet: the earlier "psi4" rows were copies of Gaussian's.

### Computing ORCA atom energies

`scripts/orca_atom_multistart.py` computes the ORCA atomic reference energies.
An atom's reference energy for a method is defined as **the lowest-energy SCF
solution of that method at the experimental ground-state spin multiplicity**.
That is well defined even where the method orders the atom's states differently
from experiment. There, though, the atomization route to formation energies is
itself becoming unreliable, and for careful work reaction energies to
well-known species (e.g. H2O -> H2 + 1/2 O2) are the better tool. Formation
energies remain far more meaningful than raw total energies, for users and for
machine-learning training data alike.

A single SCF start is not enough to find that solution for open-shell atoms: Fe
from ORCA's PModel guess lands in 3d7 4s1, 120 kJ/mol above 3d6 4s2, and small
setting changes flip Nb and Ho between solutions over 80 kJ/mol apart. So the
script runs several starts per element, method and basis (PBE0/def2-SV(P)
orbitals, PModel, HCore, Hueckel, PAtom), shares every state found in one basis
with all the others, and keeps the lowest SCF energy. <S**2> is recorded and
contamination over 20% flagged, but it does not disqualify a solution. It uses
exact exchange (`NoCOSX`), because ORCA's default COSX gets some lone atoms
wrong (He ~1, Na and Mg ~5 kJ/mol), and full orbital convergence
(`ConvCheckMode 0`) for the MP2 part of double hybrids. It writes a
`Results.csv` that `seamm-thermochemistry-import-orca` reads, plus a log of
every run and of each choice.

One known simplification: the `settings` column is a single free-form
string (e.g. `"encut=700eV"`) rather than normalized basis/cutoff columns
-- fine across today's producers, may want normalizing if that stops being
true.

### Quick start

```python
from seamm_thermochemistry import ThermoDB, formation_energy

with ThermoDB("my_reference.db") as db:
    db.add_element(1, "H", dfH0_0K=216.034, standard_state="1/2 H2(g)")
    db.add_element(8, "O", dfH0_0K=246.79, standard_state="1/2 O2(g)")
    db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
    db.add_atom_energy("O", "gaussian", "CBS-QB3", -197400.0)

    dfE = formation_energy({"H": 2, "O": 1}, system_energy, db, "gaussian", "CBS-QB3")
```

### Copyright

Copyright (c) 2026, MolSSI SEAMM
