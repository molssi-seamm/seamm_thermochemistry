seamm_thermochemistry
======================

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
  the gaussian_step/psi4_step wide CSVs) into a `ThermoDB`. Needs the
  `import` extra (`pandas`, `openpyxl`).

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

Prototype. Schema and helper API are meant to be stable enough to build
`gaussian_step`/`psi4_step`/`vasp_step` adapters against, but:

- No bundled, pre-populated `data/thermochemistry.db` yet -- build one with
  the importers (see `scripts/build_prototype_db.py`).
- Packaging is a static version, not yet wired to `versioningit` like its
  sibling packages.
- The `settings` column is a single free-form string (e.g. `"encut=700eV"`)
  rather than normalized basis/cutoff columns -- fine for today's two
  producers, may want normalizing once ORCA joins.

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
