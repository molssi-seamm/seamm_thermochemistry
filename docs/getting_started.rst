Getting Started
===============

``seamm_thermochemistry`` is the shared atomic reference-energy database and
formation-energy arithmetic for SEAMM. A raw total energy from Gaussian,
ORCA, Psi4, or VASP has an arbitrary, code-dependent zero -- meaningless to
a non-expert user, and not comparable across codes. Re-expressing it as an
**energy** or **enthalpy of formation** (relative to the elements in their
standard states) fixes both problems at once, and this package is the
single place that logic and its underlying reference data live, replacing
what used to be an independent, largely-duplicated copy in each of
``gaussian_step``, ``psi4_step``, and ``vasp_step``.

Installing
----------

The core library only depends on ``seamm-util``::

    pip install seamm_thermochemistry

Two optional extras add more:

* ``seamm_thermochemistry[import]`` -- pulls in ``pandas``/``openpyxl`` so
  the legacy Excel/CSV master-file importers in :mod:`importers` are usable.
* ``seamm_thermochemistry[installer]`` -- pulls in ``seamm-installer`` so
  the ``seamm-thermochemistry-installer`` console script can fetch the
  published reference database from Zenodo.

Fetching the reference database
--------------------------------

The reference database itself (``thermochemistry.db``) is not shipped
inside the Python package -- it is a growing dataset published separately
on Zenodo, with a version DOI for every release. Fetch (or update) the
local copy with::

    seamm-thermochemistry-installer install

which downloads the latest version, verifies it against Zenodo's own
checksum, and registers its location in ``~/.seamm.d/seamm.ini``'s
``[thermochemistry]`` section. :data:`seamm_thermochemistry.DEFAULT_DB_PATH`
then resolves to that location automatically.

A quick example
----------------

.. code-block:: python

    from seamm_thermochemistry import ThermoDB, formation_energy

    with ThermoDB("my_reference.db") as db:
        db.add_element(1, "H", dfH0_0K=216.034, standard_state="1/2 H2(g)")
        db.add_element(8, "O", dfH0_0K=246.79, standard_state="1/2 O2(g)")
        db.add_atom_energy("H", "gaussian", "CBS-QB3", -1312.0)
        db.add_atom_energy("O", "gaussian", "CBS-QB3", -197400.0)

        dfE = formation_energy(
            {"H": 2, "O": 1}, system_energy, db, "gaussian", "CBS-QB3"
        )

Two reference conventions, one schema
--------------------------------------

The ``ref_type`` column on ``atom_energy`` distinguishes:

* ``"atom"`` -- the isolated gas-phase atom (Gaussian/Psi4/ORCA's
  convention). Pairs with an experimental anchor
  (``element.dfH0_0K``/``dfH0_298K``) to give a true, cross-code-comparable
  energy/enthalpy of formation.
* ``"element_phase"`` -- energy per atom of the element's standard-state
  phase (bulk metal, graphite, O2(g), ...) -- VASP's existing convention,
  needing no experimental anchor, and useful as a fallback reference for
  elements (e.g. Mn) where the free atom is a poor DFT target.

``formation_energy(..., anchor=True)`` matches the existing
``gaussian_step``/``psi4_step`` enthalpy-of-formation arithmetic exactly
(and gives an *energy* of formation, not an enthalpy, when the system
energy excludes the zero-point energy and the 0 K anchor is used).
``formation_energy(..., anchor=False)`` matches ``vasp_step``'s existing
``DfE0`` convention exactly.

Importing ORCA results
------------------------

Production ORCA atom-energy scans (one CSV per method, covering many
basis sets and elements) are imported with::

    seamm-thermochemistry-import-orca results.csv

Each value is vetted against its expected total spin (``S^2``) before
being added, and duplicate values already in the database are left alone
unless they conflict -- see
:func:`seamm_thermochemistry.importers.import_orca_atom_results` for the
full vetting and conflict-resolution rules, including
``--prefer-lower-energy`` for the hard-atom case where two runs converge
to different, equally spin-clean self-consistent solutions.
