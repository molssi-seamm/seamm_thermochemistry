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
on Zenodo, with a version DOI for every release. Fetch the local copy with::

    seamm-thermochemistry-installer install

which downloads the latest version, verifies it against Zenodo's own
checksum, and registers its location in ``~/.seamm.d/seamm.ini``'s
``[thermochemistry]`` section. :data:`seamm_thermochemistry.DEFAULT_DB_PATH`
then resolves to that location automatically.

``seamm-thermochemistry-installer update`` (which ``seamm-manager update`` runs)
brings the copy up to the latest published version -- but only if the local file is
exactly one of the published versions. A database you have edited, for instance with
curated or recomputed reference energies, is kept as it is and reported; to replace it
with the published version anyway, run
``seamm-thermochemistry-installer update --force``, which saves the old file as
``thermochemistry.db.bak-<date>`` first. Every replacement keeps such a backup, and a
download is checked before it replaces anything. ``install`` never overwrites an
existing database, and ``uninstall`` leaves one with local changes in place unless
given ``--force``.

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

``formation_enthalpy`` and ``formation_gibbs_energy`` take the system's total
enthalpy H(T) or Gibbs energy G(T) (electronic energy + ZPE + thermal
correction, on the same absolute scale as the atom energies) and reference it
to the elements in their standard states at T. Besides the 0 K anchor they need
each element's standard-state H(298.15 K)-H(0 K)
(``element.h298_minus_h0_std_state``, tabulated per formula unit of the
standard state and halved for the diatomic gases) and, for G, the
standard-state entropy (``element.s298_std_state``). The result is exact at
298.15 K; at other temperatures the standard-state thermal enthalpy is scaled
linearly in T, which is accurate for gaseous standard states and approximate
for solids far from 298 K. As a check, water's ΔfH(298.15 K) - ΔfH(0 K) comes
out as -2.90 kJ/mol with the JANAF data, matching experiment.

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

Active Thermochemical Tables
----------------------------

``seamm_thermochemistry.atct`` holds a local copy of one version of the Active
Thermochemical Tables (ATcT, https://atct.anl.gov), the mutually consistent
enthalpies of formation that the isodesmic (CBH) thermochemistry in SEAMM uses as
references. Each installation downloads its own copy::

    seamm-thermochemistry-atct download            # from atct.anl.gov
    seamm-thermochemistry-atct download --from all.json   # a copy saved in a browser
    seamm-thermochemistry-atct list
    seamm-thermochemistry-atct lookup "C/C=C/C"

The copy lives in ``Parameters/thermochemistry/ATcT/<version>/`` and is never
modified; a newer version goes in its own directory, so earlier results stay
reproducible. In Python::

    from seamm_thermochemistry import ATcT, reaction_uncertainty

    tables = ATcT()                      # the newest cached version
    entry = tables.lookup("CCF")         # fluoroethane, plain gas-phase entry
    entry["dfH0_298K"], entry["uncertainty_298K"], tables.version

The selection rules follow how ATcT labels its entries: gas phase with the right
charge; double-bond stereochemistry must match (cis and trans are separate
species), while a lone stereocentre is ignored (enantiomers share a value); the
plain ``*0`` entry is used over conformer entries unless ``conformer=`` names
one; a labelled spin state is used when ``multiplicity=`` matches it; and the
ortho/para nuclear-spin entries are never used. Anything still ambiguous stops
with the candidates listed; ``atct_id=`` overrides everything.

ATcT uncertainties are 95 % confidence limits at 298.15 K, and the values of
different species are correlated. ``tables.covariance(ids)`` fetches and caches
the covariance block, and ``reaction_uncertainty(coefficients, uncertainties,
covariance)`` gives the 95 % uncertainty of a reaction enthalpy. As with the
tables, a block can be saved from a browser (the response of
``https://atct.anl.gov/api/v1/covariance/species/?atctids=74-84-0*0,74-98-6*0``)
and cached with ``seamm-thermochemistry-atct covariance --from FILE``; the
diagonal is checked against the listed uncertainties, which ATcT rounds up to
0.01 kJ/mol. The correlations matter: for the CBH-2 reactions of the n-alkanes
the covariance gives reference uncertainties about 40 % smaller than
quadrature. Never convert an
ATcT value between 0 K and 298.15 K yourself, and never mix TN versions in one
reaction. The index needs RDKit (``pip install seamm-thermochemistry[cbh]``).

The fragment-energy store
-------------------------

Fragments recur between targets, so their energies are kept in a store,
``fragments.db`` next to the ATcT cache, and a species is computed once per
level of theory. A job stores its species' energies as properties on the
systems named ``CBH:<SMILES>``; afterwards, on the host (never on a compute
node, and never with the store on GPFS or NFS), promote them::

    seamm-thermochemistry-fragments harvest ~/SEAMM/Jobs/projects/x/Job_42
    seamm-thermochemistry-fragments list [SMILES]
    seamm-thermochemistry-fragments lookup CC --model-chemistry MOPAC:SQM@PM7
    seamm-thermochemistry-fragments duplicates

The harvest is insert-or-ignore with provenance (job, configuration, geometry
hash, date): the first value wins, and an attempt whose values differ by more
than 0.1 kJ/mol is logged, never overwritten, which catches another conformer
or SCF solution. A record is keyed by the species (canonical SMILES, charge,
multiplicity) and the *fingerprints* of the calculations that gave the energy
and the frequencies, opaque strings that the QC step writes next to its
results; E, ZPE and H(298.15 K) are separate fields, since the levels can
differ. A record with no fingerprint is keyed on the model chemistry and
code version, flagged unverified, and used only when nothing verified exists.
In Python::

    from seamm_thermochemistry import FragmentStore
    store = FragmentStore()                 # read-only, empty if no file
    record = store.get("CC", energy_fingerprint=fp, model_chemistry=mc)
    record["E"], record["ZPE"], record["H"], record["verified"]

Reads are best effort: a missing, locked or corrupt store gives None.

Laddered reference values (CBH-ANL)
-----------------------------------

ATcT values for every fragment are not self-consistent with one another at
the 0.5 kJ/mol level, and CBH results referenced to them drift with molecular
size (the n-alkanes, above). CBH-ANL (Elliott et al. 2023) instead takes ATcT
values only for the CBH-0 species (H2, CH4, H2O, NH3, HF, ... and the radical
CBH-0 species) and "ladders" every other reference species from its own CBH
reaction one rung below its size, computed at a high level of theory::

    seamm-thermochemistry-ladder build anl-f12 CC CCC CCCC \
        --energy "ORCA:QC@CCSD(T)-F12/bse:cc-pVTZ-F12" \
        --frequencies "ORCA:DFT@B2PLYP-D3/bse:cc-pVTZ" --zpe-scale 0.9863
    seamm-thermochemistry-ladder show anl-f12

The energies and frequencies come from the fragment store (harvest the jobs
first); the levels are settings fingerprints or model chemistry strings. The
ladder is saved as a named reference set in ``ladders/`` next to the store,
with every value's ATcT coefficients (so a reaction's ATcT uncertainty is
cᵀΣc exactly), its provenance, and, for every ladder species that ATcT has,
the difference from ATcT: their RMS is the measured method term, an upper
bound since it includes ATcT's own inconsistencies. Work at 0 K: ATcT's 0 K
values refer to the lowest conformer, while its 298.15 K values are
conformer averages. In Python ``Ladder.load(name)`` and
``ladder.entry(smiles)`` give the value in the layout of an ATcT entry.

Connectivity-based hierarchy (CBH) reactions
--------------------------------------------

``seamm_thermochemistry.cbh`` builds the error-cancelling reactions of
Ramabhadran and Raghavachari's connectivity-based hierarchy from a SMILES::

    from seamm_thermochemistry import cbh

    reactions = cbh.cbh_reactions("CSCCC(N)C(=O)O", rungs=[1, 2, 3])
    for rung, reaction in reactions.items():
        print(reaction)                 # CBH-2: ... -> ...
        reaction.coefficients           # {canonical SMILES: coefficient}
        reaction.references             # the species other than the target
        reaction.balanced, reaction.trivial

Products have positive coefficients, reactants negative and the target -1.
Every fragment keeps the bond orders, charges, unpaired electrons and
stereochemistry it has in the parent, and every cut bond is replaced by
hydrogens. By default aromatic ring systems, three- and four-membered rings and
charge-separated or hypervalent groups are kept whole (``protect="default"``;
``protect="none"`` gives the plain Kekulé hierarchy, and a sequence of
``"aromatic"``, ``"small rings"``, ``"groups"``, ``"all rings"`` or SMARTS
patterns chooses). A rung that reproduces the target is ``trivial``; a rung
the formula cannot build raises ``RungUnavailable``. ``cbh.unique_fragments``
gives the union of the reference species over several reactions, so that a
batch of targets computes each fragment once.
