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
ortho/para nuclear-spin entries of small species (H2, H2O, H2O2, CH2O, ketene: at most
three heavy atoms) are never used, while on bigger species "ortho" and "para" name
positional isomers (o-tolualdehyde, p-benzyne), which are. Anything still ambiguous stops
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
    seamm-thermochemistry-fragments superseded [SMILES]

The harvest is insert-or-ignore with provenance (job, configuration, geometry
hash, date): the first value wins, and an attempt whose values differ by more
than 0.1 kJ/mol is logged, never overwritten, which catches another conformer
or SCF solution. A record is keyed by the species (canonical SMILES, charge,
multiplicity) and the *fingerprints* of the calculations that gave the energy
and the frequencies, opaque strings that the QC step writes next to its
results; E, ZPE and H(298.15 K) are separate fields, since the levels can
differ. A record with no fingerprint is keyed on the model chemistry and
code version, flagged unverified, and used only when nothing verified exists.

The conformer is not part of the key, so a value found later to be wrong (a
fragment harvested in a higher-energy conformer, say) is replaced only on
purpose: harvest the corrected job with ``--supersede "REASON"``. The stored
values that differ are then replaced, and the old ones move to the
``superseded`` table with the date and the reason, where ``superseded`` lists
them; the record keeps its key, so lookups find the new values. In Python,
``store.supersede(record, reason)`` does the same for one record.

Values that differ by less than 0.1 kJ/mol count as duplicates and are not replaced.
When a correction must reach every record of a level of theory, however small (a
dispersion correction applied to all species, say), add ``--force`` to
``--supersede "REASON"``, or ``force=True`` in Python: every value that differs at
all is replaced, with the old ones kept in ``superseded`` as before.

The store keeps E, ZPE and H unscaled, as the calculation gave them. A ZPE
scale factor is applied when the values are used (for example
``seamm-thermochemistry-ladder build ... --zpe-scale 0.9863``), so changing a
factor never needs a new harvest.

In Python::

    from seamm_thermochemistry import FragmentStore
    store = FragmentStore()                 # read-only, empty if no file
    record = store.get("CC", energy_fingerprint=fp, model_chemistry=mc)
    record["E"], record["ZPE"], record["H"], record["verified"]

Reads are best effort: a missing, locked or corrupt store gives None.

Model chemistries are stored as the bare model and the program, e.g.
``CMP@R2SCAN-3C`` and ``ORCA``, however they were given: the ORCA sub-steps
name a level ``CMP@R2SCAN-3C`` and the Model Chemistry step
``ORCA:CMP@R2SCAN-3C``. Lookups accept either spelling, and still find rows
written before the names were normalized (``split_model_chemistry`` gives the
two parts).

Every fingerprint argument of the readers takes one fingerprint or several. A
program's settings can depend on the species, not only on the level: ORCA uses
exact exchange (NoCOSX) for atoms and for Na, Mg, Zn, B and P, so a boron
species' energy has another fingerprint than methane's at the same level. A
caller passes the fingerprints it expects for each species (for ORCA,
``seamm_exec.expected_fingerprints`` predicts them)::

    store.get("CC", energy_fingerprint=[fp_as_written, fp_cosx_moot])

Write-back from jobs, and the harvest of jobs that come back
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Thermochemistry step adds the species it computes to the store itself
(write-back), through ``harvest_configuration``, the same code the harvest runs for
each system of a job, so the rows are the same either way. It does so only when the
store is on a local disk of the machine running the job: SQLite's locking is not
reliable on network and parallel filesystems (NFS, GPFS, Lustre).
``store_is_writable()`` makes the check; ``fragment-store-writable = auto | yes |
no`` in the ``[thermochemistry]`` section of ``~/.seamm.d/seamm.ini`` overrides it.

When the step cannot write the store, it leaves ``harvest_fragments.json`` in the
job's directory. A JobServer runs the plug-ins' ``job_finished`` hooks when a job
ends and its files are back on the JobServer's host; this package's hook harvests
a job carrying that marker into the store there, once (the marker becomes
``harvest_fragments.done.json``, with the tallies). A job run without a JobServer
is harvested by hand, as before.

Results are checked before they go in, by write-back and by the harvest:

* **The structure.** The bonds are perceived again from the coordinates, and must
  be those of the species in the key (its SMILES), whose charge and radicals must
  fit the structure's charge and multiplicity; and no atoms may be too close
  together (``molsystem.vsepr.check_contacts``, which also catches the ligands of a
  hypervalent centre on top of each other). A species that rearranged while it was
  optimized is refused. ``harvest --no-structure-check`` turns this off.
* **Imaginary modes** (write-back): a geometry with a mode beyond 20i is refused,
  with every energy on it. The harvest stores such a geometry, for its provenance,
  but lookups never return it.
* **Conformers** (write-back): a geometry lower in E + ZPE than every stored one of
  the species at its level, by more than 0.1 kJ/mol, and structurally distinct
  from each (RMSD above 0.05 Å after the best, symmetry-aware alignment), is added
  as the conformer ``auto:<geometry hash>``; lookups take the lowest. The same
  minimum computed again is only logged.

Geometries and energies at different levels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thermochemistry often takes the energy from a costly level of theory (for example
revDSD-PBEP86-D4 or CCSD(T)-F12) on the geometry and frequencies of a cheaper one
(for example r2SCAN-3c), and one geometry serves several energy levels. The store
keeps these in two more tables, beside the records above:

* ``geometry``: one row per conformer of a species at a geometry/frequency level,
  keyed by the species, charge, multiplicity, the frequency fingerprint and a
  conformer label. It holds the coordinates and their hash, the energy, the
  unscaled ZPE and H(298.15 K) at that level, and, when the job stored them, the
  harmonic frequencies, the number of imaginary modes, the rotational symmetry
  number and the moments of inertia.
* ``energy``: one row per single point on a geometry, keyed by the geometry's hash
  and the calculation's fingerprint. The hash is the one QC steps record for a
  calculation (``seamm_util.fingerprint.geometry_fingerprint``), so copies of the
  store on other hosts keep working.

A combined value is never stored, only its parts: the energy at one level and the
ZPE and thermal correction of the geometry's level are combined when they are used,
as ZPE factors are applied at use.

``harvest`` fills both tables from the same properties as the records: an energy
with a ZPE or enthalpy of the same level makes a geometry row, and every energy
makes an energy row on the configuration's geometry. It guards the pairing:

* an energy whose calculation record names a different geometry is refused and
  reported;
* a geometry row has its geometry and frequencies at one level: a row whose
  optimization was at another model chemistry is refused;
* a geometry with imaginary frequencies (more than 20i cm⁻¹) is stored, but not
  offered for reuse.

Jobs whose systems are not named ``CBH:<SMILES>`` (for example by an ATcT id, with
several conformers of a species) are harvested with a species map, a CSV file with the
columns ``job``, ``system name`` and ``species`` (the store key, with the stereo of the
structure), and optionally ``configuration name``, ``charge``, ``multiplicity`` and
``conformer label``::

    seamm-thermochemistry-fragments harvest --species-map map.csv --tables geometry JOB ...

``--tables geometry`` fills only the geometry and energy tables, which can hold several
conformers of a species; the records cannot. Systems the map does not list are skipped
and reported.

For jobs that stored only the energies, the map can also give each structure's thermal
corrections, in the columns ``E (kJ/mol)``, ``ZPE (kJ/mol)`` (unscaled), ``H298
(kJ/mol)`` (H at 298.15 K) and ``symmetry number``. They make the geometry row, with
H - E added to the job's energy, provided the map's E agrees with the job's within
0.1 kJ/mol; otherwise the structure is reported and gets no geometry row.

Duplicates and ``--supersede`` work as for the records, with the old values in a
``row_history`` table. To see what the store holds for a species::

    seamm-thermochemistry-fragments geometries CCO --model-chemistry ORCA:DFT@R2SCAN-3C

In Python::

    rows = store.species_at(
        "CCO",
        frequency_fingerprint=freq_fp,   # or frequency_model_chemistry=...
        energy_fingerprint=energy_fp,    # or energy_model_chemistry=...
    )
    rows[0]["E"], rows[0]["ZPE"], rows[0]["H - E"], rows[0]["geometry"]["conformer"]

The ladder (below) uses these tables first, and falls back on the records.

ZPE scale factors
-----------------

Harmonic zero-point energies are too large; scaled by a factor fitted for the level of
theory of the frequencies they come close to the anharmonic ones. The factors, with their
sources, are in ``seamm_thermochemistry.zpe``, looked up from the geometry and frequency
level's model chemistry::

    from seamm_thermochemistry import zpe_scale_factor
    f = zpe_scale_factor("ORCA:CMP@R2SCAN-3C")
    f.factor, f.source, f.flag    # 0.9899, the fit, "fitted on molecules of up to 6 atoms..."

=================================  ======  =========================================
Level                              Factor  Source and caution
=================================  ======  =========================================
r2SCAN-3c                          0.9899  fitted to 22 W4-2006 anharmonic ZPVEs;
                                           molecules of up to 6 atoms, so an
                                           extrapolation for larger ones
PBE0(-D4)/def2-TZVP                0.975   Truhlar database v5, PBE0/MG3S:
                                           basis-mismatched
B2PLYP(-D3)/cc-pVTZ                0.9863  the CBH-ANL protocol
=================================  ======  =========================================

A dispersion correction in the method's name is ignored in the lookup; a level that is
not in the table has no factor, and the Thermochemistry step then warns and uses 1. The
stored ZPEs stay unscaled.

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

``cbh.lowest_conformer(smiles, n=20)`` gives the lowest MMFF94 conformer of
``n`` RDKit embeddings (symbols, coordinates, energies), the right start for a
species whose lowest rotamer is wanted: a single embedding of a long chain is
often gauche.

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
