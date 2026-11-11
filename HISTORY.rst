=======
History
=======

2026.11.11.9 -- Write-back: quarantined rows give way; a log of every write
    * A valid new row now replaces a row quarantined by hand (model chemistry
      ``UNUSABLE: <reason>``) that holds its key, in the fragment, geometry and energy
      tables, with the old values kept as superseded ("replaced a quarantined row:
      <reason>"). Before, the quarantined row kept the key, so a species recomputed
      correctly (e.g. [CH3-], cyclic C3) could never get a usable record.
    * Quarantined rows are no longer counted as conformers when write-back decides
      whether a new geometry is a lower conformer.
    * The store keeps a log of what every harvest and every job's write-back did
      with each row -- added, a duplicate, logged, superseded or refused, and why --
      in the new ``writeback_log`` table (schema 4), shown by
      ``seamm-thermochemistry-fragments writeback-log [--job J] [--species S]``.

2026.11.11.8 -- Official stores: triplet ground states checked too
    * ``official build``'s spin check also refuses a singlet row of a species whose
      ground state is a triplet that a SMILES writes closed-shell
      (``cbh.TRIPLET_GROUND_STATES``: O2, S2, SO, OCCO, cyclic C3), as it already did
      for the singlets of ``cbh.SINGLET_GROUND_STATES``.

2026.11.11.7 -- CBH: aryl neighbours of hypervalent centres; more ground states
    * A hypervalent centre's neighbour in a protected ring stays with its ring and
      is capped like a methyl in the centre's fragments, so aryl sulfones,
      sulfoxides, sulfonamides, sulfonates and phosphonates keep CBH reactions
      (methyl phenyl sulfone: methane + PhSO2Me -> benzene + dimethyl sulfone)
      rather than becoming one trivial group with their ring.
    * Halogen cations bonded to two or more heavy atoms (the [Cl+3] of FClO3 written
      with charge separation) are hypervalent centres too.
    * The ground-state table gains the singlets HCF, HCCl, GeH2, AlH, BF and BCl,
      and the triplets that a SMILES writes closed-shell: O2, S2, SO, OCCO and
      cyclic C3.

2026.11.11.6 -- Quarantined rows left alone and left out
    * Rows quarantined by hand in a fragment store -- unverified, with the model
      chemistry ``UNUSABLE: <reason>`` -- are never found by a lookup. Now ``migrate``
      also leaves their text alone (it no longer reads "UNUSABLE" as a program),
      ``export`` skips them, and ``official build`` drops them and lists them with
      their reasons in its report.

2026.11.11.5 -- Checks for official stores; the structure guard's false refusals
    * ``official build`` checks every species state before writing anything: no
      saddle points, the structure guard, and the ground-state spin where it is known
      (the singlet silylenes, dihalocarbenes and BH of ``cbh.SINGLET_GROUND_STATES``,
      and ATcT's labelled spin states). A failing state stops the build unless excluded
      with a reason (``--exclude``). ``--frequency-fingerprint`` and
      ``--energy-fingerprint`` keep only the rows reachable at those levels. A report
      (``<file>.report.json``) records the sources, levels, rows, exclusions,
      SHA-256 and the seamm-thermochemistry that built it. An exclusion can name one
      conformer (a ``conformer`` column), keeping the species' good conformers.
    * The structure guard accepts a structure whose bonds are its species' at
      RDKit's default perception or at a covalent-radius factor from 1.15 to 1.5:
      the diagonals of strained four-membered rings and cages (azete, oxetene,
      N2O3) no longer read as bonds, nor do long real bonds go missing. Structures
      that rearranged or dissociated still fail at every factor.
    * Bugfix: the structure guard of 2026.11.11.1 refused two or more fluorines on
      one atom (CF2, CF3 and NF3 groups, so every fluorinated carbonate): it held them
      to the 2.2 Å meant for fluorines on different atoms, but such fluorines are
      2.12-2.17 Å apart. They now need only 2.0 Å.
    * Bugfix: the structure guard refused the singlet ground states of species
      whose SMILES has two radical electrons on one atom (BH, SiCl2, CF2, ...): it
      now checks only that the spin's parity fits.

2026.11.11.4 -- CBH: hypervalent and charge-separated groups kept whole
    * A hypervalent centre (S, Se, Te, P, As, Cl, Br, I with more bonds than its
      normal valence: sulfonyl, sulfite, phosphoryl, SF5-R, ClF3, ...) is kept
      together with all its heavy neighbours, so CBH-0 gives real molecules
      (sulfuric and phosphoric acid, methyl hydrogen sulfate) rather than
      hydrides such as H2SO2, H3PO, SH4 or SH6. A molecule that is one such group
      (SF6, dimethyl sulfone) is trivial at every rung.
    * Two bonded atoms with opposite formal charges are kept together: CO written
      [C-]#[O+] no longer gives [CH3-] + [OH3+], and aromatic N-oxides are now
      protected too.
    * Hypervalent halogens (ClF3, BrF5), which RDKit's valence check rejects,
      are read and canonicalized; other over-valent atoms (O=[N]=O) still are not.
    * ``cbh.ground_state_multiplicity`` gives a species' spin multiplicity, with the
      singlet silylenes, dihalocarbenes and BH that a SMILES writes as triplets.
    * Over the 1,694 PPE1694 molecules only the 32 with such groups change their
      reactions; the others are unchanged.

2026.11.11.2 -- The official and local fragment stores
    * The fragment store has two layers: the local ``fragments.db``, which this
      installation writes, and the official ``fragments-official.db``, published by
      SEAMM and only ever read (without locking, so it may be on a network
      filesystem). ``LayeredStore`` reads them as one, the local store first; an
      official correction (a superseded value) beats the same old value in the
      local store; differences over 0.1 kJ/mol are noted; the lowest conformer of
      either is used; and an energy may come from the other store when it is on the
      same geometry. Every row says which store it came from.
    * New commands: ``migrate`` (brings an older store up to date after a dated
      backup: the meta table, model chemistries stored bare with their program),
      ``merge`` (adds the rows of other copies), ``export`` (the local rows the
      official store lacks, as a contribution), ``official build`` (a new official
      version from the previous one and contributions) and ``duplicates
      --layers``. ``lookup`` and ``geometries`` read both stores.
    * Stores have a ``meta`` table (schema 3) recording their layer and, for an
      official store, its version. Older stores are read as local ones.

2026.11.11.1 -- Write-back from jobs; harvest of returning jobs; a structure guard
    * ``harvest_configuration`` puts one computed structure into the fragment store:
      what the harvest does for each system of a job, and what the Thermochemistry
      step now does for each species it computes. With ``write_back`` a geometry
      with an imaginary mode is refused, and a lower, structurally distinct
      conformer (RMSD > 0.05 Å) is added as ``auto:<hash>`` instead of only logged.
    * The harvest refuses a structure whose bonds, perceived again from its
      coordinates, are not those of its species, or whose charge or multiplicity
      does not fit, or that has atoms too close together (e.g. the fluorines of a
      PF6- on top of each other). ``harvest --no-structure-check`` turns it off.
    * ``store_is_writable()``: a job writes the store only on a local filesystem
      (SQLite is not safe on NFS or GPFS); ``fragment-store-writable`` in seamm.ini's
      ``[thermochemistry]`` section overrides it.
    * A JobServer ``job_finished`` hook harvests a job that the step could not write
      back (it leaves ``harvest_fragments.json``) when the job comes back to the
      JobServer's host. Needs seamm-jobserver 2026.11.11.

2026.11.11 -- Model chemistries stored with their program; lookups by several fingerprints
    * The fragment store keeps a model chemistry as the bare model and the program
      (``CMP@R2SCAN-3C`` and ``ORCA``), however it was given. The ORCA sub-steps
      name a level ``CMP@R2SCAN-3C`` and the Model Chemistry step
      ``ORCA:CMP@R2SCAN-3C``, so the two did not match by name before. Lookups
      accept either spelling and also find rows written before this change. The
      harvest now fills in the program of the fragment rows too.
    * The readers (``get``, ``geometries``, ``energy``, ``species_at``) accept
      several fingerprints, so that a caller can look each species up with the
      fingerprints expected for it: ORCA uses exact exchange (NoCOSX) for atoms and
      for Na, Mg, Zn, B and P, which changes their fingerprint at the same level.
    * The command line shows a level as ``PROGRAM:model``, as before.
    * The design of the next phases of the fragment store (write-back from jobs,
      official and local stores, automatic ATcT download) is in the developer
      guide.

2026.11.8.3 -- ZPE scale factors with their sources
    * ``seamm_thermochemistry.zpe_scale_factor(model_chemistry)`` gives the ZPE scale
      factor for a geometry/frequency level, with its source and any caution about
      its use: r2SCAN-3c 0.9899 (fitted to 22 W4-2006 anharmonic ZPVEs, molecules of
      up to 6 atoms), PBE0/def2-TZVP 0.975 (Truhlar's PBE0/MG3S value, basis-
      mismatched) and B2PLYP-D3/cc-pVTZ 0.9863 (CBH-ANL).

2026.11.8.2 -- Thermal corrections from the species map; job ids read again
    * The species map of ``harvest --species-map`` can give each structure's thermal
      corrections (E, ZPE, H at 298.15 K and the symmetry number), for jobs that
      stored only energies, so that their structures get geometry rows. They are used
      only when the map's energy agrees with the job's.
    * Bugfix: harvest read no job id from job_data.json files that start with SEAMM's
      "!MolSSI job_data" header line, so a species map naming the job matched
      nothing and the provenance lacked the job id.

2026.11.8.1 -- Harvesting jobs through a species map
    * ``seamm-thermochemistry-fragments harvest --species-map map.csv`` finds the
      systems of a job by name in a CSV map (job, system name, species; optionally the
      configuration, charge, multiplicity and conformer label) instead of by the
      ``CBH:`` prefix, for jobs whose systems are named otherwise, e.g. by ATcT id. The
      map gives the species with the stereo of its structure.
    * ``--tables geometry`` fills only the geometry and energy tables, so that several
      conformers of a species can be stored.

2026.11.8 -- Geometries and energies at different levels in the fragment store
    * The fragment store has two new tables: ``geometry``, one row per conformer of a
      species at a geometry/frequency level (coordinates, energy, unscaled ZPE and
      H(298.15 K), and the frequencies, symmetry number and moments of inertia when
      the job stored them), and ``energy``, one row per single point on a geometry,
      keyed by the geometry's hash and the calculation's fingerprint. An energy at a
      costly level is combined with the ZPE and thermal correction of a cheaper
      level's geometry when used, never stored combined.
    * ``harvest`` fills them from the same properties as before, refusing an energy
      computed on a different geometry and a geometry optimized at another level, and
      the ladder reads them first. Geometries with imaginary modes are kept but not
      reused.
    * ``seamm-thermochemistry-fragments geometries SMILES`` lists a species'
      geometries and the energies on each. Existing stores gain the tables the next
      time they are written; readers of an older store see none.

2026.11.4 -- Bugfix: ATcT lookup missed ortho/para positional isomers
    * ATcT labels positional isomers "ortho" and "para" as well as nuclear-spin isomers.
      The lookup dropped every such entry as a nuclear-spin isomer, so o- and
      p-tolualdehyde, o-benzyne, the o- and p-chlorophenyl series and similar species
      could not be found, and the error listed no candidates. Only species with at most
      three heavy atoms (H2, H2O, H2O2, CH2O, ketene ...) are now treated as nuclear-spin
      isomers; when a molecule has only such entries, the error says so and lists them.

2026.10.31 -- Superseding values within the duplicate tolerance
    * ``FragmentStore.supersede(record, reason, force=True)`` and
      ``seamm-thermochemistry-fragments harvest ... --supersede "REASON" --force``
      also replace values that differ by less than the 0.1 kJ/mol duplicate tolerance,
      so that a correction applied to a whole level of theory, such as a dispersion
      correction, reaches every record. The old values and the reason are kept in the
      superseded table as before; identical values are left alone, and ``--force``
      needs ``--supersede``.

2026.10.29.1 -- Superseding a stored fragment value
    * ``seamm-thermochemistry-fragments harvest --supersede "REASON"`` replaces stored
      values that differ (a fragment first harvested in the wrong conformer, say)
      instead of only logging them. The old values move to a new ``superseded`` table
      with the date and the reason, and the record keeps its key, so lookups find the
      new values. ``seamm-thermochemistry-fragments superseded [SMILES]`` lists them,
      and ``FragmentStore.supersede(record, reason)`` does it from Python.
    * Stores made by earlier versions gain the table when next opened for writing and
      read as having nothing superseded until then.
    * The documentation now says that the store keeps ZPEs unscaled; a ZPE scale
      factor is applied where the values are used.

2026.10.29 -- Bugfix: ring-closure atoms in CBH fragments
    * A cut atom bonded to two kept atoms (a ring closure, so rings opened with
      ``protect="none"`` or small-ring/aromatic protection off) capped only the first of
      them; the second lost a hydrogen. RDKit refilled it silently for an ordinary atom,
      but a stereocentre (a system's SMILES carries stereo tags) became a radical and the
      reaction failed to balance (propylene carbonate, FEC, bicyclo[1.1.0]butane). Both
      kept atoms are now capped, and a former stereocentre that ends up with two or more
      hydrogens loses its stereo tag.

2026.10.28.4 -- Lowest conformer for CBH species; harvest of every system
    * ``cbh.lowest_conformer(smiles, n)``: the lowest MMFF94 conformer of n RDKit
      embeddings, the starting structure for a species whose lowest rotamer is wanted
      (all-anti n-alkanes, anti ethanol), as CBH-ANL starts every species from its
      lowest rotamer; a single embedding of pentane to octane is gauche.
    * ``seamm-thermochemistry-fragments harvest --all-systems``: every system, its name
      (less the ``CBH:`` prefix) being the species, canonicalized; for jobs whose
      systems Read Structure named by SMILES.

2026.10.28.3 -- The fragment-energy store and CBH-ANL laddered reference values
    * ``seamm-thermochemistry-fragments harvest|list|lookup|duplicates``: a store of
      the energies of CBH species computed earlier, one SQLite file next to the ATcT
      cache, so that a species is computed once per level of theory. A record is keyed
      by the species (canonical SMILES, charge, multiplicity) and the opaque settings
      fingerprints of the energy and frequency calculations, with E, ZPE and
      H(298.15 K) as separate fields; a record without a fingerprint is keyed on the
      model chemistry and code version and flagged unverified. ``harvest`` is the single
      writer, run on the host: insert-or-ignore with provenance, and an attempt whose
      values differ by more than 0.1 kJ/mol is logged, never overwritten.
    * ``seamm-thermochemistry-ladder build|show|list``: laddered reference values as in
      CBH-ANL (Elliott et al. 2023): ATcT only for the CBH-0 species, every other
      species from its own CBH reaction one rung below its size at a high level, with
      the coefficients on the ATcT anchors (so a reaction's ATcT uncertainty is cᵀΣc
      exactly), the provenance, and the difference from ATcT for every ladder species
      ATcT has, whose RMS is the measured method term. Saved as a named reference set.
    * ``seamm-thermochemistry-atct covariance [--from FILE]`` caches a covariance block
      saved from a browser or fetched from the API. ATcT lists its uncertainties
      rounded up to 0.01 kJ/mol; the check of a block allows for that.
    * Bugfix: an entry lookup with an explicit multiplicity tripped over the
      nuclear-spin isomers of the matching state ("para singlet" formaldehyde).
    * The alkane drift golden case now checks the covariance uncertainties: the CBH-2
      drift of C4-C7 referenced to ATcT ethane and propane is outside the 95 % limits.

2026.10.28.2 -- The Active Thermochemical Tables (ATcT) and CBH reactions
    * ``seamm_thermochemistry.atct`` keeps a complete, unmodified copy of one version of
      ATcT (the ``/api/v1/all/`` download) under
      ``Parameters/thermochemistry/ATcT/<version>/``, with its download date and
      SHA-256, builds an RDKit index from canonical SMILES to the entries, and selects
      the entry for a molecule: gas phase and matching charge; double-bond stereo kept,
      a single stereocentre ignored; the plain ``*0`` entry over conformers unless a
      label is asked for; a labelled spin state when the multiplicity matches;
      nuclear-spin isomers never. An ambiguous match stops with the candidates listed,
      and an explicit ATcT ID overrides the rules.
    * Covariance blocks from ``/covariance/species/`` are cached per version, checked
      against the listed 95 % uncertainties, and used by ``reaction_uncertainty`` for
      the 95 % uncertainty of a reaction enthalpy; without a block the uncertainties
      are combined in quadrature and the caller is told.
    * ``seamm-thermochemistry-atct download|list|lookup|index|covariance``;
      ``download --from`` and ``covariance --from`` import copies saved from a
      browser. ATcT lists its uncertainties rounded up to 0.01 kJ/mol, and the check
      of a covariance block allows for that. ATcT data are never distributed with
      SEAMM.
    * ``seamm_thermochemistry.cbh`` generates the connectivity-based hierarchy (CBH)
      reactions of a molecule, CBH-0 to CBH-4, from its SMILES: induced fragments
      capped with hydrogens to the parent valence, keeping bond orders, charges,
      unpaired electrons, double-bond stereo (ring double bonds are cis) and the
      relative configuration of stereocentres; identical fragments merged and both
      sides cancelled; a strict balance of elements, charge, spin and heavy-atom bond
      types; trivial rungs detected. Aromatic ring systems, 3- and 4-membered rings and
      charge-separated or hypervalent groups (nitro, sulfonyl, phosphoryl, PF6-, ...)
      are contracted into single nodes by default ("ring-preserved" CBH); a rung that
      the formula cannot build (small unprotected rings) is reported as unavailable.
      Golden tests reproduce every reaction published by Ramabhadran & Raghavachari
      (2011) and Raghavachari et al. (2023) that we checked: cyclohexanone,
      cyclopentadiene, 2-oxazoline, heptanethiol, propyl pent-4-enoate, methionine
      CBH-0 to CBH-4, methyl 3-butenoate, methyl butanoate and its radical, and
      4-chlorotoluene both ring-preserved and Kekulé. ``unique_fragments`` gives the
      union of the reference species of several reactions, so that a batch computes
      each fragment once.
    * ``scripts/cbh_drift.py`` compares CBH-n enthalpies of formation with ATcT's own
      values for a set of targets, given their 0 K enthalpies at one level of theory.
      A golden case (run only where ATcT ver. 1.220 is cached) freezes the drift of
      the n-alkanes C4-C7 at CBH-2 with ATcT references for ethane and propane,
      +0.38, +1.18, +1.61 and +1.72 kJ/mol against ATcT's own values, from the
      CBH-ANL energies of Elliott et al. (2023). With the ATcT covariance block the
      95 % uncertainty of the drift is 0.18, 0.37, 0.47 and 0.63 kJ/mol (the alkane
      values are positively correlated), so the drift is significant for all four.
    * Found on the real ver. 1.220 data: H2, HD and D2 share the SMILES ``[H][H]``
      (resolved by mass), and the 20 exact element entries have no unit (kJ/mol).
    * RDKit is an optional dependency: ``pip install seamm-thermochemistry[cbh]``.
    * CI now installs the package with uv from its declared dependencies; the conda
      test environment file is gone.

2026.10.28.1 -- Bugfix: enthalpy and Gibbs energy of formation at T
    * ``formation_enthalpy`` and ``formation_gibbs_energy`` subtracted (5/2)RT per
      atom, the gas atom's own thermal enthalpy, where the Hess cycle needs the
      standard-state element's H(T)-H(0 K). The gas atom's thermal enthalpy
      cancels exactly between the two halves of the cycle. The result was too
      negative by 1.96 kJ/mol per H, 5.15 per C and 1.86 per O at 298 K: 5.8 kJ/mol
      for water, and about 5 kJ/mol per carbon atom for organic molecules. Every
      ``DfHT`` and ``DfGT`` reported by orca_step so far has this error; ``DfE0`` and
      the atomization energies were never affected.
    * The correction now uses the tabulated standard-state H(298.15 K)-H(0 K) per
      atom (halved for the diatomic gases, as ``derive_dfH0_0K`` already did). It is
      exact at 298.15 K and scaled linearly in T elsewhere, which is good to about
      1 kJ/mol per atom for gaseous standard states up to 1000 K but poor for
      solids such as graphite far from 298 K. The report states the values used.
    * The tests that encoded the wrong physics were replaced, including a check
      against JANAF data that ΔfH(298 K)-ΔfH(0 K) of water is -2.90 kJ/mol.

2026.10.28 -- Bugfix: ORCA atom runs use node-local scratch
    * ``scripts/orca_atom_multistart.py`` now runs ORCA in a scratch directory under
      ``$TMPDIR`` (node-local on a cluster) and keeps only ``orca.out`` and
      ``orca.gbw`` in the run directories. With hundreds of runs at once, ORCA's
      temporary files on the shared filesystem made the runs I/O-bound, about 18
      times slower, so many hit the 4 h limit per run.
    * A rerun of an unfinished run no longer starts from the orbitals the earlier
      attempt left behind (ORCA's AutoStart picked up the stale ``orca.gbw``): each
      run starts from its own guess, as its start name says.
    * The smeared seed used from La on (PBE/def2-SV(P)) now has SmearTemp 10000,
      VerySlowConv, a permanent level shift of 0.5 Eh and up to 1500 cycles. With the
      old settings it never converged for Tm and Lu, so those searches had no
      smeared starts; now it converges in under a minute, and for Tm it found a
      solution 1.3 kJ/mol lower in def2-TZVP. Seeds that already ran keep their
      input.

2026.9.28 -- Bugfix: updating never replaces a database with local changes
    * ``seamm-thermochemistry-installer update`` -- run by every
      ``seamm-manager update`` -- downloaded the published database over the local
      one, so curated changes were lost. It now replaces the local database only if
      it is exactly one of the published versions (or the copy the installer
      installed). A database with local changes is kept and reported;
      ``update --force`` replaces it anyway.
    * Every replacement first saves the old file as ``thermochemistry.db.bak-<date>``,
      and the download is verified before it replaces anything. Before, a download
      with the wrong checksum deleted the database itself.
    * ``install`` no longer overwrites an existing database, and ``uninstall`` leaves
      one with local changes in place unless given ``--force``.

2026.9.27.1 -- Bugfix: a trial installation cannot replace the shared database
    * In an installation that shares the default installation's codes and data (the
      default for any SEAMM installation other than ``~/SEAMM``, from seamm-manager
      2026.9.27.6), the installer no longer downloads the reference database over the
      shared file, nor deletes it on uninstall; it reports which database the
      installation uses, or that the default installation has none.

2026.9.27 -- The reference database can belong to the installation
    * The default database is now the SEAMM installation's own
      ``Parameters/thermochemistry/thermochemistry.db`` if it has one, then the
      installer-registered database (shared by every installation of the user), then
      the bundled snapshot. A second installation such as ``~/SEAMM_DEV`` therefore
      uses the shared database unless it has its own. Requires seamm-util 2026.9.27.1.
    * Internal: the ORCA atom-energy script uses a smeared start for Z >= 57, resumes
      from compressed outputs, stops its ORCA runs when it is stopped, and redoes runs
      that were killed.

2026.9.25 -- Internal: robust ORCA atom-energy script; no copied Psi4 data
    * Added ``scripts/orca_atom_multistart.py``, which computes ORCA atomic
      reference energies as the lowest-energy SCF solution at the experimental
      multiplicity, from several SCF starts per element, method and basis. A
      single start was found to land in higher solutions (e.g. Bi def2-TZVPD
      24.8 kJ/mol too high, Fe 3d7 4s1 instead of 3d6 4s2). Spin contamination
      is flagged rather than rejected, and ``--rechoose`` re-derives the
      choices from a finished run.
    * ``psi4_step``'s atom-energy CSV is no longer imported when building the
      database: it is ``gaussian_step``'s file with two columns removed, so its
      "psi4" rows were Gaussian numbers. Psi4 atom energies will be added once
      they are computed with Psi4.
    * The README documents the atom reference-energy convention and how the
      ORCA energies are computed.

2026.7.30.1 -- Temperature-dependent formation energies, a shared report, and DOI provenance
    * Added ``formation_enthalpy``/``formation_gibbs_energy``, generalizing
      ``formation_energy`` from the electronic-only, 0 K quantity to a real
      temperature T -- ``DfH(T)`` and ``DfG(T)``, needed to report a proper
      enthalpy/Gibbs energy of formation from a frequency calculation instead of
      the previous ad hoc "0 K + ZPE" shortcut.
    * Added ``ThermoDB.s298_std_state``/``set_s298_std_state`` (schema v2): each
      element's standard-state-phase entropy, needed for ``formation_gibbs_energy``
      -- populated for H and O so far (cited to NIST-JANAF); other elements still
      need the same curation as the rest of the experimental reference data.
    * Added ``format_report``, a shared, citable, per-atom-breakdown text report
      (atomization energy, DfE0, DfH(T), DfG(T), each with the atomic reference
      energies/citations used) -- the one shared replacement for the detailed
      report every consuming code used to duplicate on its own.
    * Added ``ThermoDB.doi``/``set_doi``: a database file can now carry the Zenodo
      DOI of the snapshot it is, so ``format_report`` cites a permanent, citable
      reference instead of just a schema version. ``scripts/publish_to_zenodo.py``
      stamps this in automatically, right before uploading.

2026.7.30 -- Initial release: reference-energy library and Zenodo-backed database
    * Added ``ThermoDB``, a shared SQLite-backed atomic reference-energy database,
      and ``formation_energy``/``atomization_energy`` helpers that turn a raw
      Gaussian, ORCA, Psi4, or VASP total energy into a physically meaningful
      energy or enthalpy of formation, replacing the separate, largely-duplicated
      copies previously carried by ``gaussian_step``, ``psi4_step``, and
      ``vasp_step``.
    * Added the ``seamm-thermochemistry-installer`` console script, which fetches
      the published reference database from Zenodo (checksum-verified) instead of
      bundling it in the package.
    * Added the ``seamm-thermochemistry-import-orca`` console script for importing
      ORCA atom-energy scans, with automatic spin (S^2) sanity checks and a
      configurable policy for resolving duplicate/conflicting entries.
