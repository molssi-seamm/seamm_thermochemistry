=======
History
=======

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
