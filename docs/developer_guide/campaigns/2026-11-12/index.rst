2026-11-12 -- The supplementary reference library (R1)
=======================================================

Status: **R1 implemented 2026-11-12** (branch ``reference-library``): the library
file, the curated-CSV reader and the NIST WebBook importer. No data is bundled or
published. The step does not use the library yet (R2). Publishing (R3) waits for
Paul's decision on the NIST terms (SRD 69 for the WebBook, SRD 13 for JANAF).

Background:

- the CBH reference-library design of 2026-10-29 (§2 resolution order, §3 the
  library, §6 user references, §8 licensing), and its addendum of 2026-11-11;
- the reference-data pilot of 2026-11-12, Paul's decision 8. Thermo curated the
  references ATcT lacks for the PPE1694 set: 62 species and 70 of its 123 targets
  can be closed from existing data.

The library holds enthalpies of formation that ATcT lacks, with every value a
source gives, the one chosen and why.

The files
---------

``references-official.db``
    Published on Zenodo, installed beside the ATcT cache and the fragment stores
    (``<...>/Parameters/thermochemistry/``), and never written locally.
``references.db``
    The local layer: user references and computed ones not yet published.

:class:`~seamm_thermochemistry.references.LayeredReferences` reads both. Within a
rank the local layer comes first. Each file carries a schema version, and a file
from newer code is never written, as for the fragment store.

The tables:

``source``
    Name, kind, rank (the resolution order), citation, licence and when it was
    checked, and the offset against ATcT (rms, number of shared species, ATcT
    version).
``species``
    Canonical SMILES (the CBH key), charge, multiplicity, Hill formula, CAS,
    InChI, and the non-stereo InChI when that was needed to find the species.
``measurement``
    Every value a source gives, chosen or not:

    - ΔfH at 298.15 K and/or 0 K;
    - the uncertainty **as stated**, with its level;
    - status (measured, review, estimate, derived, held, computed, user);
    - year, H298 − H0, standard state, ion convention;
    - the original citation, method, comment, and where it was retrieved.
``reference``
    The value chosen per species and source, with the rule that chose it, the
    spread, a flag, the value it overrides and why, and any human accept (who, when,
    why).
``derivation``
    For a derived value: the reaction, its partners and their sources, and ΔrH.

The rules
---------

- **Burcat is never a source.** A source named Burcat, or a value citing it, is
  refused.
- **No P species from the online JANAF tables.** Their PH3 and PH2 sit ~17.5 kJ/mol
  off the printed edition: early JANAF used red P. P species come from the printed
  edition.
- **Estimates are never chosen automatically.** They are stored, and used only after
  an accept, which records who, when and why.
- **The replacement rule.** A later source replaces an earlier one only when its
  95 % uncertainty is below a third of the earlier one's **and** the two agree within
  the earlier one's uncertainty. Both values are kept, with the reason. In the
  pilot, D0 replaces JANAF for NaH and KH.
- **Derived values are checked.** Each derived value has its derivation, and
  ``check`` recomputes it from its partners' current values. JANAF's tabulated PCl5
  (−360.18) fails; the derived −376.08 passes.
- **Ions are stored in ATcT's stationary-electron convention.** A value in the
  thermal-electron convention is converted on entry, by −6.197 kJ/mol per unit of
  charge at 298.15 K (JANAF H⁺ 1536.25 → 1530.05, as ATcT). An ion without its
  convention is refused.
- **Choosing among one source's measurements** (thermo's rule):

  1. a review value;
  2. else the most recent with a stated uncertainty, ties going to the smaller;
  3. else the most recent.

  Values derived from another on the page (WebBook "N/A") come last, and
  evaluations the WebBook shows ("cal" rows) are held.
- **Spreads** (science's Q3):

  - flagged when another measurement differs from the chosen one by more than
    their combined uncertainties;
  - **held** (not chosen) when that exceeds three times the combined uncertainties
    and there is no review value, until a human accepts one.
- **0 K-only values** (W4-17) are used only when no source gives a value at
  298.15 K.

Uncertainties (science's Q2):

- 95 % and 2σ values as stated; 1σ and unstated levels doubled;
- a measured value with no uncertainty gets 4 kJ/mol;
- W4-17 combines the atoms (2σ) with twice its 1 kJ/mol method error.

The source's offset, the rms of (source − ATcT) over at least 5 shared species, is
added in quadrature and kept in its own column.

The command line
----------------

::

    seamm-thermochemistry-references sources | list | check [--library FILE]
    seamm-thermochemistry-references show SMILES [--local FILE] [--official FILE]
    seamm-thermochemistry-references accept MEASUREMENT --by WHO --reason TEXT
    seamm-thermochemistry-references add-curated DIR --library FILE
    seamm-thermochemistry-references webbook fetch SPECIES.csv --cache DIR
    seamm-thermochemistry-references webbook import --cache DIR --library FILE
    seamm-thermochemistry-references offsets --library FILE [--atct-version V]

The importers (``add-curated``, ``webbook``, ``offsets``) are maintainer tools that
build a library file for publication. Users only install the published file, and
never contact NIST.

``webbook fetch`` looks species up by InChI and saves the pages, 3 s apart, logging
each request. A stereo InChI the WebBook does not find is retried without its stereo
layers; the retry is recorded and noted as "stereoisomer not verified".
``webbook import`` then reads the saved pages offline.

Acceptance (R1d)
----------------

A library was built from thermo's curated CSVs and the pilot's cached WebBook pages
and compared with the pilot's choices. The result is 58 species and 65 targets.

Four species are held under the spread rule: imidazole, isoxazole, and 3- and
5-methylisoxazole. Their measurements differ by 3.2–6.2 times their combined
uncertainties, and none has a review value. Accepting the pilot's choice for them
gives exactly the pilot's **62 species and 70 targets**, with identical values.

Science's decision (2026-11-12): the two methylisoxazoles are accepted (their
conflict is only with Tappi 1940). Isoxazole stays held: Steele et al. 1996, 82.0 ±
0.6, against McCormick & Hamilton 1978, 78.58 ± 0.54, both modern. Imidazole stays
held and goes on the tier-A compute list. The accepts are kept in a CSV and
reapplied with ``accept-file`` whenever the library is rebuilt. Result: **60 species
and 65 targets**. Isoxazole blocks 4 targets and imidazole 1.

The one other difference is HCO⁺: 826.80 kJ/mol against the pilot's 833.0. That is
the conversion of JANAF's thermal-electron value to ATcT's convention, as intended.

The tests contain no WebBook data. They use a synthetic page with the WebBook's
markup and made-up values.

Offsets: a finding
------------------

The pilot's species are, by construction, the ones ATcT lacks, so the library shares
almost nothing with ATcT:

- CODATA shares 10 atoms, with an rms of 2.39 kJ/mol, dominated by B (−5.4) and
  Si (−5.0), whose atoms ATcT has revised;
- every other source shares none.

Measuring the offsets as designed needs a **calibration set**: values from each
source for species ATcT does have. These are not used as references, since ATcT
wins for them, but they measure each source. That is data work for thermo, and the
WebBook part needs more NIST requests, so it is science's call.

Next
----

- **R2:** the resolver in seamm_thermochemistry (resolution order, one source per
  element, the uncertainties) and its use in thermochemistry_step, coordinated with
  Design. Science's Q1: a species whose override gives an element a second source is
  labelled "mixed: several sources", with the sources named, and the step prints one
  line about it, never blocking.
- **R3:** publishing ``references-official.db`` and installing it, after Paul's
  decision on the NIST terms.
