2026-11-12 -- The supplementary reference library (R1)
=======================================================

Status: **R1, R2 and R3 done 2026-11-12** (seamm-thermochemistry 2026.11.12.1,
thermochemistry-step 2026.11.12.2): the library file, the curated-CSV reader, the
NIST WebBook importer and the resolver. Paul ruled on 2026-11-12 that the
WebBook, JANAF and NBS-82 values may be used and shipped, each citing its
original source and crediting the compilation. The official library is on Zenodo
(10.5281/zenodo.24175074), installed by seamm-thermochemistry 2026.11.12.2.

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

**The sources, in rank order** (science, 2026-11-12): printed JANAF, derived
values, CODATA, the NBS tables of 1982 (NBS-82), D0, W4-17, the WebBook. NBS-82
(Wagman et al., JPCRD 11, Suppl. 2) was added the same day for six Ge, As and Se
species. It states no uncertainty, so its values get the 4 kJ/mol floor, and its
H298 − H0 is kept with each value.

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
held and goes on the tier-A compute list. The decisions (accepts and holds, with the reasons) are kept in a CSV and
reapplied with ``accept-file`` whenever the library is rebuilt. Result: **60 species
and 65 targets**. Isoxazole blocks 4 targets and imidazole 1.

Later the same day thermo added NBS-82 (six Ge, As and Se species) and science held
POF3, whose JANAF (−1254.25) and NBS-82 (−1211.3) values disagree by 43 kJ/mol,
pending a computed arbitration. A hold may name several sources ("JANAF-1998 /
NBS-82"); a source with no value, or whose file already holds it, is passed over.
Result: **65 species and 70 targets**, as thermo's own count.

The one other difference is HCO⁺: 826.80 kJ/mol against the pilot's 833.0. That is
the conversion of JANAF's thermal-electron value to ATcT's convention, as intended.

The tests contain no WebBook data. They use a synthetic page with the WebBook's
markup and made-up values.

Offsets against ATcT
--------------------

The pilot's species are, by construction, ones ATcT lacks, so the library itself
shares almost nothing with ATcT. Thermo therefore supplied a **calibration set**
(``add_calibration``; ``offsets --calibration``). It holds the values each source
gives for species ATcT 1.228 has, and these are used only to measure the source,
never as references.

Offsets with ATcT 1.228:

- JANAF (printed): n = 10, rms 0.38 kJ/mol;
- WebBook: n = 9, rms 1.16 kJ/mol. The largest difference is methyl acetate,
  +2.98, from a value with no stated uncertainty;
- CODATA: n = 8 atoms from the library itself, rms 0.69 kJ/mol. B and Si are
  excluded, as a recorded exclusion ("CODATA key value known superseded: B, Si";
  ``offsets --exclude``), since ATcT has revised both. Over all 10 atoms it is
  2.39, which is kept beside it;
- NBS-82: n = 10, rms 2.81 kJ/mol (mean +1.02), over the same species as the
  JANAF set. CF4 (+8.6) is a noted outlier, not an exclusion (science): NBS-82
  states no uncertainties and is the oldest evaluation used, and its obscure
  inorganics are where an unrevised value like CF4's would hide. With the
  4 kJ/mol floor its values carry 4.89 kJ/mol (GeH4 90.80 ± 4.89);
- D0, W4-17 and the derived values: under 5 shared species, so no offset.

Each offset is added in quadrature to the source's 95 % uncertainties, and is kept
in its own column so it can be shown and changed.

A value derived from others inherits its partners' offsets through the
derivation, sqrt(Σ (ν_i offset_i)²) (science, 2026-11-12). PCl5, derived from
JANAF's PCl3 and JANAF's ΔrH, inherits JANAF's 0.38. A D0 value with its atoms
recorded as partners inherits CODATA's offset for the metal atom only. Partners
resolve as the step would: the H atom comes from ATcT and adds nothing. NaH:
138.59 ± 1.55 (95 %, with 0.69 inherited); KH: 132.54 ± 1.06.

**Caveat (science):** the calibration species are well studied, so an offset is a
*lower bound* for the obscure species the library actually serves. It does not
measure their reliability; the stated uncertainties and the spread rule do that.

The resolver (R2a)
------------------

:class:`~seamm_thermochemistry.resolver.ReferenceResolver` resolves a target's
reference species.

**The order.** ATcT comes first. Then the user's own references, if enabled.
Then the library's sources by rank, and computed references, if enabled.
"Supplementary references = none" reproduces today's ATcT-only behaviour and
reads no library.

**One source per element.** For each element ATcT lacks, i.e. whose gas-phase atom
is not in ATcT:

- the references containing it are taken from the first source, by rank, that
  has all of them;
- a replacement (override) inside that source still applies to a single
  species. The element then has two sources: the result is labelled "mixed:
  several sources" with the sources named, and a one-line note says so;
- when no single source has them all, each reference takes its own best value,
  under the same label.

**The values** come back in the layout of an ATcT entry:

- the uncertainty is the 95 % value with the offset added;
- ``atct_coefficients`` is empty, since these values have no covariance with
  ATcT, so they are combined in quadrature;
- the entry carries its source, layer, citation, rule, flag, and any override or
  decision.

**0 K values.** A reference given only at 298.15 K is converted to 0 K with
:func:`~seamm_thermochemistry.resolver.dfh_0k_from_298k`, using the species' own
H298 − H0 (the step computes it) and CODATA's values for the elements.

**The result labels:** ATcT-only, mixed, mixed: several sources, computed
references, user references.

Next
----

- **R2b** is in thermochemistry-step 2026.11.12.2. It covers:

  - the two parameters;
  - the species table's source column;
  - the result label and the ``references`` property;
  - the non-ATcT uncertainty in quadrature;
  - the 0 K conversion;
  - the missing-reference message.

- **R3, done 2026-11-12:** the official library, version 2026.11.12, is published
  on Zenodo: https://zenodo.org/records/24175074 (DOI 10.5281/zenodo.24175074,
  concept 10.5281/zenodo.24175073), under CC-BY-4.0. Each value cites its original
  source, and NIST is credited for JANAF, the NBS tables and the WebBook. The
  record was built with the maintainer commands and published with
  ``scripts/publish_references_to_zenodo.py``, which stamps the DOI and version
  into the file. The installer (seamm-thermochemistry 2026.11.12.2) installs and
  updates it as it does the official fragment store.
