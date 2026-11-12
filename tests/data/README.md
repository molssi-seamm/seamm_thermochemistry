# Test fixtures

`atct_fixture.json` is a **synthetic** file in the layout of the ATcT API's
`/api/v1/all/` download, labelled TN version "9.999". It exists to test the
cache, the index and the entry-selection rules (plain entries, conformers,
cis/trans, mixtures, spin states, nuclear-spin isomers, ions, condensed phases,
an "exact" uncertainty, a SMILES RDKit cannot parse, and a single stereocentre).
The numbers are only roughly those of ATcT ver. 1.220 and must not be used as
reference data. Real ATcT data are never distributed with SEAMM; each
installation downloads its own copy.

`curated/` holds a few **real** rows of the hand-curated reference data
(2026-11-12), so that the tests record what the library's rules decide for real
cases. Each row cites its original source, with table or page:

- printed JANAF: M. W. Chase, NIST-JANAF Thermochemical Tables, 4th ed., J. Phys.
  Chem. Ref. Data Monograph 9 (1998), NIST Standard Reference Database 13;
- CODATA: J. D. Cox, D. D. Wagman and V. A. Medvedev, CODATA Key Values for
  Thermodynamics (Hemisphere, 1989);
- spectroscopic D0 values and the W4-17 set (A. Karton, N. Sylvetsky and
  J. M. L. Martin, J. Comput. Chem. 38, 2063 (2017)), as cited in each row;
- the derived PCl5 value, from JANAF's PCl3 and reaction enthalpy.

`webbook/make_page.py` writes a **synthetic** page in the NIST Chemistry
WebBook's markup with made-up values; no WebBook data are in the tests.
