# Test fixtures

`atct_fixture.json` is a **synthetic** file in the layout of the ATcT API's
`/api/v1/all/` download, labelled TN version "9.999". It exists to test the
cache, the index and the entry-selection rules (plain entries, conformers,
cis/trans, mixtures, spin states, nuclear-spin isomers, ions, condensed phases,
an "exact" uncertainty, a SMILES RDKit cannot parse, and a single stereocentre).
The numbers are only roughly those of ATcT ver. 1.220 and must not be used as
reference data. Real ATcT data are never distributed with SEAMM; each
installation downloads its own copy.
